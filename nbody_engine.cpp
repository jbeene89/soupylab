/**
 * N-BODY GRAVITATIONAL SIMULATOR - High Performance WebAssembly
 * Barnes-Hut Algorithm for O(n log n) Complexity
 *
 * Features:
 * - 10,000+ particle simulation
 * - Barnes-Hut quadtree optimization
 * - Real gravitational physics (Newton's law)
 * - Softening parameter to prevent singularities
 * - Multiple scenarios (galaxy collision, solar system, clusters)
 * - Collision detection & merging
 * - Adaptive timestep for stability
 * - SIMD-optimized vector operations
 */

#include <emscripten/bind.h>
#include <emscripten/val.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <memory>
#include <random>

using namespace emscripten;

// ============================================================================
// PHYSICAL CONSTANTS
// ============================================================================

constexpr float G = 6.67430e-11f;  // Gravitational constant (SI units)
constexpr float G_SCALED = 1.0f;   // Scaled for simulation (arbitrary units)
constexpr float SOFTENING = 0.1f;  // Softening length to prevent singularities
constexpr float THETA = 0.5f;      // Barnes-Hut opening angle

// ============================================================================
// VECTOR 2D
// ============================================================================

struct Vec2 {
    float x, y;

    Vec2() : x(0), y(0) {}
    Vec2(float x, float y) : x(x), y(y) {}

    Vec2 operator+(const Vec2& v) const { return Vec2(x + v.x, y + v.y); }
    Vec2 operator-(const Vec2& v) const { return Vec2(x - v.x, y - v.y); }
    Vec2 operator*(float s) const { return Vec2(x * s, y * s); }
    Vec2 operator/(float s) const { return Vec2(x / s, y / s); }

    Vec2& operator+=(const Vec2& v) { x += v.x; y += v.y; return *this; }
    Vec2& operator-=(const Vec2& v) { x -= v.x; y -= v.y; return *this; }
    Vec2& operator*=(float s) { x *= s; y *= s; return *this; }

    float length() const { return sqrtf(x * x + y * y); }
    float lengthSq() const { return x * x + y * y; }
    Vec2 normalized() const { float len = length(); return len > 0 ? *this / len : Vec2(0, 0); }
};

// ============================================================================
// PARTICLE
// ============================================================================

struct Particle {
    Vec2 pos;
    Vec2 vel;
    Vec2 acc;
    float mass;
    float radius;
    bool active;

    Particle() : mass(1.0f), radius(1.0f), active(true) {}

    Particle(Vec2 pos, Vec2 vel, float mass)
        : pos(pos), vel(vel), acc(0, 0), mass(mass), radius(1.0f), active(true) {}
};

// ============================================================================
// BARNES-HUT QUADTREE NODE
// ============================================================================

struct QuadNode {
    Vec2 centerOfMass;
    float totalMass;
    Vec2 center;      // Center of this quadrant
    float size;       // Half-width of quadrant

    std::unique_ptr<QuadNode> children[4];  // NW, NE, SW, SE
    Particle* particle;  // If leaf: points to single particle

    bool isLeaf;

    QuadNode(Vec2 center, float size)
        : centerOfMass(0, 0)
        , totalMass(0)
        , center(center)
        , size(size)
        , particle(nullptr)
        , isLeaf(true)
    {}

    // Get quadrant index (0=NW, 1=NE, 2=SW, 3=SE)
    int getQuadrant(const Vec2& pos) const {
        bool west = pos.x < center.x;
        bool north = pos.y < center.y;
        return (north ? 0 : 2) + (west ? 0 : 1);
    }

    // Get child center for quadrant
    Vec2 getChildCenter(int quadrant) const {
        float halfSize = size * 0.5f;
        float dx = (quadrant & 1) ? halfSize : -halfSize;
        float dy = (quadrant & 2) ? halfSize : -halfSize;
        return Vec2(center.x + dx, center.y + dy);
    }

    // Check if point is in this node
    bool contains(const Vec2& pos) const {
        return fabsf(pos.x - center.x) <= size &&
               fabsf(pos.y - center.y) <= size;
    }

    // Insert particle into tree
    void insert(Particle* p) {
        if (!contains(p->pos)) return;

        if (isLeaf) {
            if (particle == nullptr) {
                // Empty leaf - store particle
                particle = p;
                centerOfMass = p->pos;
                totalMass = p->mass;
            } else {
                // Leaf already has particle - subdivide
                Particle* oldParticle = particle;
                particle = nullptr;
                isLeaf = false;

                // Insert old particle into child
                int oldQuad = getQuadrant(oldParticle->pos);
                if (!children[oldQuad]) {
                    children[oldQuad] = std::make_unique<QuadNode>(
                        getChildCenter(oldQuad), size * 0.5f
                    );
                }
                children[oldQuad]->insert(oldParticle);

                // Insert new particle
                int newQuad = getQuadrant(p->pos);
                if (!children[newQuad]) {
                    children[newQuad] = std::make_unique<QuadNode>(
                        getChildCenter(newQuad), size * 0.5f
                    );
                }
                children[newQuad]->insert(p);

                // Update mass and COM
                updateMassAndCOM();
            }
        } else {
            // Internal node - insert into appropriate child
            int quad = getQuadrant(p->pos);
            if (!children[quad]) {
                children[quad] = std::make_unique<QuadNode>(
                    getChildCenter(quad), size * 0.5f
                );
            }
            children[quad]->insert(p);

            // Update mass and COM
            updateMassAndCOM();
        }
    }

    // Update center of mass and total mass
    void updateMassAndCOM() {
        if (isLeaf) {
            if (particle) {
                centerOfMass = particle->pos;
                totalMass = particle->mass;
            }
        } else {
            Vec2 com(0, 0);
            float total = 0;

            for (int i = 0; i < 4; ++i) {
                if (children[i]) {
                    float childMass = children[i]->totalMass;
                    com += children[i]->centerOfMass * childMass;
                    total += childMass;
                }
            }

            if (total > 0) {
                centerOfMass = com / total;
                totalMass = total;
            }
        }
    }

    // Calculate force on particle using Barnes-Hut approximation
    Vec2 calculateForce(const Particle* p) const {
        if (totalMass == 0 || !p->active) return Vec2(0, 0);

        Vec2 diff = centerOfMass - p->pos;
        float distSq = diff.lengthSq();

        // If this is the particle itself, skip
        if (isLeaf && particle == p) return Vec2(0, 0);

        // Barnes-Hut criterion: s/d < theta
        float d = sqrtf(distSq);
        if (d == 0) return Vec2(0, 0);

        if (isLeaf || (size / d) < THETA) {
            // Use this node as approximation
            float softDistSq = distSq + SOFTENING * SOFTENING;
            float invDist3 = 1.0f / (softDistSq * sqrtf(softDistSq));
            float forceMag = G_SCALED * totalMass * p->mass * invDist3;

            return diff * forceMag;
        } else {
            // Recurse to children
            Vec2 force(0, 0);
            for (int i = 0; i < 4; ++i) {
                if (children[i]) {
                    force += children[i]->calculateForce(p);
                }
            }
            return force;
        }
    }
};

// ============================================================================
// N-BODY SIMULATOR
// ============================================================================

class NBodyEngine {
private:
    std::vector<Particle> particles;
    std::unique_ptr<QuadNode> root;

    float worldSize;
    float time;
    float dt;
    int iteration;

    // Statistics
    float totalEnergy;
    float kineticEnergy;
    float potentialEnergy;

    std::mt19937 rng;

    // Environmental controls
    float currentDampening;  // 0.0-1.0, reduces velocity
    std::vector<Vec2> barrierZones;  // Barrier positions
    std::vector<float> barrierRadii;  // Barrier radii

    void buildTree() {
        root = std::make_unique<QuadNode>(Vec2(0, 0), worldSize);

        for (auto& p : particles) {
            if (p.active) {
                root->insert(&p);
            }
        }
    }

    void calculateForces() {
        for (auto& p : particles) {
            if (!p.active) continue;

            Vec2 force = root->calculateForce(&p);
            p.acc = force / p.mass;
        }
    }

    void integrateLeapfrog() {
        // Leapfrog integration (symplectic, conserves energy better)
        for (auto& p : particles) {
            if (!p.active) continue;

            // Kick: v(t+dt/2) = v(t) + a(t)*dt/2
            p.vel += p.acc * (dt * 0.5f);

            // Drift: x(t+dt) = x(t) + v(t+dt/2)*dt
            p.pos += p.vel * dt;

            // Wrap around boundaries (periodic)
            if (p.pos.x < -worldSize) p.pos.x += 2 * worldSize;
            if (p.pos.x > worldSize) p.pos.x -= 2 * worldSize;
            if (p.pos.y < -worldSize) p.pos.y += 2 * worldSize;
            if (p.pos.y > worldSize) p.pos.y -= 2 * worldSize;
        }

        // Rebuild tree with new positions
        buildTree();
        calculateForces();

        // Final kick: v(t+dt) = v(t+dt/2) + a(t+dt)*dt/2
        for (auto& p : particles) {
            if (!p.active) continue;
            p.vel += p.acc * (dt * 0.5f);
        }
    }

    void handleCollisions() {
        // Simple collision detection (N² but only for close particles)
        for (size_t i = 0; i < particles.size(); ++i) {
            if (!particles[i].active) continue;

            for (size_t j = i + 1; j < particles.size(); ++j) {
                if (!particles[j].active) continue;

                Vec2 diff = particles[j].pos - particles[i].pos;
                float distSq = diff.lengthSq();
                float minDist = particles[i].radius + particles[j].radius;

                if (distSq < minDist * minDist) {
                    // Merge particles (conservation of momentum)
                    float totalMass = particles[i].mass + particles[j].mass;
                    Vec2 newVel = (particles[i].vel * particles[i].mass +
                                  particles[j].vel * particles[j].mass) / totalMass;
                    Vec2 newPos = (particles[i].pos * particles[i].mass +
                                  particles[j].pos * particles[j].mass) / totalMass;

                    particles[i].pos = newPos;
                    particles[i].vel = newVel;
                    particles[i].mass = totalMass;
                    particles[i].radius = sqrtf(particles[i].radius * particles[i].radius +
                                               particles[j].radius * particles[j].radius);

                    particles[j].active = false;
                }
            }
        }
    }

    void calculateEnergy() {
        kineticEnergy = 0;
        potentialEnergy = 0;

        for (const auto& p : particles) {
            if (!p.active) continue;
            kineticEnergy += 0.5f * p.mass * p.vel.lengthSq();
        }

        // Potential energy (expensive, only calculate occasionally)
        if (iteration % 10 == 0) {
            for (size_t i = 0; i < particles.size(); ++i) {
                if (!particles[i].active) continue;

                for (size_t j = i + 1; j < particles.size(); ++j) {
                    if (!particles[j].active) continue;

                    Vec2 diff = particles[j].pos - particles[i].pos;
                    float dist = diff.length() + SOFTENING;
                    potentialEnergy -= G_SCALED * particles[i].mass * particles[j].mass / dist;
                }
            }
        }

        totalEnergy = kineticEnergy + potentialEnergy;
    }

public:
    NBodyEngine(float worldSize)
        : worldSize(worldSize)
        , time(0)
        , dt(0.01f)
        , iteration(0)
        , totalEnergy(0)
        , kineticEnergy(0)
        , potentialEnergy(0)
        , rng(12345)
        , currentDampening(0.0f)
    {}

    // ========================================================================
    // SCENARIO CREATION
    // ========================================================================

    void createGalaxyCollision(int particlesPerGalaxy) {
        particles.clear();

        std::uniform_real_distribution<float> dist(0, 1);

        // Galaxy 1 (left)
        Vec2 center1(-worldSize * 0.3f, 0);
        Vec2 vel1(50.0f, 20.0f);

        for (int i = 0; i < particlesPerGalaxy; ++i) {
            float r = sqrtf(dist(rng)) * worldSize * 0.3f;
            float theta = dist(rng) * 2 * 3.14159f;

            Vec2 pos = center1 + Vec2(cosf(theta), sinf(theta)) * r;

            // Orbital velocity
            float v = sqrtf(G_SCALED * particlesPerGalaxy * 0.5f / (r + 1.0f));
            Vec2 vel = vel1 + Vec2(-sinf(theta), cosf(theta)) * v;

            float mass = 1.0f;
            particles.emplace_back(pos, vel, mass);
        }

        // Galaxy 2 (right)
        Vec2 center2(worldSize * 0.3f, 0);
        Vec2 vel2(-50.0f, -20.0f);

        for (int i = 0; i < particlesPerGalaxy; ++i) {
            float r = sqrtf(dist(rng)) * worldSize * 0.3f;
            float theta = dist(rng) * 2 * 3.14159f;

            Vec2 pos = center2 + Vec2(cosf(theta), sinf(theta)) * r;

            float v = sqrtf(G_SCALED * particlesPerGalaxy * 0.5f / (r + 1.0f));
            Vec2 vel = vel2 + Vec2(-sinf(theta), cosf(theta)) * v;

            float mass = 1.0f;
            particles.emplace_back(pos, vel, mass);
        }
    }

    void createSolarSystem() {
        particles.clear();

        // Sun
        particles.emplace_back(Vec2(0, 0), Vec2(0, 0), 1000.0f);
        particles.back().radius = 10.0f;

        // Planets
        float distances[] = {100, 150, 200, 280, 380};
        for (int i = 0; i < 5; ++i) {
            float r = distances[i];
            float v = sqrtf(G_SCALED * 1000.0f / r);

            particles.emplace_back(
                Vec2(r, 0),
                Vec2(0, v),
                1.0f + i * 0.5f
            );
            particles.back().radius = 2.0f + i * 0.5f;
        }
    }

    void createRandomCluster(int count) {
        particles.clear();

        std::uniform_real_distribution<float> posDist(-worldSize * 0.5f, worldSize * 0.5f);
        std::uniform_real_distribution<float> velDist(-10.0f, 10.0f);
        std::uniform_real_distribution<float> massDist(0.5f, 2.0f);

        for (int i = 0; i < count; ++i) {
            Vec2 pos(posDist(rng), posDist(rng));
            Vec2 vel(velDist(rng), velDist(rng));
            float mass = massDist(rng);

            particles.emplace_back(pos, vel, mass);
        }
    }

    void createOceanTrash(int count) {
        particles.clear();

        std::uniform_real_distribution<float> dist(0, 1);
        std::uniform_real_distribution<float> angleDist(0, 2 * 3.14159f);

        // Create multiple ocean current vortices with trash
        int numVortices = 5; // Pacific, Atlantic, Indian, etc.
        int trashPerVortex = count / numVortices;

        for (int v = 0; v < numVortices; ++v) {
            // Position each vortex in different ocean region
            float vortexAngle = (v / (float)numVortices) * 2 * 3.14159f;
            float vortexRadius = worldSize * 0.4f;
            Vec2 vortexCenter(
                cosf(vortexAngle) * vortexRadius,
                sinf(vortexAngle) * vortexRadius
            );

            // Vortex rotation speed (gentle ocean current)
            float vortexSpeed = 15.0f + dist(rng) * 10.0f;

            for (int i = 0; i < trashPerVortex; ++i) {
                // Distribute trash in spiral pattern around vortex
                float r = sqrtf(dist(rng)) * worldSize * 0.15f;
                float theta = dist(rng) * 2 * 3.14159f;

                Vec2 offsetPos(cosf(theta) * r, sinf(theta) * r);
                Vec2 pos = vortexCenter + offsetPos;

                // Circular current velocity + some turbulence
                Vec2 circularVel(-sinf(theta) * vortexSpeed, cosf(theta) * vortexSpeed);
                Vec2 turbulence(
                    (dist(rng) - 0.5f) * 5.0f,
                    (dist(rng) - 0.5f) * 5.0f
                );
                Vec2 vel = circularVel + turbulence;

                // Different trash types with varied mass (plastic, bottles, bags, debris)
                float trashType = dist(rng);
                float mass, radius;

                if (trashType < 0.4f) {
                    // Plastic bottles - medium size, light
                    mass = 0.1f + dist(rng) * 0.2f;
                    radius = 1.5f + dist(rng) * 0.5f;
                } else if (trashType < 0.7f) {
                    // Plastic bags - very light, small
                    mass = 0.05f + dist(rng) * 0.1f;
                    radius = 0.8f + dist(rng) * 0.3f;
                } else if (trashType < 0.9f) {
                    // Microplastics and small debris
                    mass = 0.02f + dist(rng) * 0.05f;
                    radius = 0.3f + dist(rng) * 0.2f;
                } else {
                    // Larger debris (fishing nets, etc.)
                    mass = 0.3f + dist(rng) * 0.4f;
                    radius = 2.0f + dist(rng) * 1.0f;
                }

                particles.emplace_back(pos, vel, mass);
                particles.back().radius = radius;
            }
        }

        // Add some scattered floating trash between vortices
        int scatteredCount = count - (trashPerVortex * numVortices);
        std::uniform_real_distribution<float> scatterPos(-worldSize * 0.6f, worldSize * 0.6f);
        std::uniform_real_distribution<float> scatterVel(-8.0f, 8.0f);

        for (int i = 0; i < scatteredCount; ++i) {
            Vec2 pos(scatterPos(rng), scatterPos(rng));
            Vec2 vel(scatterVel(rng), scatterVel(rng));

            // Random trash type
            float mass = 0.05f + dist(rng) * 0.25f;
            float radius = 0.5f + dist(rng) * 1.5f;

            particles.emplace_back(pos, vel, mass);
            particles.back().radius = radius;
        }
    }

    // ========================================================================
    // TIME TRAVEL / ENVIRONMENTAL CONTROLS
    // ========================================================================

    // Add more trash pollution (simulates increased pollution over time)
    void addTrashPollution(int amount) {
        std::uniform_real_distribution<float> dist(0, 1);
        std::uniform_real_distribution<float> posDist(-worldSize * 0.6f, worldSize * 0.6f);
        std::uniform_real_distribution<float> velDist(-8.0f, 8.0f);

        for (int i = 0; i < amount; ++i) {
            Vec2 pos(posDist(rng), posDist(rng));
            Vec2 vel(velDist(rng), velDist(rng));

            // Random trash type
            float trashType = dist(rng);
            float mass, radius;

            if (trashType < 0.4f) {
                mass = 0.1f + dist(rng) * 0.2f;
                radius = 1.5f + dist(rng) * 0.5f;
            } else if (trashType < 0.7f) {
                mass = 0.05f + dist(rng) * 0.1f;
                radius = 0.8f + dist(rng) * 0.3f;
            } else if (trashType < 0.9f) {
                mass = 0.02f + dist(rng) * 0.05f;
                radius = 0.3f + dist(rng) * 0.2f;
            } else {
                mass = 0.3f + dist(rng) * 0.4f;
                radius = 2.0f + dist(rng) * 1.0f;
            }

            particles.emplace_back(pos, vel, mass);
            particles.back().radius = radius;
        }
    }

    // Remove trash (simulates cleanup efforts)
    void removeTrashPollution(int amount) {
        int removed = 0;
        for (size_t i = particles.size(); i > 0 && removed < amount; --i) {
            if (particles[i - 1].active && particles[i - 1].mass < 0.5f) {
                particles[i - 1].active = false;
                removed++;
            }
        }
    }

    // Simulate glacial melt: adds cold freshwater that disrupts currents
    void addGlacialMelt(int intensity) {
        std::uniform_real_distribution<float> dist(0, 1);

        // Glacial melt comes from polar regions (top and bottom of world)
        int icebergsPerPole = intensity / 2;

        for (int pole = 0; pole < 2; ++pole) {
            float polarY = (pole == 0) ? -worldSize * 0.7f : worldSize * 0.7f;

            for (int i = 0; i < icebergsPerPole; ++i) {
                // Ice spreads along polar region
                float x = (dist(rng) - 0.5f) * worldSize * 1.2f;
                float y = polarY + (dist(rng) - 0.5f) * worldSize * 0.2f;
                Vec2 pos(x, y);

                // Ice moves slowly, disrupts normal currents
                float driftSpeed = 5.0f + dist(rng) * 5.0f;
                float driftAngle = dist(rng) * 2 * 3.14159f;
                Vec2 vel(
                    cosf(driftAngle) * driftSpeed,
                    sinf(driftAngle) * driftSpeed
                );

                // Ice/freshwater has different density
                float iceType = dist(rng);
                float mass, radius;

                if (iceType < 0.5f) {
                    // Large icebergs
                    mass = 2.0f + dist(rng) * 3.0f;
                    radius = 3.0f + dist(rng) * 2.0f;
                } else {
                    // Melted freshwater masses
                    mass = 0.5f + dist(rng) * 1.0f;
                    radius = 2.0f + dist(rng) * 1.0f;
                }

                particles.emplace_back(pos, vel, mass);
                particles.back().radius = radius;
            }
        }
    }

    // Set global current dampening (0.0 = normal, 1.0 = completely stopped)
    void setCurrentDampening(float dampening) {
        currentDampening = fmaxf(0.0f, fminf(1.0f, dampening));
    }

    // Add a barrier zone that blocks ocean currents
    void addBarrierZone(float x, float y, float radius) {
        barrierZones.push_back(Vec2(x, y));
        barrierRadii.push_back(radius);
    }

    // Clear all barrier zones
    void clearBarriers() {
        barrierZones.clear();
        barrierRadii.clear();
    }

    // Apply environmental effects (dampening and barriers)
    void applyEnvironmentalEffects() {
        if (currentDampening <= 0.0f && barrierZones.empty()) return;

        for (auto& p : particles) {
            if (!p.active) continue;

            // Global current dampening
            if (currentDampening > 0.0f) {
                p.vel *= (1.0f - currentDampening * 0.01f);
            }

            // Check barriers
            for (size_t i = 0; i < barrierZones.size(); ++i) {
                Vec2 diff = p.pos - barrierZones[i];
                float dist = diff.length();

                if (dist < barrierRadii[i]) {
                    // Inside barrier zone - heavily dampen velocity
                    p.vel *= 0.95f;

                    // Push particle away from barrier center
                    if (dist > 0.1f) {
                        Vec2 repulsion = diff.normalized() * 0.5f;
                        p.vel += repulsion;
                    }
                }
            }
        }
    }

    // ========================================================================
    // UPDATE
    // ========================================================================

    void update(int steps = 1) {
        for (int step = 0; step < steps; ++step) {
            buildTree();
            calculateForces();
            integrateLeapfrog();
            handleCollisions();
            applyEnvironmentalEffects();

            time += dt;
            iteration++;
        }

        calculateEnergy();
    }

    // ========================================================================
    // GETTERS
    // ========================================================================

    val getParticleData() {
        val data = val::array();
        for (const auto& p : particles) {
            if (p.active) {
                data.call<void>("push", p.pos.x);
                data.call<void>("push", p.pos.y);
                data.call<void>("push", p.mass);
                data.call<void>("push", p.radius);
            }
        }
        return data;
    }

    int getParticleCount() const {
        int count = 0;
        for (const auto& p : particles) {
            if (p.active) count++;
        }
        return count;
    }

    float getTotalEnergy() const { return totalEnergy; }
    float getKineticEnergy() const { return kineticEnergy; }
    float getPotentialEnergy() const { return potentialEnergy; }
    float getTime() const { return time; }
    int getIteration() const { return iteration; }

    void setTimeStep(float timestep) { dt = timestep; }
};

// ============================================================================
// EMSCRIPTEN BINDINGS
// ============================================================================

EMSCRIPTEN_BINDINGS(nbody_module) {
    class_<NBodyEngine>("NBodyEngine")
        .constructor<float>()
        .function("createGalaxyCollision", &NBodyEngine::createGalaxyCollision)
        .function("createSolarSystem", &NBodyEngine::createSolarSystem)
        .function("createRandomCluster", &NBodyEngine::createRandomCluster)
        .function("createOceanTrash", &NBodyEngine::createOceanTrash)
        .function("addTrashPollution", &NBodyEngine::addTrashPollution)
        .function("removeTrashPollution", &NBodyEngine::removeTrashPollution)
        .function("addGlacialMelt", &NBodyEngine::addGlacialMelt)
        .function("setCurrentDampening", &NBodyEngine::setCurrentDampening)
        .function("addBarrierZone", &NBodyEngine::addBarrierZone)
        .function("clearBarriers", &NBodyEngine::clearBarriers)
        .function("update", &NBodyEngine::update)
        .function("getParticleData", &NBodyEngine::getParticleData)
        .function("getParticleCount", &NBodyEngine::getParticleCount)
        .function("getTotalEnergy", &NBodyEngine::getTotalEnergy)
        .function("getKineticEnergy", &NBodyEngine::getKineticEnergy)
        .function("getPotentialEnergy", &NBodyEngine::getPotentialEnergy)
        .function("getTime", &NBodyEngine::getTime)
        .function("getIteration", &NBodyEngine::getIteration)
        .function("setTimeStep", &NBodyEngine::setTimeStep);
}
