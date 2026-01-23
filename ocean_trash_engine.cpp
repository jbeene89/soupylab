/**
 * OCEAN TRASH SIMULATOR - Fluid Dynamics & Pollution Modeling
 * Models plastic pollution in ocean currents with environmental controls
 *
 * Features:
 * - Ocean gyre vortex simulation
 * - Different trash types (bottles, bags, microplastics, debris)
 * - Dynamic pollution addition/removal
 * - Glacial melt effects on currents
 * - Current dampening and barriers
 * - Fluid dynamics (no gravity)
 */

#include <emscripten/bind.h>
#include <emscripten/val.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>

using namespace emscripten;

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
    Vec2 normalized() const {
        float len = length();
        return len > 0 ? *this / len : Vec2(0, 0);
    }
};

// ============================================================================
// TRASH PARTICLE
// ============================================================================

struct TrashParticle {
    Vec2 pos;
    Vec2 vel;
    float mass;
    float radius;
    int type;  // 0=bottle, 1=bag, 2=microplastic, 3=debris
    bool active;

    TrashParticle() : mass(1.0f), radius(1.0f), type(0), active(true) {}

    TrashParticle(Vec2 pos, Vec2 vel, float mass, float radius, int type)
        : pos(pos), vel(vel), mass(mass), radius(radius), type(type), active(true) {}
};

// ============================================================================
// OCEAN VORTEX (GYRE)
// ============================================================================

struct OceanVortex {
    Vec2 center;
    float radius;
    float strength;     // Rotation speed
    float temperature;  // Water temperature (affects current strength)
    bool active;

    OceanVortex(Vec2 c, float r, float s, float temp = 20.0f)
        : center(c), radius(r), strength(s), temperature(temp), active(true) {}

    // Calculate current velocity at a position
    Vec2 getCurrentAt(const Vec2& pos) const {
        if (!active) return Vec2(0, 0);

        Vec2 diff = pos - center;
        float dist = diff.length();

        if (dist < 0.1f || dist > radius) return Vec2(0, 0);

        // Circular current (tangential velocity)
        float angle = atan2f(diff.y, diff.x);
        float falloff = 1.0f - (dist / radius);  // Stronger near center
        float speed = strength * falloff;

        return Vec2(-sinf(angle) * speed, cosf(angle) * speed);
    }
};

// ============================================================================
// ICE/FRESHWATER MASS
// ============================================================================

struct IceMass {
    Vec2 pos;
    Vec2 vel;
    float mass;
    float radius;
    float meltRate;  // How fast it's melting
    bool active;

    IceMass(Vec2 p, Vec2 v, float m, float r)
        : pos(p), vel(v), mass(m), radius(r), meltRate(0.01f), active(true) {}
};

// ============================================================================
// BARRIER ZONE
// ============================================================================

struct Barrier {
    Vec2 center;
    float radius;
    float dampeningFactor;  // How much it slows currents

    Barrier(Vec2 c, float r, float d = 0.9f)
        : center(c), radius(r), dampeningFactor(d) {}
};

// ============================================================================
// OCEAN TRASH ENGINE
// ============================================================================

class OceanTrashEngine {
private:
    std::vector<TrashParticle> particles;
    std::vector<OceanVortex> vortices;
    std::vector<IceMass> iceMasses;
    std::vector<Barrier> barriers;

    float worldSize;
    float time;
    float dt;
    int iteration;

    // Environmental parameters
    float currentDampening;   // Global current reduction (0-1)
    float globalTemp;         // Ocean temperature

    std::mt19937 rng;

    // Apply ocean currents to particles
    void applyOceanCurrents() {
        for (auto& p : particles) {
            if (!p.active) continue;

            Vec2 currentForce(0, 0);

            // Sum forces from all vortices
            for (const auto& vortex : vortices) {
                Vec2 vortexCurrent = vortex.getCurrentAt(p.pos);
                currentForce += vortexCurrent;
            }

            // Apply force (lighter particles affected more)
            float responsiveness = 1.0f / (p.mass + 0.1f);
            p.vel += currentForce * responsiveness * dt;

            // Drag (water resistance)
            float drag = 0.98f;
            p.vel *= drag;
        }
    }

    // Apply ice/freshwater effects
    void applyIceEffects() {
        for (auto& ice : iceMasses) {
            if (!ice.active) continue;

            // Ice disrupts nearby currents
            for (auto& vortex : vortices) {
                Vec2 diff = ice.pos - vortex.center;
                float dist = diff.length();

                if (dist < vortex.radius * 0.5f) {
                    // Ice weakens vortex
                    vortex.strength *= 0.99f;
                }
            }

            // Ice slowly melts
            ice.mass -= ice.meltRate;
            if (ice.mass <= 0.1f) {
                ice.active = false;
            }

            // Ice drifts slowly
            ice.vel *= 0.95f;
            ice.pos += ice.vel * dt;
        }
    }

    // Apply barriers
    void applyBarriers() {
        if (barriers.empty()) return;

        for (auto& p : particles) {
            if (!p.active) continue;

            for (const auto& barrier : barriers) {
                Vec2 diff = p.pos - barrier.center;
                float dist = diff.length();

                if (dist < barrier.radius) {
                    // Dampen velocity inside barrier
                    p.vel *= barrier.dampeningFactor;

                    // Push away from center
                    if (dist > 0.1f) {
                        Vec2 repulsion = diff.normalized() * 2.0f;
                        p.vel += repulsion;
                    }
                }
            }
        }
    }

    // Integrate particle positions
    void integrateParticles() {
        for (auto& p : particles) {
            if (!p.active) continue;

            // Update position
            p.pos += p.vel * dt;

            // Wrap around boundaries
            if (p.pos.x < -worldSize) p.pos.x += 2 * worldSize;
            if (p.pos.x > worldSize) p.pos.x -= 2 * worldSize;
            if (p.pos.y < -worldSize) p.pos.y += 2 * worldSize;
            if (p.pos.y > worldSize) p.pos.y -= 2 * worldSize;

            // Apply global dampening
            if (currentDampening > 0.0f) {
                p.vel *= (1.0f - currentDampening * 0.01f);
            }
        }
    }

public:
    OceanTrashEngine(float worldSize)
        : worldSize(worldSize)
        , time(0)
        , dt(0.5f)
        , iteration(0)
        , currentDampening(0.0f)
        , globalTemp(20.0f)
        , rng(12345)
    {}

    // ========================================================================
    // SCENARIO CREATION
    // ========================================================================

    void createOceanTrash(int count) {
        particles.clear();
        vortices.clear();

        std::uniform_real_distribution<float> dist(0, 1);

        // Create ocean vortices (major garbage patches)
        int numVortices = 5;
        for (int v = 0; v < numVortices; ++v) {
            float angle = (v / (float)numVortices) * 2 * 3.14159f;
            float vortexRadius = worldSize * 0.4f;
            Vec2 center(
                cosf(angle) * vortexRadius,
                sinf(angle) * vortexRadius
            );

            float radius = worldSize * 0.25f;
            float strength = 15.0f + dist(rng) * 10.0f;
            float temp = 18.0f + dist(rng) * 10.0f;

            vortices.emplace_back(center, radius, strength, temp);
        }

        // Add trash to vortices
        int trashPerVortex = count / numVortices;

        for (const auto& vortex : vortices) {
            for (int i = 0; i < trashPerVortex; ++i) {
                // Distribute in spiral around vortex
                float r = sqrtf(dist(rng)) * vortex.radius * 0.6f;
                float theta = dist(rng) * 2 * 3.14159f;

                Vec2 pos = vortex.center + Vec2(cosf(theta) * r, sinf(theta) * r);

                // Initial velocity from current
                Vec2 vel = vortex.getCurrentAt(pos);
                vel += Vec2(
                    (dist(rng) - 0.5f) * 5.0f,
                    (dist(rng) - 0.5f) * 5.0f
                );

                // Random trash type
                float typeRoll = dist(rng);
                int type;
                float mass, radius;

                if (typeRoll < 0.4f) {
                    // Plastic bottles
                    type = 0;
                    mass = 0.1f + dist(rng) * 0.2f;
                    radius = 1.5f + dist(rng) * 0.5f;
                } else if (typeRoll < 0.7f) {
                    // Plastic bags
                    type = 1;
                    mass = 0.05f + dist(rng) * 0.1f;
                    radius = 0.8f + dist(rng) * 0.3f;
                } else if (typeRoll < 0.9f) {
                    // Microplastics
                    type = 2;
                    mass = 0.02f + dist(rng) * 0.05f;
                    radius = 0.3f + dist(rng) * 0.2f;
                } else {
                    // Larger debris
                    type = 3;
                    mass = 0.3f + dist(rng) * 0.4f;
                    radius = 2.0f + dist(rng) * 1.0f;
                }

                particles.emplace_back(pos, vel, mass, radius, type);
            }
        }

        // Add scattered trash
        int scatteredCount = count - (trashPerVortex * numVortices);
        std::uniform_real_distribution<float> scatterPos(-worldSize * 0.6f, worldSize * 0.6f);
        std::uniform_real_distribution<float> scatterVel(-8.0f, 8.0f);

        for (int i = 0; i < scatteredCount; ++i) {
            Vec2 pos(scatterPos(rng), scatterPos(rng));
            Vec2 vel(scatterVel(rng), scatterVel(rng));

            float mass = 0.05f + dist(rng) * 0.25f;
            float radius = 0.5f + dist(rng) * 1.5f;
            int type = static_cast<int>(dist(rng) * 4);

            particles.emplace_back(pos, vel, mass, radius, type);
        }
    }

    // ========================================================================
    // ENVIRONMENTAL CONTROLS
    // ========================================================================

    void addTrashPollution(int amount) {
        std::uniform_real_distribution<float> dist(0, 1);
        std::uniform_real_distribution<float> posDist(-worldSize * 0.6f, worldSize * 0.6f);
        std::uniform_real_distribution<float> velDist(-8.0f, 8.0f);

        for (int i = 0; i < amount; ++i) {
            Vec2 pos(posDist(rng), posDist(rng));
            Vec2 vel(velDist(rng), velDist(rng));

            float typeRoll = dist(rng);
            int type;
            float mass, radius;

            if (typeRoll < 0.4f) {
                type = 0;
                mass = 0.1f + dist(rng) * 0.2f;
                radius = 1.5f + dist(rng) * 0.5f;
            } else if (typeRoll < 0.7f) {
                type = 1;
                mass = 0.05f + dist(rng) * 0.1f;
                radius = 0.8f + dist(rng) * 0.3f;
            } else if (typeRoll < 0.9f) {
                type = 2;
                mass = 0.02f + dist(rng) * 0.05f;
                radius = 0.3f + dist(rng) * 0.2f;
            } else {
                type = 3;
                mass = 0.3f + dist(rng) * 0.4f;
                radius = 2.0f + dist(rng) * 1.0f;
            }

            particles.emplace_back(pos, vel, mass, radius, type);
        }
    }

    void removeTrashPollution(int amount) {
        int removed = 0;
        for (size_t i = particles.size(); i > 0 && removed < amount; --i) {
            if (particles[i - 1].active && particles[i - 1].mass < 0.5f) {
                particles[i - 1].active = false;
                removed++;
            }
        }
    }

    void addGlacialMelt(int intensity) {
        std::uniform_real_distribution<float> dist(0, 1);

        int icebergsPerPole = intensity / 2;

        for (int pole = 0; pole < 2; ++pole) {
            float polarY = (pole == 0) ? -worldSize * 0.7f : worldSize * 0.7f;

            for (int i = 0; i < icebergsPerPole; ++i) {
                float x = (dist(rng) - 0.5f) * worldSize * 1.2f;
                float y = polarY + (dist(rng) - 0.5f) * worldSize * 0.2f;
                Vec2 pos(x, y);

                float driftSpeed = 5.0f + dist(rng) * 5.0f;
                float driftAngle = dist(rng) * 2 * 3.14159f;
                Vec2 vel(
                    cosf(driftAngle) * driftSpeed,
                    sinf(driftAngle) * driftSpeed
                );

                float mass = 2.0f + dist(rng) * 3.0f;
                float radius = 50.0f + dist(rng) * 30.0f;

                iceMasses.emplace_back(pos, vel, mass, radius);
            }
        }
    }

    void setCurrentDampening(float dampening) {
        currentDampening = fmaxf(0.0f, fminf(1.0f, dampening));
    }

    void addBarrierZone(float x, float y, float radius) {
        barriers.emplace_back(Vec2(x, y), radius);
    }

    void clearBarriers() {
        barriers.clear();
    }

    void setGlobalTemp(float temp) {
        globalTemp = temp;

        // Adjust vortex strength based on temperature
        for (auto& vortex : vortices) {
            float tempDiff = fabs(temp - vortex.temperature);
            if (tempDiff > 10.0f) {
                vortex.strength *= 0.9f;  // Different temps disrupt currents
            }
        }
    }

    // ========================================================================
    // UPDATE
    // ========================================================================

    void update(int steps = 1) {
        for (int step = 0; step < steps; ++step) {
            applyOceanCurrents();
            applyIceEffects();
            applyBarriers();
            integrateParticles();

            time += dt;
            iteration++;
        }
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
                data.call<void>("push", p.type);
            }
        }
        return data;
    }

    val getIceData() {
        val data = val::array();
        for (const auto& ice : iceMasses) {
            if (ice.active) {
                data.call<void>("push", ice.pos.x);
                data.call<void>("push", ice.pos.y);
                data.call<void>("push", ice.radius);
            }
        }
        return data;
    }

    val getVortexData() {
        val data = val::array();
        for (const auto& v : vortices) {
            if (v.active) {
                data.call<void>("push", v.center.x);
                data.call<void>("push", v.center.y);
                data.call<void>("push", v.radius);
                data.call<void>("push", v.strength);
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

    int getIceCount() const {
        int count = 0;
        for (const auto& ice : iceMasses) {
            if (ice.active) count++;
        }
        return count;
    }

    float getTime() const { return time; }
    int getIteration() const { return iteration; }

    void setTimeStep(float timestep) { dt = timestep; }
};

// ============================================================================
// EMSCRIPTEN BINDINGS
// ============================================================================

EMSCRIPTEN_BINDINGS(ocean_module) {
    class_<OceanTrashEngine>("OceanTrashEngine")
        .constructor<float>()
        .function("createOceanTrash", &OceanTrashEngine::createOceanTrash)
        .function("addTrashPollution", &OceanTrashEngine::addTrashPollution)
        .function("removeTrashPollution", &OceanTrashEngine::removeTrashPollution)
        .function("addGlacialMelt", &OceanTrashEngine::addGlacialMelt)
        .function("setCurrentDampening", &OceanTrashEngine::setCurrentDampening)
        .function("addBarrierZone", &OceanTrashEngine::addBarrierZone)
        .function("clearBarriers", &OceanTrashEngine::clearBarriers)
        .function("setGlobalTemp", &OceanTrashEngine::setGlobalTemp)
        .function("update", &OceanTrashEngine::update)
        .function("getParticleData", &OceanTrashEngine::getParticleData)
        .function("getIceData", &OceanTrashEngine::getIceData)
        .function("getVortexData", &OceanTrashEngine::getVortexData)
        .function("getParticleCount", &OceanTrashEngine::getParticleCount)
        .function("getIceCount", &OceanTrashEngine::getIceCount)
        .function("getTime", &OceanTrashEngine::getTime)
        .function("getIteration", &OceanTrashEngine::getIteration)
        .function("setTimeStep", &OceanTrashEngine::setTimeStep);
}
