#include <emscripten/bind.h>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace emscripten;

const float PI = 3.14159265359f;
const float EARTH_RADIUS = 6371.0f; // km
const float G = 6.67430e-11f; // gravitational constant
const float EARTH_MASS = 5.972e24f; // kg
const float MU = 398600.4418f; // Earth gravitational parameter km³/s²

enum class ObjectType {
    SATELLITE_ACTIVE,
    SATELLITE_DEAD,
    ROCKET_BODY,
    DEBRIS_LARGE,      // >10cm
    DEBRIS_SMALL,      // 1-10cm
    DEBRIS_MICRO       // <1cm
};

struct SpaceObject {
    float semiMajorAxis;  // km
    float eccentricity;
    float inclination;    // radians
    float longitude;      // radians (RAAN - Right Ascension of Ascending Node)
    float argument;       // radians (argument of periapsis)
    float trueAnomaly;    // radians

    float mass;           // kg
    float crossSection;   // m²
    ObjectType type;
    bool active;

    float x, y, z;        // Position km
    float vx, vy, vz;     // Velocity km/s

    int collisionCount;   // Number of collisions this object caused
};

struct CollisionEvent {
    int object1Index;
    int object2Index;
    float time;
    float energy;         // Joules
    int debrisGenerated;
};

class KesslerSyndrome {
private:
    std::vector<SpaceObject> objects;
    std::vector<CollisionEvent> collisionHistory;

    float time;                    // seconds
    int totalCollisions;
    int debrisCount;
    int activeS atellites;
    int deadSatellites;

    // Debris density per orbital shell (objects per km³)
    float densityLEO;      // 200-2000 km
    float densityMEO;      // 2000-35786 km
    float densityGEO;      // ~35786 km

    // Critical threshold for cascade
    float cascadeThreshold;
    bool cascadeStarted;

    void updateOrbitalPositions(float dt) {
        for (auto& obj : objects) {
            // Simplified orbital propagation
            float r = sqrtf(obj.x * obj.x + obj.y * obj.y + obj.z * obj.z);
            float v = sqrtf(obj.vx * obj.vx + obj.vy * obj.vy + obj.vz * obj.vz);

            // Gravitational acceleration
            float ax = -MU * obj.x / (r * r * r);
            float ay = -MU * obj.y / (r * r * r);
            float az = -MU * obj.z / (r * r * r);

            // Simple Euler integration (for speed)
            obj.vx += ax * dt;
            obj.vy += ay * dt;
            obj.vz += az * dt;

            obj.x += obj.vx * dt;
            obj.y += obj.vy * dt;
            obj.z += obj.vz * dt;
        }
    }

    void checkCollisions() {
        for (size_t i = 0; i < objects.size(); i++) {
            for (size_t j = i + 1; j < objects.size(); j++) {
                float dx = objects[i].x - objects[j].x;
                float dy = objects[i].y - objects[j].y;
                float dz = objects[i].z - objects[j].z;
                float distance = sqrtf(dx * dx + dy * dy + dz * dz);

                // Collision threshold (very simplified)
                float threshold = 0.01f; // 10 meters

                if (distance < threshold) {
                    handleCollision(i, j);
                }
            }
        }
    }

    void handleCollision(int idx1, int idx2) {
        SpaceObject& obj1 = objects[idx1];
        SpaceObject& obj2 = objects[idx2];

        // Calculate relative velocity
        float dvx = obj1.vx - obj2.vx;
        float dvy = obj1.vy - obj2.vy;
        float dvz = obj1.vz - obj2.vz;
        float relVel = sqrtf(dvx * dvx + dvy * dvy + dvz * dvz);

        // Kinetic energy (simplified)
        float energy = 0.5f * obj1.mass * relVel * relVel * 1000.0f; // Joules

        // Generate debris based on NASA standard breakup model
        int debrisGenerated = generateDebris(obj1, obj2, energy);

        totalCollisions++;
        obj1.collisionCount++;
        obj2.collisionCount++;

        // Record collision
        CollisionEvent event;
        event.object1Index = idx1;
        event.object2Index = idx2;
        event.time = time;
        event.energy = energy;
        event.debrisGenerated = debrisGenerated;
        collisionHistory.push_back(event);

        // Check cascade condition
        updateDensity();
        if (densityLEO > cascadeThreshold && !cascadeStarted) {
            cascadeStarted = true;
        }
    }

    int generateDebris(const SpaceObject& obj1, const SpaceObject& obj2, float energy) {
        // NASA standard breakup model (simplified)
        float totalMass = obj1.mass + obj2.mass;
        float specificEnergy = energy / totalMass; // J/kg

        // Number of debris pieces >10cm
        int debrisCount = static_cast<int>(0.1f * powf(specificEnergy / 1000.0f, 0.75f));
        debrisCount = std::min(debrisCount, 1000); // Cap for performance

        // Generate debris objects
        for (int i = 0; i < debrisCount; i++) {
            SpaceObject debris;

            // Position near collision point
            debris.x = (obj1.x + obj2.x) / 2.0f + ((rand() % 100) - 50) * 0.001f;
            debris.y = (obj1.y + obj2.y) / 2.0f + ((rand() % 100) - 50) * 0.001f;
            debris.z = (obj1.z + obj2.z) / 2.0f + ((rand() % 100) - 50) * 0.001f;

            // Velocity with random perturbation
            float dvx = ((rand() % 200) - 100) * 0.001f;
            float dvy = ((rand() % 200) - 100) * 0.001f;
            float dvz = ((rand() % 200) - 100) * 0.001f;

            debris.vx = (obj1.vx + obj2.vx) / 2.0f + dvx;
            debris.vy = (obj1.vy + obj2.vy) / 2.0f + dvy;
            debris.vz = (obj1.vz + obj2.vz) / 2.0f + dvz;

            debris.mass = totalMass / debrisCount;
            debris.crossSection = 0.01f + (rand() % 10) * 0.01f; // 1-10 cm²
            debris.type = ObjectType::DEBRIS_LARGE;
            debris.active = false;
            debris.collisionCount = 0;

            objects.push_back(debris);
        }

        this->debrisCount += debrisCount;
        return debrisCount;
    }

    void updateDensity() {
        int leoCount = 0, meoCount = 0, geoCount = 0;

        for (const auto& obj : objects) {
            float r = sqrtf(obj.x * obj.x + obj.y * obj.y + obj.z * obj.z);
            float altitude = r - EARTH_RADIUS;

            if (altitude < 2000.0f) leoCount++;
            else if (altitude < 35786.0f) meoCount++;
            else geoCount++;
        }

        // Simplified density calculation (objects per km³)
        float leoVolume = (4.0f / 3.0f) * PI * (powf(EARTH_RADIUS + 2000.0f, 3) - powf(EARTH_RADIUS + 200.0f, 3));
        float meoVolume = (4.0f / 3.0f) * PI * (powf(EARTH_RADIUS + 35786.0f, 3) - powf(EARTH_RADIUS + 2000.0f, 3));
        float geoVolume = (4.0f / 3.0f) * PI * powf(100.0f, 3); // Thin shell

        densityLEO = leoCount / leoVolume * 1e9f; // Normalize
        densityMEO = meoCount / meoVolume * 1e9f;
        densityGEO = geoCount / geoVolume * 1e9f;
    }

public:
    KesslerSyndrome() {
        reset();
    }

    void reset() {
        objects.clear();
        collisionHistory.clear();
        time = 0.0f;
        totalCollisions = 0;
        debrisCount = 0;
        activeSatellites = 0;
        deadSatellites = 0;
        densityLEO = 0.0f;
        densityMEO = 0.0f;
        densityGEO = 0.0f;
        cascadeThreshold = 100.0f; // Arbitrary threshold
        cascadeStarted = false;
    }

    void addSatellite(float altitude, float inclination, int satType) {
        SpaceObject sat;

        sat.semiMajorAxis = EARTH_RADIUS + altitude;
        sat.eccentricity = 0.0f;
        sat.inclination = inclination * PI / 180.0f;
        sat.longitude = (rand() % 360) * PI / 180.0f;
        sat.argument = 0.0f;
        sat.trueAnomaly = (rand() % 360) * PI / 180.0f;

        // Convert to Cartesian coordinates
        float r = sat.semiMajorAxis;
        sat.x = r * cosf(sat.trueAnomaly) * cosf(sat.longitude);
        sat.y = r * cosf(sat.trueAnomaly) * sinf(sat.longitude);
        sat.z = r * sinf(sat.trueAnomaly) * sinf(sat.inclination);

        // Circular orbit velocity
        float v = sqrtf(MU / r);
        sat.vx = -v * sinf(sat.trueAnomaly) * cosf(sat.longitude);
        sat.vy = -v * sinf(sat.trueAnomaly) * sinf(sat.longitude);
        sat.vz = v * cosf(sat.trueAnomaly) * sinf(sat.inclination);

        sat.mass = 1000.0f + (rand() % 5000); // 1-6 tons
        sat.crossSection = 10.0f; // 10 m²
        sat.type = static_cast<ObjectType>(satType);
        sat.active = (satType == 0);
        sat.collisionCount = 0;

        objects.push_back(sat);

        if (sat.active) activeSatellites++;
        else if (satType == 1) deadSatellites++;
    }

    void addDebrisCloud(float altitude, int count) {
        for (int i = 0; i < count; i++) {
            addSatellite(altitude, rand() % 180, 3); // Debris type
        }
    }

    void triggerCollision() {
        // Manually trigger a collision in LEO
        if (objects.size() >= 2) {
            handleCollision(0, 1);
        }
    }

    // Preset scenarios
    void setCurrentState() {
        reset();
        // Approximate current satellite population (2024)
        // LEO
        for (int i = 0; i < 3000; i++) {
            addSatellite(400.0f + (rand() % 1600), rand() % 180, rand() % 2);
        }
        // MEO
        for (int i = 0; i < 100; i++) {
            addSatellite(20000.0f, 55.0f, rand() % 2);
        }
        // GEO
        for (int i = 0; i < 500; i++) {
            addSatellite(35786.0f, 0.0f, rand() % 2);
        }
        updateDensity();
    }

    void setFuture2050() {
        reset();
        // Projected 2050 with mega-constellations
        for (int i = 0; i < 50000; i++) {
            addSatellite(400.0f + (rand() % 1000), rand() % 180, rand() % 3);
        }
        updateDensity();
    }

    void setCascadeScenario() {
        reset();
        // Dense LEO population primed for cascade
        for (int i = 0; i < 5000; i++) {
            addSatellite(800.0f, rand() % 180, rand() % 4);
        }
        // Trigger initial collision
        triggerCollision();
    }

    void update(float dt) {
        time += dt;

        updateOrbitalPositions(dt);

        // Check collisions (expensive, limit frequency)
        if (static_cast<int>(time * 10) % 5 == 0) {
            checkCollisions();
        }

        updateDensity();
    }

    // Getters
    int getObjectCount() const { return objects.size(); }
    int getActiveSatellites() const {
        int count = 0;
        for (const auto& obj : objects) {
            if (obj.type == ObjectType::SATELLITE_ACTIVE) count++;
        }
        return count;
    }
    int getDebrisCount() const { return debrisCount; }
    int getTotalCollisions() const { return totalCollisions; }
    float getTime() const { return time; }

    float getDensityLEO() const { return densityLEO; }
    float getDensityMEO() const { return densityMEO; }
    float getDensityGEO() const { return densityGEO; }

    bool isCascading() const { return cascadeStarted; }

    std::vector<float> getObjectPositions() const {
        std::vector<float> positions;
        for (const auto& obj : objects) {
            positions.push_back(obj.x);
            positions.push_back(obj.y);
            positions.push_back(obj.z);
            positions.push_back(static_cast<float>(obj.type));
        }
        return positions;
    }

    int getCollisionCount() const { return collisionHistory.size(); }

    std::vector<float> getRecentCollisions() const {
        std::vector<float> collisions;
        int count = std::min(10, static_cast<int>(collisionHistory.size()));
        for (int i = collisionHistory.size() - count; i < collisionHistory.size(); i++) {
            collisions.push_back(collisionHistory[i].time);
            collisions.push_back(collisionHistory[i].energy);
            collisions.push_back(collisionHistory[i].debrisGenerated);
        }
        return collisions;
    }
};

EMSCRIPTEN_BINDINGS(kessler_module) {
    class_<KesslerSyndrome>("KesslerSyndrome")
        .constructor<>()
        .function("reset", &KesslerSyndrome::reset)
        .function("update", &KesslerSyndrome::update)
        .function("addSatellite", &KesslerSyndrome::addSatellite)
        .function("addDebrisCloud", &KesslerSyndrome::addDebrisCloud)
        .function("triggerCollision", &KesslerSyndrome::triggerCollision)
        .function("setCurrentState", &KesslerSyndrome::setCurrentState)
        .function("setFuture2050", &KesslerSyndrome::setFuture2050)
        .function("setCascadeScenario", &KesslerSyndrome::setCascadeScenario)
        .function("getObjectCount", &KesslerSyndrome::getObjectCount)
        .function("getActiveSatellites", &KesslerSyndrome::getActiveSatellites)
        .function("getDebrisCount", &KesslerSyndrome::getDebrisCount)
        .function("getTotalCollisions", &KesslerSyndrome::getTotalCollisions)
        .function("getTime", &KesslerSyndrome::getTime)
        .function("getDensityLEO", &KesslerSyndrome::getDensityLEO)
        .function("getDensityMEO", &KesslerSyndrome::getDensityMEO)
        .function("getDensityGEO", &KesslerSyndrome::getDensityGEO)
        .function("isCascading", &KesslerSyndrome::isCascading)
        .function("getObjectPositions", &KesslerSyndrome::getObjectPositions)
        .function("getCollisionCount", &KesslerSyndrome::getCollisionCount)
        .function("getRecentCollisions", &KesslerSyndrome::getRecentCollisions);

    register_vector<float>("VectorFloat");
}
