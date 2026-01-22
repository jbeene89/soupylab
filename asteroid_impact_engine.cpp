#include <emscripten/bind.h>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace emscripten;

const float PI = 3.14159265359f;
const float EARTH_RADIUS = 6371000.0f; // meters
const float G = 6.67430e-11f; // gravitational constant
const float EARTH_MASS = 5.972e24f; // kg
const float ATMOSPHERE_SCALE_HEIGHT = 8500.0f; // meters
const float SEA_LEVEL_DENSITY = 1.225f; // kg/m³
const float SPEED_OF_SOUND = 340.0f; // m/s

enum class CompositionType {
    ROCKY,      // Density ~3000 kg/m³
    IRON,       // Density ~7800 kg/m³
    ICY,        // Density ~1000 kg/m³
    CARBONACEOUS // Density ~2000 kg/m³
};

struct ImpactStats {
    float energy; // Joules
    float megatons; // TNT equivalent
    float craterDiameter; // meters
    float craterDepth; // meters
    float thermalRadius; // meters (severe burns)
    float blastRadius; // meters (building destruction)
    float seismicMagnitude; // Richter scale
    float ejectaRadius; // meters
    float tsunamiHeight; // meters (if ocean impact)
    bool extinctionLevel;
    float globalDustMass; // kg
    float temperatureDrop; // °C
};

class AsteroidImpact {
private:
    // Asteroid properties
    float diameter; // meters
    float density; // kg/m³
    float velocity; // m/s
    float entryAngle; // degrees (0 = horizontal, 90 = vertical)
    CompositionType composition;

    // Position and state
    float altitude; // meters
    float x, y; // position on Earth surface (lat/lon in radians)
    float velocityCurrent; // current velocity during descent
    bool hasImpacted;
    bool oceanImpact;

    // Simulation state
    float time;
    ImpactStats stats;
    std::vector<float> trail; // altitude history

    // Impact effects
    float shockwaveFront; // expanding shockwave distance
    float thermalPulseIntensity;
    float seismicWaves;

    float getMass() const {
        float radius = diameter / 2.0f;
        return (4.0f / 3.0f) * PI * radius * radius * radius * density;
    }

    float getAtmosphericDensity(float alt) const {
        return SEA_LEVEL_DENSITY * expf(-alt / ATMOSPHERE_SCALE_HEIGHT);
    }

    float getDragCoefficient() const {
        // Simplified drag coefficient
        return 0.47f; // sphere approximation
    }

    void calculateImpactEffects() {
        float mass = getMass();
        float impactVelocity = velocityCurrent;

        // Kinetic energy (E = 0.5 * m * v²)
        stats.energy = 0.5f * mass * impactVelocity * impactVelocity;

        // Convert to megatons TNT (1 megaton = 4.184e15 J)
        stats.megatons = stats.energy / 4.184e15f;

        // Crater formation (using scaling laws)
        // Simple crater: D = 1.161 * (E^0.22) * (ρ_t^-0.33) * g^-0.22
        float targetDensity = oceanImpact ? 1000.0f : 2500.0f; // water or rock
        float g = 9.81f;

        // Crater diameter (meters) - empirical scaling
        stats.craterDiameter = 1.161f * powf(stats.energy, 0.22f) *
                               powf(targetDensity, -0.33f) * powf(g, -0.22f);

        // Crater depth (typically 1/5 to 1/10 of diameter)
        stats.craterDepth = stats.craterDiameter / 7.0f;

        // Thermal radiation radius (3rd degree burns)
        // Q = σ * T^4 for ~3000K fireball
        stats.thermalRadius = 2.5f * powf(stats.megatons, 0.41f) * 1000.0f; // meters

        // Blast radius (overpressure > 20 psi, building destruction)
        stats.blastRadius = 3.5f * powf(stats.megatons, 0.33f) * 1000.0f; // meters

        // Seismic magnitude (Richter scale)
        // M = 0.67 * log10(E) - 5.87
        stats.seismicMagnitude = 0.67f * log10f(stats.energy) - 5.87f;

        // Ejecta radius (material thrown out)
        stats.ejectaRadius = stats.craterDiameter * 2.5f;

        // Tsunami height (for ocean impacts)
        if (oceanImpact) {
            // Simplified tsunami scaling: h ∝ E^0.25
            stats.tsunamiHeight = 0.1f * powf(stats.megatons, 0.25f) * 100.0f; // meters
        } else {
            stats.tsunamiHeight = 0.0f;
        }

        // Global effects
        // Dust mass ejected (fraction of crater volume)
        float craterVolume = PI * powf(stats.craterDiameter / 2.0f, 2.0f) * stats.craterDepth;
        stats.globalDustMass = craterVolume * targetDensity * 0.5f; // 50% becomes dust

        // Temperature drop (simplified climate model)
        // Larger impacts inject more dust into stratosphere
        if (stats.megatons > 100000.0f) {
            stats.temperatureDrop = 10.0f + (stats.megatons / 1000000.0f); // Severe nuclear winter
            stats.extinctionLevel = true;
        } else if (stats.megatons > 10000.0f) {
            stats.temperatureDrop = 5.0f + (stats.megatons / 100000.0f);
            stats.extinctionLevel = true;
        } else if (stats.megatons > 1000.0f) {
            stats.temperatureDrop = 2.0f + (stats.megatons / 10000.0f);
            stats.extinctionLevel = false;
        } else {
            stats.temperatureDrop = 0.5f * (stats.megatons / 1000.0f);
            stats.extinctionLevel = false;
        }
    }

public:
    AsteroidImpact() {
        // Default: Tunguska-sized event
        diameter = 60.0f; // 60 meters
        composition = CompositionType::ROCKY;
        density = 3000.0f;
        velocity = 30000.0f; // 30 km/s
        entryAngle = 45.0f;

        reset();
    }

    void reset() {
        altitude = 100000.0f; // Start at 100 km
        velocityCurrent = velocity;
        x = 0.0f;
        y = 0.0f;
        hasImpacted = false;
        oceanImpact = false; // Default to land
        time = 0.0f;
        shockwaveFront = 0.0f;
        thermalPulseIntensity = 0.0f;
        seismicWaves = 0.0f;
        trail.clear();

        stats = ImpactStats{};
    }

    void setDiameter(float d) {
        diameter = std::max(1.0f, std::min(d, 15000.0f)); // 1m to 15km
        reset();
    }

    void setVelocity(float v) {
        velocity = std::max(5000.0f, std::min(v, 72000.0f)); // 5-72 km/s
        velocityCurrent = velocity;
    }

    void setEntryAngle(float angle) {
        entryAngle = std::max(0.0f, std::min(angle, 90.0f));
    }

    void setComposition(int type) {
        composition = static_cast<CompositionType>(type);
        switch(composition) {
            case CompositionType::ROCKY:
                density = 3000.0f;
                break;
            case CompositionType::IRON:
                density = 7800.0f;
                break;
            case CompositionType::ICY:
                density = 1000.0f;
                break;
            case CompositionType::CARBONACEOUS:
                density = 2000.0f;
                break;
        }
        reset();
    }

    void setOceanImpact(bool ocean) {
        oceanImpact = ocean;
    }

    // Preset scenarios
    void setTunguska() {
        diameter = 60.0f;
        velocity = 30000.0f;
        entryAngle = 30.0f;
        composition = CompositionType::ROCKY;
        density = 3000.0f;
        oceanImpact = false;
        reset();
    }

    void setBarringerCrater() {
        diameter = 50.0f;
        velocity = 12800.0f;
        entryAngle = 80.0f;
        composition = CompositionType::IRON;
        density = 7800.0f;
        oceanImpact = false;
        reset();
    }

    void setChicxulub() {
        diameter = 10000.0f; // 10 km
        velocity = 20000.0f;
        entryAngle = 60.0f;
        composition = CompositionType::CARBONACEOUS;
        density = 2000.0f;
        oceanImpact = false;
        reset();
    }

    void setChelyabinsk() {
        diameter = 20.0f;
        velocity = 19000.0f;
        entryAngle = 18.0f;
        composition = CompositionType::ROCKY;
        density = 3000.0f;
        oceanImpact = false;
        reset();
    }

    void setApophis() {
        diameter = 370.0f;
        velocity = 12600.0f;
        entryAngle = 45.0f;
        composition = CompositionType::ROCKY;
        density = 3000.0f;
        oceanImpact = false;
        reset();
    }

    void setPlanetKiller() {
        diameter = 15000.0f; // 15 km
        velocity = 72000.0f; // Maximum impact velocity
        entryAngle = 90.0f;
        composition = CompositionType::IRON;
        density = 7800.0f;
        oceanImpact = false;
        reset();
    }

    void update(float dt) {
        if (hasImpacted) {
            // Expand shockwave
            shockwaveFront += SPEED_OF_SOUND * dt * 100.0f; // Faster visualization

            // Decay thermal pulse
            thermalPulseIntensity = std::max(0.0f, thermalPulseIntensity - dt * 0.5f);

            // Seismic wave propagation
            seismicWaves += 5000.0f * dt; // 5 km/s seismic velocity

            time += dt;
            return;
        }

        if (altitude <= 0.0f) {
            // Impact!
            hasImpacted = true;
            calculateImpactEffects();
            thermalPulseIntensity = 1.0f;
            return;
        }

        // Atmospheric descent
        float rho = getAtmosphericDensity(altitude);
        float crossSection = PI * (diameter / 2.0f) * (diameter / 2.0f);
        float dragForce = 0.5f * rho * velocityCurrent * velocityCurrent * crossSection * getDragCoefficient();
        float mass = getMass();

        // Deceleration due to drag
        float deceleration = dragForce / mass;
        velocityCurrent -= deceleration * dt;
        velocityCurrent = std::max(velocityCurrent, 1000.0f); // Minimum velocity

        // Vertical descent (simplified)
        float verticalVelocity = velocityCurrent * sinf(entryAngle * PI / 180.0f);
        altitude -= verticalVelocity * dt;

        // Atmospheric breakup for small/weak objects
        float dynamicPressure = 0.5f * rho * velocityCurrent * velocityCurrent;
        if (dynamicPressure > 1e6f && composition == CompositionType::ICY && diameter < 100.0f) {
            // Object breaks up (simplified)
            velocityCurrent *= 0.95f; // Rapid deceleration
        }

        // Record trail
        if (static_cast<int>(time * 10) % 5 == 0 && trail.size() < 1000) {
            trail.push_back(altitude);
        }

        time += dt;
    }

    // Getters
    float getAltitude() const { return altitude; }
    float getVelocity() const { return velocityCurrent; }
    float getDiameter() const { return diameter; }
    float getMach() const { return velocityCurrent / SPEED_OF_SOUND; }
    bool impacted() const { return hasImpacted; }

    float getImpactEnergy() const { return stats.energy; }
    float getMegatons() const { return stats.megatons; }
    float getCraterDiameter() const { return stats.craterDiameter; }
    float getCraterDepth() const { return stats.craterDepth; }
    float getThermalRadius() const { return stats.thermalRadius; }
    float getBlastRadius() const { return stats.blastRadius; }
    float getSeismicMagnitude() const { return stats.seismicMagnitude; }
    float getEjectaRadius() const { return stats.ejectaRadius; }
    float getTsunamiHeight() const { return stats.tsunamiHeight; }
    bool isExtinctionLevel() const { return stats.extinctionLevel; }
    float getGlobalDustMass() const { return stats.globalDustMass; }
    float getTemperatureDrop() const { return stats.temperatureDrop; }

    float getShockwaveFront() const { return shockwaveFront; }
    float getThermalPulse() const { return thermalPulseIntensity; }
    float getSeismicWaves() const { return seismicWaves; }

    std::vector<float> getTrail() const { return trail; }

    float getTimeToImpact() const {
        if (hasImpacted) return 0.0f;
        float verticalVelocity = velocityCurrent * sinf(entryAngle * PI / 180.0f);
        return altitude / verticalVelocity;
    }
};

EMSCRIPTEN_BINDINGS(asteroid_module) {
    class_<AsteroidImpact>("AsteroidImpact")
        .constructor<>()
        .function("reset", &AsteroidImpact::reset)
        .function("update", &AsteroidImpact::update)
        .function("setDiameter", &AsteroidImpact::setDiameter)
        .function("setVelocity", &AsteroidImpact::setVelocity)
        .function("setEntryAngle", &AsteroidImpact::setEntryAngle)
        .function("setComposition", &AsteroidImpact::setComposition)
        .function("setOceanImpact", &AsteroidImpact::setOceanImpact)
        .function("setTunguska", &AsteroidImpact::setTunguska)
        .function("setBarringerCrater", &AsteroidImpact::setBarringerCrater)
        .function("setChicxulub", &AsteroidImpact::setChicxulub)
        .function("setChelyabinsk", &AsteroidImpact::setChelyabinsk)
        .function("setApophis", &AsteroidImpact::setApophis)
        .function("setPlanetKiller", &AsteroidImpact::setPlanetKiller)
        .function("getAltitude", &AsteroidImpact::getAltitude)
        .function("getVelocity", &AsteroidImpact::getVelocity)
        .function("getDiameter", &AsteroidImpact::getDiameter)
        .function("getMach", &AsteroidImpact::getMach)
        .function("impacted", &AsteroidImpact::impacted)
        .function("getImpactEnergy", &AsteroidImpact::getImpactEnergy)
        .function("getMegatons", &AsteroidImpact::getMegatons)
        .function("getCraterDiameter", &AsteroidImpact::getCraterDiameter)
        .function("getCraterDepth", &AsteroidImpact::getCraterDepth)
        .function("getThermalRadius", &AsteroidImpact::getThermalRadius)
        .function("getBlastRadius", &AsteroidImpact::getBlastRadius)
        .function("getSeismicMagnitude", &AsteroidImpact::getSeismicMagnitude)
        .function("getEjectaRadius", &AsteroidImpact::getEjectaRadius)
        .function("getTsunamiHeight", &AsteroidImpact::getTsunamiHeight)
        .function("isExtinctionLevel", &AsteroidImpact::isExtinctionLevel)
        .function("getGlobalDustMass", &AsteroidImpact::getGlobalDustMass)
        .function("getTemperatureDrop", &AsteroidImpact::getTemperatureDrop)
        .function("getShockwaveFront", &AsteroidImpact::getShockwaveFront)
        .function("getThermalPulse", &AsteroidImpact::getThermalPulse)
        .function("getSeismicWaves", &AsteroidImpact::getSeismicWaves)
        .function("getTrail", &AsteroidImpact::getTrail)
        .function("getTimeToImpact", &AsteroidImpact::getTimeToImpact);

    register_vector<float>("VectorFloat");
}
