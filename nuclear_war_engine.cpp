#include <emscripten/bind.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>

using namespace emscripten;

const float PI = 3.14159265359f;
const float EARTH_RADIUS = 6371000.0f; // meters
const float G = 9.81f; // gravity

enum class WarheadType {
    TACTICAL_10KT,      // 10 kilotons (battlefield)
    STRATEGIC_100KT,    // 100 kilotons (city target)
    STRATEGIC_500KT,    // 500 kilotons (hardened target)
    MIRV_350KT,         // 350 kilotons (multiple warheads)
    TSAR_BOMBA_50MT     // 50 megatons (largest ever detonated)
};

struct Warhead {
    float yield; // kilotons
    float lat, lon; // target coordinates
    float altitude; // detonation altitude (0 = ground burst)
    bool detonated;
    float time; // time since launch
    float flightTime; // total flight time
    WarheadType type;

    // Blast effects
    float blastRadius_5psi;  // Building damage
    float blastRadius_20psi; // Complete destruction
    float thermalRadius; // 3rd degree burns
    float radiationRadius; // Lethal radiation
    float falloutRadius; // Dangerous fallout (ground burst only)
    float fireballRadius;
};

struct NuclearWinterStats {
    float smokeMass; // Megatons of smoke into stratosphere
    float temperatureDrop; // Global average °C drop
    float ozoneLoss; // Percentage ozone layer destroyed
    int yearsOfDarkness; // Years of reduced sunlight
    float agricultureCollapse; // Percentage of crops lost
};

class NuclearWar {
private:
    std::vector<Warhead> warheads;
    float time;
    int detonatedCount;
    float totalYield;

    // Global effects
    NuclearWinterStats winterStats;
    float globalFallout; // Total radioactive material
    float empEffect; // EMP coverage

    // Casualties (billions)
    float immediateCasualties;
    float radiationCasualties;
    float starvationCasualties;

    void calculateBlastEffects(Warhead& w) {
        float yieldMT = w.yield / 1000.0f; // Convert to megatons

        // Fireball radius (meters)
        w.fireballRadius = 150.0f * powf(yieldMT, 0.4f);

        // Overpressure radii (scaling laws from nuclear weapons effects)
        // 20 psi: Complete building destruction
        w.blastRadius_20psi = 1600.0f * powf(yieldMT, 0.33f);

        // 5 psi: Moderate building damage
        w.blastRadius_5psi = 3500.0f * powf(yieldMT, 0.33f);

        // Thermal radiation (3rd degree burns)
        w.thermalRadius = 4000.0f * powf(yieldMT, 0.41f);

        // Prompt radiation (lethal dose 500 rem)
        w.radiationRadius = 2000.0f * powf(yieldMT, 0.19f);

        // Fallout (ground burst only)
        if (w.altitude == 0) {
            w.falloutRadius = 25000.0f * powf(yieldMT, 0.25f); // Up to 25km for large weapons
        } else {
            w.falloutRadius = 0.0f; // Air bursts produce little local fallout
        }
    }

    void calculateNuclearWinter() {
        // Smoke generation from fires (primarily from cities)
        // ~5 megatons of smoke per 100 warheads on urban targets
        int urbanWarheads = 0;
        for (const auto& w : warheads) {
            if (w.detonated && w.yield > 50.0f) { // Large enough to start firestorms
                urbanWarheads++;
            }
        }

        winterStats.smokeMass = urbanWarheads * 0.05f; // MT of smoke

        // Temperature drop (simplified climate model)
        // Based on Robock et al. nuclear winter studies
        if (totalYield > 1000000.0f) {
            // Full-scale thermonuclear war (5000+ warheads)
            winterStats.temperatureDrop = 20.0f; // 20°C drop for years
            winterStats.yearsOfDarkness = 10;
            winterStats.ozoneLoss = 75.0f; // 75% ozone destruction
            winterStats.agricultureCollapse = 99.0f; // Total collapse
        } else if (totalYield > 100000.0f) {
            // Major nuclear war (500-1000 warheads)
            winterStats.temperatureDrop = 12.0f;
            winterStats.yearsOfDarkness = 7;
            winterStats.ozoneLoss = 50.0f;
            winterStats.agricultureCollapse = 90.0f;
        } else if (totalYield > 10000.0f) {
            // Regional nuclear war (50-100 warheads)
            winterStats.temperatureDrop = 5.0f;
            winterStats.yearsOfDarkness = 3;
            winterStats.ozoneLoss = 25.0f;
            winterStats.agricultureCollapse = 50.0f;
        } else if (totalYield > 1000.0f) {
            // Limited exchange (10-20 warheads)
            winterStats.temperatureDrop = 2.0f;
            winterStats.yearsOfDarkness = 1;
            winterStats.ozoneLoss = 10.0f;
            winterStats.agricultureCollapse = 20.0f;
        } else {
            // Tactical use
            winterStats.temperatureDrop = 0.5f;
            winterStats.yearsOfDarkness = 0;
            winterStats.ozoneLoss = 2.0f;
            winterStats.agricultureCollapse = 5.0f;
        }

        // Casualty estimates (simplified)
        // Immediate: blast, thermal, radiation
        immediateCasualties = detonatedCount * 0.5f; // 500k per warhead average

        // Radiation sickness (days to weeks)
        radiationCasualties = immediateCasualties * 0.5f;

        // Starvation from nuclear winter
        if (winterStats.temperatureDrop > 10.0f) {
            starvationCasualties = 5.0f; // Billions
        } else if (winterStats.temperatureDrop > 5.0f) {
            starvationCasualties = 2.0f;
        } else {
            starvationCasualties = winterStats.agricultureCollapse / 100.0f;
        }
    }

public:
    NuclearWar() {
        reset();
    }

    void reset() {
        warheads.clear();
        time = 0.0f;
        detonatedCount = 0;
        totalYield = 0.0f;
        immediateCasualties = 0.0f;
        radiationCasualties = 0.0f;
        starvationCasualties = 0.0f;
        globalFallout = 0.0f;
        empEffect = 0.0f;
        winterStats = NuclearWinterStats{};
    }

    void addWarhead(float yieldKT, float latitude, float longitude, float altitude, int type) {
        Warhead w;
        w.yield = yieldKT;
        w.lat = latitude;
        w.lon = longitude;
        w.altitude = altitude;
        w.detonated = false;
        w.time = 0.0f;
        w.flightTime = 25.0f + (rand() % 10); // ICBM flight time 25-35 minutes
        w.type = static_cast<WarheadType>(type);

        warheads.push_back(w);
    }

    // Preset scenarios
    void setLimitedStrike() {
        reset();
        // 10 warheads on major cities
        addWarhead(500, 40.7128f, -74.0060f, 0, 2); // New York
        addWarhead(500, 34.0522f, -118.2437f, 0, 2); // Los Angeles
        addWarhead(500, 41.8781f, -87.6298f, 0, 2); // Chicago
        addWarhead(500, 51.5074f, -0.1278f, 0, 2); // London
        addWarhead(500, 48.8566f, 2.3522f, 0, 2); // Paris
        addWarhead(500, 52.5200f, 13.4050f, 0, 2); // Berlin
        addWarhead(500, 35.6762f, 139.6503f, 0, 2); // Tokyo
        addWarhead(500, 37.5665f, 126.9780f, 0, 2); // Seoul
        addWarhead(500, 39.9042f, 116.4074f, 0, 2); // Beijing
        addWarhead(500, 55.7558f, 37.6173f, 0, 2); // Moscow
    }

    void setRegionalWar() {
        reset();
        // India-Pakistan scenario (100 warheads)
        for (int i = 0; i < 50; i++) {
            float lat = 20.0f + (rand() % 20);
            float lon = 70.0f + (rand() % 20);
            addWarhead(100, lat, lon, 0, 1);
        }
        for (int i = 0; i < 50; i++) {
            float lat = 25.0f + (rand() % 15);
            float lon = 60.0f + (rand() % 20);
            addWarhead(100, lat, lon, 0, 1);
        }
    }

    void setFullScaleWar() {
        reset();
        // 1000 warheads - USA vs Russia strategic exchange

        // US strikes on Russia (500 warheads)
        for (int i = 0; i < 500; i++) {
            float lat = 45.0f + (rand() % 30);
            float lon = 30.0f + (rand() % 120);
            int type = (i % 3 == 0) ? 3 : 2; // Mix of MIRV and strategic
            addWarhead(350, lat, lon, 0, type);
        }

        // Russian strikes on USA/Europe (500 warheads)
        for (int i = 0; i < 250; i++) {
            float lat = 30.0f + (rand() % 20);
            float lon = -125.0f + (rand() % 50);
            addWarhead(500, lat, lon, 0, 2);
        }
        for (int i = 0; i < 250; i++) {
            float lat = 40.0f + (rand() % 25);
            float lon = -10.0f + (rand() % 40);
            addWarhead(500, lat, lon, 0, 2);
        }
    }

    void setArmageddon() {
        reset();
        // 5000+ warheads - total global thermonuclear war
        // All nuclear powers: USA, Russia, China, UK, France, India, Pakistan, Israel, NK

        // Massive strikes on all major population centers
        for (int i = 0; i < 5000; i++) {
            float lat = -60.0f + (rand() % 120);
            float lon = -180.0f + (rand() % 360);
            int type = rand() % 4;
            float yieldKT = 100 + (rand() % 400);
            addWarhead(yieldKT, lat, lon, 0, type);
        }
    }

    void update(float dt) {
        time += dt;

        for (auto& w : warheads) {
            if (!w.detonated) {
                w.time += dt;

                // Detonate when flight time reached
                if (w.time >= w.flightTime) {
                    w.detonated = true;
                    detonatedCount++;
                    totalYield += w.yield;
                    calculateBlastEffects(w);
                }
            }
        }

        // Calculate global effects after detonations
        if (detonatedCount > 0) {
            calculateNuclearWinter();
        }
    }

    // Getters
    int getWarheadCount() const { return warheads.size(); }
    int getDetonatedCount() const { return detonatedCount; }
    float getTotalYield() const { return totalYield; }
    float getTime() const { return time; }

    // Get warhead data for visualization
    std::vector<float> getWarheadPositions() const {
        std::vector<float> positions;
        for (const auto& w : warheads) {
            positions.push_back(w.lat);
            positions.push_back(w.lon);
            positions.push_back(w.detonated ? 1.0f : 0.0f);
            positions.push_back(w.time / w.flightTime); // Launch progress
        }
        return positions;
    }

    std::vector<float> getBlastRadii() const {
        std::vector<float> radii;
        for (const auto& w : warheads) {
            if (w.detonated) {
                radii.push_back(w.lat);
                radii.push_back(w.lon);
                radii.push_back(w.blastRadius_20psi);
                radii.push_back(w.blastRadius_5psi);
                radii.push_back(w.thermalRadius);
                radii.push_back(w.fireballRadius);
            }
        }
        return radii;
    }

    // Casualty estimates
    float getImmediateCasualties() const { return immediateCasualties; }
    float getRadiationCasualties() const { return radiationCasualties; }
    float getStarvationCasualties() const { return starvationCasualties; }
    float getTotalCasualties() const {
        return immediateCasualties + radiationCasualties + starvationCasualties;
    }

    // Nuclear winter stats
    float getTemperatureDrop() const { return winterStats.temperatureDrop; }
    int getYearsOfDarkness() const { return winterStats.yearsOfDarkness; }
    float getOzoneLoss() const { return winterStats.ozoneLoss; }
    float getAgricultureCollapse() const { return winterStats.agricultureCollapse; }
    float getSmokeMass() const { return winterStats.smokeMass; }
};

EMSCRIPTEN_BINDINGS(nuclear_module) {
    class_<NuclearWar>("NuclearWar")
        .constructor<>()
        .function("reset", &NuclearWar::reset)
        .function("update", &NuclearWar::update)
        .function("addWarhead", &NuclearWar::addWarhead)
        .function("setLimitedStrike", &NuclearWar::setLimitedStrike)
        .function("setRegionalWar", &NuclearWar::setRegionalWar)
        .function("setFullScaleWar", &NuclearWar::setFullScaleWar)
        .function("setArmageddon", &NuclearWar::setArmageddon)
        .function("getWarheadCount", &NuclearWar::getWarheadCount)
        .function("getDetonatedCount", &NuclearWar::getDetonatedCount)
        .function("getTotalYield", &NuclearWar::getTotalYield)
        .function("getTime", &NuclearWar::getTime)
        .function("getWarheadPositions", &NuclearWar::getWarheadPositions)
        .function("getBlastRadii", &NuclearWar::getBlastRadii)
        .function("getImmediateCasualties", &NuclearWar::getImmediateCasualties)
        .function("getRadiationCasualties", &NuclearWar::getRadiationCasualties)
        .function("getStarvationCasualties", &NuclearWar::getStarvationCasualties)
        .function("getTotalCasualties", &NuclearWar::getTotalCasualties)
        .function("getTemperatureDrop", &NuclearWar::getTemperatureDrop)
        .function("getYearsOfDarkness", &NuclearWar::getYearsOfDarkness)
        .function("getOzoneLoss", &NuclearWar::getOzoneLoss)
        .function("getAgricultureCollapse", &NuclearWar::getAgricultureCollapse)
        .function("getSmokeMass", &NuclearWar::getSmokeMass);

    register_vector<float>("VectorFloat");
}
