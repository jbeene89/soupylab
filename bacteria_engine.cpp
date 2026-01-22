/**
 * PLASTIC-EATING BACTERIA SIMULATOR
 * Models bacterial degradation of ocean plastic pollution
 *
 * Features:
 * - Temperature-dependent bacteria growth
 * - Weather pattern effects on bacterial activity
 * - Plastic degradation mechanics
 * - Colony spreading and competition
 * - Multiple bacteria species with different properties
 */

#include <emscripten/bind.h>
#include <emscripten/val.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>

using namespace emscripten;

// ============================================================================
// CONSTANTS
// ============================================================================

constexpr float OPTIMAL_TEMP = 25.0f;      // Celsius - optimal growth temp
constexpr float MIN_TEMP = 0.0f;           // Below this, bacteria dormant
constexpr float MAX_TEMP = 45.0f;          // Above this, bacteria die
constexpr float BASE_GROWTH_RATE = 0.05f;  // Population growth per update
constexpr float DEGRADATION_RATE = 0.01f;  // Plastic consumed per bacteria

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

    float length() const { return sqrtf(x * x + y * y); }
    float lengthSq() const { return x * x + y * y; }
};

// ============================================================================
// WEATHER SYSTEM
// ============================================================================

enum class Weather {
    SUNNY,      // High UV, optimal for some bacteria
    CLOUDY,     // Moderate conditions
    STORMY,     // Turbulent, spreads bacteria
    POLLUTED    // Toxic conditions, inhibits growth
};

struct WeatherZone {
    Vec2 center;
    float radius;
    Weather type;
    float intensity;  // 0.0 - 1.0

    WeatherZone(Vec2 c, float r, Weather t, float i = 1.0f)
        : center(c), radius(r), type(t), intensity(i) {}
};

// ============================================================================
// TEMPERATURE SYSTEM
// ============================================================================

enum class TempZone {
    TROPICAL,   // 25-30°C - fastest growth
    TEMPERATE,  // 15-25°C - moderate growth
    POLAR       // 0-15°C - slow/dormant
};

struct TemperatureField {
    Vec2 center;
    float radius;
    float temperature;  // Celsius
    TempZone zone;

    TemperatureField(Vec2 c, float r, float temp, TempZone z)
        : center(c), radius(r), temperature(temp), zone(z) {}
};

// ============================================================================
// BACTERIA COLONY
// ============================================================================

enum class BacteriaSpecies {
    IDEONELLA_SAKAIENSIS,  // Real plastic-eater, slow but efficient
    PSEUDOMONAS_PUTIDA,    // Fast grower, moderate efficiency
    ENGINEERED_SUPERBUG    // Fictional fast plastic degrader
};

struct BacteriaColony {
    Vec2 pos;
    float population;       // Bacterial count (millions)
    BacteriaSpecies species;
    float plasticConsumed;  // Total plastic degraded (kg)
    float growthRate;       // Current growth multiplier
    bool active;

    BacteriaColony()
        : population(1.0f)
        , species(BacteriaSpecies::IDEONELLA_SAKAIENSIS)
        , plasticConsumed(0.0f)
        , growthRate(1.0f)
        , active(true)
    {}

    BacteriaColony(Vec2 p, float pop, BacteriaSpecies s)
        : pos(p)
        , population(pop)
        , species(s)
        , plasticConsumed(0.0f)
        , growthRate(1.0f)
        , active(true)
    {}

    // Get species degradation efficiency
    float getDegradationEfficiency() const {
        switch (species) {
            case BacteriaSpecies::IDEONELLA_SAKAIENSIS:
                return 1.0f;  // 100% efficiency (real bacteria)
            case BacteriaSpecies::PSEUDOMONAS_PUTIDA:
                return 0.7f;  // 70% efficiency
            case BacteriaSpecies::ENGINEERED_SUPERBUG:
                return 2.0f;  // 200% efficiency (engineered)
            default:
                return 1.0f;
        }
    }

    // Get species growth rate
    float getBaseGrowthRate() const {
        switch (species) {
            case BacteriaSpecies::IDEONELLA_SAKAIENSIS:
                return 0.8f;  // Slower growth
            case BacteriaSpecies::PSEUDOMONAS_PUTIDA:
                return 1.5f;  // Fast growth
            case BacteriaSpecies::ENGINEERED_SUPERBUG:
                return 1.2f;  // Balanced
            default:
                return 1.0f;
        }
    }
};

// ============================================================================
// PLASTIC PARTICLE (simplified)
// ============================================================================

struct PlasticParticle {
    Vec2 pos;
    float mass;  // kg of plastic
    bool active;

    PlasticParticle(Vec2 p, float m) : pos(p), mass(m), active(true) {}
};

// ============================================================================
// BACTERIA ENGINE
// ============================================================================

class BacteriaEngine {
private:
    std::vector<BacteriaColony> colonies;
    std::vector<PlasticParticle> plasticParticles;
    std::vector<TemperatureField> tempFields;
    std::vector<WeatherZone> weatherZones;

    float worldSize;
    float time;
    int iteration;

    // Global environmental parameters
    float globalTemp;
    Weather globalWeather;
    float uvIntensity;

    std::mt19937 rng;

    // Calculate temperature at position
    float getTemperatureAt(const Vec2& pos) const {
        float temp = globalTemp;

        for (const auto& field : tempFields) {
            Vec2 diff = pos - field.center;
            float dist = diff.length();

            if (dist < field.radius) {
                // Blend temperature based on distance
                float influence = 1.0f - (dist / field.radius);
                temp = temp * (1.0f - influence) + field.temperature * influence;
            }
        }

        return temp;
    }

    // Get weather at position
    Weather getWeatherAt(const Vec2& pos, float& intensity) const {
        intensity = 1.0f;
        Weather weather = globalWeather;

        for (const auto& zone : weatherZones) {
            Vec2 diff = pos - zone.center;
            float dist = diff.length();

            if (dist < zone.radius) {
                weather = zone.type;
                intensity = zone.intensity;
                return weather;
            }
        }

        return weather;
    }

    // Calculate growth rate based on temperature
    float getTempGrowthMultiplier(float temp) const {
        if (temp < MIN_TEMP || temp > MAX_TEMP) {
            return 0.0f;  // Dormant/dead
        }

        // Gaussian curve around optimal temperature
        float diff = temp - OPTIMAL_TEMP;
        float variance = 100.0f;  // Temperature tolerance
        return expf(-(diff * diff) / (2.0f * variance));
    }

    // Calculate growth rate based on weather
    float getWeatherGrowthMultiplier(Weather weather, float intensity) const {
        switch (weather) {
            case Weather::SUNNY:
                return 1.2f * intensity;  // UV helps some bacteria
            case Weather::CLOUDY:
                return 1.0f;  // Neutral
            case Weather::STORMY:
                return 0.7f;  // Turbulent, harder to grow
            case Weather::POLLUTED:
                return 0.3f * (1.0f - intensity * 0.5f);  // Toxic
            default:
                return 1.0f;
        }
    }

    // Update bacteria growth
    void updateBacteriaGrowth() {
        for (auto& colony : colonies) {
            if (!colony.active || colony.population <= 0.0f) {
                colony.active = false;
                continue;
            }

            // Environmental factors
            float temp = getTemperatureAt(colony.pos);
            float weatherIntensity;
            Weather weather = getWeatherAt(colony.pos, weatherIntensity);

            // Calculate growth multipliers
            float tempMult = getTempGrowthMultiplier(temp);
            float weatherMult = getWeatherGrowthMultiplier(weather, weatherIntensity);
            float speciesMult = colony.getBaseGrowthRate();

            // Combined growth rate
            colony.growthRate = tempMult * weatherMult * speciesMult;

            // Logistic growth (limited by resources)
            float carryingCapacity = 1000.0f;  // Max population (millions)
            float growthFactor = 1.0f - (colony.population / carryingCapacity);

            if (growthFactor > 0.0f) {
                float growth = colony.population * BASE_GROWTH_RATE * colony.growthRate * growthFactor;
                colony.population += growth;
            }

            // Die off if conditions too extreme
            if (colony.growthRate < 0.1f) {
                colony.population *= 0.98f;  // 2% die-off per update
            }
        }
    }

    // Degrade plastic with bacteria
    void degradePlastic() {
        for (auto& colony : colonies) {
            if (!colony.active) continue;

            // Find nearby plastic
            for (auto& plastic : plasticParticles) {
                if (!plastic.active || plastic.mass <= 0.0f) {
                    plastic.active = false;
                    continue;
                }

                Vec2 diff = plastic.pos - colony.pos;
                float dist = diff.length();

                // Bacteria can degrade plastic within range
                float degradationRange = 50.0f;

                if (dist < degradationRange) {
                    // Amount degraded depends on population and efficiency
                    float efficiency = colony.getDegradationEfficiency();
                    float degraded = colony.population * DEGRADATION_RATE * efficiency * colony.growthRate;

                    degraded = fminf(degraded, plastic.mass);
                    plastic.mass -= degraded;
                    colony.plasticConsumed += degraded;

                    // Bacteria grow when consuming plastic
                    colony.population += degraded * 10.0f;  // Nutrient boost

                    if (plastic.mass <= 0.001f) {
                        plastic.active = false;
                    }
                }
            }
        }
    }

    // Spread bacteria to new locations
    void spreadBacteria() {
        std::uniform_real_distribution<float> dist(0, 1);
        std::vector<BacteriaColony> newColonies;

        for (const auto& colony : colonies) {
            if (!colony.active) continue;

            // Large colonies have chance to spread
            if (colony.population > 500.0f && dist(rng) < 0.01f) {
                // Spread to nearby location
                float angle = dist(rng) * 2 * 3.14159f;
                float spreadDist = 20.0f + dist(rng) * 50.0f;

                Vec2 newPos = colony.pos + Vec2(
                    cosf(angle) * spreadDist,
                    sinf(angle) * spreadDist
                );

                // Create daughter colony
                float splitPop = colony.population * 0.1f;  // 10% splits off
                newColonies.emplace_back(newPos, splitPop, colony.species);
            }
        }

        // Add new colonies
        for (auto& newColony : newColonies) {
            colonies.push_back(newColony);
        }
    }

public:
    BacteriaEngine(float worldSize)
        : worldSize(worldSize)
        , time(0)
        , iteration(0)
        , globalTemp(20.0f)
        , globalWeather(Weather::CLOUDY)
        , uvIntensity(0.5f)
        , rng(54321)
    {
        // Create default temperature zones
        createDefaultTempZones();
    }

    // ========================================================================
    // SETUP METHODS
    // ========================================================================

    void createDefaultTempZones() {
        tempFields.clear();

        // Tropical zone (center)
        tempFields.emplace_back(Vec2(0, 0), worldSize * 0.3f, 28.0f, TempZone::TROPICAL);

        // Temperate zones (mid)
        tempFields.emplace_back(Vec2(0, worldSize * 0.5f), worldSize * 0.25f, 18.0f, TempZone::TEMPERATE);
        tempFields.emplace_back(Vec2(0, -worldSize * 0.5f), worldSize * 0.25f, 18.0f, TempZone::TEMPERATE);

        // Polar zones (edges)
        tempFields.emplace_back(Vec2(0, worldSize * 0.8f), worldSize * 0.2f, 5.0f, TempZone::POLAR);
        tempFields.emplace_back(Vec2(0, -worldSize * 0.8f), worldSize * 0.2f, 5.0f, TempZone::POLAR);
    }

    void addBacteriaColony(float x, float y, int species, float initialPop) {
        BacteriaSpecies sp = static_cast<BacteriaSpecies>(species);
        colonies.emplace_back(Vec2(x, y), initialPop, sp);
    }

    void addPlasticParticle(float x, float y, float mass) {
        plasticParticles.emplace_back(Vec2(x, y), mass);
    }

    void setGlobalTemp(float temp) {
        globalTemp = temp;
    }

    void setGlobalWeather(int weather) {
        globalWeather = static_cast<Weather>(weather);
    }

    void addTempZone(float x, float y, float radius, float temp) {
        TempZone zone;
        if (temp >= 25.0f) zone = TempZone::TROPICAL;
        else if (temp >= 15.0f) zone = TempZone::TEMPERATE;
        else zone = TempZone::POLAR;

        tempFields.emplace_back(Vec2(x, y), radius, temp, zone);
    }

    void addWeatherZone(float x, float y, float radius, int weather, float intensity) {
        Weather w = static_cast<Weather>(weather);
        weatherZones.emplace_back(Vec2(x, y), radius, w, intensity);
    }

    void clearWeatherZones() {
        weatherZones.clear();
    }

    // ========================================================================
    // UPDATE
    // ========================================================================

    void update() {
        updateBacteriaGrowth();
        degradePlastic();

        // Spread bacteria occasionally
        if (iteration % 10 == 0) {
            spreadBacteria();
        }

        time += 1.0f;
        iteration++;
    }

    // ========================================================================
    // GETTERS
    // ========================================================================

    val getBacteriaData() {
        val data = val::array();
        for (const auto& colony : colonies) {
            if (colony.active && colony.population > 0.1f) {
                data.call<void>("push", colony.pos.x);
                data.call<void>("push", colony.pos.y);
                data.call<void>("push", colony.population);
                data.call<void>("push", static_cast<int>(colony.species));
                data.call<void>("push", colony.growthRate);
                data.call<void>("push", colony.plasticConsumed);
            }
        }
        return data;
    }

    val getPlasticData() {
        val data = val::array();
        for (const auto& plastic : plasticParticles) {
            if (plastic.active && plastic.mass > 0.001f) {
                data.call<void>("push", plastic.pos.x);
                data.call<void>("push", plastic.pos.y);
                data.call<void>("push", plastic.mass);
            }
        }
        return data;
    }

    int getBacteriaCount() const {
        int count = 0;
        for (const auto& c : colonies) {
            if (c.active && c.population > 0.1f) count++;
        }
        return count;
    }

    int getPlasticCount() const {
        int count = 0;
        for (const auto& p : plasticParticles) {
            if (p.active && p.mass > 0.001f) count++;
        }
        return count;
    }

    float getTotalBacteriaPopulation() const {
        float total = 0;
        for (const auto& c : colonies) {
            if (c.active) total += c.population;
        }
        return total;
    }

    float getTotalPlasticMass() const {
        float total = 0;
        for (const auto& p : plasticParticles) {
            if (p.active) total += p.mass;
        }
        return total;
    }

    float getTotalPlasticDegraded() const {
        float total = 0;
        for (const auto& c : colonies) {
            total += c.plasticConsumed;
        }
        return total;
    }

    float getTime() const { return time; }
};

// ============================================================================
// EMSCRIPTEN BINDINGS
// ============================================================================

EMSCRIPTEN_BINDINGS(bacteria_module) {
    class_<BacteriaEngine>("BacteriaEngine")
        .constructor<float>()
        .function("addBacteriaColony", &BacteriaEngine::addBacteriaColony)
        .function("addPlasticParticle", &BacteriaEngine::addPlasticParticle)
        .function("setGlobalTemp", &BacteriaEngine::setGlobalTemp)
        .function("setGlobalWeather", &BacteriaEngine::setGlobalWeather)
        .function("addTempZone", &BacteriaEngine::addTempZone)
        .function("addWeatherZone", &BacteriaEngine::addWeatherZone)
        .function("clearWeatherZones", &BacteriaEngine::clearWeatherZones)
        .function("update", &BacteriaEngine::update)
        .function("getBacteriaData", &BacteriaEngine::getBacteriaData)
        .function("getPlasticData", &BacteriaEngine::getPlasticData)
        .function("getBacteriaCount", &BacteriaEngine::getBacteriaCount)
        .function("getPlasticCount", &BacteriaEngine::getPlasticCount)
        .function("getTotalBacteriaPopulation", &BacteriaEngine::getTotalBacteriaPopulation)
        .function("getTotalPlasticMass", &BacteriaEngine::getTotalPlasticMass)
        .function("getTotalPlasticDegraded", &BacteriaEngine::getTotalPlasticDegraded)
        .function("getTime", &BacteriaEngine::getTime);
}
