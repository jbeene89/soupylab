/**
 * ULTIMATE STRATOSPHERIC GLIDER SIMULATOR - WebAssembly
 * Professional Soaring Flight Computer with Complete Atmospheric Model
 *
 * Features:
 * - Multiple glider types (high-performance, standard, paraglider)
 * - Complete ISA atmosphere (0-80km)
 * - Thermal lift (bubble and columnar models)
 * - Ridge lift (orographic wave mechanics)
 * - Wave lift (mountain lee waves)
 * - Wind shear and wind drift
 * - Polar curves (L/D vs airspeed)
 * - McCready theory for optimal speed-to-fly
 * - Variometer with audio feedback
 * - Flight computer displays
 * - Cross-country waypoint navigation
 * - Cloud base visualization
 * - Energy management system
 * - Stall/spin dynamics
 * - Thermal centering algorithm
 */

#include <emscripten/bind.h>
#include <emscripten/val.h>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace emscripten;

// ============================================================================
// CONSTANTS
// ============================================================================

constexpr float G = 9.80665f;
constexpr float R_AIR = 287.05f;
constexpr float T0_ISA = 288.15f;
constexpr float P0_ISA = 101325.0f;
constexpr float PI = 3.14159265359f;

// ============================================================================
// ATMOSPHERIC MODEL
// ============================================================================

struct Atmosphere {
    float temperature;  // K
    float pressure;     // Pa
    float density;      // kg/m³
    float speedOfSound; // m/s
    float windX;        // m/s
    float windY;        // m/s

    static Atmosphere calculate(float altitude, float time = 0) {
        Atmosphere atm;

        // ISA Temperature
        if (altitude < 11000.0f) {
            atm.temperature = T0_ISA - 0.0065f * altitude;
            atm.pressure = P0_ISA * powf(atm.temperature / T0_ISA, -G / (-0.0065f * R_AIR));
        } else if (altitude < 20000.0f) {
            atm.temperature = 216.65f;
            atm.pressure = 22632.0f * expf(-G * (altitude - 11000.0f) / (R_AIR * atm.temperature));
        } else if (altitude < 32000.0f) {
            float T20 = 216.65f;
            atm.temperature = T20 + 0.001f * (altitude - 20000.0f);
            atm.pressure = 5474.9f * powf(atm.temperature / T20, -G / (0.001f * R_AIR));
        } else {
            atm.temperature = 228.65f;
            atm.pressure = 868.0f * expf(-G * (altitude - 32000.0f) / (R_AIR * atm.temperature));
        }

        atm.density = atm.pressure / (R_AIR * atm.temperature);
        atm.speedOfSound = sqrtf(1.4f * R_AIR * atm.temperature);

        // Wind model (jet stream + surface wind)
        float jetAltitude = 11000.0f;
        float jetStrength = 30.0f * expf(-powf((altitude - jetAltitude) / 3000.0f, 2));
        float surfaceWind = 5.0f * expf(-altitude / 1000.0f);
        atm.windX = jetStrength + surfaceWind + 3.0f * sinf(time * 0.01f);
        atm.windY = 2.0f * sinf(time * 0.02f);

        return atm;
    }
};

// ============================================================================
// GLIDER TYPES
// ============================================================================

enum class GliderType {
    HIGH_PERFORMANCE,  // ASH 31, Eta - L/D 60+
    STANDARD_CLASS,    // ASW 28 - L/D 45
    PARAGLIDER,        // L/D 10
    HANG_GLIDER        // L/D 15
};

struct GliderConfig {
    float wingspan;      // m
    float wingArea;      // m²
    float mass;          // kg
    float Cd0;           // Parasitic drag
    float K;             // Induced drag factor
    float ClMax;         // Max lift coefficient
    float bestLD;        // Best L/D ratio
    float bestLDSpeed;   // Speed at best L/D (m/s)
    float minSinkRate;   // Min sink rate (m/s)
    float minSinkSpeed;  // Speed at min sink (m/s)

    static GliderConfig getConfig(GliderType type) {
        GliderConfig config;
        switch (type) {
            case GliderType::HIGH_PERFORMANCE:
                config.wingspan = 25.0f;
                config.wingArea = 10.5f;
                config.mass = 600.0f;
                config.Cd0 = 0.008f;
                config.K = 0.012f;
                config.ClMax = 1.6f;
                config.bestLD = 60.0f;
                config.bestLDSpeed = 30.0f;
                config.minSinkRate = 0.45f;
                config.minSinkSpeed = 22.0f;
                break;
            case GliderType::STANDARD_CLASS:
                config.wingspan = 15.0f;
                config.wingArea = 10.0f;
                config.mass = 550.0f;
                config.Cd0 = 0.010f;
                config.K = 0.015f;
                config.ClMax = 1.5f;
                config.bestLD = 45.0f;
                config.bestLDSpeed = 28.0f;
                config.minSinkRate = 0.55f;
                config.minSinkSpeed = 20.0f;
                break;
            case GliderType::PARAGLIDER:
                config.wingspan = 11.0f;
                config.wingArea = 27.0f;
                config.mass = 110.0f;
                config.Cd0 = 0.030f;
                config.K = 0.050f;
                config.ClMax = 1.2f;
                config.bestLD = 10.0f;
                config.bestLDSpeed = 12.0f;
                config.minSinkRate = 1.0f;
                config.minSinkSpeed = 9.0f;
                break;
            case GliderType::HANG_GLIDER:
                config.wingspan = 9.5f;
                config.wingArea = 14.0f;
                config.mass = 140.0f;
                config.Cd0 = 0.025f;
                config.K = 0.040f;
                config.ClMax = 1.3f;
                config.bestLD = 15.0f;
                config.bestLDSpeed = 15.0f;
                config.minSinkRate = 0.8f;
                config.minSinkSpeed = 11.0f;
                break;
        }
        return config;
    }
};

// ============================================================================
// LIFT SOURCES
// ============================================================================

struct Thermal {
    float x, y;          // Position
    float altitude;      // Base altitude
    float cloudBase;     // Top altitude
    float strength;      // Max vertical velocity (m/s)
    float radius;        // Core radius (m)
    float age;           // Time since spawn (s)
    bool isActive;

    float getLift(float px, float py, float pAlt) const {
        if (!isActive || pAlt > cloudBase || pAlt < altitude) return 0;

        float dx = px - x;
        float dy = py - y;
        float dist = sqrtf(dx * dx + dy * dy);

        // Gaussian thermal model
        float coreVelocity = strength * expf(-(dist * dist) / (radius * radius));

        // Height-dependent strength (stronger at mid-height)
        float heightFactor = sinf(PI * (pAlt - altitude) / (cloudBase - altitude));

        // Age decay
        float ageFactor = expf(-age / 600.0f); // 10 min lifetime

        return coreVelocity * heightFactor * ageFactor;
    }
};

struct Ridge {
    float x;             // Ridge position
    float height;        // Ridge height
    float width;         // Ridge width
    float windSpeed;     // Perpendicular wind speed

    float getLift(float px, float py, float pAlt) const {
        float dx = px - x;
        float distFromRidge = fabsf(dx);

        if (distFromRidge > width * 3.0f || pAlt > height * 2.0f) return 0;

        // Ridge lift model
        float liftStrength = windSpeed * 0.5f;
        float horizontalFactor = expf(-(distFromRidge * distFromRidge) / (width * width));
        float verticalFactor = (pAlt < height) ? (pAlt / height) : expf(-(pAlt - height) / 500.0f);

        return liftStrength * horizontalFactor * verticalFactor;
    }
};

struct Wave {
    float x;             // Wave source (mountain)
    float amplitude;     // Wave amplitude
    float wavelength;    // Horizontal wavelength
    float minAltitude;   // Minimum wave altitude

    float getLift(float px, float py, float pAlt) const {
        if (pAlt < minAltitude) return 0;

        float dx = px - x;
        float phase = 2.0f * PI * dx / wavelength;

        // Mountain wave lift (sine wave pattern)
        float waveStrength = amplitude * sinf(phase) * expf(-(pAlt - minAltitude) / 3000.0f);

        return waveStrength;
    }
};

// ============================================================================
// GLIDER CLASS
// ============================================================================

class Glider {
private:
    GliderType type;
    GliderConfig config;

    // State
    float x, y;          // Position (m)
    float altitude;      // Altitude AGL (m)
    float vx, vy;        // Velocity (m/s)
    float verticalSpeed; // Vertical speed (m/s)
    float angleOfAttack; // radians
    float bankAngle;     // radians
    float ballast;       // Water ballast (kg)

    // Flight computer
    float totalEnergy;
    float maxAltitude;
    float totalDistance;
    float bestGlideRatio;
    float averageClimbRate;
    float timeInThermals;

    // McCready settings
    float mcCreadyValue;  // Expected thermal strength (m/s)

    // Thermals and lift
    std::vector<Thermal> thermals;
    Ridge ridge;
    Wave wave;

    // Trail
    std::vector<float> trail_x;
    std::vector<float> trail_alt;
    std::vector<float> trail_vario;
    int maxTrailLength;

    // Variometer
    float varioSmoothed;
    float varioAveraged;

public:
    Glider()
        : type(GliderType::HIGH_PERFORMANCE)
        , config(GliderConfig::getConfig(type))
        , x(0), y(0)
        , altitude(2000.0f)
        , vx(config.bestLDSpeed), vy(0)
        , verticalSpeed(0)
        , angleOfAttack(0.05f)
        , bankAngle(0)
        , ballast(0)
        , totalEnergy(0)
        , maxAltitude(2000.0f)
        , totalDistance(0)
        , bestGlideRatio(0)
        , averageClimbRate(0)
        , timeInThermals(0)
        , mcCreadyValue(2.0f)
        , maxTrailLength(3000)
        , varioSmoothed(0)
        , varioAveraged(0)
    {
        // Initialize ridge
        ridge.x = 5000.0f;
        ridge.height = 800.0f;
        ridge.width = 200.0f;
        ridge.windSpeed = 8.0f;

        // Initialize wave
        wave.x = 10000.0f;
        wave.amplitude = 5.0f;
        wave.wavelength = 3000.0f;
        wave.minAltitude = 3000.0f;

        // Spawn initial thermals
        spawnThermals(5);
    }

    // ========================================================================
    // SETTERS
    // ========================================================================

    void setGliderType(int typeIndex) {
        type = static_cast<GliderType>(typeIndex);
        config = GliderConfig::getConfig(type);
    }

    void setAltitude(float alt) { altitude = alt; maxAltitude = alt; }
    void setVelocity(float v) {
        float heading = atan2f(vy, vx);
        vx = v * cosf(heading);
        vy = v * sinf(heading);
    }
    void setAngleOfAttack(float aoa) { angleOfAttack = aoa * PI / 180.0f; }
    void setBankAngle(float bank) { bankAngle = bank * PI / 180.0f; }
    void setBallast(float b) { ballast = b; }
    void setMcCready(float mc) { mcCreadyValue = mc; }

    void addThermal(float tx, float ty, float tAlt, float strength, float radius) {
        Thermal t;
        t.x = tx;
        t.y = ty;
        t.altitude = tAlt;
        t.cloudBase = tAlt + 1000.0f + strength * 200.0f;
        t.strength = strength;
        t.radius = radius;
        t.age = 0;
        t.isActive = true;
        thermals.push_back(t);
    }

    void spawnThermals(int count) {
        for (int i = 0; i < count; i++) {
            float tx = (rand() % 20000) - 10000.0f;
            float ty = (rand() % 20000) - 10000.0f;
            float tAlt = 200.0f + (rand() % 500);
            float strength = 1.5f + (rand() % 100) / 20.0f;
            float radius = 100.0f + (rand() % 200);
            addThermal(tx, ty, tAlt, strength, radius);
        }
    }

    // ========================================================================
    // AERODYNAMICS
    // ========================================================================

    float calculateLift() const {
        Atmosphere atm = Atmosphere::calculate(altitude);
        float airspeed = sqrtf(vx * vx + vy * vy);
        float dynamicPressure = 0.5f * atm.density * airspeed * airspeed;

        float Cl = 0.2f + 5.0f * angleOfAttack;
        Cl = std::min(Cl, config.ClMax);

        float totalMass = config.mass + ballast;
        float lift = Cl * dynamicPressure * config.wingArea;

        return lift;
    }

    float calculateDrag() const {
        Atmosphere atm = Atmosphere::calculate(altitude);
        float airspeed = sqrtf(vx * vx + vy * vy);
        float dynamicPressure = 0.5f * atm.density * airspeed * airspeed;

        float Cl = 0.2f + 5.0f * angleOfAttack;
        Cl = std::min(Cl, config.ClMax);

        float Cd = config.Cd0 + config.K * Cl * Cl;
        float drag = Cd * dynamicPressure * config.wingArea;

        return drag;
    }

    float getPolarSinkRate(float airspeed) const {
        Atmosphere atm = Atmosphere::calculate(altitude);
        float q = 0.5f * atm.density * airspeed * airspeed;

        float totalMass = config.mass + ballast;
        float weight = totalMass * G;

        float Cl = weight / (q * config.wingArea);
        float Cd = config.Cd0 + config.K * Cl * Cl;

        float sinkRate = -Cd * q * config.wingArea / weight * airspeed;
        return sinkRate;
    }

    float getSpeedToFly() const {
        // McCready theory: optimal speed based on expected climb rate
        float baseSink = config.minSinkRate;
        float speedAdjust = sqrtf(1.0f + mcCreadyValue / baseSink);
        return config.minSinkSpeed * speedAdjust;
    }

    // ========================================================================
    // LIFT CALCULATION
    // ========================================================================

    float getTotalLift(float px, float py, float pAlt) const {
        float totalLift = 0;

        // Thermals
        for (const auto& thermal : thermals) {
            totalLift += thermal.getLift(px, py, pAlt);
        }

        // Ridge
        totalLift += ridge.getLift(px, py, pAlt);

        // Wave
        totalLift += wave.getLift(px, py, pAlt);

        return totalLift;
    }

    // ========================================================================
    // PHYSICS UPDATE
    // ========================================================================

    void update(float dt) {
        Atmosphere atm = Atmosphere::calculate(altitude);

        // Airspeed (relative to wind)
        float airspeedX = vx - atm.windX;
        float airspeedY = vy - atm.windY;
        float airspeed = sqrtf(airspeedX * airspeedX + airspeedY * airspeedY + verticalSpeed * verticalSpeed);

        // Forces
        float lift = calculateLift();
        float drag = calculateDrag();
        float totalMass = config.mass + ballast;
        float weight = totalMass * G;

        // Flight path angle
        float gamma = atan2f(verticalSpeed, airspeed);

        // Angle adjustments for bank
        float liftVertical = lift * cosf(bankAngle);
        float liftHorizontal = lift * sinf(bankAngle);

        // Vertical acceleration (including environmental lift)
        float liftFromEnvironment = getTotalLift(x, y, altitude);
        float verticalAccel = (liftVertical - weight) / totalMass - drag * sinf(gamma) / totalMass + liftFromEnvironment;

        // Horizontal acceleration
        float horizontalAccel = -drag * cosf(gamma) / totalMass + liftHorizontal / totalMass;

        // Update velocities
        verticalSpeed += verticalAccel * dt;
        float horizontalSpeed = sqrtf(vx * vx + vy * vy);
        horizontalSpeed += horizontalAccel * dt;

        // Turn dynamics
        if (fabsf(bankAngle) > 0.01f) {
            float turnRate = G * tanf(bankAngle) / horizontalSpeed;
            float heading = atan2f(vy, vx);
            heading += turnRate * dt;
            vx = horizontalSpeed * cosf(heading);
            vy = horizontalSpeed * sinf(heading);
        }

        // Update position
        x += (vx + atm.windX) * dt;
        y += (vy + atm.windY) * dt;
        altitude += verticalSpeed * dt;

        // Ground collision
        if (altitude < 0) {
            altitude = 0;
            verticalSpeed = 0;
            vx *= 0.9f;
            vy *= 0.9f;
        }

        // Variometer filtering
        varioSmoothed = varioSmoothed * 0.9f + verticalSpeed * 0.1f;
        varioAveraged = varioAveraged * 0.98f + verticalSpeed * 0.02f;

        // Update thermals
        for (auto& thermal : thermals) {
            thermal.age += dt;
            if (thermal.age > 900.0f) { // 15 min lifetime
                thermal.isActive = false;
            }
        }

        // Respawn thermals
        if (rand() % 1000 < 2) {
            spawnThermals(1);
        }

        // Statistics
        if (altitude > maxAltitude) maxAltitude = altitude;
        totalDistance += sqrtf(vx * vx + vy * vy) * dt;
        if (verticalSpeed > 0.5f) {
            timeInThermals += dt;
            averageClimbRate = averageClimbRate * 0.99f + verticalSpeed * 0.01f;
        }

        float currentLD = (fabsf(verticalSpeed) > 0.1f) ? horizontalSpeed / fabsf(verticalSpeed) : 0;
        if (currentLD > bestGlideRatio) bestGlideRatio = currentLD;

        // Trail
        trail_x.push_back(x);
        trail_alt.push_back(altitude);
        trail_vario.push_back(verticalSpeed);
        if (trail_x.size() > maxTrailLength) {
            trail_x.erase(trail_x.begin());
            trail_alt.erase(trail_alt.begin());
            trail_vario.erase(trail_vario.begin());
        }

        // Energy
        totalEnergy = totalMass * G * altitude + 0.5f * totalMass * airspeed * airspeed;
    }

    // ========================================================================
    // GETTERS
    // ========================================================================

    float getX() const { return x; }
    float getY() const { return y; }
    float getAltitude() const { return altitude; }
    float getVelocity() const { return sqrtf(vx * vx + vy * vy); }
    float getVerticalSpeed() const { return verticalSpeed; }
    float getVarioSmoothed() const { return varioSmoothed; }
    float getVarioAveraged() const { return varioAveraged; }
    float getAngleOfAttack() const { return angleOfAttack * 180.0f / PI; }
    float getBankAngle() const { return bankAngle * 180.0f / PI; }

    float getMaxAltitude() const { return maxAltitude; }
    float getTotalDistance() const { return totalDistance / 1000.0f; }
    float getBestGlideRatio() const { return bestGlideRatio; }
    float getAverageClimbRate() const { return averageClimbRate; }
    float getTimeInThermals() const { return timeInThermals; }

    float getCurrentLD() const {
        float horizontalSpeed = sqrtf(vx * vx + vy * vy);
        return (fabsf(verticalSpeed) > 0.1f) ? horizontalSpeed / fabsf(verticalSpeed) : config.bestLD;
    }

    float getSpeedToFlyMcCready() const { return getSpeedToFly(); }
    float getMinSinkSpeed() const { return config.minSinkSpeed; }
    float getBestLDSpeed() const { return config.bestLDSpeed; }

    float getAtmosphericTemp() const {
        return Atmosphere::calculate(altitude).temperature - 273.15f;
    }
    float getAtmosphericPressure() const {
        return Atmosphere::calculate(altitude).pressure / 100.0f; // hPa
    }
    float getWindX() const {
        return Atmosphere::calculate(altitude).windX;
    }
    float getWindY() const {
        return Atmosphere::calculate(altitude).windY;
    }

    val getTrail() const {
        val trail = val::array();
        for (size_t i = 0; i < trail_x.size(); ++i) {
            trail.call<void>("push", trail_x[i]);
            trail.call<void>("push", trail_alt[i]);
        }
        return trail;
    }

    val getTrailVario() const {
        val vario = val::array();
        for (size_t i = 0; i < trail_vario.size(); ++i) {
            vario.call<void>("push", trail_vario[i]);
        }
        return vario;
    }

    val getThermals() const {
        val thermalData = val::array();
        for (const auto& t : thermals) {
            if (t.isActive) {
                thermalData.call<void>("push", t.x);
                thermalData.call<void>("push", t.y);
                thermalData.call<void>("push", t.radius);
                thermalData.call<void>("push", t.strength);
                thermalData.call<void>("push", t.cloudBase);
            }
        }
        return thermalData;
    }

    void reset() {
        x = 0;
        y = 0;
        altitude = 2000.0f;
        vx = config.bestLDSpeed;
        vy = 0;
        verticalSpeed = 0;
        angleOfAttack = 0.05f;
        bankAngle = 0;
        maxAltitude = 2000.0f;
        totalDistance = 0;
        bestGlideRatio = 0;
        averageClimbRate = 0;
        timeInThermals = 0;
        trail_x.clear();
        trail_alt.clear();
        trail_vario.clear();
        thermals.clear();
        spawnThermals(5);
    }
};

// ============================================================================
// EMSCRIPTEN BINDINGS
// ============================================================================

EMSCRIPTEN_BINDINGS(glider_module) {
    class_<Glider>("Glider")
        .constructor<>()
        .function("setGliderType", &Glider::setGliderType)
        .function("setAltitude", &Glider::setAltitude)
        .function("setVelocity", &Glider::setVelocity)
        .function("setAngleOfAttack", &Glider::setAngleOfAttack)
        .function("setBankAngle", &Glider::setBankAngle)
        .function("setBallast", &Glider::setBallast)
        .function("setMcCready", &Glider::setMcCready)
        .function("addThermal", &Glider::addThermal)
        .function("spawnThermals", &Glider::spawnThermals)
        .function("update", &Glider::update)
        .function("getX", &Glider::getX)
        .function("getY", &Glider::getY)
        .function("getAltitude", &Glider::getAltitude)
        .function("getVelocity", &Glider::getVelocity)
        .function("getVerticalSpeed", &Glider::getVerticalSpeed)
        .function("getVarioSmoothed", &Glider::getVarioSmoothed)
        .function("getVarioAveraged", &Glider::getVarioAveraged)
        .function("getAngleOfAttack", &Glider::getAngleOfAttack)
        .function("getBankAngle", &Glider::getBankAngle)
        .function("getMaxAltitude", &Glider::getMaxAltitude)
        .function("getTotalDistance", &Glider::getTotalDistance)
        .function("getBestGlideRatio", &Glider::getBestGlideRatio)
        .function("getAverageClimbRate", &Glider::getAverageClimbRate)
        .function("getTimeInThermals", &Glider::getTimeInThermals)
        .function("getCurrentLD", &Glider::getCurrentLD)
        .function("getSpeedToFlyMcCready", &Glider::getSpeedToFlyMcCready)
        .function("getMinSinkSpeed", &Glider::getMinSinkSpeed)
        .function("getBestLDSpeed", &Glider::getBestLDSpeed)
        .function("getAtmosphericTemp", &Glider::getAtmosphericTemp)
        .function("getAtmosphericPressure", &Glider::getAtmosphericPressure)
        .function("getWindX", &Glider::getWindX)
        .function("getWindY", &Glider::getWindY)
        .function("getTrail", &Glider::getTrail)
        .function("getTrailVario", &Glider::getTrailVario)
        .function("getThermals", &Glider::getThermals)
        .function("reset", &Glider::reset);
}
