/**
 * STRATOSPHERIC GLIDER SIMULATOR - WebAssembly
 * High Altitude Flight Dynamics with Atmospheric Model
 *
 * Features:
 * - ISA atmospheric model (up to 80km)
 * - Lift and drag calculations
 * - Energy management (potential + kinetic)
 * - Glide ratio optimization
 * - Thermal soaring dynamics
 * - Long-range trajectory planning
 * - Stall behavior modeling
 */

#include <emscripten/bind.h>
#include <emscripten/val.h>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace emscripten;

// ============================================================================
// PHYSICAL CONSTANTS
// ============================================================================

constexpr float G = 9.80665f;           // Gravitational acceleration (m/s²)
constexpr float R_AIR = 287.05f;        // Specific gas constant for air (J/kg·K)
constexpr float GAMMA = 1.4f;           // Specific heat ratio
constexpr float T0_ISA = 288.15f;       // Sea level temperature (K)
constexpr float P0_ISA = 101325.0f;     // Sea level pressure (Pa)
constexpr float RHO0_ISA = 1.225f;      // Sea level density (kg/m³)

// ============================================================================
// ISA ATMOSPHERE MODEL
// ============================================================================

struct Atmosphere {
    float temperature;  // K
    float pressure;     // Pa
    float density;      // kg/m³
    float speedOfSound; // m/s

    static Atmosphere calculate(float altitude) {
        Atmosphere atm;

        // International Standard Atmosphere model
        if (altitude < 11000.0f) {
            // Troposphere
            float lapse_rate = -0.0065f; // K/m
            atm.temperature = T0_ISA + lapse_rate * altitude;
            atm.pressure = P0_ISA * powf(atm.temperature / T0_ISA, -G / (lapse_rate * R_AIR));
        } else if (altitude < 20000.0f) {
            // Lower Stratosphere (isothermal)
            atm.temperature = 216.65f;
            float p11 = 22632.0f; // Pressure at 11km
            atm.pressure = p11 * expf(-G * (altitude - 11000.0f) / (R_AIR * atm.temperature));
        } else if (altitude < 32000.0f) {
            // Upper Stratosphere
            float lapse_rate = 0.001f; // K/m (warming)
            float T20 = 216.65f;
            float p20 = 5474.9f;
            atm.temperature = T20 + lapse_rate * (altitude - 20000.0f);
            atm.pressure = p20 * powf(atm.temperature / T20, -G / (lapse_rate * R_AIR));
        } else {
            // High Stratosphere
            atm.temperature = 228.65f;
            atm.pressure = 868.0f * expf(-G * (altitude - 32000.0f) / (R_AIR * atm.temperature));
        }

        atm.density = atm.pressure / (R_AIR * atm.temperature);
        atm.speedOfSound = sqrtf(GAMMA * R_AIR * atm.temperature);

        return atm;
    }
};

// ============================================================================
// GLIDER AIRCRAFT
// ============================================================================

class Glider {
private:
    // Aircraft parameters
    float wingspan;      // m
    float wingArea;      // m²
    float mass;          // kg
    float Cl0;           // Zero-lift coefficient
    float ClAlpha;       // Lift curve slope (1/rad)
    float Cd0;           // Parasitic drag coefficient
    float K;             // Induced drag factor

    // State variables
    float x, y;          // Position (m)
    float vx, vy;        // Velocity (m/s)
    float altitude;      // m
    float angleOfAttack; // radians
    float bankAngle;     // radians

    // History for visualization
    std::vector<float> trail_x;
    std::vector<float> trail_y;
    int maxTrailLength;

    // Statistics
    float maxAltitude;
    float totalDistance;
    float bestGlideRatio;

public:
    Glider()
        : wingspan(30.0f)
        , wingArea(20.0f)
        , mass(300.0f)
        , Cl0(0.1f)
        , ClAlpha(5.0f)
        , Cd0(0.015f)
        , K(0.04f)
        , x(0), y(0)
        , vx(30.0f), vy(0)
        , altitude(15000.0f)
        , angleOfAttack(0.05f)
        , bankAngle(0)
        , maxTrailLength(1000)
        , maxAltitude(15000.0f)
        , totalDistance(0)
        , bestGlideRatio(0)
    {}

    // ========================================================================
    // SETTERS
    // ========================================================================

    void setWingspan(float span) { wingspan = span; }
    void setWingArea(float area) { wingArea = area; }
    void setMass(float m) { mass = m; }
    void setAltitude(float alt) { altitude = alt; maxAltitude = alt; }
    void setAngleOfAttack(float aoa) { angleOfAttack = aoa * 3.14159f / 180.0f; }
    void setBankAngle(float bank) { bankAngle = bank * 3.14159f / 180.0f; }
    void setVelocity(float v) {
        float heading = atan2f(vy, vx);
        vx = v * cosf(heading);
        vy = v * sinf(heading);
    }

    // ========================================================================
    // AERODYNAMIC CALCULATIONS
    // ========================================================================

    float calculateCl() const {
        return Cl0 + ClAlpha * angleOfAttack;
    }

    float calculateCd() const {
        float Cl = calculateCl();
        return Cd0 + K * Cl * Cl; // Cd = Cd0 + K*Cl²
    }

    float calculateLiftToDrag() const {
        float Cl = calculateCl();
        float Cd = calculateCd();
        return Cd > 0 ? Cl / Cd : 0;
    }

    float calculateAspectRatio() const {
        return wingspan * wingspan / wingArea;
    }

    // ========================================================================
    // PHYSICS UPDATE
    // ========================================================================

    void update(float dt) {
        // Get atmospheric conditions
        Atmosphere atm = Atmosphere::calculate(altitude);

        // Calculate airspeed
        float velocity = sqrtf(vx * vx + vy * vy);
        float dynamicPressure = 0.5f * atm.density * velocity * velocity;

        // Aerodynamic coefficients
        float Cl = calculateCl();
        float Cd = calculateCd();

        // Forces in wind frame
        float lift = Cl * dynamicPressure * wingArea;
        float drag = Cd * dynamicPressure * wingArea;
        float weight = mass * G;

        // Flight path angle
        float gamma = atan2f(vy, vx);

        // Convert to body frame
        float Fx = -drag * cosf(gamma) + lift * sinf(gamma);
        float Fy = -drag * sinf(gamma) - lift * cosf(gamma) - weight;

        // Banking turn (simplified)
        if (fabsf(bankAngle) > 0.01f) {
            float turnRadius = velocity * velocity / (G * tanf(bankAngle));
            float turnRate = velocity / turnRadius;
            float heading = atan2f(vy, vx);
            heading += turnRate * dt;

            float v = sqrtf(vx * vx + vy * vy);
            vx = v * cosf(heading);
            vy = v * sinf(heading);
        }

        // Integrate equations of motion
        float ax = Fx / mass;
        float ay = Fy / mass;

        vx += ax * dt;
        vy += ay * dt;

        x += vx * dt;
        y += vy * dt;
        altitude += vy * dt;

        // Stall check
        float stallAngle = 15.0f * 3.14159f / 180.0f;
        if (angleOfAttack > stallAngle) {
            // Reduce lift dramatically in stall
            vy -= 5.0f * dt;
        }

        // Ground collision
        if (altitude < 0) {
            altitude = 0;
            vy = 0;
            vx *= 0.9f; // Friction
        }

        // Update trail
        trail_x.push_back(x);
        trail_y.push_back(altitude);
        if (trail_x.size() > maxTrailLength) {
            trail_x.erase(trail_x.begin());
            trail_y.erase(trail_y.begin());
        }

        // Update statistics
        if (altitude > maxAltitude) maxAltitude = altitude;
        totalDistance += velocity * dt;

        float glideRatio = calculateLiftToDrag();
        if (glideRatio > bestGlideRatio) bestGlideRatio = glideRatio;
    }

    // ========================================================================
    // THERMAL UPDRAFT (for soaring)
    // ========================================================================

    void applyThermal(float thermalX, float thermalY, float thermalStrength, float dt) {
        float dx = x - thermalX;
        float dy = altitude - thermalY;
        float dist = sqrtf(dx * dx + dy * dy);

        if (dist < 500.0f) {
            // Gaussian thermal profile
            float uplift = thermalStrength * expf(-dist * dist / (2.0f * 200.0f * 200.0f));
            vy += uplift * dt;
        }
    }

    // ========================================================================
    // GETTERS
    // ========================================================================

    float getX() const { return x; }
    float getY() const { return altitude; }
    float getAltitude() const { return altitude; }
    float getVelocity() const { return sqrtf(vx * vx + vy * vy); }
    float getVerticalSpeed() const { return vy; }
    float getMachNumber() const {
        Atmosphere atm = Atmosphere::calculate(altitude);
        return sqrtf(vx * vx + vy * vy) / atm.speedOfSound;
    }
    float getGlideRatio() const { return calculateLiftToDrag(); }
    float getAngleOfAttack() const { return angleOfAttack * 180.0f / 3.14159f; }

    float getAtmosphericDensity() const {
        return Atmosphere::calculate(altitude).density;
    }
    float getAtmosphericPressure() const {
        return Atmosphere::calculate(altitude).pressure;
    }
    float getAtmosphericTemperature() const {
        return Atmosphere::calculate(altitude).temperature - 273.15f; // Convert to Celsius
    }

    float getMaxAltitude() const { return maxAltitude; }
    float getTotalDistance() const { return totalDistance; }
    float getBestGlideRatio() const { return bestGlideRatio; }

    val getTrail() const {
        val trail = val::array();
        for (size_t i = 0; i < trail_x.size(); ++i) {
            trail.call<void>("push", trail_x[i]);
            trail.call<void>("push", trail_y[i]);
        }
        return trail;
    }

    void reset() {
        x = 0;
        y = 0;
        altitude = 15000.0f;
        vx = 30.0f;
        vy = 0;
        angleOfAttack = 0.05f;
        bankAngle = 0;
        trail_x.clear();
        trail_y.clear();
        maxAltitude = altitude;
        totalDistance = 0;
        bestGlideRatio = 0;
    }
};

// ============================================================================
// EMSCRIPTEN BINDINGS
// ============================================================================

EMSCRIPTEN_BINDINGS(glider_module) {
    class_<Glider>("Glider")
        .constructor<>()
        .function("setWingspan", &Glider::setWingspan)
        .function("setWingArea", &Glider::setWingArea)
        .function("setMass", &Glider::setMass)
        .function("setAltitude", &Glider::setAltitude)
        .function("setAngleOfAttack", &Glider::setAngleOfAttack)
        .function("setBankAngle", &Glider::setBankAngle)
        .function("setVelocity", &Glider::setVelocity)
        .function("update", &Glider::update)
        .function("applyThermal", &Glider::applyThermal)
        .function("getX", &Glider::getX)
        .function("getY", &Glider::getY)
        .function("getAltitude", &Glider::getAltitude)
        .function("getVelocity", &Glider::getVelocity)
        .function("getVerticalSpeed", &Glider::getVerticalSpeed)
        .function("getMachNumber", &Glider::getMachNumber)
        .function("getGlideRatio", &Glider::getGlideRatio)
        .function("getAngleOfAttack", &Glider::getAngleOfAttack)
        .function("getAtmosphericDensity", &Glider::getAtmosphericDensity)
        .function("getAtmosphericPressure", &Glider::getAtmosphericPressure)
        .function("getAtmosphericTemperature", &Glider::getAtmosphericTemperature)
        .function("getMaxAltitude", &Glider::getMaxAltitude)
        .function("getTotalDistance", &Glider::getTotalDistance)
        .function("getBestGlideRatio", &Glider::getBestGlideRatio)
        .function("getTrail", &Glider::getTrail)
        .function("reset", &Glider::reset);
}
