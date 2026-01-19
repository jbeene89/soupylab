/**
 * HYPERSONIC FLIGHT SIMULATOR - WebAssembly
 * Mach 5+ Flight Dynamics with Aerothermal Effects
 *
 * Features:
 * - Hypersonic aerodynamics (Mach 5-25)
 * - Aerodynamic heating (convective + radiative)
 * - Real gas effects at high temperatures
 * - Thermal protection system (TPS)
 * - Atmospheric re-entry trajectories
 * - Plasma sheath formation
 * - Shock layer temperature
 * - Newtonian impact theory
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
constexpr float GAMMA = 1.4f;           // Specific heat ratio (low temp)
constexpr float GAMMA_HOT = 1.2f;       // Specific heat ratio (high temp, dissociated)
constexpr float T0_ISA = 288.15f;       // Sea level temperature (K)
constexpr float P0_ISA = 101325.0f;     // Sea level pressure (Pa)
constexpr float RHO0_ISA = 1.225f;      // Sea level density (kg/m³)
constexpr float EARTH_RADIUS = 6371000.0f; // m
constexpr float STEFAN_BOLTZMANN = 5.67e-8f; // W/(m²·K⁴)

// ============================================================================
// ISA ATMOSPHERE MODEL (Extended for High Altitude)
// ============================================================================

struct Atmosphere {
    float temperature;  // K
    float pressure;     // Pa
    float density;      // kg/m³
    float speedOfSound; // m/s
    float meanFreePath; // m (for rarefaction effects)

    static Atmosphere calculate(float altitude) {
        Atmosphere atm;

        if (altitude < 11000.0f) {
            // Troposphere
            float lapse_rate = -0.0065f;
            atm.temperature = T0_ISA + lapse_rate * altitude;
            atm.pressure = P0_ISA * powf(atm.temperature / T0_ISA, -G / (lapse_rate * R_AIR));
        } else if (altitude < 20000.0f) {
            // Lower Stratosphere
            atm.temperature = 216.65f;
            float p11 = 22632.0f;
            atm.pressure = p11 * expf(-G * (altitude - 11000.0f) / (R_AIR * atm.temperature));
        } else if (altitude < 32000.0f) {
            // Upper Stratosphere
            float lapse_rate = 0.001f;
            float T20 = 216.65f;
            float p20 = 5474.9f;
            atm.temperature = T20 + lapse_rate * (altitude - 20000.0f);
            atm.pressure = p20 * powf(atm.temperature / T20, -G / (lapse_rate * R_AIR));
        } else if (altitude < 47000.0f) {
            atm.temperature = 228.65f;
            atm.pressure = 868.0f * expf(-G * (altitude - 32000.0f) / (R_AIR * atm.temperature));
        } else if (altitude < 51000.0f) {
            float lapse_rate = 0.0028f;
            float T47 = 270.65f;
            float p47 = 110.9f;
            atm.temperature = T47 + lapse_rate * (altitude - 47000.0f);
            atm.pressure = p47 * powf(atm.temperature / T47, -G / (lapse_rate * R_AIR));
        } else if (altitude < 71000.0f) {
            atm.temperature = 270.65f;
            atm.pressure = 66.9f * expf(-G * (altitude - 51000.0f) / (R_AIR * atm.temperature));
        } else if (altitude < 84852.0f) {
            float lapse_rate = -0.0028f;
            float T71 = 214.65f;
            float p71 = 3.96f;
            atm.temperature = T71 + lapse_rate * (altitude - 71000.0f);
            atm.pressure = p71 * powf(atm.temperature / T71, -G / (lapse_rate * R_AIR));
        } else {
            // Very high altitude (exponential decay)
            atm.temperature = 186.87f;
            atm.pressure = 0.3734f * expf(-G * (altitude - 84852.0f) / (R_AIR * atm.temperature));
        }

        atm.density = atm.pressure / (R_AIR * atm.temperature);
        atm.speedOfSound = sqrtf(GAMMA * R_AIR * atm.temperature);

        // Mean free path (for rarefaction regime detection)
        atm.meanFreePath = (atm.density > 1e-6f) ? (1.458e-6f * powf(atm.temperature, 1.5f) / (atm.temperature + 110.4f)) / (atm.density * atm.speedOfSound) : 1.0f;

        return atm;
    }
};

// ============================================================================
// HYPERSONIC VEHICLE
// ============================================================================

class HypersonicVehicle {
private:
    // Vehicle parameters
    float mass;             // kg
    float referenceArea;    // m² (frontal area)
    float noseRadius;       // m
    float length;           // m
    float emissivity;       // TPS emissivity (0-1)
    float heatCapacity;     // J/(kg·K) - TPS thermal mass

    // State variables
    float x, y;             // Position (m)
    float vx, vy;           // Velocity (m/s)
    float altitude;         // m
    float flightPathAngle;  // radians
    float angleOfAttack;    // radians

    // Thermal state
    float surfaceTemp;      // K (TPS outer surface)
    float internalTemp;     // K (structure)
    float totalHeatLoad;    // MJ

    // History for visualization
    std::vector<float> trail_x;
    std::vector<float> trail_y;
    std::vector<float> trail_temp;
    int maxTrailLength;

    // Statistics
    float maxMach;
    float maxHeatFlux;
    float maxTemp;
    float maxDynamicPressure;

public:
    HypersonicVehicle()
        : mass(15000.0f)        // X-15 style
        , referenceArea(18.6f)  // m²
        , noseRadius(0.3f)      // Sharp nose
        , length(15.0f)
        , emissivity(0.85f)     // High-temp ceramic
        , heatCapacity(800.0f)  // J/(kg·K)
        , x(0), y(0)
        , vx(2000.0f), vy(0)
        , altitude(50000.0f)
        , flightPathAngle(0)
        , angleOfAttack(0.1f)   // ~5.7 degrees
        , surfaceTemp(300.0f)
        , internalTemp(300.0f)
        , totalHeatLoad(0)
        , maxTrailLength(1000)
        , maxMach(0)
        , maxHeatFlux(0)
        , maxTemp(300.0f)
        , maxDynamicPressure(0)
    {}

    // ========================================================================
    // SETTERS
    // ========================================================================

    void setMass(float m) { mass = m; }
    void setReferenceArea(float area) { referenceArea = area; }
    void setAltitude(float alt) { altitude = alt; }
    void setVelocity(float v) {
        float heading = atan2f(vy, vx);
        vx = v * cosf(heading);
        vy = v * sinf(heading);
    }
    void setAngleOfAttack(float aoa) { angleOfAttack = aoa * 3.14159f / 180.0f; }
    void setFlightPathAngle(float fpa) { flightPathAngle = fpa * 3.14159f / 180.0f; }

    // ========================================================================
    // HYPERSONIC AERODYNAMICS
    // ========================================================================

    float calculateMachNumber() const {
        Atmosphere atm = Atmosphere::calculate(altitude);
        float velocity = sqrtf(vx * vx + vy * vy);
        return velocity / atm.speedOfSound;
    }

    float calculateDynamicPressure() const {
        Atmosphere atm = Atmosphere::calculate(altitude);
        float velocity = sqrtf(vx * vx + vy * vy);
        return 0.5f * atm.density * velocity * velocity;
    }

    float calculateReynoldsNumber() const {
        Atmosphere atm = Atmosphere::calculate(altitude);
        float velocity = sqrtf(vx * vx + vy * vy);
        float mu = 1.458e-6f * powf(atm.temperature, 1.5f) / (atm.temperature + 110.4f);
        return atm.density * velocity * length / mu;
    }

    float calculateKnudsenNumber() const {
        Atmosphere atm = Atmosphere::calculate(altitude);
        return atm.meanFreePath / length;
    }

    // Newtonian impact theory for hypersonic flow
    float calculateCd() const {
        float mach = calculateMachNumber();

        if (mach < 5.0f) {
            // Subsonic/supersonic regime
            return 0.2f + 0.05f * mach;
        } else {
            // Hypersonic - Newtonian theory
            float sinAlpha = sinf(angleOfAttack);
            float Cd_wave = 2.0f * sinAlpha * sinAlpha * sinAlpha;
            float Cd_base = 0.12f; // Base drag
            float Cd_friction = 0.02f / sqrtf(calculateReynoldsNumber() + 1e-6f);
            return Cd_wave + Cd_base + Cd_friction;
        }
    }

    float calculateCl() const {
        float mach = calculateMachNumber();

        if (mach < 5.0f) {
            return 5.0f * angleOfAttack; // Linear subsonic
        } else {
            // Hypersonic - Newtonian theory
            float sinAlpha = sinf(angleOfAttack);
            float cosAlpha = cosf(angleOfAttack);
            return 2.0f * sinAlpha * sinAlpha * cosAlpha;
        }
    }

    // ========================================================================
    // AEROTHERMAL HEATING
    // ========================================================================

    float calculateStagnationTemperature() const {
        Atmosphere atm = Atmosphere::calculate(altitude);
        float velocity = sqrtf(vx * vx + vy * vy);
        float mach = velocity / atm.speedOfSound;

        // Total temperature behind normal shock
        float gamma = (mach > 10.0f) ? GAMMA_HOT : GAMMA; // Real gas effects
        float T_total = atm.temperature * (1.0f + 0.5f * (gamma - 1.0f) * mach * mach);

        return T_total;
    }

    float calculateConvectiveHeatFlux() const {
        // Fay-Riddell correlation for stagnation point heating
        Atmosphere atm = Atmosphere::calculate(altitude);
        float velocity = sqrtf(vx * vx + vy * vy);

        float rho_ratio = atm.density / RHO0_ISA;
        float v_ratio = velocity / 7500.0f; // Normalize to orbital velocity

        // Simplified Fay-Riddell: q ∝ √(ρ/R_n) * V³
        float q_conv = 1.83e-4f * sqrtf(atm.density / (noseRadius + 0.01f)) * velocity * velocity * velocity;

        return q_conv; // W/m²
    }

    float calculateRadiativeHeatFlux() const {
        // Stefan-Boltzmann radiation from hot shock layer
        float T_stag = calculateStagnationTemperature();

        // Radiative heating becomes significant above Mach 10
        float mach = calculateMachNumber();
        float radiativeFactor = (mach > 10.0f) ? (mach - 10.0f) / 15.0f : 0.0f;
        radiativeFactor = std::min(radiativeFactor, 1.0f);

        float q_rad = STEFAN_BOLTZMANN * radiativeFactor * powf(T_stag, 4.0f) * 0.1f;

        return q_rad; // W/m²
    }

    float calculateThermalProtection() const {
        // Radiative cooling from TPS
        float q_reradiation = emissivity * STEFAN_BOLTZMANN * powf(surfaceTemp, 4.0f);
        return q_reradiation; // W/m²
    }

    // ========================================================================
    // PHYSICS UPDATE
    // ========================================================================

    void update(float dt) {
        // Get atmospheric conditions
        Atmosphere atm = Atmosphere::calculate(altitude);

        // Calculate airspeed
        float velocity = sqrtf(vx * vx + vy * vy);
        float dynamicPressure = calculateDynamicPressure();

        // Aerodynamic coefficients
        float Cl = calculateCl();
        float Cd = calculateCd();

        // Forces
        float lift = Cl * dynamicPressure * referenceArea;
        float drag = Cd * dynamicPressure * referenceArea;
        float weight = mass * G;

        // Flight path angle
        float gamma = atan2f(vy, vx);

        // Convert to inertial frame
        float Fx = -drag * cosf(gamma) + lift * sinf(gamma);
        float Fy = -drag * sinf(gamma) - lift * cosf(gamma) - weight;

        // Integrate equations of motion
        float ax = Fx / mass;
        float ay = Fy / mass;

        vx += ax * dt;
        vy += ay * dt;

        x += vx * dt;
        y += vy * dt;
        altitude += vy * dt;

        // Heating calculations
        float q_conv = calculateConvectiveHeatFlux();
        float q_rad = calculateRadiativeHeatFlux();
        float q_out = calculateThermalProtection();
        float q_net = q_conv + q_rad - q_out;

        // Surface temperature update (simplified 1D heat transfer)
        float thermalMass = mass * 0.1f * heatCapacity; // 10% of mass is TPS
        surfaceTemp += (q_net * referenceArea / thermalMass) * dt;
        surfaceTemp = std::max(300.0f, std::min(surfaceTemp, 2500.0f)); // Limit to material constraints

        // Internal heating (conduction through TPS)
        float k_tps = 0.5f; // W/(m·K) - insulative
        float thickness = 0.05f; // m
        float q_internal = k_tps * (surfaceTemp - internalTemp) / thickness;
        internalTemp += (q_internal * referenceArea / (mass * heatCapacity * 0.9f)) * dt;

        // Total heat load integration
        totalHeatLoad += (q_conv + q_rad) * referenceArea * dt / 1e6f; // Convert to MJ

        // Ground collision
        if (altitude < 0) {
            altitude = 0;
            vy = 0;
            vx *= 0.8f; // Friction
        }

        // Update trail
        trail_x.push_back(x);
        trail_y.push_back(altitude);
        trail_temp.push_back(surfaceTemp);
        if (trail_x.size() > maxTrailLength) {
            trail_x.erase(trail_x.begin());
            trail_y.erase(trail_y.begin());
            trail_temp.erase(trail_temp.begin());
        }

        // Update statistics
        float mach = calculateMachNumber();
        if (mach > maxMach) maxMach = mach;

        float heatFlux = q_conv + q_rad;
        if (heatFlux > maxHeatFlux) maxHeatFlux = heatFlux;

        if (surfaceTemp > maxTemp) maxTemp = surfaceTemp;

        if (dynamicPressure > maxDynamicPressure) maxDynamicPressure = dynamicPressure;
    }

    // ========================================================================
    // GETTERS
    // ========================================================================

    float getX() const { return x; }
    float getY() const { return altitude; }
    float getAltitude() const { return altitude; }
    float getVelocity() const { return sqrtf(vx * vx + vy * vy); }
    float getVerticalSpeed() const { return vy; }
    float getMachNumber() const { return calculateMachNumber(); }
    float getDynamicPressure() const { return calculateDynamicPressure() / 1000.0f; } // kPa
    float getReynoldsNumber() const { return calculateReynoldsNumber(); }
    float getKnudsenNumber() const { return calculateKnudsenNumber(); }

    float getSurfaceTemperature() const { return surfaceTemp - 273.15f; } // Celsius
    float getInternalTemperature() const { return internalTemp - 273.15f; } // Celsius
    float getStagnationTemperature() const { return calculateStagnationTemperature() - 273.15f; } // Celsius
    float getConvectiveHeatFlux() const { return calculateConvectiveHeatFlux() / 1000.0f; } // kW/m²
    float getRadiativeHeatFlux() const { return calculateRadiativeHeatFlux() / 1000.0f; } // kW/m²
    float getTotalHeatLoad() const { return totalHeatLoad; } // MJ

    float getAtmosphericDensity() const {
        return Atmosphere::calculate(altitude).density;
    }
    float getAtmosphericPressure() const {
        return Atmosphere::calculate(altitude).pressure;
    }
    float getAtmosphericTemperature() const {
        return Atmosphere::calculate(altitude).temperature - 273.15f; // Celsius
    }

    float getMaxMach() const { return maxMach; }
    float getMaxHeatFlux() const { return maxHeatFlux / 1000.0f; } // kW/m²
    float getMaxTemperature() const { return maxTemp - 273.15f; } // Celsius
    float getMaxDynamicPressure() const { return maxDynamicPressure / 1000.0f; } // kPa

    val getTrail() const {
        val trail = val::array();
        for (size_t i = 0; i < trail_x.size(); ++i) {
            trail.call<void>("push", trail_x[i]);
            trail.call<void>("push", trail_y[i]);
        }
        return trail;
    }

    val getTrailTemperature() const {
        val temps = val::array();
        for (size_t i = 0; i < trail_temp.size(); ++i) {
            temps.call<void>("push", trail_temp[i] - 273.15f);
        }
        return temps;
    }

    void reset() {
        x = 0;
        y = 0;
        altitude = 50000.0f;
        vx = 2000.0f;
        vy = 0;
        angleOfAttack = 0.1f;
        flightPathAngle = 0;
        surfaceTemp = 300.0f;
        internalTemp = 300.0f;
        totalHeatLoad = 0;
        trail_x.clear();
        trail_y.clear();
        trail_temp.clear();
        maxMach = 0;
        maxHeatFlux = 0;
        maxTemp = 300.0f;
        maxDynamicPressure = 0;
    }
};

// ============================================================================
// EMSCRIPTEN BINDINGS
// ============================================================================

EMSCRIPTEN_BINDINGS(hypersonic_module) {
    class_<HypersonicVehicle>("HypersonicVehicle")
        .constructor<>()
        .function("setMass", &HypersonicVehicle::setMass)
        .function("setReferenceArea", &HypersonicVehicle::setReferenceArea)
        .function("setAltitude", &HypersonicVehicle::setAltitude)
        .function("setVelocity", &HypersonicVehicle::setVelocity)
        .function("setAngleOfAttack", &HypersonicVehicle::setAngleOfAttack)
        .function("setFlightPathAngle", &HypersonicVehicle::setFlightPathAngle)
        .function("update", &HypersonicVehicle::update)
        .function("getX", &HypersonicVehicle::getX)
        .function("getY", &HypersonicVehicle::getY)
        .function("getAltitude", &HypersonicVehicle::getAltitude)
        .function("getVelocity", &HypersonicVehicle::getVelocity)
        .function("getVerticalSpeed", &HypersonicVehicle::getVerticalSpeed)
        .function("getMachNumber", &HypersonicVehicle::getMachNumber)
        .function("getDynamicPressure", &HypersonicVehicle::getDynamicPressure)
        .function("getReynoldsNumber", &HypersonicVehicle::getReynoldsNumber)
        .function("getKnudsenNumber", &HypersonicVehicle::getKnudsenNumber)
        .function("getSurfaceTemperature", &HypersonicVehicle::getSurfaceTemperature)
        .function("getInternalTemperature", &HypersonicVehicle::getInternalTemperature)
        .function("getStagnationTemperature", &HypersonicVehicle::getStagnationTemperature)
        .function("getConvectiveHeatFlux", &HypersonicVehicle::getConvectiveHeatFlux)
        .function("getRadiativeHeatFlux", &HypersonicVehicle::getRadiativeHeatFlux)
        .function("getTotalHeatLoad", &HypersonicVehicle::getTotalHeatLoad)
        .function("getAtmosphericDensity", &HypersonicVehicle::getAtmosphericDensity)
        .function("getAtmosphericPressure", &HypersonicVehicle::getAtmosphericPressure)
        .function("getAtmosphericTemperature", &HypersonicVehicle::getAtmosphericTemperature)
        .function("getMaxMach", &HypersonicVehicle::getMaxMach)
        .function("getMaxHeatFlux", &HypersonicVehicle::getMaxHeatFlux)
        .function("getMaxTemperature", &HypersonicVehicle::getMaxTemperature)
        .function("getMaxDynamicPressure", &HypersonicVehicle::getMaxDynamicPressure)
        .function("getTrail", &HypersonicVehicle::getTrail)
        .function("getTrailTemperature", &HypersonicVehicle::getTrailTemperature)
        .function("reset", &HypersonicVehicle::reset);
}
