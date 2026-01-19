/**
 * GROUND EFFECT VEHICLE SIMULATOR - WebAssembly
 * Wing-in-Ground Effect (WIG) Aerodynamics
 *
 * Features:
 * - Ground effect lift augmentation
 * - Reduced induced drag (RAM effect)
 * - Surface wave interaction
 * - Dynamic stability in ground effect
 * - Ekranoplan flight dynamics
 * - Cushion pressure distribution
 * - Power-augmented lift
 * - Transition to/from conventional flight
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
constexpr float RHO_AIR = 1.225f;       // Air density at sea level (kg/m³)
constexpr float RHO_WATER = 1025.0f;    // Seawater density (kg/m³)
constexpr float NU_AIR = 1.5e-5f;       // Kinematic viscosity of air (m²/s)

// ============================================================================
// GROUND EFFECT AERODYNAMICS
// ============================================================================

class GroundEffectCalculator {
public:
    // Ground effect factor (h/b where h = height, b = wingspan)
    static float calculateGroundEffect(float height, float wingspan) {
        float h_over_b = height / (wingspan + 0.1f);

        // Empirical ground effect factor (0 = on ground, 1 = out of ground effect)
        // Based on wind tunnel data for wings near ground
        float k = 1.0f - expf(-4.0f * h_over_b);

        return k;
    }

    // Lift augmentation factor
    static float liftAugmentation(float height, float wingspan, float aspectRatio) {
        float h_over_b = height / wingspan;

        // Wieselsberger-Tomotika formula modified
        float sigma = (16.0f * h_over_b * h_over_b) / (1.0f + 16.0f * h_over_b * h_over_b);

        // Lift coefficient multiplier
        float K_L = 1.0f + (aspectRatio / (aspectRatio + 2.0f)) * (1.0f / sigma - 1.0f);

        return K_L;
    }

    // Induced drag reduction factor
    static float dragReduction(float height, float wingspan) {
        float h_over_b = height / wingspan;

        // Ground effect dramatically reduces induced drag
        float phi = 16.0f * h_over_b * h_over_b / (1.0f + 16.0f * h_over_b * h_over_b);

        return phi; // 0 = max reduction, 1 = no reduction
    }

    // RAM pressure (cushion effect)
    static float ramPressure(float velocity, float height, float wingspan) {
        float dynamicPressure = 0.5f * RHO_AIR * velocity * velocity;

        // Empirical RAM pressure coefficient
        float h_over_b = height / wingspan;
        float C_ram = 0.3f * expf(-3.0f * h_over_b);

        return dynamicPressure * C_ram;
    }
};

// ============================================================================
// WAVE DYNAMICS (for water surface operation)
// ============================================================================

class WaveDynamics {
private:
    float waveHeight;
    float wavelength;
    float phase;

public:
    WaveDynamics()
        : waveHeight(0.5f)
        , wavelength(20.0f)
        , phase(0)
    {}

    void setWaveHeight(float h) { waveHeight = h; }
    void setWavelength(float lambda) { wavelength = lambda; }

    float getWaveHeight(float x, float time) const {
        float k = 2.0f * 3.14159f / wavelength;
        float omega = sqrtf(G * k); // Deep water dispersion
        return waveHeight * sinf(k * x - omega * time + phase);
    }

    void update(float dt) {
        float omega = sqrtf(G * 2.0f * 3.14159f / wavelength);
        phase += omega * dt;
    }
};

// ============================================================================
// EKRANOPLAN (GROUND EFFECT VEHICLE)
// ============================================================================

class Ekranoplan {
private:
    // Vehicle parameters
    float mass;             // kg
    float wingspan;         // m
    float chord;            // m (mean aerodynamic chord)
    float length;           // m
    float wingArea;         // m²
    float aspectRatio;      // wingspan²/area
    float Cl0;              // Zero-lift coefficient
    float ClAlpha;          // Lift curve slope
    float Cd0;              // Parasitic drag
    float thrustPower;      // W (engine power)

    // State variables
    float x, y;             // Position (m)
    float vx, vy;           // Velocity (m/s)
    float height;           // Height above surface (m)
    float pitch;            // Pitch angle (radians)
    float throttle;         // 0-1

    // Waves
    WaveDynamics waves;

    // History for visualization
    std::vector<float> trail_x;
    std::vector<float> trail_y;
    int maxTrailLength;

    // Statistics
    float maxSpeed;
    float minHeight;
    float totalDistance;
    float fuelUsed; // kg

public:
    Ekranoplan()
        : mass(380000.0f)       // KM Ekranoplan mass
        , wingspan(37.6f)
        , chord(8.0f)
        , length(92.0f)
        , wingArea(660.0f)
        , aspectRatio(2.14f)
        , Cl0(0.2f)
        , ClAlpha(4.5f)
        , Cd0(0.025f)
        , thrustPower(10000000.0f) // 10 MW total
        , x(0), y(0)
        , vx(50.0f), vy(0)
        , height(3.0f)
        , pitch(0.05f)
        , throttle(0.7f)
        , maxTrailLength(2000)
        , maxSpeed(0)
        , minHeight(3.0f)
        , totalDistance(0)
        , fuelUsed(0)
    {}

    // ========================================================================
    // SETTERS
    // ========================================================================

    void setMass(float m) { mass = m; }
    void setWingspan(float span) { wingspan = span; aspectRatio = span * span / wingArea; }
    void setHeight(float h) { height = h; }
    void setVelocity(float v) { vx = v; }
    void setPitch(float p) { pitch = p * 3.14159f / 180.0f; }
    void setThrottle(float t) { throttle = std::max(0.0f, std::min(1.0f, t)); }
    void setWaveHeight(float h) { waves.setWaveHeight(h); }
    void setWavelength(float lambda) { waves.setWavelength(lambda); }

    // ========================================================================
    // AERODYNAMIC CALCULATIONS
    // ========================================================================

    float calculateLift() const {
        float velocity = sqrtf(vx * vx + vy * vy);
        float dynamicPressure = 0.5f * RHO_AIR * velocity * velocity;

        // Angle of attack (pitch - flight path angle)
        float flightPathAngle = atan2f(vy, vx);
        float angleOfAttack = pitch - flightPathAngle;

        // Base lift coefficient
        float Cl = Cl0 + ClAlpha * angleOfAttack;

        // Ground effect augmentation
        float K_L = GroundEffectCalculator::liftAugmentation(height, wingspan, aspectRatio);
        Cl *= K_L;

        // RAM pressure contribution
        float q_ram = GroundEffectCalculator::ramPressure(velocity, height, wingspan);
        float liftRAM = q_ram * wingArea * 0.5f;

        float lift = Cl * dynamicPressure * wingArea + liftRAM;

        return lift;
    }

    float calculateDrag() const {
        float velocity = sqrtf(vx * vx + vy * vy);
        float dynamicPressure = 0.5f * RHO_AIR * velocity * velocity;

        float flightPathAngle = atan2f(vy, vx);
        float angleOfAttack = pitch - flightPathAngle;

        // Base lift coefficient
        float Cl = Cl0 + ClAlpha * angleOfAttack;
        float K_L = GroundEffectCalculator::liftAugmentation(height, wingspan, aspectRatio);
        Cl *= K_L;

        // Induced drag factor
        float K_induced = 1.0f / (3.14159f * aspectRatio);

        // Ground effect drag reduction
        float phi = GroundEffectCalculator::dragReduction(height, wingspan);

        // Total drag coefficient
        float Cd = Cd0 + K_induced * Cl * Cl * phi;

        float drag = Cd * dynamicPressure * wingArea;

        return drag;
    }

    float calculateThrust() const {
        float velocity = sqrtf(vx * vx + vy * vy) + 1.0f;
        float thrust = throttle * thrustPower / velocity;
        return thrust;
    }

    float getLiftToDrag() const {
        float lift = calculateLift();
        float drag = calculateDrag();
        return (drag > 0) ? lift / drag : 0;
    }

    // ========================================================================
    // PHYSICS UPDATE
    // ========================================================================

    void update(float dt) {
        // Update waves
        waves.update(dt);

        // Get surface height at current position
        float surfaceHeight = waves.getWaveHeight(x, 0);

        // Calculate forces
        float lift = calculateLift();
        float drag = calculateDrag();
        float thrust = calculateThrust();
        float weight = mass * G;

        // Flight path angle
        float gamma = atan2f(vy, vx);

        // Forces in inertial frame
        float Fx = thrust - drag * cosf(gamma) + lift * sinf(gamma);
        float Fy = -drag * sinf(gamma) - lift * cosf(gamma) - weight;

        // Accelerations
        float ax = Fx / mass;
        float ay = Fy / mass;

        // Integrate velocity
        vx += ax * dt;
        vy += ay * dt;

        // Integrate position
        x += vx * dt;
        y += vy * dt;

        // Update height above surface
        height = y - surfaceHeight;

        // Surface collision with wave interaction
        if (height < 0.5f) {
            // Touchdown
            height = 0.5f;
            y = surfaceHeight + height;
            vy = std::max(vy, 0.0f); // Can't go below surface
            vx *= 0.98f; // Water/ground friction
        }

        // Pitch dynamics (simplified)
        float pitchRate = 0;
        if (height < wingspan * 0.5f) {
            // In ground effect - stabilizing moment
            pitchRate = -0.1f * (pitch - 0.05f);
        } else {
            // Out of ground effect - less stable
            pitchRate = -0.05f * (pitch - 0.05f);
        }
        pitch += pitchRate * dt;

        // Fuel consumption (simplified)
        float fuelFlow = throttle * 0.5f; // kg/s at full throttle
        fuelUsed += fuelFlow * dt;

        // Update trail
        trail_x.push_back(x);
        trail_y.push_back(y);
        if (trail_x.size() > maxTrailLength) {
            trail_x.erase(trail_x.begin());
            trail_y.erase(trail_y.begin());
        }

        // Update statistics
        float velocity = sqrtf(vx * vx + vy * vy);
        if (velocity > maxSpeed) maxSpeed = velocity;
        if (height < minHeight) minHeight = height;
        totalDistance += velocity * dt;
    }

    // ========================================================================
    // GETTERS
    // ========================================================================

    float getX() const { return x; }
    float getY() const { return y; }
    float getHeight() const { return height; }
    float getVelocity() const { return sqrtf(vx * vx + vy * vy); }
    float getVerticalSpeed() const { return vy; }
    float getPitch() const { return pitch * 180.0f / 3.14159f; }
    float getThrottle() const { return throttle * 100.0f; }

    float getLift() const { return calculateLift() / 1000.0f; } // kN
    float getDrag() const { return calculateDrag() / 1000.0f; } // kN
    float getThrust() const { return calculateThrust() / 1000.0f; } // kN
    float getLiftToDragRatio() const { return getLiftToDrag(); }

    float getGroundEffectFactor() const {
        return 1.0f - GroundEffectCalculator::calculateGroundEffect(height, wingspan);
    }

    float getRAMPressure() const {
        float velocity = sqrtf(vx * vx + vy * vy);
        return GroundEffectCalculator::ramPressure(velocity, height, wingspan);
    }

    float getMaxSpeed() const { return maxSpeed; }
    float getMinHeight() const { return minHeight; }
    float getTotalDistance() const { return totalDistance / 1000.0f; } // km
    float getFuelUsed() const { return fuelUsed; } // kg

    val getTrail() const {
        val trail = val::array();
        for (size_t i = 0; i < trail_x.size(); ++i) {
            trail.call<void>("push", trail_x[i]);
            trail.call<void>("push", trail_y[i]);
        }
        return trail;
    }

    val getWaveProfile(int numPoints) const {
        val profile = val::array();
        float currentX = x;
        for (int i = 0; i < numPoints; ++i) {
            float xPos = currentX - 200.0f + (400.0f * i) / numPoints;
            float waveH = waves.getWaveHeight(xPos, 0);
            profile.call<void>("push", xPos);
            profile.call<void>("push", waveH);
        }
        return profile;
    }

    void reset() {
        x = 0;
        y = 3.0f;
        vx = 50.0f;
        vy = 0;
        height = 3.0f;
        pitch = 0.05f;
        throttle = 0.7f;
        trail_x.clear();
        trail_y.clear();
        maxSpeed = 0;
        minHeight = 3.0f;
        totalDistance = 0;
        fuelUsed = 0;
    }
};

// ============================================================================
// EMSCRIPTEN BINDINGS
// ============================================================================

EMSCRIPTEN_BINDINGS(ground_effect_module) {
    class_<Ekranoplan>("Ekranoplan")
        .constructor<>()
        .function("setMass", &Ekranoplan::setMass)
        .function("setWingspan", &Ekranoplan::setWingspan)
        .function("setHeight", &Ekranoplan::setHeight)
        .function("setVelocity", &Ekranoplan::setVelocity)
        .function("setPitch", &Ekranoplan::setPitch)
        .function("setThrottle", &Ekranoplan::setThrottle)
        .function("setWaveHeight", &Ekranoplan::setWaveHeight)
        .function("setWavelength", &Ekranoplan::setWavelength)
        .function("update", &Ekranoplan::update)
        .function("getX", &Ekranoplan::getX)
        .function("getY", &Ekranoplan::getY)
        .function("getHeight", &Ekranoplan::getHeight)
        .function("getVelocity", &Ekranoplan::getVelocity)
        .function("getVerticalSpeed", &Ekranoplan::getVerticalSpeed)
        .function("getPitch", &Ekranoplan::getPitch)
        .function("getThrottle", &Ekranoplan::getThrottle)
        .function("getLift", &Ekranoplan::getLift)
        .function("getDrag", &Ekranoplan::getDrag)
        .function("getThrust", &Ekranoplan::getThrust)
        .function("getLiftToDragRatio", &Ekranoplan::getLiftToDragRatio)
        .function("getGroundEffectFactor", &Ekranoplan::getGroundEffectFactor)
        .function("getRAMPressure", &Ekranoplan::getRAMPressure)
        .function("getMaxSpeed", &Ekranoplan::getMaxSpeed)
        .function("getMinHeight", &Ekranoplan::getMinHeight)
        .function("getTotalDistance", &Ekranoplan::getTotalDistance)
        .function("getFuelUsed", &Ekranoplan::getFuelUsed)
        .function("getTrail", &Ekranoplan::getTrail)
        .function("getWaveProfile", &Ekranoplan::getWaveProfile)
        .function("reset", &Ekranoplan::reset);
}
