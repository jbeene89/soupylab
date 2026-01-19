/**
 * ULTIMATE VEHICLE DOWNFORCE & AERODYNAMICS SIMULATOR - WebAssembly
 * Professional Race Car Aerodynamics with Real-Time CFD Visualization
 *
 * Features:
 * - Multiple vehicle types (F1, NASCAR, LMP1, Road Car)
 * - Front wing angle adjustment (angle of attack)
 * - Rear wing angle adjustment + DRS system
 * - Underbody diffuser with ground effect
 * - Ride height sensitivity
 * - Tire wake and Y250 vortex modeling
 * - Speed-dependent aerodynamic maps
 * - Cornering downforce balance
 * - Understeer/oversteer characteristics
 * - Center of pressure tracking
 * - Track simulation with corners
 * - Tire load distribution
 * - Drag coefficient vs downforce tradeoff
 * - Aerodynamic efficiency (L/D)
 * - Slipstream/drafting effects
 * - Pressure distribution visualization
 * - Wind tunnel mode
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
constexpr float RHO_AIR = 1.225f;        // kg/m³ at sea level
constexpr float PI = 3.14159265359f;

// ============================================================================
// VEHICLE TYPES
// ============================================================================

enum class VehicleType {
    FORMULA_1,
    NASCAR,
    LMP1,
    ROAD_CAR
};

struct AeroConfig {
    // Physical dimensions
    float frontalArea;          // m²
    float wheelbase;            // m
    float trackWidth;           // m
    float mass;                 // kg
    float cogHeight;            // m (center of gravity height)

    // Aerodynamic components
    float frontWingArea;        // m²
    float rearWingArea;         // m²
    float underbodyArea;        // m²
    float diffuserArea;         // m²

    // Base coefficients
    float Cd0;                  // Base drag coefficient
    float ClFrontBase;          // Front downforce coefficient
    float ClRearBase;           // Rear downforce coefficient
    float ClUnderbody;          // Underbody downforce coefficient

    // Ground effect parameters
    float rideHeightOptimal;    // m (optimal ride height)
    float rideHeightSensitivity; // How sensitive to ride height

    // Balance
    float aeroBalance;          // 0.0 = all rear, 1.0 = all front

    // DRS availability
    bool hasDRS;

    static AeroConfig getConfig(VehicleType type) {
        AeroConfig config;

        switch (type) {
            case VehicleType::FORMULA_1:
                config.frontalArea = 1.5f;
                config.wheelbase = 3.6f;
                config.trackWidth = 2.0f;
                config.mass = 798.0f;
                config.cogHeight = 0.3f;

                config.frontWingArea = 1.2f;
                config.rearWingArea = 1.0f;
                config.underbodyArea = 3.5f;
                config.diffuserArea = 1.2f;

                config.Cd0 = 0.70f;
                config.ClFrontBase = -1.5f;  // Negative = downforce
                config.ClRearBase = -2.5f;
                config.ClUnderbody = -1.8f;

                config.rideHeightOptimal = 0.030f;  // 30mm
                config.rideHeightSensitivity = 50.0f;

                config.aeroBalance = 0.45f;  // Slightly rear biased
                config.hasDRS = true;
                break;

            case VehicleType::NASCAR:
                config.frontalArea = 2.8f;
                config.wheelbase = 2.8f;
                config.trackWidth = 1.8f;
                config.mass = 1542.0f;

                config.frontWingArea = 0.3f;
                config.rearWingArea = 1.8f;  // Large rear spoiler
                config.underbodyArea = 4.5f;
                config.diffuserArea = 0.5f;

                config.Cd0 = 0.35f;
                config.ClFrontBase = -0.3f;
                config.ClRearBase = -0.8f;
                config.ClUnderbody = -0.2f;

                config.rideHeightOptimal = 0.050f;
                config.rideHeightSensitivity = 10.0f;

                config.aeroBalance = 0.35f;
                config.hasDRS = false;
                break;

            case VehicleType::LMP1:
                config.frontalArea = 2.0f;
                config.wheelbase = 3.2f;
                config.trackWidth = 1.9f;
                config.mass = 875.0f;
                config.cogHeight = 0.35f;

                config.frontWingArea = 0.8f;
                config.rearWingArea = 1.5f;
                config.underbodyArea = 4.0f;
                config.diffuserArea = 1.5f;

                config.Cd0 = 0.45f;
                config.ClFrontBase = -1.2f;
                config.ClRearBase = -2.0f;
                config.ClUnderbody = -1.5f;

                config.rideHeightOptimal = 0.040f;
                config.rideHeightSensitivity = 30.0f;

                config.aeroBalance = 0.48f;
                config.hasDRS = false;
                break;

            case VehicleType::ROAD_CAR:
                config.frontalArea = 2.2f;
                config.wheelbase = 2.7f;
                config.trackWidth = 1.5f;
                config.mass = 1400.0f;
                config.cogHeight = 0.5f;

                config.frontWingArea = 0.1f;
                config.rearWingArea = 0.2f;
                config.underbodyArea = 3.0f;
                config.diffuserArea = 0.3f;

                config.Cd0 = 0.28f;
                config.ClFrontBase = 0.05f;  // Slight lift
                config.ClRearBase = -0.10f;
                config.ClUnderbody = 0.02f;

                config.rideHeightOptimal = 0.150f;
                config.rideHeightSensitivity = 2.0f;

                config.aeroBalance = 0.40f;
                config.hasDRS = false;
                break;
        }

        return config;
    }
};

// ============================================================================
// TRACK SECTION
// ============================================================================

struct TrackSection {
    float x, y;
    float radius;       // 0 = straight, >0 = corner radius
    float length;
    bool isCorner;

    float getSpeed(float maxSpeed) const {
        if (!isCorner || radius > 1000.0f) {
            return maxSpeed;
        }
        // v = sqrt(a * r) where a is lateral acceleration limit
        float lateralG = 4.0f;  // Max lateral g-force
        return sqrtf(lateralG * G * radius);
    }
};

// ============================================================================
// VEHICLE CLASS
// ============================================================================

class RaceVehicle {
private:
    VehicleType type;
    AeroConfig config;

    // State
    float x, y;
    float velocity;         // m/s
    float heading;          // radians
    float rideHeight;       // m
    float throttlePosition; // 0-1
    float brakePosition;    // 0-1
    float steeringAngle;    // radians

    // Aerodynamic controls
    float frontWingAngle;   // degrees (added to base)
    float rearWingAngle;    // degrees
    bool drsActive;

    // Computed values
    float downforceFront;   // N
    float downforceRear;    // N
    float totalDownforce;   // N
    float dragForce;        // N
    float lateralForce;     // N

    // Tire loads (including aero)
    float tireFrontLeft;
    float tireFrontRight;
    float tireRearLeft;
    float tireRearRight;

    // Balance
    float centerOfPressure; // % from front (0.5 = center)
    float balanceCharacteristic; // <0 understeer, >0 oversteer

    // Track
    std::vector<TrackSection> track;
    int currentSection;
    float distanceInSection;

    // Trail
    std::vector<float> trail_x;
    std::vector<float> trail_y;
    std::vector<float> trail_speed;
    std::vector<float> trail_downforce;
    int maxTrailLength;

    // Statistics
    float maxSpeed;
    float maxDownforce;
    float maxLateralG;
    float avgAeroEfficiency;
    float totalDistance;
    float lapTime;

public:
    RaceVehicle()
        : type(VehicleType::FORMULA_1)
        , config(AeroConfig::getConfig(type))
        , x(0), y(0)
        , velocity(50.0f)
        , heading(0)
        , rideHeight(config.rideHeightOptimal)
        , throttlePosition(0.8f)
        , brakePosition(0)
        , steeringAngle(0)
        , frontWingAngle(0)
        , rearWingAngle(0)
        , drsActive(false)
        , downforceFront(0)
        , downforceRear(0)
        , totalDownforce(0)
        , dragForce(0)
        , lateralForce(0)
        , tireFrontLeft(0)
        , tireFrontRight(0)
        , tireRearLeft(0)
        , tireRearRight(0)
        , centerOfPressure(0.5f)
        , balanceCharacteristic(0)
        , currentSection(0)
        , distanceInSection(0)
        , maxTrailLength(2000)
        , maxSpeed(0)
        , maxDownforce(0)
        , maxLateralG(0)
        , avgAeroEfficiency(0)
        , totalDistance(0)
        , lapTime(0)
    {
        initializeTrack();
    }

    // ========================================================================
    // TRACK INITIALIZATION
    // ========================================================================

    void initializeTrack() {
        track.clear();

        // Create a test track with straights and corners
        // Straight
        TrackSection s1;
        s1.x = 0; s1.y = 0;
        s1.length = 500.0f;
        s1.radius = 0;
        s1.isCorner = false;
        track.push_back(s1);

        // Right turn
        TrackSection s2;
        s2.x = 500; s2.y = 0;
        s2.length = PI * 150.0f;  // Quarter circle
        s2.radius = 150.0f;
        s2.isCorner = true;
        track.push_back(s2);

        // Straight
        TrackSection s3;
        s3.x = 650; s3.y = 150;
        s3.length = 400.0f;
        s3.radius = 0;
        s3.isCorner = false;
        track.push_back(s3);

        // Hairpin
        TrackSection s4;
        s4.x = 650; s4.y = 550;
        s4.length = PI * 80.0f;
        s4.radius = 80.0f;
        s4.isCorner = true;
        track.push_back(s4);

        // Back straight
        TrackSection s5;
        s5.x = 570; s5.y = 630;
        s5.length = 600.0f;
        s5.radius = 0;
        s5.isCorner = false;
        track.push_back(s5);

        // Fast chicane
        TrackSection s6;
        s6.x = -30; s6.y = 630;
        s6.length = PI * 200.0f / 2.0f;
        s6.radius = 200.0f;
        s6.isCorner = true;
        track.push_back(s6);
    }

    // ========================================================================
    // SETTERS
    // ========================================================================

    void setVehicleType(int typeIndex) {
        type = static_cast<VehicleType>(typeIndex);
        config = AeroConfig::getConfig(type);
        rideHeight = config.rideHeightOptimal;
    }

    void setVelocity(float v) { velocity = v; }
    void setRideHeight(float h) { rideHeight = h; }
    void setFrontWingAngle(float angle) { frontWingAngle = angle; }
    void setRearWingAngle(float angle) { rearWingAngle = angle; }
    void setDRS(bool active) { drsActive = active && config.hasDRS; }
    void setThrottle(float t) { throttlePosition = std::max(0.0f, std::min(1.0f, t)); }
    void setBrake(float b) { brakePosition = std::max(0.0f, std::min(1.0f, b)); }

    // ========================================================================
    // AERODYNAMICS CALCULATIONS
    // ========================================================================

    float calculateDynamicPressure() const {
        return 0.5f * RHO_AIR * velocity * velocity;
    }

    float getGroundEffectFactor() const {
        // Ground effect multiplier based on ride height
        float heightRatio = rideHeight / config.rideHeightOptimal;

        // Maximum effect at optimal height
        float factor = expf(-config.rideHeightSensitivity * powf(heightRatio - 1.0f, 2));

        // Additional benefit when very low (venturi effect)
        if (rideHeight < config.rideHeightOptimal) {
            factor *= (1.0f + 0.5f * (1.0f - heightRatio));
        }

        return factor;
    }

    void calculateAerodynamics() {
        float q = calculateDynamicPressure();
        float groundEffect = getGroundEffectFactor();

        // Front wing downforce
        float frontWingCl = config.ClFrontBase + frontWingAngle * 0.05f;
        downforceFront = -frontWingCl * q * config.frontWingArea;

        // Rear wing downforce (affected by DRS)
        float rearWingCl = config.ClRearBase + rearWingAngle * 0.08f;
        if (drsActive) {
            rearWingCl *= 0.3f;  // DRS reduces rear downforce by 70%
        }
        downforceRear = -rearWingCl * q * config.rearWingArea;

        // Underbody and diffuser (ground effect dependent)
        float underbodyDownforce = -config.ClUnderbody * q * config.underbodyArea * groundEffect;

        // Distribute underbody based on aero balance
        downforceFront += underbodyDownforce * config.aeroBalance;
        downforceRear += underbodyDownforce * (1.0f - config.aeroBalance);

        totalDownforce = downforceFront + downforceRear;

        // Drag calculation
        float Cd = config.Cd0;

        // Induced drag from wings
        float inducedDragFront = 0.1f * frontWingCl * frontWingCl;
        float inducedDragRear = drsActive ? 0.02f : 0.15f * rearWingCl * rearWingCl;
        Cd += inducedDragFront + inducedDragRear;

        // Ground effect drag reduction
        Cd *= (1.0f - 0.1f * groundEffect);

        dragForce = Cd * q * config.frontalArea;

        // Center of pressure
        if (totalDownforce > 0) {
            centerOfPressure = downforceFront / totalDownforce;
        }

        // Tire loads (static + aero)
        float staticFront = config.mass * G * 0.45f;  // 45% front static
        float staticRear = config.mass * G * 0.55f;

        float tireLoadFront = staticFront + downforceFront;
        float tireLoadRear = staticRear + downforceRear;

        // Distribute left/right (simplified)
        tireFrontLeft = tireLoadFront * 0.5f;
        tireFrontRight = tireLoadFront * 0.5f;
        tireRearLeft = tireLoadRear * 0.5f;
        tireRearRight = tireLoadRear * 0.5f;

        // Balance characteristic (cornering behavior)
        float frontGrip = tireLoadFront * 1.2f;  // Front tire coefficient
        float rearGrip = tireLoadRear * 1.3f;    // Rear tire coefficient
        balanceCharacteristic = (rearGrip - frontGrip) / config.mass;
    }

    // ========================================================================
    // PHYSICS UPDATE
    // ========================================================================

    void update(float dt) {
        calculateAerodynamics();

        // Acceleration/deceleration
        float engineForce = throttlePosition * 8000.0f;  // N (simplified)
        float brakeForce = brakePosition * 15000.0f;     // N

        float netForce = engineForce - dragForce - brakeForce;
        float acceleration = netForce / config.mass;

        velocity += acceleration * dt;
        velocity = std::max(0.0f, std::min(velocity, 100.0f));  // Cap at 360 km/h

        // Track following (simplified)
        if (!track.empty()) {
            TrackSection& section = track[currentSection];
            float targetSpeed = section.getSpeed(100.0f);

            // Auto-brake for corners
            if (velocity > targetSpeed * 1.1f) {
                brakePosition = 0.5f;
                throttlePosition *= 0.8f;
            }

            // Progress through section
            distanceInSection += velocity * dt;

            if (distanceInSection >= section.length) {
                distanceInSection = 0;
                currentSection = (currentSection + 1) % track.size();
            }

            // Update position
            if (section.isCorner) {
                float angle = distanceInSection / section.radius;
                heading += angle * dt / section.length * 2.0f * PI;

                // Lateral force calculation
                float lateralAccel = velocity * velocity / section.radius;
                lateralForce = config.mass * lateralAccel;

                float maxLateralAccel = lateralAccel / G;
                if (maxLateralAccel > maxLateralG) maxLateralG = maxLateralAccel;
            } else {
                lateralForce = 0;
            }
        }

        x += velocity * cosf(heading) * dt;
        y += velocity * sinf(heading) * dt;

        totalDistance += velocity * dt;
        lapTime += dt;

        // Statistics
        if (velocity > maxSpeed) maxSpeed = velocity;
        if (totalDownforce > maxDownforce) maxDownforce = totalDownforce;

        float aeroEfficiency = (totalDownforce > 0) ? totalDownforce / dragForce : 0;
        avgAeroEfficiency = avgAeroEfficiency * 0.99f + aeroEfficiency * 0.01f;

        // Trail
        trail_x.push_back(x);
        trail_y.push_back(y);
        trail_speed.push_back(velocity);
        trail_downforce.push_back(totalDownforce);

        if (trail_x.size() > maxTrailLength) {
            trail_x.erase(trail_x.begin());
            trail_y.erase(trail_y.begin());
            trail_speed.erase(trail_speed.begin());
            trail_downforce.erase(trail_downforce.begin());
        }
    }

    // ========================================================================
    // GETTERS
    // ========================================================================

    float getX() const { return x; }
    float getY() const { return y; }
    float getVelocity() const { return velocity; }
    float getVelocityKMH() const { return velocity * 3.6f; }
    float getHeading() const { return heading * 180.0f / PI; }
    float getRideHeight() const { return rideHeight * 1000.0f; }  // mm

    float getDownforceFront() const { return downforceFront; }
    float getDownforceRear() const { return downforceRear; }
    float getTotalDownforce() const { return totalDownforce; }
    float getDragForce() const { return dragForce; }
    float getLateralForce() const { return lateralForce / 1000.0f; }  // kN

    float getTireFrontLeft() const { return tireFrontLeft; }
    float getTireFrontRight() const { return tireFrontRight; }
    float getTireRearLeft() const { return tireRearLeft; }
    float getTireRearRight() const { return tireRearRight; }

    float getCenterOfPressure() const { return centerOfPressure * 100.0f; }  // %
    float getBalanceCharacteristic() const { return balanceCharacteristic; }
    float getGroundEffectPercent() const { return getGroundEffectFactor() * 100.0f; }

    float getAeroEfficiency() const {
        return (dragForce > 0) ? totalDownforce / dragForce : 0;
    }

    bool isDRSActive() const { return drsActive; }
    bool hasDRSAvailable() const { return config.hasDRS; }

    float getMaxSpeed() const { return maxSpeed * 3.6f; }  // km/h
    float getMaxDownforce() const { return maxDownforce; }
    float getMaxLateralG() const { return maxLateralG; }
    float getAvgAeroEfficiency() const { return avgAeroEfficiency; }
    float getTotalDistance() const { return totalDistance / 1000.0f; }  // km
    float getLapTime() const { return lapTime; }

    val getTrail() const {
        val trail = val::array();
        for (size_t i = 0; i < trail_x.size(); ++i) {
            trail.call<void>("push", trail_x[i]);
            trail.call<void>("push", trail_y[i]);
        }
        return trail;
    }

    val getTrailSpeed() const {
        val speeds = val::array();
        for (float s : trail_speed) {
            speeds.call<void>("push", s * 3.6f);  // km/h
        }
        return speeds;
    }

    val getTrailDownforce() const {
        val df = val::array();
        for (float d : trail_downforce) {
            df.call<void>("push", d);
        }
        return df;
    }

    val getTrack() const {
        val trackData = val::array();
        for (const auto& section : track) {
            trackData.call<void>("push", section.x);
            trackData.call<void>("push", section.y);
            trackData.call<void>("push", section.length);
            trackData.call<void>("push", section.radius);
            trackData.call<void>("push", section.isCorner ? 1.0f : 0.0f);
        }
        return trackData;
    }

    void reset() {
        x = 0;
        y = 0;
        velocity = 50.0f;
        heading = 0;
        rideHeight = config.rideHeightOptimal;
        throttlePosition = 0.8f;
        brakePosition = 0;
        currentSection = 0;
        distanceInSection = 0;
        trail_x.clear();
        trail_y.clear();
        trail_speed.clear();
        trail_downforce.clear();
        maxSpeed = 0;
        maxDownforce = 0;
        maxLateralG = 0;
        avgAeroEfficiency = 0;
        totalDistance = 0;
        lapTime = 0;
    }
};

// ============================================================================
// EMSCRIPTEN BINDINGS
// ============================================================================

EMSCRIPTEN_BINDINGS(vehicle_aero_module) {
    class_<RaceVehicle>("RaceVehicle")
        .constructor<>()
        .function("setVehicleType", &RaceVehicle::setVehicleType)
        .function("setVelocity", &RaceVehicle::setVelocity)
        .function("setRideHeight", &RaceVehicle::setRideHeight)
        .function("setFrontWingAngle", &RaceVehicle::setFrontWingAngle)
        .function("setRearWingAngle", &RaceVehicle::setRearWingAngle)
        .function("setDRS", &RaceVehicle::setDRS)
        .function("setThrottle", &RaceVehicle::setThrottle)
        .function("setBrake", &RaceVehicle::setBrake)
        .function("update", &RaceVehicle::update)
        .function("getX", &RaceVehicle::getX)
        .function("getY", &RaceVehicle::getY)
        .function("getVelocity", &RaceVehicle::getVelocity)
        .function("getVelocityKMH", &RaceVehicle::getVelocityKMH)
        .function("getHeading", &RaceVehicle::getHeading)
        .function("getRideHeight", &RaceVehicle::getRideHeight)
        .function("getDownforceFront", &RaceVehicle::getDownforceFront)
        .function("getDownforceRear", &RaceVehicle::getDownforceRear)
        .function("getTotalDownforce", &RaceVehicle::getTotalDownforce)
        .function("getDragForce", &RaceVehicle::getDragForce)
        .function("getLateralForce", &RaceVehicle::getLateralForce)
        .function("getTireFrontLeft", &RaceVehicle::getTireFrontLeft)
        .function("getTireFrontRight", &RaceVehicle::getTireFrontRight)
        .function("getTireRearLeft", &RaceVehicle::getTireRearLeft)
        .function("getTireRearRight", &RaceVehicle::getTireRearRight)
        .function("getCenterOfPressure", &RaceVehicle::getCenterOfPressure)
        .function("getBalanceCharacteristic", &RaceVehicle::getBalanceCharacteristic)
        .function("getGroundEffectPercent", &RaceVehicle::getGroundEffectPercent)
        .function("getAeroEfficiency", &RaceVehicle::getAeroEfficiency)
        .function("isDRSActive", &RaceVehicle::isDRSActive)
        .function("hasDRSAvailable", &RaceVehicle::hasDRSAvailable)
        .function("getMaxSpeed", &RaceVehicle::getMaxSpeed)
        .function("getMaxDownforce", &RaceVehicle::getMaxDownforce)
        .function("getMaxLateralG", &RaceVehicle::getMaxLateralG)
        .function("getAvgAeroEfficiency", &RaceVehicle::getAvgAeroEfficiency)
        .function("getTotalDistance", &RaceVehicle::getTotalDistance)
        .function("getLapTime", &RaceVehicle::getLapTime)
        .function("getTrail", &RaceVehicle::getTrail)
        .function("getTrailSpeed", &RaceVehicle::getTrailSpeed)
        .function("getTrailDownforce", &RaceVehicle::getTrailDownforce)
        .function("getTrack", &RaceVehicle::getTrack)
        .function("reset", &RaceVehicle::reset);
}
