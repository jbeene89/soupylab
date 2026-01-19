#include <emscripten/bind.h>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace emscripten;

const float PI = 3.14159265359f;
const float SPEED_OF_LIGHT = 299792458.0f; // m/s
const float PARSEC = 3.0857e16f; // meters
const float LIGHT_YEAR = 9.461e15f; // meters
const float STEFAN_BOLTZMANN = 5.67e-8f; // W/(m²·K⁴)

enum class GRBType {
    LONG_DURATION,    // >2s, from massive star collapse (hypernova)
    SHORT_DURATION,   // <2s, from neutron star merger
    ULTRA_LONG        // >1000s, rare collapsar events
};

struct ExtinctionLevel {
    float ozoneDepletion;        // Percentage
    float uvRadiationIncrease;   // Multiplier
    float surfaceRadiation;      // Additional radiation dose (Sv)
    float temperatureDrop;       // Global average °C
    float speciesExtinction;     // Percentage of species extinct
    bool totalExtinction;        // All complex life ends
    float recoveryTime;          // Years for biosphere recovery
};

class GammaRayBurst {
private:
    // GRB properties
    float distance;              // light-years
    float energy;                // ergs (10^51 - 10^54)
    float duration;              // seconds
    float beamingAngle;          // degrees (jet opening angle)
    GRBType type;

    // Propagation state
    float currentTime;           // seconds since GRB
    float wavefrOntDistance;     // light-years traveled
    bool hasReachedEarth;
    float impactTime;            // time when GRB hits Earth

    // Earth impact effects
    float photonFlux;            // photons/(cm²·s)
    float energyFlux;            // ergs/(cm²·s)
    float atmosphericIonization; // Percentage of atmosphere ionized

    // Biological effects
    ExtinctionLevel extinction;
    float ozoneRecoveryYears;

    // Real-time visualization
    std::vector<float> radiationSpectrum; // Energy distribution
    float currentOzoneLevel;     // Percentage remaining (100% = normal)

    void calculateExtinctionEffects() {
        // Calculate total energy received at Earth
        float distanceMeters = distance * LIGHT_YEAR;
        float solidAngle = 2.0f * PI * (1.0f - cosf(beamingAngle * PI / 180.0f));

        // Energy per unit area at Earth (assuming in beam)
        float fluence = energy / (4.0f * PI * distanceMeters * distanceMeters);

        // Ozone depletion from nitrogen oxide production
        // NOx produced by ionization, destroys ozone catalytically
        if (distance < 100.0f) {
            extinction.ozoneDepletion = 99.0f; // Near-complete destruction
            extinction.uvRadiationIncrease = 100.0f; // 100x UV
            extinction.surfaceRadiation = 10.0f; // 10 Sv (lethal)
            extinction.temperatureDrop = 5.0f;
            extinction.speciesExtinction = 95.0f;
            extinction.totalExtinction = true;
            extinction.recoveryTime = 1000000.0f; // Never recovers
        } else if (distance < 1000.0f) {
            extinction.ozoneDepletion = 90.0f - (distance / 1000.0f) * 40.0f;
            extinction.uvRadiationIncrease = 50.0f - (distance / 1000.0f) * 30.0f;
            extinction.surfaceRadiation = 5.0f - (distance / 1000.0f) * 4.0f;
            extinction.temperatureDrop = 3.0f - (distance / 1000.0f) * 2.0f;
            extinction.speciesExtinction = 80.0f - (distance / 1000.0f) * 50.0f;
            extinction.totalExtinction = (distance < 500.0f);
            extinction.recoveryTime = 10000.0f + (distance * 10.0f);
        } else if (distance < 6000.0f) {
            // Kill zone extends to about 6000 ly
            extinction.ozoneDepletion = 50.0f - ((distance - 1000.0f) / 5000.0f) * 40.0f;
            extinction.uvRadiationIncrease = 20.0f - ((distance - 1000.0f) / 5000.0f) * 15.0f;
            extinction.surfaceRadiation = 1.0f - ((distance - 1000.0f) / 5000.0f) * 0.9f;
            extinction.temperatureDrop = 1.0f - ((distance - 1000.0f) / 5000.0f) * 0.9f;
            extinction.speciesExtinction = 30.0f - ((distance - 1000.0f) / 5000.0f) * 25.0f;
            extinction.totalExtinction = false;
            extinction.recoveryTime = 5000.0f + (distance * 5.0f);
        } else {
            // Beyond kill zone - detectable but not deadly
            extinction.ozoneDepletion = 10.0f - ((distance - 6000.0f) / 10000.0f) * 10.0f;
            extinction.uvRadiationIncrease = 5.0f - ((distance - 6000.0f) / 10000.0f) * 4.0f;
            extinction.surfaceRadiation = 0.1f;
            extinction.temperatureDrop = 0.1f;
            extinction.speciesExtinction = 5.0f - ((distance - 6000.0f) / 10000.0f) * 5.0f;
            extinction.totalExtinction = false;
            extinction.recoveryTime = 1000.0f;
        }

        // Clamp values
        extinction.ozoneDepletion = std::max(0.0f, std::min(100.0f, extinction.ozoneDepletion));
        extinction.speciesExtinction = std::max(0.0f, std::min(100.0f, extinction.speciesExtinction));

        ozoneRecoveryYears = extinction.recoveryTime;
    }

    void calculateAtmosphericEffects() {
        if (!hasReachedEarth) {
            currentOzoneLevel = 100.0f;
            atmosphericIonization = 0.0f;
            return;
        }

        // Time since impact
        float timeSinceImpact = currentTime - impactTime;

        if (timeSinceImpact < duration) {
            // During GRB
            float progress = timeSinceImpact / duration;
            atmosphericIonization = progress * 80.0f; // Peak ionization
            currentOzoneLevel = 100.0f - (progress * extinction.ozoneDepletion);
        } else if (timeSinceImpact < duration + 3600.0f) {
            // First hour after - rapid ozone destruction
            float hourProgress = (timeSinceImpact - duration) / 3600.0f;
            atmosphericIonization = 80.0f * (1.0f - hourProgress);
            currentOzoneLevel = 100.0f - extinction.ozoneDepletion;
        } else {
            // Long-term recovery (years)
            float years = timeSinceImpact / (365.25f * 86400.0f);
            float recoveryProgress = years / ozoneRecoveryYears;
            currentOzoneLevel = 100.0f - extinction.ozoneDepletion * (1.0f - recoveryProgress);
            currentOzoneLevel = std::max(0.0f, std::min(100.0f, currentOzoneLevel));
            atmosphericIonization = 0.0f;
        }
    }

public:
    GammaRayBurst() {
        // Default: Ordovician mass extinction candidate
        distance = 6000.0f; // light-years
        energy = 1.0e52f; // ergs (10^52)
        duration = 10.0f; // seconds
        beamingAngle = 5.0f; // degrees
        type = GRBType::LONG_DURATION;

        reset();
    }

    void reset() {
        currentTime = 0.0f;
        wavefrontDistance = 0.0f;
        hasReachedEarth = false;
        impactTime = distance; // Time = distance in light-years (simplified)
        photonFlux = 0.0f;
        energyFlux = 0.0f;
        atmosphericIonization = 0.0f;
        currentOzoneLevel = 100.0f;

        extinction = ExtinctionLevel{};
        calculateExtinctionEffects();

        radiationSpectrum.clear();
        radiationSpectrum.resize(100, 0.0f);
    }

    void setDistance(float ly) {
        distance = std::max(10.0f, std::min(ly, 20000.0f));
        reset();
    }

    void setEnergy(float ergs) {
        energy = std::max(1e50f, std::min(ergs, 1e54f));
        reset();
    }

    void setDuration(float seconds) {
        duration = std::max(0.1f, std::min(seconds, 10000.0f));
        reset();
    }

    void setType(int t) {
        type = static_cast<GRBType>(t);
        switch(type) {
            case GRBType::LONG_DURATION:
                duration = 10.0f;
                energy = 1e52f;
                beamingAngle = 5.0f;
                break;
            case GRBType::SHORT_DURATION:
                duration = 0.5f;
                energy = 1e51f;
                beamingAngle = 10.0f;
                break;
            case GRBType::ULTRA_LONG:
                duration = 1000.0f;
                energy = 1e53f;
                beamingAngle = 3.0f;
                break;
        }
        reset();
    }

    // Preset scenarios
    void setOrdovician() {
        // Ordovician-Silurian mass extinction (443 Mya)
        distance = 6000.0f;
        energy = 1e52f;
        duration = 10.0f;
        type = GRBType::LONG_DURATION;
        reset();
    }

    void setNearby() {
        // Hypothetical nearby GRB
        distance = 1000.0f;
        energy = 1e52f;
        duration = 10.0f;
        type = GRBType::LONG_DURATION;
        reset();
    }

    void setClose() {
        // Very close GRB - total sterilization
        distance = 100.0f;
        energy = 1e53f;
        duration = 30.0f;
        type = GRBType::ULTRA_LONG;
        reset();
    }

    void setWR104() {
        // WR 104 - potential future GRB threat
        distance = 8000.0f;
        energy = 1e52f;
        duration = 20.0f;
        type = GRBType::LONG_DURATION;
        reset();
    }

    void setGalacticCore() {
        // GRB from galactic center
        distance = 26000.0f; // Distance to galactic center
        energy = 1e54f; // Extremely energetic
        duration = 100.0f;
        type = GRBType::ULTRA_LONG;
        reset();
    }

    void update(float dt) {
        currentTime += dt;

        // Wavefront propagates at speed of light
        wavefrontDistance = currentTime; // In light-years (simplified: 1 time unit = 1 ly)

        // Check if wavefront reached Earth
        if (!hasReachedEarth && wavefrontDistance >= distance) {
            hasReachedEarth = true;
            impactTime = currentTime;
        }

        // Calculate effects
        calculateAtmosphericEffects();

        // Update radiation spectrum (simplified)
        if (hasReachedEarth && (currentTime - impactTime) < duration) {
            float progress = (currentTime - impactTime) / duration;
            for (int i = 0; i < 100; i++) {
                float energy_keV = 1.0f + i * 10.0f; // 1 keV to 1 MeV
                // Peaked spectrum around 100-1000 keV
                float peak = 500.0f;
                float intensity = expf(-powf((energy_keV - peak) / 200.0f, 2.0f));
                radiationSpectrum[i] = intensity * (1.0f - progress);
            }
        } else {
            // Decay
            for (int i = 0; i < 100; i++) {
                radiationSpectrum[i] *= 0.99f;
            }
        }
    }

    // Getters
    float getDistance() const { return distance; }
    float getEnergy() const { return energy; }
    float getDuration() const { return duration; }
    float getCurrentTime() const { return currentTime; }
    float getWavefrontDistance() const { return wavefrontDistance; }
    bool reachedEarth() const { return hasReachedEarth; }
    float getTimeToImpact() const { return std::max(0.0f, distance - wavefrontDistance); }

    float getOzoneDepletion() const { return extinction.ozoneDepletion; }
    float getUVIncrease() const { return extinction.uvRadiationIncrease; }
    float getSurfaceRadiation() const { return extinction.surfaceRadiation; }
    float getTemperatureDrop() const { return extinction.temperatureDrop; }
    float getSpeciesExtinction() const { return extinction.speciesExtinction; }
    bool isTotalExtinction() const { return extinction.totalExtinction; }
    float getRecoveryTime() const { return extinction.recoveryTime; }

    float getCurrentOzoneLevel() const { return currentOzoneLevel; }
    float getAtmosphericIonization() const { return atmosphericIonization; }

    std::vector<float> getRadiationSpectrum() const { return radiationSpectrum; }
};

EMSCRIPTEN_BINDINGS(grb_module) {
    class_<GammaRayBurst>("GammaRayBurst")
        .constructor<>()
        .function("reset", &GammaRayBurst::reset)
        .function("update", &GammaRayBurst::update)
        .function("setDistance", &GammaRayBurst::setDistance)
        .function("setEnergy", &GammaRayBurst::setEnergy)
        .function("setDuration", &GammaRayBurst::setDuration)
        .function("setType", &GammaRayBurst::setType)
        .function("setOrdovician", &GammaRayBurst::setOrdovician)
        .function("setNearby", &GammaRayBurst::setNearby)
        .function("setClose", &GammaRayBurst::setClose)
        .function("setWR104", &GammaRayBurst::setWR104)
        .function("setGalacticCore", &GammaRayBurst::setGalacticCore)
        .function("getDistance", &GammaRayBurst::getDistance)
        .function("getEnergy", &GammaRayBurst::getEnergy)
        .function("getDuration", &GammaRayBurst::getDuration)
        .function("getCurrentTime", &GammaRayBurst::getCurrentTime)
        .function("getWavefrontDistance", &GammaRayBurst::getWavefrontDistance)
        .function("reachedEarth", &GammaRayBurst::reachedEarth)
        .function("getTimeToImpact", &GammaRayBurst::getTimeToImpact)
        .function("getOzoneDepletion", &GammaRayBurst::getOzoneDepletion)
        .function("getUVIncrease", &GammaRayBurst::getUVIncrease)
        .function("getSurfaceRadiation", &GammaRayBurst::getSurfaceRadiation)
        .function("getTemperatureDrop", &GammaRayBurst::getTemperatureDrop)
        .function("getSpeciesExtinction", &GammaRayBurst::getSpeciesExtinction)
        .function("isTotalExtinction", &GammaRayBurst::isTotalExtinction)
        .function("getRecoveryTime", &GammaRayBurst::getRecoveryTime)
        .function("getCurrentOzoneLevel", &GammaRayBurst::getCurrentOzoneLevel)
        .function("getAtmosphericIonization", &GammaRayBurst::getAtmosphericIonization)
        .function("getRadiationSpectrum", &GammaRayBurst::getRadiationSpectrum);

    register_vector<float>("VectorFloat");
}
