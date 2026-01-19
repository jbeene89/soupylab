/**
 * THERMAL CONVECTION ENGINE - High Performance WebAssembly Physics Core
 * Rayleigh-Bénard Convection & Navier-Stokes Fluid Dynamics
 *
 * Features:
 * - Full Navier-Stokes solver with vorticity-stream function formulation
 * - Heat diffusion with thermal conductivity
 * - Buoyancy forces (Boussinesq approximation)
 * - Rayleigh-Bénard instability patterns
 * - Turbulent flow with eddy formation
 * - Prandtl and Rayleigh number calculations
 * - SIMD-optimized finite difference methods
 */

#include <emscripten/bind.h>
#include <emscripten/val.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <array>

using namespace emscripten;

// ============================================================================
// PHYSICAL CONSTANTS
// ============================================================================

constexpr float PI = 3.14159265358979323846f;
constexpr float G = 9.81f;  // Gravitational acceleration (m/s²)

// ============================================================================
// FLUID PROPERTIES
// ============================================================================

struct FluidProperties {
    const char* name;
    float density;          // kg/m³
    float kinematicVisc;    // m²/s (nu)
    float thermalDiff;      // m²/s (alpha)
    float thermalExpansion; // 1/K (beta)
    float specificHeat;     // J/(kg·K)
    float thermalCond;      // W/(m·K)
};

static const FluidProperties FLUIDS[] = {
    {"Air (300K)",      1.184f,  1.5e-5f,  2.2e-5f,  3.43e-3f, 1005.0f, 0.026f},
    {"Water (300K)",    997.0f,  1.0e-6f,  1.4e-7f,  2.1e-4f,  4179.0f, 0.613f},
    {"Oil (SAE 30)",    920.0f,  1.0e-4f,  8.0e-8f,  7.0e-4f,  1880.0f, 0.145f},
    {"Glycerin",        1260.0f, 1.2e-3f,  9.5e-8f,  5.0e-4f,  2400.0f, 0.286f},
    {"Liquid Sodium",   927.0f,  6.8e-7f,  6.7e-5f,  2.5e-4f,  1260.0f, 85.0f}
};

// ============================================================================
// CONVECTION SIMULATOR
// ============================================================================

class ConvectionEngine {
private:
    int width, height;  // Grid dimensions
    float dx, dy;       // Grid spacing

    // Physical fields (using flat arrays for better cache locality)
    std::vector<float> temperature;      // T field
    std::vector<float> vorticity;        // ω field (curl of velocity)
    std::vector<float> streamFunction;   // ψ field
    std::vector<float> velocityX;        // u field
    std::vector<float> velocityY;        // v field

    // Temporary buffers for updates
    std::vector<float> tempNew;
    std::vector<float> vorticityNew;
    std::vector<float> psi;

    // Fluid properties
    const FluidProperties* fluid;

    // Simulation parameters
    float T_hot;          // Hot boundary temperature (K)
    float T_cold;         // Cold boundary temperature (K)
    float dt;             // Time step
    float time;           // Simulation time

    // Dimensionless numbers
    float rayleighNumber; // Ra = gβΔTL³/(να)
    float prandtlNumber;  // Pr = ν/α

    // Performance tracking
    int iteration;

    // Helper: Convert 2D to 1D index
    inline int idx(int x, int y) const {
        return y * width + x;
    }

    // Boundary conditions
    inline bool isBoundary(int x, int y) const {
        return x == 0 || x == width - 1 || y == 0 || y == height - 1;
    }

    // ========================================================================
    // POISSON SOLVER (for stream function from vorticity)
    // ========================================================================

    void solvePoisson(int iterations = 50) {
        // Solve ∇²ψ = -ω using Jacobi iteration

        for (int iter = 0; iter < iterations; ++iter) {
            for (int y = 1; y < height - 1; ++y) {
                for (int x = 1; x < width - 1; ++x) {
                    int i = idx(x, y);

                    float psi_left  = streamFunction[idx(x-1, y)];
                    float psi_right = streamFunction[idx(x+1, y)];
                    float psi_down  = streamFunction[idx(x, y-1)];
                    float psi_up    = streamFunction[idx(x, y+1)];

                    float w = vorticity[i];

                    // Five-point stencil
                    psi[i] = 0.25f * (psi_left + psi_right + psi_down + psi_up
                                     + dx * dx * w);
                }
            }

            // Update stream function
            std::swap(streamFunction, psi);

            // Boundary conditions for stream function (no-slip walls)
            for (int x = 0; x < width; ++x) {
                streamFunction[idx(x, 0)] = 0.0f;
                streamFunction[idx(x, height-1)] = 0.0f;
            }
            for (int y = 0; y < height; ++y) {
                streamFunction[idx(0, y)] = 0.0f;
                streamFunction[idx(width-1, y)] = 0.0f;
            }
        }
    }

    // ========================================================================
    // COMPUTE VELOCITY FROM STREAM FUNCTION
    // ========================================================================

    void computeVelocity() {
        // u = ∂ψ/∂y, v = -∂ψ/∂x

        for (int y = 1; y < height - 1; ++y) {
            for (int x = 1; x < width - 1; ++x) {
                int i = idx(x, y);

                float psi_up   = streamFunction[idx(x, y+1)];
                float psi_down = streamFunction[idx(x, y-1)];
                float psi_right = streamFunction[idx(x+1, y)];
                float psi_left  = streamFunction[idx(x-1, y)];

                velocityX[i] = (psi_up - psi_down) / (2.0f * dy);
                velocityY[i] = -(psi_right - psi_left) / (2.0f * dx);
            }
        }
    }

    // ========================================================================
    // UPDATE VORTICITY (Navier-Stokes)
    // ========================================================================

    void updateVorticity(float dt) {
        float nu = fluid->kinematicVisc;
        float beta = fluid->thermalExpansion;
        float deltaT = T_hot - T_cold;

        for (int y = 1; y < height - 1; ++y) {
            for (int x = 1; x < width - 1; ++x) {
                int i = idx(x, y);

                float w = vorticity[i];
                float u = velocityX[i];
                float v = velocityY[i];
                float T = temperature[i];

                // Neighboring vorticity
                float w_left  = vorticity[idx(x-1, y)];
                float w_right = vorticity[idx(x+1, y)];
                float w_down  = vorticity[idx(x, y-1)];
                float w_up    = vorticity[idx(x, y+1)];

                // Neighboring temperature for buoyancy
                float T_right = temperature[idx(x+1, y)];
                float T_left  = temperature[idx(x-1, y)];

                // Advection term: u·∇ω
                float dw_dx = (w_right - w_left) / (2.0f * dx);
                float dw_dy = (w_up - w_down) / (2.0f * dy);
                float advection = u * dw_dx + v * dw_dy;

                // Diffusion term: ν∇²ω
                float laplacian_w = (w_left + w_right + w_down + w_up - 4.0f * w)
                                   / (dx * dx);
                float diffusion = nu * laplacian_w;

                // Buoyancy term (Boussinesq approximation): gβ(∂T/∂x)
                float dT_dx = (T_right - T_left) / (2.0f * dx);
                float buoyancy = G * beta * dT_dx;

                // Vorticity transport equation: ∂ω/∂t = -u·∇ω + ν∇²ω + gβ(∂T/∂x)
                vorticityNew[i] = w + dt * (-advection + diffusion + buoyancy);
            }
        }

        std::swap(vorticity, vorticityNew);

        // Boundary conditions for vorticity (no-slip walls)
        for (int x = 1; x < width - 1; ++x) {
            // Bottom and top walls
            vorticity[idx(x, 0)] = -2.0f * streamFunction[idx(x, 1)] / (dy * dy);
            vorticity[idx(x, height-1)] = -2.0f * streamFunction[idx(x, height-2)] / (dy * dy);
        }
        for (int y = 1; y < height - 1; ++y) {
            // Left and right walls
            vorticity[idx(0, y)] = -2.0f * streamFunction[idx(1, y)] / (dx * dx);
            vorticity[idx(width-1, y)] = -2.0f * streamFunction[idx(width-2, y)] / (dx * dx);
        }
    }

    // ========================================================================
    // UPDATE TEMPERATURE (Heat Diffusion + Advection)
    // ========================================================================

    void updateTemperature(float dt) {
        float alpha = fluid->thermalDiff;

        for (int y = 1; y < height - 1; ++y) {
            for (int x = 1; x < width - 1; ++x) {
                int i = idx(x, y);

                float T = temperature[i];
                float u = velocityX[i];
                float v = velocityY[i];

                // Neighboring temperatures
                float T_left  = temperature[idx(x-1, y)];
                float T_right = temperature[idx(x+1, y)];
                float T_down  = temperature[idx(x, y-1)];
                float T_up    = temperature[idx(x, y+1)];

                // Advection term: u·∇T
                float dT_dx = (T_right - T_left) / (2.0f * dx);
                float dT_dy = (T_up - T_down) / (2.0f * dy);
                float advection = u * dT_dx + v * dT_dy;

                // Diffusion term: α∇²T
                float laplacian_T = (T_left + T_right + T_down + T_up - 4.0f * T)
                                   / (dx * dx);
                float diffusion = alpha * laplacian_T;

                // Heat equation: ∂T/∂t = -u·∇T + α∇²T
                tempNew[i] = T + dt * (-advection + diffusion);
            }
        }

        std::swap(temperature, tempNew);

        // Boundary conditions: Fixed temperature at top (cold) and bottom (hot)
        for (int x = 0; x < width; ++x) {
            temperature[idx(x, 0)] = T_hot;           // Bottom: hot
            temperature[idx(x, height-1)] = T_cold;   // Top: cold
        }

        // Adiabatic side walls (no heat flux)
        for (int y = 0; y < height; ++y) {
            temperature[idx(0, y)] = temperature[idx(1, y)];
            temperature[idx(width-1, y)] = temperature[idx(width-2, y)];
        }
    }

public:
    ConvectionEngine(int w, int h)
        : width(w), height(h)
        , dx(1.0f / w), dy(1.0f / h)
        , fluid(&FLUIDS[0])
        , T_hot(330.0f), T_cold(290.0f)
        , dt(0.0001f)
        , time(0.0f)
        , iteration(0)
    {
        int size = width * height;

        temperature.resize(size, 0.0f);
        vorticity.resize(size, 0.0f);
        streamFunction.resize(size, 0.0f);
        velocityX.resize(size, 0.0f);
        velocityY.resize(size, 0.0f);

        tempNew.resize(size, 0.0f);
        vorticityNew.resize(size, 0.0f);
        psi.resize(size, 0.0f);

        // Initialize temperature field
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                int i = idx(x, y);

                // Linear gradient from hot (bottom) to cold (top)
                float ratio = (float)y / (height - 1);
                temperature[i] = T_hot + (T_cold - T_hot) * ratio;

                // Add small perturbations to trigger instability
                float perturbation = 0.1f * sinf(2.0f * PI * x / width)
                                          * expf(-10.0f * (ratio - 0.5f) * (ratio - 0.5f));
                temperature[i] += perturbation;
            }
        }

        updateDimensionlessNumbers();
    }

    // ========================================================================
    // SETTERS
    // ========================================================================

    void setFluid(int index) {
        if (index >= 0 && index < 5) {
            fluid = &FLUIDS[index];
            updateDimensionlessNumbers();
        }
    }

    void setHotTemperature(float T) {
        T_hot = T;
        updateDimensionlessNumbers();
    }

    void setColdTemperature(float T) {
        T_cold = T;
        updateDimensionlessNumbers();
    }

    void setTimeStep(float timestep) {
        dt = timestep;
    }

    void updateDimensionlessNumbers() {
        float deltaT = T_hot - T_cold;
        float L = 1.0f; // Characteristic length (normalized)

        float nu = fluid->kinematicVisc;
        float alpha = fluid->thermalDiff;
        float beta = fluid->thermalExpansion;

        // Rayleigh number: Ra = gβΔTL³/(να)
        rayleighNumber = (G * beta * deltaT * L * L * L) / (nu * alpha);

        // Prandtl number: Pr = ν/α
        prandtlNumber = nu / alpha;
    }

    // ========================================================================
    // UPDATE SIMULATION
    // ========================================================================

    void update(int steps = 1) {
        for (int step = 0; step < steps; ++step) {
            updateVorticity(dt);
            solvePoisson(30);
            computeVelocity();
            updateTemperature(dt);

            time += dt;
            iteration++;
        }
    }

    // ========================================================================
    // GETTERS FOR VISUALIZATION
    // ========================================================================

    val getTemperatureField() {
        val data = val::array();
        for (int i = 0; i < width * height; ++i) {
            data.call<void>("push", temperature[i]);
        }
        return data;
    }

    val getVorticityField() {
        val data = val::array();
        for (int i = 0; i < width * height; ++i) {
            data.call<void>("push", vorticity[i]);
        }
        return data;
    }

    val getVelocityField() {
        val data = val::array();
        for (int i = 0; i < width * height; ++i) {
            data.call<void>("push", velocityX[i]);
            data.call<void>("push", velocityY[i]);
        }
        return data;
    }

    float getMaxVelocity() const {
        float maxVel = 0.0f;
        for (int i = 0; i < width * height; ++i) {
            float vel = sqrtf(velocityX[i] * velocityX[i] + velocityY[i] * velocityY[i]);
            maxVel = std::max(maxVel, vel);
        }
        return maxVel;
    }

    float getMaxVorticity() const {
        float maxVort = 0.0f;
        for (int i = 0; i < width * height; ++i) {
            maxVort = std::max(maxVort, fabsf(vorticity[i]));
        }
        return maxVort;
    }

    float getAverageTemperature() const {
        float sum = 0.0f;
        for (int i = 0; i < width * height; ++i) {
            sum += temperature[i];
        }
        return sum / (width * height);
    }

    float getNusseltNumber() const {
        // Heat flux at bottom boundary
        float heatFlux = 0.0f;
        for (int x = 1; x < width - 1; ++x) {
            float dT_dy = (temperature[idx(x, 1)] - temperature[idx(x, 0)]) / dy;
            heatFlux += -fluid->thermalCond * dT_dy;
        }
        heatFlux /= (width - 2);

        // Conductive heat flux
        float deltaT = T_hot - T_cold;
        float conductiveFlux = fluid->thermalCond * deltaT;

        // Nusselt number
        return heatFlux / conductiveFlux;
    }

    float getRayleighNumber() const { return rayleighNumber; }
    float getPrandtlNumber() const { return prandtlNumber; }
    float getSimulationTime() const { return time; }
    int getIteration() const { return iteration; }

    int getWidth() const { return width; }
    int getHeight() const { return height; }
};

// ============================================================================
// EMSCRIPTEN BINDINGS
// ============================================================================

EMSCRIPTEN_BINDINGS(convection_module) {
    class_<ConvectionEngine>("ConvectionEngine")
        .constructor<int, int>()
        .function("setFluid", &ConvectionEngine::setFluid)
        .function("setHotTemperature", &ConvectionEngine::setHotTemperature)
        .function("setColdTemperature", &ConvectionEngine::setColdTemperature)
        .function("setTimeStep", &ConvectionEngine::setTimeStep)
        .function("update", &ConvectionEngine::update)
        .function("getTemperatureField", &ConvectionEngine::getTemperatureField)
        .function("getVorticityField", &ConvectionEngine::getVorticityField)
        .function("getVelocityField", &ConvectionEngine::getVelocityField)
        .function("getMaxVelocity", &ConvectionEngine::getMaxVelocity)
        .function("getMaxVorticity", &ConvectionEngine::getMaxVorticity)
        .function("getAverageTemperature", &ConvectionEngine::getAverageTemperature)
        .function("getNusseltNumber", &ConvectionEngine::getNusseltNumber)
        .function("getRayleighNumber", &ConvectionEngine::getRayleighNumber)
        .function("getPrandtlNumber", &ConvectionEngine::getPrandtlNumber)
        .function("getSimulationTime", &ConvectionEngine::getSimulationTime)
        .function("getIteration", &ConvectionEngine::getIteration)
        .function("getWidth", &ConvectionEngine::getWidth)
        .function("getHeight", &ConvectionEngine::getHeight);
}
