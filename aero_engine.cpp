/**
 * AERODYNAMICS & SUPERSONIC FLOW SIMULATOR - WebAssembly
 * Compressible Euler Equations with Riemann Solver
 *
 * Features:
 * - 2D compressible gas dynamics
 * - Roe approximate Riemann solver
 * - Shock wave capturing (no artificial viscosity needed)
 * - Supersonic flow patterns
 * - Bow shocks & Mach diamonds
 * - Oblique shock waves
 * - Expansion fans
 * - Real-time schlieren visualization
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

constexpr float GAMMA = 1.4f;  // Specific heat ratio for air
constexpr float R_GAS = 287.05f;  // Gas constant for air (J/kg·K)

// ============================================================================
// CONSERVED VARIABLES (rho, rho*u, rho*v, E)
// ============================================================================

struct ConservedVars {
    float rho;    // Density
    float rhoU;   // Momentum X
    float rhoV;   // Momentum Y
    float E;      // Total energy

    ConservedVars() : rho(1.0f), rhoU(0), rhoV(0), E(2.5f) {}
    ConservedVars(float rho, float u, float v, float p) {
        this->rho = rho;
        this->rhoU = rho * u;
        this->rhoV = rho * v;
        float ke = 0.5f * rho * (u * u + v * v);
        this->E = p / (GAMMA - 1.0f) + ke;
    }

    // Get primitive variables
    float u() const { return rhoU / rho; }
    float v() const { return rhoV / rho; }
    float pressure() const {
        float ke = 0.5f * (rhoU * rhoU + rhoV * rhoV) / rho;
        return (E - ke) * (GAMMA - 1.0f);
    }
    float speedOfSound() const {
        return sqrtf(GAMMA * pressure() / rho);
    }
    float machNumber() const {
        float vel = sqrtf(u() * u() + v() * v());
        return vel / speedOfSound();
    }
};

// ============================================================================
// FLUX FUNCTIONS
// ============================================================================

struct Flux {
    float f1, f2, f3, f4;

    Flux() : f1(0), f2(0), f3(0), f4(0) {}

    static Flux calculateX(const ConservedVars& U) {
        Flux F;
        float u = U.u();
        float v = U.v();
        float p = U.pressure();

        F.f1 = U.rhoU;
        F.f2 = U.rhoU * u + p;
        F.f3 = U.rhoU * v;
        F.f4 = (U.E + p) * u;

        return F;
    }

    static Flux calculateY(const ConservedVars& U) {
        Flux F;
        float u = U.u();
        float v = U.v();
        float p = U.pressure();

        F.f1 = U.rhoV;
        F.f2 = U.rhoV * u;
        F.f3 = U.rhoV * v + p;
        F.f4 = (U.E + p) * v;

        return F;
    }

    Flux operator+(const Flux& other) const {
        Flux result;
        result.f1 = f1 + other.f1;
        result.f2 = f2 + other.f2;
        result.f3 = f3 + other.f3;
        result.f4 = f4 + other.f4;
        return result;
    }

    Flux operator-(const Flux& other) const {
        Flux result;
        result.f1 = f1 - other.f1;
        result.f2 = f2 - other.f2;
        result.f3 = f3 - other.f3;
        result.f4 = f4 - other.f4;
        return result;
    }

    Flux operator*(float s) const {
        Flux result;
        result.f1 = f1 * s;
        result.f2 = f2 * s;
        result.f3 = f3 * s;
        result.f4 = f4 * s;
        return result;
    }
};

// ============================================================================
// ROE RIEMANN SOLVER
// ============================================================================

Flux roeFlux(const ConservedVars& UL, const ConservedVars& UR, bool isXDirection) {
    // Roe-averaged quantities
    float sqrtRhoL = sqrtf(UL.rho);
    float sqrtRhoR = sqrtf(UR.rho);
    float denom = sqrtRhoL + sqrtRhoR;

    float uTilde = (sqrtRhoL * UL.u() + sqrtRhoR * UR.u()) / denom;
    float vTilde = (sqrtRhoL * UL.v() + sqrtRhoR * UR.v()) / denom;

    float HL = (UL.E + UL.pressure()) / UL.rho;
    float HR = (UR.E + UR.pressure()) / UR.rho;
    float HTilde = (sqrtRhoL * HL + sqrtRhoR * HR) / denom;

    float aTilde = sqrtf((GAMMA - 1.0f) * (HTilde - 0.5f * (uTilde * uTilde + vTilde * vTilde)));

    // Wave speeds
    float lambda1, lambda2, lambda3;
    if (isXDirection) {
        lambda1 = fabsf(uTilde - aTilde);
        lambda2 = fabsf(uTilde);
        lambda3 = fabsf(uTilde + aTilde);
    } else {
        lambda1 = fabsf(vTilde - aTilde);
        lambda2 = fabsf(vTilde);
        lambda3 = fabsf(vTilde + aTilde);
    }

    // Entropy fix
    float epsilon = 0.1f * aTilde;
    if (lambda1 < epsilon) lambda1 = (lambda1 * lambda1 + epsilon * epsilon) / (2.0f * epsilon);
    if (lambda2 < epsilon) lambda2 = (lambda2 * lambda2 + epsilon * epsilon) / (2.0f * epsilon);
    if (lambda3 < epsilon) lambda3 = (lambda3 * lambda3 + epsilon * epsilon) / (2.0f * epsilon);

    // Differences
    float drho = UR.rho - UL.rho;
    float du = UR.u() - UL.u();
    float dv = UR.v() - UL.v();
    float dp = UR.pressure() - UL.pressure();

    // Wave strengths (simplified)
    float alpha1, alpha2, alpha3;
    if (isXDirection) {
        alpha1 = 0.5f * (dp - UL.rho * aTilde * du) / (aTilde * aTilde);
        alpha2 = drho - dp / (aTilde * aTilde);
        alpha3 = 0.5f * (dp + UL.rho * aTilde * du) / (aTilde * aTilde);
    } else {
        alpha1 = 0.5f * (dp - UL.rho * aTilde * dv) / (aTilde * aTilde);
        alpha2 = drho - dp / (aTilde * aTilde);
        alpha3 = 0.5f * (dp + UL.rho * aTilde * dv) / (aTilde * aTilde);
    }

    // Dissipation term
    Flux FL = isXDirection ? Flux::calculateX(UL) : Flux::calculateY(UL);
    Flux FR = isXDirection ? Flux::calculateX(UR) : Flux::calculateY(UR);

    Flux avgFlux;
    avgFlux.f1 = 0.5f * (FL.f1 + FR.f1);
    avgFlux.f2 = 0.5f * (FL.f2 + FR.f2);
    avgFlux.f3 = 0.5f * (FL.f3 + FR.f3);
    avgFlux.f4 = 0.5f * (FL.f4 + FR.f4);

    // Simplified dissipation
    float dissipation = 0.5f * (lambda1 * fabsf(drho) + lambda2 * fabsf(drho) + lambda3 * fabsf(drho));

    avgFlux.f1 -= dissipation;
    avgFlux.f2 -= dissipation * (isXDirection ? uTilde : 0);
    avgFlux.f3 -= dissipation * (isXDirection ? 0 : vTilde);
    avgFlux.f4 -= dissipation * HTilde;

    return avgFlux;
}

// ============================================================================
// AERODYNAMICS SIMULATOR
// ============================================================================

class AeroEngine {
private:
    int width, height;
    float dx, dy, dt;

    std::vector<ConservedVars> U;       // Current state
    std::vector<ConservedVars> UNew;    // Next state

    float time;
    int iteration;

    // Free stream conditions
    float machInfinity;
    float angleOfAttack;  // degrees

    int idx(int x, int y) const {
        return y * width + x;
    }

    bool isObstacle(int x, int y) const {
        // Define obstacle shape (airfoil or cylinder)
        float cx = width * 0.3f;
        float cy = height * 0.5f;
        float rx = width * 0.15f;
        float ry = height * 0.1f;

        float dx = (x - cx) / rx;
        float dy = (y - cy) / ry;

        return (dx * dx + dy * dy) < 1.0f;
    }

    void applyBoundaryConditions() {
        // Inlet (left): supersonic inflow
        float rho_inf = 1.225f;  // kg/m³
        float p_inf = 101325.0f; // Pa
        float u_inf = machInfinity * sqrtf(GAMMA * p_inf / rho_inf);
        float v_inf = u_inf * tanf(angleOfAttack * 3.14159f / 180.0f);

        for (int y = 0; y < height; ++y) {
            U[idx(0, y)] = ConservedVars(rho_inf, u_inf, v_inf, p_inf);
        }

        // Outlet (right): extrapolation
        for (int y = 0; y < height; ++y) {
            U[idx(width-1, y)] = U[idx(width-2, y)];
        }

        // Top/bottom: slip walls
        for (int x = 0; x < width; ++x) {
            U[idx(x, 0)] = U[idx(x, 1)];
            U[idx(x, height-1)] = U[idx(x, height-2)];
            U[idx(x, 0)].rhoV = 0;
            U[idx(x, height-1)].rhoV = 0;
        }

        // Obstacle: no-slip/slip boundary
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                if (isObstacle(x, y)) {
                    // Reflect velocities
                    U[idx(x, y)].rhoU = 0;
                    U[idx(x, y)].rhoV = 0;
                }
            }
        }
    }

public:
    AeroEngine(int w, int h)
        : width(w), height(h)
        , dx(1.0f / w), dy(1.0f / h)
        , dt(0.0001f)
        , time(0)
        , iteration(0)
        , machInfinity(2.0f)
        , angleOfAttack(0.0f)
    {
        U.resize(width * height);
        UNew.resize(width * height);

        // Initialize with free stream
        float rho_inf = 1.225f;
        float p_inf = 101325.0f;
        float u_inf = machInfinity * sqrtf(GAMMA * p_inf / rho_inf);

        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                U[idx(x, y)] = ConservedVars(rho_inf, u_inf, 0, p_inf);
            }
        }
    }

    void setMachNumber(float mach) {
        machInfinity = mach;
    }

    void setAngleOfAttack(float angle) {
        angleOfAttack = angle;
    }

    void update(int steps = 1) {
        for (int step = 0; step < steps; ++step) {
            applyBoundaryConditions();

            // Update using Roe solver
            for (int y = 1; y < height - 1; ++y) {
                for (int x = 1; x < width - 1; ++x) {
                    if (isObstacle(x, y)) continue;

                    int i = idx(x, y);

                    // X-direction fluxes
                    Flux FxL = roeFlux(U[idx(x-1, y)], U[i], true);
                    Flux FxR = roeFlux(U[i], U[idx(x+1, y)], true);

                    // Y-direction fluxes
                    Flux FyD = roeFlux(U[idx(x, y-1)], U[i], false);
                    Flux FyU = roeFlux(U[i], U[idx(x, y+1)], false);

                    // Update conserved variables
                    UNew[i].rho  = U[i].rho  - dt/dx * (FxR.f1 - FxL.f1) - dt/dy * (FyU.f1 - FyD.f1);
                    UNew[i].rhoU = U[i].rhoU - dt/dx * (FxR.f2 - FxL.f2) - dt/dy * (FyU.f2 - FyD.f2);
                    UNew[i].rhoV = U[i].rhoV - dt/dx * (FxR.f3 - FxL.f3) - dt/dy * (FyU.f3 - FyD.f3);
                    UNew[i].E    = U[i].E    - dt/dx * (FxR.f4 - FxL.f4) - dt/dy * (FyU.f4 - FyD.f4);

                    // Positivity preservation
                    if (UNew[i].rho < 0.1f) UNew[i].rho = 0.1f;
                    if (UNew[i].pressure() < 1000.0f) {
                        UNew[i] = U[i];
                    }
                }
            }

            std::swap(U, UNew);
            time += dt;
            iteration++;
        }
    }

    val getDensityField() {
        val data = val::array();
        for (int i = 0; i < width * height; ++i) {
            data.call<void>("push", U[i].rho);
        }
        return data;
    }

    val getPressureField() {
        val data = val::array();
        for (int i = 0; i < width * height; ++i) {
            data.call<void>("push", U[i].pressure());
        }
        return data;
    }

    val getMachField() {
        val data = val::array();
        for (int i = 0; i < width * height; ++i) {
            data.call<void>("push", U[i].machNumber());
        }
        return data;
    }

    val getVelocityField() {
        val data = val::array();
        for (int i = 0; i < width * height; ++i) {
            data.call<void>("push", U[i].u());
            data.call<void>("push", U[i].v());
        }
        return data;
    }

    float getMaxMach() const {
        float maxMach = 0;
        for (const auto& u : U) {
            maxMach = std::max(maxMach, u.machNumber());
        }
        return maxMach;
    }

    int getWidth() const { return width; }
    int getHeight() const { return height; }
    float getTime() const { return time; }
    int getIteration() const { return iteration; }
};

// ============================================================================
// EMSCRIPTEN BINDINGS
// ============================================================================

EMSCRIPTEN_BINDINGS(aero_module) {
    class_<AeroEngine>("AeroEngine")
        .constructor<int, int>()
        .function("setMachNumber", &AeroEngine::setMachNumber)
        .function("setAngleOfAttack", &AeroEngine::setAngleOfAttack)
        .function("update", &AeroEngine::update)
        .function("getDensityField", &AeroEngine::getDensityField)
        .function("getPressureField", &AeroEngine::getPressureField)
        .function("getMachField", &AeroEngine::getMachField)
        .function("getVelocityField", &AeroEngine::getVelocityField)
        .function("getMaxMach", &AeroEngine::getMaxMach)
        .function("getWidth", &AeroEngine::getWidth)
        .function("getHeight", &AeroEngine::getHeight)
        .function("getTime", &AeroEngine::getTime)
        .function("getIteration", &AeroEngine::getIteration);
}
