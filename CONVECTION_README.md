# Thermal Convection Simulator - WebAssembly Physics Engine

## Overview

High-performance C++ Rayleigh-Bénard convection simulator compiled to WebAssembly. Solves the full Navier-Stokes equations coupled with heat transfer to simulate natural convection patterns.

## Features

### Physics Engine
- ✅ **Full Navier-Stokes Solver** (vorticity-stream function formulation)
- ✅ **Heat Diffusion Equation** with advection
- ✅ **Buoyancy Forces** (Boussinesq approximation)
- ✅ **Rayleigh-Bénard Instability** patterns
- ✅ **Thermal Plume Formation**
- ✅ **Convection Cell Development**

### Governing Equations

**Vorticity Transport:**
```
∂ω/∂t = -u·∇ω + ν∇²ω + gβ(∂T/∂x)
```

**Stream Function (Poisson):**
```
∇²ψ = -ω
```

**Velocity from Stream Function:**
```
u = ∂ψ/∂y
v = -∂ψ/∂x
```

**Heat Equation:**
```
∂T/∂t = -u·∇T + α∇²T
```

Where:
- ω = vorticity (curl of velocity)
- ψ = stream function
- u, v = velocity components
- T = temperature
- ν = kinematic viscosity
- α = thermal diffusivity
- β = thermal expansion coefficient
- g = gravitational acceleration

## Dimensionless Numbers

### Rayleigh Number
```
Ra = gβΔTL³/(να)
```
Indicates convection strength. Critical value ~1708 for onset of convection.

### Prandtl Number
```
Pr = ν/α
```
Ratio of momentum diffusivity to thermal diffusivity.

### Nusselt Number
```
Nu = Q_convective / Q_conductive
```
Ratio of convective to conductive heat transfer.

## Fluid Properties

| Fluid | Density (kg/m³) | Kinematic Viscosity (m²/s) | Thermal Diffusivity (m²/s) | Prandtl Number |
|-------|----------------|---------------------------|---------------------------|----------------|
| **Air (300K)** | 1.184 | 1.5×10⁻⁵ | 2.2×10⁻⁵ | 0.68 |
| **Water (300K)** | 997 | 1.0×10⁻⁶ | 1.4×10⁻⁷ | 7.1 |
| **Oil (SAE 30)** | 920 | 1.0×10⁻⁴ | 8.0×10⁻⁸ | 1250 |
| **Glycerin** | 1260 | 1.2×10⁻³ | 9.5×10⁻⁸ | 12631 |
| **Liquid Sodium** | 927 | 6.8×10⁻⁷ | 6.7×10⁻⁵ | 0.01 |

## Numerical Methods

### Finite Difference Discretization
- **5-point stencil** for Laplacian operator
- **Central differences** for gradients
- **Explicit time stepping** for advection terms

### Poisson Solver
- **Jacobi iteration** for stream function
- Typically 30-50 iterations per timestep
- Convergence criterion based on residual

### Boundary Conditions
- **Fixed temperature**: Hot bottom, cold top
- **Adiabatic walls**: No heat flux through sides
- **No-slip walls**: Zero velocity at boundaries

### Stability
Time step constrained by:
- CFL condition: `dt < min(dx/|u|, dy/|v|)`
- Diffusion stability: `dt < dx²/(2ν)`

## Visualization Modes

### 1. Temperature Field 🌡️
Shows heat distribution from cold (blue) to hot (red):
- **Blue/Cyan**: Cold regions
- **Green/Yellow**: Intermediate temperatures
- **Orange/Red**: Hot regions

### 2. Vorticity Field 🌀
Shows fluid rotation (curl of velocity):
- **Blue**: Counter-clockwise rotation
- **Red**: Clockwise rotation
- Highlights eddy formation and turbulence

### 3. Velocity Field 💨
Shows flow speed magnitude:
- **Black/Purple**: Low velocity
- **Cyan**: Medium velocity
- **White**: High velocity
- Reveals convection cell structure

## Building & Usage

### Compile from Source
```bash
emcc convection_engine.cpp -o convection_engine.js \
    -s WASM=1 \
    -s MODULARIZE=1 \
    -s EXPORT_ES6=1 \
    -s EXPORT_NAME="createConvectionModule" \
    -s ALLOW_MEMORY_GROWTH=1 \
    -s MAXIMUM_MEMORY=512MB \
    -lembind \
    -O3 \
    -s ENVIRONMENT=web
```

### JavaScript Integration
```javascript
// Import module
import createConvectionModule from './convection_engine.js';

// Initialize
const Module = await createConvectionModule();
const engine = new Module.ConvectionEngine(200, 150); // width, height

// Configure
engine.setFluid(0);                  // 0=Air, 1=Water, etc.
engine.setHotTemperature(330);       // Kelvin
engine.setColdTemperature(290);      // Kelvin
engine.setTimeStep(0.0001);          // seconds

// Update simulation
engine.update(10); // 10 timesteps

// Get fields for visualization
const tempField = engine.getTemperatureField();
const vortField = engine.getVorticityField();
const velField = engine.getVelocityField();

// Get statistics
const Ra = engine.getRayleighNumber();
const Pr = engine.getPrandtlNumber();
const Nu = engine.getNusseltNumber();
const maxVel = engine.getMaxVelocity();
```

## Performance

### Grid Resolution
- Default: 200×150 cells
- Memory: ~480KB (5 fields × 200×150×4 bytes)
- Update rate: 60 FPS with 10-50 timesteps/frame

### Optimization Techniques
- **Flat array layout** for cache coherency
- **In-place updates** minimize allocations
- **SIMD-friendly** memory access patterns
- **O3 optimization** with compiler

## Physical Phenomena

### Convection Onset
When Ra > Ra_critical (~1708):
- Conduction becomes unstable
- Thermal plumes begin to rise
- Organized convection cells form

### Flow Regimes

**Laminar Convection** (Ra < 10⁴):
- Regular convection cells
- Stable thermal plumes
- Predictable flow patterns

**Transitional** (10⁴ < Ra < 10⁶):
- Cell pattern fluctuations
- Plume oscillations
- Increased Nusselt number

**Turbulent** (Ra > 10⁶):
- Chaotic flow patterns
- Multiple scales of motion
- Maximum heat transfer

### Boundary Layer Theory
Near hot bottom wall:
- Thin thermal boundary layer
- Strong temperature gradients
- High vorticity generation

## Applications

This simulator demonstrates:
- **Atmospheric Convection** (thunderstorm formation)
- **Oceanic Circulation** (thermohaline currents)
- **Mantle Convection** (plate tectonics driver)
- **Industrial Heat Transfer** (cooling systems)
- **Stellar Interiors** (energy transport)

## References

### Classical Papers
- Lord Rayleigh (1916) - "On Convection Currents in a Horizontal Layer of Fluid"
- Bénard (1900) - Original experimental observations

### Textbooks
- Tritton, D.J. (1988) - "Physical Fluid Dynamics"
- Chandrasekhar, S. (1961) - "Hydrodynamic and Hydromagnetic Stability"

### Numerical Methods
- Peyret, R. & Taylor, T.D. (1983) - "Computational Methods for Fluid Flow"

## Known Limitations

1. **2D Simulation**: Real convection is inherently 3D
2. **Boussinesq Approximation**: Valid for small temperature differences
3. **Incompressible Flow**: Cannot capture density variations
4. **No Phase Change**: No boiling/condensation
5. **Fixed Boundaries**: No free surface effects

## Future Enhancements

- [ ] 3D convection with WebGL compute
- [ ] Compressible flow formulation
- [ ] Phase change (boiling/condensation)
- [ ] Multiple fluids (stratification)
- [ ] Magnetic field effects (MHD)
- [ ] Adaptive mesh refinement

---

**Part of Soupy Labs Physics Simulation Suite**
Created by John Beene - Glen St. Mary, FL

🔥 **Experience fluid dynamics at C++ speeds in your browser!** 🌊
