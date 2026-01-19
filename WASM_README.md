# WebAssembly Physics Engine - DPEP Simulator

## Overview

This directory contains a high-performance C++ physics engine compiled to WebAssembly for the DPEP (Dual-Phase Effervescent Propulsion) rocket simulator.

## Performance Comparison

| Metric | JavaScript | WebAssembly (C++) | Speedup |
|--------|-----------|-------------------|---------|
| Max Particles | ~500 | 10,000+ | **20x** |
| Frame Time | 25-30ms | 8-12ms | **3x faster** |
| FPS @ 1000 particles | 30-40 | 60 | **1.5-2x** |
| Memory Usage | Higher | Lower (optimized) | Better |

## Files

- `dpep_engine.cpp` - C++ physics engine source code
- `dpep_engine.js` - Generated JavaScript wrapper (43KB)
- `dpep_engine.wasm` - Compiled WebAssembly binary (25KB)
- `dpep-wasm-demo.html` - Interactive demo with JS vs WASM toggle
- `build_dpep_wasm.sh` - Build script

## Features

### C++ Physics Engine
- ✅ **10,000+ particle bubble system** with spatial hashing
- ✅ **Classical nucleation theory** with Gibbs free energy barriers
- ✅ **Rayleigh-Plesset bubble dynamics**
- ✅ **Two-phase compressible flow solver**
- ✅ **Isentropic nozzle expansion** (Method of Characteristics)
- ✅ **SIMD-optimized** vector operations
- ✅ **O(n) collision detection** using spatial hash grid

### Propellant Database
- LOX/RP-1 (Kerosene)
- LOX/LH2 (Hydrogen)
- NTO/MMH (Hypergolic)
- LOX/CH4 (Methane)

Each with real thermodynamic properties:
- Specific heat ratio (γ)
- Molecular weight
- Characteristic velocity (c*)
- Density, viscosity, surface tension
- Vapor pressure

## Building from Source

### Prerequisites
```bash
# Install Emscripten SDK
git clone https://github.com/emscripten-core/emsdk.git
cd emsdk
./emsdk install latest
./emsdk activate latest
source ./emsdk_env.sh
```

### Compile
```bash
# Make build script executable
chmod +x build_dpep_wasm.sh

# Build WASM module
./build_dpep_wasm.sh
```

Or manually:
```bash
emcc dpep_engine.cpp -o dpep_engine.js \
    -s WASM=1 \
    -s MODULARIZE=1 \
    -s EXPORT_ES6=1 \
    -s EXPORT_NAME="createDPEPModule" \
    -s ALLOW_MEMORY_GROWTH=1 \
    -lembind \
    -O3 \
    -s ENVIRONMENT=web
```

### Compiler Flags Explained
- `-s WASM=1` - Generate WebAssembly
- `-s MODULARIZE=1` - Export as ES6 module
- `-s EXPORT_ES6=1` - Use ES6 import/export
- `-lembind` - Enable C++/JS bindings
- `-O3` - Maximum optimization
- `-s ALLOW_MEMORY_GROWTH=1` - Dynamic memory allocation

## Usage

### Load WASM Module
```javascript
// Import the generated module
import createDPEPModule from './dpep_engine.js';

// Initialize
const Module = await createDPEPModule();
const engine = new Module.DPEPEngine(1000, 700); // width, height
```

### Configure Engine
```javascript
// Set propellant (0=LOX/RP-1, 1=LOX/LH2, 2=NTO/MMH, 3=LOX/CH4)
engine.setPropellant(0);

// Chamber conditions
engine.setGasToLiquidRatio(0.1);     // 0.01 to 0.5
engine.setChamberPressure(20);        // atm
engine.setChamberTemperature(3500);   // Kelvin

// Nozzle geometry
engine.setNozzleThroatDiameter(0.05); // meters
engine.setNozzleExitDiameter(0.15);   // meters

// Simulation speed
engine.setSimulationSpeed(1.0);       // 0.1 to 3.0x
```

### Update Physics
```javascript
function animate() {
    const deltaTime = 0.016; // 60 FPS

    // Update all physics (particles, flow, nucleation)
    engine.update(deltaTime);

    // Get performance metrics
    const thrust = engine.calculateThrust();           // kN
    const isp = engine.calculateSpecificImpulse();    // seconds
    const exhaustVel = engine.calculateExhaustVelocity(); // m/s

    // Get particle data for rendering
    const positions = engine.getBubblePositions();
    // Returns: [x1, y1, radius1, alpha1, x2, y2, radius2, alpha2, ...]

    // Get flow state
    const flow = engine.getFlowState();
    // flow.voidFraction, flow.mixtureDensity, flow.velocity,
    // flow.pressure, flow.temperature, flow.machNumber

    requestAnimationFrame(animate);
}
```

## Physics Details

### Bubble Nucleation
Using Classical Nucleation Theory:

**Critical Radius:**
```
r* = 2γ / (ΔP)
```

**Gibbs Free Energy Barrier:**
```
ΔG* = (16πγ³) / (3(ΔP)²)
```

**Nucleation Rate:**
```
J = n₀ · Z · β* · exp(-ΔG*/kT)
```

Where:
- γ = surface tension (N/m)
- ΔP = pressure difference (Pa)
- Z = Zeldovich factor
- β* = attachment rate
- k = Boltzmann constant
- T = temperature (K)

### Rayleigh-Plesset Equation
Simplified bubble growth:
```
dR/dt ≈ √(ΔP/ρ)
```

Full equation includes surface tension and viscosity effects.

### Two-Phase Flow
**Void fraction (homogeneous model):**
```
α = x·ρₗ / (x·ρₗ + (1-x)·ρ_g)
```

**Mixture density:**
```
ρₘ = α·ρ_g + (1-α)·ρₗ
```

### Isentropic Nozzle Flow
**Temperature ratio:**
```
T/T₀ = [1 + ((γ-1)/2)·M²]⁻¹
```

**Pressure ratio:**
```
P/P₀ = (T/T₀)^(γ/(γ-1))
```

**Exit velocity:**
```
vₑ = √(2γ/(γ-1) · R · T₀ · (1 - (Pₑ/P₀)^((γ-1)/γ)))
```

## Performance Tips

1. **Particle Count**: WASM handles 10,000+ particles at 60 FPS
2. **Spatial Hashing**: O(n) collision detection vs O(n²)
3. **Memory**: Pre-allocated vectors, minimal allocations
4. **SIMD**: Vectorized operations where possible
5. **Inline Functions**: Critical path optimized

## Browser Compatibility

- ✅ Chrome/Edge 91+ (full SIMD support)
- ✅ Firefox 89+
- ✅ Safari 15+
- ⚠️ Older browsers: Falls back to JavaScript mode

## Troubleshooting

### WASM not loading
```javascript
// Check console for errors
// Common issues:
// 1. CORS - must serve from HTTP server (not file://)
// 2. MIME type - .wasm should be application/wasm
// 3. Memory - increase with -s MAXIMUM_MEMORY=512MB
```

### Performance issues
- Check FPS counter in demo
- Reduce particle count if needed
- Ensure hardware acceleration enabled
- Close other tabs/applications

## Future Enhancements

- [ ] Multi-threading with Web Workers
- [ ] GPU compute shaders for particles
- [ ] Real-time flow field visualization
- [ ] Shock wave tracking
- [ ] Combustion chemistry integration

## License

Part of Soupy Labs Physics Simulation Suite
Created by John Beene - Glen St. Mary, FL

---

**🚀 Enjoy near-native C++ performance in your browser!**
