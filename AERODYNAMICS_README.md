# 🚀 Ultimate Aerodynamics Simulation Suite

Three comprehensive C++ WebAssembly simulators covering the full spectrum of atmospheric flight physics.

## 🎯 Overview

| Simulator | Atmosphere | Speed Range | Key Physics |
|-----------|------------|-------------|-------------|
| **Stratospheric Glider** | Thin (0-25km) | 20-80 m/s | Thermals, Ridge Lift, Soaring |
| **Hypersonic Flight** | All Altitudes | Mach 5-25 | Re-entry Heating, Plasma, TPS |
| **Vehicle Downforce** | Sea Level | 0-360 km/h | Ground Effect, DRS, Balance |

---

## ⛅ Stratospheric Glider Simulator

**File:** `stratospheric-glider-wasm.html` + `glider_engine.cpp`

### Physics Features

**4 Glider Types:**
- High-Performance: L/D 60 (ASH 31, Eta)
- Standard Class: L/D 45 (ASW 28)
- Paraglider: L/D 10
- Hang Glider: L/D 15

**Lift Sources:**
- **Thermals:** Gaussian bubble model with cloud base visualization
- **Ridge Lift:** Orographic wave mechanics from terrain
- **Mountain Waves:** Lee wave oscillations (stratospheric altitudes)

**Atmospheric Model:**
- ISA (International Standard Atmosphere) up to 80km
- Temperature, pressure, density profiles
- Wind shear and jet stream modeling (11km altitude)

**Flight Computer:**
- McCready Theory for speed-to-fly calculations
- Variometer with dual filtering (smoothed + averaged)
- Polar curve calculations (sink rate vs airspeed)
- Total energy management (potential + kinetic)

**Controls:**
- Bank angle (-45° to +45°)
- Water ballast (0-150 kg)
- McCready setting (0-5 m/s)

### Interactive Features

✅ **Fullscreen Mode** - Single button immersion
✅ **Pan Camera** - Mouse drag
✅ **Zoom** - Mouse wheel + buttons (0.2x - 3.0x)
✅ **Follow Mode** - Auto-track glider
✅ **Variometer Display** - Visual + numeric climb rate
✅ **Thermal Spawning** - Interactive lift source creation
✅ **Vario-Colored Trail** - Green = lift, Red = sink

---

## 🔥 Hypersonic Flight Simulator

**File:** `hypersonic-flight-wasm.html` + `hypersonic_engine.cpp`

### Physics Features

**Flight Regimes:**
- Mach 5-10: Hypersonic aerodynamics
- Mach 10-25: Plasma sheath formation
- Entry heating corridor analysis

**Aerothermal Heating:**
- **Convective Heating:** Fay-Riddell correlation
  - `q_conv ∝ √(ρ/R_n) × V³`
- **Radiative Heating:** Stefan-Boltzmann from shock layer
  - Becomes significant above Mach 10
- **TPS Modeling:** Thermal protection system with re-radiation
  - Emissivity-based cooling

**Real Gas Effects:**
- Specific heat ratio changes (γ: 1.4 → 1.2)
- O₂/N₂ dissociation at high temperatures
- Plasma formation above 10,000 K

**Vehicle Scenarios:**
1. **X-15:** Mach 6.7, 30km altitude
2. **X-43A:** Mach 9.6, 35km altitude
3. **Re-entry:** Mach 25, 100km altitude (orbital velocity)
4. **Skip Re-entry:** Mach 20, 70km altitude

### Telemetry

- Surface temperature (TPS outer layer)
- Internal temperature (structure)
- Stagnation point temperature
- Convective heat flux (kW/m²)
- Radiative heat flux (kW/m²)
- Total heat load integration (MJ)
- Dynamic pressure (kPa)
- Plasma sheath indicator

---

## 🏎️ Vehicle Downforce Simulator

**File:** `vehicle-downforce-wasm.html` + `vehicle_aero_engine.cpp`

### Physics Features

**4 Vehicle Types:**

| Vehicle | L/D | Key Feature | Mass |
|---------|-----|-------------|------|
| Formula 1 | 3.5 | DRS System | 798 kg |
| NASCAR | 1.5 | Rear Spoiler | 1542 kg |
| LMP1 | 3.0 | Efficiency | 875 kg |
| Road Car | 0.5 | Slight Lift | 1400 kg |

**Aerodynamic Components:**
- **Front Wing:** Adjustable angle (-10° to +15°)
- **Rear Wing:** Adjustable angle (-5° to +20°)
- **DRS (Drag Reduction System):** F1 only, 70% rear downforce reduction
- **Underbody Diffuser:** Ground effect dependent
- **Ride Height Sensitivity:** Exponential ground effect curves

**Ground Effect Physics:**
- Optimal ride height (F1: 30mm, NASCAR: 50mm)
- Venturi effect at very low heights
- Lift augmentation and drag reduction

**Force Calculations:**
```
Dynamic Pressure: q = 0.5 × ρ × v²
Downforce: F = Cl × q × Area
Induced Drag: Cd_induced = K × Cl²
```

**Track Simulation:**
- Straights and corners
- Speed-dependent cornering limits
- Lateral G-force tracking (up to 4g)

### Interactive Features

✅ **Fullscreen Mode**
✅ **Pan/Zoom Camera**
✅ **3 View Modes:**
- Track View: Overhead with speed-colored trail
- Closeup View: Follows vehicle with aero arrows
- Pressure View: Downforce distribution

✅ **Real-Time Telemetry:**
- Speed, downforce, drag, L/D ratio
- Front/rear downforce split
- Center of pressure position
- Tire loads (4 corners)
- Ground effect percentage

---

## 🛠️ Compilation Instructions

### Prerequisites

Install Emscripten SDK:
```bash
git clone https://github.com/emscripten-core/emsdk.git
cd emsdk
./emsdk install latest
./emsdk activate latest
source ./emsdk_env.sh
```

### Compile All Simulators

```bash
chmod +x compile_aero.sh
./compile_aero.sh
```

### Manual Compilation

**Glider:**
```bash
emcc glider_engine.cpp -o glider_engine.js \
  -s WASM=1 -s MODULARIZE=1 -s EXPORT_ES6=1 \
  -s EXPORT_NAME="createGliderModule" \
  -lembind -O3 -s ENVIRONMENT=web
```

**Hypersonic:**
```bash
emcc hypersonic_engine.cpp -o hypersonic_engine.js \
  -s WASM=1 -s MODULARIZE=1 -s EXPORT_ES6=1 \
  -s EXPORT_NAME="createHypersonicModule" \
  -lembind -O3 -s ENVIRONMENT=web
```

**Vehicle:**
```bash
emcc vehicle_aero_engine.cpp -o vehicle_aero_engine.js \
  -s WASM=1 -s MODULARIZE=1 -s EXPORT_ES6=1 \
  -s EXPORT_NAME="createVehicleModule" \
  -lembind -O3 -s ENVIRONMENT=web
```

---

## 🚀 Running the Simulators

1. **Compile WASM modules** (see above)
2. **Start local server:**
   ```bash
   python3 -m http.server 8000
   ```
3. **Open browser:**
   ```
   http://localhost:8000/stratospheric-glider-wasm.html
   http://localhost:8000/hypersonic-flight-wasm.html
   http://localhost:8000/vehicle-downforce-wasm.html
   ```

---

## 📊 Performance

**Typical File Sizes:**
- C++ Source: 15-22 KB
- Compiled JS: ~43 KB
- WASM Binary: 25-30 KB
- Total Load: ~70 KB per simulator

**Physics Update Rate:**
- 60 FPS rendering
- 2-5 physics substeps per frame
- Total: 120-300 physics updates/second

**Particle/Object Counts:**
- Glider: Up to 50 thermals tracked
- Hypersonic: 1000-point temperature trail
- Vehicle: 6 track sections, 2000-point trail

---

## 🎮 Controls Summary

### All Simulators

| Control | Action |
|---------|--------|
| **Mouse Drag** | Pan camera |
| **Mouse Wheel** | Zoom in/out |
| **Fullscreen Button** | Immersive mode |
| **ESC** | Exit fullscreen |

### Glider-Specific

| Control | Action |
|---------|--------|
| **Bank Slider** | Turn left/right |
| **Ballast Slider** | Add water weight |
| **McCready Slider** | Set speed-to-fly |
| **Spawn Thermal Button** | Create lift source |
| **Follow Button** | Auto-track glider |

### Hypersonic-Specific

| Control | Action |
|---------|--------|
| **Velocity Slider** | Set speed (1000-7500 m/s) |
| **Altitude Slider** | Set starting height |
| **AoA Slider** | Angle of attack |
| **Scenario Buttons** | Pre-set missions |

### Vehicle-Specific

| Control | Action |
|---------|--------|
| **Front Wing** | Adjust downforce balance |
| **Rear Wing** | Adjust rear grip |
| **Ride Height** | Ground effect tuning |
| **DRS Button** | Reduce drag (F1 only) |
| **View Mode** | Switch visualization |

---

## 🧪 Physics Accuracy

**Glider:**
- Real polar curves from ASH 31, ASW 28
- Thermal models from meteorology research
- McCready theory (1958) implementation

**Hypersonic:**
- Fay-Riddell equation (1958)
- Newtonian impact theory for high Mach
- Real TPS properties (emissivity 0.85)
- ISA atmosphere extended to 100km

**Vehicle:**
- F1 regulations 2023 (mass, dimensions)
- Ground effect based on wind tunnel data
- Tire load sensitivity models
- Real drag coefficients (F1: Cd ~0.7)

---

## 📚 References

1. Anderson, J.D. "Hypersonic and High-Temperature Gas Dynamics" (2006)
2. Reichmann, H. "Cross-Country Soaring" (1993)
3. Katz, J. "Race Car Aerodynamics" (1995)
4. NASA Technical Reports (X-15, X-43A programs)
5. Wieselsberger, C. "Wing Resistance Near the Ground" (1922)

---

## 🎯 Future Enhancements

- [ ] Audio feedback (variometer beeps, sonic booms)
- [ ] Multiplayer racing (vehicle simulator)
- [ ] Weather systems (convection models)
- [ ] VR support
- [ ] Machine learning for optimal flight paths
- [ ] Real-world terrain import (glider)

---

## 📜 License

Part of Soupy Labs Physics Simulation Suite
Educational and research purposes

---

**🚀 Happy Flying! 🚀**
