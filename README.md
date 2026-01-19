# SOUPY LABS - Technical Physics Simulation Suite

**Version:** 2.0
**Author:** John Beene
**Location:** Glen St. Mary, FL
**Date:** December 2024

---

## ☕ Support This Project

**Love these free physics simulations?** Help keep them running and support more awesome content!

### **[☕ Buy Me a Coffee on Ko-fi](https://ko-fi.com/soupylab)** ⚡

Your support helps fund:
- New simulator development (quantum mechanics, relativity, particle physics)
- Server costs and domain hosting
- Advanced GPU visualization features
- Educational content and tutorials
- More physics awesomeness!

Every coffee counts! Even $3 helps keep the physics flowing. 🚀

---

## Overview

Soupy Labs is a comprehensive browser-based physics simulation suite featuring:
- Novel fusion reactor concepts (Ultra Magnetic Sun, Tokamak/Stellarator comparison)
- Plasma physics fundamentals (PIC, MHD instabilities, Debye shielding)
- GPU-accelerated field visualization
- Computational fluid dynamics

All 12 modules run entirely client-side using WebGL2 for hardware acceleration. No server required.

---

## Modules

### FUSION & REACTOR CONCEPTS

#### 1. Ultra Magnetic Sun Reactor ⭐ FLAGSHIP
- **File:** `ultra-magnetic-sun.html`
- **Concept:** Novel spherical plasma confinement with spring-mounted coils
- **Key Innovation:** Coils move WITH plasma instabilities rather than fighting them
- **Features:** Spring-damped dynamics, CME capture loops, ΨΦ reactor face

#### 2. Tokamak vs Stellarator Comparison ⭐ NEW
- **File:** `fusion-compare.html`
- **Method:** Side-by-side real-time comparison
- **Features:** Grouped parameter controls (Geometry, Plasma, Magnetic, Heating, Stability)
- **Physics:** ITER scaling (tokamak), ISS04 scaling (stellarator), Lawson criterion, Q-factor

#### 3. Fusion Cross Sections ⭐ NEW
- **File:** `fusion-xsec.html`
- **Method:** Bosch-Hale parametrization for reactivity curves
- **Reactions:** D-T, D-D, D-He³, p-B¹¹ (aneutronic)
- **Features:** Interactive temperature slider, Lawson criterion, ignition regime

#### 4. Particle Accelerator
- **File:** `accelerator.html`
- **Types:** Cyclotron, synchrotron, linac
- **Features:** RF cavity acceleration, magnetic steering, relativistic effects

### PLASMA PHYSICS

#### 5. MHD Instability Sandbox ⭐ NEW
- **File:** `instability.html`
- **Modes:** Kink (m=1), Sausage (m=0), Ballooning, Tearing, Rayleigh-Taylor, ELM
- **Parameters:** Adjustable β, safety factor q, pressure gradient, magnetic shear
- **Physics:** Kruskal-Shafranov limit, disruption probability

#### 6. Particle-in-Cell (PIC) ⭐ NEW
- **File:** `pic.html`
- **Algorithm:** Boris pusher for charged particle motion
- **Presets:** Single gyration, magnetic mirror, E×B drift, thermal plasma
- **Physics:** Larmor radius, cyclotron frequency, drift velocities

#### 7. Debye Shielding ⭐ NEW
- **File:** `debye.html`
- **Concept:** Electrostatic screening in plasma
- **Features:** Drop test charges, watch shielding clouds form
- **Physics:** Debye length λ_D, plasma frequency ω_p, Debye sphere particles

#### 8. Magnetic Reconnection
- **File:** `reconnection.html`
- **Concept:** Field line topology change releasing magnetic energy
- **Features:** X-point formation, reconnection jets, energy release
- **Relevance:** Solar flares, tokamak disruptions, magnetotail substorms

### FIELDS & SIMULATION

#### 9. MagneticGPU
- **File:** `magnetic-gpu.html`
- **Method:** GPU-accelerated vector field rendering (WebGL2)
- **Sources:** Dipoles, toroids, solenoids, current-carrying wires
- **Features:** 5-point compass sampling, draggable probe, multiple render modes

#### 10. EM Wave Propagation
- **File:** `em-wave.html`
- **Method:** FDTD (Finite-Difference Time-Domain)
- **Features:** Yee grid, PML absorbing boundaries, multiple source types

#### 11. Lattice-Boltzmann Fluid
- **File:** `fluid.html`
- **Method:** D2Q9 velocity discretization, BGK collision
- **Features:** Real-time Reynolds control, velocity/vorticity/pressure views

#### 12. Verlet Physics Engine
- **File:** `physics.html`
- **Method:** Position-based dynamics, Gauss-Seidel constraint relaxation
- **Features:** Soft bodies, cloth, rope, spatial hash collision

---

## Ultra Magnetic Sun Reactor - Technical Details

### Core Architecture

Unlike conventional tokamaks (toroidal) or stellarators (twisted torus), the Ultra Magnetic Sun maintains a **spherical plasma core** surrounded by:

1. **Spring-Mounted Mega Coils**
   - Superconducting magnets on mechanical dampers
   - Move with plasma fluctuations, not against them
   - Convert MHD instability energy to mechanical motion

2. **Multi-Path Rail Network**
   - Redundant magnetic pathways
   - Plasma self-selects lowest-energy configuration
   - Reduces disruption cascade risk

3. **CME Capture Loops**
   - Outer rings intercept plasma ejections
   - Convert ejection kinetic energy to induced current
   - Turns instabilities into harvested power

### Physics Parameters

| Parameter | Description | Target |
|-----------|-------------|--------|
| nτT (Lawson) | Triple product | >3×10²¹ keV·s/m³ |
| Q Factor | Power gain | >1 for net energy |
| β | Plasma/magnetic pressure ratio | Higher = more efficient |
| τ | Confinement time | Enhanced by spring damping |

---

## MagneticGPU - Technical Details

### Field Equations

- **Dipole:** B = (μ₀/4π) · (3(m·r̂)r̂ - m) / r³
- **Toroid:** B_φ = μ₀NI / 2πr
- **Solenoid:** B = μ₀nI (uniform interior)
- **Wire (Biot-Savart):** B = μ₀I / 2πρ

### 5-Point Compass Method

Each vector glyph samples the field at 5 points:
- Center position
- +X offset, -X offset
- +Y offset, -Y offset

The averaged direction provides smoothed field orientation while capturing local gradients, similar to iron filings aligning to magnetic field lines.

---

## Installation

1. Extract all files to a directory
2. Open `index.html` in a modern browser (Chrome, Firefox, Edge)
3. Or deploy to any static hosting (Netlify, GitHub Pages, etc.)

**Requirements:**
- WebGL2-capable browser
- No server-side dependencies
- No build step required

---

## API (for developers)

### MagneticGPU Integration

```javascript
// Initialize
const magneticGPU = new MagneticGPU(canvas, {
    fieldScale: 1.0,
    vectorLength: 1.0,
    sampleDensity: 32,
    colormap: 'coolwarm',
    renderMode: 'needle'
});

// Add field sources
magneticGPU.addFieldSource({
    type: 'dipole',      // dipole, toroid, solenoid, coil, wire
    position: [0, 0, 0],
    strength: 1.0,
    moment: [0, 1, 0]    // for dipoles
});

// Custom field function
magneticGPU.setFieldSource((x, y, z) => {
    // Return [Bx, By, Bz] at position
    return [0, Math.sin(x), 0];
});

// Render loop
function animate(dt) {
    magneticGPU.update(dt);
    magneticGPU.render();
}
```

---

## License

This work is provided for educational and research collaboration purposes.

For commercial licensing or collaboration inquiries, contact John Beene.

---

## Acknowledgments

- Plasma physics concepts informed by established MHD theory
- WebGL implementation uses standard techniques for GPU compute
- PsiPhi face visualization inspired by affective computing research
