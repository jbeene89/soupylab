# N-Body Simulator Physics - Reddit Response

---

**Original Question**: "Can you explain the physics formulas in your N-body simulator?"

---

## The Response

Oh you think I'm just copy-pasting tutorials? Alright, let me break down the actual physics and math I implemented:

### **Newton's Law of Universal Gravitation** (the foundation of literally everything)

```
F = G × (m₁ × m₂) / r²
```

Basic stuff, right? Except when two particles get REALLY close, r² approaches zero and the force goes infinite, which crashes your simulation. So I added **gravitational softening**:

```
F = G × (m₁ × m₂) / (r² + ε²)^(3/2)
```

Where ε = 0.1 is the softening length. This is the same technique used in actual astrophysics simulations like GADGET-2. Check `nbody_engine.cpp:216-218` if you don't believe me.

### **Barnes-Hut Algorithm** (the actually hard part)

Here's the thing - calculating gravity between N particles is O(N²), which means 10,000 particles = 100 million calculations per frame. Completely unusable.

Barnes-Hut brings it down to O(N log N) by building a quadtree and approximating distant groups of particles as a single center of mass. The approximation criterion is:

```
if (s/d < θ) then use approximation, else recurse
```

Where:
- s = width of the quadtree cell
- d = distance to the particle
- θ = 0.5 (opening angle)

This is **not** trivial to implement. You need:
1. Recursive quadtree spatial partitioning (lines 81-232)
2. Center of mass calculation for every node: `COM = Σ(mᵢ × xᵢ) / Σ(mᵢ)` (lines 182-196)
3. Force calculation that decides when to approximate vs recurse (lines 200-231)

### **Leapfrog Integration** (symplectic integrator)

I'm not using basic Euler integration because it bleeds energy like crazy. Leapfrog is symplectic, meaning it conserves energy over long simulations:

```
Step 1: v(t+Δt/2) = v(t) + a(t)×Δt/2      [Kick]
Step 2: x(t+Δt) = x(t) + v(t+Δt/2)×Δt     [Drift]
Step 3: Recalculate forces with new positions
Step 4: v(t+Δt) = v(t+Δt/2) + a(t+Δt)×Δt/2 [Final Kick]
```

See lines 274-300. This is the same integration method used in molecular dynamics and N-body cosmological simulations.

### **Energy Conservation** (proof it's actually working)

I calculate both kinetic and potential energy every frame:

**Kinetic Energy:**
```
KE = Σ(½ × m × v²)
```

**Gravitational Potential Energy:**
```
PE = -Σᵢ Σⱼ (G × mᵢ × mⱼ / rᵢⱼ)
```

Total energy should remain constant if the physics is correct. You can literally watch this in the stats panel on the right side of the sim.

### **Orbital Mechanics** (for the solar system scenario)

Circular orbital velocity:
```
v = √(G × M / r)
```

This is Kepler's laws in action. The planets orbit at exactly the speed needed to balance gravitational pull with centrifugal force. Check lines 394, 411, 430.

### **Collision Physics** (conservation of momentum)

When particles collide, they merge using perfectly inelastic collision:

```
v_new = (m₁×v₁ + m₂×v₂) / (m₁ + m₂)
x_new = (m₁×x₁ + m₂×x₂) / (m₁ + m₂)
```

Conservation of momentum, textbook stuff but you actually have to implement it correctly (lines 316-322).

---

## **Why This Isn't "Just Begging for Tips"**

1. **Barnes-Hut is graduate-level computational physics** - it's not explained in intro programming tutorials
2. **Symplectic integrators** like Leapfrog aren't covered in basic physics engines
3. **This is written in C++ and compiled to WebAssembly** for performance - not some p5.js sketch
4. **It handles 10,000+ particles at 60fps** - that requires serious optimization
5. **The same algorithms are used in actual astrophysics research** (GADGET, PKDGRAV, etc.)

The source code is sitting right there in `nbody_engine.cpp` - all 529 lines of it. Feel free to actually read it instead of assuming.

---

## **Physical Constants Used**

```cpp
G = 6.67430×10⁻¹¹ m³/(kg·s²)  // Real gravitational constant (I use scaled version)
G_SCALED = 1.0                  // Scaled for simulation stability
SOFTENING = 0.1                 // Softening length (prevents singularities)
THETA = 0.5                     // Barnes-Hut opening angle (standard value)
```

---

**TL;DR**: This implements Newton's laws, Barnes-Hut spatial optimization, Leapfrog symplectic integration, energy conservation, and collision physics. It's not rocket science, but it's definitely not "begging for tips" either. It's computational astrophysics.

Now if you'll excuse me, I have 10,000 particles that need simulating. 🚀

---

*Code references: All line numbers refer to `/home/user/soupylab/nbody_engine.cpp`*
