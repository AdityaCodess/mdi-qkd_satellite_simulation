## 🧭 Evolution to a Mission Design Tool

This project has evolved from a **simple link-budget calculator** into a **high-fidelity QKD Mission Design Tool**. The updates below summarize the major architectural, physical, and security-driven upgrades that define the current simulator.

---

## 🖥️ 1. Dashboard & UI Overhaul

### Dynamic 4-Box Mission Dashboard

A real-time mission header now displays:

- **PEAK SKR (ZENITH)**  
  Absolute best achievable secret key rate when the satellite is directly overhead.

- **AVERAGE SKR (PASS)**  
  Average secret key rate over the full slant-range sweep of a satellite pass.  
  This metric is critical for matching values reported in experimental literature.

- **TOTAL SYSTEM LOSS (dB)**  
  Full end-to-end loss budget, including:

  - Channel loss
  - Hardware / optical losses
  - Detector efficiency penalty

- **EFFECTIVE QBER (%)**  
  Real-time indicator of link security health.

### Integrated Plotting Engine

- Transitioned from static calculations to a **dynamic Matplotlib backend**.
- Performance curves update instantly when **RUN** is pressed.

### UI Scalability Improvements

- Added a **scrollable settings sidebar** to support an expanding mission-parameter set.

---

## 🌍 2. Physics Engine Upgrades

### Spherical Earth Geometry (GTS)

- Atmospheric channel modeling now accounts for **Earth curvature**.
- Atmospheric path length increases realistically toward the horizon.

### Satellite Altitude Clamping

- Introduced an explicit **Satellite Altitude** parameter.
- Prevents non-physical configurations (e.g., 10 km range for a 500 km satellite).
- Range sweeps always begin at **Zenith**.

### Dynamic Loss Budgeting

- Implemented the **Detector Efficiency Penalty**:

  -10 · log10(η)

- Converts previously fixed losses into **hardware-dependent dynamic losses**.

---

## 🔐 3. Protocol & Security Mathematics

### Decoy-State BB84 Refinement

- Integrated:
  - Error-correction efficiency (f_ec)
  - Intrinsic QBER
- Improved alignment with analytical and experimental QKD models.

### The “QBER Cliff” Logic

- Enforced the **11% security threshold**.
- When exceeded, the simulator correctly terminates the link (**SKR = 0**).
- Essential for reproducing experimentally observed range limits (e.g., GYS).

### MDI-QKD Engine

- Added a dedicated mathematical branch for **Measurement-Device-Independent QKD**.
- Correctly captures quadratic loss scaling and enhanced detector-side security.

---

## 🛰️ 4. Mission Verification & Calibration

### Micius (LEO QKD Satellite) Calibration

- **Beam Divergence:** 0.006 mrad
- **Dark Count Rate:** 3000 Hz (background-light inclusive)

**Verification Result**

- Average SKR ≈ 2.5e-05
- Total System Loss ≈ 28 dB  
  (Gold-standard agreement with published results)

---

### GYS (Ground-to-Satellite) Calibration

- **Wavelength:** 1550 nm
- **Detector Efficiency:** 4.5%

**Verification Result**

- Correct reproduction of the **140 km security-limited range**

---

## 📊 5. Multi-Graph Analysis Suite

The simulator now supports **10 verified analysis modes**, including:

- **SKR vs. Pointing Error**  
  (Satellite vibration and tracking tolerance testing)

- **SKR vs. Atmospheric Turbulence**  
  (Weather robustness analysis)

- **Key Rate Components**  
  (Useful key vs. information leakage)

- **QBER vs. Eve’s Attack Strength**  
  (Adversarial and cyber-security analysis)
