# qkd_simulation/parameters.py

# --- Physical Constants ---
C = 299792458  # Speed of light in m/s
H = 6.62607015e-34 # Planck's constant in J*s

# --- System Parameters ---
WAVELENGTH = 1550e-9  # Operating wavelength in meters (1550 nm)

# Sender and Receiver Assembly
TX_APERTURE_DIAMETER = 0.5  # Transmitter telescope aperture in meters
RX_APERTURE_DIAMETER = 1.0  # Receiver telescope aperture in meters

# Detector characteristics
DETECTOR_EFFICIENCY = 0.80  # Quantum efficiency of the detectors
DARK_COUNT_RATE = 20      # Dark counts per second
DETECTOR_DEAD_TIME = 50e-9  # Dead time in seconds

# --- Protocol Parameters ---
REPETITION_RATE = 1e9     # System pulse rate in Hz (1 GHz)
# MDI-QKD requires signal and decoy states
DECOY_INTENSITIES = {
    'signal': 0.5,
    'decoy': 0.1,
    'vacuum': 0.0
}

# --- Simulation Parameters ---
QBER_INTRINSIC = 0.01  # 1% intrinsic error from optical imperfections
ERROR_CORRECTION_EFFICIENCY = 1.16 # Typical efficiency factor for error correction