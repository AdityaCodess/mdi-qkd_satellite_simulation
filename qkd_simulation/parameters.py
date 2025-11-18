C = 299792458  # Speed of light in m/s
H = 6.62607015e-34 # Planck's constant in J*s

# --- System Parameters ---
# WAVELENGTH is dynamic
# DETECTOR_EFFICIENCY is dynamic

# Sender and Receiver Assembly (Example values, can be adjusted)
TX_APERTURE_DIAMETER = 0.3  # Micius transmitter aperture
RX_APERTURE_DIAMETER = 1.2  # Example: Delingha OGS receiver aperture

# Detector characteristics (Fixed properties)
DARK_COUNT_RATE = 500     # Dark counts per second
DETECTOR_DEAD_TIME = 50e-9  # Typical dead time

# --- Protocol Parameters ---
REPETITION_RATE = 1e8     # System pulse rate in Hz 
DECOY_INTENSITIES = {
    'signal': 0.5, 'decoy': 0.1, 'vacuum': 0.0
}
ERROR_CORRECTION_EFFICIENCY = 1.16 # Typical efficiency factor