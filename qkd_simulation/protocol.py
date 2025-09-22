# qkd_simulation/protocol.py
import numpy as np
from . import parameters as p

def binary_entropy(p_err):
    """Calculates the binary Shannon entropy H(p)."""
    if p_err <= 0 or p_err >= 1:
        return 0
    return -p_err * np.log2(p_err) - (1 - p_err) * np.log2(1 - p_err)

class MDIQKDProtocol:
    """Implements the MDI-QKD protocol logic and SKR calculation."""
    def __init__(self):
        self.signal_intensity = p.DECOY_INTENSITIES['signal']

    def calculate_skr(self, channel_properties):
        """
        Calculates the Secure Key Rate (SKR) for MDI-QKD.
        Ref: Based on the security analysis for MDI-QKD
        """
        loss_db = channel_properties['total_loss_db']
        qber_channel = channel_properties['effective_qber']
        
        # Total channel transmittance from one user to the central station
        eta = 10**(-loss_db / 10)
        
        # --- MDI-QKD SKR Calculation ---
        # NOTE: A full decoy-state analysis is complex. This is a robust simulation model
        # that reflects the expected behavior of the protocol.

        # 1. Background event probability at the central station
        background_prob = 1 - (1 - p.DARK_COUNT_RATE / p.REPETITION_RATE)**2
        
        # 2. Approximate overall gain (Q_z) and error rate (E_z) for the signal state.
        # This is the probability of a coincidence detection event at Charlie.
        q_coincidence = (self.signal_intensity**2 * np.exp(-2*self.signal_intensity)) * (eta**2) * p.DETECTOR_EFFICIENCY**2
        q_total_signal = 1 - (1 - q_coincidence) * (1 - background_prob)
        
        e_total_signal = (qber_channel * q_coincidence + 0.5 * background_prob) / q_total_signal if q_total_signal > 0 else 0.5
        
        # 3. Approximate the single-photon pair yield (Y_11) and phase error rate (e_11_ph)
        # This is the core result from the decoy-state method.
        y_11 = q_coincidence # Simplified for simulation
        e_11_ph = qber_channel # In MDI, phase error is closely tied to the channel's bit-flip error

        # 4. Final SKR formula
        # R = P_signal^2 * [ Y_11 * (1 - H(e_phase)) - Q_signal * f * H(E_bit) ]
        # The (0.5 * P_signal)^2 factor accounts for signal state probability and basis choices
        skr = (0.5 * self.signal_intensity)**2 * (y_11 * (1 - binary_entropy(e_11_ph)) - \
              q_total_signal * p.ERROR_CORRECTION_EFFICIENCY * binary_entropy(e_total_signal))
        
        return max(0, skr)