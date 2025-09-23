# qkd_simulation/protocol.py
import numpy as np
import math
from . import parameters as p

def binary_entropy(p_err):
    """Calculates the binary Shannon entropy H(p)."""
    if p_err <= 0 or p_err >= 1:
        return 0
    return -p_err * np.log2(p_err) - (1 - p_err) * np.log2(1 - p_err)

def poisson(mu, n):
    """Returns the probability of n photons in a pulse with mean mu."""
    if n > 20: return 0
    return (np.exp(-mu) * mu**n) / math.factorial(n)

# --- PROTOCOL 1: DECOY-STATE BB84 ---
class DecoyBB84Protocol:
    """Implements the decoy-state BB84 protocol logic and SKR calculation."""
    def __init__(self):
        self.signal_intensity = p.DECOY_INTENSITIES['signal']
        self.protocol_efficiency = 0.5

    def calculate_skr(self, channel_properties):
        """Calculates the SKR for decoy-state BB84 and its components."""
        loss_db = channel_properties['total_loss_db']
        e_intrinsic = channel_properties['effective_qber']
        
        eta_channel = 10**(-loss_db / 10)
        eta_total = eta_channel * p.DETECTOR_EFFICIENCY
        
        y0 = p.DARK_COUNT_RATE / p.REPETITION_RATE
        
        y1 = 1 - (1 - y0) * (1 - eta_total)**1
        e1 = (0.5 * y0 + e_intrinsic * (y1 - y0)) / y1 if y1 > 0 else 0.5

        q_mu = 0; e_mu_numerator = 0
        for n in range(5):
            p_n = poisson(self.signal_intensity, n)
            y_n = 1 - (1 - y0) * (1 - eta_total)**n
            e_n = (0.5 * y0 + e_intrinsic * (y_n - y0)) / y_n if y_n > 0 else 0.5
            q_mu += p_n * y_n
            e_mu_numerator += p_n * y_n * e_n
        
        e_mu = e_mu_numerator / q_mu if q_mu > 0 else 0.5
        q1 = y1 * poisson(self.signal_intensity, 1)

        leakage_rate = q_mu * p.ERROR_CORRECTION_EFFICIENCY * binary_entropy(e_mu)
        useful_rate = q1 * (1 - binary_entropy(e1))
        
        skr_per_pulse = self.protocol_efficiency * (useful_rate - leakage_rate)
        
        return {
            'skr': max(0, skr_per_pulse),
            'useful_rate': self.protocol_efficiency * useful_rate,
            'leakage_rate': self.protocol_efficiency * leakage_rate
        }

# --- PROTOCOL 2: MDI-QKD ---
class MDIQKDProtocol:
    """Implements the MDI-QKD protocol logic and SKR calculation."""
    def __init__(self):
        self.signal_intensity = p.DECOY_INTENSITIES['signal']

    def calculate_skr(self, channel_properties):
        """Calculates the Secure Key Rate (SKR) for MDI-QKD and its components."""
        loss_db = channel_properties['total_loss_db']
        qber_channel = channel_properties['effective_qber']
        eta = 10**(-loss_db / 10)
        
        background_prob = 1 - (1 - p.DARK_COUNT_RATE / p.REPETITION_RATE)**2
        
        q_coincidence = (self.signal_intensity**2 * np.exp(-2*self.signal_intensity)) * (eta**2) * p.DETECTOR_EFFICIENCY**2
        q_total_signal = 1 - (1 - q_coincidence) * (1 - background_prob)
        
        e_total_signal = (qber_channel * q_coincidence + 0.5 * background_prob) / q_total_signal if q_total_signal > 0 else 0.5
        
        y_11 = q_coincidence
        e_11_ph = qber_channel

        useful_rate = y_11 * (1 - binary_entropy(e_11_ph))
        leakage_rate = q_total_signal * p.ERROR_CORRECTION_EFFICIENCY * binary_entropy(e_total_signal)
        
        skr = (0.5 * self.signal_intensity)**2 * (useful_rate - leakage_rate)
        
        return {
            'skr': max(0, skr),
            'useful_rate': (0.5 * self.signal_intensity)**2 * useful_rate,
            'leakage_rate': (0.5 * self.signal_intensity)**2 * leakage_rate
        }