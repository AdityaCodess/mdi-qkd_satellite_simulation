# qkd_simulation/protocol.py
import numpy as np
import math
try:
    from . import parameters as p # Relative import
except ImportError:
    import parameters as p # Absolute import

def binary_entropy(p_err):
    """Calculates the binary Shannon entropy H(p)."""
    if p_err <= 0 or p_err >= 1: return 0
    return -p_err * np.log2(p_err) - (1 - p_err) * np.log2(1 - p_err)

def poisson(mu, n):
    """Returns the probability of n photons in a pulse with mean mu."""
    if n > 20 or n < 0: return 0 # Ensure n is valid
    # Use try-except for factorial calculation
    try:
        fact = math.factorial(int(n)) # Ensure n is integer
        if fact == 0: return 0 # Avoid division by zero
        # Handle potential overflow with large mu or n
        term_mu_n = mu**n
        if not np.isfinite(term_mu_n): return 0.0
        exp_neg_mu = np.exp(-mu)
        if not np.isfinite(exp_neg_mu): return 0.0

        prob = exp_neg_mu * term_mu_n / fact
        return prob if np.isfinite(prob) else 0.0
    except (ValueError, OverflowError):
        return 0.0

class DecoyBB84Protocol:
    """Implements decoy-state BB84."""
    def __init__(self):
        self.signal_intensity = p.DECOY_INTENSITIES['signal']
        self.protocol_efficiency = 0.5

    # --- UPDATED: calculate_skr takes hardware efficiencies ---
    def calculate_skr(self, channel_properties, detector_efficiency, receiver_optics_loss_db):
        """Calculates the SKR for decoy-state BB84."""
        channel_loss_db = channel_properties['channel_loss_db']
        e_effective_channel = channel_properties['effective_qber'] # Includes intrinsic + turbulence

        # Calculate total loss and transmittance
        total_loss_db = channel_loss_db + receiver_optics_loss_db
        eta_total_link = 10**(-total_loss_db / 10) # Transmittance excluding detector
        eta_detector = detector_efficiency
        eta_end_to_end = eta_total_link * eta_detector # Overall efficiency

        y0 = p.DARK_COUNT_RATE / p.REPETITION_RATE

        # Yield & Error Rate of single photons (using overall efficiency)
        y1 = 1 - (1 - y0) * (1 - eta_end_to_end)**1
        # Error rate depends on channel QBER + background
        e1 = (0.5 * y0 + e_effective_channel * (y1 - y0)) / y1 if y1 > 0 else 0.5

        # Overall Gain (Q_mu) & QBER (E_mu)
        q_mu = 0.0; e_mu_numerator = 0.0
        for n in range(5):
            p_n = poisson(self.signal_intensity, n)
            # Yield yn depends on overall efficiency eta_end_to_end
            y_n = 1 - (1 - y0) * (1 - eta_end_to_end)**n
            # Error en depends on channel QBER + background
            e_n = (0.5 * y0 + e_effective_channel * (y_n - y0)) / y_n if y_n > 0 else 0.5
            q_mu += p_n * y_n
            e_mu_numerator += p_n * y_n * e_n

        e_mu = e_mu_numerator / q_mu if q_mu > 0 else 0.5
        q1 = y1 * poisson(self.signal_intensity, 1)

        leakage_rate = q_mu * p.ERROR_CORRECTION_EFFICIENCY * binary_entropy(e_mu)
        useful_rate = q1 * (1 - binary_entropy(e1))
        skr_per_pulse = self.protocol_efficiency * (useful_rate - leakage_rate)

        return {
            'skr': max(0.0, skr_per_pulse) if np.isfinite(skr_per_pulse) else 0.0,
            'useful_rate': self.protocol_efficiency * useful_rate if np.isfinite(useful_rate) else 0.0,
            'leakage_rate': self.protocol_efficiency * leakage_rate if np.isfinite(leakage_rate) else 0.0
        }

class MDIQKDProtocol:
    """Implements MDI-QKD."""
    def __init__(self):
        self.signal_intensity = p.DECOY_INTENSITIES['signal']

    # --- UPDATED: calculate_skr takes hardware efficiencies ---
    def calculate_skr(self, channel_properties, detector_efficiency, receiver_optics_loss_db):
        """Calculates the Secure Key Rate (SKR) for MDI-QKD."""
        channel_loss_db = channel_properties['channel_loss_db']
        qber_channel = channel_properties['effective_qber'] # Includes intrinsic + turbulence

        # Calculate total loss and transmittance for ONE ARM (Alice/Bob to Charlie)
        total_loss_db_one_arm = channel_loss_db + receiver_optics_loss_db
        eta_link_one_arm = 10**(-total_loss_db_one_arm / 10) # Excluding detector
        eta_detector = detector_efficiency
        eta_end_to_end_one_arm = eta_link_one_arm * eta_detector # Overall efficiency for one arm

        background_prob = 1 - (1 - p.DARK_COUNT_RATE / p.REPETITION_RATE)**2

        # Coincidence prob depends on square of single-arm efficiency
        q_coincidence = (self.signal_intensity**2 * np.exp(-2*self.signal_intensity)) * (eta_end_to_end_one_arm**2)
        q_total_signal = 1 - (1 - q_coincidence) * (1 - background_prob)

        e_total_signal = (qber_channel * q_coincidence + 0.5 * background_prob) / q_total_signal if q_total_signal > 0 else 0.5

        y_11 = q_coincidence
        e_11_ph = qber_channel

        useful_rate = y_11 * (1 - binary_entropy(e_11_ph))
        leakage_rate = q_total_signal * p.ERROR_CORRECTION_EFFICIENCY * binary_entropy(e_total_signal)
        skr = (0.5 * self.signal_intensity)**2 * (useful_rate - leakage_rate)

        return {
            'skr': max(0.0, skr) if np.isfinite(skr) else 0.0,
            'useful_rate': (0.5 * self.signal_intensity)**2 * useful_rate if np.isfinite(useful_rate) else 0.0,
            'leakage_rate': (0.5 * self.signal_intensity)**2 * leakage_rate if np.isfinite(leakage_rate) else 0.0
        }