# qkd_simulation/channels.py
import numpy as np
from . import parameters as p

class BaseChannel:
    """Base class for all optical communication links."""
    def __init__(self, range_km, divergence_mrad):
        self.range_m = range_km * 1000
        self.divergence_rad = divergence_mrad * 1e-3

    def calculate_geometric_loss_db(self):
        """Calculates the loss from beam divergence."""
        beam_diameter_at_rx = self.divergence_rad * self.range_m
        if beam_diameter_at_rx == 0: return np.inf
        ratio = (p.RX_APERTURE_DIAMETER / beam_diameter_at_rx)**2
        return -10 * np.log10(ratio)

class GroundToSatelliteLink(BaseChannel):
    """Models the ground-to-satellite channel, including atmospheric effects."""
    def __init__(self, range_km, divergence_mrad, turbulence_strength='medium'):
        super().__init__(range_km, divergence_mrad)
        self.turbulence_strength = turbulence_strength

    def get_channel_properties(self, intrinsic_qber):
        """Calculates total loss and effective QBER for the atmospheric link."""
        total_loss = self.calculate_geometric_loss_db()
        total_loss += 2.0  # Atmospheric attenuation
        
        turbulence_losses = {'low': 1.0, 'medium': 3.0, 'high': 5.0}
        total_loss += turbulence_losses.get(self.turbulence_strength, 3.0)

        turbulence_qber = {'low': 0.005, 'medium': 0.015, 'high': 0.03}
        # The total QBER is the intrinsic hardware error PLUS the error from turbulence
        effective_qber = intrinsic_qber + turbulence_qber.get(self.turbulence_strength, 0.015)

        return {'total_loss_db': total_loss, 'effective_qber': effective_qber}

class InterSatelliteLink(BaseChannel):
    """Models the inter-satellite channel (vacuum)."""
    def __init__(self, range_km, divergence_mrad, pointing_error_urad=1.0):
        super().__init__(range_km, divergence_mrad)
        self.pointing_error_urad = pointing_error_urad

    def calculate_pointing_loss_db(self):
        """Calculates the average power loss from pointing errors (jitter)."""
        beam_radius_at_rx = (self.divergence_rad * self.range_m) / 2
        pointing_error_rad = self.pointing_error_urad * 1e-6 
        jitter_std_dev_at_rx = pointing_error_rad * self.range_m
        loss_db = 4.34 * (jitter_std_dev_at_rx / beam_radius_at_rx)**2
        return loss_db

    def get_channel_properties(self, intrinsic_qber):
        """Calculates total loss and effective QBER for the vacuum link."""
        total_loss = self.calculate_geometric_loss_db()
        total_loss += self.calculate_pointing_loss_db()
        
        # In vacuum, the total QBER is just the intrinsic hardware error
        effective_qber = intrinsic_qber

        return {'total_loss_db': total_loss, 'effective_qber': effective_qber}