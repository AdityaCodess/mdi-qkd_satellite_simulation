import numpy as np
from scipy.integrate import quad
try:
    from . import parameters as p 
except ImportError:
    import parameters as p 

# --- Atmospheric Extinction Coefficient Model ---
def extinction_coefficient(h_meters, wavelength_m):
    """
    Calculates atmospheric extinction coefficient alpha(h) in m^-1.
    Uses a simplified exponential decay model based on sea-level visibility V_km.
    Ref: e.g., Free Space Optical Communication by H. Kaushal et al. (Ch 2/3) or similar FSO texts.
    Note: A more precise model would use MODTRAN or detailed aerosol/molecular scattering data.

    Args:
        h_meters (float): Altitude in meters.
        wavelength_m (float): Wavelength in meters.

    Returns:
        float: Extinction coefficient alpha in m^-1.
    """
    if h_meters < 0: h_meters = 0
    # Assume a standard sea-level visibility, e.g., V_km = 23 km (clear air)
    V_km = 23.0
    lambda_nm = wavelength_m * 1e9

    # Kruse model approximation for scattering coefficient beta_scat
    if V_km > 50: q = 1.6
    elif V_km > 6: q = 1.3
    else: q = 0.16 * V_km + 0.34

    beta_scat_sea_level = (3.912 / V_km) * (lambda_nm / 550.0)**(-q) # in km^-1

    # Scale height for aerosols/molecules (simplified, often ~1.2 km for aerosols)
    scale_height_aerosol = 1200 # meters
    # Convert beta_scat to m^-1 and apply exponential decay with altitude
    alpha_h = (beta_scat_sea_level / 1000.0) * np.exp(-h_meters / scale_height_aerosol)

    # Add a small baseline absorption (can be refined)
    alpha_absorption = 1e-6 # m^-1 (very low at typical QKD wavelengths)
    
    total_alpha = alpha_h + alpha_absorption
    return total_alpha if np.isfinite(total_alpha) else 0.0


# --- Hufnagel-Valley C_n^2 Model ---
def hufnagel_valley(h_meters, A, v_wind_ms):
    if h_meters < 0: h_meters = 0
    try:
        term1 = 5.94e-53 * (v_wind_ms / 27)**2 * (h_meters**10) * np.exp(-h_meters / 1000)
        term2 = 2.7e-16 * np.exp(-h_meters / 1500)
        term3 = A * np.exp(-h_meters / 100)
        result = term1 + term2 + term3
        return result if np.isfinite(result) else 0.0
    except (OverflowError, ValueError):
        return 0.0

TURBULENCE_PARAMS = {
    'low': {'A': 1.7e-14, 'v_wind': 21}, 'medium': {'A': 2.75e-14, 'v_wind': 21}, 'high': {'A': 2.75e-14, 'v_wind': 57}
}
R_EARTH = 6371e3; H_ATMOS = 20000

class BaseChannel:
    """Base class requiring wavelength."""
    def __init__(self, range_km, divergence_mrad, wavelength_nm):
        self.range_m = range_km * 1000
        self.divergence_rad = divergence_mrad * 1e-3
        self.wavelength_m = wavelength_nm * 1e-9
        self.satellite_h_meters = 500e3
        self.initial_w0_radius = self.wavelength_m / (np.pi * self.divergence_rad) if self.divergence_rad > 0 else 0

    def calculate_geometric_loss_db(self, effective_beam_waist_w_diameter):
        if self.range_m <= 0 or effective_beam_waist_w_diameter <= 0: return np.inf
        rx_area = np.pi * (p.RX_APERTURE_DIAMETER / 2)**2
        beam_area_at_rx = np.pi * (effective_beam_waist_w_diameter / 2)**2
        ratio = rx_area / beam_area_at_rx if beam_area_at_rx > 0 else 0
        return -10 * np.log10(ratio) if ratio > 0 else np.inf

    def get_zenith_angle(self):
        h = self.satellite_h_meters; z = self.range_m
        if z <= 0 or R_EARTH <= 0 or 2 * z * R_EARTH == 0: return np.pi / 2
        cos_theta_arg = h/z + (h**2 - z**2) / (2 * z * R_EARTH)
        cos_theta_arg = np.clip(cos_theta_arg, -1.0, 1.0)
        return np.arccos(cos_theta_arg)

    def get_altitude_along_path(self, path_distance_y, zenith_angle_theta):
        if path_distance_y < 0: path_distance_y = 0
        cos_theta = np.cos(zenith_angle_theta) if np.isfinite(zenith_angle_theta) else 0.0
        term_inside_sqrt = R_EARTH**2 + path_distance_y**2 + 2*path_distance_y*R_EARTH*cos_theta
        if term_inside_sqrt < 0: term_inside_sqrt = 0
        return np.sqrt(term_inside_sqrt) - R_EARTH


class GroundToSatelliteLink(BaseChannel):
    """GTS channel with refined atmospheric loss."""
    def __init__(self, range_km, divergence_mrad, wavelength_nm, turbulence_strength='medium'):
        super().__init__(range_km, divergence_mrad, wavelength_nm)
        self.turbulence_strength_label = turbulence_strength
        self.hv_params = TURBULENCE_PARAMS.get(turbulence_strength, TURBULENCE_PARAMS['medium'])

    # (calculate_coherence_length_rho0_uplink remains the same)
    def calculate_coherence_length_rho0_uplink(self, zenith_angle):
        k = 2 * np.pi / self.wavelength_m
        def integrand(xi):
            xi = np.clip(xi, 0, self.range_m)
            h_at_xi = self.get_altitude_along_path(xi, zenith_angle)
            cn2_at_h = hufnagel_valley(h_at_xi, self.hv_params['A'], self.hv_params['v_wind'])
            factor = (1 - xi / self.range_m)**(5/3) if self.range_m > 0 else 0
            if factor < 0: factor = 0
            res = factor * cn2_at_h
            return res if np.isfinite(res) else 0.0
        integral_limit = self.range_m if self.range_m > 0 else 0
        try: integral_I0_up, error = quad(integrand, 0, integral_limit, limit=100, epsabs=1e-9, epsrel=1e-9)
        except Exception: integral_I0_up = 0.0
        denominator = 1.46 * k**2 * integral_I0_up
        if denominator <= 0: return 0
        rho0 = denominator**(-3/5)
        return rho0

    # (calculate_rytor_variance_uplink remains the same)
    def calculate_rytor_variance_uplink(self, zenith_angle):
        k = 2 * np.pi / self.wavelength_m
        def local_rytor_integrand(xi):
            if xi <= 0 or self.range_m <= 0 or xi >= self.range_m: return 0.0
            h_at_xi = self.get_altitude_along_path(xi, zenith_angle)
            cn2_at_h = hufnagel_valley(h_at_xi, self.hv_params['A'], self.hv_params['v_wind'])
            factor = (xi / self.range_m) * (1 - xi / self.range_m)
            if factor < 0: factor = 0.0
            integrand_value = 0.0
            try:
                if factor >= 0: integrand_value = cn2_at_h * factor**(5/6)
            except ValueError: integrand_value = 0.0
            return integrand_value if np.isfinite(integrand_value) else 0.0
        integral_limit = self.range_m if self.range_m > 0 else 0
        rytor_integral_val = 0.0
        try: rytor_integral_val, error = quad(local_rytor_integrand, 0, integral_limit, limit=100, epsabs=1e-9, epsrel=1e-9)
        except Exception: rytor_integral_val = 0.0
        sigma_R_squared = 2.25 * k**(7/6) * rytor_integral_val
        return sigma_R_squared if np.isfinite(sigma_R_squared) and sigma_R_squared >= 0 else 0.0

    # (calculate_turbulence_effects_uplink remains the same)
    def calculate_turbulence_effects_uplink(self, rho0, zenith_angle):
        w0_rad = self.initial_w0_radius; z = self.range_m; lambda_ = self.wavelength_m
        w_d_diameter = 0
        if w0_rad > 0 and lambda_ > 0 and z > 0:
             zR = np.pi * w0_rad**2 / lambda_
             w_d_diameter = 2 * w0_rad * np.sqrt(1 + (z/zR)**2) if zR > 0 else 2*w0_rad
        if rho0 <= 0 or w0_rad <= 0 or z <= 0:
            return {'sigma_TB_sq': 0, 'w_st': w_d_diameter if w_d_diameter > 0 else 4*w0_rad}

        def cn2_integrand(h): return hufnagel_valley(h, self.hv_params['A'], self.hv_params['v_wind'])
        I_infinity, _ = quad(cn2_integrand, 0, 50000)
        a = 26.28 * I_infinity**(6/5); c = 7.71 * I_infinity
        cos_zenith = np.cos(zenith_angle)
        sec_theta = 1 / cos_zenith if abs(cos_zenith) > 1e-9 else 1.0 / 1e-9
        if not np.isfinite(sec_theta) or sec_theta > 10: sec_theta = 10.0

        w0_diameter = w0_rad * 2
        sigma_TB_sq = c * w0_diameter**(-1/3) * z**2 * sec_theta if w0_diameter > 0 else 0
        delta_theta = a * lambda_**(-2/5) * sec_theta**(6/5) - (c * w0_diameter**(-1/3) * sec_theta if w0_diameter > 0 else 0)
        w_st_sq_term = w_d_diameter**2 + z**2 * delta_theta
        if w_st_sq_term < 0: w_st_sq_term = w_d_diameter**2
        w_st_diameter = np.sqrt(w_st_sq_term)
        return {'sigma_TB_sq': sigma_TB_sq, 'w_st': w_st_diameter}

    # --- UPDATED Atmospheric Loss Calculation ---
    def calculate_atmospheric_loss_db(self, zenith_angle):
        """Calculates atmospheric loss by integrating extinction coeff along slant path."""
        def integrand(y): # y is distance along path
            h_at_y = self.get_altitude_along_path(y, zenith_angle)
            # Stop integrating if altitude is above effective atmosphere
            if h_at_y > H_ATMOS: return 0.0
            alpha = extinction_coefficient(h_at_y, self.wavelength_m)
            return alpha if np.isfinite(alpha) else 0.0

        # Integrate only up to the point where altitude reaches H_ATMOS,
        # or the full range if satellite is below H_ATMOS
        integration_range = self.range_m
        if self.satellite_h_meters > H_ATMOS:
             # Find range 'z_atmos' where altitude = H_ATMOS (approximate using secant)
             z_atmos_approx = H_ATMOS / np.cos(zenith_angle) if np.cos(zenith_angle) > 1e-9 else H_ATMOS * 10
             integration_range = min(self.range_m, z_atmos_approx)

        if integration_range <= 0: return 0.0 # No path through atmosphere

        try:
             # Calculate optical depth (tau)
             optical_depth, error = quad(integrand, 0, integration_range, limit=100, epsabs=1e-9, epsrel=1e-9)
             # Transmittance = exp(-tau)
             transmittance = np.exp(-optical_depth)
             # Convert transmittance to dB loss
             loss_db = -10 * np.log10(transmittance) if transmittance > 0 else np.inf
             return loss_db if np.isfinite(loss_db) else 200.0 # Return high loss if inf
        except Exception:
             return 200.0 # Return high loss if integration fails

    def get_channel_properties(self, intrinsic_qber, qber_rytor_scaling):
        """Calculates total loss and effective QBER including refined turbulence & atmos loss."""
        zenith_angle = self.get_zenith_angle()
        rho0 = self.calculate_coherence_length_rho0_uplink(zenith_angle)
        turb_effects = self.calculate_turbulence_effects_uplink(rho0, zenith_angle)
        sigma_TB_sq = turb_effects['sigma_TB_sq']
        w_st_diameter = turb_effects['w_st']

        geometric_loss_db = self.calculate_geometric_loss_db(effective_beam_waist_w_diameter=w_st_diameter)
        beam_wander_loss_db = 0
        if w_st_diameter > 0 :
             sigma_TB_sq_safe = max(0, sigma_TB_sq)
             beam_wander_loss_db = 4.34 * (sigma_TB_sq_safe / (w_st_diameter/2)**2)

        # Use the new atmospheric loss calculation
        atmospheric_loss_db = self.calculate_atmospheric_loss_db(zenith_angle)

        channel_loss_db = geometric_loss_db + atmospheric_loss_db + beam_wander_loss_db

        rytor_variance = self.calculate_rytor_variance_uplink(zenith_angle)
        qber_from_turbulence = qber_rytor_scaling * rytor_variance
        effective_qber = intrinsic_qber + qber_from_turbulence
        effective_qber = min(effective_qber, 0.5)

        return {'channel_loss_db': channel_loss_db if np.isfinite(channel_loss_db) else 200.0,
                'effective_qber': effective_qber}


# (InterSatelliteLink class remains the same)
class InterSatelliteLink(BaseChannel):
    """ISL channel, requires wavelength."""
    def __init__(self, range_km, divergence_mrad, wavelength_nm, pointing_error_urad=1.0):
        super().__init__(range_km, divergence_mrad, wavelength_nm)
        self.pointing_error_urad = pointing_error_urad

    def calculate_pointing_loss_db(self):
        if self.range_m <= 0 or self.divergence_rad <= 0: return np.inf
        w0_rad = self.initial_w0_radius; lambda_ = self.wavelength_m; z = self.range_m
        if w0_rad <= 0 or lambda_ <= 0: return np.inf
        zR = np.pi * w0_rad**2 / lambda_
        w_z_radius = w0_rad * np.sqrt(1 + (z/zR)**2) if zR > 0 else w0_rad
        pointing_error_rad = self.pointing_error_urad * 1e-6
        jitter_std_dev_at_rx = pointing_error_rad * self.range_m
        if w_z_radius <= 0: return np.inf
        ratio_sq = (jitter_std_dev_at_rx / w_z_radius)**2
        loss_db = 4.34 * ratio_sq if np.isfinite(ratio_sq) else np.inf
        return loss_db

    def get_channel_properties(self, intrinsic_qber):
        lambda_ = self.wavelength_m; w0_rad = self.initial_w0_radius; z = self.range_m
        w_d_diameter = 0
        if w0_rad > 0 and lambda_ > 0 and z > 0:
            zR = np.pi * w0_rad**2 / lambda_
            term_in_sqrt = 1 + (z/zR)**2 if zR > 0 else 1
            w_d_diameter = 2 * w0_rad * np.sqrt(term_in_sqrt)

        geometric_loss_db = self.calculate_geometric_loss_db(effective_beam_waist_w_diameter=w_d_diameter)
        pointing_loss_db = self.calculate_pointing_loss_db()
        channel_loss_db = geometric_loss_db + pointing_loss_db
        effective_qber = intrinsic_qber
        return {'channel_loss_db': channel_loss_db if np.isfinite(channel_loss_db) else 200.0,
                'effective_qber': effective_qber}