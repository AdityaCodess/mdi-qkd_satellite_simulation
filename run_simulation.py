# run_simulation.py
import numpy as np
from qkd_simulation.channels import GroundToSatelliteLink, InterSatelliteLink
from qkd_simulation.protocol import MDIQKDProtocol

def main():
    """Main function to run the QKD simulation."""
    
    # --- Initialize the Protocol ---
    mdi_protocol = MDIQKDProtocol()
    
    # --- Simulation 1: Inter-Satellite Link (ISL) Performance vs. Range ---
    print("--- Running Inter-Satellite Link (ISL) Simulation ---")
    print("Range (km) | Loss (dB) | QBER (%)  | SKR (bits/pulse)")
    print("-" * 50)
    
    # Define a range of distances to test
    isl_ranges_km = np.linspace(1000, 7000, 7)
    
    for range_km in isl_ranges_km:
        # Create an ISL channel instance for this range
        isl_channel = InterSatelliteLink(range_km=range_km, divergence_mrad=0.005)
        
        # Get the channel's physical properties
        channel_props = isl_channel.get_channel_properties()
        
        # Calculate the SKR for MDI-QKD over this channel
        skr = mdi_protocol.calculate_skr(channel_props)
        
        # Print the results
        loss = channel_props['total_loss_db']
        qber = channel_props['effective_qber'] * 100
        print(f"{range_km:<10.0f} | {loss:<9.2f} | {qber:<9.2f} | {skr:e}")
        
    # TODO: Add logic here to save results and call the visualizer.

    # --- Simulation 2: Ground-to-Satellite Link Performance vs. Turbulence ---
    # TODO: Create a similar loop for the GroundToSatelliteLink,
    # varying the turbulence_strength parameter.
    
if __name__ == "__main__":
    main()