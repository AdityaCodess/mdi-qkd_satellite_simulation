# gui.py
import tkinter as tk
from tkinter import ttk
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from qkd_simulation.channels import GroundToSatelliteLink, InterSatelliteLink
from qkd_simulation.protocol import DecoyBB84Protocol, MDIQKDProtocol

class SimulationGUI:
    def __init__(self, master):
        self.master = master
        master.title("QKD Satellite Simulation")
        master.geometry("1100x800")
        self.create_widgets()

    def create_widgets(self):
        # --- Input Frame ---
        input_frame = ttk.LabelFrame(self.master, text="Input Parameters")
        input_frame.grid(row=0, column=0, padx=10, pady=10, sticky="ns")

        # Protocol Selection
        ttk.Label(input_frame, text="Protocol:").grid(row=0, column=0, sticky="w", padx=5, pady=5)
        self.protocol_type = tk.StringVar(value="Decoy-State BB84")
        protocol_options = ["Decoy-State BB84", "MDI-QKD"]
        self.protocol_combo = ttk.Combobox(input_frame, textvariable=self.protocol_type, values=protocol_options, state="readonly")
        self.protocol_combo.grid(row=1, column=0, padx=5, pady=2, sticky="ew")

        ttk.Separator(input_frame, orient='horizontal').grid(row=2, column=0, sticky='ew', pady=10)

        # Graph Type Selection
        ttk.Label(input_frame, text="Graph Type:").grid(row=3, column=0, sticky="w", padx=5, pady=5)
        self.graph_type = tk.StringVar(value="SKR vs. Range")
        graph_options = [
            "SKR vs. Range", "Total Loss vs. Range", "QBER vs. Range", 
            "Received Power vs. Tx Power", "QBER vs. Total Loss",
            "SKR vs. Pointing Error (ISL)", "SKR vs. Turbulence (GTS)",
            "SKR vs. QBER", "Key Rate Components vs. Range"
        ]
        self.graph_combo = ttk.Combobox(input_frame, textvariable=self.graph_type, values=graph_options, state="readonly")
        self.graph_combo.grid(row=4, column=0, padx=5, pady=2, sticky="ew")
        self.graph_combo.bind("<<ComboboxSelected>>", self.toggle_inputs)

        ttk.Separator(input_frame, orient='horizontal').grid(row=5, column=0, sticky='ew', pady=10)

        # Link Type Selection
        self.link_type = tk.StringVar(value="ISL")
        ttk.Label(input_frame, text="Link Type:").grid(row=6, column=0, sticky="w", padx=5, pady=2)
        ttk.Radiobutton(input_frame, text="Inter-Satellite (ISL)", variable=self.link_type, value="ISL", command=self.toggle_inputs).grid(row=7, column=0, sticky="w", padx=10)
        ttk.Radiobutton(input_frame, text="Ground-to-Satellite", variable=self.link_type, value="GTS", command=self.toggle_inputs).grid(row=8, column=0, sticky="w", padx=10)

        # Numerical Inputs
        self.range_label = ttk.Label(input_frame, text="Range Parameter (km):")
        self.range_label.grid(row=9, column=0, sticky="w", padx=5, pady=2)
        self.range_entry = ttk.Entry(input_frame); self.range_entry.grid(row=10, column=0, padx=5, pady=2); self.range_entry.insert(0, "4000")
        
        ttk.Label(input_frame, text="Divergence (mrad):").grid(row=11, column=0, sticky="w", padx=5, pady=2)
        self.divergence_entry = ttk.Entry(input_frame); self.divergence_entry.grid(row=12, column=0, padx=5, pady=2); self.divergence_entry.insert(0, "0.005")
        
        self.turbulence_label = ttk.Label(input_frame, text="Turbulence:")
        self.turbulence_label.grid(row=13, column=0, sticky="w", padx=5, pady=2)
        self.turbulence_combo = ttk.Combobox(input_frame, values=['low', 'medium', 'high'], state="readonly"); 
        self.turbulence_combo.grid(row=14, column=0, padx=5, pady=2); self.turbulence_combo.set('medium')
        
        self.pointing_label = ttk.Label(input_frame, text="Pointing Error (µrad):")
        self.pointing_label.grid(row=15, column=0, sticky="w", padx=5, pady=2)
        self.pointing_entry = ttk.Entry(input_frame); self.pointing_entry.grid(row=16, column=0, padx=5, pady=2); self.pointing_entry.insert(0, "1.0")
        
        ttk.Label(input_frame, text="Intrinsic QBER (%):").grid(row=17, column=0, sticky="w", padx=5, pady=2)
        self.qber_entry = ttk.Entry(input_frame); self.qber_entry.grid(row=18, column=0, padx=5, pady=2); self.qber_entry.insert(0, "1.0")

        # --- NEW WIDGETS FOR TRANSMITTED POWER ---
        self.tx_power_min_label = ttk.Label(input_frame, text="Min Tx Power (dBm):")
        self.tx_power_min_label.grid(row=19, column=0, sticky="w", padx=5, pady=2)
        self.tx_power_min_entry = ttk.Entry(input_frame)
        self.tx_power_min_entry.grid(row=20, column=0, padx=5, pady=2)
        self.tx_power_min_entry.insert(0, "-20")

        self.tx_power_max_label = ttk.Label(input_frame, text="Max Tx Power (dBm):")
        self.tx_power_max_label.grid(row=21, column=0, sticky="w", padx=5, pady=2)
        self.tx_power_max_entry = ttk.Entry(input_frame)
        self.tx_power_max_entry.grid(row=22, column=0, padx=5, pady=2)
        self.tx_power_max_entry.insert(0, "10")

        # Plot Frame and Controls
        plot_frame = ttk.LabelFrame(self.master, text="Simulation Results")
        plot_frame.grid(row=0, column=1, padx=10, pady=10, sticky="nsew")
        
        self.fig = Figure(figsize=(8, 6), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        control_frame = ttk.Frame(self.master)
        control_frame.grid(row=1, column=0, columnspan=2, pady=10)
        
        run_button = ttk.Button(control_frame, text="Run Simulation", command=self.run_simulation)
        run_button.pack()
        
        self.status_label = ttk.Label(self.master, text="Status: Ready")
        self.status_label.grid(row=2, column=0, columnspan=2, pady=5)

        self.master.grid_columnconfigure(1, weight=1); self.master.grid_rowconfigure(0, weight=1)
        self.toggle_inputs()

    def toggle_inputs(self, event=None):
        is_isl = self.link_type.get() == "ISL"
        graph_choice = self.graph_type.get()
        
        # Toggle link-specific inputs
        self.turbulence_label.config(state='disabled' if is_isl else 'normal')
        self.turbulence_combo.config(state='disabled' if is_isl else 'normal')
        self.pointing_label.config(state='normal' if is_isl else 'disabled')
        self.pointing_entry.config(state='normal' if is_isl else 'disabled')
        
        # Toggle graph-specific inputs
        is_power_graph = graph_choice == "Received Power vs. Tx Power"
        self.tx_power_min_label.config(state='normal' if is_power_graph else 'disabled')
        self.tx_power_min_entry.config(state='normal' if is_power_graph else 'disabled')
        self.tx_power_max_label.config(state='normal' if is_power_graph else 'disabled')
        self.tx_power_max_entry.config(state='normal' if is_power_graph else 'disabled')

        if graph_choice in ["SKR vs. Pointing Error (ISL)", "SKR vs. Turbulence (GTS)", "SKR vs. QBER", "Received Power vs. Tx Power"]:
             self.range_label.config(text="Fixed Range (km):")
        else:
             self.range_label.config(text="Max Range (km):")

    def update_plot(self, x_data, y_data, title, xlabel, ylabel, yscale='linear', **kwargs):
        self.ax.clear()
        if kwargs.get('is_bar'):
            self.ax.bar(x_data, y_data)
        else:
            self.ax.plot(x_data, y_data, 'o-', label=kwargs.get('label1'))
            if 'y2' in kwargs:
                self.ax.plot(x_data, kwargs['y2'], 's--', label=kwargs.get('label2'))
                self.ax.legend()
        self.ax.set_title(title); self.ax.set_xlabel(xlabel); self.ax.set_ylabel(ylabel)
        self.ax.grid(True)
        if yscale == 'log' and y_data and any(y > 0 for y in y_data):
            self.ax.set_yscale('log')
        else:
            self.ax.set_yscale('linear')
        self.fig.tight_layout()
        self.canvas.draw()

    def run_simulation(self):
        self.status_label.config(text="Status: Running..."); self.master.update_idletasks()

        protocol_choice = self.protocol_type.get()
        link_type = self.link_type.get()
        range_param = float(self.range_entry.get())
        divergence = float(self.divergence_entry.get())
        graph_choice = self.graph_type.get()
        intrinsic_qber = float(self.qber_entry.get()) / 100.0

        protocol = DecoyBB84Protocol() if protocol_choice == "Decoy-State BB84" else MDIQKDProtocol()
        
        if graph_choice == "Received Power vs. Tx Power":
            channel = InterSatelliteLink(range_km=range_param, divergence_mrad=divergence, pointing_error_urad=float(self.pointing_entry.get())) if link_type == 'ISL' else \
                      GroundToSatelliteLink(range_km=range_param, divergence_mrad=divergence, turbulence_strength=self.turbulence_combo.get())
            props = channel.get_channel_properties(intrinsic_qber=intrinsic_qber)
            total_loss_db = props['total_loss_db']
            
            # --- USE VALUES FROM NEW UI WIDGETS ---
            min_tx = float(self.tx_power_min_entry.get())
            max_tx = float(self.tx_power_max_entry.get())
            tx_powers_dbm = np.linspace(min_tx, max_tx, 30)
            
            rx_powers_dbm = tx_powers_dbm - total_loss_db
            
            self.update_plot(tx_powers_dbm, rx_powers_dbm, f"{link_type} Link Budget at {range_param} km", "Transmitted Power (dBm)", "Received Power (dBm)")
            self.status_label.config(text=f"Status: Complete. Total link loss: {total_loss_db:.2f} dB")
            return

        # --- The rest of the function remains the same ---
        # It correctly handles all the other graph types.
        # ... (full logic from previous response)
        ranges_km = np.linspace(max(100, range_param/10), range_param, 30)
        results = {'loss': [], 'qber': [], 'skr': [], 'useful': [], 'leakage': []}

        for r in ranges_km:
            if link_type == 'ISL':
                channel = InterSatelliteLink(range_km=r, divergence_mrad=divergence, pointing_error_urad=float(self.pointing_entry.get()))
            else:
                channel = GroundToSatelliteLink(range_km=r, divergence_mrad=divergence, turbulence_strength=self.turbulence_combo.get())
            
            props = channel.get_channel_properties(intrinsic_qber=intrinsic_qber)
            res = protocol.calculate_skr(props)
            results['loss'].append(props['total_loss_db']); results['qber'].append(props['effective_qber'])
            results['skr'].append(res['skr']); results['useful'].append(res['useful_rate']); results['leakage'].append(res['leakage_rate'])

        if graph_choice == "SKR vs. Range":
            self.update_plot(ranges_km, results['skr'], f"{protocol_choice} on {link_type}", "Range (km)", "SKR (bits/pulse)", yscale='log')
            self.status_label.config(text=f"Status: Complete. Final SKR at {range_param} km: {results['skr'][-1]:.2e} bits/pulse")
        
        elif graph_choice == "Total Loss vs. Range":
            self.update_plot(ranges_km, results['loss'], f"{link_type} Performance", "Range (km)", "Total Loss (dB)")
            self.status_label.config(text=f"Status: Complete. Final Loss at {range_param} km: {results['loss'][-1]:.2f} dB")
        
        elif graph_choice == "QBER vs. Range":
            self.update_plot(ranges_km, results['qber'], f"{link_type} Performance", "Range (km)", "Effective QBER", yscale='linear')
            self.status_label.config(text=f"Status: Complete. Final QBER at {range_param} km: {results['qber'][-1]:.2%}")
            
        elif graph_choice == "QBER vs. Total Loss":
            self.update_plot(results['loss'], results['qber'], f"{link_type} Performance", "Total Loss (dB)", "Effective QBER")
            self.status_label.config(text=f"Status: Complete. Final QBER: {results['qber'][-1]:.2%}")

        elif graph_choice == "Key Rate Components vs. Range":
            self.update_plot(ranges_km, results['useful'], f"{protocol_choice} Components on {link_type}", "Range (km)", "Rate (bits/pulse)", yscale='log', 
                             y2=results['leakage'], label1='Useful Rate', label2='Leakage Rate')
            self.status_label.config(text=f"Status: Complete. Analysis of key rate components.")

        elif graph_choice == "SKR vs. Pointing Error (ISL)":
            pointing_errors = np.linspace(0.5, 5.0, 30)
            skr_results = []
            for pe in pointing_errors:
                channel = InterSatelliteLink(range_km=range_param, divergence_mrad=divergence, pointing_error_urad=pe)
                props = channel.get_channel_properties(intrinsic_qber=intrinsic_qber)
                skr_results.append(protocol.calculate_skr(props)['skr'])
            self.update_plot(pointing_errors, skr_results, f"{protocol_choice} at {range_param}km", "Pointing Error (µrad)", "SKR (bits/pulse)", yscale='log')
            self.status_label.config(text=f"Status: Complete. Analysis of pointing error sensitivity.")

        elif graph_choice == "SKR vs. Turbulence (GTS)":
            turbulence_levels = ['low', 'medium', 'high']
            skr_results = []
            for tl in turbulence_levels:
                channel = GroundToSatelliteLink(range_km=range_param, divergence_mrad=divergence, turbulence_strength=tl)
                props = channel.get_channel_properties(intrinsic_qber=intrinsic_qber)
                skr_results.append(protocol.calculate_skr(props)['skr'])
            self.update_plot(turbulence_levels, skr_results, f"{protocol_choice} at {range_param}km", "Turbulence Strength", "SKR (bits/pulse)", is_bar=True, yscale='log')
            self.status_label.config(text=f"Status: Complete. Analysis of turbulence resilience.")
            
        elif graph_choice == "SKR vs. QBER":
            qber_values = np.linspace(0, 0.15, 30)
            skr_results = []
            channel = InterSatelliteLink(range_km=range_param, divergence_mrad=divergence, pointing_error_urad=float(self.pointing_entry.get())) if link_type == 'ISL' else \
                      GroundToSatelliteLink(range_km=range_param, divergence_mrad=divergence, turbulence_strength=self.turbulence_combo.get())
            base_props = channel.get_channel_properties(intrinsic_qber=0) 
            
            for qber in qber_values:
                current_props = {'total_loss_db': base_props['total_loss_db'], 'effective_qber': qber}
                skr_results.append(protocol.calculate_skr(current_props)['skr'])
            self.update_plot(qber_values * 100, skr_results, f"{protocol_choice} at {range_param}km", "Total Effective QBER (%)", "SKR (bits/pulse)", yscale='log')
            self.status_label.config(text=f"Status: Complete. Analysis of error tolerance.")


if __name__ == "__main__":
    root = tk.Tk()
    app = SimulationGUI(root)
    root.mainloop()