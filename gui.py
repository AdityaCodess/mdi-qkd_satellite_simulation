# gui.py
import tkinter as tk
from tkinter import ttk
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# Import the simulation backend
from qkd_simulation.channels import GroundToSatelliteLink, InterSatelliteLink
# Import all three protocol classes
from qkd_simulation.protocol import DecoyBB84Protocol, MDIQKDProtocol, ClassicalProtocol

class SimulationGUI:
    def __init__(self, master):
        self.master = master
        master.title("QKD Satellite Simulation")
        master.geometry("1100x850") # Adjusted size
        self.create_widgets()

    def create_widgets(self):
        # --- Input Frame ---
        input_frame = ttk.LabelFrame(self.master, text="Input Parameters")
        input_frame.grid(row=0, column=0, padx=10, pady=10, sticky="ns")

        # --- Wavelength Input --- (Row 0, 1)
        ttk.Label(input_frame, text="Wavelength (nm):").grid(row=0, column=0, sticky="w", padx=5, pady=2)
        self.wavelength_entry = ttk.Entry(input_frame); self.wavelength_entry.grid(row=1, column=0, padx=5, pady=2, sticky="ew"); self.wavelength_entry.insert(0, "850")

        ttk.Separator(input_frame, orient='horizontal').grid(row=2, column=0, sticky='ew', pady=5)

        # --- Protocol Selection (UPDATED) --- (Row 3, 4)
        ttk.Label(input_frame, text="Protocol:").grid(row=3, column=0, sticky="w", padx=5, pady=2)
        self.protocol_type = tk.StringVar(value="Decoy-State BB84")
        # Add new option to the list
        protocol_options = ["Decoy-State BB84", "MDI-QKD", "Classical (Normal Channel)"]
        self.protocol_combo = ttk.Combobox(input_frame, textvariable=self.protocol_type, values=protocol_options, state="readonly")
        self.protocol_combo.grid(row=4, column=0, padx=5, pady=2, sticky="ew")

        ttk.Separator(input_frame, orient='horizontal').grid(row=5, column=0, sticky='ew', pady=5)

        # --- Graph Type Selection (UPDATED) --- (Row 6, 7)
        ttk.Label(input_frame, text="Graph Type:").grid(row=6, column=0, sticky="w", padx=5, pady=2)
        self.graph_type = tk.StringVar(value="SKR vs. Range")
        graph_options = [
            "SKR vs. Range", "Channel Loss vs. Range", "QBER vs. Range",
            "Received Power vs. Tx Power", "QBER vs. Channel Loss",
            "SKR vs. Pointing Error (ISL)", "SKR vs. Turbulence (GTS)",
            "SKR vs. QBER", "Key Rate Components vs. Range",
            "QBER vs. Eve's Attack" # <-- NEW GRAPH
        ]
        self.graph_combo = ttk.Combobox(input_frame, textvariable=self.graph_type, values=graph_options, state="readonly")
        self.graph_combo.grid(row=7, column=0, padx=5, pady=2, sticky="ew")
        self.graph_combo.bind("<<ComboboxSelected>>", self.toggle_inputs)

        ttk.Separator(input_frame, orient='horizontal').grid(row=8, column=0, sticky='ew', pady=5)

        # --- Link Type Selection --- (Row 9, 10, 11)
        self.link_type = tk.StringVar(value="ISL")
        ttk.Label(input_frame, text="Link Type:").grid(row=9, column=0, sticky="w", padx=5, pady=2)
        ttk.Radiobutton(input_frame, text="Inter-Satellite (ISL)", variable=self.link_type, value="ISL", command=self.toggle_inputs).grid(row=10, column=0, sticky="w", padx=10)
        ttk.Radiobutton(input_frame, text="Ground-to-Satellite", variable=self.link_type, value="GTS", command=self.toggle_inputs).grid(row=11, column=0, sticky="w", padx=10)

        # --- Numerical Inputs (Channel) --- (Row 12 onwards)
        self.range_label = ttk.Label(input_frame, text="Range Parameter (km):")
        self.range_label.grid(row=12, column=0, sticky="w", padx=5, pady=2)
        self.range_entry = ttk.Entry(input_frame); self.range_entry.grid(row=13, column=0, padx=5, pady=2); self.range_entry.insert(0, "500")

        ttk.Label(input_frame, text="Divergence (mrad):").grid(row=14, column=0, sticky="w", padx=5, pady=2)
        self.divergence_entry = ttk.Entry(input_frame); self.divergence_entry.grid(row=15, column=0, padx=5, pady=2); self.divergence_entry.insert(0, "0.032")

        self.turbulence_label = ttk.Label(input_frame, text="Turbulence:")
        self.turbulence_label.grid(row=16, column=0, sticky="w", padx=5, pady=2)
        self.turbulence_combo = ttk.Combobox(input_frame, values=['low', 'medium', 'high'], state="readonly");
        self.turbulence_combo.grid(row=17, column=0, padx=5, pady=2); self.turbulence_combo.set('low')

        self.pointing_label = ttk.Label(input_frame, text="Pointing Error (µrad):")
        self.pointing_label.grid(row=18, column=0, sticky="w", padx=5, pady=2)
        self.pointing_entry = ttk.Entry(input_frame); self.pointing_entry.grid(row=19, column=0, padx=5, pady=2); self.pointing_entry.insert(0, "1.0")

        ttk.Label(input_frame, text="Intrinsic QBER (%):").grid(row=20, column=0, sticky="w", padx=5, pady=2)
        self.qber_entry = ttk.Entry(input_frame); self.qber_entry.grid(row=21, column=0, padx=5, pady=2); self.qber_entry.insert(0, "0.5")

        self.qber_scale_label = ttk.Label(input_frame, text="QBER Scaling Factor:")
        self.qber_scale_label.grid(row=22, column=0, sticky="w", padx=5, pady=2)
        self.qber_scale_entry = ttk.Entry(input_frame);
        self.qber_scale_entry.grid(row=23, column=0, padx=5, pady=2); self.qber_scale_entry.insert(0, "0.1")

        # --- Hardware Inputs --- (Row 24 onwards)
        ttk.Separator(input_frame, orient='horizontal').grid(row=24, column=0, sticky='ew', pady=5)

        ttk.Label(input_frame, text="Detector Efficiency (%):").grid(row=25, column=0, sticky="w", padx=5, pady=2)
        self.det_eff_entry = ttk.Entry(input_frame)
        self.det_eff_entry.grid(row=26, column=0, padx=5, pady=2); self.det_eff_entry.insert(0, "50.0")

        ttk.Label(input_frame, text="Receiver Optics Loss (dB):").grid(row=27, column=0, sticky="w", padx=5, pady=2)
        self.rx_loss_entry = ttk.Entry(input_frame)
        self.rx_loss_entry.grid(row=28, column=0, padx=5, pady=2); self.rx_loss_entry.insert(0, "3.0")

        # --- Tx Power Widgets (Renumbered) ---
        self.tx_power_min_label = ttk.Label(input_frame, text="Min Tx Power (dBm):")
        self.tx_power_min_label.grid(row=29, column=0, sticky="w", padx=5, pady=2)
        self.tx_power_min_entry = ttk.Entry(input_frame)
        self.tx_power_min_entry.grid(row=30, column=0, padx=5, pady=2); self.tx_power_min_entry.insert(0, "-20")

        self.tx_power_max_label = ttk.Label(input_frame, text="Max Tx Power (dBm):")
        self.tx_power_max_label.grid(row=31, column=0, sticky="w", padx=5, pady=2)
        self.tx_power_max_entry = ttk.Entry(input_frame)
        self.tx_power_max_entry.grid(row=32, column=0, padx=5, pady=2); self.tx_power_max_entry.insert(0, "10")

        # --- Plot Frame ---
        plot_frame = ttk.LabelFrame(self.master, text="Simulation Results")
        plot_frame.grid(row=0, column=1, padx=10, pady=10, sticky="nsew")
        self.fig = Figure(figsize=(8, 6), dpi=100); self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # --- Control Frame ---
        control_frame = ttk.Frame(self.master)
        control_frame.grid(row=1, column=0, columnspan=2, pady=5)
        run_button = ttk.Button(control_frame, text="Run Simulation", command=self.run_simulation)
        run_button.pack()
        self.status_label = ttk.Label(self.master, text="Status: Ready")
        self.status_label.grid(row=2, column=0, columnspan=2, pady=2)

        self.master.grid_columnconfigure(1, weight=1); self.master.grid_rowconfigure(0, weight=1)
        self.toggle_inputs()

    def toggle_inputs(self, event=None):
        # (This function is updated to handle the new graph)
        is_isl = self.link_type.get() == "ISL"
        graph_choice = self.graph_type.get()
        self.turbulence_label.config(state='disabled' if is_isl else 'normal')
        self.turbulence_combo.config(state='disabled' if is_isl else 'normal')
        self.pointing_label.config(state='normal' if is_isl else 'disabled')
        self.pointing_entry.config(state='normal' if is_isl else 'disabled')
        self.qber_scale_label.config(state='disabled' if is_isl else 'normal')
        self.qber_scale_entry.config(state='disabled' if is_isl else 'normal')
        is_power_graph = graph_choice == "Received Power vs. Tx Power"
        self.tx_power_min_label.config(state='normal' if is_power_graph else 'disabled')
        self.tx_power_min_entry.config(state='normal' if is_power_graph else 'disabled')
        self.tx_power_max_label.config(state='normal' if is_power_graph else 'disabled')
        self.tx_power_max_entry.config(state='normal' if is_power_graph else 'disabled')
        
        # All these graphs require a fixed range
        if graph_choice in ["SKR vs. Pointing Error (ISL)", "SKR vs. Turbulence (GTS)", "SKR vs. QBER", "Received Power vs. Tx Power", "QBER vs. Eve's Attack"]:
             self.range_label.config(text="Fixed Range (km):")
        else:
             self.range_label.config(text="Max Range (km):")

    def update_plot(self, x_data, y_data, title, xlabel, ylabel, yscale='linear', **kwargs):
        self.ax.clear()
        valid_data = isinstance(x_data, (list, np.ndarray)) and \
                     isinstance(y_data, (list, np.ndarray)) and \
                     len(x_data) > 0 and len(x_data) == len(y_data)
        
        is_bar = kwargs.get('is_bar', False)

        if not valid_data:
            self.ax.text(0.5, 0.5, 'No valid data to plot', ha='center', va='center', transform=self.ax.transAxes)
            full_x_range_valid = False
        else:
            full_x_range_valid = True
            if not is_bar:
                try: min_x, max_x = min(x_data), max(x_data)
                except TypeError: full_x_range_valid = False
            
            if is_bar:
                 bar_labels = [str(x) for x in x_data]
                 self.ax.bar(bar_labels, y_data)
            else:
                self.ax.plot(x_data, y_data, 'o-', label=kwargs.get('label1'))
                if 'y2' in kwargs and isinstance(kwargs['y2'], (list, np.ndarray)) and len(kwargs['y2']) == len(x_data):
                    self.ax.plot(x_data, kwargs['y2'], 's--', label=kwargs.get('label2'))
                if kwargs.get('label1') or kwargs.get('label2'):
                    self.ax.legend()

        self.ax.set_title(title); self.ax.set_xlabel(xlabel); self.ax.set_ylabel(ylabel)
        self.ax.grid(True)

        if full_x_range_valid and not is_bar:
            padding = (max_x - min_x) * 0.02 if max_x > min_x else 1
            self.ax.set_xlim(min_x - padding, max_x + padding)

        try:
            has_positive_y = any(y > 0 for y in y_data) if valid_data else False
            if yscale == 'log' and has_positive_y:
                if is_bar:
                    self.ax.set_yscale('log')
                elif valid_data:
                     filtered_data = [(x, y) for x, y in zip(x_data, y_data) if y > 0]
                     if filtered_data:
                         x_filt, y_filt_main = zip(*filtered_data)
                         self.ax.clear()
                         self.ax.plot(x_filt, y_filt_main, 'o-', label=kwargs.get('label1'))
                         if 'y2' in kwargs:
                             y2_data = kwargs.get('y2')
                             if isinstance(y2_data, (list, np.ndarray)) and len(y2_data) == len(x_data):
                                 y2_filt_pairs = [(x, y2) for x, y, y2 in zip(x_data, y_data, y2_data) if y > 0 and y2 > 0]
                                 if y2_filt_pairs:
                                     x_filt_y2, y2_filt = zip(*y2_filt_pairs)
                                     self.ax.plot(x_filt_y2, y2_filt, 's--', label=kwargs.get('label2'))
                         if kwargs.get('label1') or kwargs.get('label2'): self.ax.legend()
                         self.ax.set_yscale('log')
                         self.ax.set_title(title); self.ax.set_xlabel(xlabel); self.ax.set_ylabel(ylabel); self.ax.grid(True)
                         padding = (max_x - min_x) * 0.02 if max_x > min_x else 1
                         self.ax.set_xlim(min_x - padding, max_x + padding)
                     else: self.ax.set_yscale('linear')
                else: self.ax.set_yscale('linear')
            else: self.ax.set_yscale('linear')
        except (ValueError, TypeError): self.ax.set_yscale('linear')

        self.fig.tight_layout(); self.canvas.draw()
        
    def run_simulation(self):
        self.status_label.config(text="Status: Running..."); self.master.update_idletasks()

        try:
            protocol_choice = self.protocol_type.get()
            link_type = self.link_type.get()
            range_param = float(self.range_entry.get())
            divergence = float(self.divergence_entry.get())
            wavelength_nm = float(self.wavelength_entry.get())
            graph_choice = self.graph_type.get()
            intrinsic_qber = float(self.qber_entry.get()) / 100.0
            qber_scaling_factor = float(self.qber_scale_entry.get())
            pointing_error_val = float(self.pointing_entry.get())
            turbulence_choice = self.turbulence_combo.get()
            detector_efficiency = float(self.det_eff_entry.get()) / 100.0
            receiver_optics_loss_db = float(self.rx_loss_entry.get())
        except ValueError:
            self.status_label.config(text="Status: Error - Invalid numerical input.")
            self.update_plot([], [], "Error", "Invalid Input", "")
            return

        # --- DYNAMIC PROTOCOL SELECTION (UPDATED) ---
        if protocol_choice == "Decoy-State BB84":
            protocol = DecoyBB84Protocol()
        elif protocol_choice == "MDI-QKD":
            protocol = MDIQKDProtocol()
        elif protocol_choice == "Classical (Normal Channel)":
            protocol = ClassicalProtocol()
        else:
            self.status_label.config(text="Status: Error - Unknown protocol.")
            return

        # --- Main Dispatcher Logic ---
        
        # Graphs vs. a primary variable (Range)
        if graph_choice in ["SKR vs. Range", "Channel Loss vs. Range", "QBER vs. Range", "QBER vs. Channel Loss", "Key Rate Components vs. Range"]:
            start_range_km = max(10, range_param / 10)
            ranges_km = np.linspace(start_range_km, range_param, 30)
            results = {'loss': [], 'qber': [], 'skr': [], 'useful': [], 'leakage': []}
            valid_run = True
            for r in ranges_km:
                try:
                    if link_type == 'ISL':
                        channel = InterSatelliteLink(range_km=r, divergence_mrad=divergence, wavelength_nm=wavelength_nm, pointing_error_urad=pointing_error_val)
                        props = channel.get_channel_properties(intrinsic_qber=intrinsic_qber)
                    else: # GTS
                        channel = GroundToSatelliteLink(range_km=r, divergence_mrad=divergence, wavelength_nm=wavelength_nm, turbulence_strength=turbulence_choice)
                        props = channel.get_channel_properties(intrinsic_qber=intrinsic_qber, qber_rytor_scaling=qber_scaling_factor)
                    if not all(np.isfinite(v) for v in props.values()): raise ValueError(f"Non-finite channel prop @ {r}km")
                    res = protocol.calculate_skr(props, detector_efficiency, receiver_optics_loss_db)
                    if not all(np.isfinite(v) for v in res.values()): raise ValueError(f"Non-finite SKR @ {r}km")
                    results['loss'].append(props['channel_loss_db']); results['qber'].append(props['effective_qber'])
                    results['skr'].append(res['skr']); results['useful'].append(res['useful_rate']); results['leakage'].append(res['leakage_rate'])
                except Exception as e:
                    self.status_label.config(text=f"Status: Error - {e}"); valid_run = False; break
            if valid_run:
                # Plotting logic for "vs. Range" graphs
                if graph_choice == "SKR vs. Range":
                    self.update_plot(ranges_km, results['skr'], f"{protocol_choice} on {link_type}", "Range (km)", "SKR (bits/pulse)", yscale='log')
                    self.status_label.config(text=f"Status: Complete. Final SKR @ {range_param} km: {results['skr'][-1]:.2e} bits/pulse")
                elif graph_choice == "Channel Loss vs. Range":
                    self.update_plot(ranges_km, results['loss'], f"{link_type} Channel Loss", "Range (km)", "Channel Loss (dB)", yscale='linear')
                    self.status_label.config(text=f"Status: Complete. Final Channel Loss @ {range_param} km: {results['loss'][-1]:.2f} dB")
                elif graph_choice == "QBER vs. Range":
                     qber_percent = [q*100 for q in results['qber']]
                     self.update_plot(ranges_km, qber_percent, f"{link_type} Performance", "Range (km)", "Effective QBER (%)", yscale='linear')
                     self.status_label.config(text=f"Status: Complete. Final QBER @ {range_param} km: {results['qber'][-1]:.2%}")
                elif graph_choice == "QBER vs. Channel Loss":
                    qber_percent = [q*100 for q in results['qber']]
                    valid_losses = [l for l in results['loss'] if np.isfinite(l)]
                    valid_qbers = [q for i, q in enumerate(qber_percent) if np.isfinite(results['loss'][i])]
                    if len(valid_losses) > 0 and len(valid_losses) == len(valid_qbers):
                        self.update_plot(valid_losses, valid_qbers, f"{link_type} Performance", "Channel Loss (dB)", "Effective QBER (%)", yscale='linear')
                        self.status_label.config(text=f"Status: Complete. Final QBER: {results['qber'][-1]:.2%}")
                    else:
                        self.update_plot([], [], f"{link_type} Performance", "Channel Loss (dB)", "Effective QBER (%)", status="No valid loss data.")
                elif graph_choice == "Key Rate Components vs. Range":
                    self.update_plot(ranges_km, results['useful'], f"{protocol_choice} Components on {link_type}", "Range (km)", "Rate (bits/pulse)", yscale='log',
                                     y2=results['leakage'], label1='Useful Rate', label2='Leakage Rate')
                    self.status_label.config(text=f"Status: Complete. Key component analysis.")
        
        # Graphs where range is fixed
        elif graph_choice == "Received Power vs. Tx Power":
            try:
                if link_type == 'ISL':
                    channel = InterSatelliteLink(range_km=range_param, divergence_mrad=divergence, wavelength_nm=wavelength_nm, pointing_error_urad=pointing_error_val)
                    props = channel.get_channel_properties(intrinsic_qber=intrinsic_qber)
                else:
                    channel = GroundToSatelliteLink(range_km=range_param, divergence_mrad=divergence, wavelength_nm=wavelength_nm, turbulence_strength=turbulence_choice)
                    props = channel.get_channel_properties(intrinsic_qber=intrinsic_qber, qber_rytor_scaling=qber_scaling_factor)
                if not np.isfinite(props['channel_loss_db']): raise ValueError("Channel loss calc failed")
                total_system_loss_db = props['channel_loss_db'] + receiver_optics_loss_db + (-10*np.log10(detector_efficiency) if detector_efficiency > 0 else np.inf)
                min_tx = float(self.tx_power_min_entry.get()); max_tx = float(self.tx_power_max_entry.get())
                tx_powers_dbm = np.linspace(min_tx, max_tx, 30)
                rx_powers_dbm = tx_powers_dbm - total_system_loss_db
                self.update_plot(tx_powers_dbm, rx_powers_dbm, f"{link_type} Link Budget @ {range_param} km", "Transmitted Power (dBm)", "Received Power (dBm)")
                self.status_label.config(text=f"Status: Complete. Total System Loss: {total_system_loss_db:.2f} dB")
            except Exception as e:
                 self.status_label.config(text=f"Status: Error - {e}"); self.update_plot([], [], "Error", "Calculation Error", "")

        elif graph_choice == "SKR vs. Pointing Error (ISL)":
            pointing_errors = np.linspace(0.5, 5.0, 30); skr_results = []; valid_run = True
            for pe in pointing_errors:
                try:
                    channel = InterSatelliteLink(range_km=range_param, divergence_mrad=divergence, wavelength_nm=wavelength_nm, pointing_error_urad=pe)
                    props = channel.get_channel_properties(intrinsic_qber=intrinsic_qber)
                    if not np.isfinite(props['channel_loss_db']): props['channel_loss_db'] = 200
                    res = protocol.calculate_skr(props, detector_efficiency, receiver_optics_loss_db); skr_results.append(res['skr'])
                except Exception as e: self.status_label.config(text=f"Status: Error - {e}"); valid_run = False; break
            if valid_run:
                self.update_plot(pointing_errors, skr_results, f"{protocol_choice} @ {range_param}km", "Pointing Error (µrad)", "SKR (bits/pulse)", yscale='log')
                self.status_label.config(text=f"Status: Complete. Pointing error analysis.")

        elif graph_choice == "SKR vs. Turbulence (GTS)":
            turbulence_levels = ['low', 'medium', 'high']; skr_results = []; valid_run = True
            for tl in turbulence_levels:
                 try:
                    channel = GroundToSatelliteLink(range_km=range_param, divergence_mrad=divergence, wavelength_nm=wavelength_nm, turbulence_strength=tl)
                    props = channel.get_channel_properties(intrinsic_qber=intrinsic_qber, qber_rytor_scaling=qber_scaling_factor)
                    if not np.isfinite(props['channel_loss_db']): props['channel_loss_db'] = 200
                    res = protocol.calculate_skr(props, detector_efficiency, receiver_optics_loss_db); skr_results.append(res['skr'])
                 except Exception as e: self.status_label.config(text=f"Status: Error - {e}"); valid_run = False; break
            if valid_run:
                self.update_plot(turbulence_levels, skr_results, f"{protocol_choice} @ {range_param}km", "Turbulence Strength", "SKR (bits/pulse)", is_bar=True, yscale='log')
                self.status_label.config(text=f"Status: Complete. Turbulence resilience analysis.")

        elif graph_choice == "SKR vs. QBER":
            qber_values = np.linspace(0, 0.15, 30); skr_results = []; valid_run = True
            try:
                if link_type == 'ISL':
                    channel = InterSatelliteLink(range_km=range_param, divergence_mrad=divergence, wavelength_nm=wavelength_nm, pointing_error_urad=pointing_error_val)
                    base_props = channel.get_channel_properties(intrinsic_qber=0)
                else:
                    channel = GroundToSatelliteLink(range_km=range_param, divergence_mrad=divergence, wavelength_nm=wavelength_nm, turbulence_strength=turbulence_choice)
                    base_props = channel.get_channel_properties(intrinsic_qber=0, qber_rytor_scaling=0)
                if not np.isfinite(base_props['channel_loss_db']): raise ValueError("Base loss calc failed")
                base_channel_loss = base_props['channel_loss_db']
                for qber in qber_values:
                    current_props = {'channel_loss_db': base_channel_loss, 'effective_qber': qber}
                    res = protocol.calculate_skr(current_props, detector_efficiency, receiver_optics_loss_db); skr_results.append(res['skr'])
            except Exception as e: self.status_label.config(text=f"Status: Error - {e}"); valid_run = False
            if valid_run:
                self.update_plot(qber_values * 100, skr_results, f"{protocol_choice} @ {range_param}km ({base_channel_loss:.1f} dB Channel Loss)", "Total Effective QBER (%)", "SKR (bits/pulse)", yscale='log')
                self.status_label.config(text=f"Status: Complete. Error tolerance analysis.")

        # --- NEW GRAPH LOGIC ---
        elif graph_choice == "QBER vs. Eve's Attack":
            try:
                # 1. Get the base channel QBER (no attack)
                if link_type == 'ISL':
                    channel = InterSatelliteLink(range_km=range_param, divergence_mrad=divergence, wavelength_nm=wavelength_nm, pointing_error_urad=pointing_error_val)
                    base_props = channel.get_channel_properties(intrinsic_qber=intrinsic_qber)
                else:
                    channel = GroundToSatelliteLink(range_km=range_param, divergence_mrad=divergence, wavelength_nm=wavelength_nm, turbulence_strength=turbulence_choice)
                    base_props = channel.get_channel_properties(intrinsic_qber=intrinsic_qber, qber_rytor_scaling=qber_scaling_factor)
                
                base_qber = base_props['effective_qber']
                if not np.isfinite(base_qber): raise ValueError("Base QBER calc failed")

                # 2. Create X-axis (attack intensity)
                x_attack_levels = np.linspace(0, 1, 30) # 0% to 100% attack
                
                # 3. Calculate QBER for QKD
                # Eve introduces 25% QBER on the fraction of the signal she attacks
                qkd_qber = (1 - x_attack_levels) * base_qber + (x_attack_levels * 0.25)
                
                # 4. Calculate QBER for Classical
                # QBER is unchanged, Eve's attack is undetectable
                classical_qber = np.full_like(x_attack_levels, base_qber)
                
                # 5. Plot both lines
                self.update_plot(x_attack_levels * 100, qkd_qber * 100, 
                                 f"QKD vs Classical: Eavesdropping Detection ({link_type} @ {range_param}km)", 
                                 "Eavesdropping Attack Intensity (%)", 
                                 "Measured QBER (%)", 
                                 yscale='linear',
                                 y2=classical_qber * 100, 
                                 label1="QKD (BB84/MDI)", 
                                 label2="Classical Channel")
                self.status_label.config(text=f"Status: Complete. QKD shows detectable error increase.")
            except Exception as e:
                self.status_label.config(text=f"Status: Error - {e}"); self.update_plot([], [], "Error", "Calculation Error", "")


if __name__ == "__main__":
    root = tk.Tk()
    app = SimulationGUI(root)
    root.mainloop()