# gui.py
import tkinter as tk
from tkinter import ttk
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# Import the simulation backend
from qkd_simulation.channels import GroundToSatelliteLink, InterSatelliteLink
from qkd_simulation.protocol import MDIQKDProtocol

class SimulationGUI:
    def __init__(self, master):
        self.master = master
        master.title("MDI-QKD Satellite Simulation")
        master.geometry("1000x700")
        self.create_widgets()

    def create_widgets(self):
        # --- Input Frame ---
        input_frame = ttk.LabelFrame(self.master, text="Input Parameters")
        input_frame.grid(row=0, column=0, padx=10, pady=10, sticky="ns")

        # --- Graph Type Selection (UPDATED) ---
        ttk.Label(input_frame, text="Graph Type:").grid(row=0, column=0, sticky="w", padx=5, pady=5)
        self.graph_type = tk.StringVar(value="SKR vs. Range")
        # Added the new graph option
        graph_options = [
            "SKR vs. Range", 
            "Total Loss vs. Range", 
            "QBER vs. Range", 
            "Received Power vs. Tx Power",
            "QBER vs. Total Loss" # <-- NEW OPTION
        ]
        self.graph_combo = ttk.Combobox(input_frame, textvariable=self.graph_type, values=graph_options, state="readonly")
        self.graph_combo.grid(row=1, column=0, padx=5, pady=2, sticky="ew")
        self.graph_combo.bind("<<ComboboxSelected>>", self.toggle_inputs)


        ttk.Separator(input_frame, orient='horizontal').grid(row=2, column=0, sticky='ew', pady=10)

        # Link Type Selection
        self.link_type = tk.StringVar(value="ISL")
        ttk.Label(input_frame, text="Link Type:").grid(row=3, column=0, sticky="w", padx=5, pady=2)
        ttk.Radiobutton(input_frame, text="Inter-Satellite (ISL)", variable=self.link_type, value="ISL", command=self.toggle_inputs).grid(row=4, column=0, sticky="w", padx=10)
        ttk.Radiobutton(input_frame, text="Ground-to-Satellite", variable=self.link_type, value="GTS", command=self.toggle_inputs).grid(row=5, column=0, sticky="w", padx=10)

        # --- Numerical Inputs ---
        self.range_label = ttk.Label(input_frame, text="Max Range (km):")
        self.range_label.grid(row=6, column=0, sticky="w", padx=5, pady=2)
        self.range_entry = ttk.Entry(input_frame); self.range_entry.grid(row=7, column=0, padx=5, pady=2); self.range_entry.insert(0, "4000")

        ttk.Label(input_frame, text="Divergence (mrad):").grid(row=8, column=0, sticky="w", padx=5, pady=2)
        self.divergence_entry = ttk.Entry(input_frame); self.divergence_entry.grid(row=9, column=0, padx=5, pady=2); self.divergence_entry.insert(0, "0.005")
        
        self.turbulence_label = ttk.Label(input_frame, text="Turbulence:")
        self.turbulence_label.grid(row=10, column=0, sticky="w", padx=5, pady=2)
        self.turbulence_combo = ttk.Combobox(input_frame, values=['low', 'medium', 'high'], state="readonly"); 
        self.turbulence_combo.grid(row=11, column=0, padx=5, pady=2); self.turbulence_combo.set('medium')

        self.pointing_label = ttk.Label(input_frame, text="Pointing Error (µrad):")
        self.pointing_label.grid(row=12, column=0, sticky="w", padx=5, pady=2)
        self.pointing_entry = ttk.Entry(input_frame); self.pointing_entry.grid(row=13, column=0, padx=5, pady=2); self.pointing_entry.insert(0, "1.0")

        # --- Plot Frame ---
        plot_frame = ttk.LabelFrame(self.master, text="Simulation Results")
        plot_frame.grid(row=0, column=1, padx=10, pady=10, sticky="nsew")
        
        self.fig = Figure(figsize=(8, 6), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # --- Control Frame ---
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
        self.turbulence_label.config(state='disabled' if is_isl else 'normal')
        self.turbulence_combo.config(state='disabled' if is_isl else 'normal')
        self.pointing_label.config(state='normal' if is_isl else 'disabled')
        self.pointing_entry.config(state='normal' if is_isl else 'disabled')
        
        # Update the Range label based on context
        graph_choice = self.graph_type.get()
        if graph_choice == "Received Power vs. Tx Power":
             self.range_label.config(text="Fixed Range (km):")
        else:
             self.range_label.config(text="Max Range (km):")


    def update_plot(self, x_data, y_data, title, xlabel, ylabel, yscale='linear'):
        self.ax.clear()
        self.ax.plot(x_data, y_data, 'o-')
        self.ax.set_title(title); self.ax.set_xlabel(xlabel); self.ax.set_ylabel(ylabel)
        self.ax.grid(True); 
        
        if yscale == 'log' and y_data and any(y > 0 for y in y_data):
            self.ax.set_yscale('log')
        else:
            self.ax.set_yscale('linear')
            
        self.fig.tight_layout()
        self.canvas.draw()

    def run_simulation(self):
        self.status_label.config(text="Status: Running..."); self.master.update_idletasks()

        link_type = self.link_type.get()
        analysis_range = float(self.range_entry.get())
        divergence = float(self.divergence_entry.get())
        graph_choice = self.graph_type.get()
        
        # Logic for "Received Power vs. Tx Power"
        if graph_choice == "Received Power vs. Tx Power":
            if link_type == 'ISL':
                pointing_error = float(self.pointing_entry.get())
                channel = InterSatelliteLink(range_km=analysis_range, divergence_mrad=divergence, pointing_error_urad=pointing_error)
            else:
                turbulence = self.turbulence_combo.get()
                channel = GroundToSatelliteLink(range_km=analysis_range, divergence_mrad=divergence, turbulence_strength=turbulence)
            
            props = channel.get_channel_properties()
            total_loss_db = props['total_loss_db']
            
            tx_powers_dbm = np.linspace(-20, 10, 30)
            rx_powers_dbm = tx_powers_dbm - total_loss_db
            
            self.update_plot(tx_powers_dbm, rx_powers_dbm, 
                             f"{link_type} Link Budget at {analysis_range} km", 
                             "Transmitted Power (dBm)", "Received Power (dBm)")
            self.status_label.config(text=f"Status: Complete. Total link loss: {total_loss_db:.2f} dB")
            return

        # --- Main simulation loop for graphs vs. Range ---
        protocol = MDIQKDProtocol()
        ranges_km = np.linspace(max(100, analysis_range/10), analysis_range, 30)
        
        results = {'loss': [], 'qber': [], 'skr': []}

        for r in ranges_km:
            if link_type == 'ISL':
                pointing_error = float(self.pointing_entry.get())
                channel = InterSatelliteLink(range_km=r, divergence_mrad=divergence, pointing_error_urad=pointing_error)
            else:
                turbulence = self.turbulence_combo.get()
                channel = GroundToSatelliteLink(range_km=r, divergence_mrad=divergence, turbulence_strength=turbulence)
            
            props = channel.get_channel_properties()
            skr = protocol.calculate_skr(props)
            
            results['loss'].append(props['total_loss_db'])
            results['qber'].append(props['effective_qber'])
            results['skr'].append(skr)

        # --- Updated plotting logic ---
        if graph_choice == "SKR vs. Range":
            self.update_plot(ranges_km, results['skr'], f"{link_type} Performance", "Range (km)", "Secure Key Rate (bits/pulse)", yscale='log')
        elif graph_choice == "Total Loss vs. Range":
            self.update_plot(ranges_km, results['loss'], f"{link_type} Performance", "Range (km)", "Total Loss (dB)")
        elif graph_choice == "QBER vs. Range":
            self.update_plot(ranges_km, results['qber'], f"{link_type} Performance", "Range (km)", "Effective QBER")
        elif graph_choice == "QBER vs. Total Loss": # <-- NEW LOGIC BLOCK
            self.update_plot(results['loss'], results['qber'], f"{link_type} Performance", "Total Loss (dB)", "Effective QBER")

        final_skr = results['skr'][-1]
        self.status_label.config(text=f"Status: Complete. Final SKR at {analysis_range} km: {final_skr:.2e} bits/pulse")

if __name__ == "__main__":
    root = tk.Tk()
    app = SimulationGUI(root)
    root.mainloop()