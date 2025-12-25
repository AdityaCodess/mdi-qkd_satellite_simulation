import tkinter as tk
from tkinter import ttk
import numpy as np
import math
from scipy.integrate import quad
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# =============================================================================
# PHYSICS ENGINE
# =============================================================================

R_EARTH = 6371e3
H_ATMOS = 20000

def hufnagel_valley(h, A, v):
    if h < 0: h = 0
    try:
        t1 = 5.94e-53 * (v / 27)**2 * (h**10) * np.exp(-h / 1000)
        t2 = 2.7e-16 * np.exp(-h / 1500)
        t3 = A * np.exp(-h / 100)
        return t1 + t2 + t3
    except: return 0.0

def extinction_coefficient(h, wavelength_m):
    V_km = 23.0 
    lnm = wavelength_m * 1e9
    q = 1.6 if V_km > 50 else 1.3 if V_km > 6 else 0.16 * V_km + 0.34
    beta = (3.912 / V_km) * (lnm / 550.0)**(-q)
    return (beta / 1000.0) * np.exp(-h / 1200) + 1e-6

class QuantumLink:
    def __init__(self, cfg):
        self.alt_m = cfg['alt_km'] * 1000
        # Distance cannot be less than altitude for GTS
        self.z = max(cfg['range_km'] * 1000, self.alt_m) if cfg['link_type'] == "GTS" else cfg['range_km'] * 1000
        self.lam = cfg['wavelength_nm'] * 1e-9
        self.div = cfg['divergence_mrad'] * 1e-3
        self.dtx = cfg['tx_aperture']
        self.drx = cfg['rx_aperture']
        self.link = cfg['link_type']
        self.turb = cfg.get('turb_lvl', 'low')
        self.pe = cfg.get('pointing_error', 1.0)
        self.w0 = self.lam / (np.pi * self.div) if self.div > 0 else self.dtx/2

    def get_properties(self):
        if self.link == 'ISL':
            zR = np.pi * self.w0**2 / self.lam
            wd = 2 * self.w0 * np.sqrt(1 + (self.z/zR)**2)
            geo = -10 * np.log10(np.clip((self.drx/wd)**2, 1e-20, 1.0))
            pointing = 4.34 * ((self.pe * 1e-6 * self.z) / (wd/2))**2
            return {'loss': geo + pointing, 'turb_qber': 0.0}
        else:
            cos_th = self.alt_m/self.z + (self.alt_m**2 - self.z**2)/(2*self.z*R_EARTH)
            th = np.arccos(np.clip(cos_th, -1.0, 1.0))
            k = 2 * np.pi / self.lam
            tp = {'low': (1.7e-14, 21), 'medium': (2.75e-14, 21), 'high': (2.75e-14, 57)}[self.turb]
            
            def rho_int(xi):
                h = np.sqrt(R_EARTH**2 + xi**2 + 2*xi*R_EARTH*np.cos(th)) - R_EARTH
                return (1 - xi/self.z)**(5/3) * hufnagel_valley(h, *tp)
            i0, _ = quad(rho_int, 0, self.z)
            rho0 = (1.46 * k**2 * i0)**(-3/5) if i0 > 0 else 1e9
            
            def ry_int(xi):
                h = np.sqrt(R_EARTH**2 + xi**2 + 2*xi*R_EARTH*np.cos(th)) - R_EARTH
                return (xi/self.z * (1 - xi/self.z))**(5/6) * hufnagel_valley(h, *tp)
            ry, _ = quad(ry_int, 0, self.z)
            sigR2 = 2.25 * k**(7/6) * ry
            
            wd = 2 * self.w0 * np.sqrt(1 + (self.z/(np.pi*self.w0**2/self.lam))**2)
            wst = np.sqrt(wd**2 + 2*(self.z/(k*rho0))**2) if rho0 > 0 else wd
            geo = -10 * np.log10(np.clip((self.drx/wst)**2, 1e-20, 1.0))
            
            def ex_int(y):
                h = np.sqrt(R_EARTH**2 + y**2 + 2*y*R_EARTH*np.cos(th)) - R_EARTH
                return extinction_coefficient(h, self.lam) if h < H_ATMOS else 0
            dp, _ = quad(ex_int, 0, min(self.z, 50000))
            ext = -10 * np.log10(np.exp(-dp))
            return {'loss': geo + ext, 'turb_qber': sigR2}

# =============================================================================
# PROTOCOL ENGINE
# =============================================================================

def binary_entropy(p):
    if p <= 0 or p >= 1: return 0
    return -p * np.log2(p) - (1-p) * np.log2(1-p)

class ProtocolEngine:
    def __init__(self, p_type, hw):
        self.p_type = p_type
        self.hw = hw

    def calculate_skr(self, props, intrinsic_qber, qscale):
        if self.p_type == 'Classical (Normal Channel)': return {'skr': 0, 'useful': 0, 'leakage': 0, 'qber': intrinsic_qber}
        
        loss_total_db = props['loss'] + self.hw['rx_loss']
        eta = (10**(-loss_total_db/10)) * self.hw['det_eff']
        y0 = self.hw['dark_rate'] / self.hw['rep_rate']
        e_ch = intrinsic_qber + qscale * props['turb_qber']
        mu = self.hw['mu']

        if 'BB84' in self.p_type:
            y1 = 1 - (1-y0)*(1-eta)
            e1 = (0.5 * y0 + e_ch * (y1 - y0)) / y1 if y1 > 0 else 0.5
            q_mu = 1 - (1-y0)*np.exp(-mu * eta)
            e_mu = (0.5 * y0 + e_ch * (q_mu - y0)) / q_mu if q_mu > 0 else 0.5
            useful = 0.5 * mu * np.exp(-mu) * y1 * (1 - binary_entropy(e1))
            leakage = 0.5 * q_mu * self.hw['f_ec'] * binary_entropy(e_mu)
            return {'skr': max(0, useful - leakage), 'useful': useful, 'leakage': leakage, 'qber': e_mu}
        
        elif 'MDI' in self.p_type:
            y11 = (mu**2 * np.exp(-2*mu)) * (eta**2)
            useful = 0.25 * y11 * (1 - binary_entropy(e_ch))
            q_tot = eta**2 
            e_tot = (e_ch * q_tot + 0.5 * y0) / (q_tot + y0) if q_tot > 0 else 0.5
            leakage = 0.25 * q_tot * self.hw['f_ec'] * binary_entropy(e_tot)
            return {'skr': max(0, useful - leakage), 'useful': useful, 'leakage': leakage, 'qber': e_tot}

# =============================================================================
# GUI APPLICATION
# =============================================================================

class SimulationGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("QKD Satellite Simulator Pro")
        self.root.geometry("1400x950")
        
        style = ttk.Style()
        style.configure("Header.TLabelframe.Label", font=("Segoe UI", 11, "bold"))

        main_frame = ttk.Frame(root, padding=10)
        main_frame.pack(fill=tk.BOTH, expand=True)

        settings_pane = ttk.LabelFrame(main_frame, text="Configurations", padding=10)
        settings_pane.pack(side=tk.LEFT, fill=tk.Y, padx=5)
        
        cv = tk.Canvas(settings_pane, width=340, highlightthickness=0)
        sb = ttk.Scrollbar(settings_pane, orient="vertical", command=cv.yview)
        self.f = ttk.Frame(cv)
        self.f.bind("<Configure>", lambda e: cv.configure(scrollregion=cv.bbox("all")))
        cv.create_window((0, 0), window=self.f, anchor="nw")
        cv.configure(yscrollcommand=sb.set)
        cv.pack(side="left", fill="both", expand=True)
        sb.pack(side="right", fill="y")

        self.setup_inputs()

        results_container = ttk.Frame(main_frame)
        results_container.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=5)

        self.dashboard_frame = ttk.LabelFrame(results_container, text="Real-Time System Performance Dashboard", padding=15)
        self.dashboard_frame.pack(side=tk.TOP, fill=tk.X, pady=(0, 10))

        self.skr_val = tk.StringVar(value="--")
        self.avg_skr_val = tk.StringVar(value="--")
        self.loss_val = tk.StringVar(value="--")
        self.qber_val = tk.StringVar(value="--")

        db_grid = ttk.Frame(self.dashboard_frame)
        db_grid.pack(fill=tk.X)
        self.create_db_item(db_grid, "PEAK SKR (ZENITH)", self.skr_val, 0, "bits/pulse", "#2980b9")
        self.create_db_item(db_grid, "AVERAGE SKR (PASS)", self.avg_skr_val, 1, "bits/pulse", "#8e44ad")
        self.create_db_item(db_grid, "TOTAL SYSTEM LOSS", self.loss_val, 2, "dB", "#2c3e50")
        self.create_db_item(db_grid, "EFFECTIVE QBER", self.qber_val, 3, "%", "#c0392b")

        plot_pane = ttk.LabelFrame(results_container, text="Performance Simulation Curve", padding=10)
        plot_pane.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
        self.fig = Figure(figsize=(8, 6), dpi=100); self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_pane)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def create_db_item(self, master, label, var, col, unit, color):
        frame = tk.Frame(master, bg="#fdfdfd", highlightbackground="#e1e1e1", highlightthickness=1, padx=10, pady=10)
        frame.grid(row=0, column=col, padx=5, sticky="ew")
        master.grid_columnconfigure(col, weight=1)
        tk.Label(frame, text=label, font=("Segoe UI", 8, "bold"), fg="#7f8c8d", bg="#fdfdfd").pack()
        tk.Label(frame, textvariable=var, font=("Courier New", 14, "bold"), fg=color, bg="#fdfdfd").pack()
        tk.Label(frame, text=unit, font=("Segoe UI", 7), fg="#bdc3c7", bg="#fdfdfd").pack()

    def setup_inputs(self):
        self.protocol = tk.StringVar(value="Decoy-State BB84")
        self.add_cb("Protocol:", self.protocol, ["Decoy-State BB84", "MDI-QKD", "Classical (Normal Channel)"])
        self.link = tk.StringVar(value="GTS")
        ttk.Label(self.f, text="Link Type:").pack(anchor="w")
        ttk.Radiobutton(self.f, text="ISL (Vacuum)", variable=self.link, value="ISL", command=self.toggle).pack(anchor="w", padx=10)
        ttk.Radiobutton(self.f, text="GTS (Atmosphere)", variable=self.link, value="GTS", command=self.toggle).pack(anchor="w", padx=10)
        self.graph = tk.StringVar(value="SKR vs. Range")
        g_opts = ["SKR vs. Range", "Total System Loss vs. Range", "QBER vs. Range", "Received Power vs. Tx Power", "QBER vs. Channel Loss", "SKR vs. Pointing Error (ISL)", "SKR vs. Turbulence (GTS)", "SKR vs. QBER", "Key Rate Components vs. Range", "QBER vs. Eve's Attack"]
        self.add_cb("Graph Type:", self.graph, g_opts)
        self.params = {}
        h_data = [("Wavelength (nm)", "850"), ("Satellite Altitude (km)", "500"), ("Max Slant Range (km)", "1000"), ("Divergence (mrad)", "0.006"), ("Tx Aperture (m)", "0.3"), ("Rx Aperture (m)", "1.2"), ("Dark Rate (Hz)", "500"), ("Rep Rate (MHz)", "100"), ("Det. Efficiency (%)", "50.0"), ("Rx Optics Loss (dB)", "3.0")]
        for l, d in h_data: self.add_ent(l, d)
        self.turb = tk.StringVar(value="low")
        self.turb_lbl = ttk.Label(self.f, text="Turbulence Strength:")
        self.turb_lbl.pack(anchor="w")
        self.turb_cb = ttk.Combobox(self.f, textvariable=self.turb, values=["low", "medium", "high"], state="readonly")
        self.turb_cb.pack(fill="x", padx=5, pady=2)
        p_data = [("Pointing Error (µrad)", "1.0"), ("Intrinsic QBER (%)", "1.2"), ("Signal Intensity (mu)", "0.5"), ("Error Corr. Eff. (f)", "1.16"), ("QBER Scaling Factor", "0.15"), ("Min Tx Power (dBm)", "-20"), ("Max Tx Power (dBm)", "10")]
        for l, d in p_data: self.add_ent(l, d)
        ttk.Button(self.f, text="RUN MISSION SIMULATION", command=self.run).pack(fill="x", pady=15)
        self.status = ttk.Label(self.f, text="Ready", foreground="blue"); self.status.pack(); self.toggle()

    def add_cb(self, l, v, o):
        ttk.Label(self.f, text=l).pack(anchor="w")
        ttk.Combobox(self.f, textvariable=v, values=o, state="readonly").pack(fill="x", padx=5, pady=2)

    def add_ent(self, l, d):
        ttk.Label(self.f, text=l).pack(anchor="w")
        e = ttk.Entry(self.f); e.insert(0, d); e.pack(fill="x", padx=5, pady=2); self.params[l] = e

    def toggle(self):
        st = 'normal' if self.link.get() == "GTS" else 'disabled'
        self.turb_lbl.config(state=st); self.turb_cb.config(state=st)
        self.params["Pointing Error (µrad)"].config(state='disabled' if self.link.get() == "GTS" else 'normal')

    def run(self):
        self.status.config(text="Processing...", foreground="red"); self.root.update_idletasks()
        try:
            cfg = {
                'wavelength_nm': float(self.params["Wavelength (nm)"].get()),
                'alt_km': float(self.params["Satellite Altitude (km)"].get()),
                'range_km': float(self.params["Satellite Altitude (km)"].get()),
                'divergence_mrad': float(self.params["Divergence (mrad)"].get()),
                'tx_aperture': float(self.params["Tx Aperture (m)"].get()),
                'rx_aperture': float(self.params["Rx Aperture (m)"].get()),
                'link_type': self.link.get(), 'turb_lvl': self.turb.get(),
                'pointing_error': float(self.params["Pointing Error (µrad)"].get()),
                'dark_rate': float(self.params["Dark Rate (Hz)"].get()),
                'rep_rate': float(self.params["Rep Rate (MHz)"].get()) * 1e6,
                'det_eff': float(self.params["Det. Efficiency (%)"].get()) / 100.0,
                'rx_loss': float(self.params["Rx Optics Loss (dB)"].get()),
                'mu': float(self.params["Signal Intensity (mu)"].get()),
                'f_ec': float(self.params["Error Corr. Eff. (f)"].get()),
                'intrinsic_qber': float(self.params["Intrinsic QBER (%)"].get()) / 100.0,
                'qscale': float(self.params["QBER Scaling Factor"].get())
            }

            self.ax.clear(); g = self.graph.get(); pe = ProtocolEngine(self.protocol.get(), cfg)
            
            # Helper to calculate Total System Loss (Dynamic)
            def get_total_sys_loss(channel_loss_db):
                hw_loss = cfg['rx_loss']
                det_penalty = -10 * np.log10(cfg['det_eff']) if cfg['det_eff'] > 0 else 100
                return channel_loss_db + hw_loss + det_penalty

            # Zenith Snapshot for Dashboard
            z_props = QuantumLink(cfg).get_properties()
            z_res = pe.calculate_skr(z_props, cfg['intrinsic_qber'], cfg['qscale'])
            
            self.skr_val.set(f"{z_res['skr']:.2e}")
            self.loss_val.set(f"{get_total_sys_loss(z_props['loss']):.2f}")
            self.qber_val.set(f"{z_res['qber']*100:.2f}")

            # --- GRAPH BRANCHES ---

            if g in ["SKR vs. Range", "Total System Loss vs. Range", "QBER vs. Range", "Key Rate Components vs. Range"]:
                start = cfg['alt_km'] if cfg['link_type'] == "GTS" else 10.0
                end = float(self.params["Max Slant Range (km)"].get())
                x = np.linspace(start, end, 30); y1, y2 = [], []
                for r in x:
                    c = cfg.copy(); c['range_km'] = r
                    p = QuantumLink(c).get_properties()
                    res = pe.calculate_skr(p, cfg['intrinsic_qber'], cfg['qscale'])
                    if "SKR" in g: y1.append(res['skr'])
                    elif "Loss" in g: y1.append(get_total_sys_loss(p['loss']))
                    elif "QBER" in g: y1.append(res['qber']*100)
                    elif "Components" in g: y1.append(res['useful']); y2.append(res['leakage'])
                
                if "SKR" in g: self.avg_skr_val.set(f"{np.mean(y1):.2e}")
                self.ax.plot(x, y1, 'o-', color='#2980b9', label="Primary" if y2 else None)
                if y2: self.ax.plot(x, y2, 's--', color='#e67e22', label="Leakage")
                if "SKR" in g or "Components" in g: self.ax.set_yscale('log')
                self.ax.set_xlabel("Slant Range (km)")
                self.ax.set_ylabel(g.split(' ')[0])

            elif g == "Received Power vs. Tx Power":
                tx = np.linspace(float(self.params["Min Tx Power (dBm)"].get()), float(self.params["Max Tx Power (dBm)"].get()), 20)
                rx = tx - z_props['loss'] - cfg['rx_loss']
                self.ax.plot(tx, rx, 'o-', color='#2980b9')
                self.ax.set_xlabel("Tx Power (dBm)"); self.ax.set_ylabel("Rx Power (dBm)")

            elif g == "QBER vs. Channel Loss":
                losses = np.linspace(10, 60, 20); y = []
                for l in losses:
                    eta = (10**(-l/10)) * cfg['det_eff']; y0 = cfg['dark_rate']/cfg['rep_rate']
                    q = (0.5 * y0 + cfg['intrinsic_qber'] * eta) / (eta + y0)
                    y.append(q*100)
                self.ax.plot(losses, y, 'x-', color='#c0392b'); self.ax.set_xlabel("Channel Loss (dB)"); self.ax.set_ylabel("QBER (%)")

            elif g == "SKR vs. Pointing Error (ISL)":
                x = np.linspace(0.1, 5.0, 20); y = []
                for p in x:
                    c = cfg.copy(); c['pointing_error'] = p
                    p_props = QuantumLink(c).get_properties()
                    y.append(pe.calculate_skr(p_props, cfg['intrinsic_qber'], cfg['qscale'])['skr'])
                self.ax.plot(x, y, 'o-', color='#8e44ad'); self.ax.set_yscale('log'); self.ax.set_xlabel("Pointing Error (µrad)"); self.ax.set_ylabel("SKR")
                self.avg_skr_val.set(f"{np.mean(y):.2e}")

            elif g == "SKR vs. Turbulence (GTS)":
                lvls = ['low', 'medium', 'high']; y = []
                for l in lvls:
                    c = cfg.copy(); c['turb_lvl'] = l
                    p_props = QuantumLink(c).get_properties(); y.append(pe.calculate_skr(p_props, cfg['intrinsic_qber'], cfg['qscale'])['skr'])
                self.ax.bar(lvls, y, color=['#3498db', '#e67e22', '#e74c3c']); self.ax.set_yscale('log'); self.ax.set_ylabel("SKR")
                self.avg_skr_val.set(f"{np.mean(y):.2e}")

            elif g == "SKR vs. QBER":
                x = np.linspace(0, 0.15, 20); y = []
                for q in x: y.append(pe.calculate_skr(z_props, q, 0)['skr'])
                self.ax.plot(x*100, y, 'o-', color='#27ae60'); self.ax.set_yscale('log'); self.ax.set_xlabel("Intrinsic QBER (%)"); self.ax.set_ylabel("SKR")
                self.avg_skr_val.set(f"{np.mean(y):.2e}")

            elif g == "QBER vs. Eve's Attack":
                x = np.linspace(0, 1, 20); base_q = cfg['intrinsic_qber'] + cfg['qscale'] * z_props['turb_qber']
                self.ax.plot(x*100, ((1-x)*base_q + x*0.25)*100, 'o-', label="QKD", color='#2980b9')
                self.ax.plot(x*100, np.full_like(x, base_q)*100, '--', label="Classical", color='#7f8c8d')
                self.ax.legend(); self.ax.set_xlabel("Attack Intensity (%)"); self.ax.set_ylabel("QBER (%)")

            if self.ax.get_legend_handles_labels()[0]: self.ax.legend()
            self.ax.set_title(g); self.ax.grid(True, alpha=0.3); self.fig.tight_layout(); self.canvas.draw()
            self.status.config(text="Simulation Complete", foreground="#27ae60")
        except Exception as e: self.status.config(text=f"Error: {e}", foreground="#c0392b")

if __name__ == "__main__":
    root = tk.Tk(); SimulationGUI(root); root.mainloop()