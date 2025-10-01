# QKD Satellite Performance Simulator

A Python-based simulation tool with a graphical user interface (GUI) to analyze and compare the performance of Quantum Key Distribution (QKD) protocols over satellite communication links.

This simulator models both atmospheric ground-to-satellite (GTS) links and vacuum-based inter-satellite links (ISL), allowing for a detailed comparative analysis of protocol robustness and efficiency under different physical conditions.

---

## 🛰️ Key Features

- **Multiple Protocols:** Simulate and directly compare the performance of **Decoy-State BB84** and **Measurement-Device-Independent QKD (MDI-QKD)**.
- **Dual Channel Models:**
  - **Ground-to-Satellite (GTS):** Models the effects of atmospheric loss and **turbulence**.
  - **Inter-Satellite Link (ISL):** Models extreme geometric loss and the impact of **pointing errors** in a vacuum.
- **Interactive GUI:** A user-friendly interface built with Tkinter to configure all simulation parameters without editing code.
- **Comprehensive Analysis:** Generate a wide range of analytical plots, including:
  - SKR vs. Range
  - Total Loss vs. Range
  - QBER vs. Total Loss
  - SKR vs. Pointing Error/Turbulence
  - ...and more.

---

## 📂 Folder Structure

The project is organized in a modular and scalable structure.

```
QKD_Satellite_Simulation/
│
├── gui.py # Main script to launch the UI.
├── requirements.txt # Lists Python libraries needed.
├── README.md # This file.
│
├── qkd_simulation/ # Main source code package.
│ ├── init.py
│ ├── parameters.py # Stores all fixed system parameters.
│ ├── channels.py # Defines the GTS and ISL channel models.
│ ├── protocol.py # Implements the BB84 and MDI-QKD logic.
│ └── visualizer.py # (Currently integrated into gui.py)
│
└── results/ # Directory to save outputs (optional).
├── plots/
└── data/
```

---

## 🚀 Getting Started

### Prerequisites

- Python 3.7+
- `pip` and `venv`

### Installation

1.  **Clone the repository:**

    ```bash
    git clone <your-repository-url>
    cd QKD_Satellite_Simulation
    ```

2.  **Create and activate a virtual environment (recommended):**

    - On Windows:
      ```bash
      python -m venv venv
      .\venv\Scripts\activate
      ```
    - On macOS/Linux:
      ```bash
      python3 -m venv venv
      source venv/bin/activate
      ```

3.  **Install the required libraries:**
    Create a file named `requirements.txt` with the following content:
    ```
    numpy
    matplotlib
    ```
    Then, run the installation command:
    ```bash
    pip install -r requirements.txt
    ```

### Usage

To run the simulation and launch the graphical interface, execute the `gui.py` script:

```bash
python gui.py
```

## ⚙️ Simulation Parameters Explained

**The GUI allows you to control the following parameters:**

Protocol: Choose between Decoy-State BB84 and MDI-QKD. MDI-QKD offers higher security against detector attacks but its key rate scales quadratically with loss ($O(\eta^2)$), while BB84 scales linearly ($O(\eta)$).

Graph Type: Select the relationship you want to analyze from the dropdown menu.

Link Type:

Inter-Satellite (ISL): A vacuum link where pointing error is the main variable challenge.

Ground-to-Satellite (GTS): An atmospheric link where turbulence is the main variable challenge.

Range Parameter (km): Acts as the "Max Range" for plots vs. Range, or the "Fixed Range" for other analyses.

Divergence (mrad): The beam spread angle. A smaller value means less geometric loss but requires better pointing.

Turbulence: ('low', 'medium', 'high') Simulates the effect of atmospheric distortion on GTS links.

Pointing Error (µrad): Simulates the pointing jitter for ISLs.

Intrinsic QBER (%): The baseline error rate from hardware imperfections, independent of the channel.

```

```
