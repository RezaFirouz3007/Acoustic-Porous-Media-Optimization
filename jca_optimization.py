######################################################################
############## Reminder change the thickness d #######################
######################################################################
# STEP 0 — What we are building (one sentence)
# Input: foam parameters
# Output: absorption coefficient from an impedance tube (rigid backing)
######################################################################
# STEP 1 — Import basic Python tools
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize             #add optimizer import
import matplotlib
import os
import sys
import traceback
matplotlib.use('Agg')
######################################################################
# STEP 2 — Define frequency axis, acoustics works in angular frequency
f = np.linspace(100, 4000, 500)  # Hz
omega = 2 * np.pi * f
######################################################################
# STEP 3 — Air constants (fixed problem constants)
rho0 = 1.21          # kg/m3            air density
c0 = 343             # m/s              speed of sound in air
Z0 = rho0 * c0       # air impedance    chcaracteristic impedance of air
eta = 1.84e-5        # Pa·s             dynamic viscosity of air
gamma = 1.4          # adiabatic index  ratio of specific heats
P0 = 101325          # Pa               atmospheric pressure
Pr = 0.71            # Prandtl number   ratio of momentum diffusivity to thermal diffusivity
######################################################################
# STEP 4 — Foam physical parameters (input parameters - initial guess)
phi = 0.95           # porosity
sigma = 12000        # flow resistivity (Pa·s/m2)
alpha_inf = 1.1      # tortuosity
Lambda = 100e-6      # viscous length (m)
Lambda_p = 200e-6    # thermal length (m)
d = 0.075             # thickness (m)
######################################################################
# STEP 5 — Effective density (Johnson model)
j = 1j

term = np.sqrt(
    1 + (4 * j * alpha_inf**2 * eta * rho0 * omega) /
    (sigma**2 * phi**2 * Lambda**2)
)

rho_eff = (rho0 * alpha_inf / phi) * (
    1 + (sigma * phi) / (j * omega * rho0 * alpha_inf) * term
)
# meaning of term: accounts for viscous effects in the pores (viscous losses and controls low frequencies)
#####################################################################
# STEP 6 — Effective bulk modulus (Champoux-Allard model)
thermal_term = np.sqrt(
    1 + (j * omega * rho0 * Pr * Lambda_p**2) / (16 * eta)
)

K_eff = (gamma * P0 / phi) / (
    gamma - (gamma - 1) /
    (1 + (8 * eta) / (j * omega * rho0 * Pr * Lambda_p**2) * thermal_term)
)
# meaning of thermal_term: accounts for thermal exchanges in the pores (thermal losses and controls high frequencies)
######################################################################
# STEP 7 — Wave number and characteristic impedance of the foam
k = omega * np.sqrt(rho_eff / K_eff)
#This tells how waves decay inside foam.
######################################################################
# STEP 8 — Surface impedance of the foam layer with rigid backing
Zc = np.sqrt(rho_eff * K_eff)
# Characteristic impedance of the foam, not yet what the tube measures yet
# Meaning: ratio of pressure to particle velocity in an infinitely thick layer of foam
#######################################################################
# STEP 9 — Absorption coefficient calculation (rigid backing)
Zs = -1j * Zc / np.tan(k * d)
# Surface impedance of the foam layer (with thickness d)
# Meaning: surface impedance of a finite-thickness foam layer
# This is exactly what the tube measures
# This accounts for the phase shift due to the finite thickness of the foam layer
########################################################################
# STEP 10&11 — Reflection coefficient and absorption coefficient 
R = (Zs - Z0) / (Zs + Z0)       # Reflection coefficient
alpha = 1 - np.abs(R)**2        # Absorption coefficient (final result)
########################################################################
# STEP 12 — Plot the absorption coefficient
plt.semilogx(f, alpha)
plt.xlabel("Frequency [Hz]")
plt.ylabel("Absorption coefficient")
plt.grid(True)
plt.xlim(left=100)
plt.savefig("initial_alpha.png", dpi=150, bbox_inches="tight")
plt.close()
print("Saved initial plot: initial_alpha.png")
######################################################################
######Comparing with experimental data#####
#
#
# Step 13 — Load experimental data (from an excel file named "measured_alpha.xlsx" with columns "frequency" and "absorption")
import glob
import pandas as pd

xlsx_files = glob.glob("*.xlsx")

if not xlsx_files:
    print("No .xlsx files found in", os.getcwd())
    sys.exit(1)

df = pd.read_excel(xlsx_files[0])
print(f"Loaded file: {xlsx_files[0]}")

f_meas = df["frequency"].to_numpy()
alpha_meas = df["absorption"].to_numpy()
# Step 14 — Use the SAME frequency grid (avoid interpolation errors))
f = f_meas
omega = 2 * np.pi * f
# Model and experiment must share frequencies
####################################################################
# Step 15 — Recompute JCA model on measured frequencies
# Effective density
##Wrap JCA model into a function
def jca_absorption(params, f):

    phi, sigma, alpha_inf, Lambda, Lambda_p = params
    
    omega = 2 * np.pi * f

    # --- Effective density (Johnson)
    term = np.sqrt(
        1 + (4j * alpha_inf**2 * eta * rho0 * omega) /
        (sigma**2 * phi**2 * Lambda**2)
    )

    rho_eff = (rho0 * alpha_inf / phi) * (
        1 + (sigma * phi) / (1j * omega * rho0 * alpha_inf) * term
    )

    # --- Effective bulk modulus (Champoux–Allard)
    thermal_term = np.sqrt(
        1 + (1j * omega * rho0 * Pr * Lambda_p**2) / (16 * eta)
    )

    K_eff = (gamma * P0 / phi) / (
        gamma - (gamma - 1) /
        (1 + (8 * eta) / (1j * omega * rho0 * Pr * Lambda_p**2) * thermal_term)
    )
    ###############################################################
    # Step 16 — Compute model absorption
    
    k = omega * np.sqrt(rho_eff / K_eff)
    Zc = np.sqrt(rho_eff * K_eff)
    Zs = -1j * Zc / np.tan(k * d)

    R = (Zs - Z0) / (Zs + Z0)
    alpha_model = 1 - np.abs(R)**2

    return alpha_model

##Define error function
def error_function(params):
    alpha_model = jca_absorption(params, f_meas)

    # Example weighting options:

    # OPTION A — Emphasize low frequencies
    # weights = 1 / f_meas
    # OPTION B — Emphasize peak region (example 500–1500 Hz)
    # weights = np.where((f_meas > 500) & (f_meas < 1500), 2.0, 1.0)
    # OPTION C — Emphasize low frequencies
    weights = np.ones_like(f_meas)
    # Slight emphasis around speech band
    weights[(f_meas > 500) & (f_meas < 2000)] = 1.5

    weighted_error = weights * (alpha_meas - alpha_model)**2

    return np.mean(weighted_error)

##Initial guess (material parameters)
phi = 0.95
sigma = 15000
alpha_inf = 1.2
Lambda = 120e-6
Lambda_p = 250e-6

##Run Optimization 
x0 = [phi, sigma, alpha_inf, Lambda, Lambda_p]

bounds = [
    (0.90, 0.99),     # porosity
    (5000, 30000),     # flow resistivity
    (1.0, 2.0),        # tortuosity
    (50e-6, 250e-6),   # viscous length
    (100e-6, 500e-6)    # thermal length
]

# Multi-start optimization
best_result = None
best_error = np.inf

for i in range(20):  # try 10–50 runs

    x0_random = [
        np.random.uniform(0.97, 0.999),      # porosity
        np.random.uniform(8000, 25000),      # resistivity
        np.random.uniform(1.0, 1.3),         # tortuosity
        np.random.uniform(30e-6, 150e-6),    # Lambda
        np.random.uniform(80e-6, 300e-6)     # Lambda_p
    ]

    result_try = minimize(
        error_function,
        x0_random,
        bounds=bounds,
        method='L-BFGS-B'
    )

    if result_try.fun < best_error:
        best_error = result_try.fun
        best_result = result_try

result = best_result
print("\n✅ Best solution from multistart:")
print("Error:", result.fun)


phi, sigma, alpha_inf, Lambda, Lambda_p = result.x

print("\nFitted parameters:")
print(f"Porosity: {phi:.4f}")
print(f"Flow resistivity: {sigma:.1f}")
print(f"Tortuosity: {alpha_inf:.3f}")
print(f"Viscous length: {Lambda*1e6:.1f} µm")
print(f"Thermal length: {Lambda_p*1e6:.1f} µm")

if not result.success:
    print("⚠ Optimization warning:", result.message)

alpha_model = jca_absorption(result.x, f_meas)

###############################################################
# Step 17 — Plot model vs experiment
plt.semilogx(f, alpha_meas, "r", label="Measured")
plt.semilogx(f, alpha_model, "-", label="JCA model")

plt.xlabel("Frequency [Hz]")
plt.ylabel("Absorption coefficient")
plt.legend()
plt.grid(True)
plt.xlim(left=100)
plt.savefig("JCA_fit_plot.png", dpi=150, bbox_inches="tight")
plt.close()
print("Saved fit plot: JCA_fit_plot.png")
#######################################################################
# Step 18 — Quantify the error
mse = np.mean((alpha_meas - alpha_model)**2)
rmse = np.sqrt(mse)

ss_res = np.sum((alpha_meas - alpha_model)**2)
ss_tot = np.sum((alpha_meas - np.mean(alpha_meas))**2)
r2 = 1 - (ss_res / ss_tot)

print("\nFit quality metrics:")
print(f"MSE  : {mse:.6f}")      #→ 0 is perfect fit
print(f"RMSE : {rmse:.6f}")     #small (<0.05 very good)
print(f"R²   : {r2:.4f}")       #R² > 0.95 → 1 is perfect fit

#######################################################################
# Step 19 — Sensitivity analysis (optional but insightful)
# Vary each parameter by ±10% and observe changes in the absorption curve
sensitivity_params = [
    (phi * 1.1, sigma, alpha_inf, Lambda, Lambda_p),   # +10% porosity
    (phi * 0.9, sigma, alpha_inf, Lambda, Lambda_p),   # -10% porosity
    (phi, sigma * 1.1, alpha_inf, Lambda, Lambda_p),   # +10% flow resistivity
    (phi, sigma * 0.9, alpha_inf, Lambda, Lambda_p),   # -10% flow resistivity
    (phi, sigma, alpha_inf * 1.1, Lambda, Lambda_p),   # +10% tortuosity
    (phi, sigma, alpha_inf * 0.9, Lambda, Lambda_p),   # -10% tortuosity
]

plt.figure(figsize=(8,6))
plt.semilogx(f_meas, alpha_meas,"r", label="Measured")

for i in range(len(sensitivity_params)):
    params = sensitivity_params[i]
    alpha_sens = jca_absorption(params,f_meas)
    plt.semilogx(f_meas,alpha_sens,"--", label=f"Sensitivity {i+1}")

plt.xlabel("Frequency [Hz]")
plt.ylabel("Absorption coefficient")
plt.legend()
plt.grid(True)
plt.xlim(left=100)
plt.savefig("sensitivity_plot.png", dpi=150, bbox_inches="tight")
plt.close()
print("Saved sensitivity plot: sensitivity_plot.png")

#######################################################################
# Step 20 — Save JCA results into the same Excel file

from openpyxl import load_workbook
from openpyxl.drawing.image import Image
from openpyxl.utils.dataframe import dataframe_to_rows
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Recompute model using fitted parameters
alpha_model = jca_absorption(result.x, f_meas)

# ---- Error metrics ----
mse = np.mean((alpha_meas - alpha_model)**2)
rmse = np.sqrt(mse)

ss_res = np.sum((alpha_meas - alpha_model)**2)
ss_tot = np.sum((alpha_meas - np.mean(alpha_meas))**2)
r2 = 1 - (ss_res / ss_tot)

print("\nFit quality:")
print("MSE :", mse)
print("RMSE:", rmse)
print("R²  :", r2)

# ---- Prepare dataframe ----
df_export = pd.DataFrame({
    "Frequency (Hz)": f_meas,
    "Measured absorption": alpha_meas,
    "JCA absorption": alpha_model
})

# ---- Load workbook ----
filename = "my_experiment.xlsx"
book = load_workbook(filename)

# Remove old sheet if exists (optional)
if "JCA_Fit" in book.sheetnames:
    del book["JCA_Fit"]

sheet = book.create_sheet("JCA_Fit")

# ---- Write table ----
for row in dataframe_to_rows(df_export, index=False, header=True):
    sheet.append(row)

# ---- Write fitted parameters ----
params = result.x
labels = [
    "Porosity",
    "Flow resistivity",
    "Tortuosity",
    "Viscous length (m)",
    "Thermal length (m)"
]

sheet["E2"] = "Fitted Parameters"

for i, (label, value) in enumerate(zip(labels, params), start=3):
    sheet[f"E{i}"] = label
    sheet[f"F{i}"] = float(value)

# ---- Write metrics ----
sheet["E9"]  = "Fit Quality"
sheet["E10"] = "MSE"
sheet["F10"] = float(mse)

sheet["E11"] = "RMSE"
sheet["F11"] = float(rmse)

sheet["E12"] = "R²"
sheet["F12"] = float(r2)

# ---- Create plot ----
plt.figure()
plt.semilogx(f_meas, alpha_meas, label="Measured")
plt.semilogx(f_meas, alpha_model, label="JCA model")

plt.xlabel("Frequency (Hz)")
plt.ylabel("Absorption coefficient")
plt.legend()
plt.grid(True)
plt.xlim(left=100)

plot_path = "jca_plot.png"
plt.savefig(plot_path, dpi=150, bbox_inches="tight")
plt.close()

# ---- Insert plot ----
img = Image(plot_path)
sheet.add_image(img, "H2")

# ---- Save workbook ----
book.save(filename)

print("\n✅ Results saved into sheet 'JCA_Fit' in my_experiment.xlsx")



#######################################################################
# End of the script