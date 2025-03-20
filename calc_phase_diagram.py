import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy 
from scipy.special import erf


def fit_phase_diagram(csv_path="densities.csv"):
    """ Reads a CSV file with salt, low density, and high density data and performs two regressions:
    
        reg_Y1: (High_Density - Low_Density)^3.06 = d * (1 - S/S_c)
        reg_Y2: High_Density + Low_Density = 2*d_c + 2*A*(S - S_c)
    
    It then extracts the fitted parameters:
        - d (coefficient)
        - A (coefficient)
        - S_c (critical salt)
        - d_c (critical density)
    
    Finally, the function plots the two fits side-by-side and creates a phase diagram highlighting the critical point
    
    """
    df = pd.read_csv(csv_path)

    # Calculate regression variables
    df["reg_Y1"] = (df["High_Density"] - df["Low_Density"])**3.06
    df["reg_Y2"] = df["High_Density"] + df["Low_Density"]
    
    S = df["Salt"].values
    y1 = df["reg_Y1"].values
    y2 = df["reg_Y2"].values

    # Fit function for reg_Y1: y1 = d*(1 - S/S_c)
    fit_func1 = lambda S, d, S_c: d * (1 - (S / S_c))
    params1, cov1 = curve_fit(fit_func1, S, y1, bounds=((-1e15, 0.06), (0.08)))   # the bounds set so that d can be very negative to positive (if needed) and S_c is >0.06.
    d, S_c = params1  

    # Fit function for reg_Y2: y2 = 2*d_c + 2*A*(S - S_c)
    fit_func2 = lambda S, A, d_c: (2 * d_c) + (2 * A * (S - S_c))
    params2, cov2 = curve_fit(fit_func2, S, y2)
    A, d_c = params2
    print(f"Critical salt: {S_c}, Critical density: {d_c}")

    # Regression plots with LaTeX equation titles
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))

    ax[0].plot(S, fit_func1(S, d, S_c), label="Fitted Curve", color='crimson', linewidth=2)
    ax[0].scatter(S, y1, label="Original Data", color='navy', s=50, zorder=3)
    ax[0].set_xlabel("Salt / M", fontsize=12, fontweight="bold")
    ax[0].set_ylabel(r"$\mathrm{reg\_Y1}$", fontsize=12, fontweight="bold")
    ax[0].set_title(r"$ (\rho_h-\rho_l)^{3.06}=d\left(1-\frac{S}{S_c}\right)$", fontsize=12, fontweight="bold")
    ax[0].grid(True, linestyle="--", linewidth=0.8, alpha=0.7)
    ax[0].legend(fontsize=10)
    ax[0].tick_params(axis='both', labelsize=10)

    ax[1].plot(S, fit_func2(S, A, d_c), label="Fitted Curve", color='crimson', linewidth=2)
    ax[1].scatter(S, y2, label="Original Data", color='navy', s=50, zorder=3)
    ax[1].set_xlabel("Salt / M", fontsize=12, fontweight="bold")
    ax[1].set_ylabel(r"$\mathrm{reg\_Y2}$", fontsize=12, fontweight="bold")
    ax[1].set_title(r"$ \rho_h+\rho_l=2\rho_c+2A\,(S-S_c)$", fontsize=12, fontweight="bold")
    ax[1].grid(True, linestyle="--", linewidth=0.8, alpha=0.7)
    ax[1].legend(fontsize=10)
    ax[1].tick_params(axis='both', labelsize=10)

    plt.tight_layout()
    plt.savefig("regression_plots.png", dpi=300, bbox_inches="tight")
    plt.close()

    # Phase diagram with dots only, including the critical point
    fig, ax = plt.subplots(figsize=(6, 5))
    densities = np.concatenate((df["Low_Density"].values, df["High_Density"].values))
    salts_concat = np.concatenate((df["Salt"].values, df["Salt"].values))
    ax.scatter(densities, salts_concat, color='darkgreen', s=60, alpha=0.8, label='Coexistence density values')
    ax.scatter(d_c, S_c, marker='o', color='white', edgecolor='red', s=100, linewidth=2,
            label=f"Critical point, Cc = {round(S_c, 4)}")
    ax.set_xlabel("Density / g cm$^{-3}$", fontsize=12, fontweight="bold")
    ax.set_ylabel("Salt / M", fontsize=12, fontweight="bold")
    ax.set_title("Phase diagram of chromatin\n(linker DNA length: 25 bp)", fontsize=14, fontweight="bold", pad=10)
    ax.legend(fontsize=10)
    ax.grid(True, linestyle="--", linewidth=0.8, alpha=0.7)

    plt.tight_layout()
    plt.savefig("phase_diagram.png", dpi=300, bbox_inches="tight")
    plt.close()


    return d, A, S_c, d_c


def main():
    d, A, S_c, d_c = fit_phase_diagram("densities.csv")
    print(f"Fitted parameters:\nd = {d}\nA = {A}\nCritical salt (S_c) = {S_c}\nCritical density (d_c) = {d_c}")

if __name__ == "__main__":
    main()