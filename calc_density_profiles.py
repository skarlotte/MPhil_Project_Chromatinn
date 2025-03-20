import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy 
from scipy.special import erf


def get_frames(file_path):
    """Reads a .dump file and extracts the nucleosome z-coordinates across simulation frames
    Inputs:
        - file_path (str): path to the .dump file
    Outpus:
        - frames (lis): a list of NumPy arrays, array representing the z-coordinates of nucleosomes in a simulation frame

    """
    frame = []             
    frames = []            
    reading_atoms = False  
    filtered_types = {1}  

    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith("ITEM: TIMESTEP"):
                if frame:                          
                    frames.append(np.array(frame))
                frame = []                         
                continue
            elif line.startswith("ITEM: ATOMS"):
                reading_atoms = True
                continue
            elif line.startswith("ITEM:"):
                reading_atoms = False  

            if reading_atoms:  
                columns = line.split()
                if len(columns) < 4:  
                    continue  
                try:
                    atom_type = int(columns[-1])  
                    if atom_type in filtered_types:
                        frame.append(float(columns[3]))  
                except ValueError:
                    print(f"Warning: Skipping invalid line in {file_path}: {line.strip()}")
    if frame:  
        frames.append(np.array(frame))  

    print(f"The trajectory {file_path} has been read. Total frames read: {len(frames)}")
    return frames

def apply_pbc(z, z_min, z_max):
    """ Applies periodic boundary conditions (PBC) to the z-coordinates of the nucleosomes 
        by wrapping the coordinates when they exceed the boundaries.
    Inputs:
        - z (array-like): array of the z-coordinates of the nucleosomes
        - z_min (float): minimum z coordinate of the box
        - z_max (float): maximum z coordinate of the box
    Returns:
        - z (np.ndarray): z-coordinates of the nucleosome adjusted to periodic boundary conditions
    """
    box_length = z_max - z_min  
    z = np.array(z)  
    for i in range(len(z)):
        while z[i] > z_max:
            z[i] -= box_length
        while z[i] < z_min:
            z[i] += box_length
    return z

def center_pulse(z, z_min, z_max):
    """ Shifts z-coordinates to the center of the simulation box
        using a circular mean approach.
    Inputs:
        - z (array-like): array of the z-coordinates of the nucleosomes
        - z_min (float): minimum z coordinate of the box
        - z_max (float): maximum z coordinate of the box
    Outputs:
        - z (np.ndarray): centered z-coordinates of the nucleosome
    """
    box_length = z_max - z_min  
    C = (z_max + z_min) / 2  
    zangs = z * 2 * np.pi / box_length     
    meanangz = np.arctan2(np.mean(np.sin(zangs)), np.mean(np.cos(zangs)))     
    meanz = meanangz * box_length / (2 * np.pi)     
    delta = C - meanz     
    z = z + delta
    return z

def calc_dprof(frames, z_min, z_max, nbins, equillim, save_path):
    """ Bins nucleosome positions into a histogram and averages over frames,
        computing the density profile of nucleosomes along the z-axis.
    Inputs:
        - frames (list): z-coordinates of nucleosomes per simulation frame
        - z_min (float): minimum z coordinate of the box
        - z_max (float): maximum z coordinate of the box
        - nbins (int): number of bins for the density histogram
        - equillim (int): number of frames to discard for equilibration
        - save_path (str): file path to save the computed density profile

    Returns:
        - z_bins: center positions of histogram bins
        - dprof: density profile at positions z_bins, averaged over the frames
        - dprof_std: standard deviation of the density profile
    """
    if len(frames) < equillim:
        print("Warning: Not enough frames for analysis. Proceeding with all available frames")
    else:
        frames = frames[equillim:]

    dprof_per_frame = []
    for z in frames: 
        z = center_pulse(z, z_min, z_max)
        z = apply_pbc(z, z_min, z_max)
        y, b = np.histogram(z, bins=nbins, range=(z_min, z_max))
        z_bins = (b[1:] + b[:-1]) / 2
        dprof_per_frame.append(y)
        
    if len(dprof_per_frame) == 0:
        print("Error: No valid nucleosome data found in dump file")
        exit()
    
    dprof = np.mean(np.array(dprof_per_frame), axis=0)
    dprof_std = np.std(np.array(dprof_per_frame), axis=0)
    np.savetxt(save_path, np.array([z_bins, dprof, dprof_std]).T) 

    return z_bins, dprof, dprof_std

def rho(z, rho_a, rho_b, z_DS, t):
    """Hyperbolic tangent function for fitting the density profile
    Arguments:
        - z (array-like): z-coordinates
        - rho_a (float): low-density value
        - rho_b (float): high-density value
        - z_DS (float): transition region center
        - t (float): transition width.

    Returns:
        - fitted density values in an np.array
    """
    return (rho_a + rho_b) / 2 + (rho_b - rho_a) / 2 * np.tanh((np.abs(z) - z_DS) / t)

def calc_high_low_densities(z_bins, dprof):
    """Fit a hyperbolic tangent function to extract low/high-density values.

    Inputs:
        - z_bins (array-like): binned z-coordinates
        - dprof (array-like): density profile

    Outputs:
        - ld_density: low-density value
        - hd_density: high-density value
        - popt: fitted parameters
    """
    initial_params = [dprof.min(), dprof.max(), np.median(z_bins), -1]  
    popt, _ = curve_fit(rho, z_bins, dprof, p0=initial_params)     
    rho_a, rho_b, z_DS, t = popt                         
    ld_density, hd_density = sorted([rho_a, rho_b])      
    if ld_density < 0:
        ld_density = 1e-3    

    return ld_density, hd_density, popt

def plot_dprof(z_bins, dprof, dprof_std, popt, savefig_path, salt, hd_density):
    """Plot the density profile of nucleosomes along the z-axis with enhanced aesthetics
    Inputs:
        - z_bins (array-like): binned z-coordinates
        - dprof: density profile at positions z_bins, averaged over the frames
        - dprof_std: standard deviation of the density profile
        - popt: fitted parameters
        - savefig_path: path of directory to save the plots
        - salt: the salt values of the simulations
        - hd_density: the high density value of the densitry profile  
    
    """
    
    plt.figure(figsize=(8, 4.5), dpi=150) 
    plt.plot(z_bins, dprof, color="#E63946", linewidth=2, label="Density")
    upper_bound_std = dprof + dprof_std
    lower_bound_std = np.maximum(dprof - dprof_std, 0)  # Ensure no negative values
    plt.fill_between(z_bins, lower_bound_std, upper_bound_std, color="#E63946", alpha=0.25, label="Density ± Std")
    plt.plot(z_bins, rho(z_bins, *popt), 'b--', linewidth=2, label=f"Fitted Profile, hd: {np.round(hd_density, 2)}")
    plt.xlabel("Z Coordinate (Å)", fontsize=12, fontweight="bold")
    plt.ylabel("Density", fontsize=12, fontweight="bold")
    plt.title(f"Nucleosome Density Profile | Salt: {salt} M", fontsize=14, fontweight="bold", pad=10)
    plt.grid(True, linestyle="--", linewidth=0.6, alpha=0.6)
    plt.legend(frameon=True, fontsize=10, loc="upper right")
    plt.savefig(savefig_path, dpi=200, bbox_inches="tight")

def process_dumps(dumps_dir_name):
    """Process .dump files in the given directory and extract nucleosome z-coordinates.
    Inputs:
        - dumps_dir_name (str): name of the directory containing .dump files

    Returns:
        - frames_per_salt: list of extracted frames per salt concentration
        - salts: list of corresponding salt concentrations
    """
    current_dir = os.getcwd()
    dumps_dir = os.path.join(current_dir, dumps_dir_name)
    dump_paths = sorted(glob.glob(os.path.join(dumps_dir, "*.dump")))  
    
    if not dump_paths:
        raise FileNotFoundError(f"No .dump files found in {dumps_dir_name}")

    salts, frames_per_salt = [], []
    for dump_path in dump_paths:
        filename = os.path.basename(dump_path)
        parts = filename.replace(".dump", "").split("_")
        if len(parts) > 1:
            try:
                salt = float(parts[1])  
                salts.append(salt)
                frames = get_frames(dump_path)  
                frames_per_salt.append(frames)
            except ValueError:
                print(f"Warning: Skipping invalid filename {filename}")

    return frames_per_salt, salts

def get_dprof(frames_per_salt, salts, nbins, equillim, z_min, z_max):
    """Compute and save density profiles. The salt concentration, low density and high density values are saved in a csv file, 
     and the density fits are saved as a png. """
    
    current_dir = os.getcwd()
    output_dir = os.path.join(current_dir, "density_profile_and_fits")
    os.makedirs(output_dir, exist_ok=True)

    results = {"Salt": [], "Low_Density": [], "High_Density": []}
    for frames, salt in zip(frames_per_salt, salts):
        print(f"Processing salt concentration: {salt}")
        
        savefig_path = os.path.join(output_dir, f"{salt}_density_prof.png")
        savetxt_path = os.path.join(output_dir, f"{salt}_density_prof.txt")
        
        z_bins, dprof, dprof_std = calc_dprof(frames, z_min, z_max, nbins, equillim, savetxt_path)
        ld_density, hd_density, popt = calc_high_low_densities(z_bins, dprof)
        plot_dprof(z_bins, dprof, dprof_std, popt, savefig_path, salt, hd_density)
        np.savetxt(savetxt_path, np.column_stack([z_bins, dprof, dprof_std]))
        
        results["Salt"].append(salt)
        results["Low_Density"].append(ld_density)
        results["High_Density"].append(hd_density)
    
    csv_path = os.path.join(current_dir, "densities.csv")
    pd.DataFrame(results).to_csv(csv_path, index=False)
    print(f"Low and high densities saved to {csv_path}")


# Parameters
nbins = 900      # number of bins for histogram
equillim = 10   # calculate the density profile after EQUILLIM amount of frames for equilibration
z_min = -3500    # minimum z coordinate of the simulation box
z_max = 3500     # maximum z coordinate of the simulation box
dprof_dir_name = "dumps_cores"

def main():
    coords_salt, salts = process_dumps(dprof_dir_name)
    get_dprof(coords_salt, salts, nbins, equillim, z_min, z_max)

if __name__ == "__main__":
    main()

