import os
import pandas as pd
import glob
# Define main directory
main_dir = "."  # Update this with your actual directory
#ranks = glob.glob('system_*')
#print(ranks)
mmgbsa_folders = ["mmgbsa_0_boltz"]

# Function to extract mean and standard deviation from the FINAL_RESULTS file
def extract_values(filepath):
    try:
        with open(filepath, "r") as f:
            lines = f.readlines()
            line = lines[86]  # Extracting 88th line (index 87)
            avg = float(line.split()[2])  # Extract mean value
            std = float(line.split()[3])  # Extract standard deviation value
            return avg, std
    except Exception:
        return "N/A", "N/A"

# Loop through systems in main directory
systems = sorted([d for d in os.listdir(main_dir) if d.startswith("Mpro") and os.path.isdir(os.path.join(main_dir, d))])

for system in systems:
    print(system)
    system_dir = os.path.join(main_dir, system)
    
    #for rank in ranks:
    avg_data = []  # Store avg values
    std_data = []  # Store std values
    #rank_dir = os.path.join(system_dir, rank)  # Rank directory
    os.makedirs(system_dir, exist_ok=True)  # Ensure directory exists
        
    for folder in mmgbsa_folders:
        file_path = os.path.join(system_dir, folder, "FINAL_RESULTS_MMPBSA.dat")
        avg, std = extract_values(file_path)
        avg_data.append(f"{folder} {avg}")
        std_data.append(f"{folder} {std}")
        
        # Save as .dat files inside the rank directories
    with open(os.path.join(system_dir, "avg_energies_boltz.dat"), "w") as f:
        f.write("\n".join(avg_data) + "\n")
        
    with open(os.path.join(system_dir, "std_energies_boltz.dat"), "w") as f:
        f.write("\n".join(std_data) + "\n")

