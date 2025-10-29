import os
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
from scipy import stats

rank="."
# Define main directory
main_dir='.'
bases = sorted([os.path.join(main_dir, d) for d in os.listdir(main_dir) if d.startswith("Mpro") and os.path.isdir(os.path.join(main_dir, d))])

#main_dir2='.'
#bases2 = sorted([os.path.join(main_dir2, d) for d in os.listdir(main_dir2) if d.startswith("system") and os.path.isdir(os.path.join(main_dir2, d))])

#bases=bases1+bases2

# Function to extract the first numerical value from a file, skipping comments
def extract_value(filepath, col=1):
    try:
        with open(filepath, "r") as f:
            for line in f:
                if not line.startswith("#"):  # Ignore commented header lines
                    return float(line.split()[col])  # Read first valid data row
        return np.nan  # Return NaN if no valid data is found
    except Exception:
        return np.nan

# Function to get number of lines in a file
def count_lines(filepath):
    try:
        with open(filepath, "r") as f:
            return sum(1 for _ in f)
    except Exception:
        return np.nan

# Function to compute molecular descriptors from an SDF file
def get_descriptors(sdf_file):
    descriptor_names = [desc[0] for desc in Descriptors._descList]
    supplier = Chem.SDMolSupplier(sdf_file)
    for mol in supplier:
        if mol is not None:
            return {name: getattr(Descriptors, name)(mol) for name in descriptor_names}
    return {name: np.nan for name in descriptor_names}  # Return NaN if no valid molecules

# Initialize DataFrame
columns = ["System", "Distance", "Avg_Energy", "Std_Energy", "Hbond", "RMSD_Top", "Population", "Num_Clusters", "RMSD_Mean", "RMSD_Mode", "RMSD_Std", "RMSF_Mean", "RMSF_Mode", "RMSF_Std","Water_Mean", "Water_Mode", "Water_Std","Boltz_energy","Boltz_distance","Diffdock_distance"]
columns += [desc[0] for desc in Descriptors._descList]  # Add RDKit descriptor headers
df = pd.DataFrame(columns=columns)

# Process each base directory
for base in bases:
    #for rank in ranks:
    row = {"System": base }
    base_dir = os.path.join('.', base)
    print(base_dir)    
        # Extract confidence value from filename
    sdf_file="ligand.sdf"
    #sdf_file = next((f for f in os.listdir(base_dir) if f.startswith(rank[:-4] + "_confidence") and f.endswith(".sdf")), None)
    if sdf_file:
        print(sdf_file)
         #   row["Confidence"] = float(sdf_file.split("confidence")[-1].split(".sdf")[0])
        row.update(get_descriptors(os.path.join(base_dir, sdf_file)))  # Get descriptors
    else:
         #   row["Confidence"] = np.nan
        row.update({desc: np.nan for desc in Descriptors._descList})
        #print(sdf_file,row["Confidence"]) 
        # Extract numerical values from text files
    row["Distance"] = extract_value(os.path.join(base_dir,rank, "avg_distances_boltz.txt")) *10.0
    row["Avg_Energy"] = extract_value(os.path.join(base_dir,rank, "avg_energies_boltz.dat"))
    row["Std_Energy"] = extract_value(os.path.join(base_dir, rank,"std_energies_boltz.dat"))
    row["Hbond"] = extract_value(os.path.join(base_dir,rank, "hbond_counts_boltz.txt"))
    row["RMSD_Top"] = extract_value(os.path.join(base_dir, rank,"lig_rmsd_top5_boltz.out"))
        
    # Clustering summary
    cluster_file = os.path.join(base_dir, rank,"Clustering_boltz", "summary")
    row["Population"] = extract_value(cluster_file, col=2)
    row["Num_Clusters"] = count_lines(cluster_file)-1
        
##RMSD calculations
    rmsd_values = []
    for trial in ["trial_1_boltz", "trial_2_boltz", "trial_3_boltz"]:
        rmsd_file = os.path.join(base_dir, rank, trial, "lig_rmsd_boltz.out")
        try:
            with open(rmsd_file, "r") as f:
                for line in f:
                    if line.startswith("#"):  # Skip header lines
                        continue
                    rmsd_values.append(float(line.split()[1]))
        except Exception:
            pass
    len(rmsd_values) 
    if rmsd_values:
        row["RMSD_Mean"] = np.mean(rmsd_values)
        row["RMSD_Mode"] = stats.mode(rmsd_values,keepdims=True)[0][0]
        row["RMSD_Std"] = np.std(rmsd_values)
    else:
        row["RMSD_Mean"] = row["RMSD_Mode"] = row["RMSD_Std"] = np.nan

## RMSF of lig
    rmsf_values=[]
    for trial in ["trial_1_boltz", "trial_2_boltz", "trial_3_boltz"]:
        rmsf_file = os.path.join(base_dir, rank, trial, "rmsf_boltz.out")
        try:
            with open(rmsf_file, "r") as f:
                for line in f:
                    if line.startswith("#"):  # Skip header lines
                        continue
                    rmsf_values.append(float(line.split()[1]))
        except Exception:
            pass    
    if rmsf_values:
        row["RMSF_Mean"] = np.mean(rmsf_values)
        row["RMSF_Mode"] = stats.mode(rmsf_values,keepdims=True)[0][0]
        row["RMSF_Std"] = np.std(rmsf_values)
    else:
        row["RMSF_Mean"] = row["RMSF_Mode"] = row["RMSF_Std"] = np.nan

## Number of active site water
    water_values=[]
    water_files= os.path.join(base_dir,rank, "water_3_boltz.dat")
    try:
        with open(water_files) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                water_values.append(float(line.split()[1]))
    except Exception:
        pass
    if water_values:
        row["Water_Mean"] = np.mean(water_values)
        row["Water_Mode"] = stats.mode(water_values,keepdims=True)[0][0]
        row["Water_Std"] = np.std(water_values)
    else:
        row["Water_Mean"] = row["Water_Mode"] = row["Water_Std"] = np.nan
    row["Boltz_energy"] = extract_value(os.path.join(base_dir, "boltz_energies.dat"))
    row["Boltz_distance"] = extract_value(os.path.join(base_dir,"boltz_distances.txt"))*10.0
    row["Diffdock_distance"] = extract_value(os.path.join(base_dir,rank, "diffdock_distances.txt"))*10.0
    #if base_dir.split('/')[-2]=='diff_dock_misses20':
     #   row["Class"] =0
    #else:
     #   row["Class"] =1
### Append to df
    df = pd.concat([df, pd.DataFrame([row])], ignore_index=True)

# Save final CSV
output_file = "../final_ml_training_data_mpro_boltz.csv"
df.to_csv(output_file, index=False)
print(f"Data saved to {output_file}")
