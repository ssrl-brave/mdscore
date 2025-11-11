import glob
import pandas
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA


from mdscore import state_times

# load the known binding labels 
ucsf_answer = pandas.read_csv("answer.tsv", sep="\t")

# visual inspection indicated 20 components would be sufficient
pc = PCA(n_components=20)
# mean sub, variance division
scaler = StandardScaler()

for dset_name in ("zinc", "ucsf"):
    dirname=f'data/{dset_name}/final_data/'
    for trial in ("diffdock", "boltz"):
        
        df = pandas.read_csv(f"{dirname}/final_ml_training_data_ucsf{'_boltz' if trial=='boltz' else ''}.csv", sep=",")
        fnames = glob.glob(f"{dirname}distance*{'boltz' if trial == 'boltz' else '0'}.npy")
        for f in fnames:
            d = np.load(f)[()]
            for i_sim, d_ in enumerate((d[:,:5000], d[:,5000:10000], d[:,-5000:]) ):
              sd = scaler.fit_transform(d_)
              pd = pc.fit_transform(sd)
              resname = f.split("distance_")[1].split("_457")[0]
              # store the princ components as separate cols
              for i_comp in range(pd.shape[1]):
                pname = f"{resname}_pc{i_comp}_md{i_sim}"
                df[pname] = pd[:,i_comp]
            
              min_d, max_d, ave_d, std_d = d_.min(1), d_.max(1), d_.mean(1), d_.std(1)
              for name, val in [('min', min_d),('max', max_d),('ave', ave_d),('std', std_d)]:
                df[f"{resname}_{name}_md{i_sim}"] = val
            print(resname)
        # weird fragmentation warning.. seems harmless, but this clears it up:
        df = df.copy()

        # load the clusterID vs frame
        all_state_t = state_times.get_state_t(f"{dirname}/cluster_labels_457x15000{'_boltz' if trial=='boltz' else ''}.npy")
        max_state_t = []
        for system in all_state_t:
            st = all_state_t[system]
            if st:
                # whoops, looks like sum of the total residence time
                t = sum(st[0][1])
            else:
                t = 0
            max_state_t.append(t)
        df["max_state_t"] = max_state_t # this isnt max, its sum 

        # merge the answers with the features
        if dset_name=="ucsf":
            M = pandas.merge(df, ucsf_answer, on="System")
        else:
            M = df
            M["bound"] = M["Fragment?"]

        M.to_csv(f"{dset_name}_{trial}.tsv", sep="\t", index=False)
        print(f"Done with {dset_name}_{trial} ... \n ")

