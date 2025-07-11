import pandas
import numpy as np
from mdscore import score
from mpi4py import MPI
COMM = MPI.COMM_WORLD

results = {}
for dset_name in ("zinc", "ucsf"):
    df = pandas.read_csv(f"/Users/dermen/ensemble/{dset_name}_boltz.tsv", sep="\t")
    rdk_col = df.loc[:, "MaxAbsEStateIndex":"fr_urea"].columns.tolist()
    res_cols = [c for c in df if c.startswith("resid_")]
    md_col = list(df)[1:18] + res_cols + ["max_state_t"]

    results[dset_name] = []
    res_ids = 154, 21, 153, 152, 50
    feats_to_search= []
    for res_id in res_ids:
        for metric in ("ave", "std", "min" , "max"):
            feats = [f"resid_{res_id}_{metric}_md{x}" for x in range(3)]
            new_feat = f"resid_{res_id}_{metric}"
            df[new_feat] = df[feats].values.mean(axis=1)
            feats_to_search.append(new_feat)

    feats_to_search  += md_col + rdk_col
    for i_n, name in enumerate(feats_to_search):
        if i_n % COMM.size != COMM.rank:
            continue
        vals = df[name]
        feat_results = []
        for invert in [True, False]:

            # TODO: choose better cutoff here, 50 is arbitrary
            for v in np.linspace( vals.min(), vals.max(), 50):
                if invert:
                    f1 = score.prc(df.bound, ~(vals > v), False)["F1"]
                else:
                    f1 = score.prc(df.bound, vals > v, False)["F1"]
                result = (f1, invert, name)
                feat_results.append(result)
        best = sorted(feat_results)[-1][0]
        print(f"For {dset_name} feat {name} (invert={invert}), top hit gave F1={best:.3f} ({i_n+1} / {len(feats_to_search)})")
        results[dset_name] += feat_results

final_results = {}
for dset_name in results:
    vals = COMM.reduce(results[dset_name])
    if COMM.rank==0:
        final_results[dset_name]= vals

if COMM.rank==0:
    np.save("final_results", final_results)

