import pandas
import numpy as np
#from mdscore import score
from sklearn.metrics import f1_score
from mpi4py import MPI
COMM = MPI.COMM_WORLD

results = {}
df1 = pandas.read_csv(f"ucsf_boltz.tsv", sep="\t")
df2 = pandas.read_csv(f"zinc_boltz.tsv", sep="\t")
df = pandas.concat((df1, df2)).reset_index(drop=True)
rdk_col = df.loc[:, "MaxAbsEStateIndex":"fr_urea"].columns.tolist()
res_cols = [c for c in df if c.startswith("resid_")]
md_col = list(df)[1:18] + res_cols + ["max_state_t"]

results = []
res_ids = 154, 21, 153, 152, 50
feats_to_search= []
for res_id in res_ids:
    for metric in ("ave", "std", "min" , "max"):
        feats = [f"resid_{res_id}_{metric}_md{x}" for x in range(3)]
        new_feat = f"resid_{res_id}_{metric}"
        df[new_feat] = df[feats].values.mean(axis=1)
        feats_to_search.append(new_feat)

feats_to_search += md_col + rdk_col
feat_type = np.array(["md"]*len(feats_to_search))
feat_type[-len(rdk_col):] = "rd"
for i_n, name in enumerate(feats_to_search):
    if i_n % COMM.size != COMM.rank:
        continue
    ft = feat_type[i_n]
    vals = df[name]
    feat_results = []
    for invert in [True, False]:

        # TODO: choose better cutoff here, 50 is arbitrary
        for v in np.linspace( vals.min(), vals.max(), 100):
            if invert:
                #f1 = score.prc(df.bound, ~(vals > v), False)["F1"]
                f1 = f1_score(df.bound, ~(vals > v))
            else:
                #f1 = score.prc(df.bound, vals > v, False)["F1"]
                f1 = f1_score(df.bound, vals > v)
            result = (f1, invert, name, v, ft)
            feat_results.append(result)
    best = sorted(feat_results)[-1]
    best_f1 = best[0]
    print(f"For feat {name} (invert={invert}), top hit gave F1={best_f1:.3f} ({i_n+1} / {len(feats_to_search)})")
    results.append(best)

results = COMM.reduce(results)

if COMM.rank==0:
    np.save("all_results_both", results)
    results = sorted(results)[::-1]
    Rmd = [r for r in results if r[4] == "md"]
    Rrdkit = [r for r in results if r[4] == "rd"]

    import tabulate
    print(tabulate.tabulate(Rrdkit[:5] + Rmd[:5],
                            headers=("F1-score", "feature", "inverted", "cutoff", "feauture-type")))
    print()

