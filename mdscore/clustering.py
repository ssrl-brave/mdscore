from argparse import ArgumentParser
ap = ArgumentParser()
ap.add_argument("--nclust", default=4, type=int)
args = ap.parse_args()

from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from mdscore import score
from itertools import combinations
import numpy as np
import pandas

results = {}
for dset_name in ("ucsf", "zinc"):
    df = pandas.read_csv(f"/Users/dermen/ensemble/{dset_name}_boltz.tsv", sep="\t")
    rdk_col = df.loc[:, "MaxAbsEStateIndex":"fr_urea"].columns.tolist()

    res_cols = [c for c in df if c.startswith("resid_")]
    md_col = list(df)[1:18] + res_cols + ["max_state_t"]
    results[dset_name] = {}

    S = StandardScaler()

    for name, cols in ( ("rdk", rdk_col), ("md", md_col), ("both", rdk_col+md_col)):
        results[dset_name][name] = []

        dat = df[cols].values

        sdat = S.fit_transform(dat)

        n_clust = args.nclust
        K = KMeans(n_clusters=n_clust)
        K.fit(sdat)
        Klab = K.labels_

        combos = [s for sl in [list(combinations( range(n_clust), r=r)) for r in range(1,n_clust)] for s in sl]
        for c in combos:
            pred_bound = np.zeros_like(df.bound)
            for i in c:
                pred_bound[Klab == i] = True
            result = score.prc(df.bound, pred_bound, False)
            results[dset_name][name].append( (result["F1"], c))

for dset_name in results:
    print(dset_name, ":")
    for name in results[dset_name]:
        best = sorted(results[dset_name][name])[-1][0]
        print(name, best)
    print()



