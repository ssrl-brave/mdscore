import pandas
import numpy as np
import os

from sklearn.metrics import precision_score, recall_score, f1_score


def score(test, bound_label=0):
    """
    Parameters: `test` is a pandas DataFrame with ids col and alias col, 
    Returns: pd.DataFrame of answers with predictions
    """
    assert "ids" in test
    assert "alias" in test

    # download the answers then delete
    if not os.path.exists(".answer.npy"):
        os.system("wget https://smb.slac.stanford.edu/~dermen/answer.npy -O .answer.npy")

    answer = pandas.DataFrame(np.load(".answer.npy", allow_pickle=True)[()], columns=["alias", "bound"])
    answer.bound = answer.bound.astype(bool)
    
    m = pandas.merge(test, answer, on="alias", how='inner')
    
    actual_bound = m.bound
    pred_bound = m["ids"]==bound_label
    results = prc(actual_bound, pred_bound)
    return m, results


def prc(actual_bound, pred_bound):
    p = precision_score(actual_bound, pred_bound)
    r = recall_score(actual_bound, pred_bound)
    f = f1_score(actual_bound, pred_bound)

    actual_neg = np.logical_not(actual_bound)
    pred_neg = np.logical_not(pred_bound)

    TP = (pred_bound * actual_bound).sum()
    FP = (pred_bound * actual_neg).sum()
    assert p== TP / (FP + TP)

    FN = (pred_neg * actual_bound).sum()
    assert r==TP / (FN + TP)
    TN = (pred_neg * actual_neg).sum()

    assert np.allclose(f,2*(p*r)/(p+r))
    print(f"precision: {p}")
    print(f"recall: {r}")
    print(f"F1-score: {f}")
    print(f"TP={TP}, FP={FP}, FN={FN}, TN={TN}")
    results = {"precision":p, "recall": r, "F1": f}

    return results


if __name__=="__main__":
    from argparse import ArgumentParser
    ap = ArgumentParser()
    ap.add_argument("test", type=str, help="TSV file with columns alias  (names e.g. Azure_Walrus) and ids (0 and 1)")
    args = ap.parse_args()
    test = pandas.read_csv(args.test, sep="\t")
    _,results = score(test)

