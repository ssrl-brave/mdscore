import pandas
import numpy as np
import os


def score(test):
    """
    Parameters: `test` is a pandas DataFrame with ids col and alias col, 
    Returns: pd.DataFrame of answers with predictions
    """
    if "ids" in test:
        if not set(test.ids).issubset({0,1,False,True}):
            raise ValueError ("`ids` col must contain only 0,1 or True,False")
        test.ids = test.ids.astype(bool)
    else:
        raise KeyError( "test must contain `ids` col")

    if "alias" not in test:
        raise KeyError( "test must contain `alias` col")

    # download the answers then delete
    if not os.path.exists(".answer.npy"):
        os.system("wget https://smb.slac.stanford.edu/~dermen/answer.npy -O .answer.npy")

    answer = pandas.DataFrame(np.load(".answer.npy", allow_pickle=True)[()], columns=["alias", "bound"])
    os.remove(".answer.npy")
    
    m = pandas.merge(test, answer, on="alias", how='inner')
    bound = m.bound.values.astype(bool)
    ids = m.bound.values.astype(bool)

    print("\nAssuming id=0 means bound:")
    print(pandas.crosstab(m.bound, ~m.ids, colnames=("groundTruth",), rownames=("Prediction",)) )

    print("\nAssuming id=1 means bound:")
    print(pandas.crosstab(m.bound, m.ids, colnames=("groundTruth",), rownames=("Prediction",)) )

    acc0 = (m.ids==~m.bound).sum() / len(m) * 100
    acc1 = (m.ids==m.bound).sum() / len(m) * 100
    print(f"\n{args.test} is {acc0:.1f}% accurate if 0 is bound and {acc1:.1f}% accurate if 1 is bound!")
    return m


if __name__=="__main__":
    from argparse import ArgumentParser
    ap = ArgumentParser()
    ap.add_argument("test", type=str, help="TSV file with columns alias  (names e.g. Azure_Walrus) and ids (0 and 1)")
    args = ap.parse_args()
    test = pandas.read_csv(args.test, sep="\t")
    results = score(test)

