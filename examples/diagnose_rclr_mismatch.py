from pathlib import Path
import numpy as np, pandas as pd

interop = Path("interop")
samples = pd.read_csv(interop/"samples.csv")["sample"].astype(str).tolist()

def load_pair(i):
    R = pd.read_csv(interop/f"rclr_R/view_{i}_rclr_R.csv", index_col = 0)
    P = pd.read_csv(interop/f"rclr_PY_view_{i}.csv", index_col = 0)
    rows = R.index.intersection(P.index)
    cols = R.columns.intersection(P.columns)
    #align to the same rows/cols, sorted
    R = R.loc[rows, cols].sort_index(axis = 0).sort_index(axis = 1)
    P = P.loc[rows, cols].sort_index(axis = 0).sort_index(axis = 1)
    return R, P

def summarize(name, A, B):
    #flatten, then regression: B ≈ a + b*A
    x = A.to_numpy().ravel()
    y = B.to_numpy().ravel()
    ok = np.isfinite(x) & np.isfinite(y)
    x = x[ok]; y = y[ok]
    X = np.vstack([np.ones_like(x), x]).T
    a, b = np.linalg.lstsq(X, y, rcond = None)[0]  #intercept, slope
    corr = np.corrcoef(x, y)[0, 1]
    resid = y - (a + b*x)
    print(f"{name}: corr = {corr:.4f}, slope = {b:.6f}, intercept = {a:.6f}, max|resid| = {np.max(np.abs(resid)):.3e}")

for i in (1, 2, 3):
    R, P = load_pair(i)
    summarize(f"view_{i} (raw)", R, P)
    #hypothesis 1: R used log10; Gemelli uses ln
    summarize(f"view_{i} (R*ln10)", R*np.log(10), P)
    #hypothesis 2: sample-wise robust centering mismatch
    Rc = R.sub(R.median(axis = 0), axis = 1)
    Pc = P.sub(P.median(axis = 0), axis = 1)
    summarize(f"view_{i} (centered)", Rc, Pc)
    #hypothesis 3: both log-base and centering
    Rln = R*np.log(10)
    Rc2 = Rln.sub(Rln.median(axis = 0), axis = 1)
    Pc2 = P.sub(P.median(axis = 0), axis = 1)
    summarize(f"view_{i} (R*ln10 & centered)", Rc2, Pc2)
