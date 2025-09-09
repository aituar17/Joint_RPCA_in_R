from pathlib import Path
import numpy as np, pandas as pd

interop = Path("interop")
samples = pd.read_csv(interop/"samples.csv")["sample"].astype(str).tolist()

def load_pair(i):
    R = pd.read_csv(interop/f"rclr_R/view_{i}_rclr_R.csv", index_col = 0)
    P = pd.read_csv(interop/f"rclr_PY_view_{i}.csv", index_col = 0)
    rows = R.index.intersection(P.index)
    cols = R.columns.intersection(P.columns)
    R = R.loc[rows, cols].sort_index(0).sort_index(1)
    P = P.loc[rows, cols].sort_index(0).sort_index(1)
    return R, P

def summarize(name, A, B):
    #flatten, then regression B ≈ a + b*A
    x = A.values.ravel()
    y = B.values.ravel()
    ok = np.isfinite(x) & np.isfinite(y)
    x = x[ok]; y = y[ok]
    X = np.vstack([np.ones_like(x), x]).T
    coef, *_ = np.linalg.lstsq(X, y, rcond = None)  
    a, b = coef
    corr = np.corrcoef(x, y)[0,1]
    print(f"{name}: corr = {corr:.4f}, slope = {b:.4f}, intercept = {a:.4f}, "
          f"max|diff| = {np.max(np.abs(y - (a + b*x))):.3e}")

for i in (1,2,3):
    R, P = load_pair(i)
    summarize(f"view_{i} (raw)", R, P)
    #hypothesis 1: R used log10, Gemelli uses ln → multiply R by ln(10)
    summarize(f"view_{i} (R*ln10)", R*np.log(10), P)
    #hypothesis 2: different sample-wise centering → recenter per sample by subtracting column medians
    Rc = R - R.median(axis = 0)
    Pc = P - P.median(axis = 0)
    summarize(f"view_{i} (centered)", Rc, Pc)
    #hypothesis 3: both log-base and centering
    Rc2 = (R*np.log(10)) - (R*np.log(10)).median(axis = 0)
    Pc2 = P - P.median(axis = 0)
    summarize(f"view_{i} (R*ln10 & centered)", Rc2, Pc2)
