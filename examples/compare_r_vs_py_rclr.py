from pathlib import Path
import numpy as np
import pandas as pd
import os

interop = Path("interop")
samples = pd.read_csv(interop / "samples.csv")["sample"].astype(str).tolist()

def load_pair(i):
    R = pd.read_csv(interop / f"rclr_R/view_{i}_rclr_R.csv", index_col = 0)
    P = pd.read_csv(interop / f"rclr_PY/view_{i}_rclr_PY.csv", index_col = 0)

    #align rows (features) and columns (samples)
    rows = R.index.intersection(P.index)
    cols = R.columns.intersection(P.columns)
    R = R.loc[rows, cols].sort_index(axis = 0).sort_index(axis = 1)
    P = P.loc[rows, cols].sort_index(axis = 0).sort_index(axis = 1)
    return R, P

def summarize(name, A, B):
    x = A.values.ravel()
    y = B.values.ravel()
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if x.size == 0:
        print(f"{name}: no finite overlap")
        return
    # linear fit y ≈ a + b x
    X = np.vstack([np.ones_like(x), x]).T
    a, b = np.linalg.lstsq(X, y, rcond = None)[0]
    r = np.corrcoef(x, y)[0, 1]
    resid = y - (a + b * x)
    print(f"{name}: corr = {r:.4f}, slope = {b:.4f}, intercept = {a:.4f}, max|resid| = {np.max(np.abs(resid)):.3e}, n = {x.size}")

#discover how many views we have
view_files = sorted([f for f in os.listdir(interop / "rclr_R") if f.endswith("_rclr_R.csv")])
nviews = len(view_files)

for i in range(1, nviews + 1):
    R, P = load_pair(i)
    summarize(f"view_{i}", R, P)
