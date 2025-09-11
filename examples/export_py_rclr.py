from pathlib import Path
import os
import numpy as np
import pandas as pd
from biom import Table
from gemelli.preprocessing import matrix_rclr

interop = Path("interop")
outdir = interop / "rclr_PY"
outdir.mkdir(parents = True, exist_ok = True)

samples = pd.read_csv(interop / "samples.csv")["sample"].astype(str).tolist()
view_files = sorted([f for f in os.listdir(interop) if f.endswith("_counts.csv")])

for i, vf in enumerate(view_files, start = 1):
    #load counts as features x samples
    df = pd.read_csv(interop / vf, index_col = 0)
    #enforce sample order and numeric
    df = df.loc[:, samples].apply(pd.to_numeric, errors = "coerce")

    #gemelli rCLR expects samples x features → transpose first
    X = df.values.T  #shape: samples x features

    #handle NaNs like you did in R (turn NaNs to 0 before closure/log)
    if np.isnan(X).any():
        X = np.nan_to_num(X, nan = 0.0)

    #run rCLR
    X_rclr = matrix_rclr(X)  #samples x features
    M = pd.DataFrame(X_rclr.T, index = df.index, columns = df.columns)  

    #apply the SAME row filter as in R: finite across samples + sd > 0
    vals = M.values
    ok_finite = np.isfinite(vals).all(axis = 1)
    sds = M.std(axis = 1).to_numpy()
    keep = ok_finite & np.isfinite(sds) & (sds > 0)
    M = M.loc[keep]

    #write
    M.to_csv(outdir / f"view_{i}_rclr_PY.csv")
    print(f"[ok] view_{i}: rclr_PY shape = {M.shape}")
