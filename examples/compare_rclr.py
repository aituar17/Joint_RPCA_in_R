from pathlib import Path
import pandas as pd, numpy as np, os
from biom import Table
from gemelli.preprocessing import matrix_rclr

interop = Path("interop")

view_files = sorted([f for f in os.listdir(interop) if f.endswith("_counts.csv")])
samples = pd.read_csv(interop/"samples.csv")["sample"].astype(str).tolist()

py_rclr = []
for i, vf in enumerate(view_files, start = 1):
    df = pd.read_csv(interop/vf, index_col = 0)

    #enforce sample order and subset to declared samples
    df.columns = df.columns.astype(str)
    missing = [s for s in samples if s not in df.columns]
    if missing:
        raise ValueError(f"{vf} missing samples: {missing[:5]} (and possibly more)")
    df = df.loc[:, samples]

    #numeric cleanup
    df = df.apply(pd.to_numeric, errors = "coerce")

    #DEBUG: how many NaNs before clean?
    n_nans = int(np.isnan(df.values).sum())
    if n_nans:
        print(f"[warn] {vf}: found {n_nans} NaNs -> setting to 0")

    #replace NaNs and any negatives with 0 (counts must be nonnegative)
    df = df.fillna(0.0)
    if (df.values < 0).any():
        print(f"[warn] {vf}: found negative values -> clipping to 0")
        df[df < 0] = 0.0

    #guard against all-zero features/samples (rclr cannot handle them)
    nonzero_feat = (df.sum(axis = 1) > 0)
    nonzero_samp = (df.sum(axis = 0) > 0)
    if (~nonzero_feat).any():
        print(f"[warn] {vf}: dropping {int((~nonzero_feat).sum())} all-zero features")
    if (~nonzero_samp).any():
        bad = list(np.array(samples)[~nonzero_samp.values])
        raise ValueError(f"{vf}: these samples are all-zero after cleaning: {bad}")

    df = df.loc[nonzero_feat, nonzero_samp]

    #rCLR: (features x samples) -> transpose -> rclr -> transpose back
    X = df.values.T  #samples x features
    #final sanity check: no NaNs remain
    assert np.isfinite(X).all(), f"{vf}: still has non-finite values before rclr"
    X_rclr = matrix_rclr(X)  #samples x features
    M = pd.DataFrame(X_rclr.T, index = df.index, columns = df.columns)  #features x samples

    #match the R-side post-rclr filtering (finite + sd>0 across samples)
    ok_fin = np.isfinite(M.values).all(axis = 1)
    sds = M.std(axis = 1).to_numpy()
    keep = ok_fin & np.isfinite(sds) & (sds > 0)
    if (~keep).any():
        print(f"[info] {vf}: dropping {int((~keep).sum())} rows after (finite & sd>0) filter")
    M = M.loc[keep]

    M.to_csv(interop/f"rclr_PY_view_{i}.csv")
    print(f"[ok] rclr_PY_view_{i}.csv  shape = {M.shape}")
    py_rclr.append(M)

#compare to R 
mx_abs = []
for i in range(1, len(view_files)+1):
    R = pd.read_csv(interop/f"rclr_R/view_{i}_rclr_R.csv", index_col = 0)
    P = pd.read_csv(interop/f"rclr_PY_view_{i}.csv", index_col = 0)

    #align by intersecting features and the final (possibly reduced) sample set
    common_rows = R.index.intersection(P.index)
    common_cols = R.columns.intersection(P.columns)
    R2 = R.loc[common_rows, common_cols].sort_index(axis = 0).sort_index(axis = 1)
    P2 = P.loc[common_rows, common_cols].sort_index(axis = 0).sort_index(axis = 1)

    diff = (R2 - P2).to_numpy()
    mx = float(np.nanmax(np.abs(diff))) if diff.size else np.nan
    print(f"view_{i}: rCLR max|diff| = {mx:.3e}  (rows = {R2.shape[0]}, cols = {R2.shape[1]})")
    mx_abs.append(mx)


