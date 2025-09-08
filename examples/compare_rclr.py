from pathlib import Path
import pandas as pd, numpy as np, os
from biom import Table
from gemelli.preprocessing import matrix_rclr

interop = Path("interop")

#rebuild the same BIOM tables you fed gemelli
view_files = sorted([f for f in os.listdir(interop) if f.endswith("_counts.csv")])
samples = pd.read_csv(interop/"samples.csv")["sample"].astype(str).tolist()

#python rCLR on the same counts
py_rclr = []
for i, vf in enumerate(view_files, start = 1):
    df = pd.read_csv(interop/vf, index_col = 0)
    df = df.loc[:, samples].apply(pd.to_numeric, errors = "coerce")
    #gemelli's rclr expects counts (nonnegative); apply rclr along axis = 1 over samples
    #matrix_rclr returns samples x features for numpy arrays; use transpose carefully if needed
    #We’ll do: (features x samples) -> transpose -> rclr -> transpose back
    X = df.values.T  #samples x features
    X_rclr = matrix_rclr(X)   #samples x features (float)
    M = pd.DataFrame(X_rclr.T, index = df.index, columns = df.columns)  #features x samples

    #match R's filtering (finite + sd>0 across samples)
    ok_fin = np.isfinite(M.values).all(axis = 1)
    sds = M.std(axis = 1).values
    keep = ok_fin & np.isfinite(sds) & (sds > 0)
    M = M.loc[keep]
    M.to_csv(interop/f"rclr_PY_view_{i}.csv")
    py_rclr.append(M)

#compare to R rCLR
mx_abs = []
for i in range(1, len(view_files) + 1):
    R = pd.read_csv(interop/f"rclr_R/view_{i}_rclr_R.csv", index_col = 0)
    P = pd.read_csv(interop/f"rclr_PY_view_{i}.csv", index_col = 0)
    #align rows and columns
    common_rows = R.index.intersection(P.index)
    R = R.loc[common_rows, samples]
    P = P.loc[common_rows, samples]
    diff = (R - P).to_numpy()
    mx = np.nanmax(np.abs(diff))
    print(f"view_{i}: rCLR max|diff| = {mx:.3e}, shape = {R.shape}")
    mx_abs.append(mx)

print("ALL VIEWS max|diff|:", max(mx_abs))
