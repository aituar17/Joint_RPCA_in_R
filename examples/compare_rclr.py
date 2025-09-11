from pathlib import Path
import pandas as pd, numpy as np, os
from gemelli.preprocessing import matrix_rclr

interop = Path("interop")
view_files = sorted([f for f in os.listdir(interop) if f.endswith("_counts.csv")])
samples = pd.read_csv(interop/"samples.csv")["sample"].astype(str).tolist()




for f in interop.glob("rclr_PY_view_*.csv"):
    try: f.unlink()
    except: pass

for i, vf in enumerate(view_files, start = 1):
    df = pd.read_csv(interop/vf, index_col = 0)       
    df = df.loc[:, samples].apply(pd.to_numeric, errors = "coerce")

    #report + fix NaNs
    n_nans = int(np.isnan(df.values).sum())
    if n_nans:
        print(f"[warn] {vf}: found {n_nans} NaNs -> setting to 0 before rclr")
        df = df.fillna(0)

    if USE_TRANSPOSE:
        X = df.values.T                            
        print(f"[{vf}] USING TRANSPOSE: df {df.shape} -> X {X.shape} (samples x features)")
        X_rclr = matrix_rclr(X)                      
        M = pd.DataFrame(X_rclr.T, index = df.index, columns = df.columns)  
    else:
        X = df.values                                
        print(f"[{vf}] NO TRANSPOSE: df {df.shape} -> X {X.shape} (features x samples)")
        X_rclr = matrix_rclr(X)                     
        M = pd.DataFrame(X_rclr, index = df.index, columns = df.columns)

    #match R’s downstream filtering (finite & non-zero sd per feature)
    ok_fin = np.isfinite(M.values).all(axis = 1)
    sds = M.std(axis = 1).values
    keep = ok_fin & np.isfinite(sds) & (sds > 0)
    dropped = (~keep).sum()
    if dropped:
        print(f"[{vf}] dropping {dropped} rows after (finite & sd>0) filter")
    M = M.loc[keep]

    out = interop/f"rclr_PY_view_{i}.csv"
    M.to_csv(out)
    print(f"[ok] wrote {out}  shape = {M.shape}")
