from pathlib import Path
import numpy as np, pandas as pd
from scipy.linalg import orthogonal_procrustes
from scipy.stats import pearsonr

p = Path(__file__).parent / "interop"
samples = pd.read_csv(p/"samples.csv")["sample"].astype(str).tolist()

A = pd.read_csv(p/"gemelli_samplescores_seed42.csv", index_col = 0)              #seed 42
B = pd.read_csv(p/"gemelli_samplescores_seed777.csv", index_col=  0)             #seed 777

common = [s for s in samples if s in A.index and s in B.index]
A = A.loc[common]; B = B.loc[common]

# match dimensionality
k = min(A.shape[1], B.shape[1])
A = A.iloc[:, :k].to_numpy(); B = B.iloc[:, :k].to_numpy()

# center columns
A = A - A.mean(axis=0, keepdims=True)
B = B - B.mean(axis=0, keepdims=True)

# Procrustes align: find Q so B @ Q ≈ A
Q, _ = orthogonal_procrustes(B, A)
B_aligned = B @ Q

# per-component Pearson r and overall R, R^2
per_comp = [pearsonr(A[:, j], B_aligned[:, j])[0] for j in range(k)]
num = (A * B_aligned).sum()
den = np.sqrt((A*A).sum()) * np.sqrt((B_aligned*B_aligned).sum())
R = float(num / den); R2 = R**2

print(f"Python vs Python (seeds 42 vs 777): overall R={R:.4f}  R^2={R2:.4f}")
for j, r in enumerate(per_comp, 1):
    print(f"  Comp{j}: r={r:.4f}")
