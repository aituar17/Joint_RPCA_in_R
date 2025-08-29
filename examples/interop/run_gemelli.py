from pathlib import Path
import json
import os
import numpy as np
import pandas as pd
from gemelli.factorization import joint_rpca

#read settings and data
SCRIPT_DIR = Path(__file__).parent           
p = SCRIPT_DIR / "interop"

with open(p / "settings.json", encoding = "utf-8") as fh:
    settings = json.load(fh)
n_components = int(settings["n_components"])
max_iter     = int(settings["max_iter"])
seed         = int(settings["seed"])

#read views (features x samples) and ensure identical sample order
samples = pd.read_csv(p / "samples.csv", encoding = "utf-8")["sample"].tolist()

view_files = sorted([f for f in os.listdir(p) if f.endswith("_rclr.csv")])

if not view_files:
    raise FileNotFoundError(f"No view CSVs found in {p}. Expected files like 'view_1_rclr.csv'.")
    
views = []
for vf in view_files:
    V = pd.read_csv(p / vf, index_col = 0, encoding = "utf-8")  
    #ensure all requested samples exist
    missing = [s for s in samples if s not in V.columns]
    if missing:
        raise ValueError(f"{vf} is missing {len(missing)} samples found in samples.csv. "
                         f"First few missing: {missing[:5]}")
    V = V.loc[:, samples]  #reorder to match 'samples'
    views.append(V)

#Gemelli expects samples x features per view, so transpose each
X_list = [V.T.values for V in views]  #list of (n_samples x n_features_k)

#fit
np.random.seed(seed)
res = joint_rpca(
    tables = X_list,
    n_components = n_components,
    max_iter = max_iter,
    verbose = False
)

S = pd.DataFrame(
    res["sample_loading"],
    index = samples,
    columns = [f"comp{i+1}" for i in range(n_components)]
)

#save outputs for comparison
S.to_csv(p / "gemelli_samplescores.csv", encoding = "utf-8")

for k, V in enumerate(views, start = 1):
    Fk = pd.DataFrame(
        res["feature_loading"][k-1],
        index = V.index,
        columns = [f"comp{i+1}" for i in range(n_components)]
    )
    Fk.to_csv(p / f"gemelli_loadings_view{k}.csv", encoding = "utf-8")

print(f"[OK] Wrote {p/'gemelli_samplescores.csv'} and per-view loadings.")
