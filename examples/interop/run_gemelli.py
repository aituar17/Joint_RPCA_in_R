from pathlib import Path
import json, os, numpy as np, pandas as pd
from gemelli.factorization import joint_rpca

#read settings and data
SCRIPT_DIR = Path(__file__).parent
ROOT = SCRIPT_DIR.parent
p = ROOT / "interop"

settings = json.load(open(p/"settings.json", encoding = "utf-8"))
n_components = int(settings["n_components"]); max_iter = int(settings["max_iter"]); seed = int(settings["seed"])

#read views (features x samples) and ensure identical sample order
samples = pd.read_csv(p/"samples.tsv", sep = "\t", encoding = "utf-8")["sample"].tolist()

view_files = sorted([f for f in os.listdir(p) if f.endswith("_rclr.tsv")])
views = []
for vf in view_files:
    V = pd.read_csv(p/vf, sep = "\t", index_col = 0, encoding = "utf-8")
    V = V.loc[:, samples]
    views.append(V)

#Gemelli expects samples x features per view, so transpose each
X_list = [V.T.values for V in views]  

#fit
np.random.seed(seed)
res = joint_rpca(tables = X_list, n_components = n_components, max_iter = max_iter, verbose = False)

S = pd.DataFrame(res["sample_loading"], index = samples, columns = [f"comp{i + 1}" for i in range(n_components)])

#save outputs for comparison
S.to_csv(p/"gemelli_samplescores.tsv", sep = "\t", encoding = "utf-8")


