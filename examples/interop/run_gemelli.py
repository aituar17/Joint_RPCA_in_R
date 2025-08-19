from pathlib import Path
import json, os, numpy as np, pandas as pd

#read settings and data
SCRIPT_DIR = Path(__file__).parent
p = SCRIPT_DIR / "interop"   
settings = json.load(open(p/"settings.json"))
n_components = int(settings["n_components"])
max_iter     = int(settings["max_iter"])
seed         = int(settings["seed"])

#read views (features x samples) and ensure identical sample order
samples = pd.read_csv(p/"samples.csv")["sample"].tolist()

view_files = sorted([f for f in os.listdir(p) if f.endswith("_rclr.csv")])
views = []
for vf in view_files:
    V = pd.read_csv(p/vf, index_col = 0)
    V = V.loc[:, samples] 
    views.append(V)

#Gemelli joint-RPCA
from gemelli.factorization import joint_rpca

#Gemelli expects samples x features per view, so transpose each
X_list = [V.T.values for V in views]  

#fit
np.random.seed(seed)
res = joint_rpca(
    tables = X_list,
    n_components = n_components,
    max_iter = max_iter,
    verbose = False
)

S = pd.DataFrame(res["sample_loading"], index = samples,
                 columns = [f"comp{i + 1}" for i in range(n_components)])

#save outputs for comparison
S.to_csv(p/"gemelli_samplescores.csv")

for k, V in enumerate(views):
    Fk = pd.DataFrame(res["feature_loading"][k], index = V.index,
                      columns = [f"comp{i + 1}" for i in range(n_components)])
    Fk.to_csv(p/f"gemelli_loadings_view{k + 1}.csv")
