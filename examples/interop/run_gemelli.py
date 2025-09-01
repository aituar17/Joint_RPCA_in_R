from pathlib import Path
import json
import os
import numpy as np
import pandas as pd

# --- Gemelli import: we want the BIOM-based API from gemelli.rpca ---
try:
    from gemelli.rpca import joint_rpca  # expects BIOM Tables
except Exception as e:
    raise ImportError(
        "Could not import gemelli.rpca.joint_rpca. "
        "Ensure you're using gemelli==0.0.12 inside the gemelli39 env."
    ) from e

# BIOM table
from biom import Table

# Optional: for Procrustes comparison with R loadings
from scipy.linalg import orthogonal_procrustes


def read_settings_and_samples(interop_dir: Path):
    with open(interop_dir / "settings.json", encoding="utf-8") as fh:
        settings = json.load(fh)

    # accept either key
    n_components  = int(settings.get("n_components"))
    max_iterations = int(
        settings.get("max_iterations", settings.get("max_iter", 500))
    )
    seed = int(settings.get("seed", 42))

    samples = (
        pd.read_csv(interop_dir / "samples.csv", encoding="utf-8")["sample"]
        .astype(str)
        .tolist()
    )
    if len(samples) == 0:
        raise ValueError("samples.csv has no sample IDs.")
    return settings, n_components, max_iterations, seed, samples


def load_counts_views_as_biom(interop_dir: Path, samples):
    # We expect files like: view_1_counts.csv, view_2_counts.csv, ...
    view_files = sorted([f for f in os.listdir(interop_dir) if f.endswith("_counts.csv")])
    if not view_files:
        raise FileNotFoundError(
            f"No '*_counts.csv' found in {interop_dir}. "
            "Export them from R as features×samples CSV."
        )

    biom_tables = []
    for vf in view_files:
        df = pd.read_csv(interop_dir / vf, index_col=0)  # features x samples
        # enforce given sample order and presence
        df.columns = df.columns.astype(str)
        missing = [s for s in samples if s not in df.columns]
        if missing:
            raise ValueError(f"{vf} is missing {len(missing)} samples: first few {missing[:5]}")
        df = df.loc[:, samples]

        # Coerce to numeric (counts), drop non-finite or negative values if any
        df = df.apply(pd.to_numeric, errors="coerce")
        if not np.isfinite(df.values).all():
            # drop rows with non-finite values
            good_rows = np.isfinite(df.values).all(axis=1)
            df = df.loc[good_rows]
        # optional: drop all-zero feature rows
        nonzero_rows = (df.values.sum(axis=1) > 0)
        if nonzero_rows.sum() == 0:
            raise ValueError(f"All-zero features in {vf} after filtering; nothing to model.")
        df = df.loc[nonzero_rows]

        # Build BIOM table
        T = Table(
            df.values.astype(float),
            observation_ids=df.index.astype(str).tolist(),
            sample_ids=df.columns.astype(str).tolist()
        )
        biom_tables.append(T)
        print(f"[read] {vf}: {df.shape[0]} features × {df.shape[1]} samples")

    return biom_tables, view_files


def run_joint_rpca(biom_tables, n_components, max_iterations, seed):
    np.random.seed(seed)
    return joint_rpca(
        tables=biom_tables,
        n_components=n_components,
        max_iterations=max_iterations,
        # disable internal filtering:
        min_sample_count=1,
        min_feature_count=1,
        min_feature_frequency=0.0,
        # leave rclr on (default) since you exported raw counts
        # rclr_transform_tables=True
    )


def save_outputs(interop_dir: Path, res, samples, view_files, biom_tables):
    # unpack ordination / loadings
    if len(res) == 4:
        ord_res, dist_res, feat_loads_raw, sample_loadings = res
        S = pd.DataFrame(sample_loadings, index=[str(s) for s in samples],
                         columns=[f"comp{i+1}" for i in range(sample_loadings.shape[1])])
    else:
        ord_res, dist_res, feat_loads_raw = res
        S = ord_res.samples.copy()
        S = S.loc[[str(s) for s in samples]]
        S.columns = [f"comp{i+1}" for i in range(S.shape[1])]

    S.to_csv(interop_dir / "gemelli_samplescores.csv", encoding="utf-8")
    print(f"[write] gemelli_samplescores.csv  shape={S.shape}")

    # normalize feature loadings to per-view DataFrames using feature IDs
    per_view_df = []

    if isinstance(feat_loads_raw, (list, tuple)):
        # already per-view
        for i, arr in enumerate(feat_loads_raw):
            ids = list(biom_tables[i].ids(axis="observation"))
            Fk = pd.DataFrame(arr, index=ids,
                              columns=[f"comp{i+1}" for i in range(arr.shape[1])])
            per_view_df.append(Fk)
    else:
        # single stacked DataFrame or ndarray
        if isinstance(feat_loads_raw, pd.DataFrame) and feat_loads_raw.index.is_unique:
            FL = feat_loads_raw
        else:
            # fall back to offsets only if no IDs are present
            import numpy as np
            FL = (pd.DataFrame(feat_loads_raw)
                  .set_index(pd.Index(sum([list(b.ids(axis="observation")) for b in biom_tables], []))))
        # split by IDs per view
        for i in range(len(biom_tables)):
            ids = list(biom_tables[i].ids(axis="observation"))
            Fk = FL.loc[ids]
            # ensure nice column names
            Fk.columns = [f"comp{i+1}" for i in range(Fk.shape[1])]
            per_view_df.append(Fk)

    # write per-view files
    for i, Fk in enumerate(per_view_df, start=1):
        out = interop_dir / f"gemelli_loadings_view{i}.csv"
        Fk.to_csv(out, encoding="utf-8")
        print(f"[write] {out.name}  shape={Fk.shape}")

def optional_compare_with_R(interop_dir: Path, samples):
    """Compare Gemelli sample scores with R (your jointRPCAuniversal) using orthogonal Procrustes.
    Writes a tiny report to interop/compare_r_vs_py.txt
    """
    r_path = interop_dir / "R_samplescores.csv"
    py_path = interop_dir / "gemelli_samplescores.csv"
    if not (r_path.exists() and py_path.exists()):
        print("[skip] R_samplescores.csv or gemelli_samplescores.csv not found — skipping comparison.")
        return

    R = pd.read_csv(r_path, index_col=0)
    P = pd.read_csv(py_path, index_col=0)

    # Align rows and shared columns
    common = [s for s in samples if s in R.index and s in P.index]
    if len(common) < 3:
        print("[skip] Less than 3 common samples for Procrustes — skipping comparison.")
        return

    R = R.loc[common]
    P = P.loc[common]

    # Match dimensionality: use the min(Kr, Kp)
    k = min(R.shape[1], P.shape[1])
    Rk = R.iloc[:, :k].to_numpy()
    Pk = P.iloc[:, :k].to_numpy()

    # Center columns (zero mean) to be fair
    Rk = Rk - Rk.mean(axis=0, keepdims=True)
    Pk = Pk - Pk.mean(axis=0, keepdims=True)

    # Orthogonal Procrustes (find Q so that Pk @ Q ≈ Rk)
    Q, _ = orthogonal_procrustes(Pk, Rk)
    Pk_aligned = Pk @ Q

    # Report per-component Pearson R and overall R^2
    from scipy.stats import pearsonr
    comps = []
    for j in range(k):
        r, _ = pearsonr(Rk[:, j], Pk_aligned[:, j])
        comps.append((j+1, float(r)))

    # Overall fit: squared Frobenius correlation
    num = (Rk * Pk_aligned).sum()
    denom = np.sqrt((Rk * Rk).sum()) * np.sqrt((Pk_aligned * Pk_aligned).sum())
    overall_r = float(num / denom)
    overall_r2 = overall_r ** 2

    out = interop_dir / "compare_r_vs_py.txt"
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("R vs Python (Gemelli) sample scores — Orthogonal Procrustes alignment\n")
        fh.write(f"Common samples: {len(common)} | Matched dims: {k}\n")
        fh.write(f"Overall R: {overall_r:.4f}  R^2: {overall_r2:.4f}\n")
        for j, r in comps:
            fh.write(f" Comp{j}: Pearson r (aligned) = {r:.4f}\n")
    print(f"[write] compare_r_vs_py.txt  overall R^2={overall_r2:.4f}")


def main():
    SCRIPT_DIR = Path(__file__).parent
    interop_dir = SCRIPT_DIR / "interop"

    settings, n_components, max_iterations, seed, samples = read_settings_and_samples(interop_dir)
    biom_tables, view_files = load_counts_views_as_biom(interop_dir, samples)
    res = run_joint_rpca(biom_tables, n_components, max_iterations, seed)
    save_outputs(interop_dir, res, samples, view_files, biom_tables)
    optional_compare_with_R(interop_dir, samples)

if __name__ == "__main__":
    main()
