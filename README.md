# 📐 Joint RPCA in R
This repository provides an R implementation of **Joint-RPCA** (Robust Principal Component Analysis for multi-omics data), with reproducible workflows and interoperability with the Python **Gemelli** package.

## 🚀 Quickstart for Supervisors / Authors
To directly reproduce the **R ↔ Python (Gemelli)** comparison:
```bash
#1) clone this repository
git clone https://github.com/aituar17/Joint_RPCA_in_R.git
cd Joint_RPCA_in_R

#2) run the R workflow (produces interop/ folder)
Rscript -e 'quarto::quarto_render("examples/joint_rpca_example.qmd")'

#3) run the Python workflow (requires Gemelli installed)
python examples/run_gemelli.py

#4) view comparison results
cat examples/interop/compare_r_vs_py.txt
```

Optional:
```bash
#compare two Python runs with different seeds
#(after changing the seed in settings.json)
python examples/run_gemelli_with_different_seed.py
python examples/compare_two_python_runs.py
```

👉 Outputs (sample scores, rCLR tables, comparison reports) are written to `examples/interop/`.

## 📊 Reproducible Example: Joint RPCA on HintikkaXOData

A real-world demonstration of Joint RPCA applied to microbiome data from the **mia** package is provided in:
📄 [`examples/joint_rpca_example.qmd`](https://github.com/aituar17/Joint_RPCA_in_R/blob/main/examples/joint_rpca_example.qmd)
🖥️ Rendered HTML output: [`joint_rpca_example.html`](https://github.com/aituar17/Joint_RPCA_in_R/blob/main/examples/joint_rpca_example.html)

### What this example does
- Loads the **HintikkaXOData** dataset (via MultiAssayExperiment)
- Applies **robust centered log-ratio (rCLR)** transformation
- Splits samples into manual train/test sets
- Runs `jointRPCAuniversal()` from the **R implementation**
- Visualizes ordination of sample embeddings
- Extracts and ranks top features driving PCs
- Benchmarks Joint RPCA features against other methods with a Random Forest classifier
- **Exports results for interoperability with Python (Gemelli)**

## 🔄 R ↔ Python Interoperability
The QMD example produces an `interop/` folder containing:
- `R_samplescores.csv`, `samples.csv`, `settings.json`
- Count tables: `view_1_counts.csv`, `view_2_counts.csv`, `view_3_counts.csv`
- rCLR-transformed tables: `rclr_R/view_1_rclr_R.csv`, etc.

These files can be used directly with the Python **Gemelli** implementation.

### Running the Python side
1. Install Gemelli (conda / pip).
2. Run:
```bash
python examples/run_gemelli.py
```

This will run Joint RPCA with the settings from `settings.json`, save outputs back to `interop/`, and compare the **Python sample scores** to the **R sample scores**.
Results are written to `interop/compare_r_vs_py.txt`. You can also view gemelli_loadings_view1_seed42.csv, gemelli_loadings_view2_seed42.csv, gemelli_loadings_view3_seed42.csv, and gemelli_samplescores_seed42.csv.

## 🧪 Comparing Multiple Runs
- **R vs Python (Gemelli):**
    Run `examples/run_gemelli.py`.
    → Alignment report saved to `interop/compare_r_vs_py.txt`.
- **Python vs Python (different seeds):**
    1. Edit `settings.json` to change the seed (e.g., to 777).
    2. Run:
    ```bash
    examples/run_gemelli_with_different_seed.py
    ```

    This will produce gemelli_loadings_view1_seed777.csv, gemelli_loadings_view2_seed777.csv, gemelli_loadings_view3_seed777.csv, and gemelli_samplescores_seed777.csv.

    3. Run:
    ```bash
    examples/compare_two_python_runs.py
    ```

## ▶️ How to Reproduce Everything
To run the full R → Python comparison workflow locally:
```r
#from R
quarto::quarto_render("examples/joint_rpca_example.qmd")
```

```bash
#from Python
python examples/run_gemelli.py
python examples/run_gemelli_with_different_seed.py  #optional (after changing the seed in settings.json)
python examples/compare_two_python_runs.py          #optional
```

## 📑 Expected Output Example
When you run the R ↔ Python comparison (`compare_r_vs_py.txt`), you should see output in this format:
```yaml
R vs Python (Gemelli) sample scores — Orthogonal Procrustes alignment
Common samples: 40 | Matched dims: 3
Overall R: 0.9582  R^2: 0.9182
 Comp1: Pearson r (aligned) = 0.9416
 Comp2: Pearson r (aligned) = 0.9725
 Comp3: Pearson r (aligned) = 0.9673
```
(Values may vary slightly depending on seed and max_iterations.)

## ⚙️ Parameters Used in Comparison
For the example results above:
- **max_iterations:** 2000 (R and Python)
- **n_components:** 3
- **rclr_transform_tables:** True (default)
- **Seed:** 42
- **min_sample_count:** 1
- **min_feature_count** 1
- **min_feature_frequency:** 0.0
- **Dataset:** HintikkaXOData (from the mia package)

## 🧬 Single-Omic Replication: HMP2 IBD (16S rRNA)

This reproduces the iHMP IBD **16S-only** analysis in R with Joint-RPCA and basic benchmarking.

### What it does
- Loads **HMP2 IBD 16S** counts, taxonomy, and sample metadata from `HMP2Data`
- Builds a `SummarizedExperiment` and `MultiAssayExperiment` (single view)
- Runs `jointRPCAuniversal()` with rCLR (default)  
- Produces:
  - PC1–PC2 scatter (IBD vs non-IBD), Wilcoxon tests (PC1/PC2)
  - **PERMANOVA** on V1–V3
  - **PC1 log-ratio** test (top vs bottom loadings)
  - **Weighted Random Forest** (OOB error + **OOB AUROC**)
  - Optional **NMF** baseline and AUROC bar chart

### Where it lives
- QMD: `examples/ihmp_ibd_replication.qmd`
- Rendered HTML: `examples/ihmp_ibd_replication.html`

### Quickstart

```bash
# 1) clone this repository
git clone https://github.com/aituar17/Joint_RPCA_in_R.git
cd Joint_RPCA_in_R

# 2) render the analysis
Rscript -e "quarto::quarto_render('examples/ihmp_ibd_replication.qmd')"

# 3) open the report
open examples/ihmp_ibd_replication.html   # macOS
# xdg-open examples/ihmp_ibd_replication.html  # Linux

```

## 🧩 Multi-Omic Validation on IBDMDB Data

📄 [`examples/ibdmdb_2omic_jointrpca.qmd`](https://github.com/aituar17/Joint_RPCA_in_R/blob/main/examples/ibdmdb_2omic_jointrpca.qmd)
🖥️ [`examples/ibdmdb_2omic_jointrpca.html`](https://github.com/aituar17/Joint_RPCA_in_R/blob/main/examples/ibdmdb_2omic_jointrpca.html)

### Description:
- Integrates **metagenomics (MGX)** and **metatranscriptomics (MTX)** data from IBDMDB
- Applies strict shared-sample matching
- Produces stable ordination and top-loading taxa consistent with Gemelli results
- Confirms numerical and biological equivalence between the R and Python implementations

📄 [`examples/ibdmdb_benchmarking.qmd`](https://github.com/aituar17/Joint_RPCA_in_R/blob/main/examples/ibdmdb_benchmarking.qmd)
🖥️ [`examples/ibdmdb_benchmarking.html`](https://github.com/aituar17/Joint_RPCA_in_R/blob/main/examples/ibdmdb_benchmarking.html)

### Description:
- Benchmarks Joint-RPCA against PCA and NMF on the same IBD data
- Reports Wilcoxon, PERMANOVA, and AUROC metrics
- Confirms that Joint-RPCA achieves comparable or superior variance capture and separation performance

📄 [`examples/ibdmdb_3omic_jointrpca.qmd`](https://github.com/aituar17/Joint_RPCA_in_R/blob/main/examples/ibdmdb_3omic_jointrpca.qmd)
🖥️ [`examples/ibdmdb_3omic_jointrpca.html`](https://github.com/aituar17/Joint_RPCA_in_R/blob/main/examples/ibdmdb_3omic_jointrpca.html)

### Description:
- Extends to **3-omic integration (16S + MGX + MTX)**
- Uses strict shared-sample intersection (16 samples total)
- Performs unsupervised ordination, variance explained, and component correlation plots
- Labels unavailable for non-IBD samples, so analysis focuses on **unsupervised behavior**
- Confirms algorithmic stability and expected cross-omic structure

### 📂 Data Availability
The data needed to run examples/ibdmdb_2omic_jointrpca.qmd, examples/ibdmdb_benchmarking.qmd, 
and examples/ibdmdb_3omic_jointrpca.qmd is available in the examples/data_ibdmdb_raw folder. 
It includes taxonomic_profiles_16s.tsv (taxonomic profiles for 16S from HMP2), 
taxonomic_profiles_mgx.tsv (taxonomic profiles for MGX from HMP2), hmp2_metadata_2018-08-20.csv (HMP2 Metadata),
taxonomic_profiles_mgx_new.tsv (taxonomic profiles for MGX from HMP2_Pilot), 
taxonomic_profiles_mtx_new.tsv (taxonomic profiles for MTX from HMP2_Pilot), 
and taxonomic_profiles_16s_new.tsv (taxonomic profiles for 16S from HMP2_Pilot).

Due to size limits, the file `data_ibdmdb_raw/ecs_relab.tsv`(taxonomic profiles for MTX from HMP2)
is not included in the repository.
You can download it manually from the [iHMP IBDMDB data portal](https://ibdmdb.org/downloads/html/products_MTX_2017-12-14.html)
and place it in the same folder before running:
```r
quarto::quarto_render("examples/ibdmdb_2omic_jointrpca.qmd")
quarto::quarto_render("examples/ibdmdb_benchmarking.qmd")
```

### 🧾 Summary
These IBDMDB analyses demonstrate that:
- The **R Joint-RPCA** implementation reproduces the **Gemelli (Python)** behavior on real multi-omic IBD data.
- Results remain **numerically stable**, **biologically plausible**, and **methodologically equivalent** across omic combinations.
- 2-omic integration provides robust validation, while the 3-omic case confirms algorithmic consistency under limited overlap.
