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
