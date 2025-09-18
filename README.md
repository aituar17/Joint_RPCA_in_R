# 📊 Reproducible Example: Joint RPCA on HintikkaXOData

A real-world demonstration of Joint RPCA applied to microbiome data from the **mia** package is provided in:
📄 [`examples/joint_rpca_example.qmd`](https://github.com/aituar17/Joint_RPCA_in_R/blob/main/examples/joint_rpca_example.qmd)
🖥️ Rendered HTML output: [`joint_rpca_example.html`](https://github.com/aituar17/Joint_RPCA_in_R/blob/main/examples/joint_rpca_example.html)

## What this example does
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
