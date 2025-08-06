## Reproducible Example: Joint RPCA on HintikkaXOData

A real-world demonstration of Joint RPCA applied to microbiome data from the mia package is provided in:

📄 [`examples/joint_rpca_example.qmd`](examples/joint_rpca_example.qmd)

🖥️ View the rendered HTML output:
[`joint_rpca_example.html`](examples/joint_rpca_example.html)

📋 This example includes:
- Direct usage of the HintikkaXOData dataset via a MultiAssayExperiment object

- Preprocessing with robust centered log-ratio (rclr) transformation

- Sample grouping into manual train/test sets

- Ordination of sample embeddings using the jointRPCAuniversal() wrapper

- Visualization of sample clusters colored by set assignment

- Ranking and visualization of top features driving variation in principal components

- Covariance matrix analysis of feature loadings

This example highlights seamless integration between Joint RPCA and Bioconductor data structures, showcasing the method's ability to handle real microbiome datasets with compositional preprocessing and dimensionality reduction.

To reproduce it locally, run:

```r
quarto::quarto_render("examples/joint_rpca_example.qmd")
