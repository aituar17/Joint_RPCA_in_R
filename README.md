## Reproducible Example: Joint RPCA on Synthetic Data

A complete, self-contained example of how to run Joint RPCA with synthetic data is provided in:

📄 [`examples/joint_rpca_example.qmd`](examples/joint_rpca_example.qmd)

🖥️ View the rendered HTML output:
[`joint_rpca_example.html`](examples/joint_rpca_example.html)

📋 This example includes:
- Creation of three synthetic count tables with shared samples

- Injection of missing values to simulate sparsity

- Preprocessing with robust centered log-ratio (rclr) transformation

- Definition of a metadata-based train/test split

- Joint RPCA dimensionality reduction via OptSpace

- Visualization of the ordination results using ggplot2

- Ranking of the top features contributing to each principal component

- Computation of a covariance matrix of feature loadings

To reproduce it locally, run:

```r
quarto::quarto_render("examples/joint_rpca_example.qmd")
