
<img width="3693" height="1380" alt="image" src="https://github.com/user-attachments/assets/1e0baa91-8bd1-42f6-81ee-8e333fae17b8" />





**Spatial Multiomic Analysis and Research Toolkit for eXploration**

SMART-X is a modular R framework for integrative analysis of single-cell and spatial transcriptomics data, it transforms complex analytical workflows into simple, reproducible one-line commands while maintaining publication-quality output.

Key Features

📊 Single-Cell Analysis
One-Line Visualization: Generate publication-quality UMAP, tSNE, dot plots, and violin plots
Differential Expression: Multi-method DEG analysis with flexible group comparisons
Functional Annotation: Integrated GO enrichment with volcano plots
Cell Communication: CellChat integration for ligand-receptor analysis

🔬 Spatial Transcriptomics (v0.2.0+)
Spatial pattern visualization
Spatially variable gene identification
Spatial domain detection
Spatial cell communication networks

📈 Comparative Analysis
Cross-dataset expression correlation
DEG overlap analysis (Jaccard similarity)
GO term comparison
Publication-ready correlation heatmaps

🎨 Visualization Excellence
Automatic aesthetic optimization
Consistent styling across all plots
Customizable color schemes
PDF output with adjustable dimensions


**INSTALL**

if (!require("devtools")) install.packages("devtools")

devtools::install_github("Eric-htxiang/SmartX")
