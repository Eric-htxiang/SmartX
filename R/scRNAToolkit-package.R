#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import Seurat
#' @import ggplot2
#' @importFrom dplyr %>% arrange mutate group_by summarise ungroup filter select slice_head
#' @importFrom ggrepel geom_text_repel
#' @importFrom patchwork plot_layout
#' @importFrom grid grid.rect grid.text gpar unit
#' @importFrom reshape2 melt
#' @importFrom stats cor splinefun
#' @importFrom methods as
#' @importFrom utils write.table
## usethis namespace: end
NULL

#' scRNAToolkit: A Comprehensive Toolkit for Single-Cell RNA-Seq Analysis
#'
#' @description
#' The scRNAToolkit package provides a suite of functions for streamlined
#' single-cell RNA-seq data analysis and visualization. It builds upon the
#' Seurat framework and integrates with CellChat for cell-cell communication
#' analysis.
#'
#' @section Main Features:
#' \describe{
#'   \item{Data Processing}{Functions for creating count dataframes, 
#'   calculating statistics, and organizing single-cell data}
#'   \item{Visualization}{Comprehensive plotting functions including UMAP plots,
#'   dot plots, violin plots, stacked bar charts, and heatmaps}
#'   \item{Differential Expression}{Tools for identifying differentially 
#'   expressed genes and performing Gene Ontology enrichment analysis}
#'   \item{Correlation Analysis}{Functions for computing expression correlations
#'   and analyzing overlaps between datasets}
#'   \item{Cell Communication}{Wrapper functions for CellChat analysis}
#' }
#'
#' @section Getting Started:
#' Load the package with:
#' \code{library(scRNAToolkit)}
#'
#' For a quick introduction, see:
#' \code{vignette("introduction", package = "scRNAToolkit")}
#'
#' @section Main Functions:
#' \itemize{
#'   \item \code{\link{drawUMAP}}: Create UMAP visualizations
#'   \item \code{\link{drawExprDotPlot}}: Generate expression dot plots
#'   \item \code{\link{run_DEG}}: Perform differential expression analysis
#'   \item \code{\link{run_GOanalysis}}: Run GO enrichment analysis
#'   \item \code{\link{exprCorrelation}}: Calculate expression correlations
#'   \item \code{\link{run_cellChat}}: Perform cell-cell communication analysis
#' }
#'
#' @section Dependencies:
#' This package requires Seurat (>= 4.0.0) and several other packages.
#' Optional dependencies include CellChat, ComplexHeatmap, and clusterProfiler
#' for extended functionality.
#'
#' @docType package
#' @name scRNAToolkit-package
NULL
