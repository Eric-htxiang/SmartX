# =============================================================================
# Correlation and Overlap Analysis Functions
# =============================================================================

#' Calculate expression correlation between two datasets
#'
#' @description
#' Computes Pearson correlation coefficients between gene expression profiles
#' of two datasets (e.g., query and reference cell types).
#'
#' @param expr1 Matrix or dataframe of expression values (genes x samples/cells)
#' @param expr2 Matrix or dataframe of expression values (genes x samples/cells)
#'
#' @return A correlation matrix with dimensions: ncol(expr1) x ncol(expr2)
#'
#' @details
#' Only genes present in both datasets are used for correlation calculation.
#' Expression values should be normalized (e.g., log-normalized counts).
#'
#' @export
#' @examples
#' \dontrun{
#' # Get average expression per cluster
#' expr1 <- AverageExpression(seurat_obj1, group.by = "seurat_clusters")$RNA
#' expr2 <- AverageExpression(seurat_obj2, group.by = "cell_type")$RNA
#' 
#' cor_matrix <- exprCorrelation(expr1, expr2)
#' }
exprCorrelation <- function(expr1, expr2) {
  # Input validation
  if (!is.matrix(expr1) && !is.data.frame(expr1)) {
    stop("expr1 must be a matrix or dataframe")
  }
  if (!is.matrix(expr2) && !is.data.frame(expr2)) {
    stop("expr2 must be a matrix or dataframe")
  }
  
  # Find common genes
  common_genes <- intersect(rownames(expr1), rownames(expr2))
  
  if (length(common_genes) == 0) {
    stop("No common genes found between the two expression matrices")
  }
  
  if (length(common_genes) < 100) {
    warning("Only ", length(common_genes), " common genes found. ",
            "Results may not be reliable.")
  }
  
  # Subset to common genes
  expr1 <- expr1[common_genes, ]
  expr2 <- expr2[common_genes, ]
  
  # Calculate correlation
  cor_matrix <- cor(expr1, expr2, method = "pearson")
  
  message("Correlation calculated using ", length(common_genes), " common genes")
  return(cor_matrix)
}


#' Calculate Jaccard similarity based on DEG overlap
#'
#' @description
#' Computes Jaccard similarity coefficients between groups based on overlap
#' of their top differentially expressed genes.
#'
#' @param markers1 DEG dataframe from first dataset (from run_DEG)
#' @param markers2 DEG dataframe from second dataset (from run_DEG)
#' @param groupBy1 Character string specifying grouping column in markers1
#' @param groupBy2 Character string specifying grouping column in markers2
#' @param pCutoff Numeric, adjusted p-value cutoff (default: 0.05)
#' @param FCcutoff Numeric, log2 fold change cutoff (default: 0.5)
#' @param pct1Cutoff Numeric, minimum percentage in group 1 (default: 0.1)
#' @param pct2Cutoff Numeric, minimum percentage in group 2 (default: 0.01)
#' @param topGeneNum Integer, number of top genes to use per group (default: 100)
#' @param prefix Character string for output filename prefix
#'
#' @return A similarity matrix with Jaccard indices
#'
#' @details
#' Jaccard index = |intersection| / |union|
#' 
#' Also saves a TSV file with overlapping genes for each pair of groups.
#'
#' @export
#' @examples
#' \dontrun{
#' markers1 <- run_DEG(dataset1, groupBy = "seurat_clusters")
#' markers2 <- run_DEG(dataset2, groupBy = "cell_type")
#' 
#' jaccard_matrix <- degOverlapJaccard(
#'   markers1, markers2,
#'   groupBy1 = "seurat_clusters",
#'   groupBy2 = "cell_type",
#'   topGeneNum = 100,
#'   prefix = "Dataset_Comparison"
#' )
#' }
degOverlapJaccard <- function(
    markers1, markers2, groupBy1, groupBy2, pCutoff = 0.05, FCcutoff = 0.5,
    pct1Cutoff = 0.1, pct2Cutoff = 0.01, topGeneNum = 100, prefix) {
  
  # Filter DEGs
  df1 <- markers1[markers1$p_val_adj < pCutoff & markers1$avg_log2FC > FCcutoff & 
                   markers1$pct.1 > pct1Cutoff & markers1$pct.2 > pct2Cutoff, ]
  
  df2 <- markers2[markers2$p_val_adj < pCutoff & markers2$avg_log2FC > FCcutoff & 
                   markers2$pct.1 > pct1Cutoff & markers2$pct.2 > pct2Cutoff, ]
  
  # Select top genes per group
  df1 <- df1 %>% group_by(.data[[groupBy1]]) %>% 
    arrange(desc(avg_log2FC)) %>% slice_head(n = topGeneNum)
  
  df2 <- df2 %>% group_by(.data[[groupBy2]]) %>% 
    arrange(desc(avg_log2FC)) %>% slice_head(n = topGeneNum)
  
  groupIDs1 <- as.character(sort(unique(df1[[groupBy1]])))
  groupIDs2 <- as.character(sort(unique(df2[[groupBy2]])))
  
  # Initialize similarity matrix
  similarity_matrix <- matrix(0, nrow = length(groupIDs1), ncol = length(groupIDs2))
  rownames(similarity_matrix) <- groupIDs1
  colnames(similarity_matrix) <- groupIDs2
  
  # Calculate Jaccard similarity and collect overlaps
  overlapDEG_df <- data.frame()
  for (i in seq(1, length(groupIDs1))) {
    for (j in seq(1, length(groupIDs2))) {
      topDEG1 <- df1[df1[[groupBy1]] == groupIDs1[i], ]$gene
      topDEG2 <- df2[df2[[groupBy2]] == groupIDs2[j], ]$gene
      
      intersectSet <- intersect(topDEG1, topDEG2)
      unionSet <- union(topDEG1, topDEG2)
      similarity_matrix[i, j] <- ifelse(length(unionSet) == 0, 0, 
                                       length(intersectSet) / length(unionSet))
      
      if (length(intersectSet) > 0) {
        items <- c(groupIDs1[i], groupIDs2[j], length(intersectSet), 
                  paste(intersectSet, collapse = ","))
        overlapDEG_df <- rbind(overlapDEG_df, items)
      }
    }
  }
  
  # Save overlapping genes
  names(overlapDEG_df) <- c(groupBy1, groupBy2, "overlapNum", "overlapedGenes")
  outFile <- paste0(prefix, "_", groupBy1, ".", groupBy2, "_overlapedGenes.tsv")
  write.table(overlapDEG_df, outFile, sep = "\t", quote = FALSE, row.names = FALSE)
  message("Overlapping genes saved to: ", outFile)
  
  return(similarity_matrix)
}


#' Calculate Jaccard similarity based on GO term overlap
#'
#' @description
#' Computes Jaccard similarity coefficients between groups based on overlap
#' of enriched GO terms.
#'
#' @param GOresults GO enrichment results dataframe (from run_GOanalysis)
#' @param groupBy Character string specifying grouping column
#' @param pCutoff Numeric, adjusted p-value cutoff for GO terms (default: 0.05)
#' @param topGOnum Integer, maximum number of top GO terms per group (default: 100000)
#'
#' @return A similarity matrix with Jaccard indices
#'
#' @details
#' This function is useful for comparing functional similarity between
#' groups based on their enriched biological processes.
#'
#' @export
#' @examples
#' \dontrun{
#' go_results <- run_GOanalysis(markers, groupBy = "seurat_clusters",
#'                               prefix = "MyProject")
#' 
#' jaccard_matrix <- GOtermsOverlapJaccard(go_results,
#'                                         groupBy = "seurat_clusters",
#'                                         pCutoff = 0.05,
#'                                         topGOnum = 50)
#' }
GOtermsOverlapJaccard <- function(GOresults, groupBy, pCutoff = 0.05, topGOnum = 100000) {
  # Input validation
  if (!groupBy %in% colnames(GOresults)) {
    stop(paste("Column", groupBy, "not found in GO results"))
  }
  
  if (!"Description" %in% colnames(GOresults)) {
    stop("GO results must contain 'Description' column")
  }
  
  if (!"p.adjust" %in% colnames(GOresults)) {
    stop("GO results must contain 'p.adjust' column")
  }
  
  # Filter and select top GO terms
  GOresults <- GOresults %>% 
    filter(p.adjust < pCutoff) %>%
    arrange(.data[[groupBy]], desc(Count)) %>%
    group_by(.data[[groupBy]]) %>%
    slice_head(n = topGOnum)
  
  groupIDs <- as.character(sort(unique(GOresults[[groupBy]])))
  
  if (length(groupIDs) < 2) {
    stop("Need at least 2 groups for comparison")
  }
  
  # Initialize similarity matrix
  similarity_matrix <- matrix(0, nrow = length(groupIDs), ncol = length(groupIDs))
  rownames(similarity_matrix) <- groupIDs
  colnames(similarity_matrix) <- groupIDs
  
  # Calculate Jaccard similarity
  for (i in seq(1, length(groupIDs))) {
    for (j in seq(1, length(groupIDs))) {
      terms_i <- GOresults[GOresults[[groupBy]] == groupIDs[i], ]$Description
      terms_j <- GOresults[GOresults[[groupBy]] == groupIDs[j], ]$Description
      intersectSet <- intersect(terms_i, terms_j)
      unionSet <- union(terms_i, terms_j)
      similarity_matrix[i, j] <- ifelse(length(unionSet) == 0, 0, 
                                       length(intersectSet) / length(unionSet))
    }
  }
  
  return(similarity_matrix)
}


#' Plot correlation heatmap with ComplexHeatmap
#'
#' @description
#' Creates a publication-quality correlation heatmap using ComplexHeatmap package.
#'
#' @param cor_matrix Correlation matrix (from exprCorrelation or similar)
#' @param outPDF Character string, output PDF filename
#' @param column_title Character string, title for columns (default: "query cell type")
#' @param row_title Character string, title for rows (default: "reference cell type")
#' @param legend_title Character string, legend title
#' @param label Logical, whether to show correlation values in cells (default: TRUE)
#' @param labelSize Numeric, label font size (default: 5)
#' @param labValueCutoff Numeric, minimum correlation to show label (default: 0.4)
#' @param highLabCutoff Numeric, threshold for using highLabCol (default: 0.6)
#' @param highLabCol Character string, color for high correlation labels (default: "white")
#' @param row_names_rot Numeric, row name rotation angle (default: 0)
#' @param column_names_rot Numeric, column name rotation angle (default: 45)
#' @param diagonal_color Character string, color for diagonal cells (default: "gray80")
#' @param plotHeight Numeric, plot height in inches (default: 10)
#' @param plotWidth Numeric, plot width in inches (default: 12)
#' @param plotCommonCellType Logical, legacy parameter (default: TRUE)
#'
#' @return NULL. Saves heatmap to PDF file
#'
#' @export
#' @examples
#' \dontrun{
#' cor_mat <- exprCorrelation(expr1, expr2)
#' plot_correlation_heatmap_complex(
#'   cor_mat,
#'   outPDF = "correlation_heatmap.pdf",
#'   column_title = "Dataset 1",
#'   row_title = "Dataset 2"
#' )
#' }
plot_correlation_heatmap_complex <- function(
    cor_matrix, outPDF = "correlation_analysis_results.pdf",
    column_title = "query cell type", row_title = "reference cell type",
    legend_title = "Pearson's\ncorrelation\ncoefficient",
    label = TRUE, labelSize = 5, labValueCutoff = 0.4, highLabCutoff = 0.6, 
    highLabCol = "white", row_names_rot = 0, column_names_rot = 45,
    diagonal_color = "gray80", plotHeight = 10, plotWidth = 12, 
    plotCommonCellType = TRUE) {
  
  # Check for required package
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Package 'ComplexHeatmap' is required. Install with: BiocManager::install('ComplexHeatmap')")
  }
  
  library(ComplexHeatmap)
  
  # Create color palette
  cor_values <- seq(0, 1, length.out = 100)
  power <- 2.5
  transformed_seq <- cor_values^power
  
  col_fun <- colorRampPalette(
    rev(hcl.colors(100, palette = "Blues", alpha = 0.9))
  )(100)[findInterval(transformed_seq, cor_values)]
  
  # Create heatmap with or without labels
  if (label) {
    ht <- Heatmap(
      matrix = cor_matrix,
      col = col_fun,
      na_col = "gray90",
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      heatmap_legend_param = list(
        title = legend_title,
        title_gp = gpar(fontsize = 16, fontface = "bold"),
        labels_gp = gpar(fontsize = 12),
        legend_height = unit(4, "cm")
      ),
      column_title = column_title,
      column_title_gp = gpar(fontsize = 18, fontface = "bold"),
      row_title = row_title,
      row_title_gp = gpar(fontsize = 18, fontface = "bold"),
      column_names_side = "bottom",
      row_names_side = "left",
      row_names_gp = gpar(fontsize = 10),
      row_names_rot = row_names_rot,
      column_names_gp = gpar(fontsize = 10),
      column_names_rot = column_names_rot,
      column_names_centered = TRUE,
      cell_fun = function(j, i, x, y, width, height, fill) {
        if (i == j) {
          grid.rect(x, y, width, height, 
                   gp = gpar(fill = diagonal_color, col = "white", lwd = 0.5))
        } else {
          if (!is.na(cor_matrix[i, j]) && abs(cor_matrix[i, j]) > labValueCutoff) {
            text_color <- ifelse(abs(cor_matrix[i, j]) >= highLabCutoff, 
                               highLabCol, "black")
            grid.text(sprintf("%.4f", cor_matrix[i, j]), x, y, 
                     gp = gpar(fontsize = labelSize, col = text_color))
          }
        }
      },
      border = TRUE,
      rect_gp = gpar(col = "white", lwd = 0.5)
    )
  } else {
    ht <- Heatmap(
      matrix = cor_matrix,
      col = col_fun,
      na_col = "gray90",
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      heatmap_legend_param = list(
        title = legend_title,
        title_gp = gpar(fontsize = 16, fontface = "bold"),
        labels_gp = gpar(fontsize = 12),
        legend_height = unit(4, "cm")
      ),
      column_title = column_title,
      column_title_gp = gpar(fontsize = 18, fontface = "bold"),
      row_title = row_title,
      row_title_gp = gpar(fontsize = 18, fontface = "bold"),
      column_names_side = "bottom",
      row_names_side = "left",
      row_names_rot = row_names_rot,
      row_names_gp = gpar(fontsize = 10),
      column_names_gp = gpar(fontsize = 10),
      column_names_rot = column_names_rot,
      column_names_centered = TRUE,
      cell_fun = function(j, i, x, y, width, height, fill) {
        if (i == j) {
          grid.rect(x, y, width, height, 
                   gp = gpar(fill = diagonal_color, col = "white", lwd = 0.5))
        }
      },
      border = TRUE,
      rect_gp = gpar(col = "white", lwd = 0.5)
    )
  }
  
  # Save PDF
  pdf(outPDF, height = plotHeight, width = plotWidth)
  draw(ht)
  dev.off()
  
  cat("Heatmap saved to:", outPDF, "\n")
}
