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
  expr1 <- expr1[common_genes, , drop = FALSE]
  expr2 <- expr2[common_genes, , drop = FALSE]
  
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
#' @importFrom dplyr group_by arrange desc slice_head
#' @importFrom utils write.table
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
    markers1, markers2, groupBy1, groupBy2, groupIDs1_order=NULL, groupIDs2_order=NULL,
    pCutoff = 0.05, FCcutoff = 0.5, pct1Cutoff = 0.1, pct2Cutoff = 0.01,
    topGeneNum = 100, prefix) {
  print(groupIDs1_order)
  print(groupIDs2_order)
  # Filter DEGs
  df1 <- markers1[markers1$p_val_adj <= pCutoff &
                    markers1$avg_log2FC > FCcutoff &
                    markers1$pct.1 > pct1Cutoff &
                    markers1$pct.2 > pct2Cutoff, ]
  df2 <- markers2[markers2$p_val_adj <= pCutoff &
                    markers2$avg_log2FC > FCcutoff &
                    markers2$pct.1 > pct1Cutoff &
                    markers2$pct.2 > pct2Cutoff, ]
  
  # Select top genes per group
  df1 <- df1 %>%
    group_by(.data[[groupBy1]]) %>%
    arrange(desc(avg_log2FC), desc(pct.1)) %>%
    slice_head(n = topGeneNum)
  df2 <- df2 %>%
    group_by(.data[[groupBy2]]) %>%
    arrange(desc(avg_log2FC), desc(pct.1)) %>%
    slice_head(n = topGeneNum)
  
  #groupIDs1 <- as.character(sort(unique(df1[[groupBy1]])))
  #groupIDs2 <- as.character(sort(unique(df2[[groupBy2]])))
  
  if(is.null(groupIDs1_order)) {
    groupIDs1 <- unique(as.character(df1[[groupBy1]]))
    is_numeric <- grepl("^-?\\d+$", groupIDs1)
    groupIDs1 <- c(
      groupIDs1[is_numeric] %>% as.numeric() %>% sort() %>% as.character(),
      groupIDs1[!is_numeric] %>% sort()
    )
  }else {
    groupIDs1 <- groupIDs1_order
  }
  
  if(is.null(groupIDs2_order)) {
    groupIDs2 <- unique(as.character(df2[[groupBy2]]))
    is_numeric <- grepl("^-?\\d+$", groupIDs2)
    groupIDs2 <- c(
      groupIDs2[is_numeric] %>% as.numeric() %>% sort() %>% as.character(),
      groupIDs2[!is_numeric] %>% sort()
    )
  }else {
    groupIDs2 <- groupIDs2_order
  }
  
  # Initialize similarity matrix
  similarity_matrix <- matrix(0, nrow = length(groupIDs1), ncol = length(groupIDs2))
  rownames(similarity_matrix) <- groupIDs1
  colnames(similarity_matrix) <- groupIDs2
  
  # Calculate Jaccard similarity and collect overlaps
  overlapDEG_df <- data.frame()
  for (i in seq_along(groupIDs1)) {
    for (j in seq_along(groupIDs2)) {
      topDEG1 <- df1[df1[[groupBy1]] == groupIDs1[i], ]$gene
      topDEG2 <- df2[df2[[groupBy2]] == groupIDs2[j], ]$gene
      
      intersectSet <- intersect(topDEG1, topDEG2)
      unionSet <- union(topDEG1, topDEG2)
      
      similarity_matrix[i, j] <- ifelse(length(unionSet) == 0, 0,
                                        length(intersectSet) / length(unionSet))
      
      if (length(intersectSet) > 0) {
        items <- c(groupIDs1[i], groupIDs2[j],
                   as.character(length(intersectSet)),
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
#' @importFrom dplyr filter arrange desc group_by slice_head
#'
#' @export
#' @examples
#' \dontrun{
#' go_results <- run_GOanalysis(markers, groupBy = "seurat_clusters",
#'                              prefix = "MyProject")
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
  for (i in seq_along(groupIDs)) {
    for (j in seq_along(groupIDs)) {
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


