#' Run pseudo-bulk differential expression analysis with biological replicates
#'
#' @description
#' Performs pseudo-bulk differential gene expression analysis for single-cell
#' RNA-seq data with biological replicates. This function aggregates single-cell
#' data to pseudo-bulk samples and uses DESeq2 for statistical testing, which
#' properly accounts for biological variability between replicates.
#'
#' @param seuObj A Seurat object containing single-cell RNA-seq data
#' @param cellTypeCol Character string specifying the column name in metadata
#'   that contains cell type annotations (default: "CellType")
#' @param sampleCol Character string specifying the column name in metadata
#'   that contains sample IDs (default: "sampleID")
#' @param conditionCol Character string specifying the column name in metadata
#'   that contains condition/group information. If NULL, will be automatically
#'   extracted from sampleCol by removing numeric suffixes (default: NULL)
#' @param controlGroup Character string specifying the reference/control group
#'   name (e.g., "WT", "CTRL"). This will be used as the baseline for comparison
#' @param testGroup Character string specifying the test/treatment group name
#'   (e.g., "CKO", "KO"). If NULL, all non-control groups will be compared
#' @param cellTypes Character vector of specific cell types to analyze. If NULL,
#'   all cell types will be analyzed (default: NULL)
#' @param minCells Integer, minimum number of cells required per sample for a
#'   cell type to be included in analysis (default: 10)
#' @param minSamples Integer, minimum number of samples required per condition
#'   for analysis (default: 2)
#' @param assayName Character string specifying which assay to use (default: "RNA")
#' @param slot Character string specifying which slot to use. Must be "counts"
#'   for proper DESeq2 analysis (default: "counts")
#' @param padjThreshold Numeric, adjusted p-value threshold for significance
#'   (default: 0.05)
#' @param lfcThreshold Numeric, log2 fold change threshold for significance
#'   (default: 0)
#' @param verbose Logical, whether to print progress messages (default: TRUE)
#'
#' @return A list containing:
#' \describe{
#'   \item{results}{A dataframe with differential expression results for all
#'     cell types, containing columns:
#'     \itemize{
#'       \item cellType: Cell type name
#'       \item gene: Gene symbol
#'       \item baseMean: Mean of normalized counts across all samples
#'       \item log2FoldChange: Log2 fold change (test vs control)
#'       \item lfcSE: Standard error of log2FoldChange
#'       \item stat: Wald test statistic
#'       \item pvalue: Unadjusted p-value
#'       \item padj: Benjamini-Hochberg adjusted p-value
#'       \item significant: Logical indicating if gene meets significance criteria
#'       \item regulation: "UP" or "DOWN" for significant genes
#'     }}
#'   \item{summary}{A dataframe summarizing the number of DEGs per cell type}
#'   \item{dds_objects}{A list of DESeq2 dataset objects for each cell type
#'     (useful for downstream analyses like variance stabilization)}
#'   \item{pseudobulk_counts}{Aggregated count matrix used for analysis}
#'   \item{sample_metadata}{Metadata table describing samples and conditions}
#'   \item{parameters}{List of parameters used in the analysis}
#' }
#'
#' @details
#' This function implements the pseudo-bulk approach for differential expression
#' analysis of single-cell RNA-seq data with biological replicates. The key steps are:
#'
#' \enumerate{
#'   \item Aggregate single-cell counts to pseudo-bulk samples by combining cells
#'     from the same cell type and sample
#'   \item Create sample metadata linking each pseudo-bulk sample to its condition
#'   \item For each cell type, run DESeq2 analysis accounting for biological replicates
#'   \item Filter and annotate results based on significance thresholds
#' }
#'
#' This approach is superior to standard FindMarkers for experiments with biological
#' replicates because:
#' \itemize{
#'   \item It properly treats biological replicates as the unit of replication
#'   \item It avoids pseudo-replication issues (treating cells as independent)
#'   \item It uses appropriate statistical models (negative binomial) for count data
#'   \item It provides accurate p-values and controls false discovery rates
#' }
#'
#' @note
#' \itemize{
#'   \item Requires DESeq2 package to be installed
#'   \item Input must be raw counts (not normalized or log-transformed)
#'   \item Recommended to have at least 3 biological replicates per condition
#'   \item Cell types with too few cells or samples will be automatically excluded
#' }
#'
#' @references
#' Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of fold change
#' and dispersion for RNA-seq data with DESeq2. Genome Biology, 15:550.
#'
#' Squair, J.W. et al. (2021) Confronting false discoveries in single-cell
#' differential expression. Nature Communications, 12:5692.
#'
#' @examples
#' \dontrun{
#' # Basic usage with automatic condition detection
#' results <- run_pseudobulk_DEG(
#'   seuObj = seurat_obj,
#'   cellTypeCol = "CellType",
#'   sampleCol = "sampleID",
#'   controlGroup = "WT"
#' )
#'
#' # Analyze specific cell types only
#' results <- run_pseudobulk_DEG(
#'   seuObj = seurat_obj,
#'   cellTypeCol = "CellType",
#'   sampleCol = "sampleID",
#'   controlGroup = "WT",
#'   testGroup = "CKO",
#'   cellTypes = c("T_cells", "B_cells", "Macrophages")
#' )
#'
#' # Access results
#' deg_table <- results$results
#' summary_stats <- results$summary
#'
#' # Get significant genes for a specific cell type
#' tcell_degs <- subset(results$results,
#'                      cellType == "T_cells" & significant == TRUE)
#' }
#'
#' @importFrom Seurat AggregateExpression
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results resultsNames
#' @importFrom dplyr filter select mutate group_by summarise n left_join
#' @importFrom tidyr separate
#' @importFrom stats relevel
#' @export
run_pseudobulk_DEG <- function(seuObj,
                               cellTypeCol = "cellType",
                               sampleCol = "sampleID",
                               conditionCol = "condition",
                               controlGroup = "WT",
                               testGroup = NULL,
                               cellTypes = NULL,
                               minCells = 10,
                               minSamples = 2,
                               assayName = "RNA",
                               slot = "counts",
                               lowExprFilter = FALSE,
                               padjThreshold = 0.05,
                               lfcThreshold = 0,
                               verbose = TRUE) {

  # ============================================================================
  # 1. Input validation
  # ============================================================================
  if (!inherits(seuObj, "Seurat")) {
    stop("seuObj must be a Seurat object")
  }

  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("Package 'DESeq2' is required but not installed. Please install it using:\n",
         "BiocManager::install('DESeq2')")
  }

  required_cols <- c(cellTypeCol, sampleCol, conditionCol)
  missing_cols <- setdiff(required_cols, colnames(seuObj@meta.data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in metadata: ", paste(missing_cols, collapse = ", "))
  }

  if (slot != "counts") {
    warning("slot should be 'counts' for DESeq2 analysis. Using normalized or scaled data may produce incorrect results.")
  }

  # ============================================================================
  # 2. Determine test group if not specified
  # ============================================================================
  if (is.null(testGroup)) {
    all_conditions <- unique(seuObj@meta.data[[conditionCol]])
    testGroup <- setdiff(all_conditions, controlGroup)
    if (length(testGroup) == 0) {
      stop("No test group found. All samples belong to control group.")
    }
    if (length(testGroup) > 1) {
      warning("Multiple test groups detected: ", paste(testGroup, collapse = ", "),
              "\nAnalyzing all of them against control.")
    }
  }

  # ============================================================================
  # 3. Filter cell types if specified
  # ============================================================================
  if (!is.null(cellTypes)) {
    available_types <- unique(seuObj@meta.data[[cellTypeCol]])
    missing_types <- setdiff(cellTypes, available_types)
    if (length(missing_types) > 0) {
      warning("Cell types not found in data: ", paste(missing_types, collapse = ", "))
    }
    cellTypes <- intersect(cellTypes, available_types)
  } else {
    cellTypes <- unique(seuObj@meta.data[[cellTypeCol]])
  }

  if (verbose) {
    cat("Analyzing", length(cellTypes), "cell type(s)\n")
  }

  # ============================================================================
  # 4. Create aggregation grouping variable
  # ============================================================================
  seuObj@meta.data$celltype_sample <- paste(
    seuObj@meta.data[[cellTypeCol]],
    seuObj@meta.data[[sampleCol]],
    sep = "---"  # Use distinctive separator
  )

  # ============================================================================
  # 5. Aggregate to pseudo-bulk
  # ============================================================================
  if (verbose) cat("Aggregating cells to pseudo-bulk samples...\n")

  pseudobulk <- Seurat::AggregateExpression(
    seuObj,
    group.by = "celltype_sample",
    assays = assayName,
    slot = slot,
    return.seurat = FALSE
  )

  counts_matrix <- pseudobulk[[assayName]]

  if (verbose) {
    cat("Pseudo-bulk matrix dimensions:", nrow(counts_matrix), "genes x",
        ncol(counts_matrix), "samples\n")
  }

  # ============================================================================
  # 6. Create sample metadata
  # ============================================================================
  sample_info <- data.frame(
    sample_id = colnames(counts_matrix),
    stringsAsFactors = FALSE
  )

  # Parse celltype and sample from combined ID
  sample_info$cellType <- sapply(strsplit(sample_info$sample_id, "---"),
                                   function(x) x[1])
  sample_info$sample <- sapply(strsplit(sample_info$sample_id, "---"),
                                function(x) x[2])

  # Add condition information
  condition_map <- unique(seuObj@meta.data[, c(sampleCol, conditionCol)])
  colnames(condition_map) <- c("sample", "condition")

  sample_info <- merge(sample_info, condition_map, by = "sample", all.x = TRUE)
  rownames(sample_info) <- sample_info$sample_id

  # Add cell counts per pseudo-bulk sample
  cell_counts <- table(seuObj@meta.data$celltype_sample)
  sample_info$nCells <- as.numeric(cell_counts[sample_info$sample_id])

  # ============================================================================
  # 7. Run DESeq2 for each cell type
  # ============================================================================
  all_results <- list()
  dds_list <- list()
  summary_list <- list()

  for (ct in cellTypes) {
    if (verbose) {
      cat("\n", strrep("=", 70), "\n")
      cat("Processing cell type:", ct, "\n")
      cat(strrep("=", 70), "\n")
    }

    # Filter samples for this cell type
    ct_samples <- sample_info[sample_info$cellType == ct, ]

    # Check minimum cells requirement
    if (any(ct_samples$nCells < minCells)) {
      low_samples <- ct_samples$sample[ct_samples$nCells < minCells]
      if (verbose) {
        cat("Warning: Samples with < ", minCells, " cells: ",
            paste(low_samples, collapse = ", "), "\n")
      }
      ct_samples <- ct_samples[ct_samples$nCells >= minCells, ]
    }

    # Check minimum samples per condition
    samples_per_condition <- table(ct_samples$condition)
    if (any(samples_per_condition < minSamples)) {
      if (verbose) {
        cat("Skipping", ct, "- insufficient samples per condition\n")
        print(samples_per_condition)
      }
      next
    }

    # Filter only control and test groups
    relevant_conditions <- c(controlGroup, testGroup)
    ct_samples <- ct_samples[ct_samples$condition %in% relevant_conditions, ]

    if (nrow(ct_samples) < 4) {  # At least 2 samples per group
      if (verbose) cat("Skipping", ct, "- too few samples after filtering\n")
      next
    }

    # Extract count matrix for this cell type
    ct_counts <- counts_matrix[, ct_samples$sample_id, drop = FALSE]

    # Filter low-expressed genes
    if (lowExprFilter) {
      keep_genes <- rowSums(ct_counts >= 10) >= 2
      ct_counts <- ct_counts[keep_genes, , drop = FALSE]
    }
    
    if (verbose) {
      cat("Samples:", nrow(ct_samples), "\n")
      cat("Genes (after filtering):", nrow(ct_counts), "\n")
      print(table(ct_samples$condition))
    }

    # Create DESeq2 dataset
    tryCatch({
      dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = ct_counts,
        colData = ct_samples,
        design = ~ condition
      )

      # Set reference level
      dds$condition <- relevel(factor(dds$condition), ref = controlGroup)

      # Run DESeq2
      if (verbose) cat("Running DESeq2 analysis...\n")
      dds <- DESeq2::DESeq(dds, quiet = !verbose)

      # Store DESeq2 object
      dds_list[[ct]] <- dds

      # Extract results for each test group vs control
      for (test_cond in testGroup) {
        contrast_name <- paste0(test_cond, "_vs_", controlGroup)

        res <- DESeq2::results(
          dds,
          contrast = c("condition", test_cond, controlGroup),
          alpha = padjThreshold
        )

        # Convert to dataframe
        res_df <- as.data.frame(res)
        res_df$gene <- rownames(res_df)
        res_df$cellType <- ct
        res_df$contrast <- contrast_name
        res_df$testGroup <- test_cond
        res_df$controlGroup <- controlGroup

        # Add significance and regulation direction
        res_df$significant <- !is.na(res_df$padj) &
                              res_df$padj < padjThreshold &
                              abs(res_df$log2FoldChange) > lfcThreshold

        res_df$regulation <- NA
        res_df$regulation[res_df$significant & res_df$log2FoldChange > 0] <- "UP"
        res_df$regulation[res_df$significant & res_df$log2FoldChange < 0] <- "DOWN"

        # Reorder columns
        res_df <- res_df[, c("cellType", "contrast", "gene", "baseMean",
                              "log2FoldChange", "lfcSE", "stat",
                              "pvalue", "padj", "significant", "regulation",
                              "testGroup", "controlGroup")]

        all_results[[paste0(ct, "___", contrast_name)]] <- res_df

        # Summary statistics
        n_up <- sum(res_df$regulation == "UP", na.rm = TRUE)
        n_down <- sum(res_df$regulation == "DOWN", na.rm = TRUE)
        n_total <- n_up + n_down

        if (verbose) {
          cat(sprintf("  %s: %d DEGs (%d UP, %d DOWN)\n",
                      contrast_name, n_total, n_up, n_down))
        }

        summary_list[[paste0(ct, "___", contrast_name)]] <- data.frame(
          cellType = ct,
          contrast = contrast_name,
          nGenes_tested = nrow(res_df),
          nDEGs = n_total,
          nUP = n_up,
          nDOWN = n_down,
          nSamples_control = sum(ct_samples$condition == controlGroup),
          nSamples_test = sum(ct_samples$condition == test_cond),
          stringsAsFactors = FALSE
        )
      }

    }, error = function(e) {
      if (verbose) {
        cat("Error in DESeq2 analysis for", ct, ":", conditionMessage(e), "\n")
      }
    })
  }

  # ============================================================================
  # 8. Combine results
  # ============================================================================
  if (length(all_results) == 0) {
    stop("No successful DESeq2 analyses. Check your data and parameters.")
  }

  results_df <- do.call(rbind, all_results)
  rownames(results_df) <- NULL

  summary_df <- do.call(rbind, summary_list)
  rownames(summary_df) <- NULL

  # ============================================================================
  # 9. Prepare output
  # ============================================================================
  output <- list(
    results = results_df,
    summary = summary_df,
    dds_objects = dds_list,
    pseudobulk_counts = counts_matrix,
    sample_metadata = sample_info,
    parameters = list(
      cellTypeCol = cellTypeCol,
      sampleCol = sampleCol,
      conditionCol = conditionCol,
      controlGroup = controlGroup,
      testGroup = testGroup,
      minCells = minCells,
      minSamples = minSamples,
      padjThreshold = padjThreshold,
      lfcThreshold = lfcThreshold,
      date = Sys.time()
    )
  )

  class(output) <- c("pseudobulk_DEG", "list")

  if (verbose) {
    cat("\n", strrep("=", 70), "\n")
    cat("Analysis complete!\n")
    cat("Analyzed", length(unique(results_df$cellType)), "cell types\n")
    cat("Total DEGs:", sum(results_df$significant, na.rm = TRUE), "\n")
    cat(strrep("=", 70), "\n\n")
  }

  return(output)
}


#' Print method for pseudobulk_DEG objects
#'
#' @param x A pseudobulk_DEG object
#' @param ... Additional arguments (not used)
#' @export
print.pseudobulk_DEG <- function(x, ...) {
  cat("Pseudo-bulk Differential Expression Analysis\n")
  cat(strrep("=", 50), "\n\n")

  cat("Parameters:\n")
  cat("  Control group:", x$parameters$controlGroup, "\n")
  cat("  Test group(s):", paste(x$parameters$testGroup, collapse = ", "), "\n")
  cat("  Padj threshold:", x$parameters$padjThreshold, "\n")
  cat("  LFC threshold:", x$parameters$lfcThreshold, "\n\n")

  cat("Results:\n")
  cat("  Cell types analyzed:", length(unique(x$results$cellType)), "\n")
  cat("  Total genes tested:", nrow(x$results), "\n")
  cat("  Significant DEGs:", sum(x$results$significant, na.rm = TRUE), "\n\n")

  cat("Summary by cell type:\n")
  print(x$summary[, c("cellType", "contrast", "nDEGs", "nUP", "nDOWN")])

  cat("\nAccess results with: $results, $summary, $dds_objects\n")
}


#' Summary method for pseudobulk_DEG objects
#'
#' @param object A pseudobulk_DEG object
#' @param ... Additional arguments (not used)
#' @export
summary.pseudobulk_DEG <- function(object, ...) {
  cat("Pseudo-bulk DEG Analysis Summary\n")
  cat(strrep("=", 50), "\n\n")

  print(object$summary)

  cat("\nTop 10 genes by adjusted p-value:\n")
  top_genes <- object$results[order(object$results$padj), ]
  top_genes <- top_genes[top_genes$significant == TRUE, ]
  print(head(top_genes[, c("cellType", "gene", "log2FoldChange", "padj")], 10))
}