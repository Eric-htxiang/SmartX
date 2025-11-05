# =============================================================================
# CellChat Wrapper Functions
# =============================================================================

#' Run CellChat cell-cell communication analysis
#'
#' @description
#' Performs cell-cell communication analysis using CellChat package. This wrapper
#' function streamlines the CellChat workflow and creates comprehensive visualizations.
#'
#' @param seuObj A Seurat object with normalized data
#' @param splitCellBy Character string specifying metadata column to split analysis by
#' @param groupBy Character string specifying metadata column for cell type grouping
#' @param prefix Character string for output filename prefix
#'
#' @return NULL. Saves multiple output files:
#' \itemize{
#'   \item RDS files with CellChat objects
#'   \item CSV files with interaction matrices
#'   \item PDF files with network visualizations
#' }
#'
#' @details
#' The function performs the following steps for each group defined by splitCellBy:
#' \enumerate{
#'   \item Creates CellChat object from Seurat data
#'   \item Identifies over-expressed genes and interactions
#'   \item Computes communication probabilities
#'   \item Filters low-quality communications
#'   \item Aggregates pathway-level communications
#'   \item Creates network visualizations
#'   \item Saves results in multiple formats
#' }
#'
#' Requires CellChat package to be installed:
#' \code{devtools::install_github("sqjin/CellChat")}
#'
#' @export
#' @examples
#' \dontrun{
#' # Run CellChat analysis for each sample
#' run_cellChat(seurat_obj,
#'              splitCellBy = "sample",
#'              groupBy = "cell_type",
#'              prefix = "MyProject")
#' 
#' # The function will create:
#' # - MyProject_cell_type_Sample1_cellchat_object.rds
#' # - MyProject_cell_type_Sample1_cellchat_interaction_count.csv
#' # - MyProject_cell_type_Sample1_cellchat_interaction_weight.csv
#' # - MyProject_cell_type_Sample1_cellchat_network.pdf
#' # ... and more
#' }
run_cellChat <- function(seuObj, splitCellBy, groupBy, prefix) {
  # Check for CellChat package
  if (!requireNamespace("CellChat", quietly = TRUE)) {
    stop("Package 'CellChat' is required. Install with:\n",
         "devtools::install_github('sqjin/CellChat')")
  }
  
  library(CellChat)
  
  # Input validation
  if (!inherits(seuObj, "Seurat")) {
    stop("seuObj must be a Seurat object")
  }
  
  if (!splitCellBy %in% colnames(seuObj@meta.data)) {
    stop(paste("Column", splitCellBy, "not found in object metadata"))
  }
  
  if (!groupBy %in% colnames(seuObj@meta.data)) {
    stop(paste("Column", groupBy, "not found in object metadata"))
  }
  
  # Initialize
  CellChatDB <- CellChatDB.mouse  # Change to CellChatDB.human for human data
  CellChatDB.use <- CellChatDB
  
  cellchat_list <- list()
  cellGroupIDs <- unique(as.character(seuObj@meta.data[[splitCellBy]]))
  
  # Main loop through each group
  for (i in seq_along(cellGroupIDs)) {
    currentCellGroup <- cellGroupIDs[i]
    currentPfx <- paste0(prefix, "_", groupBy, "_", currentCellGroup)
    
    # Timestamp
    formatted_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat("\n", rep("=", 80), "\n", sep = "")
    cat("[", formatted_time, "]  Processing: \"", groupBy, "\" | \"", 
        currentCellGroup, "\"\n", sep = "")
    cat(rep("=", 80), "\n", sep = "")
    
    # Subset data
    subObj <- subset(seuObj, subset = (!!sym(splitCellBy) == currentCellGroup))
    
    # Check cell numbers
    n_cells <- ncol(subObj)
    cell_types <- table(subObj@meta.data[[groupBy]])
    cat("  Total cells: ", n_cells, "\n")
    cat("  Cell types:\n")
    print(cell_types)
    
    # Quality checks
    if (n_cells < 50) {
      cat("  Warning: Too few cells (< 50), skipping...\n")
      next
    }
    
    min_cells_per_type <- min(cell_types)
    if (min_cells_per_type < 10) {
      cat("  Warning: Some cell types have < 10 cells\n")
    }
    
    # CellChat analysis
    tryCatch({
      # Create CellChat object
      cat("  Creating CellChat object...\n")
      cellchat_obj <- createCellChat(
        object = subObj, 
        group.by = groupBy,
        assay = "SCT"
      )
      
      cellchat_obj@DB <- CellChatDB.use
      
      # Data preprocessing
      cat("  Step 1: Subsetting data...\n")
      cellchat_obj <- subsetData(cellchat_obj)
      
      cat("  Step 2: Identifying over-expressed genes...\n")
      cellchat_obj <- identifyOverExpressedGenes(cellchat_obj)
      
      cat("  Step 3: Identifying over-expressed interactions...\n")
      cellchat_obj <- identifyOverExpressedInteractions(cellchat_obj)
      
      # Compute communication probability
      cat("  Step 4: Computing communication probability...\n")
      cellchat_obj <- computeCommunProb(cellchat_obj, type = "triMean")
      
      # Filter low-quality communications
      cat("  Step 5: Filtering communication...\n")
      cellchat_obj <- filterCommunication(cellchat_obj, min.cells = 10)
      
      # Pathway-level analysis
      cat("  Step 6: Computing pathway-level communication...\n")
      cellchat_obj <- computeCommunProbPathway(cellchat_obj)
      
      cat("  Step 7: Aggregating network...\n")
      cellchat_obj <- aggregateNet(cellchat_obj)
      
      # Save results
      cat("  Step 8: Saving results...\n")
      
      # Save network matrices
      write.csv(cellchat_obj@net$count, 
                paste0(currentPfx, "_cellchat_interaction_count.csv"))
      write.csv(cellchat_obj@net$weight, 
                paste0(currentPfx, "_cellchat_interaction_weight.csv"))
      
      # Save detailed interactions
      all_interactions <- subsetCommunication(cellchat_obj)
      write.csv(all_interactions, 
                paste0(currentPfx, "_cellchat_detailed_interactions.csv"), 
                row.names = FALSE)
      
      # Save CellChat object
      saveRDS(cellchat_obj, paste0(currentPfx, "_cellchat_object.rds"))
      
      # Add to list
      list_name <- paste(groupBy, currentCellGroup, sep = "_")
      cellchat_list[[list_name]] <- cellchat_obj
      
      # Visualization
      cat("  Step 9: Generating visualizations...\n")
      groupSize <- as.numeric(table(cellchat_obj@idents))
      
      # Overall network plots
      outPDF <- paste0(currentPfx, "_cellchat_network.pdf")
      pdf(outPDF, width = 12, height = 6)
      layout(matrix(c(1, 2), 1, 2))
      par(mar = c(2, 2, 2, 2))
      
      netVisual_circle(cellchat_obj@net$count, 
                       vertex.weight = groupSize, 
                       weight.scale = TRUE, 
                       label.edge = FALSE, 
                       title.name = "Number of interactions")
      
      netVisual_circle(cellchat_obj@net$weight, 
                       vertex.weight = groupSize, 
                       weight.scale = TRUE,
                       label.edge = FALSE, 
                       title.name = "Interaction weights/strength")
      dev.off()
      
      # Individual cell type network plots
      mat <- cellchat_obj@net$weight
      n_cell_types <- nrow(mat)
      
      if (n_cell_types > 0) {
        outPDF <- paste0(currentPfx, "_eachCellGroup_cellchat_network.pdf")
        
        # Auto-adjust layout
        n_cols <- min(4, ceiling(sqrt(n_cell_types)))
        n_rows <- ceiling(n_cell_types / n_cols)
        
        pdf(outPDF, width = 6 * n_cols, height = 4 * n_rows)
        par(mfrow = c(n_rows, n_cols), xpd = TRUE, mar = c(1, 1, 2, 1))
        
        for (cell_idx in 1:n_cell_types) {
          mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), 
                         dimnames = dimnames(mat))
          mat2[cell_idx, ] <- mat[cell_idx, ]
          
          netVisual_circle(mat2, 
                           vertex.weight = groupSize, 
                           weight.scale = TRUE, 
                           edge.weight.max = max(mat), 
                           title.name = rownames(mat)[cell_idx])
        }
        dev.off()
      }
      
      cat("  Analysis completed successfully!\n")
      
    }, error = function(e) {
      cat("  Error occurred: ", e$message, "\n")
      cat("  Skipping this combination...\n")
    })
  }
  
  # Save complete list
  cat("\n", rep("=", 80), "\n", sep = "")
  cat("Saving complete cellchat list...\n")
  saveRDS(cellchat_list, paste0(prefix, "_all_cellchat_objects.rds"))
  cat("Total objects created: ", length(cellchat_list), "\n")
  cat("Object names:\n")
  print(names(cellchat_list))
  cat(rep("=", 80), "\n", sep = "")
  
  message("\nCellChat analysis complete!")
  message("Results saved with prefix: ", prefix)
}
