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
#' @importFrom CellChat createCellChat subsetData identifyOverExpressedGenes
#' @importFrom CellChat identifyOverExpressedInteractions computeCommunProb
#' @importFrom CellChat filterCommunication computeCommunProbPathway aggregateNet
#' @importFrom CellChat netVisual_circle subsetCommunication
#' @importFrom grDevices pdf dev.off
#' @importFrom graphics layout par
#' @importFrom utils write.csv
#' @importFrom rlang sym
#' @importFrom stats setNames
#' @importFrom methods is
#'
#' @export
run_cellChat <- function(seuObj, splitCellBy, groupBy, prefix) {
  # Check for CellChat package
  if (!requireNamespace("CellChat", quietly = TRUE)) {
    stop("Package 'CellChat' is required. Install with:\n",
         "devtools::install_github('sqjin/CellChat')")
  }
  
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
  CellChatDB.use <- CellChat
  
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
    subObj <- subset(seuObj, subset = (!!rlang::sym(splitCellBy) == currentCellGroup))
    
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


#' Draw bubble plot for CellChat ligand-receptor interactions
#'
#' @description
#' Creates a bubble plot visualization of cell-cell communication patterns from
#' CellChat analysis results. The plot shows ligand-receptor interactions with
#' bubble size representing statistical significance and color representing 
#' interaction probability.
#'
#' @param df Dataframe containing CellChat interaction results with columns:
#'   \itemize{
#'     \item source - Source cell type
#'     \item target - Target cell type
#'     \item ligand - Ligand gene name
#'     \item receptor - Receptor gene name
#'     \item interaction_name_2 - Ligand-receptor pair name
#'     \item pathway_name - Signaling pathway name
#'     \item prob - Interaction probability
#'     \item pval - P-value of interaction
#'   }
#' @param sourceItems Character vector, filter for specific source cell types (default: NULL, use all)
#' @param targetItems Character vector, filter for specific target cell types (default: NULL, use all)
#' @param pathway Character vector, filter for specific signaling pathways (default: NULL, use all)
#' @param ligands Character vector, filter for specific ligands (default: NULL, use all)
#' @param facetBy Character string, column name for faceting the plot (default: NULL, no faceting)
#' @param probCutoff Numeric, minimum interaction probability threshold (default: 0)
#' @param pValueCutoff Numeric, maximum p-value threshold (default: 1)
#' @param drawSingleLR Logical, whether to create individual plots for each L-R pair when faceting (default: FALSE)
#' @param prefix Character string, output filename prefix
#' @param pltHeight Numeric, plot height in inches
#' @param pltWidth Numeric, plot width in inches
#'
#' @return NULL. Saves bubble plot(s) to PDF file(s)
#'
#' @details
#' The bubble plot visualizes cell-cell communication with:
#' \itemize{
#'   \item \strong{X-axis}: Source-target cell type pairs
#'   \item \strong{Y-axis}: Ligand-receptor interaction pairs
#'   \item \strong{Bubble size}: Statistical significance (p-value categories)
#'     \itemize{
#'       \item Small: 0.01 < p
#'       \item Medium: 0.001 < p <= 0.01
#'       \item Large: p <= 0.001
#'     }
#'   \item \strong{Bubble color}: Interaction probability (gradient from blue to red)
#' }
#' 
#' Filtering options allow focused analysis on:
#' \itemize{
#'   \item Specific cell type pairs (source and/or target)
#'   \item Particular signaling pathways
#'   \item Selected ligands of interest
#' }
#' 
#' Output filename: \code{<prefix>_netVisual_bubble.pdf}
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_size_manual scale_color_gradientn
#' @importFrom ggplot2 facet_wrap theme_bw theme element_text element_rect labs
#' @importFrom ggplot2 margin
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette pdf dev.off
#' @importFrom dplyr case_when
#' @importFrom stats na.omit
#' @importFrom utils write.csv
#' @importFrom graphics par
#' @importFrom methods is
#' @importFrom rlang .data
#' 
#' @export
draw_CellChat_bubble <- function(df, sourceItems = NULL, targetItems = NULL, 
                                 pathway = NULL, ligands = NULL,
                                 facetBy = NULL, probCutoff = 0, pValueCutoff = 1, 
                                 drawSingleLR = FALSE,
                                 prefix, pltHeight, pltWidth) {
  
  # Input validation
  if (!is.data.frame(df)) {
    stop("df must be a dataframe")
  }
  
  required_cols <- c("source", "target", "interaction_name_2", "prob", "pval")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  if (!is.numeric(probCutoff) || probCutoff < 0 || probCutoff > 1) {
    stop("probCutoff must be a number between 0 and 1")
  }
  
  if (!is.numeric(pValueCutoff) || pValueCutoff < 0 || pValueCutoff > 1) {
    stop("pValueCutoff must be a number between 0 and 1")
  }
  
  # Define color palette
  n.colors <- 10
  color.use <- brewer.pal(n.colors, "Spectral")
  color.use <- rev(color.use)
  col_func <- colorRampPalette(color.use)(99)
  
  # Filter data
  subDF <- df[df$prob > probCutoff & df$pval < pValueCutoff, ]
  
  if (nrow(subDF) == 0) {
    warning("No interactions pass the filtering criteria. Try relaxing probCutoff or pValueCutoff.")
    return(invisible(NULL))
  }
  
  if (!is.null(ligands)) {
    if (!"ligand" %in% colnames(subDF)) {
      warning("Column 'ligand' not found in dataframe. Skipping ligand filtering.")
    } else {
      subDF <- subDF[subDF$ligand %in% ligands, ]
    }
  }
  
  if (!is.null(sourceItems)) subDF <- subDF[subDF$source %in% sourceItems, ]
  if (!is.null(targetItems)) subDF <- subDF[subDF$target %in% targetItems, ]
  if (!is.null(pathway)) {
    if (!"pathway_name" %in% colnames(subDF)) {
      warning("Column 'pathway_name' not found in dataframe. Skipping pathway filtering.")
    } else {
      subDF <- subDF[subDF$pathway_name %in% pathway, ]
    }
  }
  
  if (nrow(subDF) == 0) {
    warning("No interactions remain after filtering. Check filter parameters.")
    return(invisible(NULL))
  }
  
  # Create source-target label
  subDF$sourceTarget <- paste0(subDF$source, " -> ", subDF$target)
  
  # Create p-value size categories
  subDF$pSize <- case_when(
    subDF$pval <= 0.001 ~ "p <= 0.001",
    subDF$pval <= 0.01  ~ "0.001 < p <= 0.01",
    TRUE                ~ "0.01 < p"
  )
  subDF$pSize <- factor(subDF$pSize, 
                       levels = c("0.01 < p", "0.001 < p <= 0.01", "p <= 0.001"))
  
  # Create plot
  p <- ggplot(subDF, aes(x = sourceTarget, y = interaction_name_2)) +
    geom_point(aes(size = pSize, color = prob)) +
    scale_size_manual(values = c(1, 2, 3), name = "p-value") +
    scale_color_gradientn(colors = col_func, na.value = "grey90", 
                         name = "Interaction\nProbability")
  
  if (!is.null(facetBy)) {
    if (!facetBy %in% colnames(subDF)) {
      warning("Column '", facetBy, "' not found. Skipping faceting.")
    } else {
      p <- p + facet_wrap(aes(.data[[facetBy]]), 
                         scales = "free_x", 
                         ncol = length(unique(subDF[[facetBy]])))
    }
  }
  
  p <- p + theme_bw() +
    theme(
      strip.text = element_text(size = 10, face = "bold", color = "white", 
                               margin = margin(4, 4, 4, 4)),
      strip.background = element_rect(fill = "#4A90E2", color = "#4A90E2", 
                                      linewidth = 0.8),
      panel.border = element_rect(color = "#333333", fill = NA, linewidth = 0.8),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold")
    ) +
    labs(x = "Cell Type Pairs (Source -> Target)",
         y = "Ligand-Receptor Interactions",
         title = "Cell-Cell Communication Network")
  
  # Save plot
  outPDF <- paste0(prefix, "_netVisual_bubble.pdf")
  pdf(outPDF, height = pltHeight, width = pltWidth)
  print(p)
  dev.off()
  
  message("Bubble plot saved to: ", outPDF)
  message("Number of interactions plotted: ", nrow(subDF))
  
  # Draw individual L-R pair plots if requested
  if (!is.null(facetBy) & drawSingleLR) {
    if (!facetBy %in% colnames(subDF)) {
      warning("Cannot create single L-R plots: facetBy column not found")
      return(invisible(NULL))
    }
    
    LRpairs <- sort(unique(subDF$interaction_name_2))
    message("Creating ", length(LRpairs), " individual L-R pair plots...")
    
    for (i in seq_along(LRpairs)) {
      lr_data <- subDF[subDF$interaction_name_2 == LRpairs[i], ]
      
      p_single <- ggplot(lr_data, aes(x = sourceTarget, y = .data[[facetBy]])) +
        geom_point(aes(size = pSize, color = prob)) +
        scale_size_manual(values = c(1, 2, 3), name = "p-value") +
        scale_color_gradientn(colors = col_func, na.value = "grey90",
                            name = "Interaction\nProbability") +
        theme_bw() +
        theme(
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 12, face = "bold")
        ) +
        labs(x = "Cell Type Pairs (Source -> Target)",
             y = facetBy,
             title = LRpairs[i])
      
      # Clean filename
      safe_filename <- gsub("[^A-Za-z0-9_-]", "_", LRpairs[i])
      out_single <- paste0(prefix, "_LR_", safe_filename, ".pdf")
      
      pdf(out_single, height = pltHeight * 0.6, width = pltWidth * 0.8)
      print(p_single)
      dev.off()
    }
    
    message("Individual plots saved with prefix: ", prefix, "_LR_")
  }
  
  return(invisible(NULL))
}
