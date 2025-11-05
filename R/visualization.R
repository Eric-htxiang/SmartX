# =============================================================================
# Visualization Functions
# =============================================================================

#' Draw UMAP plots for Seurat object
#'
#' @description
#' Creates comprehensive UMAP visualizations including overall clustering,
#' individual cluster highlighting, and sample-specific plots.
#'
#' @param seuObj A Seurat object with dimensional reduction computed
#' @param groupBy Character string specifying metadata column for grouping
#' @param reducName Character string specifying reduction name (e.g., "umap", "tsne")
#' @param sampleIDs Optional character vector of sample IDs to plot separately
#' @param colorCodes Named vector of colors for each group
#' @param drawSingleGroup Logical, whether to create individual plots for each group (default: FALSE)
#' @param prefix Character string for output filename prefix
#' @param pltHeight Numeric, plot height in inches (default: 8)
#' @param pltWidth Numeric, plot width in inches (default: 8)
#' @param pointSize Numeric, point size for cells (default: 0.1)
#' @param legCol Integer, number of legend columns (default: auto-calculated)
#' @param labBox Logical, whether to draw boxes around labels (default: TRUE)
#' @param labColor Character string, label color (default: "white")
#' @param labSize Numeric, label size (default: 3)
#'
#' @return NULL. Saves PDF file with UMAP plots
#'
#' @details
#' The function creates three types of plots:
#' \enumerate{
#'   \item Overall UMAP with all groups colored
#'   \item Individual plots highlighting each group (if drawSingleGroup = TRUE)
#'   \item Sample-specific plots (if sampleIDs provided)
#' }
#'
#' Output filename: \code{<prefix>_<groupBy>_<reducName>_UMAP.pdf}
#'
#' @export
#' @examples
#' \dontrun{
#' # Basic usage
#' colors <- c("0" = "red", "1" = "blue", "2" = "green")
#' drawUMAP(seurat_obj, 
#'          groupBy = "seurat_clusters",
#'          reducName = "umap",
#'          colorCodes = colors,
#'          prefix = "MyProject")
#' 
#' # With individual group plots
#' drawUMAP(seurat_obj,
#'          groupBy = "cell_type",
#'          reducName = "umap",
#'          colorCodes = colors,
#'          drawSingleGroup = TRUE,
#'          prefix = "CellTypes")
#' }
drawUMAP <- function(seuObj, groupBy, reducName, sampleIDs = NULL, colorCodes,
                     drawSingleGroup = FALSE, prefix, pltHeight = 8, pltWidth = 8,
                     pointSize = 0.1, legCol = NULL, labBox = TRUE, labColor = "white", 
                     labSize = 3) {
  
  # Input validation
  if (!inherits(seuObj, "Seurat")) {
    stop("seuObj must be a Seurat object")
  }
  
  if (!groupBy %in% colnames(seuObj@meta.data)) {
    stop(paste("Column", groupBy, "not found in object metadata"))
  }
  
  if (!reducName %in% names(seuObj@reductions)) {
    stop(paste("Reduction", reducName, "not found. Available:", 
               paste(names(seuObj@reductions), collapse = ", ")))
  }
  
  groupIDs <- sort(unique(as.character(seuObj@meta.data[[groupBy]])))
  if (is.null(names(colorCodes))) names(colorCodes) <- groupIDs
  colorCodes[!is.na(names(colorCodes))]
  
  if (is.null(legCol)) legCol <- ceiling(length(groupIDs) / 40)
  
  outPDF <- paste0(prefix, "_", groupBy, "_", reducName, "_UMAP.pdf")
  pdf(outPDF, height = pltHeight, width = pltWidth + legCol)
  
  # Main UMAP plot
  pltTitle <- paste0(prefix, "  ", groupBy, "  (", ncol(seuObj), " cells)")
  legTitle <- paste0(groupBy, " (", length(groupIDs), " groups)")
  p <- DimPlot(seuObj, reduction = reducName, pt.size = pointSize, group.by = groupBy, 
               cols = colorCodes, raster = FALSE)
  p <- p + labs(title = pltTitle, color = legTitle)
  p <- p + guides(color = guide_legend(ncol = legCol, override.aes = list(size = 3)))
  p <- p + theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
                 legend.title = element_text(size = 18, face = "bold"))
  
  p <- LabelClusters(p, id = groupBy, color = labColor, size = labSize, box = labBox, 
                     max.overlaps = 50, box.padding = 0.1, segment.color = NA)
  print(p)
  
  # Individual group plots
  if (drawSingleGroup) {
    for (i in seq(1, length(groupIDs))) {
      seuObj@meta.data$groupTmp <- seuObj@meta.data[[groupBy]]
      seuObj@meta.data[seuObj@meta.data[[groupBy]] != groupIDs[i], ]$groupTmp <- "Others"
      cellNum <- nrow(seuObj@meta.data[seuObj@meta.data[[groupBy]] == groupIDs[i], ])
      pltTitle <- paste0(prefix, "  ", groupBy, "  \"",  groupIDs[i], "\"  (", cellNum , " cells)")
      p <- DimPlot(seuObj, reduction = reducName, pt.size = pointSize, 
                   group.by = "groupTmp", cols = c("red", "grey"), raster = FALSE)
      p <- p + labs(title = pltTitle)
      p <- p + guides(color = guide_legend(ncol = legCol, override.aes = list(size = 3)))
      p <- p + theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
      print(p)
    }
  }
  
  # Sample-specific plots
  if (!is.null(sampleIDs)) {
    for (sampleID in sampleIDs) {
      cells <- rownames(seuObj@meta.data[seuObj@meta.data$sample == sampleID, ])
      p <- DimPlot(seuObj, cells = cells, reduction = reducName, pt.size = pointSize, 
                   group.by = groupBy, cols = colorCodes, raster = FALSE)
      p <- p + labs(title = paste0(prefix, "  ", groupBy, "  ", sampleID, "  (", length(cells), " cells)"))
      p <- p + guides(color = guide_legend(ncol = legCol, override.aes = list(size = 3)))
      p <- p + theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
      p <- LabelClusters(p, id = groupBy, color = labColor, size = labSize, box = labBox, 
                         max.overlaps = 50, box.padding = 0.1, segment.color = NA)
      print(p)
    }
  }
  dev.off()
  
  message("UMAP plots saved to: ", outPDF)
}


#' Draw gene expression dot plot
#'
#' @description
#' Creates a dot plot showing gene expression across groups with cell count
#' annotations on y-axis labels.
#'
#' @param seuObj A Seurat object
#' @param groupBy Character string specifying metadata column for grouping
#' @param geneLists Character vector of gene names to plot
#' @param prefix Character string for output filename prefix
#'
#' @return NULL. Saves PDF file with dot plot
#'
#' @details
#' The function creates a dot plot where:
#' \itemize{
#'   \item Dot size represents percentage of cells expressing the gene
#'   \item Dot color represents average expression level
#'   \item Y-axis labels show percentage and cell count for each group
#' }
#'
#' Output filename: \code{<prefix>_<groupBy>_givenGeneExpDotPlot.pdf}
#'
#' @export
#' @examples
#' \dontrun{
#' genes <- c("CD3D", "CD4", "CD8A", "MS4A1")
#' drawExprDotPlot(seurat_obj,
#'                 groupBy = "seurat_clusters",
#'                 geneLists = genes,
#'                 prefix = "MarkerGenes")
#' }
drawExprDotPlot <- function(seuObj, groupBy, geneLists, prefix) {
  # Input validation
  if (!inherits(seuObj, "Seurat")) {
    stop("seuObj must be a Seurat object")
  }
  
  if (!groupBy %in% colnames(seuObj@meta.data)) {
    stop(paste("Column", groupBy, "not found in object metadata"))
  }
  
  # Check for missing genes
  available_genes <- rownames(seuObj)
  missing_genes <- setdiff(geneLists, available_genes)
  if (length(missing_genes) > 0) {
    warning("The following genes are not in the dataset: ", 
            paste(missing_genes, collapse = ", "))
    geneLists <- intersect(geneLists, available_genes)
  }
  
  if (length(geneLists) == 0) {
    stop("No valid genes to plot")
  }
  
  counts_df <- createCountsDf(seuObj, groupBy)
  seuObj@meta.data[[groupBy]] <- factor(seuObj@meta.data[[groupBy]], levels = counts_df$groupID)
  
  p <- DotPlot(seuObj, features = geneLists, group.by = groupBy)
  p <- p + scale_y_discrete(labels = counts_df$label)
  p <- p + labs(title = paste0(prefix, " (", groupBy, ") (", ncol(seuObj), " cells)"))
  p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                 plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))
  
  outPDF <- paste0(prefix, "_", groupBy, "_givenGeneExpDotPlot.pdf")
  pdf(outPDF, height = set_height(nrow(counts_df)), width = set_width(length(geneLists)))
  print(p)
  dev.off()
  
  message("Dot plot saved to: ", outPDF)
}


#' Draw stacked bar chart for cell composition
#'
#' @description
#' Creates stacked bar charts showing cell counts and proportions across groups.
#'
#' @param dataframe A dataframe containing the data
#' @param xItem Character string specifying column for x-axis (fill colors)
#' @param yItem Character string specifying column for y-axis (bars)
#' @param orderedYimtes Optional character vector specifying order of y-axis items
#' @param colorCodes Named vector of colors
#' @param prefix Character string for output filename prefix
#' @param pltHeight Numeric, plot height in inches
#' @param pltWidth Numeric, plot width in inches
#' @param legPos Character string, legend position ("right" or "bottom", default: "right")
#' @param legCol Integer, number of legend columns (default: auto-calculated)
#'
#' @return NULL. Saves PDF file with two stacked bar plots (counts and proportions)
#'
#' @export
#' @examples
#' \dontrun{
#' colors <- setNames(rainbow(5), c("A", "B", "C", "D", "E"))
#' drawStackBar(meta.data,
#'              xItem = "cell_type",
#'              yItem = "sample",
#'              colorCodes = colors,
#'              prefix = "CellComposition",
#'              pltHeight = 6,
#'              pltWidth = 10)
#' }
drawStackBar <- function(dataframe, xItem, yItem, orderedYimtes = NULL, colorCodes, 
                         prefix, pltHeight, pltWidth, legPos = "right", legCol = NULL) {
  
  df_summary <- dataframe %>%
    group_by(.data[[yItem]], .data[[xItem]]) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(.data[[yItem]]) %>%
    mutate(proportion = count / sum(count)) %>%
    ungroup()
  
  df_summary[[xItem]] <- factor(df_summary[[xItem]], 
                                levels = as.character(sort(unique(df_summary[[xItem]]))))
  if (!is.null(orderedYimtes)) {
    df_summary[[yItem]] <- factor(df_summary[[yItem]], levels = orderedYimtes)
  }
  
  if (is.null(legCol)) {
    groupNum <- length(unique(as.character(df_summary[[xItem]])))
    if (legPos == "right") {
      legCol <- ceiling(groupNum / 20)
    } else if (legPos == "bottom") {
      legCol <- ceiling(groupNum / 10)
    }
  }
  
  # Count plot
  p1 <- ggplot(df_summary, aes(x = !!sym(yItem), y = count, fill = !!sym(xItem))) +
    geom_col() +
    coord_flip() +
    scale_y_continuous(expand = c(0.005, 0.005)) +
    scale_fill_manual(values = colorCodes) +
    labs(title = paste0("Stacked Bar Plot of Cell Counts by ", yItem, " and ", xItem),
         x = yItem, y = "Number of Cells", fill = xItem) +
    guides(fill = guide_legend(ncol = legCol, override.aes = list(size = 3))) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          legend.position = legPos,
          plot.title = element_text(size = 15, face = "bold", hjust = 0.5))
  
  # Proportion plot
  p2 <- ggplot(df_summary, aes(x = !!sym(yItem), y = proportion, fill = !!sym(xItem))) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = colorCodes) +
    scale_y_continuous(expand = c(0.005, 0.005)) +
    labs(title = paste0("Stacked Bar Plot of Cell Proportion by ", yItem, " and ", xItem),
         x = yItem, y = "Proportion of Cells (%)", fill = xItem) +
    guides(fill = guide_legend(ncol = legCol, override.aes = list(size = 3))) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          legend.position = legPos,
          plot.title = element_text(size = 15, face = "bold", hjust = 0.5))
  
  outPDF <- paste0(prefix, "_", xItem, "_", yItem, "_cellComposition.pdf")
  pdf(outPDF, height = pltHeight, width = pltWidth)
  print(p1)
  print(p2)
  dev.off()
  
  message("Stacked bar plots saved to: ", outPDF)
}


#' Draw gene expression violin plots
#'
#' @description
#' Creates violin plots showing gene expression distributions across groups.
#'
#' @param seuObj A Seurat object
#' @param layer Character string specifying data layer (default: "data")
#' @param geneLists Character vector of gene names to plot
#' @param groupBy Character string specifying metadata column for grouping
#' @param prefix Character string for output filename prefix
#' @param pltHeight Numeric, plot height in inches
#' @param pltWidth Numeric, plot width in inches
#' @param scale Character string, scaling method ("width", "area", or "count", default: "width")
#' @param plotBox Logical, whether to overlay box plots (default: FALSE)
#'
#' @return NULL. Saves PDF file with violin plots
#'
#' @export
#' @examples
#' \dontrun{
#' genes <- c("CD3D", "CD4", "CD8A")
#' drawFeatureExprVlnPlot(seurat_obj,
#'                        geneLists = genes,
#'                        groupBy = "seurat_clusters",
#'                        prefix = "MarkerExpression",
#'                        pltHeight = 8,
#'                        pltWidth = 12)
#' }
drawFeatureExprVlnPlot <- function(seuObj, layer = "data", geneLists, groupBy, 
                                   prefix, pltHeight, pltWidth, scale = "width", 
                                   plotBox = FALSE) {
  
  if (!requireNamespace("pals", quietly = TRUE)) {
    stop("Package 'pals' is required. Please install it with: install.packages('pals')")
  }
  
  library(pals)
  colorCodes <- as.vector(c(alphabet(), alphabet2()))
  colorCodes <- colorCodes[1:length(unique(seuObj@meta.data[[groupBy]]))]
  
  geneExpr <- FetchData(seuObj, layer = layer, vars = geneLists)
  geneExpr[[groupBy]] <- seuObj@meta.data[[groupBy]]
  expDF <- melt(geneExpr, id.var = groupBy)
  names(expDF) <- c(groupBy, "gene", "expression")
  
  p <- ggplot(expDF, aes(x = expression, y = !!sym(groupBy), fill = !!sym(groupBy))) +
    geom_violin(alpha = 0.6, scale = scale) +
    labs(title = prefix, x = "Expression Level", y = groupBy, fill = groupBy)
  
  if (plotBox) {
    p <- p + geom_boxplot(width = 0.15, alpha = 0.9, fill = "white",
                          position = position_dodge(width = 0.9),
                          outlier.size = 0.3)
  }
  
  p <- p + facet_wrap(~ gene, scales = "free_x", ncol = length(levels(expDF$gene))) +
    scale_fill_manual(values = colorCodes) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 10, face = "bold", color = "white", 
                               margin = margin(4, 4, 4, 4)),
      strip.background = element_rect(fill = "#4A90E2", color = "#4A90E2", linewidth = 0.8),
      panel.border = element_rect(color = "#333333", fill = NA, linewidth = 0.8),
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 25, face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(size = 22, face = "bold"),
      legend.text = element_text(size = 15)
    )
  
  outPDF <- paste0(prefix, "_geneExprVlnPlot.pdf")
  pdf(outPDF, height = pltHeight, width = pltWidth)
  print(p)
  dev.off()
  
  message("Violin plot saved to: ", outPDF)
}
