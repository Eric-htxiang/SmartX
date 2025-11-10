# =============================================================================#
# Visualization Functions
# =============================================================================#

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
#' @importFrom Seurat DimPlot LabelClusters
#' @importFrom ggplot2 labs theme element_text guide_legend
#' @importFrom grDevices pdf dev.off
#' @export
drawUMAP <- function(seuObj, groupBy, reducName, sampleIDs = NULL, colorCodes, drawSingleGroup = FALSE, prefix, pltHeight = 8, pltWidth = 8, pointSize = 0.1, legCol = NULL, labBox = TRUE, labColor = "white", labSize = 3) {
  # Input validation
  if (!inherits(seuObj, "Seurat")) {
    stop("seuObj must be a Seurat object")
  }
  if (!groupBy %in% colnames(seuObj@meta.data)) {
    stop(paste("Column", groupBy, "not found in object metadata"))
  }
  if (!reducName %in% names(seuObj@reductions)) {
    stop(paste("Reduction", reducName, "not found. Available:", paste(names(seuObj@reductions), collapse = ", ")))
  }

  groupIDs <- sort(unique(as.character(seuObj@meta.data[[groupBy]])))
  if (is.null(names(colorCodes))) names(colorCodes) <- groupIDs
  colorCodes[!is.na(names(colorCodes))]

  if (is.null(legCol)) legCol <- ceiling(length(groupIDs) / 40)

  outPDF <- paste0(prefix, "_", groupBy, "_", reducName, "_UMAP.pdf")
  grDevices::pdf(outPDF, height = pltHeight, width = pltWidth + legCol)

  # Main UMAP plot
  pltTitle <- paste0(prefix, " ", groupBy, " (", ncol(seuObj), " cells)")
  legTitle <- paste0(groupBy, " (", length(groupIDs), " groups)")
  p <- DimPlot(seuObj, reduction = reducName, pt.size = pointSize, group.by = groupBy, cols = colorCodes, raster = FALSE)
  p <- p + labs(title = pltTitle, color = legTitle)
  p <- p + guides(color = guide_legend(ncol = legCol, override.aes = list(size = 3)))
  p <- p + theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 18, face = "bold")
  )

  # Add cluster labels
  p <- LabelClusters(p, id = groupBy, color = labColor, size = labSize, box = labBox,
                     max.overlaps = 50, box.padding = 0.1, segment.color = NA)
  print(p)

  # Individual group plots
  if (drawSingleGroup) {
    for (i in seq(1, length(groupIDs))) {
      seuObj@meta.data$groupTmp <- as.character(seuObj@meta.data[[groupBy]])
      seuObj@meta.data[seuObj@meta.data[[groupBy]] != groupIDs[i], ]$groupTmp <- "Others"
      cellNum <- nrow(seuObj@meta.data[seuObj@meta.data[[groupBy]] == groupIDs[i], ])
      pltTitle <- paste0(prefix, " ", groupBy, " \"", groupIDs[i], "\" (", cellNum , " cells)")
      p <- DimPlot(seuObj, reduction = reducName, pt.size = pointSize, group.by = "groupTmp", cols = c("red", "grey"), raster = FALSE)
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
      p <- DimPlot(seuObj, cells = cells, reduction = reducName, pt.size = pointSize, group.by = groupBy, cols = colorCodes, raster = FALSE)
      p <- p + labs(title = paste0(prefix, " ", groupBy, " ", sampleID, " (", length(cells), " cells)"))
      p <- p + guides(color = guide_legend(ncol = legCol, override.aes = list(size = 3)))
      p <- p + theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
      p <- LabelClusters(p, id = groupBy, color = labColor, size = labSize, box = labBox,
                         max.overlaps = 50, box.padding = 0.1, segment.color = NA)
      print(p)
    }
  }

  grDevices::dev.off()
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
#' @importFrom Seurat DotPlot
#' @importFrom ggplot2 scale_y_discrete labs theme element_text
#' @importFrom grDevices pdf dev.off
#' @importFrom dplyr mutate
#' @export
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
    warning("The following genes are not in the dataset: ", paste(missing_genes, collapse = ", "))
  }
  geneLists <- intersect(geneLists, available_genes)
  if (length(geneLists) == 0) {
    stop("No valid genes to plot")
  }

  counts_df <- createCountsDf(seuObj, groupBy)
  seuObj@meta.data[[groupBy]] <- factor(seuObj@meta.data[[groupBy]], levels = counts_df$groupID)

  p <- DotPlot(seuObj, features = geneLists, group.by = groupBy)
  p <- p + scale_y_discrete(labels = counts_df$label)
  p <- p + labs(title = paste0(prefix, " (", groupBy, ") (", ncol(seuObj), " cells)"))
  p <- p + theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold")
  )

  outPDF <- paste0(prefix, "_", groupBy, "_givenGeneExpDotPlot.pdf")
  grDevices::pdf(outPDF, height = set_height(nrow(counts_df)), width = set_width(length(geneLists)))
  print(p)
  grDevices::dev.off()
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
#' @importFrom dplyr group_by summarise mutate ungroup n
#' @importFrom ggplot2 ggplot aes geom_col coord_flip scale_y_continuous
#' @importFrom ggplot2 scale_fill_manual labs guides guide_legend theme_bw
#' @importFrom ggplot2 element_blank element_text
#' @importFrom rlang sym
#' @importFrom grDevices pdf dev.off
#' @export
drawStackBar <- function(dataframe, xItem, yItem, orderedYimtes = NULL, colorCodes, 
                         prefix, pltHeight, pltWidth, legPos = "right", legCol = NULL) {
  df_summary <- dataframe %>%
    group_by(.data[[yItem]], .data[[xItem]]) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(.data[[yItem]]) %>%
    mutate(proportion = .data$count / sum(.data$count)) %>%
    ungroup()

  df_summary[[xItem]] <- factor(df_summary[[xItem]], levels = as.character(sort(unique(df_summary[[xItem]]))))

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
  p1 <- ggplot(df_summary, aes(x = !!sym(yItem), y = .data$count, fill = !!sym(xItem))) +
    geom_col() +
    coord_flip() +
    scale_y_continuous(expand = c(0.005, 0.005)) +
    scale_fill_manual(values = colorCodes) +
    labs(
      title = paste0("Stacked Bar Plot of Cell Counts by ", yItem, " and ", xItem),
      x = yItem, y = "Number of Cells", fill = xItem
    ) +
    guides(fill = guide_legend(ncol = legCol, override.aes = list(size = 3))) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = legPos,
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5)
    )

  # Proportion plot
  p2 <- ggplot(df_summary, aes(x = !!sym(yItem), y = .data$proportion, fill = !!sym(xItem))) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = colorCodes) +
    scale_y_continuous(expand = c(0.005, 0.005)) +
    labs(
      title = paste0("Stacked Bar Plot of Cell Proportion by ", yItem, " and ", xItem),
      x = yItem, y = "Proportion of Cells (%)", fill = xItem
    ) +
    guides(fill = guide_legend(ncol = legCol, override.aes = list(size = 3))) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = legPos,
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5)
    )

  outPDF <- paste0(prefix, "_", xItem, "_", yItem, "_cellComposition.pdf")
  grDevices::pdf(outPDF, height = pltHeight, width = pltWidth)
  print(p1)
  print(p2)
  grDevices::dev.off()
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
#' @param legCol Integer, number of legend columns (default: auto-calculated)
#' @param scale Character string, scaling method ("width", "area", or "count", default: "width")
#' @param plotBox Logical, whether to overlay box plots (default: FALSE)
#'
#' @return NULL. Saves PDF file with violin plots
#'
#' @importFrom Seurat FetchData
#' @importFrom ggplot2 ggplot aes geom_violin labs facet_wrap scale_fill_manual
#' @importFrom ggplot2 theme_minimal theme element_text element_rect
#' @importFrom ggplot2 geom_boxplot position_dodge margin
#' @importFrom reshape2 melt
#' @importFrom grDevices pdf dev.off
#' @importFrom pals alphabet alphabet2
#' @export
drawFeatureExprVlnPlot <- function(seuObj, layer = "data", geneLists, groupBy, prefix,
                                   pltHeight, pltWidth, legCol=NULL, scale = "width", 
                                   plotBox = FALSE) {
  if (!requireNamespace("pals", quietly = TRUE)) {
    stop("Package 'pals' is required. Please install it with: install.packages('pals')")
  }

  colorCodes <- as.vector(c(alphabet(), alphabet2()))
  colorCodes <- colorCodes[1:length(unique(seuObj@meta.data[[groupBy]]))]

  geneLists <- intersect(geneLists, rownames(seuObj))
  geneExpr <- FetchData(seuObj, layer = layer, vars = geneLists)
  geneExpr[[groupBy]] <- seuObj@meta.data[[groupBy]]
  expDF <- melt(geneExpr, id.var = groupBy)
  names(expDF) <- c(groupBy, "gene", "expression")

  p <- ggplot(expDF, aes(x = .data$expression, y = !!rlang::sym(groupBy), fill = !!rlang::sym(groupBy))) +
    geom_violin(alpha = 0.6, scale = scale) +
    labs(
      title = paste0(prefix, " (scale=", scale, ")"),
      x = "Expression Level", y = groupBy, fill = groupBy
    )

  if (plotBox) {
    p <- p + geom_boxplot(width = 0.15, alpha = 0.9, fill = "white", position = position_dodge(width = 0.9), outlier.size = 0.3)
  }
  if(!is.null(legCol)) {
    p <- p + guides(fill = guide_legend(ncol = legCol, override.aes = list(size = 3)))
  }

  p <- p + facet_wrap(~ .data$gene, scales = "free_x", ncol = length(levels(expDF$gene))) +
    scale_fill_manual(values = colorCodes) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 10, face = "bold", color = "white", margin = margin(4, 4, 4, 4)),
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

  outPDF <- paste0(prefix, "_", groupBy, "_geneExprVlnPlot.", scale, ".pdf")
  grDevices::pdf(outPDF, height = pltHeight, width = pltWidth)
  print(p)
  grDevices::dev.off()
  message("Violin plot saved to: ", outPDF)
}


#' Plot heatmap with ComplexHeatmap
#'
#' @description
#' Creates a publication-quality heatmap using the ComplexHeatmap
#' package. Supports extensive customization including diagonal masking,
#' automatic label sizing, and intersection group plotting.
#'
#' @param cor_matrix Numeric matrix, correlation matrix (from exprCorrelation or similar)
#' @param outPDF Character string, output PDF filename (default: "correlation_analysis_results.pdf")
#' @param column_title Character string, title for columns (default: "query cell type")
#' @param row_title Character string, title for rows (default: "reference cell type")
#' @param legend_title Character string, legend title (default: "Pearson's\\ncorrelation\\ncoefficient")
#' @param label Logical, whether to show correlation values in cells (default: TRUE)
#' @param labelSize Numeric, font size for cell labels (default: 5)
#' @param labValueCutoff Numeric, minimum correlation value to display label (default: 0.4)
#' @param highLabCutoff Numeric, threshold for using highLabCol color (default: 0.6)
#' @param highLabCol Character string, color for high correlation labels (default: "white")
#' @param row_names_rot Numeric, row name rotation angle in degrees (default: 0)
#' @param column_names_rot Numeric, column name rotation angle in degrees (default: 45)
#' @param row_names_size Numeric, font size for row names (default: 10)
#' @param column_names_size Numeric, font size for column names (default: 10)
#' @param column_names_max_height Unit object or character string:
#'   \itemize{
#'     \item NULL (default): uses unit(6, "cm")
#'     \item "auto": automatically calculates based on label lengths
#'     \item unit object: custom height specification
#'   }
#' @param plotIntersectGroup Logical, whether to plot only intersecting row/column names (default: TRUE)
#' @param hideDiagonal Logical, whether to mask diagonal cells (default: FALSE)
#' @param diagonal_color Character string, color for diagonal cells if hidden (default: "gray80")
#' @param plotHeight Numeric, plot height in inches (default: 10)
#' @param plotWidth Numeric, plot width in inches (default: 12)
#'
#' @return NULL. Saves heatmap to PDF file
#'
#' @importFrom ComplexHeatmap Heatmap draw
#' @importFrom grid gpar unit grid.rect grid.text
#' @importFrom grDevices colorRampPalette hcl.colors pdf dev.off
#' @export
drawHeatmap <- function(
    cor_matrix,
    outPDF = "correlation_analysis_results.pdf",
    column_title = "query cell type",
    row_title = "reference cell type",
    legend_title = "Pearson's\ncorrelation\ncoefficient",
    label = TRUE,
    labelSize = 5,
    labValueCutoff = 0.4,
    highLabCutoff = 0.6,
    highLabCol = "white",
    row_names_rot = 0,
    column_names_rot = 45,
    row_names_size = 10,
    column_names_size = 10,
    column_names_max_height = NULL,
    plotIntersectGroup = TRUE,
    hideDiagonal = FALSE,
    diagonal_color = "gray80",
    plotHeight = 10,
    plotWidth = 12
) {
  # Input validation
  if (!is.matrix(cor_matrix) && !is.data.frame(cor_matrix)) {
    stop("cor_matrix must be a matrix or dataframe")
  }
  if (!is.numeric(as.matrix(cor_matrix))) {
    stop("cor_matrix must contain numeric values")
  }
  if (any(!is.finite(as.matrix(cor_matrix)), na.rm = TRUE)) {
    warning("cor_matrix contains infinite values")
  }

  # Create color palette with power transformation
  cor_values <- seq(0, 1, length.out = 100)
  power <- 2.5
  transformed_seq <- cor_values^power
  col_fun <- colorRampPalette(rev(hcl.colors(100, palette = "Blues", alpha = 0.9)))(100)[findInterval(transformed_seq, cor_values)]

  # Set column names centering based on rotation
  column_names_centered <- if (column_names_rot == 0) TRUE else FALSE

  # Handle column names max height
  if (is.null(column_names_max_height)) {
    column_names_max_height <- unit(6, "cm")
  } else if (is.character(column_names_max_height) && column_names_max_height == "auto") {
    colLabels <- colnames(cor_matrix)
    column_names_max_height <- set_max_text_width(colLabels, fontsize = column_names_size)
    message("Auto-calculated column names max height: ", as.numeric(column_names_max_height), " mm")
  } else if (!inherits(column_names_max_height, "unit")) {
    stop("column_names_max_height must be NULL, 'auto', or a unit object")
  }

  # Plot only intersecting groups if requested
  if (plotIntersectGroup) {
    commonIDs <- intersect(rownames(cor_matrix), colnames(cor_matrix))
    if (length(commonIDs) == 0) {
      stop("No common IDs found between rows and columns")
    }
    if (length(commonIDs) < nrow(cor_matrix) || length(commonIDs) < ncol(cor_matrix)) {
      message("Plotting ", length(commonIDs), " common groups (intersection mode)")
    }
    cor_matrix <- cor_matrix[commonIDs, commonIDs, drop = FALSE]
  }

  # Create heatmap
  if (label) {
    ht <- Heatmap(
      matrix = cor_matrix,
      col = col_fun,
      na_col = "gray90",
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      # Legend Settings
      heatmap_legend_param = list(
        title = legend_title,
        title_gp = gpar(fontsize = 16, fontface = "bold"),
        labels_gp = gpar(fontsize = 12),
        legend_height = unit(4, "cm")
      ),
      # Title Setting
      column_title = column_title,
      column_title_gp = gpar(fontsize = 18, fontface = "bold"),
      row_title = row_title,
      row_title_gp = gpar(fontsize = 18, fontface = "bold"),
      # Label Setting
      row_names_side = "left",
      row_names_rot = row_names_rot,
      row_names_gp = gpar(fontsize = row_names_size),
      column_names_side = "bottom",
      column_names_rot = column_names_rot,
      column_names_centered = column_names_centered,
      column_names_gp = gpar(fontsize = column_names_size),
      column_names_max_height = column_names_max_height,
      cell_fun = function(j, i, x, y, width, height, fill) {
        # Handle diagonal cells
        if (i == j && hideDiagonal) {
          grid.rect(x, y, width, height, gp = gpar(fill = diagonal_color, col = "white", lwd = 0.5))
        } else {
          # Display values if above cutoff
          val <- cor_matrix[i, j]
          if (!is.na(val) && abs(val) > labValueCutoff) {
            text_color <- ifelse(abs(val) >= highLabCutoff, highLabCol, "black")
            grid.text(sprintf("%.4f", val), x, y, gp = gpar(fontsize = labelSize, col = text_color))
          }
        }
      },
      # Border settings
      border = TRUE,
      rect_gp = gpar(col = "white", lwd = 0.5)
    )
  } else {
    # Heatmap without labels
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
      row_names_side = "left",
      row_names_rot = row_names_rot,
      row_names_gp = gpar(fontsize = row_names_size),
      column_names_side = "bottom",
      column_names_rot = column_names_rot,
      column_names_centered = column_names_centered,
      column_names_gp = gpar(fontsize = column_names_size),
      column_names_max_height = column_names_max_height,
      cell_fun = function(j, i, x, y, width, height, fill) {
        if (i == j && hideDiagonal) {
          grid.rect(x, y, width, height, gp = gpar(fill = diagonal_color, col = "white", lwd = 0.5))
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
  message("Heatmap saved to: ", outPDF)
  message("Matrix dimensions: ", nrow(cor_matrix), " x ", ncol(cor_matrix))
}
