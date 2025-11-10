# =============================================================================#
# Differential Expression and GO Analysis Functions
# =============================================================================#

#' Run differential expression analysis
#'
#' @description
#' Performs differential gene expression analysis using Seurat's FindMarkers
#' or FindAllMarkers functions. Supports both simple cluster comparisons and
#' complex pairwise comparisons within groups.
#'
#' @param seuObj A Seurat object
#' @param groupBy Character string specifying the primary grouping variable
#' @param groupIDs Optional character vector of specific groups to analyze
#' @param compareBy Optional character string for secondary grouping (e.g., "treatment")
#' @param compareDF Optional dataframe with two columns specifying pairwise comparisons
#' @param testUse Character string specifying statistical test (default: "wilcox")
#' @param recorrectUMI Logical, whether to recorrect UMI counts (default: FALSE)
#' @param onlyPos Logical, return only positive markers (default: TRUE)
#' @param assayName Character string specifying assay to use (default: "SCT")
#'
#' @return A dataframe containing differential expression results with columns:
#' \describe{
#'   \item{gene}{Gene name}
#'   \item{p_val}{Unadjusted p-value}
#'   \item{p_val_adj}{Adjusted p-value}
#'   \item{avg_log2FC}{Average log2 fold change}
#'   \item{pct.1}{Percentage of cells expressing in group 1}
#'   \item{pct.2}{Percentage of cells expressing in group 2}
#'   \item{log2pctDiff}{Log2 ratio of pct.1 to pct.2}
#' }
#'
#' @details
#' Two modes of operation:
#' \enumerate{
#'   \item Simple mode: When compareBy and compareDF are NULL, runs FindAllMarkers
#'   \item Complex mode: Performs pairwise comparisons specified in compareDF within each group
#' }
#'
#' @importFrom Seurat FindMarkers FindAllMarkers Idents
#' @importFrom dplyr select
#' @importFrom rlang sym
#' @importFrom stats na.omit
#' @export
run_DEG <- function(seuObj, groupBy, groupIDs = NULL, compareBy = NULL, compareDF = NULL, 
                    testUse = "wilcox", recorrectUMI = FALSE, onlyPos = TRUE, assayName = "SCT") {
  # Input validation
  if (!inherits(seuObj, "Seurat")) {
    stop("seuObj must be a Seurat object")
  }
  if (!groupBy %in% colnames(seuObj@meta.data)) {
    stop(paste("Column", groupBy, "not found in object metadata"))
  }

  Idents(seuObj) <- seuObj@meta.data[[groupBy]]

  if (is.null(groupIDs)) groupIDs <- sort(unique(seuObj@meta.data[[groupBy]]))
  allMarkers <- data.frame()

  # Complex comparison mode
  if (!is.null(compareBy) & !is.null(compareDF)) {
    if (!compareBy %in% colnames(seuObj@meta.data)) {
      stop(paste("Column", compareBy, "not found in object metadata"))
    }
    if (ncol(compareDF) != 2) {
      stop("compareDF must have exactly 2 columns")
    }

    for (i in seq(1, length(groupIDs))) {
      subObj <- subset(seuObj, subset = !!sym(groupBy) == groupIDs[i])
      cat("Processing group:", groupIDs[i], "(", ncol(subObj), "cells )\n")

      for (j in seq(1, nrow(compareDF))) {
        compareA <- compareDF[, 1][j]
        compareB <- compareDF[, 2][j]

        if (!all(c(compareA, compareB) %in% subObj@meta.data[[compareBy]])) next

        formatted_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        cat("[", formatted_time, "] Starting FindMarkers: ", groupIDs[i], " - ", compareA, " vs ", compareB, "\n", sep = "")

        cellNumA <- nrow(subObj@meta.data[subObj@meta.data[[compareBy]] == compareA, ])
        cellNumB <- nrow(subObj@meta.data[subObj@meta.data[[compareBy]] == compareB, ])

        if (cellNumA >= 3 & cellNumB >= 3) {
          df <- FindMarkers(subObj, ident.1 = compareA, ident.2 = compareB, group.by = compareBy, test.use = testUse, recorrect_umi = recorrectUMI)
          df$log2pctDiff <- log2(df$pct.1 / df$pct.2)
          df$gene <- row.names(df)
          df[[groupBy]] <- groupIDs[i]
          df$sampleA <- compareA
          df$sampleB <- compareB
          df$cellNumA <- cellNumA
          df$cellNumB <- cellNumB
          df <- df %>% select(all_of(groupBy), sampleA, sampleB, cellNumA, cellNumB, gene, p_val, p_val_adj, avg_log2FC, pct.1, pct.2, log2pctDiff)
          df$compare <- paste0(df$sampleA, ".vs.", df$sampleB)
          df$direct <- "UP"
          df[df$avg_log2FC < 0, ]$direct <- "DOWN"
          df$group <- paste0(df$compare, "_", df$direct)
          allMarkers <- rbind(allMarkers, df)
        }
      }
    }
  } else {
    # Simple mode: FindAllMarkers
    cat("Running FindAllMarkers for", groupBy, "\n")
    df <- FindAllMarkers(seuObj, only.pos = onlyPos, assay = assayName, recorrect_umi = recorrectUMI)
    df$log2pctDiff <- log2(df$pct.1 / df$pct.2)
    df <- df %>% dplyr::select(cluster, gene, p_val, p_val_adj, avg_log2FC, pct.1, pct.2, log2pctDiff)
    
    colnames(df)[colnames(df) == "cluster"] <- groupBy
    allMarkers <- df
  }

  return(allMarkers)
}


#' Run Gene Ontology enrichment analysis
#'
#' @description
#' Performs GO enrichment analysis on differentially expressed genes and creates
#' integrated volcano plots with GO dot plots.
#'
#' @param df Dataframe of DEG results from run_DEG
#' @param groupBy Character string specifying grouping column
#' @param pvalCol Character string specifying p-value column (default: 'p_val_adj')
#' @param pCutoff Numeric, adjusted p-value cutoff (default: 0.05)
#' @param FCcutoff Numeric, log2 fold change cutoff (default: 0.5)
#' @param pct1Cutoff Numeric, minimum percentage in group 1 (default: 0.15)
#' @param pct2Cutoff Numeric, minimum percentage in group 2 (default: 0.01)
#' @param relaxFilter Logical, use relaxed filtering if few genes pass (default: FALSE)
#' @param topGene2Display Integer, number of top genes to label on volcano plot (default: 20)
#' @param topGenes4GO Integer, maximum genes to use for GO analysis (default: NULL)
#' @param prefix Character string for output filename prefix
#' @param pltHeight Numeric, plot height in inches (default: 8)
#' @param pltWidth Numeric, plot width in inches (default: 20)
#'
#' @return A dataframe containing GO enrichment results
#'
#' @details
#' For each group, creates a two-panel plot:
#' \enumerate{
#'   \item Left panel: Volcano plot of differential expression
#'   \item Right panel: GO term enrichment dot plot
#' }
#'
#' Requires clusterProfiler and org.Mm.eg.db packages.
#'
#' @importFrom clusterProfiler enrichGO
#' @importFrom enrichplot dotplot
#' @importFrom org.Mm.eg.db org.Mm.eg.db
#' @importFrom ggplot2 ggplot geom_hline geom_vline geom_point geom_text labs scale_colour_manual scale_x_continuous scale_y_continuous theme_minimal theme element_blank element_rect
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr filter arrange slice_head
#' @importFrom grDevices pdf dev.off
#' @importFrom scales alpha
#' @importFrom patchwork plot_layout
#' @export
run_GOanalysis <- function(df, groupBy, pvalCol="p_val_adj", pCutoff = 0.05, FCcutoff = 0.5, pct1Cutoff = 0.15,
                           pct2Cutoff = 0.01, relaxFilter = FALSE, topGene2Display = 20, topGenes4GO = NULL, 
                           prefix, pltHeight = 8, pltWidth = 20) {
  # Check for required packages
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("Package 'clusterProfiler' is required. Install with: BiocManager::install('clusterProfiler')")
  }
  if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
    stop("Package 'org.Mm.eg.db' is required. Install with: BiocManager::install('org.Mm.eg.db')")
  }
  if (!requireNamespace("enrichplot", quietly = TRUE)) {
    stop("Package 'enrichplot' is required. Install with: BiocManager::install('enrichplot')")
  }

  # Filter DEGs
  filDF <- df[df[[pvalCol]] <= pCutoff & df$avg_log2FC > FCcutoff & df$pct.1 > pct1Cutoff & df$pct.2 > pct2Cutoff, ]
  outPDF <- paste0(prefix, "_", groupBy, "_DEG.GO_plot.pdf")
  message("Saving GO plots to: ", outPDF)

  pdf(outPDF, height = pltHeight, width = pltWidth)

  groupIDs <- unique(as.character(df[[groupBy]]))
  is_numeric <- grepl("^-?\\d+$", groupIDs)
  groupIDs <- c(
    groupIDs[is_numeric] %>% as.numeric() %>% sort() %>% as.character(),
    groupIDs[!is_numeric] %>% sort()
  )
  
  allGOresults <- data.frame()

  for (i in seq(1, length(groupIDs))) {
    subDF <- df[df[[groupBy]] == groupIDs[i], ]
    filSubDF <- filDF[filDF[[groupBy]] == groupIDs[i], ]

    # Relaxed filtering if needed
    if (relaxFilter) {
      if (nrow(filSubDF) == 0) {
        cat("No genes passed filter for", groupIDs[i], "- using relaxed criteria\n")
        filSubDF <- subDF %>%
          filter(p_val < pCutoff) %>%
          arrange(desc(avg_log2FC)) %>%
          slice_head(n = 100)
      } else if (nrow(filSubDF) > 0 & nrow(filSubDF) < 20) {
        cat("Only", nrow(filSubDF), "genes for", groupIDs[i], "- adding more\n")
        tmpDF <- subDF %>%
          filter(!gene %in% filSubDF$gene) %>%
          arrange(desc(avg_log2FC)) %>%
          slice_head(n = 100 - nrow(filSubDF))
        filSubDF <- rbind(filSubDF, tmpDF)
      }
    }

    cat("Number of DEGs for", groupBy, "_", groupIDs[i], ":", nrow(filSubDF), "\n")
    if (nrow(filSubDF) == 0) next

    formatted_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat("[", formatted_time, "] Processing", groupBy, "_", groupIDs[i], "(", i, "of", length(groupIDs), ")\n", sep="")

    # Create volcano plot
    subDF$label <- "Un_significant"
    subDF[subDF$avg_log2FC > FCcutoff, ]$label <- "Up_regulated"
    if (nrow(subDF[subDF$avg_log2FC < -FCcutoff, ]) > 0) {
      subDF[subDF$avg_log2FC < -FCcutoff, ]$label <- "Down_regulated"
    }

    topDF <- filSubDF %>% slice_head(n = as.numeric(topGene2Display))

    myColors <- c("Up_regulated" = "#8400ff", "Down_regulated" = "#dab3ff", "Un_significant" = "#b9b9b9")
    xRange <- range(subDF$avg_log2FC, na.rm = TRUE)
    yRange <- range(-log10(subDF$p_val_adj)[is.finite(-log10(subDF$p_val_adj))], na.rm = TRUE)
    xBreak <- pretty(c(-ceiling(max(abs(xRange))), 0, ceiling(max(abs(xRange)))), n = 5)
    yBreak <- pretty(c(0, ceiling(max(yRange))), n = 5)

    volcano_plot <- ggplot() +
      geom_hline(yintercept = c(-log10(pCutoff)), linewidth = 0.7, color = "#0059FF", lty = "dashed") +
      geom_vline(xintercept = c(-FCcutoff, FCcutoff), linewidth = 0.7, color = "#0059FF", lty = "dashed") +
      geom_point(data = subDF, aes(x = avg_log2FC, y = -log10(p_val_adj), color = label), size = 2) +
      geom_text(aes(x = xBreak[1] * 0.95, y = -log10(pCutoff), label = pCutoff), color = "red", vjust = 1.2, hjust = 0, size = 3) +
      geom_text(aes(x = -FCcutoff, y = max(-log10(subDF$p_val_adj), na.rm = TRUE), label = -FCcutoff), angle = 90, vjust = -0.5, hjust = 1, color = "red", size = 3) +
      geom_text(aes(x = FCcutoff, y = max(-log10(subDF$p_val_adj), na.rm = TRUE), label = FCcutoff), angle = 90, vjust = 1.5, hjust = 1, color = "red", size = 3) +
      labs(x = 'avg_log2FC', y = '-log10(p_val_adj)') +
      scale_colour_manual(name = '', values = alpha(myColors, 0.7)) +
      scale_x_continuous(limits = range(xBreak), breaks = xBreak) +
      scale_y_continuous(limits = range(yBreak), breaks = yBreak) +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
        axis.text = element_text(color = "black", size = 10)
      )

    if (nrow(topDF) > 0) {
      volcano_plot <- volcano_plot +
        geom_text_repel(
          data = topDF,
          aes(x = avg_log2FC, y = -log10(p_val_adj), label = gene),
          color = "black",
          max.overlaps = 50,
          force = 2,
          size = 3
        )
    }

    # GO enrichment
    if (!is.null(topGenes4GO)) filSubDF <- filSubDF %>% slice_head(n = as.numeric(topGenes4GO))
    geneSymbols <- filSubDF$gene

    ego <- enrichGO(
      gene = geneSymbols,
      OrgDb = org.Mm.eg.db,
      keyType = 'SYMBOL',
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      readable = TRUE
    )

    egoDF <- as.data.frame(ego)
    if (nrow(egoDF) == 0) {
      cat("No GO terms enriched for", groupIDs[i], "\n")
      next
    }

    egoDF[[groupBy]] <- groupIDs[i]
    allGOresults <- rbind(allGOresults, egoDF)

    plotTitle <- paste0(prefix, " ", groupBy, "_", groupIDs[i])
    GO_dotplot <- dotplot(ego, showCategory = 30, label_format = 50) +
      ggtitle(plotTitle) +
      theme(plot.title = element_text(color = "#d63031", face = "bold", hjust = 0.5))

    p <- volcano_plot + GO_dotplot + plot_layout(ncol = 2, widths = c(1, 1))
    print(p)
  }

  dev.off()

  # Save GO results
  outFile <- paste0(prefix, "_", groupBy, "_GOresults.tsv")
  write.table(allGOresults, outFile, sep = "\t", quote = FALSE, row.names = FALSE)
  message("GO results saved to: ", outFile)

  return(allGOresults)
}
