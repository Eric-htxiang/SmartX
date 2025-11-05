# =============================================================================
# Data Processing Functions
# =============================================================================

#' Create counts dataframe from Seurat object
#'
#' @description
#' Creates a summary dataframe of cell counts grouped by a specified metadata 
#' column. The function automatically sorts group IDs (numeric first, then 
#' alphabetic) and places "Undefined" and "Unknown" categories at the end.
#'
#' @param obj A Seurat object
#' @param groupBy Character string specifying the metadata column name to group by
#'
#' @return A dataframe with the following columns:
#' \describe{
#'   \item{groupID}{Factor of group identifiers, properly sorted}
#'   \item{cellCounts}{Numeric count of cells in each group}
#'   \item{percent}{Percentage of total cells in each group (rounded to 3 decimals)}
#'   \item{label}{Formatted label combining percentage, count, and group ID}
#' }
#'
#' @details
#' The sorting logic is:
#' \enumerate{
#'   \item Numeric group IDs (sorted numerically)
#'   \item Alphabetic group IDs (sorted alphabetically)
#'   \item "Undefined" and "Unknown" categories (if present)
#' }
#'
#' @export
#' @examples
#' \dontrun{
#' library(Seurat)
#' # Assuming you have a Seurat object named 'seurat_obj'
#' counts_df <- createCountsDf(seurat_obj, "seurat_clusters")
#' head(counts_df)
#' 
#' # Example output:
#' #   groupID cellCounts percent              label
#' # 1       0       1234   15.4  (15.4% 1234) 0
#' # 2       1       2345   29.3  (29.3% 2345) 1
#' }
createCountsDf <- function(obj, groupBy) {
  # Input validation
  if (!inherits(obj, "Seurat")) {
    stop("obj must be a Seurat object")
  }
  
  if (!groupBy %in% colnames(obj@meta.data)) {
    stop(paste("Column", groupBy, "not found in object metadata"))
  }
  
  # Create counts dataframe
  counts_df <- as.data.frame(table(obj@meta.data[[groupBy]]))
  names(counts_df) <- c("groupID", "cellCounts")
  allGroupIDs <- unique(as.character(counts_df$groupID))
  
  # Sort numeric IDs
  sortedNum <- allGroupIDs %>%
    .[grepl("^-?\\d+$", .)] %>%
    as.numeric() %>%
    sort() %>%
    as.character()
  
  # Sort character IDs
  sortedChar <- allGroupIDs %>%
    .[!grepl("^-?\\d+$", .)] %>%
    sort()
  
  # Combine: numeric first, then character
  allGroupIDs <- c(sortedNum, sortedChar)
  
  # Move Undefined/Unknown to end
  undef <- c("Undefined", "Unknown")
  part1 <- setdiff(allGroupIDs, undef)
  part2 <- intersect(allGroupIDs, undef)
  allGroupIDs <- c(part1, part2)
  
  # Convert to factor with proper levels
  counts_df$groupID <- factor(counts_df$groupID, levels = allGroupIDs)
  counts_df <- counts_df %>% arrange(groupID)
  
  # Add percentage and label
  counts_df <- counts_df %>%
    mutate(percent = round(cellCounts / sum(cellCounts) * 100, digits = 3)) %>%
    mutate(label = paste0("(", percent, "% ", cellCounts, ") ", groupID))
  
  return(counts_df)
}


#' Set plot height based on number of clusters
#'
#' @description
#' Calculates appropriate plot height in inches based on the number of 
#' cluster IDs to ensure readable visualizations.
#'
#' @param n_clstIDs Integer, number of cluster IDs
#'
#' @return Numeric value representing recommended plot height in inches
#'
#' @details
#' Height scaling:
#' \itemize{
#'   \item <= 10 clusters: 6 inches
#'   \item 11-15 clusters: 7 inches
#'   \item 16-20 clusters: 8 inches
#'   \item > 20 clusters: 8 + 2 * (steps after 10) inches
#' }
#'
#' @export
#' @examples
#' set_height(8)   # Returns 6
#' set_height(12)  # Returns 7
#' set_height(25)  # Returns 10
set_height <- function(n_clstIDs) {
  if (!is.numeric(n_clstIDs) || n_clstIDs < 1) {
    stop("n_clstIDs must be a positive number")
  }
  
  if (n_clstIDs <= 10) {
    return(6)
  } else if (n_clstIDs > 10 & n_clstIDs <= 15) {
    return(7)
  } else if (n_clstIDs > 15 & n_clstIDs <= 20) {
    return(8)
  } else {
    base_height <- 8 
    step_height <- 2 
    steps <- ceiling((n_clstIDs - 10) / 10)
    return(base_height + (steps - 1) * step_height)
  }
}


#' Set plot width based on number of genes
#'
#' @description
#' Calculates appropriate plot width in inches based on the number of genes
#' using monotonic spline interpolation for smooth scaling.
#'
#' @param num_genes Integer, number of genes to display
#'
#' @return Numeric value representing recommended plot width in inches
#'
#' @details
#' Uses predefined control points and monotonic Hermite spline interpolation
#' to ensure smooth, monotonic width scaling. Minimum width is 7 inches.
#'
#' @export
#' @examples
#' set_width(5)    # Returns 7
#' set_width(10)   # Returns ~9
#' set_width(50)   # Returns ~18
set_width <- function(num_genes) {
  if (!is.numeric(num_genes) || num_genes < 1) {
    stop("num_genes must be a positive number")
  }
  
  geneNum_points <- c(5, 6, 8, 10, 15, 20, 30, 40, 50, 60, 70, 80, 100, 150, 200, 250, 300)
  width_points <- c(7, 8, 8, 9, 11, 12, 13, 15, 18, 20, 21, 22, 25, 35, 40, 55, 60)
  
  if (num_genes <= min(geneNum_points)) {
    return(min(width_points))
  } else {
    spline_fun <- splinefun(geneNum_points, width_points, method = "monoH.FC")
    return(spline_fun(num_genes))
  }
}
