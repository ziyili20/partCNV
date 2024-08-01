#' Extract and normalize gene expression counts for a SingleCellExperiment object
#'
#' This function helps normalize the gene expression count matrix for a SingleCellExperiment object.
#'
#' @param obj          The SingleCellExperiment object.
#' @param scale_factor Feature counts for each cell are divided by the total counts for 
#' that cell and multiplied by the scale.factor, and then natural-log transformed using log1p.
#'
#' @return A normalized gene expression counts matrix. 
#' 
#' @importClassesFrom  SingleCellExperiment SingleCellExperiment
#' @export
#' @examples
#' data(SimDataSce)
#' counts_mat <- NormalizeCounts(SimDataSce)
#' 
#'
NormalizeCounts <- function(obj, scale_factor=10000) {
  counts <- obj@assays@data@listData[[1]]
  counts <- as.matrix(counts)
  counts_normalized <- log1p(scale_factor * sweep(counts, 2, colSums(counts), FUN="/"))
  return(counts_normalized)
}
