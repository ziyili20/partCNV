#' Get normalized gene expression counts for selected genes
#'
#' This function helps normalize the gene expression count matrix if needed
#' and select the genes that are located in the interested region. This procedure
#' happens after applying GetCytoLocation().
#'
#' @param cytoloc_output The output from the function GetCytoLocation().
#' The function needs to be run with a complete cytogenetics feature input, e.g.,
#' chr20(q11.1-11.1), or providing the chr/start/end loction of the interested
#' region.
#' @param Counts The single cell expression matrix for the whole genome of the
#' sample. Rows are genes and columns are cell IDs.
#' @param normalization Specify whether the data need to be normalized. Default is TRUE.
#' @param qt_cutoff A quantile cut-off to remove genes that are almost all zeros.
#' If the cut-off is 0.99, then all the genes expressed in less than 0.01 percent
#' of cells will be eliminated for further analysis.
#'
#' @return A list with normalized and ordered gene expression for the interested
#' cytogenetics region.
#' @importFrom magrittr %>%
#' @export
#' @examples
#' res <- GetCytoLocation(cyto_feature = "chr20(q11.1-q13.1)")
#' data(SimData)
#' GetExprCountCyto(cytoloc_output = res, Counts = as.matrix(SimData), normalization = TRUE, 
#' qt_cutoff = 0.99)
#'
GetExprCountCyto <- function(cytoloc_output,
                             Counts = NULL,
                             normalization = TRUE,
                             qt_cutoff = 0.99) {

    Hg38_gtf <- NULL

    data("Hg38_gtf", envir=environment())

    if(cytoloc_output$Downstream_index == 0) {
        stop("Please re-run GetCytoLocation with a complete cyto input, e.g., chr20(q11.1-q13.1), or provide chr/start/end information.")
    } else if(cytoloc_output$Downstream_index == 1) {

        if(!normalization) {
            normalizedCounts <- Counts
        } else {
            normFactor <- colMeans(Counts)/stats::median(colMeans(Counts))
            normalizedCounts <- sweep(Counts, 2, normFactor, "/")
            rownames(normalizedCounts) <- rownames(Counts)
        }

        fCounts <- normalizedCounts[which(rownames(normalizedCounts) %in% cytoloc_output$overGeneName), ]
        tmpGL <- Hg38_gtf[which(Hg38_gtf$gene_name %in% rownames(fCounts)), ]
        GeneLocation <- tmpGL$start
        fCounts_ord <- fCounts[order(GeneLocation),]
        GeneLocation_ord <- GeneLocation[order(GeneLocation)]

        fout <- Fselect(fCounts_ord, cutoff = qt_cutoff, status = GeneLocation_ord)

        return(list(ProcessedCount = fout$count,
                    status = fout$status))
    }

}
