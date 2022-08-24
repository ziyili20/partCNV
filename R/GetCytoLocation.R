#' Get exact location of the interested cytogenetics feature
#'
#' This function helps you identify the location of the cytogenetics feature.
#' For example, if the region of interest is chr20(q11.1-q13.1), this function
#' greps the start and end location of this
#' region. Additionally, you can just put in "chr20", and it provides you all
#' the available cytogenetics locations on chromosome 20.
#' It also report the number of genes within the region. If the number of genes
#' is too few, we recommend to include neighboring regions to provide more
#' stable results.
#'
#' @param cyto_feature the cytogenetics location you are interested. It can be
#' of two format: chr20(q11.1-q13.1) or chr20. For the first format, the start
#' and end regions need to be separated by "-". If you are interested in one
#' region for example, chr20(q11.1), put it as chr20(q11.1-11.1). For the second
#' format, all the available regions will be printed for selection.
#' @param chr chromosome location of the interested region. This is only used
#' when cyto_feature is null.
#' @param start starting location of the interested region. This is only used
#' when cyto_feature is null.
#' @param end ending location of the interested region. This is only used
#' when cyto_feature is null.
#'
#' @return If the first format of cyto_feature is provided, the starting and
#' ending location as well as the number of genes overlapped with be provided.
#' If the second format of cyto_feature is provided, all the cytogenetics locations
#' will be displayed for review. If the region location (chr, start, end) is
#' provided, the number of genes overlapped will be the output.
#' @export
#' @examples
#' ### example 1
#' GetCytoLocation(cyto_feature = "chr20(q11.1-q13.1)")
#' ### example 2
#' GetCytoLocation(cyto_feature = "chr20")
#' ### example 3
#' GetCytoLocation(chr = "chr20", start = 25600000, end = 49800000)
#'
GetCytoLocation <- function(cyto_feature = NULL,
                            chr = NULL,
                            start = NULL,
                            end = NULL) {

    Hg38_gtf <- NULL
    cytoBand_Hg38 <- NULL

    data("cytoBand_Hg38", envir=environment())
    data("Hg38_gtf", envir=environment())

    GeneNum_index <- 0
    Downstream_index <- 0

    if(is.null(cyto_feature) & is.null(chr)) {
        stop("If cyto_feature is null, please provide the region location through chr, start, and end!")
    } else if (!is.null(cyto_feature)) {

        if(length(grep("\\(", cyto_feature)) == 0) {

            if(length(grep(cyto_feature, cytoBand_Hg38$chr)) == 0) {
                stop("Please provide a valid cyto_feature input!")
            } else {
                message(paste0("Print out all cytogenetics features on ", cyto_feature, ":"))
                tmpout <- cytoBand_Hg38[grep(cyto_feature, cytoBand_Hg38$chr), ]
                rownames(tmpout) <- NULL
                print(tmpout)
                return(list(cytogeneticsInfo = tmpout,
                            Downstream_index = Downstream_index))
            }

        } else if (length(grep("\\(", cyto_feature)) > 0) {

            tmpout2 <- unlist(strsplit(cyto_feature, split = "\\("))
            chr_location <- tmpout2[1]
            if(length(grep(chr_location, cytoBand_Hg38$chr)) == 0) {
                stop("Please provide a valid chromosome information in cyto_feature input!")
            } else {
                filtered_cyto <- cytoBand_Hg38[grep(chr_location, cytoBand_Hg38$chr), ]
            }

            tmpout3 <- gsub(")", "", tmpout2[2])
            tmpout4 <- unlist(strsplit(tmpout3, split = "-"))

            idx1 <- grep(tmpout4[1], filtered_cyto$cytoband)
            idx2 <- grep(tmpout4[2], filtered_cyto$cytoband)

            if(length(idx1) == 0 | length(idx2) == 0) {
                stop("The cytoband information couldn't be recognized. Please try using the chr info only.")
            } else {
                cyto_start <- min(filtered_cyto$start[idx1])
                cyto_end <- max(filtered_cyto$end[idx2])
            }
            GeneNum_index <- 1

        }
    } else if(is.null(cyto_feature) & !is.null(chr) & !is.null(start) & !is.null(end)) {

        chr_location <- chr
        cyto_start <- start
        cyto_end <- end
        GeneNum_index <- 1

    } else {
        stop("The chr/start/end information couldn't be recognized.")
    }

    if(GeneNum_index == 1) {

        gene.gr <- GenomicRanges::GRanges(seqnames = Hg38_gtf$seqnames,
                                          ranges = IRanges::IRanges(start = as.numeric(Hg38_gtf$start),
                                                                    end = as.numeric(Hg38_gtf$end)))
        cyto.gr <- GenomicRanges::GRanges(seqnames = chr_location,
                                          ranges = IRanges::IRanges(start = cyto_start,
                                                                    end = cyto_end))
        FOout <- as.data.frame(GenomicRanges::findOverlaps(cyto.gr, gene.gr))
        overGene <- unique(Hg38_gtf$gene_id[unique(FOout$subjectHits)])
        overGeneName <- unique(Hg38_gtf$gene_name[unique(FOout$subjectHits)])
        message(paste0("Interested region: ", chr_location, ":", cyto_start, "-", cyto_end, "."))
        message(paste0("A total of ", length(overGeneName), " genes are located in this region."))
        Downstream_index <- 1

        return(list(chr_location = chr_location,
                    cyto_start = cyto_start,
                    cyto_end = cyto_end,
                    overGeneID = overGene,
                    overGeneName = overGeneName,
                    Downstream_index = Downstream_index))
    }
}


