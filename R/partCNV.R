#' Infer cells that are locally aneuploid using partCNV
#'
#' This function uses EM algorithm to cluster the cells with a Poisson Mixture model.
#' Cells will be grouped into two groups, locally aneuploid (status = 1) and
#' diploid (status = 0).
#'
#' @param int_counts Normalized gene expression counts for the genes in the interested
#' region, e.g., the ProcessedCount variable from the output of GetExprCountCyto().
#' @param cyto_type The type of the cytogenetics alteration. It can only be "del" or "amp"
#' @param cyto_p The percentage of cells with the cytogenetics alteration, e.g., 0.2.
#' @param tau The variance of the prior information. Default is 0.1. If you have less confidence,
#' specify a larger tau, e.g., 10.
#' @param maxniter The maximum number of iterations of the EM algorithm.
#' 
#' @importFrom Seurat CreateSeuratObject NormalizeData ScaleData RunPCA FindNeighbors Idents FindClusters
#'
#' @return A vector with the cell status inferred by the method, 1 is aneuploid and 0 is diploid.
#' @export
#' @examples
#' ### example 1
#' cytoloc <- GetCytoLocation(cyto_feature = "chr20(q11.1-q13.1)")
#' data(SimData)
#' exprout <- GetExprCountCyto(cytoloc_output = cytoloc, Counts = as.matrix(SimData), normalization = TRUE, qt_cutoff = 0.99)
#' status <- partCNV(int_counts = exprout$ProcessedCount, cyto_type = "del", cyto_p = 0.2)
#'
partCNV <- function(int_counts,
                  cyto_type,
                  cyto_p,
                  tau = 0.1,
                  maxniter = 1000) {

    ngene <- nrow(int_counts)
    ncell <- ncol(int_counts)

    tmp <- int_counts
    int_counts <- round(as.matrix(int_counts))

    ### initialization with PCA+louvian
    input <- int_counts
    if(is.null(rownames(int_counts))) {
        rownames(input) <- paste0("gene", seq_len(nrow(input)))
    }
    if(is.null(colnames(int_counts))) {
        colnames(input) <- paste0("cell", seq_len(ncol(input)))
    }
    sim <- CreateSeuratObject(counts = input, project = "sim", min.cells = 0, min.features = 0)
    sim <- NormalizeData(sim)
    sim <- ScaleData(sim, features = rownames(input))
    sim <- RunPCA(sim, features = rownames(input))
    sim <- FindNeighbors(sim, dims = seq_len(10))
    reso_vec <- seq(0.001, 0.5, by = 0.01)
    K <- 1
    j <- 1
    while(K == 1) {
        sim <- FindClusters(sim, resolution = reso_vec[j], algorithm = 1)
        j <- j + 1
        K <- length(unique(Idents(sim)))
    }
    cellZ <- Idents(sim)

    mres <- tapply(colMeans(int_counts), cellZ, mean)
    qi <- sum(cellZ == names(which.min(mres)))/length(cellZ)
    choice1 <- cutfunc(abs(stats::rnorm(length(cellZ), mean = qi, sd = 0.1)))
    choice0 <- cutfunc(abs(stats::rnorm(length(cellZ), mean = 1-qi, sd = 0.1)))
    pi <- ifelse(cellZ == names(which.min(mres)), choice1, choice0)
    mu_g1 <- rowSums(int_counts * getweightmat(pi, ngene, ncell))/sum(pi)
    mu_g0 <- rowSums(int_counts * getweightmat(1-pi, ngene, ncell))/sum(1-pi)
    mu_g1[mu_g1<0] <- 0.01
    mu_g0[mu_g0<0] <- 0.01

    ### q0 is always the probability of group with lower gene expressions
    if (tolower(cyto_type) == "del") {
        q0 <- cyto_p
    } else if (tolower(cyto_type) == "amp") {
        q0 <- 1 - cyto_p
    } else {
        stop("cyto_type must be del or amp!")
    }

    pi_old <- rep(0, length(cellZ))
    diff_p <- 10
    niter <- 1
    Ci <- rep(0, length(cellZ))
    diff_prp <- 10
    ### update the parameters
    while(diff_p > 10^-5 & diff_prp > 10^-5 & niter < maxniter) {
        pi_old <- pi
        Cold <- Ci

        # ## mu_g1 group should be the group with smaller mean
        if(mean(mu_g1) > mean(mu_g0)) {
            tmp <- mu_g0
            mu_g0 <- mu_g1
            mu_g1 <- tmp
            qi <- 1-qi
        }

        # ## update qi
        Ci <- ifelse(pi>=0.5, 1, 0)
        Ctotal <- sum(Ci)
        qi <- Ctotal/ncell
        fnx <- function(x) {
            Ctotal/x - (ncell - Ctotal)/(1-x) - ncell*(x - q0)/tau^2
        }
        rrout <- stats::uniroot(fnx, interval = c(0.001,0.99))
        qi <- rrout$root

        ## update qi
        logP1 <- stats::dpois(as.matrix(int_counts), lambda = mu_g1, log=TRUE)
        logP1sum <- colSums(logP1)
        logP0 <- stats::dpois(as.matrix(int_counts), lambda = mu_g0, log=TRUE)
        logP0sum <- colSums(logP0)
        pi <- qi/(qi + exp(logP0sum - logP1sum)*(1-qi))

        ## update mu
        mu_g1 <- rowSums(int_counts * getweightmat(pi, ngene, ncell))/sum(pi)
        mu_g0 <- rowSums(int_counts * getweightmat(1-pi, ngene, ncell))/sum(1-pi)

        Ci <- ifelse(pi>=0.5, 1, 0)
        Ctotal <- sum(Ci)

        if(Ctotal == 0 | Ctotal == ncol(int_counts)) {
            cutoff <- stats::quantile(pi, qi)
            Ci <- ifelse(pi >= cutoff, 1, 0)
            Ctotal <- sum(Ci)
        }

        niter <- niter + 1
        diff_p <- abs(sum(pi - pi_old))
        diff_prp <- sum(abs(Ci - Cold))/ncell
        message("qi = ", qi)
        message("diff_p = ", diff_p)
        message(Ctotal/ncell)
    }

    if (tolower(cyto_type) == "amp") {
        Ci <- 1 - Ci
        message(sum(Ci)/ncell)
    }

    return(Ci)
}
