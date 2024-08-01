#' Infer cells that are locally aneuploid using partCNVH
#'
#' This function uses EM algorithm to cluster the cells with a Poisson Mixture model in the first step.
#' With the results, it applies hidden markov model to improve feature selection.
#' After that, another round of EM algorithm is applied to obtain the final cell status.
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
#' @param navg Number of genes used for rolling average.
#'
#' @importFrom data.table frollmean
#' @importFrom depmixS4 depmix fit posterior
#'  
#' @return A vector with the cell status inferred by the method, 1 is aneuploid and 0 is diploid.
#' @export
#' @examples
#' cytoloc <- GetCytoLocation(cyto_feature = "chr20(q11.1-q13.1)")
#' data(SimData)
#' exprout <- GetExprCountCyto(cytoloc_output = cytoloc, Counts = as.matrix(SimData), normalization = TRUE, qt_cutoff = 0.99)
#' status <- partCNVH(int_counts = exprout$ProcessedCount, cyto_type = "del", cyto_p = 0.2, navg = 50)
#'
partCNVH <- function(int_counts,
                     cyto_type,
                     cyto_p,
                     tau = 0.1,
                     maxniter = 1000,
                     navg = 50) {

    if (nrow(int_counts) <= 20) {
        warning("The number of genes is only ", nrow(int_counts), ". It maybe better to use partCNV only.")
    }

    EMlabel <- partCNV(int_counts = int_counts,
                     cyto_type = cyto_type,
                     cyto_p = cyto_p,
                     tau = tau,
                     maxniter = maxniter)

    if(cyto_type == "del") {
        meanratio <- rowMeans(int_counts[, EMlabel == 0])/rowMeans(int_counts[, EMlabel == 1])
        meanratio2 <- frollmean(meanratio, n = navg, na.rm = TRUE, align = "center")
        mysumdata <- data.frame(rowmean = meanratio2)
        initStatus <- rep(1, length(mysumdata$rowmean))
        initStatus[mysumdata$rowmean > stats::median(mysumdata$rowmean)] <- 2
        mod <- depmix(rowmean ~ 1, data = mysumdata, nstates = 2, initdata = initStatus, trstart = c(0.9,0.1,0.1,0.9)) # use gaussian() for normally distributed data
        fit.mod <- depmixS4::fit(mod)
        est.states <- posterior(fit.mod)

        if(mean(mysumdata$rowmean[est.states$state == 2], na.rm = TRUE) < mean(mysumdata$rowmean[est.states$state == 1], na.rm = TRUE)) {
            myidealstate <- 1
        } else {
            myidealstate <- 2
        }
    } else if(cyto_type == "amp") {
        meanratio <- rowMeans(int_counts[, EMlabel == 1])/rowMeans(int_counts[, EMlabel == 0])
        meanratio2 <- frollmean(meanratio, n = navg, na.rm = TRUE, align = "center")
        mysumdata <- data.frame(rowmean = meanratio2)
        initStatus <- rep(1, length(mysumdata$rowmean))
        initStatus[mysumdata$rowmean > stats::median(mysumdata$rowmean)] <- 2
        mod <- depmix(rowmean ~ 1, data = mysumdata, nstates = 2, initdata = initStatus, trstart = c(0.9,0.1,0.1,0.9)) # use gaussian() for normally distributed data
        fit.mod <- depmixS4::fit(mod)
        est.states <- posterior(fit.mod)

        if(mean(mysumdata$rowmean[est.states$state == 2], na.rm = TRUE) < mean(mysumdata$rowmean[est.states$state == 1], na.rm = TRUE)) {
            myidealstate <- 1
        } else {
            myidealstate <- 2
        }
    }

    EMHMMlabel <- partCNV(int_counts[est.states$state == myidealstate, ],
                        cyto_type = cyto_type,
                        cyto_p = cyto_p,
                        tau = tau,
                        maxniter = maxniter)

    return(list(EMlabel = EMlabel,
                EMHMMlabel = EMHMMlabel))
}


