Fselect <- function(fCounts, cutoff = 0.95, status = NULL) {
    expRatio <- rowSums(fCounts>0)/ncol(fCounts)
    fCounts2 <- fCounts[expRatio>1 - cutoff, ]
    if(length(status) == 0) {
        return(list(count = fCounts2))
    } else {
        status2 <- status[expRatio>1 - cutoff]
        return(list(count = fCounts2,
                    status = status2))
    }
}
