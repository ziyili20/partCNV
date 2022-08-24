getweightmat <- function(pi, ngene, ncell) {
    return(matrix(rep(pi, each = ngene), ngene, ncell))
}
