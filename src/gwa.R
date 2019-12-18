library(plinkBED)

#' linear kernel weighted by GWAS report
gwk <- function(gwa, gno, out=NULL)
{
    gwa <- read.delim(gwa)
    gmx <- readBED(gno)

    ## w <- sqrt(abs(gwa$STAT))
    w <- gwa$STAT
    x <- t(t(gmx) * w)
    x <- as.matrix(scale(x))
        
    a <- is.na(x)
    M <- tcrossprod(1 - a)              # pairwise non-NA

    ## set NA to zero
    x[a] <- 0.0
    rm(a)

    k <- tcrossprod(x) / M

    if(!is.null(out))
        saveRDS(k, out)
    invisible(k)
}
