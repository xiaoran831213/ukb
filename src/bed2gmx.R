library(BEDMatrix)
gmx <- function(pfx, imp=TRUE, sav=NULL)
{
    ## get file names
    bed <- paste(pfx, ".bed", sep="")
    fam <- paste(pfx, ".fam", sep="")
    bim <- paste(pfx, ".bim", sep="")

    ## individuals
    fam <- read.table(
        fam, col.names=c('FID', 'IID', 'PID', 'MID', 'SEX', 'PHE'), as.is=TRUE)
    iid <- fam$IID
    
    ## genomic variants
    bim <- read.table(
        bim, col.names=c('CHR', 'SNP', 'NULL', 'POS', 'REF', 'ALT'),
        colClasses=c('integer', 'character', 'NULL', 'integer', 'character', 'character'),
        as.is=TRUE)
    snp <- bim$SNP

    ## genomic matrix
    bed <- BEDMatrix(bed, n=nrow(fam), p=nrow(bim))
    gmx <- bed[,]
    dimnames(gmx) <- list(iid=iid, snp=snp)

    ## naive imputation
    if(imp)
    {
        gmx <- apply(gmx, 2L, function(g)
        {
            frq <- table(g)
            lvl <- as.integer(names(frq))
            mrk <- is.na(g)
            g[mrk] <- sample(lvl, sum(mrk), TRUE, frq)
            g
        })
    }

    ## save and return
    ret <- list(gmx=gmx, map=bim, fam=fam)
    if(is.null(sav))
        sav <- paste0(pfx, '.rds')
    if(!is.na(sav))
        saveRDS(ret, sav)
    ret
}
