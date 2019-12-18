## read a list of IID
.get.id <- function(fn, only.IID=TRUE)
{
    . <- read.table(fn, col.names=c('FID', 'IID'), as.is=TRUE)
    if(only.IID)
        .$IID
    else
        with(., paste(FID, IID, sep='.'))
}

## read binary matrix, guess the shape of recorded entries.
.get.bm <- function(fn, id, dg=1)
{
    M <- length(id)
    S <- file.size(fn)                  # file size
    R <- matrix(.0, M, M, dimnames=list(id, id))

    ## try: lower triangle with diagonal
    L <- M * (M + 1.0) / 2.0
    U <- S / L
    if(U == 4 || U == 8)
    {
        R[upper.tri(R, 1)] <- readBin(fn, .0, L, U)
        R[lower.tri(R, 0)] <- t(R)[lower.tri(R, 0)]
        print(data.frame(file.size=S, num.entries=L, unit.bytes=U, shape='LWD'))
        return(R)
    }

    ## try: lower triangle without diagonal
    L <- M * (M - 1.0) / 2.0              # number of entries
    U <- S / L                          # unit size
    if(U == 4 || U == 8)
    {
        R[upper.tri(R, 0)] <- readBin(fn, .0, L, U)
        R[lower.tri(R, 0)] <- t(R)[lower.tri(R, 0)]
        diag(R) <- dg                 # assigned diagnal
        print(data.frame(file.size=S, num.entries=L, unit.bytes=U, shape='LND'))
        return(R)
    }

    ## try: a squre
    L <- 1.0 * M * M
    U <- S / L
    if(U == 4 || U == 8)
    {
        R[ , ] <- readBin(fn, .0, L, U)
        print(data.frame(file.size=S, num.entries=L, unit.bytes=U, shape='SQR'))
        return(R)
    }

    ## fail: return
    stop("can not figure out the shape of recoded entries in the relatedness matrix.")
}

#' read IBS binary of plink
readIBS <- function(pfx)
{
    .get.bm(paste0(pfx, ".mibs.bin"), .get.id(paste0(pfx, ".mibs.id")))
}

readREL <- function(pfx)
{
    .get.bm(paste0(pfx, ".rel.bin"), .get.id(paste0(pfx, ".rel.id")))
}

#' read GRM binary of GCTA
#'
#' GRM is equivalent to standardized REL
readGRM <- function(pfx)
{
    .get.bm(paste0(pfx, ".grm.bin"), .get.id(paste0(pfx, ".grm.id")))
}

saveGRM <- function(pfx, grm)
{
    ## get file names
    fn.rmx <- paste0(pfx, ".grm.bin")
    fn.N <- paste0(pfx, ".grm.N.bin")
    fn.id <- paste0(pfx, ".grm.id")

    ## complete id and N
    if(is.matrix(grm))
    {
        grm <- list(rmx=grm, id=makeID(grm), N=1.0)
    }

    with(grm,
    {
        ## upper.tri of col major = lower.tri of row major
        idx <- upper.tri(diag(nrow(id)), T)
        
        ## genomic relatedness matrix
        rmx <- rmx[idx]
        writeBin(rmx, fn.rmx, 4L)

        ## genomic variant count matrix
        N <- N[idx]
        writeBin(N, fn.N, 4L)

        ## subject IDs
        write(t(id), fn.id, 2, sep='\t')
    })
}

## mean of relationship
meanGRM <- function(x)
{
    ## the first GRM
    rmx <- .0
    N <- 0L
    id <- NA

    for(e in x)
    {
        if(is.character(e))
            e <- readGRM(e)
        rmx <- rmx + e$rmx * e$N
        N <- N + e$N
        id <- e$id
    }
    rmx <- rmx/N
    
    list(rmx=rmx, N=N, id=id)
}

makeID <- function(x)
{
    N <- nrow(x)
    fid <- sprintf('F%04X', seq(N))
    iid <- sprintf('I%04X', seq(N))
    data.frame(FID=fid, IID=iid, stringsAsFactors=FALSE)
}
