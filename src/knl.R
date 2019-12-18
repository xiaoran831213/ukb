## Read Plink genotype or relatedness file, produce kernels and save
## them in RDS format
## Make sure the PLINK BED reader is installed.

library(plinkBED)
main <- function(pfx)
{
    ## linear kernel pre-calcuated by PLINK
    rds <- paste0(pfx, ".rel")
    if(file.exists(rds))
    {
        cat("read ", rds, "\n", sep="")
        rel <- readRDS(rds)
    }
    else
    {
        cat("read ", pfx, ".???\n", sep="")
        rel <- readREL(pfx)
        saveRDS(rel, paste0(pfx, ".rel"))
    }

    ## exponential kernel
    rds <- paste0(pfx, ".exp")
    if(file.exists(rds))
        cat("skip ", rds, "\n", sep="")
    else
    {
        knl <- exp(rel)
        cat("save ", rds, "\n", sep="")
        saveRDS(knl, rds)
        rm(knl)
    }

    ## sigmoid kernel, L=2 (steepness)
    rds <- paste0(pfx, ".sgm")
    if(file.exists(rds))
        cat("skip ", rds, "\n", sep="")
    else
    {
        knl <- 1 / (1 + exp(-rel * 2))
        cat("save ", rds, "\n", sep="")
        saveRDS(knl, rds)
        rm(knl, rel)
    }

    ## IBS kernel pre-calculated by PLINK
    rds <- paste0(pfx, ".ibs")
    if(file.exists(rds))
    {
        cat("read ", rds, "\n", sep="")
        ibs <- readRDS(rds)
    }
    else
    {
        cat("readIBS ", pfx, "\n", sep="")
        ibs <- readIBS(pfx)
        saveRDS(ibs, rds)
    }

    ## Laplacian kernel
    rds <- paste0(pfx, ".lap")
    if(file.exists(rds))
        cat("skip ", rds, "\n", sep="")
    else
    {
        knl <- exp(ibs - 1)
        paste0("save ", pfx, ".lap")
        saveRDS(knl, rds)
        rm(ibs, knl)
    }
    gc()

    gmx <- NULL
    ## Gaussian kernel
    rds <- paste0(pfx, ".gsu")
    if(file.exists(rds))
        cat("skip ", rds, "\n", sep="")
    else
    {
        ## read genotype
        if(is.null(gmx))
            gmx <- readBED(pfx, 1)
        knl <- gsu(gmx)
        cat("save ", rds, "\n", sep="")
        saveRDS(knl, rds)
        rm(knl)
    }

    ## Heterozogosity kernel
    rds <- paste0(pfx, ".het")
    if(file.exists(rds))
        cat("skip ", rds, "\n", sep="")
    else
    {
        ## read genotype
        if(is.null(gmx))
            gmx <- readBED(pfx, 1)
        knl <- add(gmx)
        cat("save ", rds, "\n", sep="")
        saveRDS(knl, rds)
        rm(knl)
    }

    ## Additive kernel
    rds <- paste0(pfx, ".add")
    if(file.exists(rds))
        cat("skip ", rds, "\n", sep="")
    else
    {
        ## read genotype
        if(is.null(gmx))
            gmx <- readBED(pfx, 1)
        knl <- add(gmx)
        cat("save ", rds, "\n", sep="")
        saveRDS(knl, rds)
        rm(knl)
    }
    
    ## Dominence kernel
    rds <- paste0(pfx, ".dom")
    if(file.exists(rds))
        cat("skip ", rds, "\n", sep="")
    else
    {
        ## read genotype
        if(is.null(gmx))
            gmx <- readBED(pfx, 1)
        knl <- dom(gmx)
        cat("save ", rds, "\n", sep="")
        saveRDS(knl, rds)
        rm(knl)
    }
    
    invisible(NULL)
}

## fast (squared) Standardized Euclidean distance
gsu <- function(x, v=1)
{
    a <- is.na(x)
    nn <- tcrossprod(1 - a)              # pairwise non-NA
    x[a] <- 0                            # set NA to 0
        
    ## squared Euclidian distance
    x2 <- rowSums(x^2)
    ## standardized by the number of non-na SNPs
    k <- (outer(x2, x2, `+`) - 2 * tcrossprod(x)) / nn

    ## Gaussian kernel
    k <- exp(-k / (2 * v))
    k
}

add <- function(x)
{
    ## count pairwise complete SNPs
    a <- is.na(x)
    nn <- tcrossprod(1 - a)             # pairwise non-NA
    x[a] <- 0                           # set NA to 0

    ## product kernel
    k <- tcrossprod(x) / nn
    k
}

het <- function(x)
{
    ## recode to heterozygosity
    x <- x == 1                        # 1 -> 1, 0 and 2 -> 0
    ## count pairwise complete SNPs
    a <- is.na(x)
    nn <- tcrossprod(1 - a)             # pairwise non-NA
    x[a] <- 0                           # set NA to 0

    ## product kernel
    k <- tcrossprod(x) / nn
    k
}


dom <- function(x)
{
    ## recode to dominence
    x <- x > 1                          # 0 -> 0, 1 and 2 -> 1

    ## count pairwise complete SNPs
    a <- is.na(x)
    nn <- tcrossprod(1 - a)             # pairwise non-NA
    x[a] <- 0                           # set NA to 0

    ## product kernel
    k <- tcrossprod(x) / nn
    k
}

enc.knl <- function(kns, enc=0, ...)
{
    nms <- names(kns)
    if(enc >= 1)
    {
        for(i in seq_along(nms))
        {
            r <- list(exp(kns[[i]]))
            names(r) <- paste0("EX", i)
            kns <- c(kns, r)
        }
    }
    if(enc >= 2)
    {
        sgm <- function(x) 1 / (1 + exp(-x))
        for(i in seq_along(nms))
        {
            r <- list(sgm(kns[[i]]))
            names(r) <- paste0("SG", i)
            kns <- c(kns, r)
        }
    }
    kns
}
