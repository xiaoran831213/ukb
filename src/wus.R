library(CompQuadForm)

#' weight U statistics
#'
#' @param Y design matrix of one or more response variables
#' @param w weight kernel
#' @param X design matrix of covariate
#'
#' 
ust <- function(Y, w, X=NULL, acc=1e-9, ...)
{
    ## response and covariate
    Y <- as.matrix(Y)
    N <- nrow(Y)

    ## normal rank quantile standardization
    Y <- scale(apply(Y, 2, rank))

    ## regression residual matrix, R = I - X(X'X)^X'
    if(is.null(X))
        X = matrix(1, N, 1L)
    else
        X = cbind(1, as.matrix(X))

    hat <- diag(1, N, N) - tcrossprod(X %*% solve(crossprod(X)), X)
    
    ## exclude liner covariant effect on Y, leave residual of Y
    Y <- hat %*% Y;
    Y <- apply(Y, 2, function(u)
    {
        u/sqrt(sum(u^2)/(N-ncol(X)))
    })
    
    ## the response kernel is the pair wise similarity between phenotypes
    ## now restricted to product kernel
    f <- tcrossprod(Y);
    
    diag(w) <- 0; # ??
    
    ## compute u score
    u <- sum(w * f);

    ## exclude coveriant effect on the weight kernel
    ## W* = hat' W hat
    w <- tcrossprod(hat %*% w, hat);

    ## calculate p-value of u.
    pval <- tryCatch(
    {
        coef <- eigen(w, symmetric=TRUE, only.values=TRUE)$values;
        p = davies(u, coef, acc=acc)$Qq
        p
    }, warning = function(wa)
    {
        print(wa)
        p
    }, error = function(e)
    {
        print(e)
        NA
    })
    pval
}

#' centralize kernel
ckn <- function(K)
{
    K - outer(rowMeans(K), colMeans(K), '+') + mean(K)
}

#' Gussian Similarity Kernel
GUS <- function(x)
{
    ## normalize features
    ## .map.std.norm <- function(y) qnorm((rank(y)-0.5)/length(y))
    ## x <- apply(x, 2L, .map.std.norm);
    
    ## exp(- weighted gaussian distance) = weight gaussian similiarity
    ## centralize the similarity.

    ## squared euclidian distance
    ec2 <- as.matrix(dist(x, method='euclidean')^2)
    exp(-ec2 * (.5 / NCOL(x)))
}

#' Lineaer Kernel
LNR <- function(x)
{
    ## NA mask, and pairwise non-NA counts
    na.msk <- is.na(x)
    na.cnt <- tcrossprod(1 - na.msk)

    ## pairwise products / pairwise non-NA counts
    x[na.msk] <- 0.0
    tcrossprod(x) / na.cnt
}

#' Identity by State Kernel
#'
#' @param x matrix of N row samples and P feature columns.
#' @param m maximum possible distance between two samples.
#'
#' for dosage coded genome data, the maximum distance this is 2.
IBS <- function(x, m=NULL)
{
    if(is.null(m))
        m <- max(x, na.rm=TRUE) - min(x, na.rm=TRUE)
    1 - as.matrix(dist(x, method='manhattan') / (m * NCOL(x)))
}

#' Genetic Relatedness matrix
GRM <- function(x)
{
    LNR(gsc(x))
}

#' z-score standardization
zsc <- function(x) as.matrix(scale(x))

#' g-score standardization
#'
#' Let q be the allele frequency of a genomic variant {g}, treat
#' that variant as an a binomial random variable of 2 trails, a mean
#' of 2q, and consequently a variance is 2pq. The standardized
#' variant is 
#'
#' z = (g - 2q) / sqrt(2pq)
#' 
#' @param x matrix N row samples and P column variants.
#' @return N by P matrix  of standardized variants.
gsc <- function(x, q=NULL)
{
    if(is.null(q))                      # allele frequency
        q <- colMeans(x, TRUE) / 2
    s <- sqrt((2 * q * (1 - q)))        # weights
    as.matrix(scale(x, q * 2, s))
}
