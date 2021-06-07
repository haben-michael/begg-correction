## adapted from src/library/stats/R/cor.test.R

begg.test <- function(y,v,adjusted=TRUE,...) {
    
    DNAME <- paste(deparse(substitute(y)), "and", deparse(substitute(v)))
    if(length(y) != length(v))
        stop("'y' and 'v' must have the same length")
    if(!is.numeric(y)) stop("'y' must be a numeric vector")
    if(!is.numeric(v)) stop("'v' must be a numeric vector")
    OK <- complete.cases(y, v)
    y <- y[OK]
    v <- v[OK]
    n <- length(y)
    NVAL <- 0    
    if(n < 2)
        stop("not enough finite observations")
    PARAMETER <- NULL
    TIES <- (min(length(unique(y)), length(unique(v))) < n)
    if(TIES) warning("Cannot compute exact p-value with ties")
    if(adjusted) {
        method <- "Begg's test of publication bias (bias-adjusted)"
    } else { method <- "Begg's test of publication bias (standard)"}
    names(NVAL) <- "tau"
    tau.hat <- tau(y,v,...)    
    ESTIMATE <- c(tau = tau.hat)

    if(!is.finite(ESTIMATE)) {
        ESTIMATE[] <- NA
        STATISTIC <- c(T = NA)
        PVAL <- NA
    }
    else {
        ## if(exact && !TIES) {
        ## q <- round((tau.hat + 1) * n * (n - 1) / 4)
        STATISTIC <- c(T = round((tau.hat + 1) * n * (n - 1) / 4))
        ## pkendall2*pnorm(sqrt(n/(4/9-bias.hat))*abs(tau.hat) <- function(q, n) .Call(C <- pKendall, q, n)
        bias.hat <- if(adjusted) {
                        s <- sqrt(1/v)
                        s.md.hat <- mean(abs(outer(s,s,`-`)))
                        s2.hat <- mean(s^2)
                        bias.hat <- s.md.hat^2/s2.hat/pi
                    } else {
                        bias.hat <- 0
                    }
        PVAL <-2*pnorm(sqrt(n/(4/9-bias.hat))*abs(tau.hat),lower.tail=FALSE)
    }
    
    RVAL <- list(statistic = STATISTIC,
                 parameter = PARAMETER,
                 p.value = as.numeric(PVAL),
                 estimate = ESTIMATE,
                 null.value = c(tau=0),
                 alternative = 'two.sided',
                 method = method,
                 data.name = DNAME)
    class(RVAL) <- "htest"
    RVAL
}
