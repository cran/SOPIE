adgofteststat <-
function (data, distr.fun, lower ,upper, tol=0.001, ...) 
{
    DNAME <- deparse(substitute(data))
    data <- sort(data)
    if (missing(distr.fun) && (data[1] < 0 || data[length(data)] > 1)) 
        stop(paste("Data ", DNAME, " is not in the [0,1] range."))
    if (!missing(distr.fun)) {
        data <- distr.fun(data, lower-tol,upper+tol)
        DNAME <- paste(DNAME, " and ", deparse(substitute(distr.fun)))
    }
    tmp <- data * (1 - rev(data))
    tmp <- (2 * seq(data) - 1) * log(tmp)
    STATISTIC <- -mean(tmp) - length(data)   
}
