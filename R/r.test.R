r.test <-
function (data) 
{
    n <- length(data)
    ss <- sum(sin(data))
    cc <- sum(cos(data))
    rbar <- (sqrt(ss^2 + cc^2))/n
    z <- (n * rbar^2)
    p.value <- exp(-z)
    if (n < 50) 
        temp <- (1 + (2 * z - z^2)/(4 * n) - (24 * z - 132 * 
            z^2 + 76 * z^3 - 9 * z^4)/(288 * n^2))
    else temp <- 1
    result <- list(r.bar = rbar, p.value = p.value * temp)
    return (result)
}
