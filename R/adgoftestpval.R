adgoftestpval <-
function (data, n) 
{
    if (data < 2) 
        data <- exp(-1.2337141/data)/sqrt(data) * (2.00012 + (0.247105 - 
            (0.0649821 - (0.0347962 - (0.011672 - 0.00168691 * 
                data) * data) * data) * data) * data)
    else data <- exp(-exp(1.0776 - (2.30695 - (0.43424 - (0.082433 - 
        (0.008056 - 0.0003146 * data) * data) * data) * data) * data))
    if (data > 0.8) 
        return(data + (-130.2137 + (745.2337 - (1705.091 - (1950.646 - 
            (1116.36 - 255.7844 * data) * data) * data) * data) * data)/n)
    z <- 0.01265 + 0.1757/n
    if (data < z) {
        v <- data/z
        v <- sqrt(v) * (1 - v) * (49 * v - 102)
        return(data + v * (0.0037/(n * n) + 0.00078/n + 6e-05)/n)
    }
    v <- (data - z)/(0.8 - z)
    v <- -0.00022633 + (6.54034 - (14.6538 - (14.458 - (8.259 - 
        1.91864 * v) * v) * v) * v) * v
    data + v * (0.04213 + 0.01365/n)/n
}
