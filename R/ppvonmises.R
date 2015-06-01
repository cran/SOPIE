ppvonmises <-
function (q, mu, kappa, from = NULL, tol = 1e-20)
{
    q <- q%%(2 * pi)
    n <- length(q)
    mu <- mu%%(2 * pi)
    pvm.mu0 <- function(q, kappa, tol) {
        flag <- TRUE
        p <- 1
        sum <- 0
        while (flag) {
            term <- (besselI(x = kappa, nu = p, expon.scaled = FALSE) * 
                sin(p * q))/p
            sum <- sum + term
            p <- p + 1
            if (abs(term) < tol) 
                flag <- FALSE
        }
        return(q/(2 * pi) + sum/(pi * besselI(x = kappa, nu = 0, 
            expon.scaled = FALSE)))
    }
    result <- rep(NA, n)
    if (mu == 0) {
        for (i in 1:n) {
            result[i] <- pvm.mu0(q[i], kappa, tol)
        }
    }
    else {
        for (i in 1:n) {
            if (q[i] <= mu) {
                upper <- (q[i] - mu)%%(2 * pi)
                if (upper == 0) 
                  upper <- 2 * pi
                lower <- (-mu)%%(2 * pi)
                result[i] <- pvm.mu0(upper, kappa, tol) - pvm.mu0(lower, 
                  kappa, tol)
            }
            else {
                upper <- q[i] - mu
                lower <- mu%%(2 * pi)
                result[i] <- pvm.mu0(upper, kappa, tol) + pvm.mu0(lower, 
                  kappa, tol)
            }
        }
    }
    return(result)
}
