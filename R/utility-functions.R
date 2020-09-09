#' Calculate R0 implied from the model
#'
#' @param tbar Average duration of the communicable window.
#' @param lambda Poisson rate
#' @param mu Expectation of the logarithmic distribution giving the
#'     number of individuals infected per event
#' @param q Probability of node having a communicable window of the
#'     second kind.
#' @param mbar Average duration for communicable windows of the second kind.
#'     
#'
#' @return The value of \eqn{R_0}.
#' @export
calc_R0 <- function(tbar, lambda, mu, q = 0, mbar = NA) {
    if (q == 0) {
        a <- 1
    } else {
        a <- 1 - q * (1 - mbar / tbar)
    }
    return(a * mu * lambda * tbar)
}


#' Find Gamma rate parameter
#'
#' Given the average duration of the communicable window, and an upper
#' bound corresponding to the \eqn{p}-th quantile (default \code{p =
#' 0.95}), find the corresponding Gamma rate parameter.
#'
#' @param tbar A scalar greater than 0. Average duration of the
#'     communicable window.
#' @param p_val A scalar greater than 0. Upper quantile of the
#'     communicable window.
#' @param p A scalar greater than 0 and smaller than 1. Quantile level
#'     (default 0.95)
#'
#' @return A shape parameter higher than zero
#' @export
find_kappa <- function(tbar, p_val, p = .95) {
    uniroot(function(kappa) {p_val - qgamma(p, tbar * kappa, kappa)},
            interval = c(.01, 200))[["root"]]
}


#' Compute average epidemic path
#'
#' @param lambda A non-negative scalar. Rate parameter for the Poisson
#'     distribution that governs the number of infectious events per
#'     parent node.
#' @param mu A scalar greater than 1. Average number of infections
#'     resulting from each infectious interaction.
#' @param A The shape parameter of the gamma life time distribution (default: 5.5).
#' @param B The rate parameter of the gamma life time distribution (default: 0.85).
#' @param tmax A scalar greater than 0. Calculate path for times between 0 and this value.
#' @param nstep Grid size.
#'
#' @return A \code{data.frame} with columns \code{time},
#'     \code{n_infected} (cumulative number infected -- strictly
#'     increasing), \code{n_infectious} (number infectious at time \code{time}).
#' 
#' @export
calc_avg_path <- function(lambda = .11,
                          mu = 1.5,
                          A = 5.5,
                          B = 0.85,
                          tmax = 100,
                          nstep = 10000) {
    
    get_p <- Vectorize(find_p)
    p <- get_p(mu)

    dmu_fun <- dMu(A=A, B=B, Lambda = lambda, P = p)

    # Calculate the theoretical number of infected / infectious
    # individuals vs time
    eval_tibble <- left_join(renewal_function(dmu_fun,
                                              # n_infected
                                              G = 1,
                                              Time_limit = tmax,
                                              nstep = nstep) %>%
                             rename(n_infected = solution),
                             #####
                             renewal_function(dmu_fun,
                                              # n_infectious
                                              G = G(A = A, B = B),
                                              Time_limit = tmax,
                                              nstep = nstep) %>%
                             rename(n_infectious = solution),
                             #####
                             by = "time")

    return(eval_tibble)
}
