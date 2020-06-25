#' Calculate R0 implied from the model
#'
#' @param tbar Average duration of the communicable window.
#' @param lambda Poisson rate
#' @param p Logarithmic distribution parameter
#' @param q Probability of interruption
#' @param mbar Average communicable period elapsed at the time of interuption
#'
#' @return The value of \eqn{R_0}.
#' @export
calc_R0 <- function(tbar, lambda, mu, q = 0, mbar = NA) {
    if (q == 0) {
        a <- 1
    } else {
        a <- q * (mbar / tbar + 1)
    }
    return(a * mu * lambda * tbar)
}


#' Find Gamma rate parameter
#'
#' @param tbar
#' @param q
#'
#' @return A shape parameter higher than zero
#'
#' @export
find_kappa <- function(tbar, p_val, p = .95) {
    uniroot(function(kappa) {p_val - qgamma(p, tbar * kappa, kappa)},
            interval = c(.01, 200))[["root"]]
}


#' Compute average epidemic path
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
