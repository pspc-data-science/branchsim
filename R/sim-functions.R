#' Initialize the tree
#'
#' Creates the first layer of a tree.
#'
#' @param n_start An integer greater or equal to 1. Number of nodes in
#'     the initial layer.
#' @param tbar A non-negative scalar. Average duration of the
#'     communicable windows of the new nodes. Together with
#'     \code{kappa}, it determines a shape parameter \code{tbar *
#'     kappa} for the the Gamma distribution from which we draw a
#'     communicable period for each new node.
#' @param kappa A non-negative scalar (including \code{Inf}). Rate
#'     parameter for the Gamma distribution from which we draw a
#'     communicable period duration for each new node. If \code{kappa
#'     = Inf}, then the draw is deterministic, and all communicable
#'     periods are exactly equal to \code{tbar}.
#' 
#' @return A \code{tibble} with one row for each new node, and three
#'     columns: \code{id_parent}, \code{t_infect}, and \code{t_comm}.
#'
make_first_layer <- function(n_start, tbar, kappa) {
      first_layer <- tibble(id = seq_len(n_start),
                            id_parent = 0,
                            id_layer = 1,
                            t_infect = 0,
                            t_comm = rgamma(n_start, kappa * tbar, kappa))
      return(first_layer)
}


#' Builds the next generation of nodes in a tree
#'
#' Given a set of parent nodes and their properties, this function
#' simulates the next generation of infected nodes.
#'
#' @param parent_layer A \code{data.frame}. Contains the properties of
#'     each parent node. The columns are \code{id_parent} (the index
#'     of the parent's parent node), \code{t_infect} (the clock time
#'     at which that parent was infected), \code{t_comm} (the duration
#'     of the parent's communicable period).
#' @param tmax A non-negative scalar. Cutoff time for the simulated
#'     branches.
#' @param tbar A non-negative scalar. Average duration of the
#'     communicable windows of the new nodes. Together with
#'     \code{kappa}, it determines a shape parameter \code{tbar *
#'     kappa} for the the Gamma distribution from which we draw a
#'     communicable period for each new node.
#' @param p A scalar between 0 and 1. Parameter for a logarithmic
#'     distribution, giving the expected number of new infections per
#'     infection event.
#' @param lambda A non-negative scalar. Parameter for a Poisson
#'     distribution, giving the expected number of infection events
#'     per parent node.
#' @param kappa A non-negative scalar (including \code{Inf}). Rate
#'     parameter for the Gamma distribution from which we draw a
#'     communicable period duration for each new node. If \code{kappa
#'     = Inf}, then the draw is deterministic, and all communicable
#'     periods are exactly equal to \code{tbar}.
#' @param q A scalar between 0 and 1. The probability that a parent
#'     node is intercepted (e.g. through contact tracing).
#' @param mbar A non-negative scalar. The value replacing \code{tbar}
#'     for intercepted parent nodes. This parameter has no role if
#'     \code{q = 0}.
#' @param kappaq A non-negative scalar. The value replacing
#'     \code{kappa} for intercepted parent nodes. This parameter has
#'     no role if \code{q = 0}.
#'
#' @return A \code{tibble} with one row for each new node, and three
#'     columns: \code{id_parent}, \code{t_infect}, and \code{t_comm}.
#'
#' @export
add_layer <- function(parent_layer,
                      tmax,
                      tbar,
                      p,
                      lambda,
                      kappa,
                      q,
                      mbar,
                      kappaq) {

    # temporary index variable: keep only parents infected before tmax
    p_id <- which(parent_layer$t_infect < tmax)
    n_parents <- length(p_id)
    # number of infection events for each parent
    n_events <- rpois(n_parents, lambda * parent_layer$t_comm[p_id])
    # If no events, stop and return NULL
    if (sum(n_events) == 0) {
        return(NULL)
    }
    # Compute the timestamp of each infection event (t_infect) by
    # adding the parent's t_infect, to uniformly distributed random
    # draws from each parent's t_comm (one draw per infection event).
    t_infect_parents <- rep(parent_layer$t_infect[p_id], times = n_events)
    t_comm_parents <- rep(parent_layer$t_comm[p_id], times = n_events)
    t_infect <-  t_infect_parents + t_comm_parents * runif(sum(n_events))
    # Random draws for the number of new infections per infection
    # event, and the communicable period, for each new infected node.
    # No need to keep track of parents here.
    n_infect <- extraDistr::rlgser(sum(n_events), p)
    t_comm <- rgamma(sum(n_infect), kappa * tbar, kappa)
    # When q > 0, the natural communicable period can be interrupted.
    if (q > 0) {
        n_catch <- rbinom(1, sum(n_infect), q)
        if (is.infinite(kappaq)) {
            t_stop <- rep(mbar, n_catch)
        } else {
            t_stop <- rgamma(n_catch, kappaq * mbar, kappaq)
        }
        # Pad t_stop with Inf; these Inf's won't be picked by pmin
        # below, and t_comm will default to the uninterrupted case for
        # those.
        t_stop <- c(t_stop, rep(Inf, sum(n_infect) - n_catch))
        # For caught individuals, pick the minimum between t_comm and
        # t_stop, then shuffle with `sample`
        if (sum(n_infect) > 1) {
            # shuffle multiple values with `sample`
            t_comm <- sample(pmin(t_comm, t_stop))
        } else {
            # `sample` will not work with a single value
            t_comm <- min(t_comm, t_stop)
        }
    }
    new_layer <- tibble(id = seq_len(sum(n_infect)) + max(parent_layer$id),
                        id_parent =
                            rep(parent_layer$id[p_id], times = n_events) %>%
                            rep(times = n_infect),
                        t_infect = rep(t_infect, times = n_infect),
                        t_comm)
    return(new_layer)
}

#' Simulate a single epidemic path
#'
#' For information about input parameters, see the documentation for
#' \code{\link{add_layer}}. The function starts by generating a single
#' parent and its communicable period of random duration, and starting
#' at \code{t_infect = 0}. The function then builds a tree starting at
#' the initial node by iteratively calling \code{\link{add_layer}},
#' until it receives an empty tibble, which happens either when the
#' epidemic extinguishes, or all branches have extended beyond
#' \code{t_max}.
#'
#' @return A list of tibbles, each one containing information about
#'     new nodes generated from the previous ones. Effectively this is
#'     a tree.
#'
#' @export
build_tree <- function(nstart,
                       tmax,
                       tbar,
                       p,
                       lambda,
                       kappa,
                       q,
                       mbar,
                       kappaq) {
    # This is the first parent layer: one parent only, with a
    # communicable period drawn randomly from a Gamma distribution.
    old_layer <- make_first_layer(nstart, tbar, kappa)
    # The tree starts with the initial layer
    layer_count <- 1
    tree <- list(old_layer)
    # Call add_layer iteratively, until add_layer returns NULL
    while (!is.null(old_layer)) {
        layer_count <- layer_count + 1
        new_layer <- add_layer(old_layer, tmax, tbar, p, lambda, kappa,
                               q, mbar, kappaq)
        if (!is.null(new_layer)) {
            tree[[length(tree) + 1]] <-
                new_layer %>%
                mutate(id_layer = layer_count)
            old_layer <- new_layer
        } else {
            break
        }
    }
    return(tree)
}

#' Multiple path simulations over a range of input parameters
#'
#' This function calls \code{\link{build_tree}} over a range of input
#' parameters. It computes \code{nsim} replications for each unique
#' combination of input parameters.
#'
#' @param nsim A non-negative integer. The number of replications for
#'     each set of input parameters.
#' @param ... Input parameters for \code{\link{add_layer}}. Each input
#'     parameter can receive a vector with multiple values. If that is
#'     the case, then all input parameters should have the same
#'     length, or one.
#'
#' @return A tibble with a column for each input parameter, and a
#'     list-column \code{treelist}. This is a list of lists of lists.
#'     Each sublist corresponds to a unique combination of input
#'     parameters, and contains \code{nsim} trees (i.e. simulated
#'     paths). Each individual tree is itself a list of tibbles, with
#'     one tibble for each depth level.
#'
#' @export
run_sims <- function(nsim = 10,
                     nstart = 1,
                     tmax = 100,
                     tbar = 10,
                     p = .5,
                     lambda = .11,
                     kappa = 1,
                     q = .6,
                     mbar = 5,
                     kappaq = 3,
                     keep_trees = TRUE) {
    # input arguments
    args <- tibble(nstart,
                   tmax,
                   tbar,
                   p,
                   lambda,
                   kappa,
                   q,
                   mbar,
                   kappaq)
    # small functions to run nsim replications using purrr::rerun
    run_reps <- function(...) {
        reps <- rerun(nsim, build_tree(...))
        return(reps)
    }
    # run nsim replications for each row of input args
    treelist <- pmap(args, run_reps)
    # put results in a tbl, along with sim parameters
    sim_data <- bind_cols(args, tibble(treelist))
    return(sim_data)
}
