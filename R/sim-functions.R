
add_layer <- function(old_layer,
                      tmax,
                      tbar,
                      p,
                      lambda,
                      kappa,
                      q,
                      mbar,
                      kappaq) {

    p_id <- which(old_layer$t_infect < tmax)
    n_parents <- length(p_id)
    n_events <- rpois(n_parents, lambda * old_layer$t_comm[p_id])
    t_infect_parents <- rep(old_layer$t_infect[p_id], times = n_events)
    t_comm_parents <- rep(old_layer$t_comm[p_id], times = n_events)
    t_infect <-  t_infect_parents + t_comm_parents * runif(sum(n_events))
    n_infect <- extraDistr::rlgser(sum(n_events), p)
    t_comm <- rgamma(sum(n_infect), kappa * tbar, kappa)
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
        # t_stop, then shuffle
        t_comm <- sample(pmin(t_comm, t_stop))
    }
    tibble(t_comm, t_infect = rep(t_infect, times = n_infect))
}

# @export
run_model <- function(tmax,
                      tbar,
                      p,
                      lambda,
                      kappa,
                      q,
                      mbar,
                      kappaq) {
    old_layer <- tibble(t_infect = 0, t_comm = rgamma(1, kappa * tbar, kappa))
    layers <- list(old_layer)
    while (nrow(old_layer) > 0) {
        new_layer <- add_layer(old_layer, tmax, tbar, p, lambda, kappa,
                               q, mbar, kappaq)
        layers[[length(layers) + 1]] <- new_layer
        old_layer <- new_layer
    }
    return(layers)
}


# @export
process_tree <- function(tree_raw, tmax = Inf) {
    tree <-
        tree_raw %>%
        bind_rows %>%
        mutate(t_heal = t_infect + t_comm) %>%
        select(-t_comm) %>%
        pivot_longer(everything(), names_to = "label", values_to = "time") %>%
        mutate(delta = ifelse(label == "t_infect", 1, -1)) %>%
        select(-label) %>%
        arrange(time) %>%
        mutate(date = floor(time)) %>%
        filter(date < tmax)

    if (!is.infinite(tmax)) {
        all_dates <- seq(0, round(tmax) - 1, by = 1)
        extra_rows <- tibble(date = setdiff(all_dates, tree[["date"]]),
                             delta = 0)
        tree <-
            bind_rows(tree, extra_rows) %>%
            arrange(date)
    }

    tree <-
        tree %>%
        group_by(date) %>%
        summarize(n_infected = sum(delta == 1),
                  delta = sum(delta)) %>%
        mutate(n_infected = cumsum(n_infected),
               n_infectious = cumsum(delta)) %>%
        select(-delta)

    return(tree)
}


# @export
process_sims <- function(tree_list,
                         tmax = Inf,
                         keep_trees = TRUE) {
    sim_data <- tibble(sim = 1:length(tree_list),
                       trees = tree_list,
                       paths = map(tree_list, ~ process_tree(., tmax = tmax)))
    if (!keep_trees) {
        sim_data <- sim_data %>% select(-trees)
    }
    return(sim_data)
}


# @export
run_sims <- function(nsim = 10,
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
    args <- tibble(tmax = tmax,
                   tbar = tbar,
                   p = p,
                   lambda = lambda,
                   kappa = kappa,
                   q = q,
                   mbar = mbar,
                   kappaq = kappaq)
    # small functions to run nsim replications using purrr::rerun
    run_reps <- function(...) {
        reps <-
            rerun(nsim, run_model(...)) %>%
            process_sims(keep_trees = keep_trees)
        return(reps)
    }
    # run nsim replications for each row of input args
    sims <- pmap(args, run_reps)
    # put results in a tbl, along with sim parameters
    sims_data <- bind_cols(args, tibble(sims = sims))
    return(sims_data)
}
