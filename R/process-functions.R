#' Leaving this function here for now, but will probably remove soon.
#' \code{\link{treelist_to_paths}} is much faster.
#' @export
make_path <- function(tree_raw, tmax = Inf) {
    tree <-
        # bind all layers
        tree_raw %>%
        bind_rows %>%
        select(t_infect, t_comm) %>%
        mutate(t_heal = t_infect + t_comm) %>%
        select(-t_comm) %>%
        pivot_longer(everything(), names_to = "label", values_to = "time") %>%
        mutate(delta = ifelse(label == "t_infect", 1, -1)) %>%
        select(-label) %>%
        arrange(time) %>%
        mutate(date = floor(time)) %>%
        filter(date < tmax)

    # For all paths that extinguished before tmax, pad with zeros
    # until tmax
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



#' Compute paths from a list of trees
#'
#' This function computes paths from all trees at once; it is much
#' faster than iterating over trees (e.g. as one would do using
#' \code{purrr::map}).
#'
#' @param treelist A list of lists. Each sublist corresponds to a
#'     single tree (i.e. a simulated path). Each tree is itself a list
#'     of tibbles, with one tibble for each depth level.
#' @export
treelist_to_paths <- function(treelist,
                              tmax = Inf) {
    tibble(id_sim = seq_along(treelist),
           tree = treelist) %>%
        # first, unnest the sublists
        unnest(tree) %>%
        # second, unnest the tibbles
        unnest(tree) %>%
        # compute a timestamp for the end of the communicable period
        mutate(t_heal = t_infect + t_comm) %>%
        # all we need from hereon are id_sim (for grouping), and t_infect,
        # t_heal
        select(id_sim, t_infect, t_heal) %>%
        # below we build n_infectious as a cumulative sum of births (+1)
        # and deaths (-1), and n_infected as the cumulative sum of births
        # only.
        pivot_longer(cols = c(t_infect, t_heal),
                     names_to = "label",
                     values_to = "time") %>%
        mutate(date = floor(time),
               delta_pos = (label == "t_infect"),
               delta = ifelse(delta_pos, 1, -1)) %>%
        select(-time, -label) %>%
        group_by(id_sim, date) %>%
        summarize_at(vars(delta, delta_pos), sum) %>%
        mutate(n_infected = cumsum(delta_pos),
               n_infectious = cumsum(delta)) %>%
        # cleanup
        ungroup %>%
        select(-starts_with("delta")) %>%
        filter(date <= tmax)
}


#' Make all paths end at a fixed time
#'
#' @param df_paths A \code{data.frame}. Typically, the output of
#'     \code{\link{treelist_to_paths}}. The following columns are
#'     required: \code{date}, \code{id_sim}.
#' @param tmax A non-negative scalar.
#'
#' @return A \code{tbl} with the same format as \code{df_paths}, but
#'     with the last point for each path fixed at \code{tmax}.
#'
#' @export
equalize_paths <- function(df_paths, tmax) {
    # Truncate points beyond tmax
    df_paths <- df_paths %>% filter(date <= tmax)
    # Extend paths that need it
    tails <-
        df_paths %>%
        arrange(date) %>%
        group_by(id_sim) %>%
        summarize_all(~ tail(., 1)) %>%
        filter(date < tmax) %>%
        mutate(date = tmax)
    # Add extensions to truncated paths
    df_paths_equal <-
        bind_rows(df_paths, tails) %>%
        arrange (id_sim, date)
    return(df_paths_equal)
}
