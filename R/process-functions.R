#' Compute paths from a list of trees
#'
#' This function computes paths from all trees at once; it is much
#' faster than iterating over trees (e.g. as one would do using
#' \code{purrr::map}).
#'
#' @param treelist A list of lists. Each sublist corresponds to a
#'     single tree (i.e. a simulated path). Each tree is itself a list
#'     of tibbles, with one tibble for each depth level.
#' @param tmax A positive scalar. The cut point for paths.
#' @param equalize A logical scalar. If TRUE, make all paths end at
#'     t_max by calling \code{\link{equalize_paths}}.
#' 
#' @export
treelist_to_paths <- function(treelist,
                              tmax,
                              equalize = TRUE) {
    paths <-
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
        mutate(delta_pos = (label == "t_infect"),
               delta = ifelse(delta_pos, 1, -1)) %>%
        select(-label) %>%
        group_by(id_sim, time) %>%
        summarize_at(vars(delta, delta_pos), sum) %>%
        mutate(n_infected = cumsum(delta_pos),
               n_infectious = cumsum(delta)) %>%
        # cleanup
        ungroup %>%
        select(-starts_with("delta")) %>%
        filter(time <= tmax)
    # make all paths end at t_max
    if (equalize) {
        paths <- equalize_paths(paths, tmax)
    }
    return(paths)
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
equalize_paths <- function(df_paths, tmax) {
    # Truncate points beyond tmax
    df_paths <- df_paths %>% filter(time <= tmax)
    # Extend paths that need it
    tails <-
        df_paths %>%
        arrange(time) %>%
        group_by(id_sim) %>%
        summarize_all(~ tail(., 1)) %>%
        filter(time < tmax) %>%
        mutate(time = tmax)
    # Add extensions to truncated paths
    df_paths_equal <-
        bind_rows(df_paths, tails) %>%
        arrange (id_sim, time)
    return(df_paths_equal)
}


get_extinction_data <- function(treelist, tmax = Inf, ceil = Inf) {
    tibble(tree = map(treelist, bind_rows),
           id_sim = seq_along(tree)) %>%
        unnest(tree) %>%
        group_by(id_sim) %>%
        summarize(id_ext = which.max(t_infect + t_comm),
                  t_ext = t_infect[id_ext] + t_comm[id_ext],
                  id_max = max(id),
                  gen_ext = id_layer[id_ext],
                  gen_max = max(id_layer)) %>%
        filter(t_ext < tmax & id_max < ceil)
}



#' Convert a tree to a \code{tidygraph} object.
#'
#' @param tree A list of tibbles, each corresponding to a generation
#'     in the infection propagation
#'
#' @return a \code{tidygraph} object
#'
#' @export
as_tidygraph <- function(tree) {
    # If tree is passed as a list, this will convert it to a tibble
    tree <- tree %>% bind_rows
    # Generate `nodes`, `edges`, and `graph`
    nodes <- tree
    edges <-
        tree %>%
        select(from = id_parent, to = id) %>%
        tail(-1)
    graph <- tidygraph::tbl_graph(nodes = nodes, edges = edges)
    return(graph)
}
