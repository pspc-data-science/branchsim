#' Plot epidemic paths
#'
#' @param paths_data A \code{data.frame}, with at least the columns
#'     \code{id_sim}, \code{date}, and a quoted variable passed as
#'     \code{yvar}.
#' @param y_var Column in \code{paths_data} to use as y-variable.
#' @param n_max An integer. Maximum number of paths to display. Useful
#'     to avoid graphs that are too crowded.
#' @export
plot_paths <- function(paths_data,
                       y_var = n_infected,
                       n_max = 200) {
    y_var = enquo(y_var)
    # Sample down to n_max paths if required
    if (length(unique(paths_data[["id_sim"]])) > n_max) {
        paths_data <-
            paths_data %>%
            filter(id_sim %in% sample(unique(paths_data[["id_sim"]]), n_max))
    }
    ggplot(paths_data, aes(date, !!y_var, group = id_sim)) +
        geom_step(alpha = .3) +
        labs(x = "Time") +
        theme_bw() +
        scale_y_log10()
}
