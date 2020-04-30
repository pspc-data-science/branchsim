#' Plot epidemic paths
#'
#' @param paths_data A \code{data.frame}, with at least the columns
#'     \code{id_sim}, \code{date}, and a quoted variable passed as
#'     \code{yvar}.
#' @param y_var Column in \code{paths_data} to use as y-variable.
#' @param n_max An integer. Maximum number of paths to display. Useful
#'     to avoid graphs that are too crowded.
#' @param show_hist A logical scalar. If \code{TRUE}, show the
#'     distribution of paths at the end.
#' @param show_1_mass A logical scalar. Show the probability mass at
#'     \code{y_var = 1}. When \code{y_var} is \code{n_infected}, this
#'     corresponds to epidemics that did not continue past the first
#'     infected individual. NOTE: this parameter is meaningful only
#'     when \code{show_hist} is \code{TRUE}.
#'
#' @return A plot.
#'
#' @export
plot_paths <- function(paths_data,
                       y_var = n_infected,
                       n_max = 200,
                       show_hist = TRUE,
                       show_1_mass = FALSE) {
    y_var <- enquo(y_var)
    # Sample down to n_max paths if required
    if (length(unique(paths_data[["id_sim"]])) > n_max) {
        paths_plt <-
            paths_data %>%
            filter(id_sim %in% sample(unique(paths_data[["id_sim"]]), n_max))
    } else {
        paths_plt <- paths_data
    }
    # Pick y_var limits
    max_y <-
        paths_data %>%
        pull(!!y_var) %>%
        max
    min_y <-
        paths_data %>%
        pull(!!y_var) %>%
        min
    if (min_y == 1 & show_1_mass) {
        min_y <- 0.5
    }
    max_date <- max(paths_data[["date"]])
    # Plot paths
    plt1 <-
        paths_plt %>%
        ggplot(aes(date, !!y_var, group = id_sim)) +
        geom_step(alpha = .3) +
        labs(x = "Time") +
        theme_bw() +
        scale_x_continuous(expand = c(NA, max_date)) +
        scale_y_log10(limits = c(1, max_y))
    # Show distribution of paths position at the end
    if (show_hist) {
        margins <- theme_get()[["plot.margin"]]
        margins[4] <- unit(-10, "pt")
        plt2 <-
            paths_data %>%
            filter(!!y_var > min_y) %>%
            filter(date == max_date) %>%
            ggplot(aes(!!y_var, ..density..)) +
            geom_histogram() +
            scale_x_log10(limits = c(1, max_y)) +
            theme_bw() +
            coord_flip() +
            theme(axis.title.y = element_blank(),
                  axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  plot.margin = margins)
        return(egg::ggarrange(plt1 + coord_cartesian(xlim = c(0, max_date)), plt2, widths = c(5, 1)))
    } else {
        return(plt1)
    }
}
