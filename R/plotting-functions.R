#' Plot epidemic paths
#'
#' @param paths_data A \code{data.frame}, with at least the columns
#'     \code{id_sim}, \code{time}, and a quoted variable passed as
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
                       show_extinct = TRUE) {
    # Sample down to n_max paths if required
    if (length(unique(paths_data[["id_sim"]])) > n_max) {
        paths_plt <-
            paths_data %>%
            filter(id_sim %in% sample(unique(paths_data[["id_sim"]]), n_max))
    } else {
        paths_plt <- paths_data
    }
    # Tweak extinct paths y-value for plotting
    min_infectious = exp(- .05 * max(log(paths_data[["n_infected"]])))
    paths_plt <-
        paths_data %>%
        mutate(is_extinct = n_infectious == 0) %>%
        mutate(n_infectious = ifelse(n_infectious == 0,
                                     min_infectious,
                                     n_infectious))
    if (!show_extinct) {
        paths_plt <- paths_plt %>% filter(!is_extinct)
    }
    y_var <- enquo(y_var)
    # Pick y_var limits
    max_y <-
        paths_plt %>%
        pull(!!y_var) %>%
        max
    min_y <-
        paths_plt %>%
        pull(!!y_var) %>%
        min
    bin_width <- (log(max_y, 10) - log(min_y, 10)) / 100
    max_y <- 10^(log(max_y, 10) + bin_width)
    min_y <- 10^(log(min_y, 10) - bin_width)
    # Plot paths
    plt1 <-
        paths_plt  %>%
        filter(id_sim %in% sample(unique(id_sim), min(max(id_sim), n_max))) %>%
        ggplot(aes(time, !!y_var, group = id_sim, color = is_extinct)) +
        geom_step(alpha = .3, show.legend = FALSE) +
        labs(x = "Time") +
        theme_bw() +
        scale_x_continuous(expand = expansion(c(0, 0))) +
        scale_y_log10(limits = c(min_y, max_y),
                      expand = expansion(c(0, .05))) +
        scale_color_manual(values = c("TRUE" = "red",
                                         "FALSE" = "black"))
    # Show distribution of paths position at the end
    if (show_hist) {
        # Adjust margins so the plots are joined perfectly.
        max_time <- max(paths_plt[["time"]])
        margins <- theme_get()[["plot.margin"]]
        margins[4] <- unit(-8, "pt")
        plt2 <-
            paths_plt %>%
            filter(time == max_time) %>%
            ggplot(aes(!!y_var, ..density..)) +
            geom_histogram(aes(fill = is_extinct),
                           binwidth = bin_width,
                           alpha = .5,
                           show.legend = show_extinct) +
            scale_fill_manual(values = c("TRUE" = "red",
                                         "FALSE" = "black")) +
            scale_x_log10(limits = c(min_y, max_y),
                          expand = expansion(c(0, .05))) +
            scale_y_continuous(expand = expansion(mult = c(0, .05))) +
            labs(fill = "extinct") +
            theme_bw() +
            coord_flip() +
            theme(axis.title.y = element_blank(),
                  axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.text.x = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  plot.margin = margins)
        return(egg::ggarrange(plt1, plt2, widths = c(5, 1)))
    } else {
        return(plt1)
    }
}
