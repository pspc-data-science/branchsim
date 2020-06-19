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



#' Plot epidemic paths
#'
#' @param paths_data A \code{data.frame}, with at least the columns
#'     \code{id_sim}, \code{time}, and a quoted variable passed as
#'     \code{yvar}.
#' @param y_var Column in \code{paths_data} to use as y-variable.
#' @param n_max An integer. Maximum number of paths to display. Useful
#'     to avoid graphs that are too crowded.
#' @param tlim A scalar. 
#' @param show_hist A logical scalar. If \code{TRUE}, show the
#'     distribution of paths at the end.
#' @param show_extinct A logical scalar. Show extinct paths along with
#'     non-extinct ones.
#' @param smooth_data A \code{data.frame}, with at least the columns
#'     \code{time}, \code{n_infected}, and \code{n_infectious}.
#' @param extinct_color A color value understood by ggplot.
#' @param n_bins An integer scalar. Number of bins in the histogram
#'     panel.
#' @return A plot.
#'
#' @export
plot_paths <- function(paths_data,
                       y_var = n_infected,
                       n_max = 200,
                       tlim = NA,
                       highlight_path_id = NA,
                       show_hist = TRUE,
                       show_extinct = TRUE,
                       smooth_data = NULL,
                       extinct_color = "red",
                       n_bins = 100) {
    y_var <- enquo(y_var)
    # Sample down to n_max paths if necessary
    id_sim_all <- unique(paths_data[["id_sim"]])
    id_sim_subset <- id_sim_all[1:min(n_max, length(id_sim_all))]
    paths_data <- paths_data %>% filter(id_sim %in% id_sim_subset)
    # Tweak extinct paths y-value for plotting
    min_infectious <- exp(- .05 * max(log(paths_data[["n_infected"]])))
    paths_data <-
        paths_data %>%
        mutate(is_extinct = n_infectious == 0) %>%
        mutate(n_infectious = ifelse(n_infectious == 0,
                                     min_infectious,
                                     n_infectious))
    if (!show_extinct) {
        paths_data <- paths_data %>% filter(!is_extinct)
    }
    # Make an extra dataset just with the path(s) to highlight
    id_sim_highlight <- id_sim_subset[highlight_path_id]
    highlight_path_data <- paths_data %>% filter(id_sim %in% id_sim_highlight)
    # Pick y_var limits
    max_y <-
        paths_data %>%
        pull(!!y_var) %>%
        max
    min_y <-
        paths_data %>%
        pull(!!y_var) %>%
        min
    bin_width <- (log(max_y, 10) - log(min_y, 10)) / n_bins
    max_y <- 10^(log(max_y, 10) + bin_width)
    min_y <- 10^(log(min_y, 10) - bin_width)
    # Plot paths
    plt1 <-
        paths_data  %>%
        filter(id_sim %in% sample(unique(id_sim), min(max(id_sim), n_max))) %>%
        ggplot(aes(time, !!y_var, group = id_sim, color = is_extinct)) +
        geom_step(alpha = .3, show.legend = FALSE) +
        labs(x = "Time") +
        theme_bw() +
        scale_x_continuous(limits = c(0, tlim),
                           expand = expansion(c(0, 0))) +
        scale_y_log10(limits = c(min_y, max_y),
                      expand = expansion(c(0, .05))) +
        scale_color_manual(values = c("TRUE" = extinct_color,
                                      "FALSE" = "black"))
    if (nrow(highlight_path_data) > 0) {
        plt1 <-
            plt1 +
            geom_step(data = highlight_path_data,
                      alpha = 1,
                      size = 1,
                      show.legend = FALSE)
    }
    if (!is.null(smooth_data)) {
        plt1 <-
            plt1 +
            geom_line(data = smooth_data,
                      mapping = aes(time, !!y_var),
                      inherit.aes = FALSE,
                      color = "blue",
                      size = 1)
    }
    # Show distribution of paths position at the end
    if (show_hist) {
        # Adjust margins so the plots are joined perfectly.
        max_time <- max(paths_data[["time"]])
        margins <- theme_get()[["plot.margin"]]
        margins[4] <- unit(-8, "pt")
        plt2 <-
            paths_data %>%
            filter(time == max_time) %>%
            ggplot(aes(!!y_var, ..density..)) +
            geom_histogram(aes(fill = is_extinct),
                           binwidth = bin_width,
                           ## alpha = .5,
                           show.legend = show_extinct) +
            scale_fill_manual(values = c("TRUE" = extinct_color,
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


plot_single_path <- function(path_data,
                             tlim = NA,
                             ylim = NA,
                             scale = "log",
                             smooth_data = NULL) {
    # Tweak extinct paths y-value for plotting
    path_data <-
        path_data %>%
        mutate(kind = "path")
    if (!is.null(smooth_data)) {
        path_data <- bind_rows(path_data,
                               smooth_data %>% mutate(kind = "expectation"))
    }
    # n_infected, and n_infectious will be in a single column
    path_data <- path_data %>% pivot_longer(starts_with("n_"))
    # Plot paths
    plt <-
        path_data %>%
        ggplot(aes(time, value, linetype = name, color = kind)) +
        geom_step() +
        labs(x = "Time") +
        theme_bw() +
        scale_x_continuous(limits = c(0, tlim)) +
        scale_color_manual(values = c(path = "black", expectation = "blue"))
    if (scale == "log") {
        plt <- plt + scale_y_log10(limits = c(0.95, ylim))
    } else {
        plt <- plt + scale_y_continuous(limits = c(0, ylim))
    }
    return(plt)
}


#' Plot dendrogram corresponding to a single path
#'
#'
#' @export
plot_dendrogram <- function(single_path_data, tmax = NA) {
    # Plot as dendogram
    single_path_data %>%
        as_tidygraph %>%
        ggraph::ggraph(layout = "dendrogram",
                       height = t_infect) +
        ggraph::geom_edge_elbow2() +
        theme(axis.line.x = element_line(),
              axis.text.x = element_text(),
              axis.ticks.x = element_line(),
              axis.title.x = element_text(),
              panel.background = element_blank()) +
        labs(y = "time") +
        coord_flip(ylim = c(0, tmax))

}
