#' Plot epidemic paths
#'
#' @param paths_data A \code{data.frame}, with at least the columns
#'     \code{id_sim}, \code{time}, and a quoted variable passed as
#'     \code{yvar}.
#' @param y_var Column in \code{paths_data} to use as y-variable.
#' @param n_max An integer. Maximum number of paths to display. Useful
#'     to avoid graphs that are too crowded.
#' @param tlim A scalar.
#' @param highlight_path_id A subset of values in column \code{id_sim}
#'     of input \code{paths_data} that will be highlighted (black and
#'     thicker). If the values provided are not found in column
#'     \code{id_sim}, no path will be highlighted.
#' @param show_hist A logical scalar. If \code{TRUE}, show the
#'     distribution of paths at the right-end of the plot. In this
#'     case, the output will be a \code{ggarrange} object, rather than
#'     a \code{ggplot} object.
#' @param show_extinct A logical scalar. Show extinct paths along with
#'     non-extinct ones.
#' @param smooth_data A \code{data.frame}, with at least the columns
#'     \code{time}, \code{n_infected}, and \code{n_infectious}. If
#'     given, will be superposed on top of the paths plot.
#' @param extinct_color A color value understood by ggplot. This will
#'     be the colour of paths that went extinct before the end of the
#'     simulation.
#' @param n_bins An integer scalar. Number of bins in the histogram
#'     panel.
#' 
#' @return A plot. If \code{show_hist = TRUE} (default), the output
#'     will be a \code{ggarrange} object showing both the paths, and a
#'     histogram of the level at which they terminate. Otherwise the
#'     output will be a \code{ggplot} object showing only the paths.
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
                       n_bins = 20) {
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
        pull({{y_var}}) %>%
        max
    min_y <-
        paths_data %>%
        pull({{y_var}}) %>%
        min
    bin_width <- (log(max_y, 10) - log(min_y, 10)) / n_bins
    max_y <- 10^(log(max_y, 10) + bin_width)
    min_y <- 10^(log(min_y, 10) - bin_width)
    # Plot paths
    plt1 <-
        paths_data  %>%
        filter(id_sim %in% sample(unique(id_sim), min(max(id_sim), n_max))) %>%
        ggplot(aes(time, {{y_var}}, group = id_sim, color = is_extinct)) +
        geom_step(alpha = .3, show.legend = FALSE) +
        # Fancy y label for the two most common cases (n_infected, n_infectious)
        labs(x = "Time",
             y = rlang::as_name(y_var) %>%
                 stringr::str_replace("n_infected", "Cumulative infection count") %>%
                 stringr::str_replace("n_infectious", "Ongoing infection count")) +
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
                      mapping = aes(time, {{y_var}}),
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
            # Do not use `..density..` as y aesthetic below because
            # it shows the wrong relative proportions between extinct /
            # non-extinct cases
            ggplot(aes({{y_var}})) +
            geom_histogram(aes(fill = is_extinct),
                           binwidth = bin_width,
                           show.legend = show_extinct) +
            scale_fill_manual(values = c("TRUE" = extinct_color,
                                         "FALSE" = "black")) +
            scale_x_log10(limits = c(min_y, max_y),
                          expand = expansion(c(0, .05))) +
            scale_y_continuous(expand = expansion(mult = c(0, .05))) +
            labs(fill = "extinct") +
            theme_bw() +
            coord_flip() +
            labs(y = "density") +
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

#' Show infection counts (cumulative and ongoing) for a single path
#'
#' @param path_data A \code{data.frame}, with at least the columns
#'     \code{time}, \code{n_infected}, and \code{n_infectious}
#' @param tlim Upper time limit to display. If \code{NA}, show the
#'     entire path.
#' @param ylim Explicit limit for the top of the graph. Useful when
#'     showing the detail of one path in comparison with the output of
#'     \code{plot_paths}, so that the y-scales are the same between
#'     the two plots.
#' @param scale If equal to "log", show in log scale (default). If any
#'     other value, show in linear scale. Like argument \code{ylim},
#'     this is useful when comparing to the output of
#'     \code{plot_paths}.
#' @param smooth_data  A \code{data.frame}, with at least the columns
#'     \code{time}, \code{n_infected}, and \code{n_infectious}. If
#'     given, will be superposed on top of the paths plot.
#'
#' @return A \code{ggplot} plot showing the cumulative infection
#'     count, the ongoing infection count, and, if \code{smooth_data}
#'     is given, the expectation values of the stochastic process.
#' 
#' @export
plot_single_path <- function(path_data,
                             tlim = NA,
                             ylim = NA,
                             scale = "log",
                             smooth_data = NULL) {
    # n_infected, and n_infectious will be in a single column
    path_data <-
        path_data %>%
        pivot_longer(starts_with("n_"),
                     names_to = "count_type",
                     values_to = "count") %>%
        mutate(count_type = count_type %>%
                   stringr::str_replace("n_infected", "Cumulative infection count") %>%
                   stringr::str_replace("n_infectious", "Ongoing infection count"))
    # If `smooth_data` is non-null, we want to differentiate paths
    # from the expected value coming from the renewal equation.
    if (!is.null(smooth_data)) {
        path_data <- bind_rows(path_data %>% mutate(kind = "path"),
                               smooth_data %>% mutate(kind = "expectation"))
        plt_proto <- path_data %>% ggplot(aes(time, count, linetype = count_type,
                                              color = kind))
    } else {
        plt_proto <- path_data %>% ggplot(aes(time, count, linetype = count_type))
    }    
    # Plot paths
    plt <-
        plt_proto +
        geom_step() +
        labs(x = "Time",
             linetype = "Count type",
             color = "") +
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


#' Plot dendrogram for a single path
#'
#' This function is in construction. The main inconvenience is that
#' the communicable window of each node is only shown if it had at
#' least one offspring.
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

#' Contour plot of \eqn{R_eff} for different values of \eqn{\lambda} and \eqn{\mu}.
#'
#' The plot shows a colour map of \eqn{R_eff} across a range of values
#' for the rate at which infectious interactions occur, and the 0.95th
#' quantile for the number of individuals infected at each infectious
#' event. The plot also highlights a specific value of \eqn{R_eff},
#' corresponding to \eqn{\lambda =} \code{my_lambda}, and 0.95th
#' quantile \code{my_q95}. The value of \eqn{R_eff} is shown in shades
#' of red for values greater than one, and blue for values smaller
#' than one.
#'
#' @param my_lambda The Poisson rate of infectious interactions.
#' @param my_q95 The 0.95th quantile of the logarithmic distribution
#'     governing the
#' @param tbar The average duration of the communicable period.
#' @param lambda_lim A vector of two values covering the range of
#'     \eqn{lambda}.
#' @param q95_lim A vector of two values covering the range of the
#'     0.95th quantile of the logarithmic distribution governing the
#'     number of new infections occurring in each infectious
#'     interaction.
#' @param q The probability of an infection belonging to a second
#'     class of infections, with average duration \code{mbar}, instead
#'     of \code{tbar} (e.g. when modelling contact tracing).
#' @param mbar Average value of the communicable period for infections
#'     of the second kind. This parameter is unused when \code{q = 0}.
#' @param n Grid size for calculating \eqn{R_eff} values.
#'
#' @return A ggplot object.
#'
#' @export
plot_R0 <- function(my_lambda,
                    my_q95,
                    tbar,
                    lambda_lim = c(.02, .2),
                    q95_lim = c(2, 10),
                    q = 0,
                    mbar = NA,
                    n = 32) {
    my_p <- find_p_from_quantile(my_q95 - 1)
    my_mu <- calc_mu(my_p)
    my_R0 <- branchsim::calc_R0(tbar, my_lambda, my_mu, q, mbar)
    lambda_vec <- seq(lambda_lim[1],
                      lambda_lim[2],
                      length.out = n)
    q95_vec <- seq(q95_lim[1], q95_lim[2])
    p_vec <- map_dbl(q95_vec - 1, find_p_from_quantile)
    mu_vec <- calc_mu(p_vec)
    dat_2d <-
        expand_grid(lambda = lambda_vec,
                    q95 = q95_vec) %>%
        mutate(p = map_dbl(q95, find_p_from_quantile),
               mu = calc_mu(p),
               R0 = map2_dbl(lambda, mu,
                             ~ branchsim::calc_R0(tbar, .x, .y, q, mbar)))
    # Find (lambda, mu) coordinates corresponding to contour R0 == 1
    # For each lambda, find mu for which R0 == 1
    one_contour_lambda <-
        dat_2d %>%
        select(lambda, q95, R0) %>%
        group_by(lambda) %>%
        group_nest %>%
        filter(map_lgl(data, ~ min(.[["R0"]]) < 1 & max(.[["R0"]]) > 1)) %>%
        mutate(q95 = map_dbl(data, ~ approx(.[["R0"]],
                                            .[["q95"]], 1)[["y"]])) %>%
        # cleanup
        select(-data)
    # For each mu, find lambda for which R0 == 1
    one_contour_q95 <-
        dat_2d %>%
        select(lambda, q95, R0) %>%
        group_by(q95) %>%
        group_nest %>%
        filter(map_lgl(data, ~ min(.[["R0"]]) < 1 & max(.[["R0"]]) > 1)) %>%
        mutate(lambda = map_dbl(data, ~ approx(.[["R0"]],
                                               .[["lambda"]], 1)[["y"]])) %>%
        # cleanup
        select(-data)
    # Join the two tibbles
    one_contour <-
        bind_rows(one_contour_lambda, one_contour_q95) %>%
        arrange(lambda) %>%
        unique
    # Label position for R0 value
    label_pos <- one_contour %>% filter(q95 == quantile(q95, .5, type = 1))
    label_nudge_y <- 0.05 * abs(q95_lim[2] - q95_lim[1])
    ggplot(dat_2d, aes(lambda, q95, fill = R0)) +
        ## geom_contour_filled() +
        geom_raster(interpolate = TRUE) +
        geom_line(aes(lambda, q95),
                  data = one_contour,
                  inherit.aes = FALSE,
                  color = "black",
                  linetype = 3,
                  size = 2) +
        annotate("point",
                 x = my_lambda,
                 y = my_q95,
                 size = 8,
                 color = "black",
                 shape = 8) +
        annotate("text",
                 x = my_lambda,
                 y = my_q95,
                 hjust = "inward",
                 vjust = "inward",
                 label = bquote(R[eff] == .(format(my_R0, digits = 2))),
                 color = "black",
                 size = 10) +
        theme_bw() +
        scale_x_continuous(breaks = seq(lambda_lim[1], lambda_lim[2])) +
        scale_y_continuous(breaks = q95_vec) +
        scale_fill_gradient2(low = scales::muted("blue"),
                             high = scales::muted("red"),
                             trans = "log",
                             midpoint = 0,
                             breaks = c(.1, .3, 1, 3, 10)) +
        # For math notation and bquote, see
        # https://trinkerrstuff.wordpress.com/2018/03/15/2246/
        labs(x = "Average number of infectious interactions per day",
             y = "95% of infectious interactions involve at most",
             fill = bquote(R[eff]))
}


#' Plot the Gamma density for the communicable window
#'
#' @param tbar Average duration of the communicable window.
#' @param kappa Rate parameter of the Gamma distribution.
#' @param pmax A scalar between 0 and 1. Upper limit for the time axis.
#' @param n Number of points to plot.
#'
#' @return A ggplot object
#'
#' @export
plot_gamma <- function(tbar, kappa, pmax = 0.999, n = 1000) {
    # shape
    alpha <- tbar * kappa
    beta <- kappa
    xmax <- qgamma(pmax, alpha, beta)
    dat <- tibble(time = seq(from = 0, to = xmax, length.out = 1000),
                  f = dgamma(time, alpha, beta))
    ggplot(dat, aes(time, f)) +
        geom_line() +
        theme_bw() +
        geom_vline(xintercept = tbar, linetype = 2) +
        annotate("text",
                 x = tbar * 1.02,
                 y = dgamma(tbar, alpha, beta) * 1.02,
                 hjust = 0,
                 vjust = 1,
                 label = str_c("mean: ", tbar)) +
        labs(x = "time",
             y = "probability")
}


#' Plot the time between infectious events
#'
#' @param lambda Poisson rate
#' @param pmax A scalar between 0 and 1. Upper limit for the time axis.
#' @param n Number of points to plot.
#'
#' @return A ggplot object
#'
#' @export
plot_exp <- function(lambda, pmax = 0.995, n = 1000) {
    xmax <- qexp(pmax, lambda)
    xavg <- 1 / lambda
    dat <- tibble(time = seq(from = 0, to = xmax, length.out = n),
                  f = dexp(time, lambda))
    ggplot(dat, aes(time, f)) +
        geom_line() +
        theme_bw() +
        geom_vline(xintercept = 1 / lambda,
                   linetype = 2) +
        annotate("text",
                 x = xavg,
                 y = mean(c(dexp(0, lambda), dexp(xavg, lambda))),
                 hjust = 0,
                 label = str_c("mean time between events: ",
                               format(xavg, digits = 2))) +
        labs(x = "time between infectious events",
             y = "probability")
}


#' Plot the number of individuals infected per infectious event
#'
#' @param p Parameter of the logarithmic distribution.
#' @param q_high Upper quantile to bound the plot.
#'
#' @return A ggplot object
#'
#' @export
plot_log <- function(p, q_high = 0.995) {
    mu <- -p / ((1 - p) * log(1 - p))
    x_high <- extraDistr::qlgser(q_high, p)
    dat <- tibble(n = seq(1, x_high, by = 1),
                  f = extraDistr::dlgser(n, p))
    ggplot(dat, aes(n, f)) +
        geom_segment(aes(x = n, xend = n,
                         y = 0, yend = f)) +
        geom_point(aes(x = n, y = f)) +
        geom_vline(xintercept = mu,
                   linetype = 2) +
        annotate("text",
                 x = mu * 1.02,
                 y = extraDistr::dlgser(floor(mu), p),
                 hjust = 0,
                 label = str_c("mean : ", format(mu, digits = 2))) +
        theme_bw() +
        scale_x_continuous(breaks = seq(1, x_high, by = 1)) +
        labs(x = "infections per event",
             y = "probability")
}
