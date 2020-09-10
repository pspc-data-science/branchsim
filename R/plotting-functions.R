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

#' @export
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


#' Contour plot of \eqn{R_0} for different values of \eqn{\lambda} and \eqn{\mu}.
#'
#' @param lambda A vector of two values covering the range of \eqn{lambda}.
#' @param mu A vector of two values covering the range of \eqn{mu}.
#'
#' @return A ggplot object
#'
#' @export
plot_R0 <- function(my_lambda, my_mu, tbar,
                    lambda_lim = c(.01, 0.2),
                    mu_lim = c(1.01, 2.),
                    q = 0, mbar = NA, n = 32) {
    my_R0 <- calc_R0(tbar, my_lambda, my_mu, q, mbar)
    lambda_vec <- seq(lambda_lim[1], lambda_lim[2], length.out = n)
    mu_vec <- seq(mu_lim[1], mu_lim[2], length.out = n)
    dat_2d <-
        expand_grid(lambda = lambda_vec, mu = mu_vec) %>%
        mutate(R0 = calc_R0(tbar, lambda, mu, q, mbar))
    # Find (lambda, mu) coordinates corresponding to contour R0 == 1
    # For each lambda, find mu for which R0 == 1
    one_contour_lambda <-
        dat_2d %>%
        group_by(lambda) %>%
        group_nest %>%
        filter(map_lgl(data, ~ min(.[["R0"]]) < 1 & max(.[["R0"]]) > 1)) %>%
        mutate(mu = map_dbl(data, ~ approx(.[["R0"]],
                                           .[["mu"]], 1)[["y"]])) %>%
        # cleanup
        select(-data)
    # For each mu, find lambda for which R0 == 1
    one_contour_mu <-
        dat_2d %>%
        group_by(mu) %>%
        group_nest %>%
        filter(map_lgl(data, ~ min(.[["R0"]]) < 1 & max(.[["R0"]]) > 1)) %>%
        mutate(lambda = map_dbl(data, ~ approx(.[["R0"]],
                                               .[["lambda"]], 1)[["y"]])) %>%
        # cleanup
        select(-data)
    # Join the two tibbles
    one_contour <-
        bind_rows(one_contour_lambda, one_contour_mu) %>%
        arrange(lambda) %>%
        unique
    # Label position for R0 value
    label_pos <- one_contour %>% filter(mu == quantile(mu, .5, type = 1))
    label_nudge_y <- 0.05 * abs(mu_lim[2] - mu_lim[1])
    ggplot(dat_2d, aes(lambda, mu, fill = R0)) +
        ## geom_contour_filled() +
        geom_raster(interpolate = TRUE) +
        geom_line(aes(lambda, mu),
                  data = one_contour,
                  inherit.aes = FALSE,
                  color = "black",
                  linetype = 3,
                  size = 2) +
        ## annotate("text",
        ##          x = label_pos[["lambda"]],
        ##          y = label_pos[["mu"]],
        ##          label = bquote(R[0] == 1),
        ##          color = "white",
        ##          vjust = 0,
        ##          hjust = 0,
        ##          size = 10) +
        annotate("point",
                 x = my_lambda,
                 y = my_mu,
                 size = 8,
                 color = "black",
                 shape = 8) +
        annotate("text",
                 y = ifelse(my_R0 <= 1,
                            my_mu - label_nudge_y,
                            my_mu + label_nudge_y),
                 x = my_lambda,
                 hjust = 0.5,
                 vjust = ifelse(my_R0 <= 1, 1, 0),
                 label = bquote(R[eff] == .(my_R0)),
                 color = "black",
                 size = 10) +
        theme_bw() +
        scale_fill_gradient2(low = muted("blue"),
                             high = muted("red"),
                             midpoint = 1) +
        # For math notation and bquote, see
        # https://trinkerrstuff.wordpress.com/2018/03/15/2246/
        labs(x = "Average number of infectious events per day",
             y = "Average number infected per infectious event",
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
#' @param mu Average number of individuals infected
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
