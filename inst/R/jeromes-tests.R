library(tidyverse)
library(tidygraph)
library(ggraph)
## ggraph::set_graph_style(plot_margin = margin(1,1,1,1))

devtools::load_all()

set.seed(123)
sims <-
    branchsim::run_sims(nsim = 200, tmax = 25, q = 0) %>%
    mutate(paths = map2(treelist, tmax, branchsim::treelist_to_paths)) %>%
    mutate(smooth = pmap(list(lambda,
                              -p / ((1 - p) * log(1 - p)),
                              tbar * kappa,
                              kappa,
                              tmax), branchsim::calc_avg_path))


plot_single_path(sims[["paths"]][[1]] %>% filter(id_sim == 2) %>% branchsim::equalize_paths(10), smooth_data = sims[["smooth"]][[1]] %>% filter(time <= 10), ylim = max(sims[["paths"]][[1]][["n_infected"]]), tlim = 25)


## # Plot with energy minimization algo
## sims[["treelist"]][[1]][[6]] %>%
##     as_tidygraph %>%
##     mutate(special = id %in% igraph::subcomponent(., 4, mode = "out")) %>% print %>%
##     ggraph::ggraph(layout = "stress") + #, height = t_infect) +
##     ggraph::geom_edge_link() +
##     ggraph::geom_node_point() + 
##     ## ggraph::geom_edge_diagonal(aes(colour = node.special), show.legend = FALSE) +
##     coord_cartesian(ylim = c(0, 20)) +
##     scale_edge_colour_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
##     coord_flip()
