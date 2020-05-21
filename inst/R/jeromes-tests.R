library(tidygraph)
library(ggraph)
ggraph::set_graph_style(plot_margin = margin(1,1,1,1))

set.seed(123)
branchsim::run_sims(nsim = 10, tmax = 20)

# Plot as dendogram
sims[["treelist"]][[1]][[6]] %>%
    as_tidygraph %>%
    mutate(special = id %in% igraph::subcomponent(., 4, mode = "out")) %>% print %>%
    ggraph::ggraph(layout = "dendrogram", height = t_infect) +
    ggraph::geom_edge_elbow2(aes(colour = node.special), show.legend = FALSE) +
    coord_cartesian(ylim = c(0, 20)) +
    scale_edge_colour_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
    coord_flip()

# Plot with energy minimization algo
sims[["treelist"]][[1]][[6]] %>%
    as_tidygraph %>%
    mutate(special = id %in% igraph::subcomponent(., 4, mode = "out")) %>% print %>%
    ggraph::ggraph(layout = "stress") + #, height = t_infect) +
    ggraph::geom_edge_link() +
    ggraph::geom_node_point() + 
    ## ggraph::geom_edge_diagonal(aes(colour = node.special), show.legend = FALSE) +
    coord_cartesian(ylim = c(0, 20)) +
    scale_edge_colour_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
    coord_flip()
