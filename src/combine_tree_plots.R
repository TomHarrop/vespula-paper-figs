#!/usr/bin/env Rscript

library(cowplot)
library(ggplot2)
library(ggtree)

p450_tree <- readRDS("output/plot_data/p450_tree.Rds")
tree_plots <- readRDS("output/plot_data/tree_plots.Rds")

# get the legend for the gene trees
legend <- get_legend(tree_plots[[1]] +
                         theme(legend.direction = "horizontal",
                               legend.justification = "centre",
                               legend.key.width = unit(25, "mm"),
                               legend.key.height = unit(0, "mm")))

# top row: the species tree
p_row1 <- p450_tree +  xlim(c(0, 1.7))

# row two: compound plot of trees
p_col1 <- plot_grid(
    tree_plots[[1]] + xlim(c(0,1.3)) + theme(legend.position = "none"),
    labels = "B"
)
p_col2 <- plot_grid(plotlist = list(
    tree_plots[[2]] + xlim(c(0, 2.9)) + theme(legend.position = "none"),
    tree_plots[[3]] + xlim(c(0,1.3)) + theme(legend.position = "none")),
    align = "hv",
    ncol = 1,
    labels = c("C", "D"))

p_row2 <- plot_grid(
    p_col1, p_col2,
    align = "hv",
    ncol = 2,
    rel_widths = c(4, 3)
)

cowplot <- plot_grid(
    p_row1,
    p_row2,
    legend,
    ncol = 1,
    labels = c("A"),
    # align = "h",
    axis = "b",
    rel_heights = c(1, 0.9, 0.1))

ggsave("figs/P450_supplementary_figure.pdf",
       cowplot,
       device = cairo_pdf,
       width = 170,
       height = 225,
       units = "mm")
