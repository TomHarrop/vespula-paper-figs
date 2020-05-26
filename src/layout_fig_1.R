#!/usr/bin/env Rscript

library(cowplot)
library(ggmap)
library(ggplot2)
library(ggtree)

# species
germanica <- draw_image("data/img/germanica3_cropped.jpg",
                        height = 1,
                        clip = "on")
vulgaris <- draw_image("data/img/vulgaris_cropped2.jpg",
                       height = 1,
                       clip = "on")
pensylvanica <- draw_image("data/img/pensylvanica_rotated.jpg",
                           height = 1,
                           clip = "on")

# map panel
map_plot <- readRDS("output/plot_data/map_plot.Rds")
mappoly <- readRDS("output/plot_data/map_poly.Rds")
vesp_map <- map_plot +
    theme_minimal(base_size = 8,
                  base_family = "Helvetica") +
    theme(legend.text = element_text(face = "italic"),
          legend.position = "bottom",
          legend.direction = "vertical",
          legend.key.size = unit(2/3, "lines"),
          axis.text = element_blank(),
          panel.grid = element_blank()) +
    xlab(NULL) + ylab(NULL) +
    facet_wrap(~ Range, ncol = 1, strip.position = "left") +
    scale_colour_viridis_d() +
    scale_shape_manual(values = c(16, 17), guide = FALSE) +
    guides(colour = guide_legend(title = NULL,
                                 override.aes = list(size = 2,
                                                     shape = 16))) +
    mappoly + 
    geom_point(alpha = 0.8,
               size = 1)


# tree panel
spec_tree <- readRDS("output/plot_data/all_species_tree.Rds")

# to view nodes
# gt + geom_nodelab(aes(label = node))

MyCladelab <- function(my_lab, my_node) {
    geom_cladelabel(
        node = my_node,
        label = my_lab,
        offset = 0.7,
        offset.text = 0.025,
        align = TRUE,
        angle = -90,
        hjust = 0.5,
        fontsize = 2
    )
}

# set up tip colours
vcols <- structure(viridisLite::viridis(3),
                   names = c("Vespula germanica",
                             "Vespula pensylvanica",
                             "Vespula vulgaris"))

tree_names <- spec_tree$data$spec_name
missing_names <- tree_names[!is.na(tree_names) & !tree_names %in% names(vcols)]

all_cols <- c(vcols,
              structure(rep(NA, length(missing_names)),
                        names = missing_names))

spec_plot <- spec_tree +
    theme(text = element_text(family = "Helvetica")) +
    geom_tree(size = 0.5) +
    scale_colour_manual(values = all_cols,
                        na.value = "black",
                        guide = FALSE) +
    geom_tiplab(aes(label = spec_name,
                    colour = spec_name),
                size = 2.5,
                fontface = "italic") + 
    geom_treescale(fontsize = 2,
                   linesize = 0.5,
                   x = 0,
                   y = 6,
                   offset = 0.1) +
    MyCladelab("Formicidae", 36) +
    MyCladelab("Apoidea", 31) +
    MyCladelab("Vespidae", 46)

# layout
my_ggdraw <- function(draw_offset) {
    ggdraw(xlim = c(0 - draw_offset, 1 + draw_offset),
           ylim = c(0 - draw_offset, 1 + draw_offset),
           clip = "on")
}

draw_offset <- 1/100

row1 <- plot_grid(
    my_ggdraw(draw_offset) + germanica, 
    my_ggdraw(draw_offset) + pensylvanica, 
    my_ggdraw(draw_offset) + vulgaris,
    ncol = 3,
    labels = c("A", "B", "C"),
    label_fontfamily = "Helvetica",
    label_colour = "white",
    label_size = 10)

row2 <- plot_grid(
    vesp_map,
    spec_plot,
    ncol = 2,
    labels = c("D", "E"),
    label_fontfamily = "Helvetica",
    label_size = 10,
    # align = "h",
    axis = "b")

cowplot <- plot_grid(row1,
                     row2,
                     ncol = 1,
                     rel_heights = c(4, 10))

ggsave("figs/figure_1.pdf",
       cowplot,
       device = cairo_pdf,
       width = 170,
       height = 225 * 2/3,
       units = "mm")

