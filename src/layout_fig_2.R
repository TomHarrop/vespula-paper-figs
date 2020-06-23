#!/usr/bin/env Rscript

library(data.table)
library(cowplot)
library(ggplot2)
library(grid)
library(ggplotify)
library(ggupset)
library(tidyverse)

# chromosome ideogram
assembly_pd <- readRDS("output/plot_data/chr_plot.Rds")
end_coords <- assembly_pd[, .(starts,
                              stops = max(stops)),
                          by = .(assembly, name)]
end_offset <- 1/1e6 # in MB
ideogram <- ggplot(assembly_pd, aes(x = name,
                                    xend = name,
                                    y = starts / 1e6,
                                    yend = stops / 1e6,
                                    colour = gc)) +
    theme_minimal(base_size = 8, base_family = "Helvetica") +
    ylim(c(-5, 25)) +
    scale_y_continuous(expand = expand_scale(add = 1)) +
    theme(panel.grid.major.x = element_blank(),
          strip.text = element_text(face = "italic"),
          axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     vjust = 0.5),
          legend.key.size = unit(3, "mm"),
          legend.title = element_text(size = 5),
          legend.text = element_text(size = 4)) +
    facet_wrap(~ assembly, ncol = 1) +
    xlab(NULL) + ylab("Megabases") +
    scale_colour_viridis_c(
        guide = guide_colourbar(title = "GC\ncontent")) +
    geom_segment(data = end_coords,
                 mapping = aes(x = name,
                               xend = name,
                               y = (starts / 1e6)  - end_offset,
                               yend = (stops / 1e6) + end_offset),
                 colour = "grey90",
                 lineend = "round",
                 size = 3,
                 alpha = 0.25) +
    geom_segment(size = 3)


# synteny plot
syn_file <- "data/img/Figure2.jpg"
# syn_file <- magick::image_read_svg("data/img/Fig 2 size.svg")
synteny <- draw_image(syn_file,
                      x = 0.5,
                      y = 0.5,
                      hjust = 0.5,
                      vjust = 0.5)
my_ggdraw <- function(draw_offset) {
    ggdraw(xlim = c(0 - draw_offset, 1 + draw_offset),
           ylim = c(0 - draw_offset, 1 + draw_offset),
           clip = "on")
}
draw_offset <- 0
synteny_draw <- my_ggdraw(draw_offset) + synteny
    

# venn diagram
vesp_ortho_list <- readRDS("output/plot_data/vesp_ortho_list.Rds")
og_upset <- ggplot(as_tibble(vesp_ortho_list),
                   aes(x = Genome)) +
    theme_minimal(base_size = 8, base_family = "Helvetica") +
    theme(axis.title.y = element_text(margin = margin(0, 2, 0, 10, "pt"))) +
    xlab(NULL) + ylab("Shared\northogroups") +
    scale_x_upset() +
    geom_bar(fill = viridisLite::viridis(1)) +
    theme_combmatrix(combmatrix.label.text = element_text(size = 6,
                                                          face = "italic"),
                     combmatrix.label.make_space = FALSE, 
                     combmatrix.panel.point.color.fill = viridisLite::viridis(1),
                     combmatrix.panel.line.size = 0)

# pan / core genome
pan_core <- readRDS("output/plot_data/pan_core_pd.Rds")
pc_plot <- ggplot(pan_core, aes(x = number_of_species,
                                y = value / 1e3,
                                colour = Genome)) +
    theme_minimal(base_size = 8, base_family = "Helvetica") +
    theme(legend.key.size = unit(2/3, "lines")) +
    xlab("Number of genomes") + ylab("Orthogroups (thousands)") +
    scale_colour_viridis_d(guide = guide_legend(override.aes = list(alpha = 1))) +
    geom_point(shape =16, alpha = 0.5)


# layout
col_r <- plot_grid(
    synteny_draw,
    og_upset,
    ncol = 1,
    labels = c("B", "C"),
    label_fontfamily = "Helvetica",
    label_size = 10)
row1 <- plot_grid(
    ideogram,
    col_r,
    ncol = 2,
    rel_widths = c(4, 2),
    labels = c("A", ""),
    label_fontfamily = "Helvetica",
    label_size = 10)

cowplot <- plot_grid(
    row1,
    pc_plot,
    ncol = 1,
    rel_heights = c(3, 2),
    labels = c("", "D"),
    label_fontfamily = "Helvetica",
    label_size = 10)

ggsave("figs/figure_2.pdf",
       cowplot,
       device = cairo_pdf,
       width = 170,
       height = 225 * 3/4,
       units = "mm")

# ggsave("figs/figure_2.png",
#        cowplot,
#        device = "png",
#        dpi = 300,
#        width = 170,
#        height = 225 * 3/4,
#        units = "mm")

