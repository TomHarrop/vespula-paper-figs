library(data.table)
library(ggmap)
library(readxl)

# read data
germloc <- data.table(
    read_xlsx("data/Distribution data for Vespula wasps.xlsx",
              skip = 2,
              sheet = "Vespula germanica"))
vulgloc <- data.table(
    read_xlsx("data/Distribution data for Vespula wasps.xlsx",
              skip = 2,
              sheet = "Vespula vulgaris"))
pensloc <- data.table(
    read_xlsx("data/Distribution data for Vespula wasps.xlsx",
              skip = 2,
              sheet = "Vespula pensylvanica"))

wasploc <- rbindlist(list("Vespula germanica" = germloc,
                          "Vespula vulgaris" = vulgloc,
                          "Vespula pensylvanica" = pensloc),
                     idcol = "species")

# order the ranges
wasploc[, Range := factor(Range, levels = c("Native", "Invaded"))]

# decide limits
min_long <- wasploc[, min(Long)]
max_long <- wasploc[, max(Long)]
min_lat <- wasploc[, min(Lat)]
max_lat <- wasploc[, max(Lat)]

nudge <- 5

map_xlim <- c(min_long - nudge, max_long + nudge)
map_ylim <- c(min_lat - nudge, 60)

# choose a background colour
pal <- RColorBrewer::brewer.pal(4, "Dark2")
base_colour <- "grey80"
mappoly <- borders("world",
                   colour = alpha(base_colour, 0.75),
                   fill = alpha(base_colour, 0.5),
                   xlim = map_xlim,
                   ylim = map_ylim,
                   size = 0.1)

# draw the plt
gp <- ggplot(wasploc[!is.na(Range)],
             aes(x = Long, y = Lat, colour = species, shape = Range)) +
    theme_minimal(base_size = 8,
                  base_family = "Helvetica") +
    theme(legend.text = element_text(face = "italic"),
          legend.position = "bottom",
          legend.key.size = unit(0.5, "lines"),
          axis.text = element_blank(),
          panel.grid = element_blank()) +
    xlab(NULL) + ylab(NULL) +
    coord_quickmap() +
    facet_wrap(~ Range, ncol = 1, strip.position = "left") +
    scale_colour_brewer(palette = "Dark2") +
    scale_shape_manual(values = c(16, 17), guide = FALSE) +
    guides(colour = guide_legend(title = NULL,
                                 override.aes = list(size = 2,
                                                     shape = 16))) +
    mappoly + 
    geom_point(alpha = 0.8,
               size = 1)


saveRDS(gp, "gp.Rds")

PtToMm <- function(x){
    grid::convertUnit(grid::unit(x, "pt"), "mm", valueOnly = TRUE)
}



ggsave("vespula_worldmap.pdf",
       gp,
       device = cairo_pdf(),
       width = PtToMm(398),
       height = PtToMm(227),
       units = "mm")



