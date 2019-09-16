library(seqinr)
library(data.table)
library(ggplot2)
library(future.apply)

CalculateGCWindows <- function(seq_object, window_size) {
    print(getName(seq_object))
    my_seqlength <- length(seq_object)
    
    # from https://github.com/avrilcoghlan/LittleBookofRBioinformatics/blob/master/src/chapter2_answers.rst
    my_starts <- seq(1, my_seqlength - window_size, by = window_size)
    my_stops <- seq(1 + window_size, my_seqlength, by = window_size) - 1
    my_stops[length(my_stops)] <- my_seqlength
    seq_ranges <- data.table(name = getName(seq_object),
                             starts = my_starts,
                             stops = my_stops)
    pd <- seq_ranges[, .(gc = GC(seq_object[starts:stops])),
                     by = .(name, starts, stops)]
    return(pd)
}

GetGCPlotData <- function(seq_object, window_size = 1e5){
    my_chr <- seq_object[grepl("^Chr", getName(seq_object))]
    gc_list <- future_lapply(my_chr, function(x)
        CalculateGCWindows(x, window_size))
    pd <- rbindlist(gc_list)
    chr_only <- pd[grepl("^Chr", name)]
    chr_only[, chr_num := as.integer(gsub("[^[:digit:]]+", "", name))]
    setorder(chr_only, chr_num)
    chr_only[, name := factor(name, levels = unique(name))]
    return(chr_only)
}

PtToMm <- function(x){
    grid::convertUnit(grid::unit(x, "pt"), "mm", valueOnly = TRUE)
}

###########
# GLOBALS #
###########

assembly_files <- list.files("data/assemblies",
                             pattern = "assembly",
                             full.names = TRUE)

########
# MAIN #
########

plan(multiprocess)
options(future.globals.maxSize = +Inf)

names(assembly_files) <- sub("\\..+", "", basename(assembly_files))

# read the fasta files
assembly_fastas <- future_lapply(assembly_files, read.fasta)
assembly_pd_list <- future_lapply(assembly_fastas,
                                  GetGCPlotData,
                                  window_size = 1e5)
assembly_pd <- rbindlist(assembly_pd_list, idcol = "assembly")

# rename the assemblies
spec_order <- c("Vgerm" = "V. germanica",
  "Vpens" = "V. pensylvanica",
  "Vvulg" = "V. vulgaris")

assembly_pd[, assembly := factor(plyr::revalue(assembly, spec_order),
                     levels = spec_order)]

# get the ends
end_coords <- assembly_pd[, .(starts,
                              stops = max(stops)),
                          by = .(assembly, name)]
end_offset <- 0.005 # in MB

saveRDS(assembly_pd, "output/plot_data/chr_plot.Rds")

gp <- ggplot(assembly_pd, aes(x = name,
                        xend = name,
                        y = starts / 1e6,
                        yend = stops / 1e6,
                        colour = gc)) +
    theme_minimal(base_size = 8) +
    theme(panel.grid.major.x = element_blank(),
          strip.text = element_text(face = "italic"),
          axis.text.x = element_text(angle = 30,
                                     hjust = 1,
                                     vjust = 1),
          legend.key.size = unit(3, "mm"),
          legend.title = element_text(size = 5),
          legend.text = element_text(size = 4)) +
    facet_wrap(~ assembly, ncol = 1) +
    xlab(NULL) + ylab("Megabases") +
    scale_colour_viridis_c(
        guide = guide_colourbar(title = "GC content")) +
    geom_segment(data = end_coords,
                 mapping = aes(x = name,
                               xend = name,
                               y = starts / 1e6  - end_offset,
                               yend = stops / 1e6 + end_offset),
                 colour = "grey90",
                 lineend = "round",
                 size = 4,
                 alpha = 0.25) +
    geom_segment(size = 4)

ggsave("figs/assembly_ideogram.pdf",
       gp,
       device = cairo_pdf(),
       width = 150,
       height = 100,
       units = "mm")

