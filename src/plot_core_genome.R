#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

set.seed(14)
pan_core_counts <- readRDS("output/plot_data/pan_core_counts.Rds")
var_order <- c("Single" = "Single", "n_core" = "Core", "n_pan" = "Pan")

counts_long <- melt(pan_core_counts, id.vars = c("number_of_species", "species"))

counts_no_singletons <- counts_long[!(number_of_species == 1 &
                                          variable == "n_core")]
counts_no_singletons[number_of_species == 1,
                     variable := "Single"]
counts_no_singletons[, Genome := factor(plyr::revalue(variable, var_order),
                                        levels = var_order)]

subset <- counts_no_singletons[, .I[sample(1:.N,
                                           min(1e3, .N),
                                           replace = FALSE)],
                               by = number_of_species][, c(V1)]

plot_data <- counts_no_singletons[subset]
saveRDS(plot_data, "output/plot_data/pan_core_pd.Rds")
