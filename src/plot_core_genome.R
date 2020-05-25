#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

set.seed(14)

gene_counts_wide <- fread("data/Orthogroups.GeneCount.tsv")
pan_core_counts <- readRDS("output/plot_data/pan_core_counts.Rds")

# count the number in all species
gene_counts_wide[, Total := NULL]
gene_counts_wide[, Dmel := NULL]
gene_counts <- melt(gene_counts_wide,
                    id.vars = "Orthogroup",
                    variable.name = "spec_code",
                    value.name = "gene_count")
gene_counts[, has_gene := gene_count > 0]
setkey(gene_counts, spec_code)
all_spec <- gene_counts[, unique(spec_code)]
count_all <- gene_counts[,
            .(core = all(has_gene),
              pan = any(has_gene)),
            by = Orthogroup][, .(
                number_of_species = length(all_spec),
                species = paste(all_spec, collapse = "_"),
                n_core = sum(core),
                n_pan = sum(pan))]


var_order <- c("Single" = "Single", "n_core" = "Core", "n_pan" = "Pan")

counts_long <- melt(rbind(pan_core_counts, count_all),
                    id.vars = c("number_of_species", "species"))

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
