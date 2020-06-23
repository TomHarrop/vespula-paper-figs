#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(ggupset)
library(tidyverse)

spec_order <- c("V. germanica" = "Vespula_germanica",
                "V. pensylvanica" = "Vespula_pensylvanica",
                "V. vulgaris" = "Vespula_vulgaris")


gene_counts_wide <- fread("data/Orthogroups.GeneCount.final.tsv")
gene_counts_wide[, Total := NULL]
gene_counts <- melt(gene_counts_wide,
                    id.vars = "Orthogroup",
                    variable.name = "spec_code",
                    value.name = "gene_count")[
                        spec_code %in% spec_order
                        ]
gene_counts[, has_gene := gene_count > 0]
gene_counts[, spec_code := droplevels(spec_code)]
gene_counts[, Genome := factor(plyr::mapvalues(spec_code, spec_order, names(spec_order)),
                               levels = names(spec_order))]
setkey(gene_counts, Genome)

og_lists <- dcast(gene_counts[has_gene == TRUE], Orthogroup ~ .,
      value.var = "Genome",
      fun.aggregate = list)
setnames(og_lists, ".", "Genome")

saveRDS(og_lists, "output/plot_data/vesp_ortho_list.Rds")

# counts
gene_counts_all <- melt(gene_counts_wide,
                    id.vars = "Orthogroup",
                    variable.name = "spec_code",
                    value.name = "gene_count")[spec_code != "Drosophila_melanogaster"]
gene_counts_all[, spec_code := droplevels(spec_code)]
gene_counts_all[, has_gene := gene_count > 0]

# number of genomes
gene_counts_all[, length(unique(as.character(spec_code)))]

# number of orthogroups in pan-genome
gene_counts_all[, all(has_gene), by = Orthogroup][V1 == TRUE, length(unique(Orthogroup))]

# number in vespula species
gene_counts[, sum(has_gene), by = .(spec_code)]

# number in non-vespula hymenoptera
gene_counts_all[!spec_code %in% spec_order, sum(has_gene), by = spec_code][, summary(V1)]


