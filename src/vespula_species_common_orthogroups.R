#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(ggupset)
library(tidyverse)

spec_order <- c("V. germanica" = "Vgerm",
                "V. pensylvanica" = "Vpens",
                "V. vulgaris" = "Vvulg")


gene_counts_wide <- fread("data/Orthogroups.GeneCount.tsv")
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

