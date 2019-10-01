#!/usr/bin/env Rscript

library(data.table)
library(gtools)
library(future.apply)

plan(multiprocess(workers = 4))
options(future.globals.maxSize = +Inf)

gene_counts_wide <- fread("data/Orthogroups.GeneCount.tsv")
gene_counts_wide[, Total := NULL]
gene_counts_wide[, Dmel := NULL]
gene_counts <- melt(gene_counts_wide,
                    id.vars = "Orthogroup",
                    variable.name = "spec_code",
                    value.name = "gene_count")
gene_counts[, has_gene := gene_count > 0]
setkey(gene_counts, spec_code)

GetComboCounts <- function(gene_counts, spec_no) {
    print(paste("spec_no", spec_no))
    comb_res <- gene_counts[, combinations(
        length(as.character(unique(spec_code))),
        spec_no,
        as.character(unique(spec_code)),
        set = TRUE,
        repeats.allowed = FALSE)]
    my_count_list <- future_apply(
        X = comb_res,
        MARGIN = 1,
        future.packages = "data.table",
        FUN = function(x)
        {
            gene_counts[spec_code %in% x,
                        .(core = all(has_gene),
                          pan = any(has_gene)),
                        by = Orthogroup][, .(
                            number_of_species = length(x),
                            species = paste(x, collapse = "_"),
                            n_core = sum(core),
                            n_pan = sum(pan))]})
    
    my_counts <- rbindlist(my_count_list, fill = TRUE)
    return(my_counts)}

number_of_specs <- gene_counts[, length(unique(spec_code))]

all_count_list <- lapply(rev(1:(number_of_specs - 1)), function(x)
    GetComboCounts(gene_counts, x))
all_counts <- rbindlist(all_count_list, fill = TRUE)

saveRDS(all_counts, "output/plot_data/pan_core_counts.Rds")
