#!/usr/bin/env Rscript

library(data.table)
library(gtools)
library(future.apply)

plan(multiprocess(workers = 100))
#plan(multiprocess)
options(future.globals.maxSize = +Inf)

gene_counts_wide <- fread("data/Orthogroups.GeneCount.final.tsv")
gene_counts_wide[, Total := NULL]
gene_counts_wide[, Drosophila_melanogaster := NULL]
gene_counts <- melt(gene_counts_wide,
                    id.vars = "Orthogroup",
                    variable.name = "spec_code",
                    value.name = "gene_count")
gene_counts[, has_gene := gene_count > 0]

# convert spec_code to numeric (save mem?)
gene_counts[, spec_id := as.numeric(as.factor(spec_code))]
spec_code_to_id <- unique(gene_counts[, .(spec_code, spec_id)])
fwrite(spec_code_to_id, "output/plot_data/spec_code_to_id.csv")

gene_counts[, spec_code := NULL]
setnames(gene_counts, "spec_id", "spec_code")

setkey(gene_counts, spec_code)

GetComboCounts <- function(gene_counts, spec_no, out_csv) {
    print(date())
    print(paste("spec_no", spec_no))
    comb_res <- gene_counts[, combinations(
        length(unique(spec_code)),
        spec_no,
        unique(spec_code),
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
    fwrite(my_counts, out_csv, append = TRUE, compress = "gzip")
    gc(TRUE)
}

out_csv <- "output/plot_data/pan_core_counts.csv.gz"
number_of_specs <- gene_counts[, length(unique(spec_code))]
lapply(rev(1:(number_of_specs - 1)), function(x)
    GetComboCounts(gene_counts, x, out_csv))

