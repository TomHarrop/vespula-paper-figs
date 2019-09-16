#!/usr/bin/env Rscript

library(data.table)

ReadAlnNames <- function(x) {
    fread(cmd = paste("grep '^>'",
                      x,
                      "| cut -f1 | sed 's/>//g'"),
          fill = TRUE,
          sep = ',',
          header = FALSE,
          col.names = "gene_name")
}

# read the gene names out of the alignment files
aln_files <- list.files("data/MultipleSequenceAlignments",
                        full.names = TRUE,
                        pattern = ".fa")
names(aln_files) <- sub(".fa", "", basename(aln_files))
aln_list <- lapply(aln_files, ReadAlnNames)
aln_names <- rbindlist(aln_list, idcol = "orthogroup", fill = TRUE)
aln_names[, spec_code := substr(gene_name, 1, 4)]

# identify our genomes
wasp_prefixes <- c("Vpens", "Vgerm", "Vvulg")
aln_names[, wasp_hit := any(sapply(wasp_prefixes, function(y)
    grepl(y, gene_name))), by = gene_name]

# count P450s
p450_counts <- aln_names[, .(p450_count = length(unique(gene_name))),
                         by = .(orthogroup, spec_code)]

# count backround
total_genes <- aln_names[, length(unique(gene_name))]
total_wasp_genes <- aln_names[wasp_hit == TRUE, length(unique(gene_name))]

# look for enrichment / depletion
wasp_counts <- merge(aln_names[, .(orthogroup_size = .N), by = orthogroup],
                     aln_names[wasp_hit == TRUE, .(wasp_members = .N), by = orthogroup])

wasp_counts[, wasp_enrichment := (wasp_members/orthogroup_size) / (total_wasp_genes/total_genes)]

wasp_counts[wasp_enrichment > 1 &
                wasp_members > 10,
            p_enrich := phyper(
                q = wasp_members - 1,
                m = total_wasp_genes,
                n = total_genes - total_wasp_genes,
                k = orthogroup_size,
                lower.tail = FALSE
            )]
wasp_counts[wasp_enrichment < 0.5 &
                wasp_members <= 3,
            p_deplete := phyper(
                q = wasp_members,
                m = total_wasp_genes,
                n = total_genes - total_wasp_genes,
                k = orthogroup_size,
                lower.tail = TRUE
            )]

LowestP <- function(x){
    my_c <- x[!is.na(x)]
    if(length(my_c) == 0) {
        return(as.numeric(NA))
    }
    return(min(my_c))
}

wasp_counts[, p_lowest := LowestP(c(p_enrich, p_deplete)),
            by = orthogroup]
wasp_counts[!is.na(p_lowest), p_lowest_adj := p.adjust(p_lowest, "fdr")]

setorder(wasp_counts, p_lowest_adj, na.last = TRUE)
wasp_counts
