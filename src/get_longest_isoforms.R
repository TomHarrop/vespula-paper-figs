#!/usr/bin/env Rscript

library(data.table)
library(rtracklayer)

# gff_file <- "data/annotations/Amel_v4.5_refseq_withRefSeqChrIds.gff3.gz"
# isoform_file <- ""
# feature_types <- c("mRNA")

gff_file <- snakemake@input[[1]]
isoform_file <- snakemake@output[[1]]
feature_types <- c("mRNA")

my_gff <- import.gff(gff_file,
                     feature.type = feature_types)
my_gff_dt <- as.data.table(my_gff)[, .(ID,
                                       Parent = as.character(Parent),
                                       Name,
                                       width,
                                       Dbxref)]
# get the protein id from the DBxref field
GetRefSeqProtein <- function(dbxref){
    flat_xref <- unlist(dbxref)
    refseq_prot <- flat_xref[sapply(flat_xref, function(x)
        grepl("^RefSeq_Prot", x))]
    if(length(refseq_prot) == 0){
        return(as.character(NA))
    }
    return(sub("RefSeq_Prot:", "", refseq_prot))
}
my_gff_dt[, refseq_protein := GetRefSeqProtein(Dbxref), by = ID]
my_gff_dt[, Dbxref := NULL]

setorder(my_gff_dt, -width) # just grabbing the longest
keep_rows <- my_gff_dt[, .I[which.max(width)], by = Parent][, V1]
unique_by_parent <- my_gff_dt[keep_rows]
fwrite(unique_by_parent, isoform_file)
