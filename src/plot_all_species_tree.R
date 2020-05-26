#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(ggtree)
library(treeio)

ReadAlnNames <- function(x) {
    fread(cmd = paste("grep '^>'",
                      x,
                      "| cut -f1 | sed 's/>//g'"),
          fill = TRUE,
          sep = ',',
          header = FALSE,
          col.names = "fasta_name")
}


wasp_specs <- c("Pdom", "Pcan", "Vger", "Vpen", "Vvul")

# species tree
species_data <- fread("data/hymenoptera_genomes.txt",
                      header = FALSE,
                      col.names = c("genus", "species"))
species_data[, spec_name := paste(genus, species)]
species_data[, spec_code := paste0(toupper(substring(genus, 1, 1)),
                                   tolower(substring(species, 1, 3)))]

# read the tree
tree_data <- read.newick("data/SpeciesTree_rooted.txt")

# join species names
tip_dt <- data.table(tip_label = tree_data$tip.label)
tip_dt[, spec_code := gsub("^([[:alpha:]])[[:alpha:]]+_([[:alpha:]]{3}).*", "\\1\\2", tip_label)]
new_labels <- merge(tip_dt,
                    species_data[, .(spec_code, spec_name, species, genus)],
                    all.x = TRUE)
setcolorder(new_labels, neworder = c("tip_label"))
new_labels[genus == "Vespula", is_vespula := TRUE]
new_labels[genus != "Vespula", is_vespula := FALSE]
setcolorder(new_labels, neworder = c("tip_label"))

# make the tree
gt <- ggtree(tree_data)
spec_tree <- gt %<+% new_labels

saveRDS(spec_tree, "output/plot_data/all_species_tree.Rds")
