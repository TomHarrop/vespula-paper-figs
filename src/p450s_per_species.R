#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(ggtree)

ReadAlnNames <- function(x) {
    fread(cmd = paste("grep '^>'",
                      x,
                      "| cut -f1 | sed 's/>//g'"),
          fill = TRUE,
          sep = ',',
          header = FALSE,
          col.names = "fasta_name")
}


wasp_specs <- c("Pdom", "Vger", "Vpen", "Vvul")

# species tree
species_data <- fread("data/hymenoptera_genomes.txt",
                      header = FALSE,
                      col.names = c("genus", "species"))
species_data[, spec_name := paste(genus, species)]
species_data[, spec_code := paste0(toupper(substring(genus, 1, 1)),
                                   tolower(substring(species, 1, 3)))]

y <- read.newick("data/trees/wg.txt")

# join species names
tip_dt <- data.table(tip_label = y$tip.label)
tip_dt[, spec_code := substr(tip_label, 1, 4)]
new_labels <- merge(tip_dt,
                    species_data[, .(spec_code, spec_name, species, genus)],
                    all.x = TRUE)
setcolorder(new_labels, neworder = c("tip_label"))
new_labels[genus == "Vespula", is_vespula := TRUE]
new_labels[genus != "Vespula", is_vespula := FALSE]


# annotate with number of P450s
aln_files <- list.files("data/MultipleSequenceAlignments",
                        full.names = TRUE,
                        pattern = ".fa")
names(aln_files) <- sub(".fa", "", basename(aln_files))
aln_list <- lapply(aln_files, ReadAlnNames)
aln_names <- rbindlist(aln_list, idcol = "orthogroup", fill = TRUE)
aln_names[, spec_code := substr(fasta_name, 1, 4)]

# mung the labels
aln_names[!spec_code %in% wasp_specs, gene_name := fasta_name]
aln_names[, gene_name := sub("Pdomi_r3_proteins_", "", gene_name)]
aln_names[, gene_name := sub("-mRNA-1-PROTEIN", "", gene_name)]
aln_names[, gene_name := gsub("\\|.*", "", gene_name)]

# try to remove isoforms
isoform_list <- fread("data/annotations/all_isoforms.csv")
to_keep <- isoform_list[, unique(refseq_protein)]
annotated_species <- isoform_list[, unique(spec_code)]
aln_names[!is.na(gene_name), extracted_prot_id := gsub(
    paste0(spec_code, "_?"),
    "",
    gene_name),
    by = gene_name]
drop1 <- aln_names[(!spec_code %in% c(wasp_specs, "Dmel")) &
              (spec_code %in% annotated_species) &
              (!extracted_prot_id %in% to_keep),
          unique(fasta_name)]


# get Dmel annotations
dmel_ids <- fread("data/dmel_fbpp_to_name.csv")

# dmel isoforms
dmel_droptable <- merge(aln_names[spec_code == "Dmel"],
      dmel_ids,
      by.x = "extracted_prot_id",
      by.y = "ID",
      all.x = TRUE,
      all.y = FALSE)
setorder(dmel_droptable, Name)           
drop2 <- dmel_droptable[duplicated(dmel_droptable, by = "Name"), unique(fasta_name)]

# filter the table
filtered_p450s <- aln_names[!fasta_name %in% c(drop1, drop2)]

# make a note of species we couldn't filter
filtered_p450s[spec_code %in% c(annotated_species, wasp_specs),
               spec_filtered := TRUE]
filtered_p450s[!spec_code %in% c(annotated_species, wasp_specs),
               spec_filtered := FALSE]

p450_counts <- filtered_p450s[, .(p450_count = length(unique(fasta_name)),
                                  spec_filtered = unique(spec_filtered)),
                         by = spec_code]

# make an annotation table
annot_table <- merge(new_labels, p450_counts, by = "spec_code")
annot_table[, pm_label := paste0(
    'italic("',
    spec_name,
    '")~"(',
    p450_count,
    ifelse(spec_filtered, "", "*"),
    ')"'
)]
# don't colour the un-filtered ones
annot_table[spec_filtered == FALSE, p450_count := NA]
setcolorder(annot_table, neworder = c("tip_label"))

gt <- ggtree(y)

gt + geom_nodelab(aes(label = node))

MyCladelab <- function(my_lab, my_node) {
    geom_cladelabel(
        node = my_node,
        label = my_lab,
        offset = 0.55,
        offset.text = 0.01,
        align = TRUE,
        angle = -90,
        hjust = 0.5,
        fontsize = 2
    )
}

p450_tree <- gt %<+% annot_table +
    theme(legend.position = c(1/4, 3/4),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.direction = "horizontal",
          legend.key.size = unit(0.75, "lines")) +
    geom_tree(aes(colour = p450_count),
              size = 1) +
    scale_colour_viridis_c(guide = guide_colourbar(
        title = "P450 genes",
        title.position = "top",
        title.hjust = 0.5,
        title.theme = element_text(size = 8),
        label.theme = element_text(size = 6)
    )) +
    geom_tiplab(aes(label = pm_label,
                    colour = p450_count),
                size = 3,
                parse = TRUE) + 
    geom_treescale(fontsize = 3,
                   linesize = 1,
                   x = 1.25,
                   y = 0) +
    MyCladelab("Formicidae", 35) +
    MyCladelab("Anthophila", 29) +
    MyCladelab("Vespidae", 45)

saveRDS(p450_tree, "output/plot_data/p450_tree.Rds")

# ggsave("figs/p450_tree.pdf",
#        p450_tree,
#        device = cairo_pdf,
#        width = 10,
#        height = 7.5,
#        units = "in")

