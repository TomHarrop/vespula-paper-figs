#!/usr/bin/env Rscript

library(data.table)
library(cowplot)
library(ggplot2)
library(ggtree)
library(treeio)

# plot Cyp4e1 and Cyp12 orthogroups, only showing Dmel, Amel and Vespidae genes.
# get trees as list and lapply

# get full species names
species_data <- fread("data/hymenoptera_genomes.txt",
                      header = FALSE,
                      col.names = c("genus", "species"))
species_data[, spec_name := paste(genus, species)]
species_data[, spec_code := paste0(toupper(substring(genus, 1, 1)),
                                   tolower(substring(species, 1, 3)))]

# get Dmel annotations
dmel_ids <- fread("data/dmel_fbpp_to_name.csv")
setkey(dmel_ids, "ID")

# get longest isoforms
longest_isoforms <- fread("data/annotations/all_isoforms.csv")
plottable <- longest_isoforms[, unique(refseq_protein)]

# Cyp4e1 orthogroup
# orthogroup <- "OG0000273"

# Cyp12d1 orthogroup
# orthogroup <- "OG0002213"

wasp_specs <- c("Pdom", "Vger", "Vpen", "Vvul")
specs_to_plot <- c("Amel" , "Dmel", "Vger", "Vpen", "Vvul")

GetTree <- function(orthogroup){
  orthogroup_file <- paste0("data/trees/", orthogroup, "_tree.txt")
  tree_data <- read.newick(orthogroup_file)
  return(tree_data)}
# annot <- data.table(tip_label = tree_data$tip.label,
#                     spec_code = substr(tree_data$tip.label, 1, 4))

GetAnnot <- function(tree_list) {
  annot_list <- lapply(tree_list, function(x)
    data.table(tip_label = x$tip.label,
               spec_code = substr(x$tip.label, 1, 4)))
  annot <- rbindlist(annot_list)
  
  # get the tip labels from the tree data
  annot[spec_code %in% wasp_specs,
        gene_name := gsub("V.{4}_r3_6_", "", tip_label)]
  
  # mung the tip labels
  annot[!spec_code %in% wasp_specs, gene_name := tip_label]
  annot[, gene_name := sub("Pdomi_r3_proteins_", "", gene_name)]
  annot[, gene_name := sub("-mRNA-1-PROTEIN", "", gene_name)]
  annot[, gene_name := gsub("\\|.*", "", gene_name)]
  
  # add the full species names
  annot2 <- merge(annot, species_data, by = "spec_code")
  
  # use Dmel gene names
  annot2[spec_code == "Dmel", ID := sub("Dmel_", "", gene_name)]
  annot3 <- merge(annot2,
                  dmel_ids,
                  by.x = "ID",
                  by.y = "ID",
                  all.x = TRUE,
                  all.y = FALSE)
  # annot3[spec_code == "Dmel", gene_name := gsub("-P[[:upper:]]$", "", Name)]
  annot3[spec_code == "Dmel",
         gene_name := Name]
  
  # tip labels have to be column 1
  setcolorder(annot3, "tip_label")
  return(annot3)}


PlotTree <- function(tree_object, annotation, palette) {
  my_tree <- ggtree(tree_object, layout="circular")
  # add the annotation
  my_annotated_tree <- my_tree %<+% annotation +
    theme(legend.position="right") +
    scale_fill_manual(
      values = palette,
      guide = guide_legend(
        title = NULL,
        # label.hjust = -0.1,
        label.theme = element_text(face = "italic",
                                   size = 6),
        label.position = "bottom"
      )) +
    scale_colour_manual(
      values = palette,) +
    geom_tippoint(aes(fill = spec_name),
                  shape = 21,
                  size = 1.5,
                  stroke = 0.001) +
    geom_tiplab2(aes(label = gene_name,
                     colour = spec_name),
                 size = 1.5,
                 offset = 0.05,
                 show.legend = FALSE)
  
  return(my_annotated_tree)}


tree_list <- lapply(list("OG0000016", "OG0000273", "OG0002213"),
                    GetTree)
tree_annot <- GetAnnot(tree_list)
# try to squish isoforms
tree_annot[, extracted_prot_id := gsub(
  paste0(spec_code, "_?"),
  "",
  gene_name),
  by = gene_name]
drop1 <- tree_annot[(!spec_code %in% c(wasp_specs, "Dmel")) &
                      (!extracted_prot_id %in% plottable),
                    unique(tip_label)]

# ok, but need to deal with Dros (not in refseq)
dmel_droptable <- tree_annot[spec_code == "Dmel", .(
  tip_label,
  Name,
  short_name = gsub("-P[[:upper:]]$", "", Name))]
setorder(dmel_droptable, short_name)           
drop3 <- dmel_droptable[duplicated(dmel_droptable, by = "short_name"),
                        unique(tip_label)]

#drop the non-species of interest
drop2 <- tree_annot[!spec_code %in% specs_to_plot, unique(tip_label)]
pruned_trees <- lapply(tree_list, drop.tip, tip = c(drop1, drop2, drop3))
class(pruned_trees) <- "multiPhylo"

# set up a palette
manual_labels <- tree_annot[spec_code %in% specs_to_plot, unique(spec_name)]
manual_pal <- viridis::viridis_pal()(length(manual_labels))
names(manual_pal) <- manual_labels

tree_plots <- lapply(pruned_trees,
                     PlotTree,
                     annotation = tree_annot,
                     palette = manual_pal)


# save plot data
saveRDS(tree_plots, "output/plot_data/tree_plots.Rds")



##################
# PLOT ALL TREES #
##################

# og_files <- list.files("data/trees", pattern = "OG.*_tree.txt")
# all_ogs <- sub("_tree.txt", "", og_files)
# names(all_ogs) <- all_ogs
# all_trees <- lapply(all_ogs, GetTree)
# all_annot <- GetAnnot(all_trees)
# 
# # prune
# all_annot[, extracted_prot_id := gsub(
#   paste0(spec_code, "_?"),
#   "",
#   gene_name),
#   by = gene_name]
# 
# drop1 <- all_annot[(!spec_code %in% c(wasp_specs, "Dmel")) &
#                      (!extracted_prot_id %in% plottable),
#                    unique(tip_label)]
# 
# # ok, but need to deal with Dros (not in refseq)
# dmel_droptable <- all_annot[spec_code == "Dmel", .(
#   tip_label,
#   Name,
#   short_name = gsub("-P[[:upper:]]$", "", Name))]
# setorder(dmel_droptable, short_name)           
# drop3 <- dmel_droptable[duplicated(dmel_droptable, by = "short_name"),
#                         unique(tip_label)]
# 
# #drop the non-species of interest
# drop2 <- all_annot[!spec_code %in% specs_to_plot, unique(tip_label)]
# 
# all_pruned <- lapply(all_trees, drop.tip, tip = c(drop1, drop2, drop3))
# class(all_pruned) <- "multiPhylo"
# 
# all_plots <- lapply(all_pruned, PlotTree, annotation = all_annot)
# 
# PrintPlot <- function(plot, plot_name) {
#   plot_file <- paste0("output/trees/", plot_name, ".pdf")
#   ggsave(plot_file,
#          plot,
#          device = cairo_pdf,
#          width = 10, 
#          height = 7.5,
#          units = "in")
# }
# 
# lapply(names(all_plots), function(x)
#   PrintPlot(all_plots[[x]], x))
# 
# 
# print(all_trees[["OG0000273"]] + xlim(c(0, 2.5))) 
