library(data.table)

isoform_files <- list.files("data/annotations",
                            pattern = "_longest_isoforms.csv",
                            full.names = TRUE)
names(isoform_files) <- sapply(basename(isoform_files),
                               substr,
                               start = 1,
                               stop = 4)
all_isoforms <- rbindlist(lapply(isoform_files, fread), idcol = "spec_code")
fwrite(all_isoforms, "data/annotations/all_isoforms.csv")
