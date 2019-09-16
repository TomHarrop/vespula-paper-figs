#!/usr/bin/env Rscript

library(data.table)

# get Dmel annotations
dmel_ids <- fread("data/dmel_fbpp_to_name.csv")
search_result_files <- list.files("output/vespula_p450s/search_results",
                                  pattern = ".txt",
                                  full.names = TRUE)
names(search_result_files) <- sub(".txt", "", basename(search_result_files))

FreadHmmerOutput <- function(x){
    my_empty_dt <- data.table("target_name" = character(),
                      "target_accession" = character(),
                      "query_name" = character(),
                      "query_accession" = character(),
                      "e_value" = numeric())
    my_dt <- fread(cmd = paste0('grep -v "^#" ',
                                x),
                   header = FALSE,
                   sep = ' ',
                   select = 1:5,
                   fill = TRUE,
                   col.names = names(my_empty_dt))
    if(dim(my_dt)[[1]] == 0) {
        return(my_empty_dt)
    }
    my_dt[, e_value := as.numeric(e_value)]
    return(my_dt)
}

search_result_list <- lapply(search_result_files, FreadHmmerOutput)
search_results <- rbindlist(search_result_list)

best_hit_rows <- search_results[, .I[which.min(e_value)], by = query_name][, V1]
best_hits <- search_results[best_hit_rows]
hits_with_name <- merge(best_hits,
      dmel_ids,
      by.x = "target_name",
      by.y = "ID",
      all.x = TRUE,
      all.y = FALSE)

wasp_specs <- c("Vger" = "Vespula germanica",
                "Vpen" = "Vespula pensylvanica",
                "Vvul" = "Vespula vulgaris")

hits_with_name[, spec_code := substr(query_name, 1, 4)]

vespula_p450s <- hits_with_name[, .(
    annotation = plyr::revalue(spec_code, wasp_specs),
    accession = query_name,
    best_Dmel_hit = target_name,
    best_Dmel_hit_name = Name,
    phmmer_e_value = e_value
)]
setorder(vespula_p450s, annotation, accession)
fwrite(vespula_p450s, "output/vespula_p450s/vespula_p450s.csv")

