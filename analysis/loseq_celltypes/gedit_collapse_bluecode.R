#!/usr/bin/env Rscript
#-------------------------------------------------------
# Collapse cell types - BlueCode
# Author: Molly Martorella
#-------------------------------------------------------

set.seed(1)

library(tidyverse)
library(data.table)


# Inputs ----------------------

#res_path <- "analysis/loseq_celltypes/results/loseq5mil_bluecode.tsv"
res_path <- "analysis/loseq_celltypes/results/loseq2_5mil_thresh_bluecode.tsv"

out_path <- "analysis/loseq_celltypes/results/"
out_name <- paste0(gsub(basename(res_path), pattern = ".tsv", replacement = ""), "_collapsed_cells.txt")


# Functions ----------------------

source(here("src/remove_txn_id.R"))


# Read in data ----------------------

res <- fread(res_path, data.table = FALSE) %>% dplyr::rename("group_name" = "V1",
                                                             "T.cell.CD4" = "T.Cell.CD4")

# Fix group_names ----------------------

res$group_name <- gsub(res$group_name, pattern = "\\.", replacement = "-")

# Cell type collapse categories ----------------------

cell_types <- c("Chondrocyte", "B.cell", "Monocyte", "NK.cell", "Endothelial", "Epithelial", "Fibroblast", "keratinocyte", "Macrophage", "neutrophil", "Melanocyte", "Smooth.muscle", "T.cell", "Osteoblast", "MSC.like.pluripotent.cell")

out_df <- data.frame("group_name" = res[, c("group_name")]) # set up results df

for (i in 1:length(cell_types)) {

  cols <- as.data.frame(res[,str_detect(colnames(res), pattern = cell_types[i])]) # columns matching cell type pattern

  if (ncol(cols) == 1) {
    out_df[, cell_types[i]] <- cols[, 1] # otherwise get error from rowSums
  } else {
    out_df[, cell_types[i]] <- rowSums(cols)
  }

}

out_df[, "Myocyte"] <- rowSums(res[, c("Myometrial.cell", "Myocyte.regular.cardiac")])

out_df <- out_df %>% mutate("Myocyte" = Smooth.muscle + Myocyte,
                            "Epithelial" = Epithelial + keratinocyte) %>% select(-Smooth.muscle, -keratinocyte) %>%
  dplyr::rename("MSC" = "MSC.like.pluripotent.cell")


# Write to file ----------------------

fwrite(out_df, file = paste0(out_path, out_name), sep = "\t", col.names = TRUE)

















