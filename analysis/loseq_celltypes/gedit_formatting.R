#!/usr/bin/env Rscript
#-------------------------------------------------------
# GEDIT cts processing - raw cts, rownames should be gene names
# Author: Molly Martorella
#-------------------------------------------------------

set.seed(1)

library(tidyverse)
library(data.table)


# Inputs ----------------------

#cts_path <- "data/downsampled_counts/loseq_allpreps_threshold5000000_downsampled5000000_raw.txt"
cts_path <- "data/downsampled_counts/loseq_allpreps_threshold2500000_raw.txt"

annot_path <- "data/gencode.v33.annotation.GRCh38.gene_id_name_type.txt"
meta_path <- "data/samplekey_485.txt"

out_path <- "analysis/loseq_celltypes/inputs/"
#out_name <- "loseq5mil.txt"
out_name <- "loseq2_5mil_thresh.txt"

# Functions ----------------------

source(here("src/remove_txn_id.R"))


# Read in data ----------------------

cts <- fread(cts_path, data.table = FALSE) %>% rename("gene_id" = "V1")
annot <- fread(annot_path, data.table = FALSE) %>% select("gene_id", "gene_name")
meta <- fread(meta_path, data.table = FALSE)


# Select samples/filter hek ----------------------

samps <- meta %>% filter(group_name %in% colnames(cts),
                         tissue != "hek") %>%
  pull(group_name)

cts <- cts[, colnames(cts) %in% c("gene_id", samps)]


# Format ----------------------

cts <- merge(cts, annot, by = "gene_id")

cts <- cts[!duplicated(cts$gene_name),]
rownames(cts) <- cts$gene_name

cts <- cts %>% select(-gene_id, -gene_name)


# Write to file ----------------------

fwrite(cts, file = paste0(out_path, out_name), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

