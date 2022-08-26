#-------------------------------------------------------
# xCell formatting
# Author: Molly Martorella
#------------------------------------------------------

set.seed(1)

library(tidyverse)
library(data.table)
library(here)


# Functions -----------

source(here::here("src/base_functions.R"))


# Inputs -----------

# Need normalized counts for xCell, and rownames should be hgnc gene symbols
# Use downsampled gtex cts: analyses/data_processing/gtex_sampling_downsampling.R

gtex_cts_path <- "data/gtex/gtex_downsampled_75samps_5e+06_tpms.txt"
gtex_key_path <- "data/gtex/GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt"

loseq_cts_path <- "data/downsampled_counts/loseq_allpreps_threshold5000000_downsampled5000000_tpms.txt"
loseq_key_path <- "data/samplekey_485.txt"

out_path <- "analysis/gtex_comparison/xcell_formatted/"

gene_annot_path <- "data/gencode.v26.ensemblv103.gene_annot.txt"


# Read in data -----------

gtex_cts <- fread(gtex_cts_path, sep = "\t", data.table = FALSE) %>% column_to_rownames("V1")
gtex_key <- fread(gtex_key_path, data.table = FALSE)

loseq_cts <- fread(loseq_cts_path, sep = "\t", data.table = FALSE) %>% column_to_rownames("V1")
loseq_key <- fread(loseq_key_path, data.table = FALSE)

gene_annot <- fread(gene_annot_path, data.table = FALSE)


# Process loseq -----------

loseq_key <- loseq_key %>% filter(group_name %in% colnames(loseq_cts),
                                  tissue != "hek",
                                  prep == "loseq")

loseq_key <- loseq_key[!duplicated(loseq_key$donor_col_tis),]

loseq_key <- loseq_key %>% select(group_name, tissue)

loseq_cts <- loseq_cts[, colnames(loseq_cts) %in% loseq_key$group_name]


# Process gtex -----------

gtex_key <- gtex_key %>% filter(SAMPID %in% colnames(gtex_cts)) %>%
  rename("SAMPID" = "group_name",
         "SMTSD" = "tissue") %>%
  select(group_name, tissue)


# rbind keys -----------

key <- rbind(loseq_key, gtex_key)


# Merge cts -----------

all_tpms <- merge(loseq_cts, gtex_cts, by = 0)
all_tpms <- all_tpms %>% rename("Row.names" = "gene_id")


# Format gene names -----------

all_tpms <- merge(all_tpms, gene_annot[, c("gene_id", "gene_name")], by = c("gene_id"))

# collapse duplicate gene names, take the mean of the counts: TAKES A WHILE - 5MIN OR SO
all_tpms <- all_tpms %>% select(-c("gene_id")) %>% group_by(gene_name) %>% summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))

all_tpms <- all_tpms %>% column_to_rownames("gene_name")


# Write out to file

fwrite(all_tpms, file = file.path(out_path, "xcell_formatted_loseq_gtex_5e6_tpms.txt"), sep = "\t", col.names = TRUE, row.names = TRUE)

fwrite(key, file = file.path(out_path, "xcell_formatted_loseq_gtex_5e6_key.txt"), sep = "\t")






