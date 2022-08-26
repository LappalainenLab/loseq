#-------------------------------------------------------
# GTEx sampling and downsampling
# Author: Molly Martorella
#------------------------------------------------------

set.seed(1)

library(tidyverse)
library(data.table)


# Functions -----------

source(here::here("src/base_functions.R"))
source(here::here("src/convert_tpms.R"))
source(here::here("src/remove_txn_id.R"))
source(here::here("src/downsampling_fxns.R"))

# Inputs -----------

gtex_cts_path <- "data/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz"
gtex_key_path <- "data/gtex/GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt"

gtex_n <- 75
gtex_downsamp_cts <- 5e6 # change and rerun if want different downsample

out_path <- "data/gtex/"


# Read in data -----------

gtex_cts <- fread(gtex_cts_path, sep = "\t", data.table = FALSE)
gtex_key <- fread(gtex_key_path, data.table = FALSE)

gene_annot <- fread("data/gencode.v26.ensemblv103.gene_annot.txt", data.table = FALSE)
gene_lengths <- setNames(object = gene_annot$gene_length, nm = gene_annot$gene_id)


# Process data -----------

gtex_cts$Name <- remove_txn_id(x = gtex_cts$Name)

# make ensembl ids rownames for gtex
gtex_cts <- gtex_cts %>% column_to_rownames("Name")

# filter for samples in gtex key
gtex_cts <- gtex_cts[, colnames(gtex_cts) %in% gtex_key$SAMPID]
gtex_key <- gtex_key %>% filter(SAMPID %in% colnames(gtex_cts))

# downsample each gtex tissue of interest to 19
low_samp_tis <- gtex_key %>% group_by(SMTSD) %>% summarise(n = n()) %>% filter(n < gtex_n) %>% pull(SMTSD)
gtex_key <- gtex_key %>% filter(SMTSD %notin% low_samp_tis) %>% group_by(SMTSD) %>% sample_n(gtex_n)
gtex_cts <- gtex_cts[, colnames(gtex_cts) %in% gtex_key$SAMPID]


# Downsample GTEx counts -----------

# get binomial fractions
fracs <- get_fracs(sample_depths = colSums(gtex_cts), desired_depths = gtex_downsamp_cts)
# binomial sample counts 5x - takes a few minutes
gtex_cts <- downsample(cts_df = gtex_cts, fractions = fracs, sample_size = 5)

# convert gtex_cts to tpms
gtex_tpms <- convert_tpms_nocols(gtex_cts, lengths = gene_lengths)


# Write out files -----------

fwrite(gtex_cts, file = paste0(out_path, "gtex_downsampled_", gtex_n, "samps_", gtex_downsamp_cts, "_cts.txt"), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

fwrite(gtex_tpms, file = paste0(out_path, "gtex_downsampled_", gtex_n, "samps_", gtex_downsamp_cts, "_tpms.txt"), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)



