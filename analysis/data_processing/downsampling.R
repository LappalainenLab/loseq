# Thresholding and Downsampling
# Molly Martorella
# 10/19/21

# Iteratively threshold and downsample counts files
# Output sample number remaining per group and downsampled, thresholded file

library(tidyverse)

set.seed(1)



#### INPUTS ####

cts_df <- data.table::fread("sequencing/200914_novaseq_processing/featureCounts/loseq485_raw_nofilt_featcounts.txt", data.table = FALSE)
key_df <- read.table("data/samplekey_508.txt", sep = "\t", header = TRUE)

sample_id_col <- "group_name"
depths <- c(1e6, 2.5e6, 5e6, 7.5e6, 10e6) # also used as threshold
groups <- c("prep", "tissue") # colnames of key_df

out_dir <- "data/downsampled_counts/" # path for output, make fig path dir.create(file.path(mainDir, subDir))
out_file_prefix <- "loseq_allpreps_"

# optional changes:
protein_coding <- read.table("data/ensembl_gene_type_grch38_v103.txt", sep = "\t", header = TRUE) # optional file with gene annotations, used to filter for protein coding/lncRNA before thresholding.
gene_col <- "Geneid" # optional name of gene name column
length_col <- "Length" # optional
output_threshold_file <- TRUE # boolean, optional, set to FALSE by default
output_norm_counts <- TRUE # boolean, optional, set to FALSE by default; REQUIRES length column in input
samp_size=5

#### CONSTANTS ####

samples <- key_df %>% select(all_of(sample_id_col)) %>% pull()
samples <- as.character(samples)

# gene vector for subsetting, so calculate depth using this subset and that is used for binomial success prob
if (length(protein_coding) != 0) {
  genes <- protein_coding %>% filter(type == "protein_coding" | type == "lncRNA") %>% pull(gene_id)
  genes <- as.character(genes)
} else if (length(gene_col) != 0) {
  genes <- cts_df %>% pull(all_of(gene_col))
  genes <- as.character(genes)
} else {
  genes <- rownames(cts_df)
}

if (length(length_col) != 0 & length(gene_col) != 0) {
  length <- c(cts_df[, length_col])
  names(length) <- c(cts_df[, gene_col])
} else {
  print("No length_col and/or gene_col provided.")
}



#### FUNCTIONS ####

# outputs fraction of total reads (binomial probability) to achieve target depth per sample
get_fracs <- function(sample_depths, desired_depths) {

  temp <- lapply(sample_depths, function(x) {

    ifelse(desired_depths/x >= 1, yes = 1, no =  desired_depths/x)

  })

}


# convert to tpms
convert_tpms <- function(r, geneid_colname, length_colname) {
  temp <- r %>% select(-geneid_colname, -length_colname)
  mat <- as.matrix(temp)
  len <- r %>% select(all_of(length_colname))
  len <- unlist(len)

  # normalize for gene length
  mat <- mat/len

  # normalize for sequencing depth
  depth <- colSums(mat)
  scaled <- depth/(1e6)
  tpms <- as.data.frame(t(t(mat)/scaled))

  # finalize df/add back geneids
  temp <- r %>% select(all_of(geneid_colname))
  tpms <- cbind(temp, tpms)

  return(tpms)
}


# convert to tpms df without cols
convert_tpms_nocols <- function(r, geneids, lengths) {
  mat <- as.matrix(r)

  # normalize for gene length
  mat <- mat/lengths

  # normalize for sequencing depth
  depth <- colSums(mat)
  scaled <- depth/(1e6)
  tpms <- as.data.frame(t(t(mat)/scaled))

  # finalize df/add back geneids
  rownames(tpms) <- geneids

  return(tpms)
}


# downsample
downsample <- function(df_list, fractions, sample_size = 5) {

  df_down_list <- vector(mode = "list", length = length(df_list))

  for (i in 1:length(df_list)) {

    df_down <- data.frame(row.names = rownames(df_list[[i]]))

    for (k in 1:ncol(df_list[[i]])) {

      df_down[,as.character(names(df_list[[i]])[k])] <- rowMeans(as.data.frame(replicate(sample_size, rbinom(n = nrow(df_list[[i]]), size = df_list[[i]][,k], prob = unlist(fractions[[i]][k])), simplify = FALSE)))

    }

    df_down_list[[i]] <- df_down

  }

  return(df_down_list)

}





#### Main ####

#PLAN
# - get vector of protein coding depths per sample
# - convert total raw counts to tpms/cpms
# - implement sample filter based on protein coding depth minimum
# - output thresholded raw, cpm, and tpm files
# - calculate binomial sampling success rate using total depth and only samples that passed protein depth threshold
# - downsample original raw counts to different depths - iterate over thresholded raw counts
# - normalize downsampled raw counts
# - output thresholded, downsamples raw and norm counts



# NORMALIZE COUNTS

if (output_norm_counts == TRUE & length(length_col) != 0 & length(gene_col) != 0) {
  cts_tpms <- convert_tpms(r = cts_df, geneid_colname = gene_col, length_colname = length_col)
  cts_cpms <- as.data.frame(edgeR::cpm(cts_df[, colnames(cts_df) %in% samples]))
} else {
  cts_tpms <- NULL
  cts_cpms <- NULL
  print("Skipping norm counts. Either no gene_col or length_col was provided, or output_norm_counts set to false.")
}



# FORMAT

if (length(gene_col) != 0){
  cts_df <- cts_df %>% column_to_rownames(gene_col)
} else {
  print("Assuming rownames are gene ids.")
}

if (length(cts_tpms) != 0 & length(cts_cpms) != 0){
  cts_tpms <- cts_tpms %>% column_to_rownames(gene_col)
  rownames(cts_cpms) <- rownames(cts_df)
} else {
  print("Output_norm_counts set to false or missing necessary columns.")
}

cts_df <- cts_df[, colnames(cts_df) %in% samples] # removes length col

df_list <- list(cts_df, cts_cpms, cts_tpms)
df_list <- df_list[lengths(df_list) != 0]



# df_list <- lapply(df_list, function(x) { # filters for protein coding/lncRNA genes, TPM columns will not add to 1e6 because of this but that should be fine.
#   x[rownames(x) %in% genes,]
# })


# THRESHOLD FOR PASSING SAMPLES

thresh_samp_list <- lapply(depths, function(x) { # generates list of samples passing p/lnc thresholds
  temp <- colSums(cts_df[rownames(cts_df) %in% genes,]) >= x # filters based on protein coding depth if gene list is protein coding genes (if read in file)
  filt <- cts_df[, temp]
  return(colnames(filt))
})

options(scipen = 999)
names(thresh_samp_list) <- as.character(depths)


df_filts <- lapply(df_list, function(x) {
  lapply(thresh_samp_list, function(y) {x[,colnames(x) %in% y]})
})


if (output_norm_counts == TRUE) {
  names(df_filts) <- c("raw", "cpms", "tpms")
} else {
  names(df_filts) <- c("raw")
}



# WRITE OUT THRESHOLD FILES

for (i in 1:length(df_filts)) {
  for (k in 1:length(df_filts[[i]])) {

    data.table::fwrite(df_filts[[i]][[k]], file = paste0(out_dir, out_file_prefix, "threshold", names(df_filts[[i]][k]), "_", names(df_filts[i]), ".txt"), sep = "\t", row.names = TRUE, col.names = TRUE)

  }

}



# DOWNSAMPLE COUNTS THEN NORMALIZE

# get binomial fractions

fracs <- vector(mode = "list", length = length(df_raw_thresh))

df_raw_thresh <- df_filts[[1]]

for (i in 1:length(df_raw_thresh)) {

  fracs[[i]] <- get_fracs(sample_depths = colSums(df_raw_thresh[[i]]), desired_depths = depths[i])

}

names(fracs) <- as.character(depths)


# binomial sample counts 5x

df_down <- downsample(df_list = df_raw_thresh, fractions = fracs) ## here

options(scipen = 999)
names(df_down) <- as.character(depths)

# s <- df_raw_thresh[[1]][,1]
# num <- nrow(df_raw_thresh[[1]])
# p <- unlist(fracs[[1]][1])
# test <- rbinom(n = num, size = s, prob = p)
# test <- downsample(df_list = df_raw_thresh, fractions = fracs)


# create norm dfs

if (output_norm_counts == TRUE & length(length_col) != 0 & length(gene_col) != 0) {
  down_tpms <- lapply(df_down, function(x) {convert_tpms_nocols(r = x, geneids = names(length), lengths = length)})
  down_cpms <- lapply(df_down, function(x) {edgeR::cpm(x)})
} else {
  down_tpms <- NULL
  down_cpms <- NULL
  print("Skipping downsampled norm counts. Either no gene_col or length_col was provided, or output_norm_counts set to false.")
}


# WRITE OUT DOWNSAMPLED COUNTS

df_down_list <- list(df_down, down_cpms, down_tpms)
df_down_list <- df_down_list[lengths(df_down_list) != 0]


if (output_norm_counts == TRUE) {
  names(df_down_list) <- c("raw", "cpms", "tpms")
} else {
  names(df_down_list) <- c("raw")
}


for (i in 1:length(df_down_list)) {
  for (k in 1:length(df_down_list[[i]])) {

    data.table::fwrite(df_down_list[[i]][[k]], file = paste0(out_dir, out_file_prefix, "threshold", names(df_down_list[[i]][k]), "_downsampled", names(df_down_list[[i]][k]), "_", names(df_down_list[i]), ".txt"), sep = "\t", row.names = TRUE, col.names = TRUE)

  }

}















