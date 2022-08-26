# Creating comprehensive gene reference
# Molly Martorella
# 12/12/21

# Libraries

library(tidyverse)

### INPUTS

gtf <- rtracklayer::readGFF("data/gencode.v26.annotation.collapsed.gtf")
fc <- data.table::fread(file = "sequencing/200914_novaseq_processing/featureCounts/loseq485_raw_nofilt_featcounts.txt", data.table = FALSE) %>% select(Geneid, Length)
gene_types <- read.table("data/ensembl_gene_type_grch38_v103.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

out_path <- "data/"


### PROCESSING

# Filter for genes only

gtf <- gtf %>% filter(type == "gene") %>% rename(annot_type = type)

# reformat gene_ids

gtf$gene_id <- gsub(pattern = "\\.\\d+", replacement = "", gtf$gene_id)

# add gene length

fc <- fc %>% rename(gene_id = Geneid,
                   gene_length = Length)

g_fc <- merge(gtf, fc, by = "gene_id")

ref_out <- left_join(g_fc, gene_types, by = "gene_id")

# add gene length and gc content --> need fucking newer version of R to do this, love NYGC cluster
#info <-  EDASeq::getGeneLengthAndGCContent(gtf$gene_id, org = "hsa") # outputs matrix with ensembl id rownames and columns length, gc
# info <- as.data.frame(info)
# info <- info %>% rownames_to_column("gene_id")

# Fill NAs

ref_out <- ref_out %>% mutate(type = ifelse(is.na(type), gene_type, type))


### WRITE TO FILE

write.table(ref_out, file = paste0(out_path, "gencode.v26.ensemblv103.gene_annot.txt"), sep = "\t", col.names = TRUE, row.names = FALSE)


