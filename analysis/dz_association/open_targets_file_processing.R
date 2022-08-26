#!/usr/bin/env Rscript
#-------------------------------------------------------
# Open Targets File Processing
# Author: Molly Martorella
#-------------------------------------------------------

set.seed(1)

library(dplyr, lib.loc = "/gpfs/commons/home/mmartorella/R/x86_64-pc-linux-gnu-library/3.6/")
library(sparklyr, lib.loc = "/gpfs/commons/home/mmartorella/R/x86_64-pc-linux-gnu-library/3.6/")
library(sparklyr.nested, lib.loc = "/gpfs/commons/home/mmartorella/R/x86_64-pc-linux-gnu-library/3.6/")

# Inputs -----------

dzPath <- "data/open_targets/diseases"

assocPath <- "data/open_targets/associationByOverallDirect"

out_path <- "data/open_targets/220209_20dz_genetargets_"

# select diseases of interest - 20 dz
id_list <- c("EFO_0003086", "EFO_0003818", "EFO_0000270", "EFO_0000537", "EFO_0003882", "HP_0003124", "EFO_0005856", "EFO_0003914", "EFO_0000712", "EFO_0000400", "EFO_0000249", "EFO_0004247", "EFO_0000253", "EFO_0003885", "EFO_0002690", "EFO_0003767", "EFO_0000676", "EFO_0000274", "EFO_0001421", "EFO_0003777")


# Establish connection -----------

sc <- spark_connect(master = "local")


# Process disease data -----------

dz <- spark_read_parquet(sc, path = dzPath) # read dz file

schemacolumns <- dz %>%  sdf_schema() %>%
  lapply(function(x) do.call(tibble, x)) %>%
  bind_rows() # View the disease columns

dzSelect <- dz %>%
  select(id, description, name, therapeuticAreas, ontology) %>%
  sdf_explode(therapeuticAreas) # select columns of interest

dz_df <- dzSelect %>% collect() # convert to tibble


# Process association data -----------

assoc <- spark_read_parquet(sc, path = assocPath) # read association file

#schemacolumns <- assoc %>% sdf_schema() %>% lapply(function(x) do.call(tibble, x)) %>% bind_rows() # if want to check columns

assoc_df <- assoc %>% select(diseaseId, targetId, score, evidenceCount) %>% collect() # select columns, convert to tibble


# Select dzs, merge dfs -----------

dz_df <- dz_df %>% filter(id %in% id_list) %>% select(id, name) %>% distinct() %>% rename("diseaseId" = "id")

# all evidence
assoc_df <- assoc_df %>% filter(diseaseId %in% id_list)

# merge
df_out <- merge(dz_df, assoc_df)


# Write out to files -----------

write.table(df_out, file = out, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)






