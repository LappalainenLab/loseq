#!/usr/bin/env Rscript
#-------------------------------------------------------
# Samples passing QC per prep
# Author: Molly Martorella
#-------------------------------------------------------

set.seed(1)

library(tidyverse)
library(data.table)

# Inputs ----------------------

meta_full_path <- "data/samplekey_508.txt"
meta_seq_path <- "data/samplekey_485.txt"

cts_path <- "data/downsampled_counts/loseq_allpreps_threshold1000000_raw.txt" #1e6 threshold

out_path <- "analysis/qc/"
fig_path <- "analysis/qc/figs/samples_passing_qc_" # with script prefix


# Functions ----------------------

source(here::here("src/base_functions.R"))


# Read in data ----------------------

meta_full <- fread(meta_full_path, data.table = FALSE)
meta_seq <- fread(meta_seq_path, data.table = FALSE)

cts <- fread(cts_path, data.table = FALSE)

source(here::here("src/gtex_loseq_colors.R"))
prep_colors3 <- c(prep_colors, "seagreen3")
names(prep_colors3) <- c("illumina", "loseq_replicate_1", "smartseq",  "loseq_replicate_2") # named character vector, names are preps, values are colors
legendord <- prep_colors3[c("loseq_replicate_1", "loseq_replicate_2", "illumina", "smartseq")]


# Clean data ----------------------

# set up replicate naming, remove unreplicated loseq, remove HEK, add pass indicator
meta_reps <- meta_full %>%
  mutate(prep = ifelse(prep == "loseq" & replicate %in% c(0,1), "loseq_replicate_1", prep),
         prep = ifelse(prep == "loseq" & replicate == 2, "loseq_replicate_2", prep)) %>%
  filter(prep != "loseq", # does not remove unreplicated loseq, not a totally necessary line
         replicate != 3) %>%  # removes 3rd HEK replicate
  filter(tissue != "hek") %>%
  select(group_name, donor_id, tissue, prep, collection_number) %>%
  mutate(pass = ifelse(group_name %in% colnames(cts), "", "X"))

# add counts column for stacked barplot plot
temp <- meta_reps %>%
  group_by(donor_id, tissue, prep, collection_number) %>%
  summarise(n = n()) %>%
  mutate(prep = factor(prep, levels = c("smartseq", "illumina", "loseq_replicate_2", "loseq_replicate_1")))

meta_reps <- merge(temp, meta_reps, by = c("donor_id", "tissue", "prep", "collection_number"))

# put collection number/prep in correct order
meta_reps <- meta_reps %>%
  select(donor_id, tissue, prep, collection_number, n, pass) %>%
  group_by(donor_id, tissue) %>%
  arrange(desc(prep), .by_group = TRUE)



# Stacked barplot with prep, coll#, and pass/fail indicator ----------------------

meta_reps %>%
  ggplot(aes(x = donor_id, y = n)) +
  geom_col(aes(fill = prep), colour = 'black') +
  geom_text(aes(label = as.factor(collection_number)), position = position_stack(vjust = 0.5), color = "black") +
  geom_text(aes(label = pass), position = position_stack(vjust = 0.5), color = "white", alpha = 0.75, lwd = 5) +
  scale_y_continuous(breaks = seq(0,8,1), limits = c(0,8)) +
  facet_grid(~tissue) +
  coord_flip() +
  # scale_pattern_type_discrete(choices = c('vertical', 'horizontal', 'circles', "hs_cross")) +
  scale_fill_manual(name = "Prep", values = legendord, labels = c("loseq rep1", "loseq rep2", "illumina", "smartseq"), drop = TRUE) +
  ggtitle("Samples per donor",
          subtitle = "X's indicate failed QC, Numbers indicate collection #") +
  xlab("Donor") +
  ylab("Samples") +
  theme_bw() +
  scale_x_discrete(labels = c(seq(1,19,1))) +
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        plot.subtitle = element_text(size = 13),
        axis.title = element_text(size = 13),
        axis.title.y = element_text(vjust = 2),
        axis.title.x = element_text(vjust = -0.75),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())

ggsave(filename = paste0(fig_path, "per_donor_grid.png"), height = 4.5, width = 8)

write.table(meta_reps, file = "analysis/final_figs/inputs/fig1_MAIN_per_donor_grid.txt", sep = "\t")

# Proportion samples passing ----------------------

meta_full %>%
  mutate(passed = recode(passed, cross_prep_QC = "seq_QC", loseq_QC = "seq_QC", wetlab_prep = "wetlab_QC")) %>%
  mutate(passed = fct_relevel(passed, c("failed", "wetlab_QC", "seq_QC"))) %>%
  group_by(tissue, passed, prep) %>%
  summarise(n = n()) %>%
  mutate(proportion = case_when(

    prep == "loseq" & tissue != "hek" ~ n/90,
    prep == "loseq" & tissue == "hek" ~ n/24,
    prep == "illumina" & tissue != "hek" ~ n/16,
    prep == "illumina" & tissue == "hek" ~ n/6,
    prep == "smartseq" & tissue != "hek" ~ n/12,
    prep == "smartseq" & tissue == "hek" ~ n/6

  )) %>%
  ggplot(aes(x = prep, y = proportion, fill = passed)) +
  geom_bar(position = "stack", stat = "identity") +
  geom_text(aes(label = round(proportion, digits = 2)), position = position_stack(vjust = 0.5), color = "white") +
  facet_wrap(~tissue, nrow = 1) +
  ylab("Proportion") +
  #xlab("Library Prep") +
  ggtitle("Sample outcomes following each QC step") +
  # scale_fill_manual(values = prep_colors) +
  # scale_alpha("passed") +
  scale_fill_manual(labels = c("Failed Prep QC", "Failed Seq QC", "Passed"),
                    values = c("tomato3", "orange", "green4")) + # "orangered", "orangered3", "orangered4"
  labs(fill = "Outcome") +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        plot.subtitle = element_text(size = 13),
        axis.title = element_text(size = 13),
        axis.title.y = element_text(vjust = 2),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 12),
        strip.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.position = "top",
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey92", colour = NA))

ggsave(file = paste0(fig_path, "proportion_passing.png"), width = 6.5, height = 3.5)

write.table(meta_full, file = "analysis/final_figs/inputs/fig1_SUPP_prop_pass.txt", sep = "\t")


# FOR AI OUTPUT

meta_full %>%
  mutate(passed = recode(passed, cross_prep_QC = "seq_QC", loseq_QC = "seq_QC", wetlab_prep = "wetlab_QC")) %>%
  mutate(passed = fct_relevel(passed, c("failed", "wetlab_QC", "seq_QC"))) %>%
  group_by(tissue, passed, prep) %>%
  summarise(n = n()) %>%
  mutate(proportion = case_when(

    prep == "loseq" & tissue != "hek" ~ n/90,
    prep == "loseq" & tissue == "hek" ~ n/24,
    prep == "illumina" & tissue != "hek" ~ n/16,
    prep == "illumina" & tissue == "hek" ~ n/6,
    prep == "smartseq" & tissue != "hek" ~ n/12,
    prep == "smartseq" & tissue == "hek" ~ n/6

  )) %>%
  ggplot(aes(x = prep, y = proportion, fill = passed)) +
  geom_bar(position = "stack", stat = "identity") +
  geom_text(aes(label = round(proportion, digits = 2)), position = position_stack(vjust = 0.5), color = "white") +
  facet_wrap(~tissue, nrow = 1) +
  ylab("Proportion") +
  #xlab("Library Prep") +
  #ggtitle("Sample outcomes following each QC step") +
  # scale_fill_manual(values = prep_colors) +
  # scale_alpha("passed") +
  scale_fill_manual(labels = c("Failed Prep QC", "Failed Seq QC", "Passed"),
                    values = c("tomato3", "orange", "green4")) + # "orangered", "orangered3", "orangered4"
  labs(fill = "Outcome:") +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  theme_bw() +
  theme(plot.title = element_text(size = 10),
        plot.subtitle = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.title.y = element_text(vjust = 2),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 6),
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 7),
        strip.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.position = "bottom",
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA),
        strip.text.x = element_text(margin = margin(0.34,0,0.34,0, "cm")))



ggsave(file = paste0(fig_path, "proportion_passing_ai.pdf"), width = 5.5, height = 5, useDingbats = FALSE)




