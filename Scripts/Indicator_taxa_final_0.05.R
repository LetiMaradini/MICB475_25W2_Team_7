############################################################
# Load packages
############################################################

library(phyloseq)
library(dplyr)
library(tibble)
library(stringr)
library(indicspecies)

############################################################
# STEP 1 — Load phyloseq object and extract tables
############################################################

load("micro_rare.RData")  # loads micro_rare

# Extract OTU table (taxa as rows, samples as columns)
taxa_table <- as.data.frame(otu_table(micro_rare))
if (!taxa_are_rows(micro_rare)) {
  taxa_table <- t(taxa_table)
}

# Extract metadata
metadata <- as.data.frame(sample_data(micro_rare))

############################################################
# STEP 2 — Run ISA
############################################################

indicator_multipatt <- multipatt(t(taxa_table), metadata$Condition, duleg = TRUE)

### Look at output
summary(indicator_multipatt)

### Write raw summary to file
# Note: if grouping by higher-order taxa, names may wrap across lines —
#   zoom out in RStudio to keep each row on one line.
indicator_output <- capture.output(summary(indicator_multipatt, indvalcomp = TRUE))
write.table(indicator_output, file = "indicator_values.txt", row.names = FALSE, quote = FALSE)

############################################################
# STEP 3 — Extract ISA results
############################################################

isa_results <- as.data.frame(indicator_multipatt$sign)
isa_results <- rownames_to_column(isa_results, "FeatureID")

############################################################
# STEP 4 — Extract taxonomy
############################################################

taxonomy <- as.data.frame(tax_table(micro_rare))
taxonomy <- rownames_to_column(taxonomy, "FeatureID")

############################################################
# STEP 5 — Merge ISA + taxonomy
############################################################

isa_tax <- left_join(isa_results, taxonomy, by = "FeatureID")

############################################################
# STEP 6 — Detect group columns (multipatt format: s.groupname)
############################################################

group_cols <- grep("^s\\.", colnames(isa_tax), value = TRUE)

############################################################
# STEP 7 — Assign group labels
############################################################

isa_tax$Group <- apply(
  isa_tax[, group_cols],
  1,
  function(x) {
    grp <- group_cols[which(x == 1)]
    if (length(grp) == 0) NA else grp[1]
  }
)

isa_tax$Group <- sub("^s\\.", "", as.character(isa_tax$Group))

############################################################
# STEP 8 — Create Taxon column (last non-NA rank)
############################################################

taxonomy_ranks <- colnames(tax_table(micro_rare))

isa_tax$Taxon <- apply(
  isa_tax[, taxonomy_ranks],
  1,
  function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) NA else tail(x, 1)
  }
)

isa_tax$Taxon <- as.character(isa_tax$Taxon)

############################################################
# STEP 9 — Filter to significant indicators and export
############################################################

isa_table <- isa_tax %>%
  filter(p.value <= 0.05) %>%
  select(Taxon, Group, stat, p.value) %>%
  arrange(Group, desc(stat))

write.csv(isa_table, "Indicator_species_results_unfiltered.csv", row.names = FALSE)

############################################################
# STEP 10 — Filter, clean, and export final table
############################################################

final_table <- isa_table %>%
  filter(stat >= 0.7, Group != "Basal Control") %>%
  group_by(Taxon, Group) %>%
  slice_max(stat, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(Group, desc(stat)) %>%
  mutate(
    stat = round(stat, 3),
    Rank = case_when(
      str_starts(Taxon, "g__") ~ "genus",
      str_starts(Taxon, "f__") ~ "family",
      TRUE ~ NA_character_
    ),
    Taxon = str_remove(Taxon, "^[gf]__")
  ) %>%
  select(Taxon, Rank, stat, p.value, Group)

# Split by group and recombine with section headers
basal_control  <- final_table %>% filter(Group == "Basal Control")  %>% select(-Group)
ground_control <- final_table %>% filter(Group == "Ground Control") %>% select(-Group)
space_flight   <- final_table %>% filter(Group == "Space Flight")   %>% select(-Group)

blank <- data.frame(Taxon = NA, Rank = NA, stat = NA, p.value = NA)

final_df <- bind_rows(
  data.frame(Taxon = "=== Ground Control ===", Rank = NA, stat = NA, p.value = NA),
  ground_control,
  blank,
  data.frame(Taxon = "=== Space Flight ===",   Rank = NA, stat = NA, p.value = NA),
  space_flight
)

write.csv(final_df, "Indicator_species_filtered.csv", row.names = FALSE, na = "")
