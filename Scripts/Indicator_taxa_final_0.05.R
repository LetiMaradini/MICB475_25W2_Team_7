############################################################
# Load packages
############################################################

library(phyloseq)
library(dplyr)
library(tibble)

############################################################
# STEP 1 — Extract ISA results
############################################################

isa_results <- as.data.frame(indicator_multipatt_1$sign)
isa_results <- rownames_to_column(isa_results, "FeatureID")

############################################################
# STEP 2 — Extract taxonomy
############################################################

taxonomy <- as.data.frame(tax_table(ps))
taxonomy <- rownames_to_column(taxonomy, "FeatureID")

############################################################
# STEP 3 — Merge ISA + taxonomy
############################################################

isa_tax <- left_join(isa_results, taxonomy, by = "FeatureID")

############################################################
# STEP 4 — Detect group columns (multipatt format)
############################################################

group_cols <- grep("^s\\.", colnames(isa_tax), value = TRUE)

############################################################
# STEP 5 — Assign group labels correctly
############################################################

isa_tax$Group <- apply(
  isa_tax[, group_cols],
  1,
  function(x) {
    grp <- group_cols[which(x == 1)]
    if(length(grp) == 0) NA else grp[1]
  }
)

isa_tax$Group <- as.character(isa_tax$Group)
isa_tax$Group <- sub("^s\\.", "", isa_tax$Group)

############################################################
# STEP 6 — Create Taxon column safely
############################################################

taxonomy_ranks <- colnames(tax_table(ps))

isa_tax$Taxon <- apply(
  isa_tax[, taxonomy_ranks],
  1,
  function(x) {
    x <- x[!is.na(x)]
    if(length(x) == 0) NA else tail(x, 1)
  }
)

isa_tax$Taxon <- as.character(isa_tax$Taxon)

############################################################
# STEP 7 — Filter significance threshold
############################################################

isa_table <- isa_tax %>%
  filter(p.value <= 0.05) %>%
  select(Taxon, Group, stat, p.value) %>%
  arrange(Group, desc(stat))

############################################################
# STEP 8 — Export table
############################################################

write.csv(
  isa_table,
  "Indicator_species_results_unfiltered.csv",
  row.names = FALSE
)

############################################################
# STEP 9 — Filter table

library(dplyr)
library(stringr)

# Load data
df <- read.csv("Indicator_species_results_unfiltered.csv")

# Filter, deduplicate, round, extract rank, clean taxon names
df <- df %>%
  filter(stat >= 0.7, Group != "control") %>%
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

# Split into Exposed and Unexposed, then combine with a section header
exposed   <- df %>% filter(Group == "exposed")   %>% select(-Group)
unexposed <- df %>% filter(Group == "unexposed") %>% select(-Group)

blank <- data.frame(Taxon = NA, Rank = NA, stat = NA, p.value = NA)

exposed_labeled <- bind_rows(
  data.frame(Taxon = "=== Exposed ===", Rank = NA, stat = NA, p.value = NA),
  exposed
)
unexposed_labeled <- bind_rows(
  data.frame(Taxon = "=== Unexposed ===", Rank = NA, stat = NA, p.value = NA),
  unexposed
)

final_df <- bind_rows(exposed_labeled, blank, unexposed_labeled)

# Write to CSV
write.csv(final_df, "Indicator_species_filtered.csv", row.names = FALSE, na = "")