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
  filter(p.value <= 0.07) %>%
  select(Taxon, Group, stat, p.value) %>%
  arrange(Group, desc(stat))

############################################################
# STEP 8 — Export table
############################################################

write.csv(
  isa_table,
  "Indicator_species_results_p0.07.csv",
  row.names = FALSE
)

############################################################
# STEP 9 — Confirm success
############################################################

str(isa_table)