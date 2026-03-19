library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)
library(indicspecies)
library(stringr)
setwd("~/Indicator Taxa")
#### Load data ####
metafp <- "microgravity_export/microgravity_metadata.tsv"
meta <- read_delim(metafp, delim = "\t")

otufp <- "microgravity_export/feature-table.txt"
otu <- read_delim(file = otufp, delim = "\t", skip = 1)

taxfp <- "microgravity_export/taxonomy.tsv"
tax <- read_delim(taxfp, delim = "\t")

phylotreefp <- "microgravity_export/tree.nwk"
phylotree <- read.tree(phylotreefp)

#### Format OTU table ####
otu_mat <- as.matrix(otu[, -1])
rownames(otu_mat) <- otu$`#OTU ID`
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)

#### Format sample metadata ####
meta <- meta %>%
  mutate(
    # sample_type from sample-id
    sample_type = str_extract(`sample-id`, "(?i)(Necropsy|Fresh)"),
    # treatment_group from host_id: remove any digits
    treatment_group = str_remove_all(host_id, "\\d+")
  )

samp_df <- as.data.frame(meta[, -1])
rownames(samp_df) <- meta$`sample-id`
SAMP <- sample_data(samp_df)

#### Formatting taxonomy ####
tax_mat <- tax %>%
  select(-Confidence) %>%
  separate(
    col  = Taxon,
    sep  = "; ",
    into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
    fill = "right"
  ) %>%
  as.matrix()

tax_mat <- tax_mat[, -1]
rownames(tax_mat) <- tax$`Feature ID`

common_taxa <- intersect(rownames(otu_mat), tax$`Feature ID`)
otu_mat <- otu_mat[common_taxa, ]
tax_mat <- tax_mat[common_taxa, ]
TAX <- tax_table(tax_mat)

#### Create phyloseq object ####
micro <- phyloseq(OTU, SAMP, TAX, phylotree)

######### FILTER / RAREFY ##########
micro_filt_nolow <- filter_taxa(micro, function(x) sum(x) > 5, prune = TRUE)
micro_final <- prune_samples(sample_sums(micro_filt_nolow) > 100, micro_filt_nolow)

# rarecurve(t(as.data.frame(otu_table(micro_final))), cex = 0.1)
micro_rare <- rarefy_even_depth(micro_final, rngseed = 1, sample.size = 64411)

save(micro_final, file = "micro_final.RData")
save(micro_rare,  file = "micro_rare.RData")

#### INDICATOR SPECIES ANALYSIS ####

# Work on the rarefied object
ps <- micro_rare

# Extract sample data as data frame for easier manipulation
sdat <- as(sample_data(ps), "data.frame")

## 1) Parse condition from treatment_group (already created above as letters-only host_id)
##    Example: "F9-LAR" -> treatment_group "F-LAR", "GC3-ISST" -> "GC-ISST", "BL1" -> "BL"
##    If you've already done this when building SAMP, treatment_group is available.

# (If you had not already made treatment_group, it would be something like:)
# sdat$treatment_group <- str_remove_all(sdat$host_id, "\\d+")

## 2) Parse week from sample-id: "…_Week 9_…" -> "9", "…_Week 4.5_…" -> "4.5"

library(stringr)
sdat <- sdat %>%
  mutate(
    Week = case_when(
      Time == "Launch plus 0"     ~ "0",
      Time == "Launch plus 4.5"   ~ "4.5", 
      Time == "Launch plus 9"     ~ "9",
      TRUE ~ NA_character_
    )
  )


# Optional checks:
#table(sdat$treatment_group)
#table(sdat$Week)

## 3) ISA‑1: exposed vs unexposed vs control

isa1_group <- with(sdat, case_when(
  # exposed to microgravity:
  (treatment_group == "F-ISST" & Week == "9")              ~ "exposed",
  (treatment_group == "F-LAR"  & Week %in% c("4.5", "9"))  ~ "exposed",
  
  # no exposure to microgravity:
  (treatment_group == "GC-ISST" & Week == "9")             ~ "unexposed",
  (treatment_group == "GC-LAR"  & Week %in% c("0", "4.5", "9")) ~ "unexposed",
  (treatment_group == "F-LAR"   & Week == "0")             ~ "unexposed",
  
  # control:
  (treatment_group == "BL"      & Week == "0")             ~ "control",
  
  TRUE ~ NA_character_
))

isa1_group <- factor(isa1_group,
                     levels = c("control", "unexposed", "exposed"))

keep1 <- !is.na(isa1_group)
ps_isa1 <- prune_samples(keep1, ps)
isa1_group <- droplevels(isa1_group[keep1])

# OTU table for ISA‑1 (taxa in columns)
otu_isa1 <- as(otu_table(ps_isa1), "matrix")
if (taxa_are_rows(ps_isa1)) {
  otu_isa1 <- t(otu_isa1)
}

indicator_multipatt_1 <- multipatt(otu_isa1, isa1_group, duleg = TRUE)

summary(indicator_multipatt_1)

indicator_output_1 <- capture.output(summary(indicator_multipatt_1, indvalcomp = TRUE))
writeLines(indicator_output_1, "indicator_values_ISA1.txt")

## 4) ISA‑2: nine specific condition × time categories

sdat2 <- sdat[rownames(sdat) %in% sample_names(ps), ]

sdat2 <- sdat2 %>%
  mutate(CondWeek = paste(treatment_group, Week, sep = "_"))

isa2_levels <- c(
  "F-ISST_9",
  "GC-ISST_9",
  "F-LAR_0",
  "F-LAR_4.5",
  "F-LAR_9",
  "GC-LAR_0",
  "GC-LAR_4.5",
  "GC-LAR_9",
  "BL_0"
)

isa2_group <- ifelse(sdat2$CondWeek %in% isa2_levels, sdat2$CondWeek, NA)
isa2_group <- factor(isa2_group, levels = isa2_levels)

keep2 <- !is.na(isa2_group)
ps_isa2 <- prune_samples(keep2, ps)
isa2_group <- droplevels(isa2_group[keep2])

otu_isa2 <- as(otu_table(ps_isa2), "matrix")
if (taxa_are_rows(ps_isa2)) {
  otu_isa2 <- t(otu_isa2)
}

indicator_multipatt_2 <- multipatt(otu_isa2, isa2_group, duleg = TRUE)

summary(indicator_multipatt_2)

indicator_output_2 <- capture.output(summary(indicator_multipatt_2, indvalcomp = TRUE))
writeLines(indicator_output_2, "indicator_values_ISA2.txt")

# =============================================================================
# FIXED ISA VISUALIZATION SCRIPT
# =============================================================================

library(ggplot2)
library(dplyr)
library(viridis)
library(stringr)
library(phyloseq)

# ===== 1. CORRECTED EXTRACTION (access $sign directly on the multipatt object) =====
extract_isa_results <- function(mp_obj, ps_obj, alpha = 0.05) {
  
  # $sign lives on the multipatt object itself, not on summary()
  sign_df <- as.data.frame(mp_obj$sign)
  sign_df$Taxon <- rownames(sign_df)
  
  # Filter by p-value
  sign_df <- sign_df[sign_df$p.value <= alpha, ]
  
  if (nrow(sign_df) == 0) return(NULL)
  
  # Identify which group columns are present (everything except stat, p.value, Taxon)
  group_cols <- setdiff(colnames(sign_df), c("stat", "p.value", "Taxon", "index"))
  
  # Decode the numeric group index to a label
  # mp_obj$sign$index is an integer index into the combination list
  cluster_names <- colnames(mp_obj$A)  # group labels in order
  
  sign_df$Group <- cluster_names[sign_df$index]
  sign_df$IndVal <- sign_df$stat
  
  # Attach taxonomy if available
  if (!is.null(ps_obj)) {
    tax_df <- as.data.frame(tax_table(ps_obj))
    tax_df$Taxon <- rownames(tax_df)
    sign_df <- left_join(sign_df, tax_df[, c("Taxon", "Genus", "Family", "Phylum")],
                         by = "Taxon")
    # Build a readable label: prefer Genus, fall back to Family, then Phylum, then raw ID
    sign_df$Label <- coalesce(
      na_if(sign_df$Genus,  "NA"), na_if(sign_df$Genus,  ""),
      na_if(sign_df$Family, "NA"), na_if(sign_df$Family, ""),
      na_if(sign_df$Phylum, "NA"), na_if(sign_df$Phylum, ""),
      sign_df$Taxon
    )
  } else {
    sign_df$Label <- sign_df$Taxon
  }
  
  sign_df <- sign_df[, c("Taxon", "Label", "Group", "IndVal", "p.value")]
  sign_df[order(sign_df$IndVal, decreasing = TRUE), ]
}

# ===== 2. HEATMAP — capped at top 50 taxa so the file stays renderable =====
plot_isa_heatmap <- function(df, title, filename, top_n = 50) {
  if (is.null(df) || nrow(df) == 0) { message("No data for: ", filename); return() }
  
  plot_df <- df %>% slice_max(IndVal, n = top_n)
  
  p <- ggplot(plot_df, aes(x = Group, y = reorder(Label, IndVal), fill = IndVal)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_viridis_c(name = "IndVal", option = "plasma", limits = c(0, 1)) +
    geom_text(aes(label = sprintf("%.2f", IndVal)), size = 2.5, color = "white") +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.y  = element_text(size = 8),
      axis.text.x  = element_text(angle = 45, hjust = 1),
      plot.title   = element_text(hjust = 0.5, face = "bold"),
      panel.grid   = element_blank()
    ) +
    labs(title = paste0(title, "\n(top ", top_n, " by IndVal)"),
         x = "Group", y = "Taxon")
  
  # Fixed, safe dimensions
  ggsave(filename, p, width = 10, height = 14, dpi = 300, limitsize = FALSE)
  message("Saved: ", filename)
}

# ===== 3. TOP 15 BARPLOT =====
plot_top_indicators <- function(df, filename) {
  if (is.null(df) || nrow(df) == 0) { message("No data for: ", filename); return() }
  
  top_df <- df %>%
    slice_max(IndVal, n = 15) %>%
    mutate(Label_short = str_trunc(Label, 35))
  
  p <- ggplot(top_df, aes(x = reorder(Label_short, IndVal), y = IndVal, fill = Group)) +
    geom_col(alpha = 0.85) +
    coord_flip() +
    scale_fill_brewer(type = "qual", palette = "Set2") +
    theme_minimal(base_size = 12) +
    labs(title = "Top 15 Indicator Taxa (Highest IndVal)",
         x = "Taxon", y = "Indicator Value") +
    theme(legend.position = "bottom")
  
  ggsave(filename, p, width = 12, height = 8, dpi = 300)
  message("Saved: ", filename)
}

# ===== 4. ALPHA DIVERSITY =====
plot_diversity <- function() {
  ps_plot <- ps
  sample_data(ps_plot)$isa1_group <- NA
  sample_data(ps_plot)$isa1_group[keep1] <- as.character(isa1_group)
  
  p <- plot_richness(ps_plot, x = "isa1_group",
                     measures = c("Shannon", "Simpson")) +
    theme_bw() +
    labs(title = "Alpha Diversity by Microgravity Exposure",
         x = "Exposure Group", y = "Diversity Index") +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 30, hjust = 1))
  
  ggsave("diversity_exposure.png", p, width = 8, height = 6, dpi = 300)
  message("Saved: diversity_exposure.png")
}

# ===== 5. PCoA =====
plot_pcoa <- function() {
  ps_plot <- ps
  sample_data(ps_plot)$isa1_group <- NA
  sample_data(ps_plot)$isa1_group[keep1] <- as.character(isa1_group)
  
  ord <- ordinate(ps_plot, "PCoA", "bray")
  p <- plot_ordination(ps_plot, ord, color = "isa1_group") +
    theme_bw() +
    labs(title = "Beta Diversity (PCoA, Bray-Curtis) by Exposure Group") +
    scale_color_brewer(type = "qual", palette = "Set1", na.value = "grey70")
  
  ggsave("pcoa_exposure.png", p, width = 8, height = 6, dpi = 300)
  message("Saved: pcoa_exposure.png")
}

# ===== 6. RUN =====

isa1_df <- extract_isa_results(indicator_multipatt_1, ps_isa1)
isa2_df <- extract_isa_results(indicator_multipatt_2, ps_isa2)

cat("ISA-1 significant taxa:", if (!is.null(isa1_df)) nrow(isa1_df) else 0, "\n")
cat("ISA-2 significant taxa:", if (!is.null(isa2_df)) nrow(isa2_df) else 0, "\n")

plot_isa_heatmap(isa1_df, "ISA-1: Microgravity Exposure Indicators", "ISA1_heatmap.png")
plot_top_indicators(isa1_df, "ISA1_top15.png")
plot_isa_heatmap(isa2_df, "ISA-2: Condition × Time Indicators",     "ISA2_heatmap.png")

plot_diversity()
plot_pcoa()

if (!is.null(isa1_df)) write.csv(isa1_df, "ISA1_significant_157taxa.csv", row.names = FALSE)
if (!is.null(isa2_df)) write.csv(isa2_df, "ISA2_significant_177taxa.csv", row.names = FALSE)
