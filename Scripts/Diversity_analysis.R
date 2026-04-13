# --- 1. SETUP & LIBRARIES ---
# Required: install.packages(c("phyloseq", "ape", "tidyverse", "vegan", "picante", "FSA"))
library(phyloseq); library(ape); library(tidyverse); library(vegan); library(picante); library(FSA)

# --- 2. DATA LOADING & ENGINEERING ---
meta_raw <- read_delim("microgravity_metadata.tsv", delim="\t")
otu_raw  <- read_delim("feature-table.txt", delim="\t", skip=1)
tax_raw  <- read_delim("taxonomy.tsv", delim="\t")
tree     <- read.tree("tree.nwk")

# Format OTU Table
otu_mat <- as.matrix(otu_raw[,-1]); rownames(otu_mat) <- otu_raw$`#OTU ID`
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)

# Advanced Metadata Construction
meta <- meta_raw %>%
  mutate(
    Sample_State = ifelse(grepl("Necropsy", `sample-id`, ignore.case = TRUE), "Post", "Pre"),
    Group = str_remove_all(host_id, "\\d+"),
    Timeline = case_when(
      Time == "Launch plus 0" ~ "W0_Pre",
      Time == "Launch plus 4.5" ~ "W4.5_Pre",
      Time == "Launch plus 9" & Sample_State == "Pre" ~ "W9_Pre",
      Time == "Launch plus 9" & Sample_State == "Post" ~ "W9_Post",
      TRUE ~ "Baseline"
    )
  )
SAMP <- sample_data(meta); rownames(SAMP) <- meta$`sample-id`

# Format Taxonomy
tax_mat <- tax_raw %>% select(-Confidence) %>%
  separate(Taxon, sep="; ", into = c("D","P","C","O","F","G","S"), fill="right") %>%
  as.matrix(); tax_mat <- tax_mat[,-1]; rownames(tax_mat) <- tax_raw$`Feature ID`
TAX <- tax_table(tax_mat)

# Build Phyloseq & Pre-process
ps <- phyloseq(OTU, SAMP, TAX, tree)
ps_rare <- rarefy_even_depth(filter_taxa(ps, function(x) sum(x) > 5, prune=TRUE), rngseed=1, verbose=FALSE)

# Calculate Alpha Metrics
sample_data(ps_rare)$Shannon <- estimate_richness(ps_rare, measures="Shannon")$Shannon
sample_data(ps_rare)$FaithPD <- pd(t(as(otu_table(ps_rare), "matrix")), phy_tree(ps_rare), include.root=TRUE)$PD

# --- 3. MASTER PLOTTING & STATS FUNCTIONS ---

# Function for 2-Group Comparisons (Stats + Plots)
run_test_2group <- function(ps_obj, targets, group_col, title, use_wunifrac = TRUE) {
  cat("\nTEST:", title, "\n")
  ps_sub <- prune_samples(sample_data(ps_obj)[[group_col]] %in% targets, ps_obj)
  df <- data.frame(sample_data(ps_sub))
  df[[group_col]] <- factor(df[[group_col]], levels = targets)

  # Stats: Alpha
  p_shannon <- wilcox.test(Shannon ~ get(group_col), data=df)$p.value
  p_faith   <- wilcox.test(FaithPD ~ get(group_col), data=df)$p.value
  cat("alpha shannon: p =", round(p_shannon, 5), "\n")
  cat("alpha faith PD: p =", round(p_faith, 5), "\n")

  # Alpha Plot
  p1 <- ggplot(df, aes(x=.data[[group_col]], y=Shannon, fill=.data[[group_col]])) +
    geom_boxplot() + theme_bw() + labs(title=paste("Shannon:", title))
  print(p1)

  # Stats: Beta Bray-Curtis
  dist_bray <- phyloseq::distance(ps_sub, "bray")
  perm_bray <- adonis2(dist_bray ~ df[[group_col]])
  cat("beta bray curtis: R^2 =", round(perm_bray[1,4], 4), ", p =", perm_bray[1,5], "\n")
  
  # Stats: Beta Weighted UniFrac
  if(use_wunifrac) {
    dist_wuni <- phyloseq::distance(ps_sub, "wunifrac")
    perm_wuni <- adonis2(dist_wuni ~ df[[group_col]])
    disper_p  <- permutest(betadisper(dist_wuni, df[[group_col]]))$tab[1,6]
    cat("beta weighted unifrac: R^2 =", round(perm_wuni[1,4], 4), ", p =", perm_wuni[1,5], "\n")
    cat("betadisper (wunifrac): p =", round(disper_p, 4), "\n")
    
    # Beta Plot
    p2 <- plot_ordination(ps_sub, ordinate(ps_sub, "PCoA", "wunifrac"), color=group_col) +
      geom_point(size=3) + stat_ellipse() + theme_minimal() + labs(title=paste("W-UniFrac:", title))
    print(p2)
  }
}

# --- 4. EXECUTION ---
sink("Statistical_Analysis_Results.txt")
pdf("Biological_Visualizations.pdf", width=11, height=8.5)

# 1. Microgravity vs Control (F-ISST x GC-ISST)
run_test_2group(subset_samples(ps_rare, Sample_State == "Post"), c("F-ISST", "GC-ISST"), "Group", "F-ISST x GC-ISST")

# 2. Longitudinal within Cohorts (GC-LAR and F-LAR)
for(cohort in c("GC-LAR", "F-LAR")) {
  cat("\nTEST: Longitudinal", cohort, "\n")
  ps_sub <- prune_samples(sample_data(ps_rare)$Group == cohort, ps_rare)
  df <- data.frame(sample_data(ps_sub))
  
  # Alpha (Kruskal + Dunn)
  k_shan <- kruskal.test(Shannon ~ Timeline, data=df)
  cat("alpha shannon: p =", round(k_shan$p.value, 5), "\n")
  if(k_shan$p.value < 0.05) print(dunnTest(Shannon ~ Timeline, data=df, method="bh")$res)
  
  # Beta
  perm_bray <- adonis2(phyloseq::distance(ps_sub, "bray") ~ df$Timeline)
  perm_wuni <- adonis2(phyloseq::distance(ps_sub, "wunifrac") ~ df$Timeline)
  cat("beta bray curtis: R^2 =", round(perm_bray[1,4], 4), ", p =", perm_bray[1,5], "\n")
  cat("beta weighted unifrac: R^2 =", round(perm_wuni[1,4], 4), ", p =", perm_wuni[1,5], "\n")
  
  p_long <- ggplot(df, aes(x=Timeline, y=Shannon, fill=Timeline)) + geom_boxplot() + labs(title=paste("Timeline:", cohort))
  print(p_long)
}

# 3. Pairwise Time-points (GC-LAR x F-LAR)
for(tp in c("W0_Pre", "W4.5_Pre", "W9_Pre", "W9_Post")) {
  run_test_2group(subset_samples(ps_rare, Timeline == tp), c("GC-LAR", "F-LAR"), "Group", paste("GC-LAR vs F-LAR at", tp))
}

# 4. Recovery vs Flight (F-ISST x F-LAR)
run_test_2group(subset_samples(ps_rare, Timeline == "W9_Post"), c("F-ISST", "F-LAR"), "Group", "Recovery: F-ISST x F-LAR")

# 5. Control Comparison (GC-LAR x GC-ISST)
run_test_2group(subset_samples(ps_rare, Timeline == "W9_Post"), c("GC-LAR", "GC-ISST"), "Group", "Control: GC-LAR x GC-ISST")

# 6. Control (W9) vs Basal (Combined BL + W0)
sample_data(ps_rare)$Basal_Logic <- ifelse(sample_data(ps_rare)$Group == "BL" | sample_data(ps_rare)$Timeline == "W0_Pre", "Basal", "Other")
run_test_2group(subset_samples(ps_rare, (Group == "GC-LAR" & Timeline == "W9_Post") | Basal_Logic == "Basal"), 
                c("Basal", "Other"), "Basal_Logic", "GC-LAR (W9) x Basal", use_wunifrac = FALSE)

dev.off(); sink()
