library(tidyverse)
library(phyloseq)
library(DESeq2)
library(dplyr)
library(purrr)
library(tibble)
library(ggplot2)
library(readr)

# Load data and setup
load("micro_final.RData")

# Dividing by treatment group
sample_data(micro_final)$treatment_group <- factor(
  sample_data(micro_final)$treatment_group,
  levels = c("BL", "F-LAR", "GC-LAR", "F-ISST", "GC-ISST")
)

# Dividing by sample time
sample_data(micro_final)$Time <- factor(
  sample_data(micro_final)$Time,
  levels = c("Launch plus 0", "Launch plus 4.5", "Launch plus 9")
)

# Dividing by sample type
sample_data(micro_final)$sample_type_group <- sample_data(micro_final)$treatment_group %>%
  as.character() %>%  
  ifelse(
    sample_data(micro_final)$Time == "Launch plus 9" & 
      sample_data(micro_final)$treatment_group %in% c("F-LAR", "GC-LAR"),
    paste0(
      sample_data(micro_final)$treatment_group,
      case_when(
        sample_data(micro_final)$sample_type == "Fresh" ~ "-Fresh",
        sample_data(micro_final)$sample_type == "Necropsy" ~ "-Necropsy",
        TRUE ~ "-Other"
      )
    ),
    .
  ) %>%
  factor(levels = c("BL", "F-LAR", "F-LAR-Fresh", "F-LAR-Necropsy", "GC-LAR", 
                    "GC-LAR-Fresh", "GC-LAR-Necropsy", "F-ISST", "GC-ISST"))

# Add 1 to 0s in the data to allow analysis
micro_plus1 <- transform_sample_counts(micro_final, function(x) x + 1)

#### DESeq and Volcano plots ####

# Helper function to use sample_type_group
run_deseq_treatment <- function(micro_sub, group1, group2) {
  groups <- c(group1, group2)
  keep <- sample_data(micro_sub)$sample_type_group %in% groups  
  micro_sub <- prune_samples(keep, micro_sub)
  dds <- phyloseq_to_deseq2(micro_sub, ~ sample_type_group) 
  dds <- DESeq(dds)
  results(dds, tidy = TRUE, contrast = c("sample_type_group", group1, group2))
}

# Helper function to make volcano plots
volcano_plot <- function(res, filename) {
  p <- res %>%
    mutate(significant = padj < 0.05 & abs(log2FoldChange) > 2) %>%
    ggplot() + 
    geom_point(aes(x = log2FoldChange, y = -log10(padj), col = significant)) +
    theme_minimal()
  print(p)
  ggsave(filename, p)
}


# Group comparisons:
# 1. Microgravity vs control: F-ISST x GC-ISST (Launch plus 9) 
micro_iss_9 <- subset_samples(micro_plus1, treatment_group %in% c("F-ISST", "GC-ISST") & Time == "Launch plus 9")
res_microgravity <- run_deseq_treatment(micro_iss_9, "F-ISST", "GC-ISST")

# Save volcano plot
volcano_plot(res_microgravity, "vol_plot_F-ISST_vs_GC-ISST_9.png")

# 2. Timepoints within cohorts: GC-LAR and F-LAR
# F-LAR
micro_FLAR <- subset_samples(micro_plus1, treatment_group == "F-LAR")
sample_data(micro_FLAR)$host_id <- factor(sample_data(micro_FLAR)$host_id)

sample_data(micro_FLAR)$time_group <- factor(
  case_when(
    sample_data(micro_FLAR)$Time == "Launch plus 0" ~ "Week0",
    sample_data(micro_FLAR)$Time == "Launch plus 4.5" ~ "Week4.5", 
    sample_data(micro_FLAR)$Time == "Launch plus 9" & 
      sample_data(micro_FLAR)$sample_type == "Fresh" ~ "Week9_Fresh",
    sample_data(micro_FLAR)$Time == "Launch plus 9" & 
      sample_data(micro_FLAR)$sample_type == "Necropsy" ~ "Week9_Necropsy",
    TRUE ~ "Other"
  ),
  levels = c("Week0", "Week4.5", "Week9_Fresh", "Week9_Necropsy")
)

dds_FLAR <- phyloseq_to_deseq2(micro_FLAR, ~ host_id + time_group)
dds_FLAR <- DESeq(dds_FLAR)

# Save volcano plots
res_FLAR_45_vs_0 <- results(dds_FLAR, contrast = c("time_group", "Week4.5", "Week0"), tidy = TRUE)
volcano_plot(res_FLAR_45_vs_0, "vol_plot_F-LAR_4.5_vs_0.png")
res_FLAR_9fresh_vs_0 <- results(dds_FLAR, contrast = c("time_group", "Week9_Fresh", "Week0"), tidy = TRUE)
volcano_plot(res_FLAR_9fresh_vs_0, "vol_plot_F-LAR_9fresh_vs_0.png")
res_FLAR_9necro_vs_0 <- results(dds_FLAR, contrast = c("time_group", "Week9_Necropsy", "Week0"), tidy = TRUE)
volcano_plot(res_FLAR_9necro_vs_0, "vol_plot_F-LAR_9necro_vs_0.png")
res_FLAR_9fresh_vs_45 <- results(dds_FLAR, contrast = c("time_group", "Week9_Fresh", "Week4.5"), tidy = TRUE)
volcano_plot(res_FLAR_9fresh_vs_45, "vol_plot_F-LAR_9fresh_vs_4.5.png")
res_FLAR_9necro_vs_45 <- results(dds_FLAR, contrast = c("time_group", "Week9_Necropsy", "Week4.5"), tidy = TRUE)
volcano_plot(res_FLAR_9necro_vs_45, "vol_plot_F-LAR_9necro_vs_4.5.png")
res_FLAR_9fresh_vs_9necro <- results(dds_FLAR, contrast = c("time_group", "Week9_Fresh", "Week9_Necropsy"), tidy = TRUE)
volcano_plot(res_FLAR_9fresh_vs_9necro, "vol_plot_F-LAR_9fresh_vs_9necro.png")

# GC-LAR
micro_GCLAR <- subset_samples(micro_plus1, treatment_group == "GC-LAR")
sample_data(micro_GCLAR)$host_id <- factor(sample_data(micro_GCLAR)$host_id)

sample_data(micro_GCLAR)$time_group <- factor(
  case_when(
    sample_data(micro_GCLAR)$Time == "Launch plus 0" ~ "Week0",
    sample_data(micro_GCLAR)$Time == "Launch plus 4.5" ~ "Week4.5", 
    sample_data(micro_GCLAR)$Time == "Launch plus 9" & 
      sample_data(micro_GCLAR)$sample_type == "Fresh" ~ "Week9_Fresh",
    sample_data(micro_GCLAR)$Time == "Launch plus 9" & 
      sample_data(micro_GCLAR)$sample_type == "Necropsy" ~ "Week9_Necropsy",
    TRUE ~ "Other"
  ),
  levels = c("Week0", "Week4.5", "Week9_Fresh", "Week9_Necropsy")
)

dds_GCLAR <- phyloseq_to_deseq2(micro_GCLAR, ~ host_id + time_group)
dds_GCLAR <- DESeq(dds_GCLAR)

# Save volcano plots
res_GCLAR_45_vs_0 <- results(dds_GCLAR, contrast = c("time_group", "Week4.5", "Week0"), tidy = TRUE)
volcano_plot(res_GCLAR_45_vs_0, "vol_plot_GC-LAR_4.5_vs_0.png")
res_GCLAR_9fresh_vs_0 <- results(dds_GCLAR, contrast = c("time_group", "Week9_Fresh", "Week0"), tidy = TRUE)
volcano_plot(res_GCLAR_9fresh_vs_0, "vol_plot_GC-LAR_9fresh_vs_0.png")
res_GCLAR_9necro_vs_0 <- results(dds_GCLAR, contrast = c("time_group", "Week9_Necropsy", "Week0"), tidy = TRUE)
volcano_plot(res_GCLAR_9necro_vs_0, "vol_plot_GC-LAR_9necro_vs_0.png")
res_GCLAR_9fresh_vs_45 <- results(dds_GCLAR, contrast = c("time_group", "Week9_Fresh", "Week4.5"), tidy = TRUE)
volcano_plot(res_GCLAR_9fresh_vs_45, "vol_plot_GC-LAR_9fresh_vs_4.5.png")
res_GCLAR_9necro_vs_45 <- results(dds_GCLAR, contrast = c("time_group", "Week9_Necropsy", "Week4.5"), tidy = TRUE)
volcano_plot(res_GCLAR_9necro_vs_45, "vol_plot_GC-LAR_9necro_vs_4.5.png")
res_GCLAR_9fresh_vs_9necro <- results(dds_GCLAR, contrast = c("time_group", "Week9_Fresh", "Week9_Necropsy"), tidy = TRUE)
volcano_plot(res_GCLAR_9fresh_vs_9necro, "vol_plot_GC-LAR_9fresh_vs_9necro.png")

# 3. Same time-points between cohorts: GC-LAR x F-LAR
time_comparisons <- list(
  list(time = "Launch plus 0",   groups = c("F-LAR", "GC-LAR"), use_time = TRUE),
  list(time = "Launch plus 4.5", groups = c("F-LAR", "GC-LAR"), use_time = TRUE),
  list(groups = c("F-LAR-Fresh", "GC-LAR-Fresh"),         use_time = FALSE),
  list(groups = c("F-LAR-Necropsy", "GC-LAR-Necropsy"),   use_time = FALSE)
)

res_LAR_pairs <- purrr::map_dfr(seq_along(time_comparisons), function(i) {
  cmp <- time_comparisons[[i]]
  
  if (cmp$use_time) {
    keep <- sample_data(micro_plus1)$Time == cmp$time &
      sample_data(micro_plus1)$treatment_group %in% cmp$groups
    comp_label <- paste(cmp$time, paste(cmp$groups, collapse = " vs "), sep = " | ")
  } else {
    keep <- sample_data(micro_plus1)$sample_type_group %in% cmp$groups
    comp_label <- paste(cmp$groups, collapse = " vs ")
  }
  
  if (sum(keep, na.rm = TRUE) < 4) {
    warning(paste("Skipping", comp_label, "- only", sum(keep, na.rm = TRUE), "samples"))
    return(tibble())
  }
  
  micro_sub <- prune_samples(keep, micro_plus1)
  if (nsamples(micro_sub) == 0) return(tibble())
  
  res_tmp <- run_deseq_treatment(micro_sub, cmp$groups[1], cmp$groups[2])
  res_tmp$comparison <- comp_label
  res_tmp
})

# Save volcano plots
res_LAR_pairs %>%
  filter(comparison == "Launch plus 0 | F-LAR vs GC-LAR") %>%
  volcano_plot("vol_plot_F-LAR_vs_GC-LAR_0.png")
res_LAR_pairs %>%
  filter(comparison == "Launch plus 4.5 | F-LAR vs GC-LAR") %>%
  volcano_plot("vol_plot_F-LAR_vs_GC-LAR_4.5.png")
res_LAR_pairs %>%
  filter(comparison == "F-LAR-Fresh vs GC-LAR-Fresh") %>%
  volcano_plot("vol_plot_F-LAR-Fresh_vs_GC-LAR-Fresh.png")
res_LAR_pairs %>%
  filter(comparison == "F-LAR-Necropsy vs GC-LAR-Necropsy") %>%
  volcano_plot("vol_plot_F-LAR-Necropsy_vs_GC-LAR-Necropsy.png")

# 4. Return vs flight at 9 weeks after euthanasia: F-ISST x F-LAR
micro_return_9 <- subset_samples(micro_plus1, 
                                 sample_type_group %in% c("F-ISST", "F-LAR-Necropsy") & 
                                   Time == "Launch plus 9")
res_return_flight <- run_deseq_treatment(micro_return_9, "F-ISST", "F-LAR-Necropsy")
volcano_plot(res_return_flight, "vol_plot_F-ISST_vs_F-LAR-Necropsy_9.png")

# 5. Controls week 9 after euthanasia: GC-LAR x GC-ISST
micro_controls_9 <- subset_samples(micro_plus1, 
                                   sample_type_group %in% c("GC-LAR-Necropsy", "GC-ISST") & 
                                     Time == "Launch plus 9")
res_controls_9 <- run_deseq_treatment(micro_controls_9, "GC-LAR-Necropsy", "GC-ISST")
volcano_plot(res_controls_9, "vol_plot_GC-LAR-Necropsy_vs_GC-ISST_9.png")

# 6. Control vs basal week 0: GC-LAR x BL
micro_ctrl_basal_0 <- subset_samples(micro_plus1, 
                                     treatment_group %in% c("GC-LAR", "BL") & 
                                       Time == "Launch plus 0")
res_ctrl_basal <- run_deseq_treatment(micro_ctrl_basal_0, "GC-LAR", "BL")  # Week 0 = "GC-LAR" (no split)
volcano_plot(res_ctrl_basal, "vol_plot_GC-LAR_vs_BL_0_all_samples.png")


#### Count up/down ASVs for each comparison #### 
lar_0w  <- res_LAR_pairs %>% filter(comparison == "Launch plus 0 | F-LAR vs GC-LAR")
lar_45w  <- res_LAR_pairs %>% filter(comparison == "Launch plus 4.5 | F-LAR vs GC-LAR")

sig_summary <- tibble(
  comparison = c(
    "F-ISST vs GC-ISST (9w)",
    "F-LAR 4.5vs0", "F-LAR 9fresh vs0", "F-LAR 9necro vs0", "F-LAR 9fresh vs4.5",
    "F-LAR 9necro vs4.5", "F-LAR 9fresh vs9necro",
    "GC-LAR 4.5vs0", "GC-LAR 9fresh vs0", "GC-LAR 9necro vs0", "GC-LAR 9fresh vs4.5",
    "GC-LAR 9necro vs4.5", "GC-LAR 9fresh vs9necro",
    "F-LAR vs GC-LAR (0w)", "F-LAR vs GC-LAR (4.5w)",
    "F-LAR-Fresh vs GC-LAR-Fresh (9w)", "F-LAR-Necropsy vs GC-LAR-Necropsy (9w)",
    "F-ISST vs F-LAR-Necropsy (9w)", "GC-LAR-Necropsy vs GC-ISST (9w)", "GC-LAR vs BL (0w)"
  ),
  res_obj = list(
    res_microgravity,
    res_FLAR_45_vs_0, res_FLAR_9fresh_vs_0, res_FLAR_9necro_vs_0,
    res_FLAR_9fresh_vs_45, res_FLAR_9necro_vs_45, res_FLAR_9fresh_vs_9necro,
    res_GCLAR_45_vs_0, res_GCLAR_9fresh_vs_0, res_GCLAR_9necro_vs_0,
    res_GCLAR_9fresh_vs_45, res_GCLAR_9necro_vs_45, res_GCLAR_9fresh_vs_9necro,
    lar_0w,
    lar_45w,
    res_LAR_pairs %>% filter(comparison == "F-LAR-Fresh vs GC-LAR-Fresh"),
    res_LAR_pairs %>% filter(comparison == "F-LAR-Necropsy vs GC-LAR-Necropsy"),
    res_return_flight, res_controls_9, res_ctrl_basal
  )
) %>%
  rowwise() %>%
  mutate(
    tested = sum(!is.na(res_obj$padj)),
    up     = sum(res_obj$padj < 0.05 & res_obj$log2FoldChange > 2, na.rm = TRUE),
    down   = sum(res_obj$padj < 0.05 & res_obj$log2FoldChange < -2, na.rm = TRUE),
    ns     = tested - up - down
  ) %>%
  ungroup()

write_csv(sig_summary %>% select(comparison, tested, up, down, ns), "ASVs_count_summary.csv")
print(sig_summary %>% select(comparison, tested, up, down, ns))

#### Generate sig ASV barplots + save tables for each comparison ####
# Explicit LAR subsets
lar_0w  <- res_LAR_pairs %>% filter(comparison == "Launch plus 0 | F-LAR vs GC-LAR")
lar_45w <- res_LAR_pairs %>% filter(comparison == "Launch plus 4.5 | F-LAR vs GC-LAR")
# Make sure the two subsets exist
stopifnot(nrow(lar_0w) > 0, nrow(lar_45w) > 0)

sig_plots <- tibble(
  comparison = c(
    "F-ISST vs GC-ISST (9w)",
    "F-LAR 4.5vs0", "F-LAR 9fresh vs0", "F-LAR 9necro vs0", "F-LAR 9fresh vs4.5",
    "F-LAR 9necro vs4.5", "F-LAR 9fresh vs9necro",
    "GC-LAR 4.5vs0", "GC-LAR 9fresh vs0", "GC-LAR 9necro vs0", "GC-LAR 9fresh vs4.5",
    "GC-LAR 9necro vs4.5", "GC-LAR 9fresh vs9necro",
    "F-LAR vs GC-LAR (0w)", "F-LAR vs GC-LAR (4.5w)",
    "F-LAR-Fresh vs GC-LAR-Fresh (9w)", "F-LAR-Necropsy vs GC-LAR-Necropsy (9w)",
    "F-ISST vs F-LAR-Necropsy (9w)", "GC-LAR-Necropsy vs GC-ISST (9w)", "GC-LAR vs BL (0w)"
  ),
  res_obj = list(
    res_microgravity,
    res_FLAR_45_vs_0, res_FLAR_9fresh_vs_0, res_FLAR_9necro_vs_0,
    res_FLAR_9fresh_vs_45, res_FLAR_9necro_vs_45, res_FLAR_9fresh_vs_9necro,
    res_GCLAR_45_vs_0, res_GCLAR_9fresh_vs_0, res_GCLAR_9necro_vs_0,
    res_GCLAR_9fresh_vs_45, res_GCLAR_9necro_vs_45, res_GCLAR_9fresh_vs_9necro,
    lar_0w,
    lar_45w,
    res_LAR_pairs %>% filter(comparison == "F-LAR-Fresh vs GC-LAR-Fresh"),
    res_LAR_pairs %>% filter(comparison == "F-LAR-Necropsy vs GC-LAR-Necropsy"),
    res_return_flight, res_controls_9, res_ctrl_basal
  )
) %>%
  rowwise() %>%
  mutate(
    tested = sum(!is.na(res_obj$padj)),
    up     = sum(res_obj$padj < 0.05 & res_obj$log2FoldChange > 2, na.rm = TRUE),
    down   = sum(res_obj$padj < 0.05 & res_obj$log2FoldChange < -2, na.rm = TRUE),
    ns     = tested - up - down,
    sig_table = list({
      strict_sig <- res_obj %>%
        filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 2)
      
      if (nrow(strict_sig) == 0) {
        tibble()
      } else {
        sigASVs_vec <- strict_sig %>% pull(row)
        tax_df <- tax_table(prune_taxa(sigASVs_vec, micro_final)) %>%
          as.data.frame() %>%
          rownames_to_column("ASV")
        
        strict_sig_clean <- strict_sig %>%
          mutate(ASV = row) %>%
          filter(!is.na(ASV))
        
        tax_df %>%
          right_join(strict_sig_clean, by = "ASV") %>%
          arrange(log2FoldChange) %>%
          mutate(Genus = make.unique(Genus))
      }
    }),
    plot = list({
      if (nrow(sig_table) == 0) {
        ggplot() +
          theme_void() +
          labs(title = comparison, subtitle = "No strict-significant ASVs")
      } else {
        ggplot(sig_table) +
          geom_bar(aes(x = Genus, y = log2FoldChange), stat = "identity") +
          geom_errorbar(aes(
            x = Genus,
            ymin = log2FoldChange - lfcSE,
            ymax = log2FoldChange + lfcSE
          )) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
          labs(title = comparison)
      }
    }),
    filename = paste0("sig_ASVs_", gsub("[^A-Za-z0-9]", "_", comparison), ".csv")
  ) %>%
  ungroup()

# Save tables
purrr::walk2(sig_plots$sig_table, sig_plots$filename, readr::write_csv)

# Save plots for all 20 comparisons
purrr::walk2(
  sig_plots$plot,
  paste0("sig_ASVs_", gsub("[^A-Za-z0-9]", "_", sig_plots$comparison), ".png"),
  ~ggsave(filename = .y, plot = .x, width = 12, height = 6)
)
