# =========================
# 1. Load libraries
# =========================
library(tidyverse)
library(caret)
library(randomForest)
library(ranger)
library(pROC)
library(boot)
library(pheatmap)

# =========================
# 2. Read and parse ISA table
# =========================
isa_lines <- readLines("indicator_values_ISA1.txt")

current_group <- NA_character_
isa_list <- list()

for (ln in isa_lines) {
  ln_trim <- trimws(ln)
  
  if (grepl("^Group\\s+", ln_trim)) {
    current_group <- sub("^Group\\s+([^ ]+).*$", "\\1", ln_trim)
  } else if (grepl("^[0-9a-f]{32}", ln_trim)) {
    parts <- strsplit(ln_trim, "\\s+")[[1]]
    
    isa_list[[length(isa_list) + 1]] <- data.frame(
      FeatureID = parts[1],
      isa_group = current_group,
      A = as.numeric(parts[2]),
      B = as.numeric(parts[3]),
      stat = as.numeric(parts[4]),
      p.value = as.numeric(parts[5]),
      stringsAsFactors = FALSE
    )
  }
}

isa_df <- bind_rows(isa_list)

# Keep only exposed/unexposed ISA features with stat > 0.7
selected_ids <- isa_df %>%
  filter(isa_group %in% c("exposed", "unexposed"),
         stat > 0.7) %>%
  pull(FeatureID) %>%
  unique()

if (length(selected_ids) == 0) {
  stop("No ISA features passed the filter (exposed/unexposed with stat > 0.7).")
}

# =========================
# 3. Load files
# =========================
otu <- read.csv("micro_final_otu.csv", row.names = 1, check.names = FALSE)
tax <- read.delim("taxonomy.tsv", sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
meta <- read.delim("microgravity_metadata.tsv", sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

# =========================
# 4. Keep only exposed and unexposed samples
# =========================
meta <- meta %>%
  mutate(treatment_group = case_when(
    Condition == "Space Flight" ~ "exposed",
    Condition == "Ground Control" ~ "unexposed",
    Condition == "Basal Control" ~ "control"
  ))

meta_sub <- meta %>%
  filter(treatment_group %in% c("exposed", "unexposed"))

if (nrow(meta_sub) == 0) {
  stop("No metadata rows found for treatment_group = exposed/unexposed.")
}

# Keep only ISA-selected OTU rows and exposed/unexposed sample columns
otu_sub <- otu[
  rownames(otu) %in% selected_ids,
  colnames(otu) %in% meta_sub$`sample-id`,
  drop = FALSE
]

if (nrow(otu_sub) == 0) {
  stop("No OTU features left after applying ISA feature filter.")
}

# Make metadata order match OTU column order
meta_sub <- meta_sub %>%
  filter(`sample-id` %in% colnames(otu_sub)) %>%
  arrange(match(`sample-id`, colnames(otu_sub)))

# Reorder OTU columns to match metadata
otu_sub <- otu_sub[, meta_sub$`sample-id`, drop = FALSE]

# =========================
# 5. Make genus-level count table
# =========================
tax2 <- tax %>%
  rename(FeatureID = `Feature ID`) %>%
  select(FeatureID, Taxon)

# Extract genus from taxonomy string
tax2$Genus <- stringr::str_extract(tax2$Taxon, "g__[^;]+")
tax2$Genus <- gsub("g__", "", tax2$Genus)
tax2$Genus[is.na(tax2$Genus) | tax2$Genus == ""] <- "Unassigned"

# Keep only selected ISA features
tax2 <- tax2 %>% filter(FeatureID %in% rownames(otu_sub))

# Add FeatureID as a column so we can join
otu_sub$FeatureID <- rownames(otu_sub)

# Join counts to taxonomy
otu_genus <- otu_sub %>%
  left_join(tax2[, c("FeatureID", "Genus")], by = "FeatureID")

# Sum selected ASVs within each genus
otu_genus <- otu_genus %>%
  select(-FeatureID) %>%
  group_by(Genus) %>%
  summarise(across(everything(), sum), .groups = "drop")

# Convert to matrix with genus as row names
otu_genus_mat <- as.data.frame(otu_genus)
rownames(otu_genus_mat) <- otu_genus_mat$Genus
otu_genus_mat$Genus <- NULL

# =========================
# 6. Convert to relative abundance
# =========================
otu_rel <- sweep(as.matrix(otu_genus_mat), 2, colSums(otu_genus_mat), "/")

# Samples should be rows, genera should be columns
df <- as.data.frame(t(otu_rel))
df$sample_id <- rownames(df)

# =========================
# 7. Add outcome labels
# =========================
df <- df %>%
  left_join(meta_sub %>% select(`sample-id`, treatment_group),
            by = c("sample_id" = "sample-id"))

df <- df %>%
  mutate(group = case_when(
    treatment_group == "unexposed" ~ "Unexposed",
    treatment_group == "exposed" ~ "Exposed"
  ))

# =========================
# 8. Final RF input table
# =========================
df_final <- df %>%
  select(-sample_id, -treatment_group)

# Predictors
predictors <- df_final %>% select(-group)

# Outcome
outcome <- df_final$group %>%
  factor(levels = c("Unexposed", "Exposed"))

if (ncol(predictors) == 0) {
  stop("No predictor columns left after preprocessing.")
}

# =========================
# 9. Cross-validation folds
# =========================
set.seed(421)
folds <- createFolds(outcome, k = 5, list = TRUE)

# =========================
# 10. Hyperparameter grid
# =========================
mtry_vals <- unique(pmin(c(1, 3, 6, 10), ncol(predictors)))

tune_grid <- expand.grid(
  mtry = mtry_vals,
  splitrule = c("gini", "extratrees"),
  min.node.size = c(2, 3, 4)
)

# =========================
# 11. Run RF
# =========================
source("randomforest_functions.R")

rf_model <- run_rf(
  X = predictors,
  y = outcome,
  fold_list = folds,
  hyper = tune_grid,
  rngseed = 421
)

# =========================
# 12. Inspect outputs
# =========================
table(outcome)
dim(predictors)
rf_model$auc_train
rf_model$auc_test
head(rf_model$importance)

# =========================
# 13. Generate plots
# =========================

# -------------------------
# 1. ROC Curve with test only
# -------------------------
roc_obj <- roc(
  response = rf_model$test_labels$true_labels,
  predictor = rf_model$test_labels$predicted_probabilities,
  levels = c("Unexposed", "Exposed")
)

png("RF_ISA_exposed_unexposed_ROC_curve.png", width = 1800, height = 1400, res = 220)
plot(
  roc_obj,
  main = paste0("RF ROC Curve (AUC = ", round(rf_model$auc_test, 3), ")"),
  col = "blue",
  lwd = 3
)
abline(a = 0, b = 1, lty = 2, col = "gray50")
dev.off()

# -------------------------
# 1.2 ROC Curve with test and training
# -------------------------
roc_test <- roc(
  response = rf_model$test_labels$true_labels,
  predictor = rf_model$test_labels$predicted_probabilities,
  levels = c("Unexposed", "Exposed")
)

roc_train <- roc(
  response = rf_model$train_labels$true_labels,
  predictor = rf_model$train_labels$predicted_probabilities,
  levels = c("Unexposed", "Exposed")
)

png("RF_ISA_exposed_unexposed_ROC_curve_both.png", width = 1800, height = 1400, res = 220)

plot(
  roc_test,
  main = paste0(
    "Random Forest ROC Curve\n",
    "Test AUC = ", round(rf_model$auc_test, 3),
    " | Train AUC = ", round(rf_model$auc_train, 3)
  ),
  col = "blue",
  lwd = 3
)

plot(
  roc_train,
  col = "red",
  lwd = 3,
  add = TRUE
)

abline(a = 0, b = 1, lty = 2, col = "gray50")

legend(
  "bottomright",
  legend = c("Test", "Train"),
  col = c("blue", "red"),
  lwd = 3,
  bty = "n"
)

dev.off()

# -------------------------
# 2. Top 15 Important Genera
# -------------------------
top15 <- rf_model$importance %>%
  slice_max(order_by = MeanDecreaseGini, n = min(15, nrow(rf_model$importance)))

p1 <- ggplot(top15, aes(x = reorder(Feature, MeanDecreaseGini), y = MeanDecreaseGini)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Top Important Genera (ISA-filtered RF)",
    x = "Genus",
    y = "Mean Decrease Gini"
  ) +
  theme_minimal(base_size = 14)

ggsave("RF_ISA_exposed_unexposed_top15_genera.png", p1, width = 9, height = 7, dpi = 300)

# -------------------------
# 3. Confusion Matrix
# -------------------------
pred_class <- ifelse(rf_model$test_labels$predicted_probabilities > 0.5, "Exposed", "Unexposed")
pred_class <- factor(pred_class, levels = c("Unexposed", "Exposed"))

cm <- confusionMatrix(pred_class, rf_model$test_labels$true_labels, positive = "Exposed")
print(cm)

cm_df <- as.data.frame(cm$table)

p2 <- ggplot(cm_df, aes(x = Prediction, y = Reference, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = Freq), size = 7) +
  scale_fill_gradient(low = "white", high = "tomato") +
  labs(title = "Confusion Matrix") +
  theme_minimal(base_size = 14)

ggsave("RF_ISA_exposed_unexposed_confusion_matrix.png", p2, width = 6, height = 5, dpi = 300)

# -------------------------
# 4. Predicted Probability Distribution
# -------------------------
p3 <- ggplot(
  rf_model$test_labels,
  aes(x = predicted_probabilities, fill = true_labels)
) +
  geom_histogram(alpha = 0.6, bins = 10, position = "identity") +
  labs(
    title = "Predicted Probability of Exposed",
    x = "Predicted probability",
    y = "Count",
    fill = "True group"
  ) +
  theme_minimal(base_size = 14)

ggsave("RF_ISA_exposed_unexposed_probability_distribution.png", p3, width = 7, height = 5, dpi = 300)

# -------------------------
# 5. Sample Count by Group
# -------------------------
group_counts <- data.frame(
  Group = names(table(outcome)),
  Count = as.numeric(table(outcome))
)

p4 <- ggplot(group_counts, aes(x = Group, y = Count)) +
  geom_col() +
  labs(
    title = "Sample Count by Group",
    x = "Group",
    y = "Number of samples"
  ) +
  theme_minimal(base_size = 14)

ggsave("RF_ISA_exposed_unexposed_group_counts.png", p4, width = 5, height = 4, dpi = 300)

# -------------------------
# 6. Clustered heatmap of top important genera
# -------------------------

# Take top 12 important genera from RF
top_features <- rf_model$importance %>%
  slice_max(order_by = MeanDecreaseGini, n = 12) %>%
  pull(Feature)

# Get those features from predictors
heatmap_mat <- predictors[, top_features, drop = FALSE]

# rows = genera, columns = samples
heatmap_mat <- t(as.matrix(heatmap_mat))

# log transform to reduce extreme sparsity
heatmap_mat <- log10(heatmap_mat + 1e-4)

# row-scale for better contrast
heatmap_mat <- t(scale(t(heatmap_mat)))
heatmap_mat[is.na(heatmap_mat)] <- 0

# sample colors
col_side <- ifelse(outcome == "Exposed", "#E76F51", "#4C78A8")

png("RF_clustered_heatmap_top12_clean.png", width = 2200, height = 1600, res = 240)

heatmap(
  heatmap_mat,
  ColSideColors = col_side,
  scale = "none",
  col = colorRampPalette(c("#313695", "#74ADD1", "white", "#FDAE61", "#A50026"))(120),
  margins = c(8, 12),
  labCol = FALSE,
  cexRow = 1.2,
  main = "Clustered Heatmap of Top 12 Important Genera"
)

legend(
  "topright",
  legend = c("Unexposed", "Exposed"),
  fill = c("#4C78A8", "#E76F51"),
  border = NA,
  bty = "n",
  cex = 1.1
)

dev.off()


# -------------------------
# 7. Grouped heatmap of top important genera
# -------------------------

# Take top 12 important genera from RF
top_features <- rf_model$importance %>%
  slice_max(order_by = MeanDecreaseGini, n = 12) %>%
  pull(Feature)

# Get those features from predictors
heatmap_mat <- predictors[, top_features, drop = FALSE]

# rows = genera, columns = samples
heatmap_mat <- t(as.matrix(heatmap_mat))

# log transform + row scaling
heatmap_mat <- log10(heatmap_mat + 1e-4)
heatmap_mat <- t(scale(t(heatmap_mat)))
heatmap_mat[is.na(heatmap_mat)] <- 0

# Order columns by group
ord <- order(outcome)
heatmap_mat <- heatmap_mat[, ord, drop = FALSE]
outcome_ord <- outcome[ord]

# Top color bar for groups
col_side <- ifelse(outcome_ord == "Exposed", "#E76F51", "#4C78A8")

png("RF_grouped_heatmap_top12.png", width = 2200, height = 1600, res = 240)

par(mar = c(8, 12, 0, 4))

heatmap(
  heatmap_mat,
  Colv = NA,
  Rowv = TRUE,
  ColSideColors = col_side,
  scale = "none",
  col = colorRampPalette(c("#313695", "#74ADD1", "white", "#FDAE61", "#A50026"))(120),
  margins = c(8, 12),
  labCol = FALSE,
  cexRow = 1.1,
  main = "Grouped Heatmap of Top 12 Important Genera"
)

legend(
  "topright",
  inset = c(-0.01, 0.00),
  legend = c("Unexposed", "Exposed"),
  fill = c("#4C78A8", "#E76F51"),
  border = NA,
  bty = "n",
  cex = 0.95
)

dev.off()

