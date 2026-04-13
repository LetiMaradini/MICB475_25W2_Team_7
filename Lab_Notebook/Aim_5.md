# Aim_5 - Random Forest Analysis

## Purpose

The purpose of this analysis was to determine whether gut microbiome composition could be used to classify samples as exposed or unexposed using a Random Forest (RF) model. This analysis was designed to test whether microbiome profiles contain a reproducible multivariate signature associated with exposure status. Because the dataset contains many ASVs, Indicator Species Analysis (ISA1) results were first used as a feature selection step to reduce dimensionality and retain only the most informative taxa for classification.

## Rationale

Random Forest was chosen because it is well-suited for microbiome data, where relationships between taxa and outcomes are likely nonlinear and involve interactions among many features. Compared with simple univariate approaches, RF can better capture multivariable patterns across taxa. To avoid including too many weakly informative predictors, ISA1 results were used to pre-filter the feature set before model training.

## Input files and software

This analysis was performed in R/RStudio using the following files:
- micro_final_otu.csv
- taxonomy.tsv
- microgravity_metadata.tsv
- indicator_values_ISA1.txt
- randomforest_functions.R

Main R packages used:
- tidyverse
- caret
- randomForest
- ranger
- pROC
- boot

## Methods

The ISA1 output was first parsed to identify taxa associated with the three groups: control, unexposed, and exposed. To focus the RF analysis specifically on the exposed vs unexposed comparison, taxa assigned to the control group were excluded. Among the remaining taxa, only those with an IndVal stat greater than 0.7 were retained. This threshold was used to keep relatively strong indicator taxa while reducing noise and feature dimensionality.

Metadata did not contain a pre-existing exposed/unexposed grouping variable, so one was created from the Condition column:
- Space Flight = exposed
- Ground Control = unexposed
- Basal Control = control

Only exposed and unexposed samples were retained for RF analysis. The OTU table was then filtered to keep only the ASVs selected from the ISA1 table, and sample order was aligned between the OTU table and metadata.

Taxonomy strings were processed to extract the genus assignment for each ASV. Selected ASVs were then collapsed to the genus level by summing counts across ASVs assigned to the same genus. Any missing genus assignment was labeled as Unassigned. The genus-level counts were converted to relative abundance, and the resulting table was transposed so that rows represented samples and columns represented genera.

The final RF input dataset contained:
- 100 samples total
- 50 exposed
- 50 unexposed
- 20 genus-level predictors

The response variable was a binary factor with levels:
- Unexposed
- Exposed

Five-fold cross-validation was used for model evaluation. A hyperparameter grid was tested with:
- mtry = 1, 3, 6, 10 (capped to available predictor number)
- splitrule = gini, extratrees
- min.node.size = 2, 3, 4

The Random Forest model was trained using the helper functions in randomforest_functions.R. Model performance was assessed using:
- training AUC
- test AUC
- confusion matrix metrics
- feature importance

## Code
[RF code](../Scripts/RFresult_exposed_unexposed.R)

## Results

The ISA-filtered RF model showed strong classification performance for exposed vs unexposed samples. The model achieved:
- Training AUC = 0.871
- Test AUC = 0.924

These values indicate that the microbiome profiles contained a strong signal allowing separation of exposed and unexposed samples. The confusion matrix showed:
- Accuracy = 0.84
- Sensitivity = 0.82
- Specificity = 0.86
- Balanced accuracy = 0.84

These results indicate that the model performed well for both classes and did not appear strongly biased toward one group.

The most important genera contributing to classification included:
- Escherichia-Shigella
- Lachnoclostridium
- Akkermansia
- Ligilactobacillus
- Turicibacter
- Roseburia

Additional visualizations were generated to summarize the model, including:
- ROC curve
- confusion matrix
- top important genera bar plot
- grouped heatmap of top important genera
- predicted probability distribution
- sample count by group

The grouped heatmap showed that the top important genera displayed different abundance patterns across exposed and unexposed samples, supporting the interpretation that the RF model was using biologically meaningful microbial differences rather than random variation.

## Interpretation

Overall, the RF model suggests that exposed and unexposed samples can be distinguished with relatively high accuracy based on microbiome composition after ISA-based feature selection. This supports the idea that exposure status is associated with a detectable microbial signature. The taxa highlighted by the feature importance analysis are candidate contributors to that signature and may be useful for interpretation alongside the ISA and differential abundance results.

## Discussion

The test AUC was slightly higher than the training AUC. This can happen in small datasets due to variation in train-test splits and differences in subset difficulty. In this case, the pattern does not suggest obvious overfitting, since overfitting would more commonly appear as a much higher training AUC than test AUC. However, because the sample size is limited, the result should still be interpreted with caution.

A major strength of this approach was the use of ISA1 as a feature selection step before RF. This reduced dimensionality and restricted the model to taxa already identified as associated with exposed or unexposed groups, likely improving interpretability. Collapsing features to the genus level also reduced sparsity and made the final model easier to explain biologically.

One limitation is that the exposed and unexposed categories combine multiple biological contexts under broad labels. This makes the model useful for overall exposed vs unexposed classification, but less specific for identifying which particular cohort or timepoint differences are driving the separation. In addition, the RF model identifies predictive taxa, but this does not establish causation.

## Conclusion

The ISA-filtered Random Forest analysis successfully classified exposed and unexposed microbiome samples with high performance. The results support the presence of a reproducible microbial signature associated with exposure status and identified several genera as important contributors to this classification. This analysis provides a predictive complement to the diversity, ISA, and differential abundance analyses performed elsewhere in the project.
	
