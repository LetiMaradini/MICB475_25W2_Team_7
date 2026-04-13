# Aim_4 -  Differential Abundance (DESeq2)

Mar 17th, 2026

## Purpose: 
To identify differentially abundant ASVs between flight (F-LAR/F-ISST) and ground control (GC-LAR/GC-ISST) microbiome samples across mission timepoints (Week 0, 4.5, 9) and sample types (fresh vs necropsy on week 9 for LAR) using DESeq2 analysis. Generate volcano plots, significant ASV tables (padj<0.05, |LFC|>2), barplots of log2 fold changes, and summary statistics of significant ASVs per comparison for figures and tables.
### Comparison groups:
   1. Difference between microgravity and control: F-ISST x GC-ISST
   2. Difference between different time-points within the same cohort (0 x 4.5 x 9 weeks before euthanasia x 9 weeks after euthanasia): GC-LAR and F-LAR
   3. Difference between the same time-points in the cohorts (0 x 0, 4.5 x 4.5, 9 weeks before euthanasia x 9 weeks before euthanasia, and 9 weeks after euthanasia x 9 weeks after euthanasia): GC-LAR x F-LAR
   4. Difference between return (time for recovery) and flight (both at 9 weeks after euthanasia): F-ISST x F-LAR
   5. Compare control and basal (week 0 not dividing these between pre/post euthanasia because both of them have these samples so it would already be consistent): GC-LAR x Basal

## Material: 
1. R & Rstudio (packages: phyloseq, tidyverse, DESeq2, dplyr, purrr, tibble, ggplot2, readr, patchwork)
2. [micro_final.RData (phyloseq object unrarefied from P03)](../Backup_Files/RStudio/Phyloseq_object/micro_final.RData)

## Method:
### Data Loading & Preprocessing
1. Loaded phyloseq object micro_final from micro_final.RData.
2. Factored treatment_group into levels BL, F-LAR, GC-LAR, F-ISST, and GC-ISST.
3. Factored Time into levels Launch plus 0, Launch plus 4.5, and Launch plus 9.
4. Created sample_type_group by appending -Fresh/-Necropsy to F-LAR/GC-LAR samples at Launch plus 9; factored with appropriate levels.
5. Applied a +1 pseudocount transformation to all sample counts, generating the micro_plus1 phyloseq object.
### DESeq2 Differential Abundance Analysis
1. Defined the run_deseq_treatment helper to subset by sample_type_group, convert to a DESeq2 dataset with a ~ sample_type_group design, run DESeq(), and extract tidy contrast results.
2. For microgravity at Launch plus 9, analyzed F-ISST versus GC-ISST using the helper.
3. For F-LAR longitudinal analysis, subset F-LAR samples, factored host_id, created time_group (Week0, Week4.5, Week9_Fresh, Week9_Necropsy), fit DESeq2 with ~ host_id + time_group, and extracted six contrasts (4.5vs0, 9freshvs0, 9necrovs0, 9freshvs4.5, 9necrovs4.5, 9freshvs9necro).
4. The same workflow was used for GC-LAR longitudinal analysis.
5. For cross-cohort LAR comparisons, mapped four contrasts (Launch plus 0 and Launch plus 4.5 for F-LAR vs GC-LAR, plus Launch plus 9 Fresh and Necropsy), filtered samples with a minimum of four samples, ran DESeq2 via the helper, and combined the results as res_LAR_pairs.
6. Additional contrasts included F-ISST versus F-LAR-Necropsy at week 9, GC-LAR-Necropsy versus GC-ISST at week 9, and GC-LAR versus BL at week 0.

### Visualization & Output Generation
1. Defined volcano_plot to flag significant ASVs using padj < 0.05 and abs(log2FoldChange) > 2, plot log2 fold change versus -log10(padj), and save PNGs for all 20 comparisons.
2. Compiled a 20-row significance summary from the result objects, including the two LAR time-point subsets, and calculated tested, upregulated, downregulated, and non-significant ASVs using the same strict cutoff; saved the table as ASVs_count_summary.csv.
3. Generated significant ASV barplots and tables for all 20 comparisons by extracting strict-significant ASVs, joining taxonomy from micro_final, sorting by log2 fold change, making genus names unique, and plotting genus-level log2 fold change with lfcSE error bars and rotated x-axis labels.
4. Saved the barplots as PNGs and the per-comparison tables as CSV files.
5. Created a publication ready figure using:
   - vol_plot_F-ISST_vs_GC-ISST_9.png
   - vol_plot_F-LAR_vs_GC-LAR_0.png
   - vol_plot_F-LAR_vs_GC-LAR_4.5.png
   - vol_plot_F-LAR-Fresh_vs_GC-LAR-Fresh.png
   - vol_plot_F-LAR-Necropsy_vs_GC-LAR-Necropsy.png

## Code: 
[DESeq code](../Scripts/DESeq.R)

## Output files:
1. [Volcano plots for comparison groups](../Images/RStudio/Aim_4/Volcano_plots) 
2. [ASVs_count_summary.csv (table of 20 comparisons with columns tested, up, down, and ns ASVs)](../Images/RStudio/Aim_4/ASVs_count_summary.csv)
3. [Significant ASVs in barplots per comparison group](../Images/RStudio/Aim_4/Bar_plots_sig_ASVs+Lists)
4. [Tables of significant ASVs per comparison group](../Images/RStudio/Aim_4/Bar_plots_sig_ASVs+Lists)

## Results: 
Check output files

## Discussion:
N/A

## Future direction:
Write figure legends, results and methods for this section
