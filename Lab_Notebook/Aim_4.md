# P03 -  Differential Abundance (DESeq2)

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
1. R & Rstudio (packages: phyloseq, tidyverse, DESeq2)
2. micro_final.RData (phyloseq object unrarefied from P03)

## Method:
### Data Loading & Preprocessing
1. Loaded phyloseq object "micro_final" from "micro_final.RData".
2. Factored treatment_group into levels: BL, F-LAR, GC-LAR, F-ISST, GC-ISST.
3. Factored Time into levels: Launch plus 0, Launch plus 4.5, Launch plus 9.
4. Created sample_type_group by appending -Fresh/-Necropsy to F-LAR/GC-LAR samples at Launch plus 9; factored with appropriate levels.
5. Applied +1 pseudocount transformation to all sample counts, generating micro_plus1 phyloseq object.
### DESeq2 Differential Abundance Analysis
1. Defined run_deseq_treatment helper: subset by sample_type_group, convert to DESeqDataSet (~sample_type_group design), run DESeq(), extract tidy contrast results.
2. Microgravity (Launch plus 9): Subset F-ISST vs GC-ISST; DESeq2 contrast via helper.
3. F-LAR longitudinal: Subset F-LAR samples; factored host_id; created time_group (Week0, Week4.5, Week9_Fresh, Week9_Necropsy); DESeq2 (~host_id + time_group); extracted 6 contrasts (4.5vs0, 9freshvs0, 9necrovs0, 9freshvs4.5, 9necrovs4.5, 9freshvs9necro).
4. GC-LAR longitudinal: Identical workflow to F-LAR; extracted same 6 contrasts.
5. LAR cross-cohort timepoints: Mapped 4 comparisons (0w/4.5w F-LAR vs GC-LAR; 9w Fresh/Necropsy); filtered samples (skip <4 samples); DESeq2 via helper; combined results as res_LAR_pairs.
6. Return vs flight (Week 9): Subset F-ISST vs F-LAR-Necropsy; DESeq2 contrast.
7. Controls Week 9: Subset GC-LAR-Necropsy vs GC-ISST; DESeq2 contrast.
8. Basal control (Week 0): Subset GC-LAR vs BL; DESeq2 contrast.
### Visualization & Output Generation
1. Defined volcano_plot helper: flag significant ASVs (padj < 0.05, |log2FC| > 2); plot log2FC vs -log10(padj); save PNGs for all 20 comparisons.
2. Significance summary: Compiled 20 result objects (slicing res_LAR_pairs for 0w/4.5w); rowwise computed tested/up/down/non-sig ASVs (padj < 0.05, |log2FC| > 2); saved "ASVs_count_summary.csv".
3. Significant ASV barplots/tables: Filtered 20 comparisons for >0 strict significant ASVs; per comparison: extracted sig ASVs, joined with tax_table from micro_final, sorted by log2FC, unique-ified Genus names; generated log2FC barplots by Genus with lfcSE error bars (90° x-axis rotation); saved PNGs and per-comparison CSV tables.

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
