# Aim 3 – Indicator Taxa Analysis
**Mar 17th, 2026**

---

## Purpose
To identify microbial taxa significantly associated with microgravity exposure and experimental conditions in a murine gut microbiome dataset using Indicator Species Analysis. This approach detects taxa that are both **specific** and **frequent** within defined groups, making it well-suited for ecological interpretation and biomarker discovery.

Two complementary analyses were conducted:
- **ISA-1 (Exposure model):** simplifies samples into exposure categories to identify robust indicator taxa and support machine learning applications  
- **ISA-2 (Condition × Time model):** retains experimental resolution to detect fine-scale ecological and temporal patterns  

---

## Comparison groups

### ISA-1: Microgravity Exposure Model

Samples were grouped into three biologically meaningful categories:

**Exposed (microgravity or recovery exposure):**
- F-ISST Week 9  
- F-LAR Week 4.5  
- F-LAR Week 9  

**Unexposed (no microgravity exposure):**
- GC-ISST Week 9  
- GC-LAR Week 0  
- GC-LAR Week 4.5  
- GC-LAR Week 9  
- F-LAR Week 0  

**Control (baseline):**
- BL Week 0  

---

### ISA-2: Condition × Time Model (9 groups)

Samples were grouped into nine treatment × timepoint categories:

- F-ISST Week 9  
- GC-ISST Week 9  
- F-LAR Week 0  
- F-LAR Week 4.5  
- F-LAR Week 9  
- GC-LAR Week 0  
- GC-LAR Week 4.5  
- GC-LAR Week 9  
- BL Week 0  

---

## Materials

- **Software:** R, RStudio  
- **R Packages:**
  - `phyloseq`
  - `tidyverse`
  - `vegan`
  - `indicspecies`
  - `ape`
  - `stringr`

- **Input Files:**
  - `microgravity_metadata.tsv`
  - `feature-table.txt`
  - `taxonomy.tsv`
  - `tree.nwk`

---

## Method
Uses phyloseq object created in [PO3](../Lab_Notebook/P03.md)
### Data Import and Phyloseq Construction

- Imported OTU table, taxonomy, metadata, and phylogenetic tree  
- Converted OTU table into matrix format and assigned taxa IDs  
- Parsed taxonomy into hierarchical ranks (Domain → Species)  
- Created sample metadata including:
  - `treatment_group` (derived from host ID)
  - `Week` (0, 4.5, 9)

- Matched taxa between OTU and taxonomy tables  
- Constructed a unified `phyloseq` object (`micro`)

---

### Filtering and Rarefaction

- Removed low-abundance taxa (total counts ≤ 5)  
- Removed low-depth samples (< 100 reads)  
- Rarefied all samples to equal depth (64,411 reads) using:

```r
rarefy_even_depth()
- Generated:
  - micro_final (filtered dataset)
  - micro_rare (rarefied dataset used for ISA)

# Indicator Species Analysis

## ISA-1: Exposure-Based Grouping

Created grouping variable `isa1_group` with levels:
- `control`
- `unexposed`
- `exposed`

Assigned samples based on treatment group and timepoint, and removed samples not matching defined categories.

Subset phyloseq object:

```r
ps_isa1 <- prune_samples(keep1, ps)
```

Extracted OTU table (samples × taxa format) and ran Indicator Species Analysis:

```r
indicator_multipatt_1 <- multipatt(otu_isa1, isa1_group, duleg = TRUE)
```

Saved output:

```r
indicator_output_1 <- capture.output(summary(indicator_multipatt_1, indvalcomp = TRUE))
writeLines(indicator_output_1, "indicator_values_ISA1.txt")
```

---

## ISA-2: Condition × Time Grouping

Created combined grouping variable:

```r
CondWeek = paste(treatment_group, Week, sep = "_")
```

Defined a nine-category grouping factor (`isa2_group`) and filtered samples to retain only defined groups.

Subset phyloseq object:

```r
ps_isa2 <- prune_samples(keep2, ps)
```

Extracted OTU table (samples × taxa) and ran Indicator Species Analysis:

```r
indicator_multipatt_2 <- multipatt(otu_isa2, isa2_group, duleg = TRUE)
```

Saved output:

```r
indicator_output_2 <- capture.output(summary(indicator_multipatt_2, indvalcomp = TRUE))
writeLines(indicator_output_2, "indicator_values_ISA2.txt")
```

---

## Visualization and Post-processing

- Extracted significant taxa (`p < 0.05`)
- Annotated taxa using taxonomy (Genus, Family, Phylum)
- Generated:
  - Heatmaps of top indicator taxa
  - Barplots of top taxa (IndVal ranking)
  - Alpha diversity plots (Shannon, Simpson)
  - Beta diversity (PCoA, Bray-Curtis)

---
## Code

[Indicator_species code](../Scripts/Indicator_taxa_final.R)

## Output files

### Indicator analysis outputs
- `indicator_values_ISA1.txt`  
  - Significant indicator taxa for **microgravity exposure model (ISA-1)**  
  - Includes IndVal scores, associated group, and p-values  

- `indicator_values_ISA2.txt`  
  - Significant indicator taxa for **condition × time model (ISA-2)**  
  - Higher-resolution associations across all 9 experimental groups  

---

### Visualization outputs
- `ISA1_heatmap.png`  
  - Heatmap of top indicator taxa for ISA-1 (exposure groups)

- `ISA2_heatmap.png`  
  - Heatmap of top indicator taxa for ISA-2 (9 condition × time groups)

- `ISA1_top15.png`  
  - Barplot of top 15 indicator taxa ranked by IndVal (ISA-1)

- `diversity_exposure.png`  
  - Alpha diversity (Shannon, Simpson) across exposure groups

- `pcoa_exposure.png`  
  - Beta diversity (Bray-Curtis PCoA) colored by exposure group  

---

### Tabular outputs
- `ISA1_significant_157taxa.csv`  
  - Table of significant taxa (p < 0.05) from ISA-1  
  - Includes taxonomy (Genus/Family/Phylum) and IndVal scores  

- `ISA2_significant_177taxa.csv`  
  - Table of significant taxa (p < 0.05) from ISA-2  
  - Provides higher-resolution taxon-group associations  

---

## Results
See [Indicator Species
### ISA-1 (Microgravity Exposure Model)
- Indicator taxa analysis identified microbial taxa significantly associated with:
  - **Exposed** (microgravity or recovery)
  - **Unexposed** (terrestrial/no exposure)
  - **Control** (baseline)

- Taxa identified in this model represent **strong ecological indicators of microgravity exposure status**

---

### ISA-2 (Condition × Time Model)
- Indicator taxa were identified for all **9 treatment × timepoint combinations**

- This model captured:
  - **Temporal dynamics** of microbiome shifts
  - Differences between:
    - Flight vs ground control
    - Early vs late exposure
    - Recovery vs active exposure

- ISA-2 revealed **more nuanced and condition-specific microbial associations** compared to ISA-1

- However:
  - Increased group resolution results in **lower statistical power per group**
  - More variability due to smaller sample sizes per category

---

## Discussion

- Indicator Species Analysis (ISA) differs from differential abundance methods (e.g., DESeq2) in that it:
  - Focuses on **specificity and consistency of association**
  - Identifies taxa that are **characteristic of particular environments or conditions**

---

### ISA-1
- Produces **strong candidate biomarkers**
- Well-suited for:
  - Machine learning classification
  - Predictive modeling of microgravity exposure

---

### ISA-2
- Captures:
  - Time-dependent microbial shifts
  - Subtle differences between experimental cohorts
- Enables ecological interpretation of:
  - Exposure duration
  - Recovery effects
  - Environmental vs microgravity influences

---

### Limitations
- ISA assumes taxa are independent, which may not reflect ecological interactions
- Rare taxa may be excluded due to filtering and rarefaction
- Uneven sample sizes across groups (especially in ISA-2) can reduce statistical power
- Taxonomic resolution limitations (many taxa unresolved beyond higher ranks)

---

### Integration with other analyses
- ISA complements:
  - **DESeq2** (magnitude of change)
  - **Alpha diversity** (within-sample diversity)
  - **Beta diversity** (community structure differences)

- Combining these approaches provides a **more complete understanding of microbiome responses**

---

## Future direction

- **Machine learning applications**
  - Use ISA-1 taxa as features for classification models (Random Forest)
  - Predict exposure status based on microbiome composition

---

- **Functional analysis**
  - Map indicator taxa to:
    - Metabolic pathways
    - Functional gene categories
  - Determine ecological roles of key taxa

---

- **Cross-method validation**
  - Compare ISA results with DESeq2 outputs
  - Identify taxa that are both:
    - Differentially abundant
    - Strong ecological indicators

---

- **Ecological interpretation**
  - Investigate how microgravity influences:
    - Microbial community assembly
    - Stability and resilience of the microbiome



