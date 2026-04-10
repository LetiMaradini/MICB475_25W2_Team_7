# Aim_3 – Indicator Taxa Analysis
**Mar 17th, 2026**

---

## Purpose
To identify microbial taxa significantly associated with microgravity exposure and experimental conditions in a murine gut microbiome dataset using Indicator Species Analysis. This approach detects taxa that are both **specific** and **frequent** within defined groups, making it well-suited for ecological interpretation and biomarker discovery.

ISA Exposure Model: simplifies samples into exposure categories to identify robust indicator taxa and support machine learning applications  

---

## Comparison groups

Samples were grouped into two categories:

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

## Visualization and Post-processing

- Extracted significant taxa (`p < 0.05`)
- Annotated taxa using taxonomy (Genus, Family, Phylum)
- Filtered for stat values > 0.7
- Removed duplicate genus/families
- Generated:
  - Table of ISA values 

---
## Code

[Indicator_species code](../Scripts/Indicator_taxa_final.R)

## Output files

### Indicator analysis outputs
- `indicator_values_ISA1.txt`  
  - Significant indicator taxa for **microgravity exposure model (ISA-1)**  
  - Includes IndVal scores, associated group, and p-values
 
---

## Results
- Indicator taxa analysis identified microbial taxa significantly associated with:
  - **Exposed** (microgravity or recovery)
  - **Unexposed** (terrestrial/no exposure)
- Taxa identified in this model represent **strong ecological indicators of microgravity exposure status**
- Top 3 stat values: Lachnoclostridium, Colidextribacter, Turicibacter (genus)

---

## Discussion

- Indicator Species Analysis (ISA) differs from differential abundance methods (e.g., DESeq2) in that it:
  - Focuses on **specificity and consistency of association**
  - Identifies taxa that are **characteristic of particular environments or conditions**

- Produces **strong candidate biomarkers**
- Well-suited for:
  - Machine learning classification
  - Predictive modeling of microgravity exposure

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
  - Use ISA taxa as features for classification models (Random Forest)
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



