# Aim_1 – Diversity Analysis
**Mar 10th, 2026**

---

## Purpose
To investigate alpha & beta diversity metrics both within and between mice gut microbiome samples. Preliminary analysis was used to see whether statistical significance warrented the need for additional investigation between samples. Diversity metrics were used to see whether sampling conditions, i.e exposure to microgravity and/or spatial isolation, would impact mice gut microbiome communities in terms of abundance, composition, and environment-dependent fluctuations. 

Two main analyses were conducted:
- **Sample alpha diversity calculations (to check for statistical significance):** looked at general changes in diversity to see whether significance would justify further investigation
- **Alpha & Beta diversity (sampling group comparisons):** compared diversity metrics across different groups and time points in order to extract relevant information

---

## Comparison groups (CONTINUE EDIT FROM HERE)

### Microgravity Exposure Model

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

### Non-exposed to Microgravity (Condition × Time Model i.e 9 groups)

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
  - `picante`
  - `ape`
  - `FSA`

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
- Construction of sample metadata including:
  - `Time` (i.e "Launch plus 0" ~ "W0_Pre)
  - `Sample_State` ("Pre" or "Post")
- Format taxonomy via OTU and taxonomy tables  
- Build Phyloseq & pre-process (rarify) to construct `phyloseq` object

---

### Filtering and Rarefaction

- Removed low-abundance taxa (total counts ≤ 5)  
- Removed low-depth samples (< 100 reads)  
- Rarefied all samples to equal depth (64,411 reads) using:

### Group comparison statistical analysis & plotting

#### Build function for comparison between 2 sample groups

```r
run_test_2group()
run_test_2group <- function(ps_obj, targets, group_col, title, use_wunifrac = TRUE)
```

#### Microgravity vs Control

- comparing F-ISST vs GC-ISST sample groups
- function: ```run_test_2group()```
- inputs: ps_rare, Sample_state, F-ISST, GC-ISST, Group, F-ISST x GC-ISST

#### Longitudinal 

- comparison within sample group cohort for GC-LAR & F-LAR
- create dataframe for cohort
- run Alpha tests (Kruskal & Dunn)
- run Beta tests (bray curtis & weighted unifrac)
- plot using: ```ggplot()```

#### Pairwise Time-points

- comparison for time-points in GC-LAR & F-LAR
- function: ```run_test_2group()```
- for timepoints tp in W0_Pre, W4.5_Pre, W9_Pre, W9_Post

#### Recovery vs Flight

- comparison for samples F-ISST vs F-LAR
- function: ```run_test_2group()```
- input F-ISST & F-LAR for W9_Post tp

#### Ground Controls 

- comparison for ground control in both samples i.e GC-LAR vs GC-ISST
- function: ```run_test_2group()```
- input GC-ISST & GC-LAR for W9_Post tp

#### Basal vs Control

- comparison for basal (combined BL + W0) vs control (W9)
- function: ```run_test_2group()```
- input GC-LAR W9_Post vs Basal (BL + W0_Pre)

### Saved outputs

2 files generated
- analysis metric values in `Satistical_Analysis_Results.txt`
- diversity plots in `Biological_Visualizations.pdf`



---
## Code

[Diversity analysis](../Scripts/Diversity analysis) 

 

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



---

## Discussion

- Indicator Species Analysis (ISA) differs from differential abundance methods (e.g., DESeq2) in that it:
  - Focuses on **specificity and consistency of association**
  - Identifies taxa that are **characteristic of particular environments or conditions**

---




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
  

---

- **Ecological interpretation**
