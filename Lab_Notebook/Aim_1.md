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
- Alpha test: Wilcoxon rank-sum
- Alpha metrics: Shannon, Faith's PD
- Beta test: PERMANOVA & betadisper
- Beta distances: Bray-Curtis, wUniFrac

#### Longitudinal 

- comparison between different time-points within same cohort for GC-LAR & F-LAR
- create dataframe for cohort
- Alpha tests: Kruskal & Dunn
- Alpha metrics: Shannon, Faith's PD
- Beta tests: PERMANOVA
- Beta distances: Bray-Curtis & wUniFrac
- plot using: ```ggplot()```

#### Pairwise Time-points

- comparison between the same time-points in GC-LAR & F-LAR
- function: ```run_test_2group()```
- for timepoints tp in W0_Pre, W4.5_Pre, W9_Pre, W9_Post
- Alpha test: 3 Wilcoxon rank-sum
- Alpha metrics: Shannon, Faith's PD
- Beta test: 3 pairwise PERMANOVA
- Beta distances: Bray-Curtis, wUniFrac

#### Recovery vs Flight

- comparison for samples F-ISST vs F-LAR both at 9 weeks after euthanasia
- function: ```run_test_2group()```
- input F-ISST & F-LAR for W9_Post tp
- Alpha test: Wilcoxon rank-sum
- Alpha metrics: Shannon, Faith's PD
- Beta test: PERMANOVA + betadisper
- Beta distances: Bray-Curtis, wUniFrac

#### Ground Controls 

- comparison for ground controls GC-LAR vs GC-ISST at 9 weeks after euthanasia
- function: ```run_test_2group()```
- input GC-ISST & GC-LAR for W9_Post tp
- Alpha test: Wilcoxon rank-sum
- Alpha metrics: Shannon, Faith's PD
- Beta test: PERMANOVA + betadisper
- Beta distances: Bray-Curtis, wUniFrac

#### Basal vs Control

- comparison for basal (combined BL + W0) vs control (W9)
- function: ```run_test_2group()```
- input GC-LAR W9_Post vs Basal (BL + W0_Pre)
- Alpha test: Wilcoxon rank-sum
- Alpha metrics: Shannon, Faith's PD
- Beta test: PERMANOVA + betadisper
- Beta distances: Bray-Curtis

### Saved outputs

2 files generated
- analysis metric values in `Satistical_Analysis_Results.txt`
- diversity plots in `Biological_Visualizations.pdf`



---
## Code

[Diversity analysis](../Scripts/Diversity_analysis) 

 

---

## Results

### Microgravity (F-ISST) vs Control (GC-ISST)

- compared the mice remaining on ISS against ground control at the end of the study
- Alpha diversity: significant difference in *richness* (Faith's PD, p = 0.00684) eventhough the difference in *evenness* appeared insignificant (Shannon, p = 0.06301)
- Beta diversity: significant difference in the types of bacteria present (Bray-Curtis, p = 0.005), but remains insignificant when considering evolutionary relationships (wUniFrac, p = 0.181)
- dispersion: relatively similar levels of variability (betadisper, p = 0.055), further supporting the PERMANOVA result

### Longitudinal (GC-LAR & F-LAR)

- compared how the microbiome changed over time within the same group of mice
- GC-LAR: significant shifts in alpha diversity (p = 0.01484) most notably between week 0 and week 4.5 (p-adj = 0.014). beta diversity revealed significant shift in community composition over time (p = 0.001)
- F-LAR: significant alpha diversity shift (p = 0.0119), most notably between week 0 and week 4.5 also (p-adj = 0.0109). beta diversity revealed a significant shift in commmunity composition also over time (Bray-Curtis R^2 = 32.4%, and wUniFrac R^2 = 28.3%)

### Pairwise Time-points

- compared differences at identical time intervals between F-LAR (flight) & GC-LAR (ground) cohorts
- weeks 0 & 4.5: no significant differences, both initially responded siml=ilarly to their respective environments
- week 9 (pre-euthanasia): no significant differences, but the significant `betadisper` result (p = 0.014) revealed that F-LAR started to show more individual variation
- week 9 (post-euthanasia): simila alpha diversity, but significant beta diversity (Bray-Curtis p = 0.032, wUniFrac p = 0.05) suggesting that the composition of the F-LAR group significantly diverged from the GC-LAR ground group

### Recovery vs Flight

- compared mice euthanized whilst on the ISS (F-ISST) vs those euthanized later after returned to Earth and recovered (F-LAR)
- diversity metrics: every metric significantly different (Shannon p = 0.00004, Beta p = 0.001). significant R^2 values (27.6% and 39.3%) suggest that return to Earth and subsequent recovery time created diverging microbiome profiles

### Ground controls

- compared 2 different ground control groups to see if they remained similar independent of microgravity factors
- diversity metrics: significant differences in alpha (Shannon p = 0.00032) and beta (p = 0.001) suggesting that different housing/handling protocol conditions led to microbiome differences independent of flight/microgravity conditions

### Basal vs Control 

- compared the final state of the ground control group (GC-LAR at W9) vs original baseline (Basal)
- alpha diversity: significant differences in *richness* (Faith's PD, p = 0.025)
- beta diversity: significant differences in *composition* (Bray-Curtis, p = 0.001)
- diversity results suggest a natural drift in the microbiome over the 9 week period even without the influence of microgravity

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
