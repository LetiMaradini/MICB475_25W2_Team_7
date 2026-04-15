library(phyloseq)
library(ape) # importing trees
library(tidyverse)
library(vegan)
library("microbiome")
library("ggVennDiagram")
library("ggplot2")

#### Load data ####
metafp <- "microgravity_export/microgravity_metadata.tsv"
meta <- read_delim(metafp, delim="\t")

otufp <- "microgravity_export/feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

taxfp <- "microgravity_export/taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "microgravity_export/tree.nwk"
phylotree <- read.tree(phylotreefp)

#### Format OTU table ####
# Save everything except first column (OTU ID) into a matrix
otu_mat <- as.matrix(otu[,-1])
# Make first column (#OTU ID) the rownames of the new matrix
rownames(otu_mat) <- otu$`#OTU ID`
# Use the "otu_table" function to make an OTU table
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 

#### Format sample metadata ####
# Creating sample_type column
meta <- meta %>%
  mutate(sample_type = str_extract(`sample-id`, "(?i)(Necropsy|Fresh)"))
# Creating  treatment_group column
meta <- meta %>%
  mutate(
    treatment_group = str_remove_all(host_id, "\\d+")
  )
# Save everything except sampleid as new data frame
samp_df <- as.data.frame(meta[,-1])
# Make sample names the rownames
rownames(samp_df)<- meta$'sample-id'
# Make phyloseq sample data with sample_data() function
SAMP <- sample_data(samp_df)

#### Formatting taxonomy ####
# Convert taxon strings to a table with separate taxa rank columns
tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() # Saving as a matrix
# Save everything except feature IDs 
tax_mat <- tax_mat[,-1]
# Make sampleids the rownames
rownames(tax_mat) <- tax$`Feature ID`
# Make taxa table
TAX <- tax_table(tax_mat)

#### Create phyloseq object ####
# Merge all into a phyloseq object
micro <- phyloseq(OTU, SAMP, TAX, phylotree)
# View the components of phyloseq object with the following
otu_table(micro)
sample_data(micro)
tax_table(micro)
phy_tree(micro)

######### ANALYZE ##########
# Did not remove non-bacterial sequences, because did it before in QIIME (chloroplast & mitochondria removed used table-no-mitochondria-no-chloroplast.qza) and wanted to keep other microbes (archaea etc)
# Remove ASVs that have less than 5 counts total
micro_filt_nolow <- filter_taxa(micro, function(x) sum(x)>5, prune = TRUE)
# Remove samples with less than 100 reads
micro_final <- prune_samples(sample_sums(micro_filt_nolow)>100, micro_filt_nolow)

# Convert phyloseq object to relative abundance
mt_rela <- transform_sample_counts(micro, fun=function(x) x/sum(x)) # Determine relative abundance

# Subset dataset into treatment and control groups
basal <- subset_samples(mt_rela, treatment_group == "BL") # Filter data for fresh and necropsy
GCLAR <- subset_samples(mt_rela, treatment_group == "GC-LAR")
FLAR <- subset_samples(mt_rela, treatment_group == "F-LAR")
GCISST <- subset_samples(mt_rela, treatment_group == "GC-ISST")
FISST <- subset_samples(mt_rela, treatment_group == "F-ISST")

basal_fresh <- subset_samples(basal, sample_type == "Fresh")
basal_necropsy <- subset_samples(basal, sample_type == "Necropsy")
GCLAR_fresh <- subset_samples(GCLAR, sample_type == "Fresh")
GCLAR_necropsy <- subset_samples(GCLAR, sample_type == "Necropsy")
FLAR_fresh <- subset_samples(FLAR, sample_type == "Fresh")
FLAR_necropsy <- subset_samples(FLAR, sample_type == "Necropsy")
GCISST_necropsy <- subset_samples(GCISST, sample_type == "Necropsy")
FISST_necropsy <- subset_samples(FISST, sample_type == "Necropsy")

#week 0
GCLAR_week0 <- subset_samples(GCLAR, Time == "Launch plus 0")
FLAR_week0 <- subset_samples(FLAR, Time == "Launch plus 0")


#week 4.5 
GCLAR_week4_5 <- subset_samples(GCLAR, Time == "Launch plus 4.5")
FLAR_week4_5 <- subset_samples(FLAR, Time == "Launch plus 4.5")


#WEEK 9
GCLAR_week9 <- subset_samples(GCLAR, Time == "Launch plus 9")
GCLAR_week9_fresh <- subset_samples(GCLAR_week9, sample_type == "Fresh")
GCLAR_week9_necropsy <- subset_samples(GCLAR_week9, sample_type == "Necropsy")
FLAR_week9 <- subset_samples(FLAR, Time == "Launch plus 9")
FLAR_week9_fresh <- subset_samples(FLAR_week9, sample_type == "Fresh")
FLAR_week9_necropsy <- subset_samples(FLAR_week9, sample_type == "Necropsy")


# Set the prevalence threshold and abundance threshold
basal_ASVs <- core_members(basal, detection = 0.001, prevalence = 0.10)
GCLAR_ASVs <- core_members(GCLAR, detection = 0.001, prevalence = 0.10)
FLAR_ASVs <- core_members(FLAR, detection = 0.001, prevalence = 0.10)
GCISST_ASVs <- core_members(GCISST, detection = 0.001, prevalence = 0.10)
FISST_ASVs <- core_members(FISST, detection = 0.001, prevalence = 0.10)

basal_fresh_ASVs <- core_members(basal_fresh, detection = 0.001, prevalence = 0.10)
basal_necropsy_ASVs <- core_members(basal_necropsy, detection = 0.001, prevalence = 0.10)
GCLAR_fresh_ASVs <- core_members(GCLAR_fresh, detection = 0.001, prevalence = 0.10)
GCLAR_necropsy_ASVs <- core_members(GCLAR_necropsy, detection = 0.001, prevalence = 0.10)
FLAR_fresh_ASVs <- core_members(FLAR_fresh, detection = 0.001, prevalence = 0.10)
FLAR_necropsy_ASVs <- core_members(FLAR_necropsy, detection = 0.001, prevalence = 0.10)
GCISST_necropsy_ASVs <- core_members(GCISST_necropsy, detection = 0.001, prevalence = 0.10)
FISST_necropsy_ASVs <- core_members(FISST_necropsy, detection = 0.001, prevalence = 0.10)

#week 0
GCLAR_week0_ASVs <- core_members(GCLAR_week0, detection = 0.001, prevalence = 0.10)
FLAR_week0_ASVs <- core_members(FLAR_week0, detection = 0.001, prevalence = 0.10)

#week 4.5
GCLAR_week4_5_ASVs <- core_members(GCLAR_week4_5, detection = 0.001, prevalence = 0.10)
FLAR_week4_5_ASVs <- core_members(FLAR_week4_5, detection = 0.001, prevalence = 0.10)


#week 9
GCLAR_week9_ASVs <- core_members(GCLAR_week9, detection = 0.001, prevalence = 0.10)
GCLAR_week9_fresh_ASVs <- core_members(GCLAR_week9_fresh, detection = 0.001, prevalence = 0.10)
GCLAR_week9_necropsy_ASVs <- core_members(GCLAR_week9_necropsy, detection = 0.001, prevalence = 0.10)
FLAR_week9_ASVs <- core_members(FLAR_week9, detection = 0.001, prevalence = 0.10)
FLAR_week9_fresh_ASVs <- core_members(FLAR_week9_fresh, detection = 0.001, prevalence = 0.10)
FLAR_week9_necropsy_ASVs <- core_members(FLAR_week9_necropsy, detection = 0.001, prevalence = 0.10)

# The detection is set to 0.001 (0.1% relative abundance) in order to filter out rare things
# and only consider those that are somewhat abundant 
# The prevalence threshold is set to 0.1 in order to see differences between groups
# which means it only has to be present in 10% of samples in that group to be considered

# basal vs gclar no split
venn_plot_basal_GCLAR_no_split <- ggVennDiagram(x = list(Basal = basal_ASVs,
                                                         GCLAR = GCLAR_ASVs)) +
  scale_fill_gradient(low = "white", high = "blue", limits = c(0, 100)) + # Plot aesthetics
  scale_x_continuous(expand = expansion(mult = .1))

# Save the Venn diagram as a png
ggsave(filename = "core_microbiome_basal_GCLAR_no_split.png", plot = venn_plot_basal_GCLAR_no_split, + scale_x_discrete(expand = c(0, 1), width = 10, height = 10)

# GCLAR week 9 vs GCISST
venn_plot_GCLARweek9_GCISST_no_split <- ggVennDiagram(x = list(GCLAR = GCLAR_week9_ASVs,
                                                         GCISST = GCISST_ASVs), label_size = 100) +
  scale_fill_gradient(low = "white", high = "blue", limits = c(0, 100)) + # Plot aesthetics
  scale_x_continuous(expand = expansion(mult = .1))

# Save the Venn diagram as a png
ggsave(filename = "core_microbiome_GCISST_GCLAR_no_split.png", plot = venn_plot_GCLARweek9_GCISST_no_split, width = 11, height = 11)

#F-ISST x F-LAR 9 WEEKS NECROPSY
venn_plot_FISST_FLAR_no_split_week9 <- ggVennDiagram(x = list(FISST = FISST_ASVs,
                                                        FLAR = FLAR_week9_ASVs), label_size = 100) +
  scale_fill_gradient(low = "white", high = "blue", limits = c(0, 100)) + # Plot aesthetics
  scale_x_continuous(expand = expansion(mult = .1))

# Save the Venn diagram as a png
ggsave(filename = "core_microbiome_FISST_FLAR_no_split_week9.png", plot = venn_plot_FISST_FLAR_no_split_week9, width = 11, height = 11)

#GCLAR FLAR WEEK 0 X 0
venn_plot_GCLAR_FLAR_week0 <- ggVennDiagram(x = list(GCLAR_Week_0 = GCLAR_week0_ASVs,
                                                              FLAR_Week_0 = FLAR_week0_ASVs), label_size = 100) +
  scale_fill_gradient(low = "white", high = "blue", limits = c(0, 100)) + # Plot aesthetics
  scale_x_continuous(expand = expansion(mult = .1))

# Save the Venn diagram as a png
ggsave(filename = "core_microbiome_GCLAR_FLAR_week0.png", plot = venn_plot_GCLAR_FLAR_week0, width = 11, height = 11)

#GCLAR FLAR WEEK 4.5
venn_plot_GCLAR_FLAR_week4_5 <- ggVennDiagram(x = list(GCLAR_Week_4_5 = GCLAR_week4_5_ASVs,
                                                     FLAR_Week_4_5 = FLAR_week4_5_ASVs), label_size = 100) +
  scale_fill_gradient(low = "white", high = "blue", limits = c(0, 100)) + # Plot aesthetics
  scale_x_continuous(expand = expansion(mult = .1))

# Save the Venn diagram as a png
ggsave(filename = "core_microbiome_GCLAR_FLAR_week4_5.png", plot = venn_plot_GCLAR_FLAR_week4_5 , width = 11, height = 11)

#GCLAR FLAR WEEK 9 FRESH
venn_plot_GCLAR_FLAR_week9_fresh <- ggVennDiagram(x = list(A = GCLAR_week9_fresh_ASVs,
                                                       B = FLAR_week9_fresh_ASVs), label_size = 100) +
  scale_fill_gradient(low = "white", high = "blue", limits = c(0, 100)) + # Plot aesthetics
  scale_x_continuous(expand = expansion(mult = .1))

# Save the Venn diagram as a png
ggsave(filename = "core_microbiome_GCLAR_FLAR_week9_fresh.png", plot = venn_plot_GCLAR_FLAR_week9_fresh , width = 11, height = 11)

#GCLAR FLAR WEEK 9 necropsy
venn_plot_GCLAR_FLAR_week9_Necropsy <- ggVennDiagram(x = list(GCLAR_Week_9_N = GCLAR_week9_necropsy_ASVs,
                                                           FLAR_week9_N = FLAR_week9_necropsy_ASVs), label_size = 100) +
  scale_fill_gradient(low = "white", high = "blue", limits = c(0, 100)) + # Plot aesthetics
  scale_x_continuous(expand = expansion(mult = .1))

# Save the Venn diagram as a png
ggsave(filename = "core_microbiome_GCLAR_FLAR_week9_necropsy.png", plot = venn_plot_GCLAR_FLAR_week9_Necropsy, width = 11, height = 11)

#GCLAR WEEK 0, 4.5, 9 fresh
venn_plot_GCLAR_week0_9_fresh <- ggVennDiagram(x = list(GCLAR_Week_9_F = GCLAR_week9_fresh_ASVs,
                                                        GCLAR_Week_0 = GCLAR_week0_ASVs,
                                                        GCLAR_Week_4_5 = GCLAR_week4_5_ASVs), label_size = 100) +
  scale_fill_gradient(low = "white", high = "blue", limits = c(0, 100)) + # Plot aesthetics
  scale_x_continuous(expand = expansion(mult = 0.01))

# Save the Venn diagram as a png
ggsave(filename = "core_microbiome_GCLAR_week0_9_fresh.png", plot = venn_plot_GCLAR_week0_9_fresh, width = 11, height = 11)

#GCLAR WEEK 0, 4.5, 9 NECROPSY
venn_plot_GCLAR_week0_9_necropsy <- ggVennDiagram(x = list(GCLAR_Week_9_N = GCLAR_week9_necropsy_ASVs,
                                                        GCLAR_Week_0 = GCLAR_week0_ASVs,
                                                        GCLAR_Week_4_5 = GCLAR_week4_5_ASVs)) +
  scale_fill_gradient(low = "white", high = "blue", limits = c(0, 100)) + # Plot aesthetics
  scale_x_continuous(expand = expansion(mult = .1))

# Save the Venn diagram as a png
ggsave(filename = "core_microbiome_GCLAR_week0_9_necropsy.png", plot = venn_plot_GCLAR_week0_9_necropsy, width = 11, height = 11)

#FISST VS GISST 
venn_plot_FISST_GCISST <- ggVennDiagram(x = list(FISST = FISST_ASVs,
                                                GCISST = GCISST_ASVs), label_size = 10) +
  scale_fill_gradient(low = "white", high = "blue", limits = c(0, 100)) + # Plot aesthetics
  scale_x_continuous(expand = expansion(mult = .1))

# Save the Venn diagram as a png
ggsave(filename = "core_microbiome_FISST_GCISST.png", plot = venn_plot_FISST_GCISST, width = 11, height = 11)

# Rarefy samples
# rngseed sets a random number to sample. If you want to reproduce this exact analysis, you need
# to set rngseed the same number each time
# t transposes the table to use rarecurve function
# cex decreases font size
# rarecurve(t(as.data.frame(otu_table(micro_final))), cex=0.1) #this took long time and it is just a visualization, dont need to run
# micro_rare <- rarefy_even_depth(micro_final, rngseed = 1, sample.size = 64411) #chosen rarefaction depth

##### Saving #####
# save(micro_final, file="micro_final.RData") #rarefaction not always needed
# save(micro_rare, file="micro_rare.RData") #needed for diversity metrics
