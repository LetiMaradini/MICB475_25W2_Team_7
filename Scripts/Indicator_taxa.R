library(phyloseq)
library(ape) # importing trees
library(tidyverse)
library(vegan)

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
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species"),
           fill = "right") %>%
  as.matrix() # Saving as a matrix
# Save everything except feature IDs 
tax_mat <- tax_mat[,-1]
# Make sampleids the rownames
rownames(tax_mat) <- tax$`Feature ID`
# Keep only taxa present in both OTU and taxonomy
common_taxa <- intersect(rownames(otu_mat), tax$`Feature ID`)
otu_mat <- otu_mat[common_taxa, ]
tax_mat <- tax_mat[common_taxa, ]
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

# Rarefy samples
# rngseed sets a random number to sample. If you want to reproduce this exact analysis, you need
# to set rngseed the same number each time
# t transposes the table to use rarecurve function
# cex decreases font size
rarecurve(t(as.data.frame(otu_table(micro_final))), cex=0.1) #this took long time and it is just a visualization, dont need to run
micro_rare <- rarefy_even_depth(micro_final, rngseed = 1, sample.size = 64411) #chosen rarefaction depth

##### Saving #####
save(micro_final, file="micro_final.RData") #rarefaction not always needed
save(micro_rare, file="micro_rare.RData") #needed for diversity metrics


library(indicspecies)

# Extract OTU table
otu_mat <- as.data.frame(otu_table(micro_rare))

# Ensure taxa are columns
otu_mat <- t(otu_mat)

# Get grouping variable from metadata
group <- sample_data(micro_rare)$treatment_group

# Run indicator species analysis
indicator_multipatt <- multipatt(otu_mat, group, duleg = TRUE)

# View results
summary(indicator_multipatt)

# Save output
indicator_output <- capture.output(summary(indicator_multipatt, indvalcomp = TRUE))
writeLines(indicator_output, "indicator_values.txt")