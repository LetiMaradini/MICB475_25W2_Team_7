# February 3rd pre-meeting notes

## NASA datasets that work

- Space environmental factor impacts upon murine colon microbiota and mucosal homeostasis (OSD-72)
  https://osdr.nasa.gov/bio/repo/data/studies/OSD-72
- Temporal dynamics of the gut microbiota in people sharing a confined environment, a 520-day ground-based space simulation, MARS500 (OSD-191)
  https://osdr.nasa.gov/bio/repo/data/studies/OSD-191
- Reproducible changes in gut microbiome reveal a shift in microbial and host metabolism during spaceflight (OSD-212)
  https://osdr.nasa.gov/bio/repo/data/studies/OSD-212
- Specific Host Metabolite and Gut Microbiome Alterations Are Associated with Bone-loss During Spaceflight (OSD-417)
  https://osdr.nasa.gov/bio/repo/data/studies/OSD-417
- Salivary microbiome sequencing of astronauts (OSD-280)
  https://osdr.nasa.gov/bio/repo/data/studies/OSD-280



## Lila novel research idea

Research question: How does the ISS’ environmental microbiome affect the observed dysbiotic shifts in human oral microbiome during spaceflight?

Datasets: use OSD-280 (astronaut salivary microbiome) + OSD-72 (ISS environmental surface microbes)

Analysis: 
- quantify % of astronaut salivary microbiome that originates from ISS surfaces
- identify shared ASVs, see if certain ISS microbes are found in astronaut’s oral microbiome (pre-flight vs during vs post-flight)
- compare the diversity of human-associated microbes in the ISS environment surfaces to the expected decrease in diversity in astronaut oral microbiome (i.e do ISS surfaces containing more human-associated microbes stabilize the astronaut’s salivary microbiome?)

## Meeting suggestions

use ISS server data to look at abiotic microbiome vs astronaut oral microbiome (OSD-280)
- combine datasets (would need to normalize)
=> would be top heavy

or

use OSD-417 looking at microgravity affect on rodent gut microbiome
- functional analysis
- build machine learning model
=> would be end heavy


## Selected topic

OSD-417 dataset: Specific Host Metabolite and Gut Microbiome Alterations Are Associated with Bone-loss During Spaceflight

Novel question: looking at microgravity affect on rodent gut microbiome

Workflow
- process through QIIME
- compare diversity metrics
- core microbiome (how many shared vs unique microbes)
- annotated taxa (presence + abundance microbes in samples, and association with parameters)
- DESeq (abundance difference)
  
if not significant results
- functional analysis
  
if significant results
- machine learning model (prediction)


notes:
- Data should be processed to table.qza in Qiime2 for the proposal
- Always check with Ritu before trimming data etc
- Upload files after each step in Qiime2 to GitHub for backup

## Anson's Notes

effect of microgravity on gut microbiome in space document
- label is bone loss (but isn’t bone loss)

combining 280 and osd 272 
surface microbiome and oral microbiome on the ISS and see whether iss acts as buffer or amplifier to oral microobiome (both osds) compare and see what percentage (beta diversity metrics) to see how distant and close the microbiomes are and distance between pre and post flights. 

- idea is nice but can’t do so because you can’t compare mice and human (one on server called ISS is human (the one where they looked at different surfaces) 
- osd 191
- compare oral microbiome vs abiotic surfaces 
- the one on the server is very small so yeah. you have to e careful
- control (group and take fecal sample in week 45) 
- space and euthanize mice but both still ahve ground controls)
- 1 (diversity metrics)
- 2 (core microbiome
- 3 indicator taxa
- 4 deseq
- if data not significant than functional analysis (elevates project a little) —> piecrust 2 (kind of like kegg pathways)
- if data is signficiant then machine learning (random forced) there is code on the modules so 
- diversity metrics will give us a sense of if there is anu difference on a community level
- core microbiome see if there is any venn diagram any difference
- indicator taxa anything indicative conditions 
- could bin all ground samples so see which groups available 
- deseq —> to see if there is a expression difference
- if taxanomically not different than functionally it might be different 
- will have to have this data processed in qiime before the proposal

