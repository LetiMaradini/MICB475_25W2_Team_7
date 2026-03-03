## Meeting agenda ##
- Discuss and go over proposal:
  * Add more background to introduction
  * Further expanding and clarifying the comparison groups (experimental aims)
      * Decide groups used for comparision (planning to use model to predict all of the complex scenarios, or aiming to predict gravity versus no gravity?)
      * Further discussion of individual factors (space versus gravity, time, gravity stimulation, etc) is needed in your experimental aims
  * Include comparison groups used to proposed approach
      * Determine which statistical tests will be used
 - Make sure everyone is ready to go/understands next analysis (diversity metrics, core microbiome, indicator taxa, deseq) due on March 8th so there is enough time for machine learning (due last week of March)


## Meeting notes ##
* Groups to be compared (add to the aims):
  * Difference between microgravity and control: F-ISST x GC-ISST
  * Difference between different time-points within the same cohort (0 x 4.5 x 9 weeks after euthanasia): GC-LAR and F-LAR
  * Difference between the same time-points in the cohorts (0 x 0, 4.5 x 4.5, and 9 weeks after euthanasia x 9 weeks after euthanasia): GC-LAR x F-LAR
  * Difference between return (time for recovery) and flight (both at 9 weeks after euthanasia): F-ISST x F-LAR
  * Compare controls (at week 9 after euthanasia): GC-LAR x GC-ISST
  * Compare control and basal (week 0 after euthanasia): GC-LAR x Basal
 
* For DESeq, blast the ones that do not have identification
* Remove the pre-euthasia data from the phyloseq object
