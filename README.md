<img 
  src="documentation_materials/header_logos.png" 
  alt="University of Cambridge and Wellcome - MRC Cambridge Stem Cell Institute logos" 
  title="University of Cambridge and Wellcome - MRC Cambridge Stem Cell Institute logos" 
/>

# A single-cell comparison of adult and foetal human epicardium defines the age-associated changes in epicardial activity
#### Vincent R. K night-Schrijver

# Introduction
This report documents the full analysis of scRNA-seq data underlying the project entitled 
"*A single-cell comparison of adult and foetal human epicardium defines the age-associated changes in epicardial activity*".
. This project was a combination of other projects worked on simultaneously within the Sinha Lab initialised by a single
 Foetal heart sample, following by the slow acquisition of futher foetal samples andn finally catalysed by 
 the sudden appearance of an adult dataset with identifiable scRNA-seq epicardium. This led to the natural questions:

1. How does foetal epicardium differ from adult epicardium in normal conditions and otherwise healthy ageing?
2. Is there a signal here responsible for the perceived regenerative potential of hearts in more embryonic and neonatal stages of mammals?
3. How does the adult epicardium compare with our in vitro hESC-derived epicardium?

There are multiple sub-questions as well but these were the initial drivers of the project. As to the structure of this report, I will try to keep 
the reports separated into differenti stages of the analysis. We will begin with loading in some of the key packages and occasionally load an 
extra package where required within the specific stage of the analysis.

# Analysis workflow:
The overall structure of the analysis is as follows:

1. Read both adult and foetal data
2. Process data to similar QC standards
3. Integrate datasets for a single object
4. Clustering and cell-type annotation
5. Sub-cluster epicardium
6. Parallel stage-separate DEGs analysis for epicardial markers
7. Secreted epicardial factors and CellPhoneDB
8. Comparison with hESC-EPI and classification model



