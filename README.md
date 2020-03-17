## Github Repository for
# Taxonomically-linked growth phenotypes during arsenic stress among arsenic resistant bacteria isolated from soils overlying the Centralia coal seam fire
## by Taylor K. Dunivin, Justine Miller, Ashley Shade 

<i>This work is published.</i>

### Data
All sequences >200 bp were submitted to NCBI, and sequences can be accessed from GenBank with the following accession numbers: 16S rRNA KX825887- KX825911, arsC KY405022- KY405029, ACR3(2) KY405030- KY405032, and arsB KY405033- KY405040.

### To cite this work or code
Dunivin TK, Miller J, Shade A (2018) Taxonomically-linked growth phenotypes during arsenic stress among arsenic resistant bacteria isolated from soils overlying the Centralia coal seam fire. PLoS ONE 13(1): e0191893. https://doi.org/10.1371/journal.pone.0191893

### Abstract
Arsenic (As), a toxic element, has impacted life since early Earth. Thus, microorganisms have evolved many As resistance and tolerance mechanisms to improve their survival outcomes given As exposure. We isolated As resistant bacteria from Centralia, PA, the site of an underground coal seam fire that has been burning since 1962. From a 57.4°C soil collected from a vent above the fire, we isolated 25 unique aerobic As resistant bacterial strains spanning seven genera. We examined their diversity, resistance gene content, transformation abilities, inhibitory concentrations, and growth phenotypes. Although As concentrations were low at the time of soil collection (2.58 ppm), isolates had high minimum inhibitory concentrations (MICs) of arsenate and arsenite (>300 mM and 20 mM respectively), and most isolates were capable of arsenate reduction. We screened isolates (PCR and sequencing) using 12 published primer sets for six As resistance genes (AsRGs). Genes encoding arsenate reductase (arsC) and arsenite efflux pumps (arsB, ACR3(2)) were present, and phylogenetic incongruence between 16S rRNA genes and AsRGs provided evidence for horizontal gene transfer. A detailed investigation of differences in isolate growth phenotypes across As concentrations (lag time to exponential growth, maximum growth rate, and maximum OD590) showed a relationship with taxonomy, providing information that could help to predict an isolate’s performance given As exposure in situ. Our results suggest that microbiological management and remediation of environmental As could be informed by taxonomically-linked As tolerance, potential for resistance gene transferability, and the rare biosphere.

### Funding
The project was supported by start-up funds from Michigan State University to AS. This work was supported in part by Michigan State University through computational resources provided by the [Institute for Cyber-Enabled Research](https://icer.msu.edu/). 

### More info
[ShadeLab](http://ashley17061.wixsite.com/shadelab/home)


# Arsenic_Growth_Analysis
This folder contains all raw data, contextual data, individual analysis, and meta analysis for growth curves of arsenic resistant isolates cultivated from surface soil of an active vent in Centralia, PA in 2014. 

## File organization
### /Arsenic_Growth_Analysis/
* __Raw data file__: Date_AsType_MIC_IsolateAbbreviation.csv
* __Plate organization__ (created manually): Date_platemap.csv
* __Modified plate organization for grofit__ (created manually): Date_platemap_grofit.csv
* __Annotated raw data__ (combined raw data and plate organization): Date_MIC_annotated.csv
* __Results__ (raw results from grofit): Date_results.csv
* __Header for data__: Time,wells.csv
* __Compiled grofit data__: orig_grofit.csv
* __Compiled grofit data with binary column__ (Quality) describing whether to use model (0) or spline (1): orig_grofit_model.csv
* __Isolate and genus key__: isolate_genus.csv

### /Arsenic_Growth_Analysis/R_scripts/Grofit/
* __Grofit run for individual experiment__: Date_Data_Analysis.R
* __Grofit analysis for all isolates__: Meta_analysis_FigureS2.R
* __Grofit analysis for all isolates by genus__: Meta_analysis_Figure4.R

### /Arsenic_Growth_Analysis/R_scripts/Growth_curves
* __Growth curves for individual experiment__: Date_Growth_analysis.R
* __Final growth curves for all isolates__: MIC_curves_FigureS1.R

##Workflow
1. Individual experiment (plate) analysis
    * Growth curves were made for each individual experiment (plate)
      * __Script__: /Arsenic_Growth_Analysis/R_scripts/Growth_curves/Date_Growth_analysis.R
      * __Input__: Date_AsType_MIC_IsolateAbbreviation.csv, Date_platemap.csv, Time,wells.csv
      * __Output__: Date_MIC_annotated.csv
    * Each experiment (plate) was analyzed using Grofit
      * __Script__: /Arsenic_Growth_Analysis/R_scripts/Grofit/Date_Data_Analysis.R
      * __Input__: Date_AsType_MIC_IsolateAbbreviation.csv, Date_platemap.csv
      * __Optional input__: Date_platemap_grofit.csv (used when certain wells were removed from analysis (less than 5 pts >0)
      * __Output__: Date_results.csv
      
2. Meta analysis
    * All final growth curves from experiment were compiled for Figure S1
      * __Script__: /Arsenic_Growth_Analysis/R_scripts/Growth_curves/MIC_curves_FigureS1.R
      * __Input__: each individual Date_MIC_annotated.csv
      * __Output__: Figure S1
    * All grofit data was compiled to make individual growth parameter curves for each isolate (Figure S2)
      * __Script__: /Arsenic_Growth_Analysis/R_scripts/Grofit/Meta_analysis_FigureS2.R
      * __Input__: each individual Date_results.csv, orig_grofit_model.csv
      * __Output__: orig_grofit.csv, Figure S2
    * Grofit data for max growth rate and lag time were combined by genus for all isolates (Figure 4)
      * __Script__: /Arsenic_Growth_Analysis/R_scripts/Grofit/Meta_analysis_Figure4.R
      * __Input__: Date_results.csv, isolate_genus.csv, orig_grofit_model.csv
      * __Output__: orig_grofit.csv, Figure 4
