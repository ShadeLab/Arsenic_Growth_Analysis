# Arsenic_Growth_Analysis
This folder contains all raw data, contextual data, individual analysis, and meta analysis for growth curves of arsenic resistant isolates cultivated from surface soil of an active vent in Centralia, PA in 2014. 

##File organization
### /Arsenic_Growth_Analysis/
* __Raw data file__: Date_AsType_MIC_IsolateAbbreviation.csv
* __Plate organization__ (created manually): Date_platemap.csv
* __Modified plate organization for grofit__ (created manually): Date_platemap_grofit.csv
* __Annotated raw data__ (combined raw data and plate organization): Date_MIC_annotated.csv
* __Results__ (raw results from grofit): Date_results.csv
* __Header for data__: Time,wells.csv
* __Compiled grofit data__: orig_grofit.csv
* __Compiled grofit data with binary column__ (Quality) describing whether to use model (0) or spline (1): orig_grofit_model.csv

### /Arsenic_Growth_Analysis/R_scripts/Grofit/
* __Grofit run for individual experiment__: Date_Data_Analysis.R
* __Grofit analysis for all isolates__: Meta_analysis_FigureS2.R
* __Grofit analysis for all isolates by genus__: Meta_analysis_Figure4.R

### /Arsenic_Growth_Analysis/R_scripts/Growth_curves
* __Growth curves for individual experiment__: Date_Growth_analysis.R
* __Final growth curves for all isolates__: MIC_curves_FigureS1.R

##Workflow
__1. Individual experiment (plate) analysis__
    * Growth curves were made for each individual experiment (plate)
      * Script: /Arsenic_Growth_Analysis/R_scripts/Growth_curves/Date_Growth_analysis.R
      * Input: Date_AsType_MIC_IsolateAbbreviation.csv, Date_platemap.csv, Time,wells.csv
      * Output: Date_MIC_annotated.csv
    * Each experiment (plate) was analyzed using Grofit
      * Script: /Arsenic_Growth_Analysis/R_scripts/Grofit/Date_Data_Analysis.R
      * Input: Date_AsType_MIC_IsolateAbbreviation.csv, Date_platemap.csv
      * Optional input: Date_platemap_grofit.csv (used when certain wells were removed from analysis (less than 5 pts >0) 
      * Output: Date_results.csv
      
__2. Meta analysis__
    * All final growth curves from experiment were compiled for Figure S1
      * Script: /Arsenic_Growth_Analysis/R_scripts/Growth_curves/MIC_curves_FigureS1.R
      * Input: each individual Date_MIC_annotated.csv
      * Output: Figure S1
    * All grofit data was compiled to make individual growth parameter curves for each isolate (Figure S2)
      * Script: /Arsenic_Growth_Analysis/R_scripts/Grofit/Meta_analysis_FigureS2.R
      * Input: each individual Date_results.csv, orig_grofit_model.csv
      * Output: orig_grofit.csv, Figure S2
    * Grofit data for max growth rate and lag time were combined by genus for all isolates (Figure 4)
      * Script: /Arsenic_Growth_Analysis/R_scripts/Grofit/Meta_analysis_Figure4.R
      * Input: Date_results.csv, orig_grofit_model.csv
      * Output: orig_grofit.csv, Figure 4
