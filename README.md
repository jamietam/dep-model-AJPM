# dep-model-AJPM
**Lifetime Prevalence of Major Depressive Episodes and Adjustment of Recall Error Through Simulation Modeling**
Code written by Jamie Tam

This repository contains code to run a simulation model of major depressive episodes that adjusts for recall error in estimates of lifetime major depressive episodes (MDEs). The model is run in R and calibrated to data from the National Survey on Drug Use and Health (NSDUH) 2005-2017. 

**Download, clean, and harmonize NSDUH data**
1. Download the full combined NSDUH dataset for 2002-2017 <a href="https://www.datafiles.samhsa.gov/study-dataset/nsduh-2002-2017-ds0001-nsduh-2002-2017-ds0001-nid18471"> here.</a>
2. Load the NSDUH data in R
3. Run `data_clean.R`. This produces the `depprevs_2005-2017.Rda` file.

**Perform model calibration and parameter estimation**
1. Specify the parameters to be estimated by adjusting the code in `main.R` (see code comments). This can also be done manually in `parameters.xlsx`.
2. Run `calibration.R` to generate parameter estimates that minimize the sum of squared distances between model and NSDUH data.

**Run the model and generate figures:**
1. Run `main.R` with the calibrated parameter estimates. 
2. Run `manuscript_figures.R`.
