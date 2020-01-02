# dep-model-AJPM
<bold> Lifetime Prevalence of Major Depressive Episodes and Adjustment of Recall Error Through Simulation Modeling</bold>

This repository contains code to run a simulation model of major depressive episodes that adjusts for recall error in estimates of lifetime MD episodes. The model is run in R and calibrated to data from the National Survey on Drug Use and Health (NSDUH) 2005-2017. 

<u>Download, clean, and harmonize NSDUH data:</u>
<bold>Step 1)</bold> Download the full combined NSDUH dataset for 2002-2017 <a href="https://www.datafiles.samhsa.gov/study-dataset/nsduh-2002-2017-ds0001-nsduh-2002-2017-ds0001-nid18471"> here.</a>
<bold>Step 2)</bold> Load the NSDUH data in R
<bold>Step 3)</bold> Run data_clean.R
This produces the depprevs_2005-2017.Rda file.

<u>Perform model calibration and parameter estimation:</u>
<bold>Step 1)</bold> Specify the parameters to be estimated by adjusting the code in 'main.R' (see code comments). This can also be done manually in 'parameters.xlsx'.
<bold>Step 2)</bold> Run 'calibration.R' to generate parameter estimates that minimize the sum of squared distances between model and NSDUH data.

<u>Run the model and generate figures:</u>
<bold>Step 1)</bold> Run 'main.R' with the calibrated parameter estimates. 
<bold>Step 2)</bold> Run 'manuscript_figures.R'.
