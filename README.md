# Location-scale models and cross validation to advance quantitative evidence synthesis

This repository contains data and code to accompany 'Known unknowns and 
model selection in ecological evidence synthesis'
https://doi.org/10.1101/2024.12.18.629303

There are two folders - one for each case study presented in the ms.  
Briefly, each case study represents an example of two common forms  
of quantitative evidence synthesis in ecology: meta-analysis and an analysis of  
a primary data compilation. Each case study is described in detail in the accompanying paper. 

Use the init-dir.R file to set up the working directory, load packages and   
source custom functions for calculating and plotting the results.   

Each case study folder contains multiple subfolders:   
- one with the data needed to reproduce results,  
- one with the figures, and 
- one with the model fits and results of cross-validation.  
The main folder for each case study has files to: 
- fit models, 
- do simulation-based calibration, 
- do cross-validation (cross validation for the primary data case study  
was done a scientific computing cluster), 
- wrangle results and plot figures.

## Fragment size richness relationships

This is the primary data case study. It extends the work of  
Chase et al. 2020 *Nature* to examine heteroscedasticity.

## Native exotic richness relationships

This is the meta-analysis case study. It extends the work of  
Peng et al. 2019 *Ecology* to quantify heteroscedasticity.