# Modelling heterogeneity for quantitative evidence synthesis in ecology

Repo for models exploring heterogeneity in quantitative evidence synthesis.

There are two folders - one for each case study. Use the init-dir.R file
to set up the working directory, load packages and source a custom function
for plotting the results. Each folder contains multiple subfolders: one with
the data needed to reproduce results, one with the figures, and one with the 
model fits and results of cross-validation. The main folder for each case study
has the scripts to fit models, multiple files for simulation-based calibration 
for the series of models fit in each case study, do cross-valiation, 
and to wrangle results and create results figures.

Note that some scripts were written to be executed on a scientific computing 
cluster.

Each case study is described in detail in the accompanying paper. 

## Fragment size richness relationships

Extends the work of Chase et al. 2020 *Nature* to examine heteroscedasticity.

## Native exotic richness relationships

Extends the work of Peng et al. 2019 *Ecology* to quantify heterogeneity.