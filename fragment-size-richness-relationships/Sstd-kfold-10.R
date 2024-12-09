# stratified kfold cross-validation; executed on scientific computing cluster

library(tidyverse)
library(brms)
library(future)

load('/data/idiv_chase/sablowes/evid-synth/results/sstd_frag-4338431.Rdata')
load('/data/idiv_chase/sablowes/evid-synth/results/sstd_frag_fs-4427402.Rdata')

cpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
future::plan(future::multicore(workers = cpus))

Sstd_fragSize_kf10 <- kfold(x = Sstd_lognorm_fragSize,
                            folds = 'stratified', 
                            group = 'dataset_label',
                            k = 10, chains = 4)

Sstd_fragSize_sigma_kf10 <- kfold(x = Sstd_lognorm_fragSize_sigma,
                            folds = 'stratified', 
                            group = 'dataset_label',
                            k = 10, chains = 4)

Sstd_fragSize_sigma_cor_kf10 <- kfold(x = Sstd_lognorm_fragSize_sigma_cor,
                            folds = 'stratified', 
                            group = 'dataset_label',
                            k = 10, chains = 4)

Sstd_fragSize_sigma_fs_kf10 <- kfold(x = Sstd_lognorm_fragSize_sigma_fs,
                            folds = 'stratified', 
                            group = 'dataset_label',
                            k = 10, chains = 4)

Sstd_fragSize_sigma_fs_cor_kf10 <- kfold(x = Sstd_lognorm_fragSize_sigma_fs_cor,
                            folds = 'stratified', 
                            group = 'dataset_label',
                            k = 10, chains = 4)

save(Sstd_fragSize_kf10,
     Sstd_fragSize_sigma_kf10,
     Sstd_fragSize_sigma_cor_kf10,
     Sstd_fragSize_sigma_fs_kf10,
     Sstd_fragSize_sigma_fs_cor_kf10,
     file = Sys.getenv('OFILE'))