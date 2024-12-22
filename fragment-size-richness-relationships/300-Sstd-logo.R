# leave-one-group-out cross-validation; executed on scientific computing cluster

library(tidyverse)
library(brms)
library(future)

load('/data/idiv_chase/sablowes/evid-synth/results/sstd_frag-4338431.Rdata')
load('/data/idiv_chase/sablowes/evid-synth/results/sstd_frag_fs-4427402.Rdata')

cpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
future::plan(future::multicore(workers = cpus))

# model 2.1
cv10g_sstd <- kfold(x = Sstd_lognorm_fragSize,
                            group = 'dataset_label')
# model 2.2
cv10g_sstd_sigma <- kfold(x = Sstd_lognorm_fragSize_sigma,
                                  group = 'dataset_label')
# model 2.3
cv10g_sstd_sigma_fs <- kfold(x = Sstd_lognorm_fragSize_sigma_fs,
                             group = 'dataset_label')
# model 2.4
cv10g_sstd_sigma_cor <- kfold(x = Sstd_lognorm_fragSize_sigma_cor,
                                      group = 'dataset_label')
# model 2.5
cv10g_sstd_sigma_fs_cor <- kfold(x = Sstd_lognorm_fragSize_sigma_fs_cor, 
                                         group = 'dataset_label')
# save results
save(cv10g_sstd,
     cv10g_sstd_sigma,
     cv10g_sstd_sigma_fs,
     cv10g_sstd_sigma_cor,
     cv10g_sstd_sigma_fs_cor,
     file = Sys.getenv('OFILE'))