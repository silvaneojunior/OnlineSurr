# ------------------------------------------------------------------------------
# This is a optional code, meant to run multiple simulation batches in parallel
# ------------------------------------------------------------------------------

library(parallel)
n.cores=detectCores()-2

file=local_path('run_single.R')

prepare_parallel=function(file,batch.idx,batch.size){
  set.seed(batch.idx)
  assign('batch.idx',batch.idx,.GlobalEnv)
  assign('batch.size',batch.size,.GlobalEnv)
  source(file)
}


set.seed(13031998)
# For the paper, we used n.cores=20
seeds=round(runif(n.cores)*1000000)

clust=makeCluster(n.cores)
clusterApply(clust,seeds,prepare_parallel,batch.size=50,file=file)
stopCluster(clust)
