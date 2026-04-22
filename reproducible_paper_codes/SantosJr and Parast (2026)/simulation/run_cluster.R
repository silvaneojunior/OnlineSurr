# ------------------------------------------------------------------------------
# This is a optional code, meant to run multiple simulation batches in parallel
# ------------------------------------------------------------------------------

library(parallel)
n.cores=detectCores()-2

file=local_path('run_single.R')

prepare_parallel=function(file,batch.idx,batch.size,cases){
  set.seed(batch.idx)
  assign('batch.idx',batch.idx,.GlobalEnv)
  assign('batch.size',batch.size,.GlobalEnv)
  assign('cases',cases,.GlobalEnv)
  source(file)
}


# Simulation 1
set.seed(11041992)
# To get the exact same results as the paper, it is necessary to use the exact same seeds, which depend on the number of cores used.
# For the paper, we used n.cores=20.
seeds=round(runif(n.cores)*1000000)

cases=paste0('Study 1 - Case ',1:3)
clust=makeCluster(n.cores)
clusterApply(clust,seeds,prepare_parallel,batch.size=2000/n.cores,file=file,cases=cases)
stopCluster(clust)

# Simulation 2
set.seed(13031998)
# To get the exact same results as the paper, it is necessary to use the exact same seeds, which depend on the number of cores used.
# For the paper, we used n.cores=20.
seeds=round(runif(n.cores)*1000000)
cases=paste0('Study 2 - Case ',2:5)
clust=makeCluster(n.cores)
clusterApply(clust,seeds,prepare_parallel,batch.size=3000/n.cores,file=file,cases=cases)
stopCluster(clust)
