# Runs a batch of replications of the simulation.

source('extra/SantosJr and Parast (2026) - codes/simulation/init.R')

file=file(local_path(paste0('data/',timestamp(),'.csv')), "w") # Create the file to store the results
run_set(((batch.idx-1)*batch.size):(batch.idx*batch.size)) # Run the simulations
close(file) # Close the file
