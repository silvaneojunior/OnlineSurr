# Runs a batch of replications of the simulation.

source('reproducible_paper_codes/SantosJr and Parast (2026)/simulation/init.R')

file=file(local_path(paste0('data/',timestamp(),'.csv')), "w") # Create the file to store the results
run_set(((batch.idx-1)*batch.size):(batch.idx*batch.size),cases) # Run the simulations
close(file) # Close the file
