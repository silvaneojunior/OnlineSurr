# ------------------------------------------------------------------------------
# Auxiliary code to define the simulation parameters and import dependencies.
# This script is also used to import the other dependencies.
# ------------------------------------------------------------------------------

library(tidyverse)
library(cubature)
library(stringr)
library(tictoc)
library(kDGLM)
library(Rfast)
library(readr)
library(clubSandwich)
library(extrafont)
library(lme4)
library(geepack)
library(latex2exp)
library(knitr)

path='reproducible_paper_codes/SantosJr and Parast (2026)/simulation/'

local_path=function(file){
  paste0(path,'functions.R')
}

ginv=kDGLM:::ginv
rmvnorm=kDGLM:::rmvnorm
source(local_path('functions.R'))

# Configuration of the hypothetical study
freq=12/3 # Number of measurements per year
V=0.05/freq # Variance for the observational noise
W=0.2*V # Variance for the evolution noise
phi1=0.95 # AR coefficient for the outcome evolution
phi2=0.95 # AR coefficient for the surrogate evolution
S.effect=1 # Scale of the surrogate effect

N.bootstrap=2000

# Coefficient for the CAR prior for the treatment effect
# This is optional, but may be useful
# Values closer to 1 imply in a treatment effect that changes smoothly over time
# Values closer to 0 imply in a treatment effect that changes drastically over time
# In a non-Bayesian perspective, you can think of this as a penalty for time heterogeneity for the treatment effect
# In practice, the prior is VERY vague (the penalty is very small), so this has barely any influence in the final result.
rho=0.9


T.list=c(10,5,3,2,1) # Set of study lengths to be tested
N.list=c(300,250,200,150,100,50) # Set of sample sizes to be tested

# Set of total treatment effects to be tested
# Different settings may have different values
# The values were chosen so that the effect of the sample size would be visible in the graphs
total.effect=list('Study 1 - Case 1'=5,
              'Study 1 - Case 2'=5,
              'Study 1 - Case 3'=5,
              'Study 2 - Case 2'=1,
              'Study 2 - Case 3'=1,
              'Study 2 - Case 4'=1,
              'Study 2 - Case 5'=1)
# Set of true PTEs to be tested
# For study 2, we use a single starting PTE, as that is the time heterogenous setting, so the PTE changes over time.
pte.list=list('Study 1 - Case 1'=c(0.25,0.5,0.75,0.9,1),
              'Study 1 - Case 2'=c(0.25,0.5,0.75,0.9,1),
              'Study 1 - Case 3'=c(0.25,0.5,0.75,0.9,1),
              'Study 2 - Case 2'=c(0.9),
              'Study 2 - Case 3'=c(0.9),
              'Study 2 - Case 4'=c(0.9),
              'Study 2 - Case 5'=c(0.9))


# ------------------------------------------------------------------------------
#### Loading fonts ####
# Run only once
# font_import()

loadfonts(device = "win")
font_size=12
family_font='serif'
par(family=family_font)

theme_custom=theme_bw()+theme(text=element_text(size=font_size,family=family_font),
                              legend.text = element_text(size=font_size,family=family_font),
                              axis.title = element_text(size=font_size,family=family_font),
                              legend.position = 'bottom',
                              panel.grid = element_blank())
# ------------------------------------------------------------------------------
