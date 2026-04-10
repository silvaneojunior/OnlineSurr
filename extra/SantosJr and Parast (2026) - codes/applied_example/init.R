library(tidyverse)
library(Rfast)
library(kDGLM)
library(knitr)
library(qs)


path='extra/SantosJr and Parast (2026) - codes/applied_example/'

local_path=function(file){
  paste0(path,'functions.R')
}

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

# Sorry, but your data is in another castle :(
raw.data=read.csv(local_path('data_long_stream.csv'))
dim(raw.data)
head(raw.data)

# Changes to the data format:
#  - Making treatment a boolean variable
#  - Renaming some variable
#  - Swithing to a longer format: each combination of patient-time has its own line
#  - Filling the missing lines with NA. Note that the measurement of the surrogate and of the outcome have different frequencies.
data=raw.data %>%
  mutate(TREAT=(GROUP=='EXPERIMENTAL')) %>%
  select(MASK_PAT,
         TREAT,
         SEX=FEMALE,SMOKER,BATCH,
         3:39) %>%
  pivot_longer(-c(1:5)) %>%
  mutate(name=str_split_i(name,'V',2)) %>%
  rename(HBAQV=value,Time=name) %>%
  full_join(
    raw.data %>%
      select(MASK_PAT,
             40:49) %>%
      pivot_longer(-1) %>%
      mutate(name=str_split_i(name,'R',2)) %>%
      rename(LAER=value,Time=name),
    by=c('MASK_PAT','Time')) %>%
  arrange(MASK_PAT,Time) %>%
  mutate(Time=as.numeric(Time)+1)

# Checks for duplicated patients. If 1, then there is no duplication.
raw.data$MASK_PAT %>% table %>% max

# Changing the center and scale of the data.
data$LAER=-(data$LAER %>% {(.-mean(.,na.rm=TRUE))/sd(.,na.rm=TRUE)})
data$HBAQV=-(data$HBAQV %>% {(.-mean(.,na.rm=TRUE))/sd(.,na.rm=TRUE)})

# Taking a look at the result
ggplot(
  data %>%
    group_by(Time,TREAT) %>%
    summarize(LAER.icl=quantile(LAER,0.025,na.rm=TRUE),
              LAER.icu=quantile(LAER,0.975,na.rm=TRUE),
              LAER=mean(LAER,na.rm=TRUE),
              count=sum(!is.na(LAER))))+
  geom_point(aes(x=Time %>% as.numeric,y=LAER,color=TREAT))+
  geom_line(aes(x=Time %>% as.numeric,y=LAER,color=TREAT))+
  geom_errorbar(aes(x=Time %>% as.numeric,
                    ymin=LAER+(LAER.icl-LAER)/sqrt(count),
                    ymax=LAER+(LAER.icu-LAER)/sqrt(count),color=TREAT))+
  theme_custom

ggplot(
  data %>%
    group_by(Time,TREAT) %>%
    summarize(HBAQV.icl=quantile(HBAQV,0.025,na.rm=TRUE),
              HBAQV.icu=quantile(HBAQV,0.975,na.rm=TRUE),
              HBAQV=mean(HBAQV,na.rm=TRUE),
              count=sum(!is.na(HBAQV))) %>% ungroup
)+
  geom_point(aes(x=Time %>% as.numeric,y=HBAQV,color=TREAT))+
  geom_line(aes(x=Time %>% as.numeric,y=HBAQV,color=TREAT))+
  theme_custom



# Defining important quantities
N=data$MASK_PAT %>% unique %>% length
T=(data %>% filter(!is.na(LAER)))$Time %>% unique %>% length # Number of outcome measurements
T.full=data$Time %>% unique %>% length # Number of measurement (i.e., including times where the outcome was NOT measured)
S=data$HBAQV %>% matrix(T.full,N)
Y=data$LAER %>% matrix(T.full,N)

# Making sure that when the surrogate stopped being measured the outcome was dropped
# Otherwise, we would have to consider a model with missing covariates, which would be much more complicated
# In practice, there are no cases like this in the dataset, this part is just to be safe.
S.index=(data %>% filter(!is.na(LAER)))$Time %>% unique
for(i in S.index){
  Y[i:T.full,is.na(S[i,])]=NA
}

# Number of observed outcome measurements per patient.
# Patients with no outcome measurements would be excluded;
# in practice, none were dropped. This check is included as a safeguard.
time.length = T.full - colSums(is.na(Y))

# Keep only the covariate values at time points included in the outcome model.
# This refers to the study design, not to patient-specific missingness. In other
# words, some time points were not scheduled for outcome measurement for any
# patient, although the surrogate was still measured at those times.
#
# The outcome model is indexed only by the time points at which the outcome was
# scheduled to be measured. Time points that were excluded by design from the
# outcome process are therefore not part of the model index. However, surrogate
# measurements collected at those times are still used as covariates.
#
# This is distinct from ordinary missing data for individual patients: a time
# point is removed only when the outcome was not measured for any patient, not
# when the outcome is missing for a subset of patients.
#
# Because the baseline covariates are time-invariant, it is sufficient to retain
# their values only at the time points used by the outcome model.
sex.base=matrix(data$SEX-1,T.full,N)[S.index,time.length>0,drop=FALSE]
treat.base=matrix(data$TREAT,T.full,N)[S.index,time.length>0,drop=FALSE]
smoker.base=matrix(data$SMOKER,T.full,N)[S.index,time.length>0,drop=FALSE]

# Removing patients that were not measured (no one was removed)
S.redux=S[,time.length>0,drop=FALSE]
Y=Y[S.index,time.length>0,drop=FALSE]


# Making sure that when the surrogate stopped being measured the outcome was dropped
for(i in 1:T){
  Y[i:T,is.na(Y[i,])]=NA
}

# Redefining the sample size in case someone was dropped (no one was removed)
N=dim(Y)[2]
