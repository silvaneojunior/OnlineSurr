source('reproducible_paper_codes/SantosJr and Parast (2026)/simulation/init.R')

# We assume that the model has been fitted
load(local_path('linear_N=1441.bkp'))


format.lag=function(lag){
  ifelse(lag==0,'No lags',
         ifelse(lag<=4,
                paste0(lag/4,' year'),
                paste0(lag/4,' years')))
}

# Point estimate for the total treatment effect
delta.est=(data.plot %>% filter(lag==-1))$y

data.total.effect=data.frame()
data.CPTE=data.frame()
data.plot4=data.frame()
data.time.test=data.frame()
data.LPTE=data.frame()

# The treatment effect should be 0 at baseline, but this is not built in the model
# Here we do a correction for that, thus avoiding the PTE being negative in the beginning do to finite-sample variations.
# We do this correction via a linear transformation defined by the matrix A
A=diag(T)
A[,1]=A[,1]-A[1,1]
delta.est=A%*%delta.est

# First, we compute the treatment effect for the marginal model, which will be used to compute the PTE for all conditional models
# Note that in the fitting script, at each bootstrap iteration, we used the same indexes to all models,
# which allows us to use the same marginal model for all conditional models
beta.marg=data.smp[[1]]$mt.smp[(structure.base$n-T):(structure.base$n-1),]
beta.marg=A%*%beta.marg[c(2:T,1),] # We need to reorder the sampled betas, similar to the reordering done for the point estimate in the fitting script.
data.total.effect=rbind(data.total.effect,
                 data.frame(time=1:T,
                            lag=-1,
                            mean=rowMeans(beta.marg),
                            median=rowQuantile(beta.marg,0.5),
                            icu=rowQuantile(beta.marg,0.95),
                            icl=rowQuantile(beta.marg,0.05)))
for(idx in 2:length(data.smp)){
  beta=data.smp[[idx]]$mt.smp[(structure.base$n-T):(structure.base$n-1),]
  beta=A%*%beta[c(2:T,1),]
  cat(paste0(idx,'                     \r'))

  # Total treatment effect
  data.total.effect=rbind(data.total.effect,
                   data.frame(time=1:T,
                              lag=data.smp[[idx]]$lag.max,
                              mean=rowMeans(beta),
                              median=rowQuantile(beta,0.5),
                              icu=rowQuantile(beta,0.95),
                              icl=rowQuantile(beta,0.05)))

  # LPTE
  pte.cur=1-beta/beta.marg
  data.LPTE=rbind(data.LPTE,
                   data.frame(time=1:T,
                              lag=data.smp[[idx]]$lag.max,
                              mean=rowMeans(pte.cur),
                              median=rowQuantile(pte.cur,0.5),
                              icu=rowQuantile(pte.cur,0.95),
                              icl=rowQuantile(pte.cur,0.05)))

  # CPTE
  pte.cur=1-colCumSums(beta)/colCumSums(beta.marg)
  data.CPTE=rbind(data.CPTE,
                   data.frame(time=1:T,
                              lag=data.smp[[idx]]$lag.max,
                              mean=rowMeans(pte.cur),
                              median=rowQuantile(pte.cur,0.5),
                              icu=rowQuantile(pte.cur,0.95),
                              icl=rowQuantile(pte.cur,0.05)))

  # Hyphotesis test
  delta.R.est=data.smp[[idx]]$delta.est-data.smp[[idx]]$delta.est[1]
  pte.est=1-sum(delta.R.est)/sum(delta.est)

  h=(delta.R.est-(1-pte.est)*delta.est)[-1]

  h.smp=(beta-(1-matrix(pte.cur,nrow=10,ncol=2000,byrow=TRUE))*beta.marg)[-1,]
  var.h=var(t(h.smp))
  sd.h=diag(1/sqrt(diag(var.h)))

  stat.smp=colQuantile(rmvnorm(50000,rep(0,T-1),var.h %>% cov2cor) %>% abs,1)
  stat.smp %>% quantile(0.95)

  stat=abs(sd.h%*%h) %>% max()


  data.time.test=rbind(data.time.test,
                   data.frame(lag=data.smp[[idx]]$lag.max,
                              crit.val=stat.smp %>% quantile(0.95),
                              stat=stat,
                              p.value=mean(stat.smp>stat)))
}

# The p-value for the time homogeneity test
ggplot(data.time.test,aes(x=lag))+
  geom_line(aes(y=p.value))+
  geom_hline(yintercept=0.05,linetype='dashed')+
  scale_y_continuous(breaks=c(seq.int(0,1,0.25),0.05))+
  ylab('p-value')+
  scale_x_continuous('Maximum surrogate lag (in years)',labels=~./4,breaks=seq.int(0,36,8))+
  coord_cartesian(ylim=c(0,1))+
  theme_custom
ggsave(local_path('figures/fig3_applied_pvalue.pdf'),
       width=3,height=2.5)

# LPTE at key points (for visualization only)
(ggplot(data.LPTE %>%
          filter(lag %in% c(0,4,12,20,28)) %>%
          mutate(lag=format.lag(lag) %>% factor(.,unique(.))) %>%
          group_by(lag),
        aes(x=(time-1)*4+1,
            color=lag,
            fill=lag))+
    geom_line(aes(y=median))+
    geom_hline(yintercept=0,linetype='dashed')+
    geom_hline(yintercept=1,linetype='dashed')+
    scale_x_continuous('Time (in years)',labels=~./4,breaks=seq.int(0,36,8))+
    ylab('CPTE')+
    scale_color_hue('K')+
    scale_fill_hue('K')+
    coord_cartesian(ylim=c(0,1.5))+
    theme_custom) %>%
  (plotly::ggplotly)

# CPTE at key points
(ggplot(data.CPTE %>%
          filter(lag %in% c(0,4,12,20,28)) %>%
          mutate(lag=format.lag(lag) %>% factor(.,unique(.))) %>%
          group_by(lag),
        aes(x=(time-1)*4+1,color=lag,fill=lag))+
    geom_line(aes(y=median))+
    scale_x_continuous('Time (in years)',labels=~./4,breaks=seq.int(0,36,8))+
    scale_y_continuous('CPTE',breaks=seq.int(0,1,0.5),expand=FALSE)+
    coord_cartesian(ylim=c(0,1))+
    scale_color_hue('K')+
    scale_fill_hue('K')+
    theme_custom+
    theme(legend.position='left',
          panel.grid.major.y=element_line(color = 'lightgray')))
ggsave(local_path('figures/fig1_applied_cpte.pdf'),
       width=3.55,height=2)
