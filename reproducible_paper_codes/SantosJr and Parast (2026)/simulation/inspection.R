source('reproducible_paper_codes/SantosJr and Parast (2026)/simulation/init.R')

# Auxiliary function to make the study names more intuitive.
convert.case=function(x){
  ifelse(x=='Study 1 - Case 1','Monotone',
         ifelse(x=='Study 1 - Case 2','Parabole',
                'Random Walk'))
}

# File list - assumes you already computed the simulations
file.list=paste0(local_path('data/'),list.files(local_path('data/')))

B=data.frame()

for(file in file.list){
  A=read_csv(file,progress = FALSE,show_col_types = FALSE) %>%
    mutate(rep=paste0(file,rep))
  B=rbind(B,A)
}

B=B %>%
  factor(.,unique(.)) %>%
  arrange(case,rep,N,time,PTE,alpha,phi,sign,method)
B$case %>% is.na %>% sum

# Evaluating Study 1: time homogeneous PTE
# In this comparisson we compare our method with some of its natural competitors
# Through different scenarios, each corresponding to different dynamics for the treatment effect
# The PTE stays constant in every scenario.
ref.case='Study 1'

# Evaluating bias. We set reference duration as 40 (the maximum duration)
# The behavior of the bias is the same for all durations, but they have different levels
# If we combine the different durations, the graph becomes awkward because we're mixing groups centered around different points
# Qualitatively, the behavior is the same, so we just show the results for one duration
# At first, using the maximum duration may look like looking at the best scenario (more measurements),
# but in reality it is the other way around: the longer study is the hardest case
# This is so because each additional time point adds more effects to estimate, while the number of patients is the same for all durations
ggplot(B %>%
         filter(time==40) %>%
         filter(grepl(ref.case,case)) %>%
         mutate(Sample_size=N*time,
                case=convert.case(case)) %>%
         group_by(N,PTE,method,case) %>%
         filter(abs(point-PTE)<quantile(abs(point-PTE),0.99)) %>%
         summarize(point=mean(abs(point-PTE),na.rm=TRUE))) +
  geom_point(aes(x=N,y=point,color=method),size=1)+
  geom_hline(yintercept=0,linetype='dashed')+
  scale_color_hue('')+
  xlab('Sample size')+
  scale_y_continuous('Average bias',
                     breaks=function(x){
                       seq.int(0,x[2],l=3)
                     },
                     limits = function(x){
                       c(0,0.1*ceiling(x[2]/0.1))
                     },expand=c(0.05,0,0,0))+
  facet_grid(case~paste0('PTE=',PTE),scale='free')+
  theme_custom+
  theme(legend.position='bottom',axis.text.x = element_text(angle=90),
        panel.grid = element_blank())
ggsave(local_path('figures/fig2_simul1_bias.pdf'),
       width=6,height=4)

# Evaluating the coverage
ggplot(B %>%
         filter(time==40) %>%
         filter(grepl(ref.case,case)) %>%
         mutate(Sample_size=N*time,
                case=convert.case(case)) %>%
         group_by(N,PTE,method,case) %>%
         filter(q025>quantile(q025,0.01),q975<quantile(q975,0.99)) %>%
         summarize(point=mean((PTE>=q025) & (PTE<=q975)))) +
  geom_point(aes(x=N,y=point,color=method),size=1)+
  geom_hline(yintercept=0.95,linetype='dashed')+
  scale_color_hue('')+
  xlab('Sample size')+
  scale_y_continuous('Empirical coverage',
                     breaks=function(x){
                       seq.int(x[1],x[2]-0.01,l=3)
                     },
                     limits = function(x){
                       c(0.05*floor(x[1]/0.05),0.05*ceiling(x[2]/0.05)+0.01)
                     },expand=c(0,0,0,0))+
  facet_grid(case~paste0('PTE=',PTE))+
  theme_custom+
  theme(legend.position='bottom',axis.text.x = element_text(angle=90),
        panel.grid = element_blank())
ggsave(local_path('figures/fig3_simul1_coverage.pdf'),
       width=6,height=4)

# Evaluating the width of the confidence interval.
ggplot(B %>%
         filter(time==40) %>%
         filter(grepl(ref.case,case)) %>%
         mutate(Sample_size=N*time,
                case=convert.case(case)) %>%
         group_by(N,PTE,method,case) %>%
         filter(q025>quantile(q025,0.01),q975<quantile(q975,0.99)) %>%
         summarize(point=mean(q975-q025))) +
  geom_point(aes(x=N,y=point,color=method),size=1)+
  geom_hline(yintercept=0,linetype='dashed')+
  scale_color_hue('')+
  xlab('Sample size')+
  scale_y_continuous('C.I. width',
                     breaks=function(x){
                       if(x[2]>1){
                         seq.int(0,ceiling(x[2])-1,l=3)
                       }else{
                         seq.int(0,x[2],l=3)
                       }
                     },
                     limits = function(x){
                       if(x[2]>1){
                         c(0,2*ceiling(x[2]/2))
                       }else{
                         c(0,ceiling(x[2]))
                       }})+
  facet_grid(case~paste0('PTE=',PTE),scales='free')+
  theme_custom+
  theme(legend.position='bottom',axis.text.x = element_text(angle=90),
        panel.grid = element_blank())
ggsave(local_path('figures/fig4_simul1_width.pdf'),
       width=6,height=4)

# Evaluating the power of the test for the validity of the PTE.
# The reference threshold is 0.75 and the significance level is 0.05
ggplot(B %>%
         filter(time==40) %>%
         filter(!is.na(point)) %>%
         filter(grepl(ref.case,case)) %>%
         mutate(Sample_size=N*time,
                case=convert.case(case)) %>%
         group_by(N,PTE,method,case) %>%
         summarize(point=mean(q050>0.75))) +
  geom_point(aes(x=N,y=point,color=method),size=1)+
  geom_hline(yintercept=0.05,linetype='dashed')+
  scale_color_hue('')+
  xlab('Sample size')+
  scale_y_continuous('Prop. of rejections',breaks=seq.int(0,1,0.5))+
  facet_grid(case~paste0('PTE=',PTE))+
  theme_custom+
  theme(legend.position='bottom',axis.text.x = element_text(angle=90),
        panel.grid = element_blank())
ggsave(local_path('figures/fig5_simul1_power.pdf'),
       width=6,height=4)

# Case study 2: time homogeneity test
# Goal for this part is to
# (1) check if the proposed hyphotesis test has a empirical Type I error close to the nominal value
# (2) compare the MSD and Wald-based tests.
# We consider the test based on the SSM estimator,
# as the difference between different models is not important here

ggplot(B %>%
         filter(method=='SSM') %>%
         filter(grepl(ref.case,'Study 1 - Case 1') | grepl('Study 2',case)) %>%
         mutate(case=str_split_i(case,' - ',2) %>%
                  str_replace_all('Case','Scenario')) %>%
         mutate(MSD=MSD>MSD.crit,
                Wald=Wald>qchisq(0.95,time),
                # N=as.factor(N),
                time=paste0(time/4,' years') %>% factor(.,levels=unique(.))) %>%
         select(case,time,N,MSD,Wald) %>%
         pivot_longer(4:5) %>%
         group_by(case,time,N,name) %>%
         summarize(mean=mean(value,na.rm=TRUE),
                   sd=sd(value,na.rm=TRUE),
                   count=length(value))
)+
  geom_point(aes(x=N,y=mean,color=name))+
  geom_hline(yintercept=0.05,linetype='dashed')+
  scale_color_manual('',values=c('#5555ff','#ff5555'))+
  scale_x_continuous(TeX('Sample size'),breaks=seq.int(0,300,100))+
  scale_y_continuous('Proportion of rejections',breaks=seq.int(0,1,0.5))+
  facet_grid(case~time)+
  theme_custom+
  theme(axis.text.x = element_text(angle=90),
        axis.ticks.x = ,
        legend.position = 'bottom')
ggsave(local_path('figures/fig6_simul2_power.pdf'),
       width=6,height=6)
