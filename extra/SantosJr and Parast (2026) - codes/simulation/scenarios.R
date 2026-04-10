source('extra/SantosJr and Parast (2026) - codes/simulation/init.R')

# First, a visualization of the treatment effects used in Study 1

T=2000
data.plot=data.frame()
data.plot=rbind(data.plot,
                data.frame(x=seq.int(1,20,l=T),
                           y=tanh((seq.int(1,20,l=T)-1)*2.64665/(20-1)),
                           color='Monotone',
                           group='Monotone'))
data.plot=rbind(data.plot,
                data.frame(x=seq.int(1,20,l=T),
                           y=-((seq.int(1,20,l=T)-10.5)/(9.5))**2+1,
                           color='Parabole',
                           group='Parabole'))

T=20
for(i in 1:1){
  y=cumsum(rnorm(T,0,0.5/sqrt(T)))
  y=y-min(y)
  data.plot=rbind(data.plot,
                  data.frame(x=seq.int(1,20,l=T),
                             y=y,
                             color='Random Walk',
                             group=i))
}

ggplot(data.plot)+
  geom_line(aes(x=x,y=y,color=color,group=group))+
  scale_y_continuous('Treatment effect' )+
  scale_x_continuous('')+
  scale_color_manual('',values=c('#dd6666','#66dd66','#6666dd99'))+
  theme_custom+
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),axis.text.y = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'bottom')

ggsave(file=local_path('figures/fig1_effect_dynamic.pdf'),
       width=6,height=3)


# Now, for Study 2

N=1000
# Scenario 1
T=10
x=seq.int(1,T,l=N+1)
alpha=0.01
direct=0.1*tanh((x-1)*2.64665/(20-1))
indirect=0.9*tanh((x-1)*2.64665/(20-1))
total=direct+indirect

direct=direct/max(total)
indirect=indirect/max(total)
total=total/max(total)

local.pte=indirect/total
cum.pte=cumsum(indirect)/cumsum(total)

ggplot()+
  geom_line(aes(x=x,y=local.pte),color='black')+
  geom_line(aes(x=x,y=cum.pte),color='black',linetype='dashed')+
  geom_line(aes(x=x,y=total),color='gray')+
  geom_hline(yintercept=0.75,linetype='dashed')+
  scale_y_continuous('',breaks=-1)+scale_x_continuous('',breaks=-1)+
  coord_cartesian(ylim=c(0,1))+
  theme_custom

ggsave(local_path('simul8_scenario_1.pdf'),
       width=2,height=2)

# Scenario 2
T=10
x=seq.int(1,T,l=N+1)
alpha=0.1
direct=(sin(3*(x)*pi/(3*T)))**10
indirect=rep(alpha,N+1)
total=direct+indirect

direct=direct/max(total)
indirect=indirect/max(total)
total=total/max(total)

local.pte=indirect/total
cum.pte=cumsum(indirect)/cumsum(total)

ggplot()+
  geom_line(aes(x=x,y=local.pte),color='black')+
  geom_line(aes(x=x,y=cum.pte),color='black',linetype='dashed')+
  geom_line(aes(x=x,y=total),color='gray')+
  geom_hline(yintercept=0.75,linetype='dashed')+
  scale_y_continuous('',breaks=-1)+scale_x_continuous('',breaks=-1)+
  coord_cartesian(ylim=c(0,1))+
  theme_custom

ggsave(local_path('simul8_scenario_2.pdf'),
       width=2,height=2)

# Scenario 3
T=10
x=seq.int(1,T,l=N+1)
alpha=0.1
direct=rep(alpha,N+1)
indirect=(sin(3*(x)*pi/(3*T)))**10
total=direct+indirect

direct=direct/max(total)
indirect=indirect/max(total)
total=total/max(total)

local.pte=indirect/total
cum.pte=cumsum(indirect)/cumsum(total)

ggplot()+
  geom_line(aes(x=x,y=local.pte),color='black')+
  geom_line(aes(x=x,y=cum.pte),color='black',linetype='dashed')+
  geom_line(aes(x=x,y=total),color='gray')+
  geom_hline(yintercept=0.75,linetype='dashed')+
  scale_y_continuous('',breaks=-1)+scale_x_continuous('',breaks=-1)+
  coord_cartesian(ylim=c(0,1))+
  theme_custom
ggsave(local_path('simul8_scenario_3.pdf'),
       width=2,height=2)


# Scenario 4
T=10
x=seq.int(1,T,l=N+1)
alpha=0.2
direct=rep(alpha,N+1)
indirect=tanh(3*seq.int(0.01,1,l=N+1))/tanh(3)
total=direct+indirect

direct=direct/max(total)
indirect=indirect/max(total)
total=total/max(total)

local.pte=indirect/total
cum.pte=cumsum(indirect)/cumsum(total)

ggplot()+
  geom_line(aes(x=x,y=local.pte),color='black')+
  geom_line(aes(x=x,y=cum.pte),color='black',linetype='dashed')+
  geom_line(aes(x=x,y=total),color='gray')+
  geom_hline(yintercept=0.75,linetype='dashed')+
  scale_y_continuous('',breaks=-1)+scale_x_continuous('',breaks=-1)+
  coord_cartesian(ylim=c(0,1))+
  theme_custom
ggsave(local_path('simul8_scenario_4.pdf'),
       width=2,height=2)


# Scenario 5
T=10
x=seq.int(1,T,l=N+1)
alpha=0.2
direct=tanh(3*seq.int(0.01,1,l=N+1))/tanh(3)
indirect=rep(alpha,N+1)
total=direct+indirect

direct=direct/max(total)
indirect=indirect/max(total)
total=total/max(total)

local.pte=indirect/total
cum.pte=cumsum(indirect)/cumsum(total)

ggplot()+
  geom_line(aes(x=x,y=local.pte),color='black')+
  geom_line(aes(x=x,y=cum.pte),color='black',linetype='dashed')+
  geom_line(aes(x=x,y=total),color='gray')+
  geom_hline(yintercept=0.75,linetype='dashed')+
  scale_y_continuous('',breaks=-1)+
  scale_x_continuous('',breaks=-1)+
  coord_cartesian(ylim=c(0,1))+
  theme_custom

ggsave(local_path('simul8_scenario_5.pdf'),
       width=2,height=2)
