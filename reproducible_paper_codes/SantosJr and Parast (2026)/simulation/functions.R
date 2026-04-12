# ------------------------------------------------------------------------------
# Auxiliary functions for the simulation study
#
# This file defines helper functions used to run the simulation experiments
# reported in the paper. Its main entry point is `run_set()`, which generates
# data under the study scenarios, fits the competing methods, computes summary
# statistics and test quantities, and writes the results to an output file.
#
# This is an auxiliary script: it is intended to be sourced by a higher-level
# driver file rather than executed on its own.
#
# Dependencies:
#   This file assumes that several objects and helper functions are defined
#   elsewhere in the project, including the simulation settings, scenario
#   definitions, model hyperparameters, bootstrap size, and output connection.
#   Those values are specified in the simulation/init.R file.
# ------------------------------------------------------------------------------

# Generates the observations.
generate_obs=function(treat.val.y,treat.val.sur,phi,df,sign,N,T){
  treat=c(rep(0,N/2),rep(0,N/2)) %>% sample(N)

  treat=c(rep(0,N/2),rep(1,N/2)) %>% sample(N)
  treat.mat=matrix(treat.val.sur,T,N)*matrix(treat,T,N,byrow=TRUE)

  S=matrix(NA,T,N)
  S[1,]=rnorm(N,0,sqrt(W/(1-phi1)))
  for(t in 2:T){
    S[t,]=phi1*S[t-1,]+rnorm(N,0,sqrt(W))
  }
  S=S+treat.mat

  treat.mat=matrix(treat.val.y,T,N)*matrix(treat,T,N,byrow=TRUE)
  Y=matrix(NA,T,N)
  V.mixed=rgamma(N,df/2,df/(2*V))
  Y[1,]=sqrt(V.mixed/(1-phi2))*sign*log(rgamma(N,phi,exp(digamma(phi))))
  for(t in 2:T){
    Y[t,]=phi2*Y[t-1,]+sqrt(V.mixed)*sign*log(rgamma(N,phi,exp(digamma(phi))))
  }
  Y=Y+matrix(S.effect,T,N)*S+treat.mat
  return(list(Y=Y,S=S,treat=treat,
              delta=treat.val.y+treat.val.sur,
              delta.R=treat.val.y,
              sign=sign))
}

# Generates data under the study scenarios, fits the competing methods, computes summary
# statistics and test quantities, and writes the results to an output file.
run_set=function(reps){
  count=0
  for(case in c(
    'Study 1 - Case 1','Study 1 - Case 2','Study 1 - Case 3',
    'Study 2 - Case 2','Study 2 - Case 3','Study 2 - Case 4','Study 2 - Case 5'
  )){
    for(PTE.true in pte.list[[case]]){
      # N is the number of patients in the study.
      for(N in N.list){
        # n.years is the duration of the study, measured in years.
        # In this hypothetical study, measurements are taken every 3 months,
        # so each patient has 4 observations per year.
        for(n.years in T.list){
          # phi controls the over dispersion of the errors (both the observational errors and the evolution errors).
          # It is analogous to the degrees of freedom of the student-t distribution.
          # We tested several values of phi and found no qualitative differences in the results.
          # For brevity, we report only the most challenging setting: the smallest tested value of phi.
          # It is important to note that, if phi is too small the outcome may have infinite variance,
          # violating a condition required by our method. We therefore use phi = 5, which remains a hard
          # case while still ensuring finite first and second moments.
          for(phi in c(5)){
            # alpha controls the asymmetry of the errors, both in the observation equation
            # and in the state evolution equation. Larger values produce a more symmetric
            # distribution, whereas smaller values imply greater asymmetry.
            # We tested several values of alpha and found no qualitative differences in the results.
            # For brevity, we report only the most challenging setting: the smallest tested value of alpha.
            for(alpha in c(15)){
              T=floor(n.years*freq) # Time series length
              labels=paste0('mu',formatC(1:N,flag = '0',width=3)) # Labels for the state-space model (for the kDGLM package)

              # Initializing the structural blocks for the state-space model (for the kDGLM package)
              loc.level=init_level(N,D=0.97)
              glob.G=init_G(labels,N,T,rho,900)
              loc.G=init_G('mu',1,T,rho,900*N)

              # Treatment index in the latent variables
              idx.G.cond=1:(T+1)
              idx.G.marg=1:(T)

              # Random effect index in the latent variables
              idx.levels.cond=(T+1)+1:N
              idx.levels.marg=(T)+1:N

              # For the OLS
              placeholder.G=rep(diag(T),N) %>% matrix(T,N*T)

              idx.ref=1:N

              # Placeholder for the bootstrap samples
              delta.smp=matrix(NA,N.bootstrap,T); delta.R.smp=matrix(NA,N.bootstrap,T);

              if(case=='Study 1 - Case 1'){
                # Study 1 - Case 1: Time homogeneous PTE; Monotonically increasing treatment effect.
                delta.true=total.effect[[case]]*sqrt(V/(1-phi1))
                delta.R.true=(1-PTE.true)*delta.true
                treat.y=delta.R.true
                treat.s=delta.true*PTE.true/S.effect

                treat.val=tanh(((1:T)-1)*2.64665/(20-1))
                treat.val.y=treat.y*treat.val
                treat.val.sur=treat.s*treat.val
              }else if(case=='Study 1 - Case 2'){
                # Study 1 - Case 2: Time homogeneous PTE; Increasing and then decreasing (parabole) treatment effect.
                delta.true=total.effect[[case]]*sqrt(V/(1-phi1))
                delta.R.true=(1-PTE.true)*delta.true
                treat.y=delta.R.true
                treat.s=delta.true*PTE.true/S.effect

                treat.val=-((1:T-10.5)/9.5)**2+1
                treat.val.y=treat.y*treat.val
                treat.val.sur=treat.s*treat.val
              }else if(case=='Study 1 - Case 3'){
                # Study 1 - Case 3: Time homogeneous PTE; Random walk treatment effect.
                delta.true=total.effect[[case]]*sqrt(V/(1-phi1))
                delta.R.true=(1-PTE.true)*delta.true
                treat.y=delta.R.true
                treat.s=delta.true*PTE.true/S.effect

                treat.val=cumsum(rnorm(T,0,0.1))
                treat.val=(treat.val-min(treat.val))
                treat.val.y=treat.y*treat.val
                treat.val.sur=treat.s*treat.val
              }else if(case=='Study 2 - Case 2'){
                # Study 2 - Case 2: Time heterogeneous PTE; Seasonal direct effect + constant indirect effect.
                treat.val.y=(sin(3*(1:T)*pi/(3*T)))**10
                treat.val.sur=0.1

                delta.true=total.effect[[case]]*sqrt(V/(1-phi1))
                treat.y=delta.true/max(treat.val.y+treat.val.sur)
                treat.val.y=treat.y*treat.val.y
                treat.val.sur=treat.y*treat.val.sur
              }else if(case=='Study 2 - Case 3'){
                # Study 2 - Case 3: Time heterogeneous PTE; Constant direct effect + seasonal indirect effect.
                treat.val.y=0.1
                treat.val.sur=(sin(3*(1:T)*pi/(3*T)))**10

                delta.true=total.effect[[case]]*sqrt(V/(1-phi1))
                treat.y=delta.true/max(treat.val.y+treat.val.sur)
                treat.val.y=treat.y*treat.val.y
                treat.val.sur=treat.y*treat.val.sur
              }else if(case=='Study 2 - Case 4'){
                # Study 2 - Case 4: Time heterogeneous PTE; Constant direct effect + Increasing indirect effect.
                treat.val.y=0.2
                treat.val.sur=tanh(2*seq.int(0.01,1,l=T))/tanh(2)

                delta.true=total.effect[[case]]*sqrt(V/(1-phi1))
                treat.y=delta.true/max(treat.val.y+treat.val.sur)
                treat.val.y=treat.y*treat.val.y
                treat.val.sur=treat.y*treat.val.sur
              }else if(case=='Study 2 - Case 5'){
                # Study 2 - Case 5: Time heterogeneous PTE; Increasing direct effect + constant indirect effect.
                treat.val.y=tanh(2*seq.int(0.01,1,l=T))/tanh(2)
                treat.val.sur=0.2

                delta.true=total.effect[[case]]*sqrt(V/(1-phi1))
                treat.y=delta.true/max(treat.val.y+treat.val.sur)
                treat.val.y=treat.y*treat.val.y
                treat.val.sur=treat.y*treat.val.sur
              }else{
                stop('Invalid treatment dynamic')
              }

              for(rep in reps){
                sign=(-1)**rep # Left or right skewed
                obs=generate_obs(treat.val.y,treat.val.sur,df=alpha,phi=phi,sign=(-1)**rep,N=N,T=T)
                Y=obs$Y;S=obs$S;treat=obs$treat;delta=obs$delta;delta.R=obs$delta.R;

                ############ DLM ############
                glob.G$FF[1,,] = matrix(treat,N,T)

                # First, fit the full model to compute W (variance of the evolution errors).
                naive.dlm=fit_full_normal(Y,S,N,T,loc.level,glob.G)

                # The variance of the observational noise. This is more of a placeholder.
                # Can be useful if one wants to give different weights to each patient,
                # Otherwise it does not affect the estimates.
                V.cond=1
                V.marg=1

                # The variance of the evolution noise.
                # This is computed using a discount strategy.
                W.cond=naive.dlm$cond$W[(1:N)+(T),(1:N)+(T),]
                W.marg=naive.dlm$marg$W[(1:N)+(T),(1:N)+(T),]

                W.cond.list=sapply(1:N,function(x){W.cond[x,x,]})
                W.marg.list=sapply(1:N,function(x){W.marg[x,x,]})

                # Prior mean and variance. This is equivalent to a L2 penalization.
                cond.mean=0;cond.var=9;
                marg.mean=0;marg.var=9;

                # Placeholders for the contribution of each patient.
                St.marg=array(NA,c(T+1,T+1,N)) # Precision matrices
                St.cond=array(NA,c(T+2,T+2,N)) # Precision matrices
                mt.cond=array(NA,c(T+2,N));mt.marg=array(NA,c(T+1,N));

                # This loop computes the contribution of each patient, which will later be combined to generate the bootstrap samples.
                for(j in 1:N){
                  out=Y[,j,drop=FALSE]
                  loc.G$FF[1,1,]=treat[j]/sqrt(V.cond)
                  models.cond=fit_single_normal(out,S[,j,drop=FALSE],
                                                V.cond,W.cond.list[,j],
                                                a1=cond.mean,R1=cond.var,N,T,
                                                struct.G = loc.G,
                                                sample=FALSE)
                  St.cond[,,j]=St=models.cond$Ct[,,T] %>% ginv
                  mt.cond[,j]=St%*%models.cond$mt[,T]

                  loc.G$FF[1,1,]=treat[j]/sqrt(V.marg)
                  models.marg=fit_single_normal(out,S=NA,
                                                V.marg,W.marg.list[,j],
                                                a1=marg.mean,R1=marg.var,N,T,
                                                struct.G = loc.G,
                                                sample=FALSE)
                  St.marg[,,j]=St=models.marg$Ct[,,T] %>% ginv
                  mt.marg[,j]=St%*%models.marg$mt[,T]
                }

                # Computing the point estimates
                idx=1:N

                # Conditional model
                vals.diag=St.cond[(T+1),(T+1),]
                A_11=rowSums(St.cond[-(T+1),-(T+1),],dims=2)
                A_12=St.cond[(T+1),-(T+1),]
                S_11=A_11-tcrossprod(sweep(A_12,2,vals.diag,FUN='/'),A_12)

                a_1=rowSums(mt.cond[-(T+1),])
                a_2=mt.cond[(T+1),]/vals.diag

                mt.full.cond=solve(S_11,a_1-A_12%*%a_2)

                # Marginal model
                vals.diag=St.marg[(T+1),(T+1),]
                A_11=rowSums(St.marg[-(T+1),-(T+1),],dims=2)
                A_12=St.marg[(T+1),-(T+1),]
                S_11=A_11-tcrossprod(sweep(A_12,2,vals.diag,FUN='/'),A_12)

                a_1=rowSums(mt.marg[-(T+1),])
                a_2=mt.marg[(T+1),]/vals.diag

                mt.full.marg=solve(S_11,a_1-A_12%*%a_2)
                for(i in 1:N.bootstrap){
                  # Making sure each bootstrap sample keeps the balance between treatment and control.
                  idx.treat=sample.int(N/2,N/2,replace=TRUE)
                  idx.control=sample.int(N/2,N/2,replace=TRUE)
                  idx=c(which(treat==1)[idx.treat],which(treat==0)[idx.control])

                  # Conditional model

                  vals.diag=St.cond[(T+1),(T+1),idx]
                  A_11=rowSums(St.cond[-(T+1),-(T+1),idx],dims=2)
                  A_12=St.cond[(T+1),-(T+1),idx]
                  S_11=A_11-tcrossprod(sweep(A_12,2,vals.diag,FUN='/'),A_12)

                  a_1=rowSums(mt.cond[-(T+1),idx])
                  a_2=mt.cond[(T+1),idx]/vals.diag

                  mt.est.cond=solve(S_11,a_1-A_12%*%a_2)

                  # Marginal model
                  vals.diag=St.marg[(T+1),(T+1),idx]
                  A_11=rowSums(St.marg[-(T+1),-(T+1),idx],dims=2)
                  A_12=St.marg[(T+1),-(T+1),idx]
                  S_11=A_11-tcrossprod(sweep(A_12,2,vals.diag,FUN='/'),A_12)

                  a_1=rowSums(mt.marg[-(T+1),idx])
                  a_2=mt.marg[(T+1),idx]/vals.diag

                  mt.est.marg=solve(S_11,a_1-A_12%*%a_2)


                  delta.smp[i,]=mt.est.marg[c(2:T,1)]
                  delta.R.smp[i,]=mt.est.cond[c(2:T,1)]
                }
                pte.smp=1-rowSums(delta.R.smp)/rowSums(delta.smp)
                delta.est=mt.full.marg[c(2:T,1)]
                delta.R.est=mt.full.cond[c(2:T,1)]

                pte=1-sum(delta.R.est)/sum(delta.est)

                q050.dlm=quantile(pte.smp,0.05)
                q025.dlm=quantile(pte.smp,0.025)
                q950.dlm=quantile(pte.smp,0.95)
                q975.dlm=quantile(pte.smp,0.975)
                point.dlm=pte
                mean.dlm=mean(pte.smp)
                sd.dlm=sd(pte.smp)
                median.dlm=quantile(pte.smp,0.5)

                test.stats=calc_tests(delta.R.smp,delta.smp,pte.smp,
                                      delta.R.est,pte,delta.est)

                Wald.dlm=test.stats$Wald
                Wald.dlm.crit=test.stats$Wald.test
                MSD.dlm=test.stats$MSD
                MSD.dlm.crit=test.stats$MSD.test

                ############ OLS ############

                X.full=rbind(1,c(S),
                             placeholder.G*matrix(rep(treat,each=T),T,N*T,byrow=TRUE))
                Y.full=c(Y)

                X.idx=X.full
                Y.idx=Y.full
                A=X.idx %>% tcrossprod(.,.) %>% chol
                out=X.idx%*%Y.idx
                beta.cond=backsolve(A,forwardsolve(A,out,upper.tri = TRUE,transpose = TRUE))

                X.idx=X.idx[-2,]
                A=X.idx %>% tcrossprod(.,.) %>% chol
                out=X.idx%*%Y.idx
                beta.marg=backsolve(A,forwardsolve(A,out,upper.tri = TRUE,transpose = TRUE))

                delta.R.est=beta.cond[1:T+2]
                delta.est=beta.marg[1:T+1]

                for(i in 1:N.bootstrap){
                  idx.treat=sample.int(N/2,N/2,replace=TRUE)
                  idx.control=sample.int(N/2,N/2,replace=TRUE)
                  idx=c(which(treat==1)[idx.treat],which(treat==0)[idx.control])

                  idx=rep((idx-1)*T,each=T)+rep(1:T,times=N)
                  X.idx=X.full[,idx]
                  Y.idx=Y.full[idx]
                  A=X.idx %>% tcrossprod(.,.) %>% chol
                  out=X.idx%*%Y.idx
                  beta.cond=backsolve(A,forwardsolve(A,out,upper.tri = TRUE,transpose = TRUE))

                  X.idx=X.idx[-2,]
                  A=X.idx %>% tcrossprod(.,.) %>% chol
                  out=X.idx%*%Y.idx
                  beta.marg=backsolve(A,forwardsolve(A,out,upper.tri = TRUE,transpose = TRUE))

                  delta.R.smp[i,]=beta.cond[1:T+2]
                  delta.smp[i,]=beta.marg[1:T+1]
                }
                pte.smp=1-rowSums(delta.R.smp)/rowSums(delta.smp)

                pte=1-sum(delta.R.est)/sum(delta.est)

                q050.ols=quantile(pte.smp,0.05)
                q025.ols=quantile(pte.smp,0.025)
                q950.ols=quantile(pte.smp,0.95)
                q975.ols=quantile(pte.smp,0.975)
                point.ols=pte
                mean.ols=mean(pte.smp)
                sd.ols=sd(pte.smp)
                median.ols=quantile(pte.smp,0.5)


                test.stats=calc_tests(delta.R.smp,delta.smp,pte.smp,
                                      delta.R.est,pte,delta.est)

                Wald.ols=test.stats$Wald
                Wald.ols.crit=test.stats$Wald.test
                MSD.ols=test.stats$MSD
                MSD.ols.crit=test.stats$MSD.test

                ############ Diff ############
                # Similar to the OLS, but using only the difference from baseline.

                X.full=rbind(1,S[T,]-S[1,],treat)
                Y.full=Y[T,]-Y[1,]

                X.idx=X.full
                Y.idx=Y.full
                A=X.idx %>% tcrossprod(.,.) %>% chol
                out=X.idx%*%Y.idx
                beta.cond=backsolve(A,forwardsolve(A,out,upper.tri = TRUE,transpose = TRUE))

                X.idx=X.idx[-2,]
                A=X.idx %>% tcrossprod(.,.) %>% chol
                out=X.idx%*%Y.idx
                beta.marg=backsolve(A,forwardsolve(A,out,upper.tri = TRUE,transpose = TRUE))

                delta.R.est=beta.cond[3]
                delta.est=beta.marg[2]

                for(i in 1:N.bootstrap){
                  idx.treat=sample.int(N/2,N/2,replace=TRUE)
                  idx.control=sample.int(N/2,N/2,replace=TRUE)
                  idx=c(which(treat==1)[idx.treat],which(treat==0)[idx.control])

                  X.idx=X.full[,idx]
                  Y.idx=Y.full[idx]
                  A=X.idx %>% tcrossprod(.,.) %>% chol
                  out=X.idx%*%Y.idx
                  beta.cond=backsolve(A,forwardsolve(A,out,upper.tri = TRUE,transpose = TRUE))

                  X.idx=X.idx[-2,]
                  A=X.idx %>% tcrossprod(.,.) %>% chol
                  out=X.idx%*%Y.idx
                  beta.marg=backsolve(A,forwardsolve(A,out,upper.tri = TRUE,transpose = TRUE))

                  delta.R.smp[i,]=beta.cond[3]
                  delta.smp[i,]=beta.marg[2]
                }
                pte.smp=1-rowSums(delta.R.smp)/rowSums(delta.smp)
                pte=1-sum(delta.R.est)/sum(delta.est)

                q050.diff=quantile(pte.smp,0.05)
                q025.diff=quantile(pte.smp,0.025)
                q950.diff=quantile(pte.smp,0.95)
                q975.diff=quantile(pte.smp,0.975)
                point.diff=pte
                mean.diff=mean(pte.smp)
                sd.diff=sd(pte.smp)
                median.diff=quantile(pte.smp,0.5)

                Wald.diff=NA
                Wald.diff.crit=NA
                MSD.diff.crit=NA
                MSD.diff=NA

                ############ LMM ############

                data.X=data.frame(outcome=c(Y),surrogate=c(S),
                                  time=as.factor(rep(1:T,N)),
                                  patient=as.factor(rep(1:N,each=T)),
                                  treatment=as.factor(rep(treat,each=T)))

                fm.marg <- lmer(outcome ~ time:treatment + time + (1 | patient),
                                data = data.X)
                fm.cond <- lmer(outcome ~ surrogate + time:treatment + time + (1 | patient),
                                data = data.X)

                delta.R.mean=coef(fm.cond)$patient[1,1:T+T+1] %>% as.numeric
                delta.mean=coef(fm.marg)$patient[1,1:T+T] %>% as.numeric
                delta.R.var=vcovCR(fm.cond, cluster = data.X$patient, type = "CR2")[1:T+T+1,1:T+T+1] %>% as.matrix
                delta.var=vcovCR(fm.marg, cluster = data.X$patient, type = "CR2")[1:T+T,1:T+T] %>% as.matrix

                delta.R.smp=rmvnorm(N.bootstrap,delta.R.mean,delta.R.var) %>% t
                delta.smp=rmvnorm(N.bootstrap,delta.mean,delta.var) %>% t
                pte.smp=1-rowSums(delta.R.smp)/rowSums(delta.smp)

                local.pte=1-(delta.R.smp)/(delta.smp)

                delta.R.est=colMeans(delta.R.smp)
                delta.est=colMeans(delta.smp)
                pte=1-sum(delta.R.est)/sum(delta.est)

                q050.lmm=quantile(pte.smp,0.05)
                q025.lmm=quantile(pte.smp,0.025)
                q950.lmm=quantile(pte.smp,0.95)
                q975.lmm=quantile(pte.smp,0.975)
                point.lmm=pte
                mean.lmm=mean(pte.smp)
                sd.lmm=sd(pte.smp)
                median.lmm=quantile(pte.smp,0.5)

                test.stats=calc_tests(delta.R.smp,delta.smp,pte.smp,
                                      delta.R.est,pte,delta.est)

                Wald.lmm=test.stats$Wald
                Wald.lmm.crit=test.stats$Wald.test
                MSD.lmm=test.stats$MSD
                MSD.lmm.crit=test.stats$MSD.test

                ############ GEE ############

                gee.marg <- geeglm(outcome ~ time:treatment + time,
                                   id=patient,
                                   corstr = "independence",
                                   data = data.X)
                gee.cond <- geeglm(outcome ~ surrogate + time:treatment + time,
                                   id=patient,
                                   corstr = "independence",
                                   data = data.X)

                delta.R.mean=coef(gee.cond)[1:T+T+1] %>% as.numeric
                delta.mean=coef(gee.marg)[1:T+T] %>% as.numeric
                delta.R.var=summary(gee.cond)$cov.scaled[1:T+T+1,1:T+T+1] %>% as.matrix
                delta.var=summary(gee.marg)$cov.scaled[1:T+T,1:T+T] %>% as.matrix

                delta.R.smp=rmvnorm(N.bootstrap,delta.R.mean,delta.R.var) %>% t
                delta.smp=rmvnorm(N.bootstrap,delta.mean,delta.var) %>% t

                pte.smp=1-rowSums(delta.R.smp)/rowSums(delta.smp)


                local.pte=1-(delta.R.smp)/(delta.smp)

                delta.R.est=colMeans(delta.R.smp)
                delta.est=colMeans(delta.smp)
                pte=1-sum(delta.R.est)/sum(delta.est)

                q050.gee=quantile(pte.smp,0.05)
                q025.gee=quantile(pte.smp,0.025)
                q950.gee=quantile(pte.smp,0.95)
                q975.gee=quantile(pte.smp,0.975)
                point.gee=pte
                mean.gee=mean(pte.smp)
                sd.gee=sd(pte.smp)
                median.gee=quantile(pte.smp,0.5)


                test.stats=calc_tests(delta.R.smp,delta.smp,pte.smp,
                                      delta.R.est,pte,delta.est)

                Wald.gee=test.stats$Wald
                Wald.gee.crit=test.stats$Wald.test
                MSD.gee=test.stats$MSD
                MSD.gee.crit=test.stats$MSD.test


                ############ Saving ############
                vals=data.frame(
                  case=rep(case,5),
                  rep=rep(rep,5), N=rep(N,5), time=rep(T,5), PTE=rep(PTE.true,5),
                  alpha=rep(alpha,5),phi=rep(phi,5),sign=rep(sign,5),
                  method=c('DLM','Diff','OLS','LMM','GEE'),
                  q050      =c(q050.dlm,q050.diff,q050.ols,q050.lmm,q050.gee),
                  q025      =c(q025.dlm,q025.diff,q025.ols,q025.lmm,q025.gee),
                  q950      =c(q950.dlm,q950.diff,q950.ols,q950.lmm,q950.gee),
                  q975      =c(q975.dlm,q975.diff,q975.ols,q975.lmm,q975.gee),
                  point     =c(point.dlm,point.diff,point.ols,point.lmm,point.gee),
                  mean      =c(mean.dlm,mean.diff,mean.ols,mean.lmm,mean.gee),
                  median    =c(median.dlm,median.diff,median.ols,median.lmm,median.gee),
                  sd        =c(sd.dlm,sd.diff,sd.ols,sd.lmm,sd.gee),
                  Wald      =c(Wald.dlm,Wald.diff,Wald.ols,Wald.lmm,Wald.gee),
                  Wald.crit =c(Wald.dlm.crit,Wald.diff.crit,Wald.ols.crit,Wald.lmm.crit,Wald.gee.crit),
                  MSD       =c(MSD.dlm,MSD.diff,MSD.ols,MSD.lmm,MSD.gee),
                  MSD.crit  =c(MSD.dlm.crit,MSD.diff.crit,MSD.ols.crit,MSD.lmm.crit,MSD.gee.crit)
                )
                if(count==0){
                  vals %>% names %>%
                    paste0(collapse=',') %>%
                    writeLines(con=file)
                }
                count=count+1
                for(i in 1:dim(vals)[1]){
                  vals[i,] %>%
                    paste0(collapse=',') %>%
                    writeLines(con=file)
                }
              }
            }
          }
        }
      }
    }
  }
}



# Fits a model for a single patient
fit_single_normal=function(Y,S,V,H1,N,T,struct.G,a1=0,R1=900,sample=FALSE){
  args1=list(smooth=FALSE,safe.mode=FALSE)
  structure=list()

  idx=length(structure)+1
  structure[[idx]]=struct.G
  idx=idx+1
  args=list()
  args$mu=1/sqrt(V)
  args$H=H1
  args$a1=a1
  args$R1=R1
  args$name='level'
  structure[[idx]]=do.call(polynomial_block,args)

  if(all(!is.na(S))){
    idx=idx+1
    args=list()
    args$mu=S/sqrt(V)
    args$D=1
    args$name='S'
    # The variance must be scaled by the sample size.
    # The penalty is applied at the patient level, so we distribute it evenly
    # across patients so that, when the contributions are combined, the total
    # penalty equals 900.
    args$R1=900*N
    structure[[idx]]=do.call(polynomial_block,args)
  }

  # Here the scaling is done on Y, such that Y/sqrt(V) has variance 1.
  # This allows for time heteroscedasticity, which the kDGLM package does not support for Gaussians with known variance.
  args1[['outcome']]=Normal(mu='mu',V=1,data=Y/sqrt(V))
  structure=do.call(block_superpos,structure)

  args1$structure=structure

  return(do.call(fit_model,args1))
}


# Fits the model with all patients.
fit_full_normal=function(Y,S,N,T,struc.level,struct.G){
  labels=paste0('mu',formatC(1:N,flag = '0',width=3))
  args1=list(smooth=TRUE,safe.mode=FALSE)
  structure=list()

  idx.G=length(structure)+1
  structure[[idx.G]]=struct.G
  structure[[idx.G+1]]=struc.level
  args2=as.list(rep(1,N))
  names(args2)=labels
  args2$D=0.9
  structure[[idx.G+2]]=do.call(polynomial_block,args2)

  args1[['outcome']]=Normal(mu=labels,V=diag(N),data=Y)
  structure.marg=do.call(block_superpos,structure)

  args1$structure=structure.marg

  model.marg=do.call(fit_model,args1)

  args=S %>% apply(2,c,simplify='list')
  names(args)=labels
  args$D=1
  args$name='S'
  args$R1=900
  structure.cond=append(structure,
                        list(do.call(polynomial_block,args)))
  structure.cond=do.call(block_superpos,structure.cond)
  args1$structure=structure.cond
  model.cond=do.call(fit_model,args1)

  return(list(marg=model.marg,cond=model.cond))
}

# Initializes the random effects
init_level=function(N,D){
  structure=list()
  labels=paste0('mu',formatC(1:N,flag = '0',width=3))
  for(i in 1:N){
    args=list()
    args[[labels[i]]]=1
    args$D=D
    args$name=labels[i]
    structure=append(structure,
                     list(do.call(polynomial_block,args)))
  }

  structure=do.call(block_superpos,structure)
  return(structure)
}

# Initializes the treatment effects
init_G=function(labels,N,T,rho,R1){
  M=diag(T)*0
  diag(M[-1,-T])=1
  diag(M[-T,-1])=1
  D=diag(rowSums(M))
  R=solve((1-rho)*diag(T)+rho*(D-M))

  args=matrix(1,T,N)
  args=lapply(1:N,function(x){args[,x]})
  names(args)=labels

  args$period=T
  args$name='G'
  structure=do.call(ffs_block,args)
  structure$R1=R*R1
  return(structure)
}

# Compute the time homogeneity test.
calc_tests=function(delta.R.smp,delta.smp,pte.smp,
                    delta.R.est,delta.est,pte){
  # _____________________________ Wald _________________________________
  err=delta.R.smp-(1-pte.smp)*delta.smp

  th=delta.R.est-(1-pte)*delta.est
  var.th=var(err)

  A=eigen(var.th)
  A$values=c(A$values[-T],0)
  var.th=A$vectors%*%diag(A$values)%*%t(A$vectors)
  A$values=c(1/A$values[-T],0)
  inv.var.th=A$vectors%*%diag(A$values)%*%t(A$vectors)

  Wald=t(th)%*%inv.var.th%*%th
  Wald.crit=rowSums((err%*%inv.var.th)*th) %>% quantile(0.95)

  # _______________________________ MSD __________________________________
  A=rmvnorm(50000,rep(0,T),var.th)
  sd.th=sqrt(diag(var.th))
  sigma=matrix(sd.th,T,50000)
  MSD.crit=colMaxs(abs(A/sigma),value=TRUE) %>% quantile(0.95)
  MSD=max(abs(th/sd.th))

  return(list(MSD=MSD,MSD.crit=MSD.crit,Wald=Wald,Wald.crit=Wald.crit))
}



timestamp=function(){
  Sys.time() %>% format.Date('%Y%m%d_%H_%M_%OS6') %>%
    str_replace_all('\\.','')
}
