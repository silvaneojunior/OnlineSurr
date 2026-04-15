source('reproducible_paper_codes/SantosJr and Parast (2026)/simulation/init.R')

set.seed(seed)

# Computing the splines basis for the surrogate
for(surr_type in c('linear')){
  if(surr_type=='linear'){
    # Linear relationship.
    # For consistency, we keep the same format as the splines case.
    P=1
    K=1
    spline=t(S.redux)
    T.ref=T.full
    S.index.ref=S.index
    dim(spline)=c(1,dim(spline))
  }else if(surr_type=='skip'){
    # Still linear, but we skip times were the outcome was not measured
    # The idea is just to save on computational power but removing redundant variables
    # In practice, we observed that this exclusion does affect the results, so we adivise against it.
    # We keep the code here just for the sake of transparency.
    P=1
    K=1
    spline=t(S.redux[seq.int(1,T.full,4),])
    T.ref=dim(spline)[2]
    S.index.ref=1:T.ref
    dim(spline)=c(1,dim(spline))
    }else{
      # B-splines, but with some possible variations latter on
    P=3
    T.ref=T.full
    S.index.ref=S.index
    K=10+P
    if(surr_type=='splines'){
      # Notes distributed uniformly within the observed surrogate range.
      t=seq.int(min(S.redux,na.rm=TRUE),max(S.redux,na.rm=TRUE),l=K-P+1)
    }else{
      # Notes distributed based on the quantiles of the data.
      t=seq.int(0,1,l=K-P+1)
      t=quantile(S.redux,t,na.rm=TRUE)
    }
    t=c(rep(min(t),P),t,rep(max(t),P))
    spline=array(S.redux,c(T.full,N,K+P))
    for(k in 1:(K+P)){
      spline[,,k]=(S.redux>=t[k]) & ((S.redux<t[k+1]) | ((k==(K+P)) & (S.redux==t[k+1])))
    }
    for(p in 1:P){
      for(k in 1:(K+P-p)){
        spline[,,k]=
          spline[,,k  ]*(if(t[k+p  ]==t[k  ]){0}else{(S.redux   -t[k] )/(t[k+p  ]-t[k  ])})+
          spline[,,k+1]*(if(t[k+p+1]==t[k+1]){0}else{(t[k+p+1]-S.redux)/(t[k+p+1]-t[k+1])})
      }
      spline=spline[,,-(K+P-p+1),drop=FALSE]
    }
    spline=aperm(spline,c(3,2,1))
  }
  spline[is.na(spline)]=0
  #####################################

  data.plot=data.frame()
  data.smp=list()
  boots.idx.mat=matrix(sample.int(N,N*N.boots,replace=TRUE),N,N.boots)

  # Computing the treatment for several possible maximum surrogate lags
  # -1 means no surrogate values (marginal model)
  # 0 means only the current surrogate value (conditional model, but with no lags)
  # k>0 includes all lags up to k.
  # Save the treatment effect at each time and compute the PTE separately later.
  for(lag.max in -1:T.ref){
    cat(paste0('lag.max=',lag.max,' - start!                                                                  \r'))

    # Random effects
    structure.base=
      (N*polynomial_block(mu=1,D=0.9,name='local.level'))

    # Labels for the kDGLM functions
    labels=structure.base$pred.names

    # Global level; that is, the trend common to all patients
    input=matrix(1,T,N) %>% apply(2,c,simplify=FALSE)
    names(input)=labels
    structure.base=structure.base+
      do.call(ffs_block,append(input,list(period=T,name='global.level')))

    # Sex effect
    input=apply(sex.base,2,c,simplify=FALSE)
    names(input)=labels
    structure.base=structure.base+
      do.call(polynomial_block,append(input,list(name='Sex')))

    # Smoker status effect
    input=apply(smoker.base,2,c,simplify=FALSE)
    names(input)=labels
    structure.base=structure.base+
      do.call(polynomial_block,append(input,list(name='Smoker')))

    # Treatment effect
    input=apply(treat.base,2,c,simplify=FALSE)
    names(input)=labels
    structure.base=structure.base+
      do.call(ffs_block,append(input,list(period=T,name='Treat')))

    # Defining the prior/regularization
    structure.base$R1=diag(structure.base$n)*9

    # If we're fitting the conditional model
    if(lag.max>=0){
      # Including the current surrogate value
      input=lapply(1:N,function(x){rep(1,T)})
      names(input)=labels
      structure.S=do.call(polynomial_block,append(input,list(order=K,D=1,name='spline')))
      structure.S$FF[,,]=spline[,,S.index.ref]
      structure.S$G[,,]=diag(K)
      structure.S$R1=diag(structure.S$n)*9

      # If we're also including surrogate lags
      if(lag.max>0){
        for(lag in 1:lag.max){
          input=lapply(1:N,function(x){rep(1,T)})
          names(input)=labels
          structure.S2=do.call(polynomial_block,append(input,list(order=K,D=1,name='spline2')))
          structure.S2$FF[,,which(S.index.ref>lag)]=spline[,,S.index.ref[S.index.ref>lag]-lag]
          structure.S2$FF[,,which(S.index.ref<=lag)]=0
          structure.S2$G[,,]=diag(K)
          structure.S2$R1=diag(structure.S2$n)*9
          structure.S=structure.S+structure.S2
        }
      }
      structure=structure.base+structure.S
    }else{
      structure=structure.base
    }

    cat(paste0('lag.max=',lag.max,' - Fitting: start!                                                                \r'))
    # Fitting the full model to get the evolution covariance
    model=fit_model(structure,outcome=Normal(mu=labels,V=diag(N),data=Y))
    cat(paste0('lag.max=',lag.max,' - Fitting: done!                                                             \r'))
    W.est=model$W

    # Saving the point estimate. Note that some reordering is required.
    delta.est=model$mt[(N+T+1+1+1:T),T][c(2:(T),1)]

    # Defining the individual models used to compute the contribution of each patient for the bootstrap
    structure.base=
      polynomial_block(mu=rep(1,T),name='local.level')+
      ffs_block(mu=1,period=T,name='global.level')+
      polynomial_block(mu=0,name='Sex')+
      polynomial_block(mu=0,name='Smoker')+
      ffs_block(mu=0,period=T,name='Treat')

    # The penalty must be divided among the patient contributions, so we multiply the prior variance by N.
    structure.base$R1=diag(structure.base$n)*9*N
    structure.base$R1[1,1]=9

    # If we're fitting the conditional model
    if(lag.max>=0){
      structure.S=polynomial_block(mu=1,order=K*(lag.max+1),D=1,name='splines')
      structure.S$FF[,,]=0
      structure.S$G[,,]=diag(K*(lag.max+1))
      structure.S$R1=diag(structure.S$n)*9*N

      structure=structure.base+structure.S
    }else{
      structure=structure.base
    }


    mt=matrix(NA,structure$n,N)
    St=array(NA,c(structure$n,structure$n,N))

    start=Sys.time()
    for(i in 1:N){
      args=list()
      args[['Outcome']]=Normal(mu='mu',V=1,data=Y[,i])
      args$smooth=FALSE
      args$safe.mode=FALSE

      # Replacing the covariantes
      structure.cur=structure
      structure.cur$FF[structure.cur$var.names$Sex,1,]=sex.base[,i]
      structure.cur$FF[structure.cur$var.names$Smoker,1,]=smoker.base[,i]

      init=structure.cur$var.names$Treat[1]
      structure.cur$FF[init,1,]=treat.base[,i]

      idx=-((1:N)[-i])
      structure.cur$H[,,]=W.est[idx,idx,]

      # Replacing the surrogate values
      if(lag.max>=0){
        structure.cur$FF[structure.base$n+1:K,1,]=spline[,i,S.index.ref]
        if(lag.max>0){
          for(lag in 1:lag.max){
            structure.cur$FF[structure.base$n+1:K+K*lag,1,which(S.index.ref>lag)]=spline[,i,S.index.ref[S.index.ref>lag]-lag]
            structure.cur$FF[structure.base$n+1:K+K*lag,1,which(S.index.ref<=lag)]=0
          }
        }
      }

      # Computing the contribution for patient i
      args$structure=structure.cur
      model.cur=do.call(fit_model,args)

      St[,,i]=S.cur=model.cur$Ct[,,T] %>% ginv
      mt[,i]=S.cur%*%model.cur$mt[,T]

      spent=difftime(Sys.time(),start,units = 'mins')
      ETA=(N-i)*spent/i
      cat(paste0('lag.max=',lag.max,' - pre-processing: ',round(100*i/N,2),'% - ETA: ',round(ETA,2),'mins                   \r'))
    }

    # Now that the patient contributions were computed, we can run the bootstrap
    N.boots=2000
    mt.smp=matrix(NA,structure$n-1,N.boots)
    rownames(mt.smp)=row.names(model.cur$mt)[-1]

    start=Sys.time()
    for(i in 1:N.boots){
      # The bootstrap indexes are computed seperately to make sure all models use the same indexes
      # This is important because we are only computing the marginal model once
      # But to account for the marginal and conditional models we must use the same indexes for both.
      boots.idx=boots.idx.mat[,i]

      mt.remix=mt[,boots.idx,drop=FALSE]
      St.remix=St[,,boots.idx,drop=FALSE]

      vals.diag=St.remix[1,1,] %>% ifelse(.==0,1,.)
      A_11=rowSums(St.remix,dims=2)[-1,-1]
      A_12=St.remix[1,-1,]
      S_11=A_11-tcrossprod(sweep(A_12,2,vals.diag,FUN='/'),A_12)

      a_1=rowSums(mt.remix[-1,])
      a_2=mt.remix[1,]/vals.diag

      mt.smp[,i]=solve(S_11,a_1-A_12%*%a_2)

      spent=difftime(Sys.time(),start,units = 'mins')
      ETA=(N.boots-i)*spent/i
      cat(paste0('lag.max=',lag.max,' - bootstraping: ',round(100*i/N.boots,2),'% - ETA: ',round(ETA,2),'mins                   \r'))
    }

    data.smp=append(data.smp,
                    list(list(lag.max=lag.max,
                              delta.est=delta.est,
                              mt.smp=mt.smp)))
    cat(paste0('lag.max=',lag.max,' - finshed!                                                                  \n'))
  }
  save.image(local_path(paste0(surr_type,'_N=',N,'.bkp')))
}
