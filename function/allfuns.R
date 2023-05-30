##################################################################################################################
##Yu(2012) and ren (2017) method
####

 #SFOD method(Yu et. al.,2012,tech)
  Stat_SFOD<-function(N,C){
    
    cbar=apply(C,2,mean)
    cpp=cov(C)*(N-1)/N
    ev_cpp=eigen(cpp)
    
    d=0
    per=0
    while(per<0.85){
      d=d+1 
      per=per+ev_cpp$values[d]/sum(ev_cpp$values)
    }
    stat=matrix(0,N,1)
    for (i in 1:N){
      temp1=0
      for (j in 1:d){
        temp2=(C[i,]-cbar)%*%ev_cpp$vectors[,j]
        temp2=temp2**2/ev_cpp$values[j]
        temp1=temp1+temp2
      }
      stat[i]=temp1
    }
    maxstat=max(stat);
    maxno=which.max(stat);
    list(maxstat=maxstat,maxno=maxno,d=d,stat=stat)
  }#SFOD statistics
  
  Und<-function(d,N,alpha){
    und=2*qGumbel(1-alpha)+2*log(N)+(d-2)*log(log(N))-2*log(gamma(d/2)) 
    und
  }#SFOD's cutoff
  
  Gnd<-function(d,N,alpha){
    ns=1000
    G=matrix(0,ns,1)
    for(i in 1:ns){
      dat=rnorm(N*d)
      dim(dat)=c(N,d)
      mdat=apply(dat,2,mean)
      sdat=apply(dat,1,function(x){sum((x-mdat)^2)})
      G[i]=max(sdat)
    }
    gnd=quantile(G,1-alpha)
    gnd
  }#SFOD's cutoff
  
  SFOD<-function(N,C,alpha){
    samNo=1:N
    sN=N;sC=C
    out=0
    for(k in 1:N){
      res=Stat_SFOD(sN,sC)
      if(sN<100){lnd=Gnd(res$d,sN,alpha)}
      if(sN>=100){lnd=Und(res$d,sN,alpha)}
      #lnd=Und(res$d,sN,alpha)
      if(res$maxstat<lnd) break;
      if(res$maxstat>=lnd){
        sN=sN-1;sC=sC[-res$maxno,]
        out=c(out,samNo[res$maxno])
        samNo=samNo[-res$maxno]
      }
      k=k+1    
    }  
    out=out[-1]
    out
  }#SFOD procedure
  
  

  
  
  
  ##MDP
  MDP<-function(m1,h,X,N,p){
    H0=matrix(0,2,m1)
    H0=apply(H0,2,function(X){sort(sample(N,2))})    #initial subset for MDP
    H0=t(H0) #dimension m1*k
    
    H0_LTS=matrix(0,m1,h )  #save the best MDP subset 
    detD_mdp=matrix(0,m1,1) #save the best MDP detD
    trR2_mdp=matrix(0,m1,1) #save the best  MDP trR2
    dis_mdp=matrix(0,m1,N)
    
    for(i in 1:m1){
      Y=X[H0[i,],]
      Ybar=apply(Y,2,mean)
      S1=cov(Y)
      D=diag(S1)
      detD=prod(D)
      dis=matrix(0,N,1)
      for (j in 1:N){
        temp2=as.matrix(X[j,]-Ybar)
        dis[j]=t(temp2/D)%*%temp2
      }
      nn=sort(order(dis)[1:h]) #
      crit=100 #
      
      k=1
      while(crit!=0 & k<10){
        Y=X[nn,]
        Ybar=apply(Y,2,mean)
        S2=cov(Y)#
        D1=diag(S2)
        
        detD=prod(D1)
        dis=matrix(0,N,1)
        for (j in 1:N){
          temp2=as.matrix(X[j,]-Ybar)
          dis[j]=t(temp2/D1)%*%temp2
        }
        nn2=sort(order(dis)[1:h])
        crit=sum(abs(nn2-nn))
        nn=nn2
        k=k+1
      }
      ER=cor(X[nn,])
      trR2=sum(diag(ER%*%ER))-p^2/h
      
      H0_LTS[i,]=nn
      detD_mdp[i]=detD
      trR2_mdp[i]=trR2
      dis_mdp[i,]=dis
      
    }
    
    loc_mdp=which.min(detD_mdp)
    list(Hmdp=H0_LTS[loc_mdp,],dis=dis_mdp[loc_mdp,],trR2=trR2_mdp[loc_mdp])
    
  }#####end function MDP  
  
 #find H_LTFS
  LTFS<-function(mdp,m1,m2,C,N,h){
    cpp=cov(C[mdp,])   
    while(det(cpp)==0){
      a=setdiff(1:N,mdp)[1]
      mdp=sort(c(mdp,a))
      cpp=cov(C[mdp,])
    }
    ev_cpp=eigen(cpp)
    
    d=0
    per=0
    while(per<0.90){
      d=d+1
      per=per+ev_cpp$values[d]/sum(ev_cpp$values)
    }
    
    H0_LTS=matrix(0,2,m1)
    H0_LTS=apply(H0_LTS,2,function(X){sort(sample(N,2))})    #initial subset for MDP
    H0_LTS=t(H0_LTS) #dimension m1*k
    
    H2_LTS=matrix(0,m1,h)
    sumDH2=matrix(0,m1,1)
    for(i in 1:m1){
      H1_LTS=Cstep(H0_LTS[i,],h,C,N,ev_cpp,d)$H1
      temp=Cstep(H1_LTS,h,C,N,ev_cpp,d)
      H2_LTS[i,]=temp$H1
      sumDH2[i]=temp$sumDH 
    }
    
    ord_H2=order(sumDH2)[1:m2]
    H3=H2_LTS[ord_H2,]
    sumDH3=sumDH2[ord_H2]
    
    for(i in 1:m2){
      sum_dh0=1000000
      sum_dh1=sumDH3[i]
      h1=H3[i,]
      
      dk=1
      while((sum_dh0-sum_dh1)>0 & dk<=10){
        sum_dh1=sum_dh0
        temp=Cstep(h1,h,C,N,ev_cpp,d)
        h1=temp$H1
        sum_dh0=temp$sumDH
        dk=dk+1
      }
      H3[i,]=h1
      sumDH3[i]=sum_dh1
    }
    tt=which.min(sumDH3)
    Hltfs=H3[tt,]
    
    cpp=cov(C[Hltfs,])
    ev_cpp=eigen(cpp)
    d=0
    per=0
    while(per<0.90){
      d=d+1
      per=per+ev_cpp$values[d]/sum(ev_cpp$values)
    }
    Cbar=apply(C[Hltfs,],2,mean)
    DH=matrix(0,N,1)
    for (i in 1:N){
      temp1=0
      for (j in 1:d){
        temp2=(C[i,]-Cbar)%*%ev_cpp$vectors[,j]
        temp2=temp2**2/ev_cpp$values[j]
        temp1=temp1+temp2
      }
      DH[i]=temp1
    }
    
    list(Hltfs=Hltfs,DH=DH,dltfs=d)
    
  }#LTFS function  
  ##Cstep function 
  Cstep<-function(H0_LTS,h,X,N,ev_cpp,d){
    C=X[H0_LTS,]
    Cbar=apply(C,2,mean)
    DH=matrix(0,N,1)
    for (i in 1:N){
      temp1=0
      for (j in 1:d){
        temp2=(X[i,]-Cbar)%*%ev_cpp$vectors[,j]
        temp2=temp2**2/ev_cpp$values[j]
        temp1=temp1+temp2
      }
      DH[i]=temp1
    }
    H1=order(DH)[1:h]   
    sumDH=sum(DH[H1])
    list(H1=H1,sumDH=sumDH,DH=DH,d=d)
  }##end function Cstep
    ##Refined LTFS step
  Refined<-function(N,od,allC){
    nod=length(od)
    C=allC[-od,]
    if(nod==0){C=allC}
    newN=N-nod
    
    cbar=apply(C,2,mean)
    cpp=cov(C)
    ev_cpp=eigen(cpp)
    
    d=0
    per=0
    while(per<0.90){
      d=d+1
      per=per+ev_cpp$values[d]/sum(ev_cpp$values)
    }
    
    stat=matrix(0,N,1)
    for (i in 1:N){
      temp1=0
      for (j in 1:d){
        temp2=(allC[i,]-cbar)%*%ev_cpp$vectors[,j]
        temp2=temp2**2/ev_cpp$values[j]
        temp1=temp1+temp2
      }
      stat[i]=temp1
    }
    list(stat=stat,d=d)
  }
 ###############################################################################################################################################



###
####our method
###
#################################################################################
  fun_statistic=function(n,p,y,y_sample) 
  {
    T_stat <- rep(0,n)
    R <- FPCA(y_sample$Ly,y_sample$Lt,list(nRegGrid=50))   
    for (i in 1:n) 
    {
      
      ##sigma Y_ihat
      l <- y$Lt[[i]]
      #cat(l,"\n")
      location <- l
      #cat(location,"\n")
      G_matrix <- R$fittedCov[location,location]
      n_i <- length(l)
      var <- R$sigma2
      sigma_hat <- G_matrix + var*diag(n_i)
      
      ##H matrix and sigma_i 
      phii <- R$phi[location,]
      hhat <- diag(R$lambda)%*%t(phii)
      sigmai <- hhat%*%ginv(sigma_hat)%*%t(hhat)
      
      ##calculate T_i
      D <- dim(phii)[2]
      xi=rep(0,D)
      for (K in 1:D) {
        lambda <- R$lambda[K]
        phiik <- phii[,K]
        yi <- y$Ly[[i]]
        mui <- R$mu[location]
        xiik <- lambda*t(phiik)%*%ginv(sigma_hat)%*%(yi - mui)
        xi[K] <- xiik
      }
      T_stat[i] <- t(xi) %*% ginv(sigmai) %*% (xi)
    }
    #cat(T_stat,"\n")
    list(T_stat=T_stat,D=D)
  }
  
  
  #老的，已经废弃，用上面那个。
  fun_statistic_old=function(n,p,y,y_sample) 
{

T_stat <- rep(0,n)
R <- FPCA(y_sample$Ly,y_sample$Lt,list(nRegGrid=p,FVEthreshold=0.9,useBinnedData="OFF"))   
mu <- R$mu
lambda <- R$lambda
phi <- R$phi 
D <- dim(phi)[2]
xi=matrix(0,n,D)
for (i in 1:n) 
{   
      for (k in 1:D) {   
        xi[i,k]<- t(matrix((y$Ly[[i]]-mu)))%*%matrix(phi[,k])}
      T_stat[i] <- sum(xi[i,]^2/lambda) 
}
  
  list(T_stat=T_stat,D=D)
}
#################################################################################
#################################################################################
  fun_h_clean <- function(n,p,y,iter,k_sub)
  {
    T_matrix <- matrix(0,iter,n)
    for (m in 1:iter) 
    {
      ind <- sample(c(1:n),floor(n/2), replace = FALSE)
      y_sample <- list()
      for (k in 1:length(ind)) {
        y_sample$Lt[[k]] <- y$Lt[[ind[k]]]
        y_sample$Ly[[k]] <- y$Ly[[ind[k]]]
      }
      T_matrix[m,]=fun_statistic(n,p,y,y_sample)$T_stat 
    }
    T_max <- apply(T_matrix, 2, max)
    T_min <- apply(T_matrix, 2, min)
    clean_max=order(T_max,decreasing=FALSE)[1:floor(n*k_sub)]
    clean_min=order(T_min,decreasing=FALSE)[1:floor(n*k_sub)]
    cleanset <- intersect(clean_max,clean_min)
    list(cleanset=cleanset)
  }

  #老的，已经废弃，用上面那个。
fun_h_clean_old <- function(n,p,y,iter,k_sub)
{

T_matrix <- matrix(0,iter,n)
for (m in 1:iter) 
{
      ind <- sample(c(1:n),floor(n/2), replace = FALSE)
      y_sample <- list()
      for (k in 1:length(ind)) {
        y_sample$Lt[[k]] <- y$Lt[[ind[k]]]
        y_sample$Ly[[k]] <- y$Ly[[ind[k]]]
      }
   T_matrix[m,]=fun_statistic(n,p,y,y_sample)$T_stat 
}
  T_max <- apply(T_matrix, 2, max)
  T_min <- apply(T_matrix, 2, min)
  clean_max=order(T_max,decreasing=FALSE)[1:floor(n*k_sub)]
  clean_min=order(T_min,decreasing=FALSE)[1:floor(n*k_sub)]
  cleanset <- intersect(clean_max,clean_min)
list(cleanset=cleanset)
}


######################################################################################
fun_test=function(n,p,y,y_clean, n_out,alpha0)
{  
  T.test=rep(0,n)
  Q_T=rep(0,n)  
  
  AA=fun_statistic(n,p,y,y_clean)
  T.test=AA$T_stat
  D=AA$D
  th_median=median(T.test)/qchisq(0.5,D)
  T.test=T.test/th_median
  
  for (i in 1:n)
  {Q_T[i]=qnorm((pchisq(T.test[i],D)+1)/2)
  #cat(Q_T[i],"\n")
  }
  
  ############ next find the threshold t0 by iteration 
  len.tt=100
  U_thr=sqrt(2*log(n))
  t_it=seq(0,U_thr,length=len.tt)
  t_fdr=rep(0,len.tt)
  for (it in 1:len.tt)
  {
    t_fdr[it]=2*n*(1-pnorm(t_it[it]))/max(length(which(abs(Q_T)>=t_it[it])),1)
  }
  # print(T.test)
  # print(Q_T)
  
  if (length(which(t_fdr<=alpha0))>0){
    t0=t_it[min(which((t_fdr-alpha0)<=0))]}
  else
  {t0=sqrt(2*log(n))}
  # print(t_fdr)
  #cat("t0: ",t0,"\n")  
  ############### next perform outlier detection
  out_set=which(Q_T>t0)
  tp <- length(which(out_set <= n_out))/(n_out)
  fp <- length(which(out_set > n_out))/(n-n_out)
  
  ###### next for estimation
  clean_curve=setdiff(c(1:n),out_set)
  y_clean_curve <- list()
  for (i in 1:length(clean_curve)) 
  {
    y_clean_curve$Lt[[i]] <- y$Lt[[clean_curve[i]]]
    y_clean_curve$Ly[[i]] <- y$Ly[[clean_curve[i]]]
  }
  
  res_curve=FPCA(y_clean_curve$Ly,y_clean_curve$Lt,list(nRegGrid=50))
  workGrid_curve=res_curve$workGrid
  mu_curve=res_curve$mu
  cov_curve=res_curve$fittedCov
  
  list(out_set=out_set,pow=tp,size=fp,mu_curve=mu_curve,cov_curve=cov_curve,workGrid_curve=workGrid_curve)
}

#老的，已经废弃，用上面那个。
fun_test_old=function(n,p,y,y_clean, alpha0)
{  
  T.test=rep(0,n)
  Q_T=rep(0,n)  
  
  
  AA=fun_statistic(n,p,y,y_clean)
  T.test=AA$T_stat
  D=AA$D
  th_median=median(T.test)/qchisq(0.5,D)  
  
  # BB=fun_statistic_refined(y,y_clean,th_median) 
    # T.test=BB$T_stat
  # D=BB$D
  T.test=T.test/th_median
  for (i in 1:n)
  {Q_T[i]=qnorm((pchisq(T.test[i],D)+1)/2)
    #cat(Q_T[i],"\n")
  }
  
  ############ next find the threshold t0 by iteration 
  len.tt=100
  U_thr=sqrt(2*log(n))
  t_it=seq(0,U_thr,length=len.tt)
  t_fdr=rep(0,len.tt)
  for (it in 1:len.tt)
  {
    t_fdr[it]=2*n*(1-pnorm(t_it[it]))/max(length(which(abs(Q_T)>=t_it[it])),1)
  }
  #print(T.test)
  #print(Q_T)
  
  if (length(which(t_fdr<=alpha0))>0){
    t0=t_it[min(which((t_fdr-alpha0)<=0))]}
  else
  {t0=sqrt(2*log(n))}
  #print(t_fdr)
  #print(t0)  
  ############### next perform outlier detection
  #cat("t0: ",t0,"\n")
  out_set=which(Q_T>t0)
  list(out_set=out_set)
}
  
 
 
