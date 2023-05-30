### Required Packages
library(bnlearn)
library(fda)
library(MASS)
library(pracma)
library(abind)
library(keep)
library(fda)
library(bnlearn)
library(graph)
library(igraph)

### Function: Learns the parameters for a known structure
### Inputs: structure of the graph, curves, number of simulations
### Outputs: mean square prediction errors: dim1=method(1-FDGM, 2-FPCR) dim2=non-root nodes dim3=# of simulations
results = function(d, dat, sim)
{
  simulation = dat
  #Param. definition 
  spb=10                          #number of functions for the b-spline basis
  m=dim(simulation)[1]/sim        #number of replicates of functions
  n=dim(simulation)[2]            #defines number of time obs.
  # Graph structure
  nodes = strtoi(nodes(d))
  arcs = matrix(strtoi(arcs(d)), ncol = 2)
  structure = matrix(0, length(nodes), length(nodes))
  for(i in 1:dim(arcs)[1])
  {
    structure[arcs[i,1],arcs[i,2]] = 1
  }
  roots = nodes[!(nodes %in% arcs[,2])]
  childs = nodes[(nodes %in% arcs[,2])]
  #Creates the coordinates for the x axis
  x=as.matrix(seq(0,1,length.out=n))
  
  ##Iterations for each simulation##
  results_test = array(0,dim = c(2,length(childs),sim))
  for(it in 1:sim)
  {
    data = simulation[((it-1)*m+1):(it*m),,]
    train.ind = sample.int(m,m*0.8,replace = FALSE) #80% of the data is used for training
    train = data[train.ind,,]
    test = data[-train.ind,,]
    
    ###############Basis estimation#############
    splinebasis_B = create.bspline.basis(c(0,1),spb)
    base_B = eval.basis(as.vector(x),splinebasis_B)
    B = array(dim=c(length(train.ind),spb,length(nodes)))
    pc = vector("list", length = length(nodes))
    nharm = rep(0,length(nodes))
    for(i in 1:length(nodes))
    {
      smooth = smooth.basis(x,t(drop(train[,,i])),splinebasis_B)
      B[,,i] = t(smooth$fd$coefs)
      obj = smooth$fd
      pca = pca.fd(obj, spb)
      var.pca = cumsum(pca$varprop)
      nharm[i] = sum(var.pca < 0.95) + 1
      pc[[i]] = pca.fd(obj, nharm[i])
    }
    
    ###############FDGM esimates####################
    beta = vector("list", length = length(nodes))
    J=t(base_B)%*%base_B/n
    for(i in 1:length(childs))
    {
      Z = NULL
      par = which(structure[,childs[i]]==1)
      beta[[childs[i]]] = vector("list", length = length(par))
      for(j in 1:length(par))
      {
        Z = cbind(Z, (train[,,par[j]]%*%base_B)*(1/n)) 
      }
      ZZ=t(Z)%*%Z
      k_JZ=kronecker(J,ZZ)
      B_h = solve(k_JZ)%*%c(t(Z)%*%drop(B[,,childs[i]])%*%J)
      B_h=matrix(B_h,ncol=spb)
      B_h=array(B_h, dim = c(spb,length(par),spb))
      for(j in 1:length(par))
      {
        beta_h = (base_B%*%B_h[,j,]%*%t(base_B))
        beta[[childs[i]]][[j]] = beta_h
      }
    }
    
    ###############FPCR esimates####################
    coef_pc = vector("list", length = length(nodes))
    for(i in 1:length(childs))
    {
      Xreg = 1
      par = which(structure[,childs[i]]==1)
      for(j in 1:length(par))
      {
        Xreg = cbind(Xreg, pc[[par[j]]]$scores) 
      }
      coef_pc[[childs[i]]] = t(solve(t(Xreg)%*%Xreg)%*%t(Xreg)%*%pc[[childs[i]]]$scores)
    } 
  
    #Prediction error - Test data set 
    pred = m-length(train.ind)
    err_te = matrix(0,2,length(childs))
    for(i in 1:pred)
    {
      x_bn = array(0, dim = c(n,length(nodes)))
      x_pc = list()
      for(j in 1:length(nodes))
      {
        x_pc[[j]] = rep(0, nharm[j])
      }
      xx_pc = array(0, dim = c(n,length(nodes)))
      for(j in 1:length(roots))
      {
        x_bn[,roots[j]] = drop(test[i,,roots[j]])
        aux_pc = drop(test[i,,roots[j]])-base_B%*%pc[[roots[j]]]$meanfd$coefs
        smooth_te=smooth.basis(x,aux_pc,splinebasis_B)
        B_te=t(smooth_te$fd$coefs)
        x_pc[[roots[j]]] = (B_te%*%(pc[[roots[j]]]$harmonics$coefs))*(1/10)
        xx_pc[,roots[j]] = drop(test[i,,roots[j]])
       }
      for(j in 1:length(childs))
      {
        par = which(structure[,childs[j]]==1)
        x_aux_pc = 1
        for(k in 1:length(par))
        {
          x_bn[,childs[j]] = x_bn[,childs[j]] + drop(drop(x_bn[,par[k]])%*%drop(beta[[childs[j]]][[k]])/n)
          x_aux_pc = c(x_aux_pc,drop(x_pc[[par[k]]]))
        }
        x_pc[[childs[j]]] = drop(drop(coef_pc[[childs[j]]])%*%x_aux_pc)
        xx_pc[,childs[j]] = base_B%*%pc[[childs[j]]]$harmonics$coefs%*%x_pc[[childs[j]]] + base_B%*%pc[[childs[j]]]$meanfd$coefs
        err_te[1,j] = err_te[1,j] + sum((drop(test[i,,childs[j]])-drop(x_bn[,childs[j]]))^2)/n
        err_te[2,j] = err_te[2,j] + sum((drop(test[i,,childs[j]])-drop(xx_pc[,childs[j]]))^2)/n
      }
    }
    err_te = err_te/pred
    ### Results
    results_test[,,it] = err_te
  }
  return(results_test)
}

### Parameter learning
set.seed(17031990)
print(Sys.time())
load("dat.Rdata")
res = results(d, dat, 2)
print(Sys.time())


