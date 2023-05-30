### This function performs funtion to function linear regression
### Input: lambda, gamma, y, D, initial b, tol1, tol2
### Return: optimal b 
FFR = function(tol1, tol2, lambda, gamma, y, D, b_0)
{
  n = ncol(b_0) #number of coordinates
  # Cycling through the coordinates
  cont = 1
  L_0 = 1000000
  diff = Inf
  while(diff > tol2)
  {
    # Parameter estimation for each coordinate
    for(j in 1:n)
    {
      b0 = theta0 = b_0[,j]
      sum = 0
      for(k in 1:n)
      {
        if(k!=j){sum = sum + D[,,k]%*%b_0[,k]}
      }
      r = y - sum
      aux_L = t(D[,,j])%*%D[,,j]+lambda*diag(dim(t(D[,,j])%*%D[,,j])[1])
      L = max(eigen(aux_L)$values)
      s=1/L
      normp = Inf
      t = 1
      while(normp > tol1)
      {
        z = theta0 +  s*(t(D[,,j])%*%(r-D[,,j]%*%theta0)-lambda*theta0)
        if(1-(s*gamma*sqrt(length(b0)))/(Norm(z,2))>0) 
        {b = (1-(s*gamma*sqrt(length(b0)))/(Norm(z,2)))*z} else {b = rep(0,length(b0))}
        theta = b +(t/(t+3))*(b-b0)
        normp = Norm(b-b0,2)
        b0 = b
        theta0 = theta
        t = t+1
        if(t>1000)
        {
          print(paste0("Max. number of iterations (t): ", ", parent ", j, ", iteration ", cont))
          break
        }
        s = s/t
      }
      b_0[,j] = b0
    }
    sum1 = rep(0,length(y))
    sum2 = 0
    sum3 = 0
    for(k in 1:n)
    {
      sum1 = sum1 + D[,,k]%*%b_0[,k]
      sum2 = sum2 + Norm(b_0[,k],2)^2
      sum3 = sum3 + sqrt(length(b_0[,k]))*Norm(b_0[,k],2)
    }
    L = (1/2)*Norm(y-sum1,2)^2+(lambda/2)*sum2+gamma*sum3
    diff = abs(L-L_0)/L_0
    L_0 = L
    cont = cont + 1
    if(cont>10000)
    {
      print(paste0("Max. number of iterations (cont): ", ", parent ", j, ", iteration ", cont))
      break
    }
  }
  return(b_0)
}


### This function performs cross-validation to tune the parameters of the FFR
### Input: lambda, gamma, Y, Z, P, tol1, tol2
### Return: optimal lambda, gamma, b 
cvFFR = function(lambda, gamma, Y, Z, P, tol1, tol2)
{
  n = nrow(Y)
  it = 5   
  nv = floor(n^0.7)  
  errcv = matrix(0,length(gamma),length(lambda))
  for(i in 1:it)
  {
    print(paste0("CV: ", i))
    # Training data
    s = sample.int(n,nv)
    Ytr = Y[-s,]
    Ztr = Z[-s,,,drop = FALSE]
    ytr = c(t(Ytr))
    Dtr = array(dim=c(dim(Ztr)[1]*dim(P)[2],dim(Ztr)[2]*dim(P)[1],dim(Ztr)[3]))
    for(j in 1:dim(Ztr)[3])
    {
      Dtr[,,j] = kronecker(Ztr[,,j],t(P))
    }
    # Testing data
    Yte = Y[s,]
    Zte = Z[s,,,drop = FALSE]
    yte = c(t(Yte))
    Dte = array(dim=c(dim(Zte)[1]*dim(P)[2],dim(Zte)[2]*dim(P)[1],dim(Zte)[3]))
    for(j in 1:dim(Zte)[3])
    {
      Dte[,,j] = kronecker(Zte[,,j],t(P))
    }
    # Parameters initialization
    b_0 = matrix(rnorm(dim(P)[1]^2*dim(Z)[3]),ncol = dim(Z)[3])
    # Testing error matrix
    err = matrix(Inf,length(gamma),length(lambda))
    for(j in 1:length(gamma))
    {
      for(k in 1:length(lambda))
      {
        # Function to function regression
        b_0 = FFR(tol1, tol2, lambda[k], gamma[j], ytr, Dtr, b_0)
        # Model retrieval
        d = which(b_0[1,]!=0)
        # Parameter estimation
        if(length(d)>=1)
        {
          Dtr.new = NULL
          for(l in 1:length(d))
          {
            Dtr.new = cbind(Dtr.new,Dtr[,,d[l]])
          }
          b = solve(t(Dtr.new)%*%Dtr.new+(lambda[k])*diag(dim(t(Dtr.new)%*%Dtr.new)[1]))%*%t(Dtr.new)%*%ytr
          b = matrix(b, ncol = length(d))
          b.new = b_0
          b.new[,d] = b
        } else
        {
          b.new = b_0
        }
        
        # Compute testing error
        sum1 = rep(0,length(yte))
        for(l in 1:dim(Zte)[3])
        {
          sum1 = sum1 + Dte[,,l]%*%b.new[,l]
        }
        err[j,k] = (1/nv)*Norm(yte-sum1,2)^2
      }
    }
    errcv = errcv + err
  }
  errcv = errcv/it
  opt = which(errcv == min(errcv), arr.ind = TRUE)
  lambda = lambda[opt[1,2]]
  gamma = gamma[opt[1,1]]
  y = c(t(Y))
  D = array(dim=c(dim(Z)[1]*dim(P)[2],dim(Z)[2]*dim(P)[1],dim(Z)[3]))
  for(j in 1:dim(Z)[3])
  {
    D[,,j] = kronecker(Z[,,j],t(P))
  }
  b_0 = FFR(tol1, tol2, lambda, gamma, y, D, b_0)
  return(list(b_0,c(lambda,gamma)))
}

### This function finds the optimal b for every node in the structure (train)
### and computes structure and prediction errors
### Input: nodes, arcs, data, lambda, gamma, tol1, tol2
### Return: structure learned, structure error, prediction error
NFDA = function(nodes, arcs, data, lambda, gamma, tol1, tol2)
{
  # true structure of the network
  structure = matrix(0, length(nodes), length(nodes))
  for(i in 1:dim(arcs)[1])
  {
    structure[arcs[i,1],arcs[i,2]] = 1
  }
  #Creates the coordinates for the x axis
  x=as.matrix(seq(0,1,length.out=dim(data)[2]))
  # Basis
  spb = 10
  splinebasis_B=create.bspline.basis(c(0,1),spb)
  base_B=eval.basis(as.vector(x),splinebasis_B)
  # Substract the mean
  B = array(dim=c(dim(data)[1],spb,length(nodes)))
  mean = matrix(0,dim(data)[2],dim(data)[3])
  for(i in 1:length(nodes))
  {
    smooth = smooth.basis(x,t(drop(data[,,i])),splinebasis_B)
    B[,,i] = t(smooth$fd$coefs)
    mean_B = colMeans(drop(B[,,i]))
    mean[,i] = base_B%*%mean_B
  }
  new_dat = array(dim=dim(data))
  for(i in 1:length(nodes))
  {
    new_dat[,,i] = data[,,i]-t(matrix(mean[,i],length(mean[,i]),dim(data)[1]))
  }
  data = new_dat
  # train and tets data sets
  train.ind = sample.int(dim(data)[1],(dim(data)[1])*0.8,replace = FALSE)
  train = data[train.ind,,]
  test = data[-train.ind,,]
  # recursion
  res = vector("list", length = length(nodes))
  opt_param = matrix(99,length(nodes),2)
  for(i in 2:length(nodes))
  {
    print(paste0("Node: ", i))
    Y = train[,,i]
    Z = array(dim=c(dim(train)[1],spb,i-1))
    for(j in 1:(i-1))
    {
      Z[,,j] = train[,,j]%*%base_B/dim(train)[2] 
    }
    P = t(base_B)
    ans = cvFFR(lambda, gamma, Y, Z, P, tol1, tol2)
    res[[i]] = ans[[1]]
    opt_param[i,] = ans[[2]]
  }
  # structure error
  str = matrix(0, length(nodes), length(nodes))
  for(i in 2:length(nodes))
  {
    for(j in 1:(i-1))
    {
      if(res[[i]][1,j]!=0)
      {
        str[j,i] = 1
      }
    }
  }
  diff = structure-str
  mat = matrix(99,2,2)
  mat[1,2] = sum(diff == 1)
  mat[2,1] = sum(diff == -1)
  mat[1,1] = sum(diff == 0 & str == 1)
  mat[2,2] = length(nodes)*(length(nodes)-1)/2-mat[1,2]-mat[2,1]-mat[1,1]
  # prediction error
  roots = NULL
  for(i in 1:length(nodes))
  {
    if(sum(str[,i])==0){roots = c(roots,i)}
  }
  childs = nodes[-roots]
  errpred = rep(0,length(nodes))
  for(i in 1:length(childs))
  {
    Y = test[,,childs[i]]
    Z = array(dim=c(dim(test)[1],spb,childs[i]-1))
    for(j in 1:(childs[i]-1))
    {
      Z[,,j] = test[,,j]%*%base_B/dim(train)[2] 
    }
    P = t(base_B)
    y = c(t(Y))
    D = array(dim=c(dim(Z)[1]*dim(P)[2],dim(Z)[2]*dim(P)[1],dim(Z)[3]))
    for(j in 1:dim(Z)[3])
    {
      D[,,j] = kronecker(Z[,,j],t(P))
    }
    b_0 = res[[childs[i]]]
    # Compute testing error
    sum1 = rep(0,length(y))
    for(l in 1:dim(Z)[3])
    {
      sum1 = sum1 + D[,,l]%*%b_0[,l]
    }
    errpred[childs[i]] = sum((y-sum1)^2)/(dim(test)[1]*dim(test)[2])
  }
  return(list(str,mat,errpred,opt_param))
}
