### Required Packages
library(bnlearn)
library(fda)

### Function: Generation of the graph
### Inputs: number of nodes, proportion of arcs in the graph (compared with a fully connected graph)
### Returns: graph structure
dagraph = function(num.nodes, den)
{
  dag = random.graph(paste(1:num.nodes), num = 1, method = "ordered", prob = den, debug = FALSE)
  ### Parameter definition
  nodes = strtoi(nodes(dag))
  arcs = matrix(strtoi(arcs(dag)), ncol = 2)
  return(list(dag,nodes,arcs))
}

### Function: Data generation
### Inputs: graph (nodes, arcs), snr, number of simulations 
### Returns: simulated curves
datgen = function(nodes, arcs, snr, sim)
{
  # Graph structure
  structure = matrix(0, length(nodes), length(nodes))
  for(i in 1:dim(arcs)[1])
  {
    structure[arcs[i,1],arcs[i,2]] = 1
  }
  roots = nodes[!(nodes %in% arcs[,2])]
  childs = nodes[(nodes %in% arcs[,2])]
  #Param. definition
  n=50            #defines number of time obs.
  m=100           #number of replicates of functions
  #Generates all the info
  simnodes = array(dim=c((m*sim),n,length(nodes)))
  #Roots generation
  r = root(roots, snr, m, n, sim)
  #Childs generation
  for(i in 1:length(roots))
  {
    simnodes[,,roots[i]] = r[[i]]
  }
  for(i in 1:length(childs))
  {
    pa_val =  list()
    for(j in 1:sum(structure[,childs[i]]))
    {
      par = which(structure[,childs[i]]==1)
      pa_val[[j]] = drop(simnodes[,,par[j]])
    }
    simnodes[,,childs[i]] = node(pa_val, snr, sim, m, n)
  }
  return(simnodes)
}

### Function: Generates the data for the root nodes
### Inputs: root nodes, snr, # of curves to generate, # of time obs. to generate, # of simulations
### Outputs: root nodes' curves
root = function(roots, snr, m, n, sim)
{
  x_1 = list()
  for(i in 1:length(roots))
  {
    x = seq(0,1,length=n)
    E = as.matrix(dist(x, diag=T, upper=T))
    Sigma = exp(-10*E^2)
    eig = eigen(Sigma)
    Sigma.sqrt = eig$vec%*%diag(sqrt(eig$val+10^(-10)))%*%t(eig$vec)
    mean1 = Sigma.sqrt%*%rnorm(n)
    S_noise = exp(-0.1*E^2)
    eig_noise = eigen(S_noise)
    S.sqrt_noise = eig_noise$vec%*%diag(sqrt(eig_noise$val+10^(-10)))%*%t(eig_noise$vec)
    noise = S.sqrt_noise%*%rnorm(n)
    signal = mean1 + noise
    var = var(signal)
    ds1 = sqrt(var/snr)
    S.sqrt_err = diag(n)*drop(ds1)
    x1 = matrix(0,m*sim,n)
    for(j in 1:(m*sim))
    {
      noise = S.sqrt_noise%*%rnorm(n)
      error = S.sqrt_err%*%rnorm(n)
      x1[j,] = mean1 + noise + error
    }
    x_1[[i]] = x1 
  }
  return(x_1)
}

### Function: Generates the data for the non-root nodes
### Inputs: parents of the node, snr, # of curves to generate, # of time obs. to generate, # of simulations
### Outputs: curves generated for the non-root node
node = function(pa_val, snr, sim, m, n)
{
  beta = list()
  for(i in 1:length(pa_val))
  {
    x = seq(0,1,length=n)
    E = as.matrix(dist(x, diag=T, upper=T))
    Sigma = exp(-10*E^2)
    eig = eigen(Sigma)
    Sigma.sqrt = eig$vec%*%diag(sqrt(eig$val+10^(-10)))%*%t(eig$vec)
    g = matrix(0,n,3)
    p = matrix(0,n,3)
    for(j in 1:3)
    {
      g[,j] = Sigma.sqrt%*%rnorm(n)
      p[,j] = Sigma.sqrt%*%rnorm(n)
    }
    beta[[i]] = (g%*%t(p)) 
  }
  xn = matrix(0,m*sim,n)
  for(i in 1:(m*sim))
  {
    xaux = matrix(0,1,n)
    for(j in 1:length(pa_val))
    {
      xaux = xaux + (pa_val[[j]][i,]%*%t(beta[[j]]))/n
    }
    ds = sqrt(var(drop(xaux))/(snr)) 
    xn[i,] = xaux + rnorm(n,0,ds)
  }
  return(xn)
}


### Simulation scenarios
sim=2        #number of simulations
########### SNR ###############
set.seed(17031990)
num.nodes = 10
den = 0.2
snr = 200
dag = dagraph(num.nodes, den)
plot(dag[[1]])
set.seed(1502)
dat = datgen(dag[[2]], dag[[3]], snr, sim)
d=dag[[1]]
save(d, dat, file = "dat.Rdata")
