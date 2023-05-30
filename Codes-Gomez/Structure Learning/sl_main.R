### Required Packages
suppressWarnings(library(bnlearn))
suppressWarnings(library(fda))
suppressWarnings(library(pracma))

source('sl_functions_mcv_int.R', echo=TRUE)

# Data for the first simulation
load("dat.Rdata")
data = dat[((1-1)*100+1):(1*100),,]

nodes = strtoi(nodes(d))
arcs = matrix(strtoi(arcs(d)), ncol = 2)

tol1 = 0.001
tol2 = 0.001

lambda = c(10^-10,10^-8,10^-6,10^-4,10^-2)
gamma = c(0,0.2,0.4,0.6,0.8,1,3,5,7,9,11,13,15)

ans = NFDA(nodes, arcs, data, lambda, gamma, tol1, tol2)


