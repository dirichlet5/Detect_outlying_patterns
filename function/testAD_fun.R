###########################################################################################################
# Gina-Maria Pomann
# 9/11/2012
#This script will test along each dimension of of the input data sets to see if they are different
#using the Anderson Darling Test
###########################################################################################################

# packageurl = "https://cran.r-project.org/src/contrib/Archive/adk/adk_1.0-2.tar.gz"
# install.packages(packageurl, repos=NULL, type="source")
library(adk)

testAD<-function(datX,datY,p){
  #datX is the matrix of coefficients for the X data  (number of subjects x number of PC scores)
  #datY is the matrix of coefficients for the Y data (number of subjects x number of PC scores)
  #p the number of coefficient vectors to test using the AD test
  #this gets overwritten, but you can easily comment out those lines if you 
  #want to specify p 
  
  #no transformation, just test 
  U<- datX
  V<- datY
  if(is.vector(U)==T || is.vector(V)==T){p=1}
  if(is.vector(U)==F && is.vector(V)==F){p=min(dim(U)[2],dim(V)[2])}
  #p=min(dim(U)[2],dim(V)[2])	
  pvalvec<-rep(NA,p)
  
  
  if(p==1){k=1;pvalvec[k]<-adk.test(U,V)$adk[1,2]}
  if(p>1){
    for(k in 1:p){
    #get the p value from the anderson darling test
      pvalvec[k]<-adk.test(U[,k],V[,k])$adk[1,2] 
    }
  }
  #now return the pvalue vector
  list(pval=pvalvec)
}
