library(fdapace)

source('function/allfuns.R')
##############  parameters for our method
MM_FOD = function(data, iter, k_sub, alpha0, n, nout, p){

  tt = seq(0, 10, length.out = p)  ## grid for each curve
  yy = data  ## list
  #### first find clean set
  clean_hat = fun_h_clean(n,p,yy,iter,k_sub)$cleanset
  y_clean = list()
  for (i in 1:length(clean_hat)) 
  {
    y_clean$Lt[[i]] <- yy$Lt[[clean_hat[i]]]
    y_clean$Ly[[i]] <- yy$Ly[[clean_hat[i]]]
  }
  ####### result: our method
  out_hat = fun_test(n, p, yy, y_clean, nout,alpha0)$out_set
}
















