### 1. 
source("./function/Qu_Zhuo_JCGS/simulation-function.R")
source("./function/Qu_Zhuo_JCGS/simulation_parameter.R")
source("./function/Qu_Zhuo_JCGS/one_sample.R")
source("./function/Qu_Zhuo_JCGS/multi_sparse.R")



### 2.Define functions used to add the four types of outliers.

add_outlier_type4 <- function(normalData, p , n_outliers, outlier_rate = 0.05,
                              mu = 4, q = 15, a = 0.1, b = 0.9, kprob = 0.5,
                              cov_alpha = 1, cov_beta = 1, cov_nu = 1,
                              deterministic = TRUE, seed = NULL,
                              plot = F,
                              plot_title = "Simulation Model 3",
                              title_cex = 1.5,
                              show_legend = T,
                              ylabel = "",
                              xlabel = "gridpoints"){
  tt <- seq(0, 1, length.out = p)
  covfun <- covfunexp(tt, cov_alpha, cov_beta, cov_nu)
  muu <- mu*tt
  L <- chol(covfun)
  
  ### Generate Data
  if(!is.null(seed)){
    set.seed(seed)
  }
  e <- matrix(rnorm(n*p), nrow = p, ncol = n)
  #y <- muu+t(L)%*%e
  y <- t(normalData) #The original y is in p*n form.
  # set seed
  
  # check deterministic or probabilistic
  #dtt <- determine(deterministic, n, outlier_rate)
  #n_outliers <- dtt$n_outliers
  true_outliers <- sample(1:dim(normalData)[1],n_outliers,replace = F)
  
  # generate outliers
  if(n_outliers > 0){
    e <- matrix(rnorm(p*n_outliers), nrow = p)
    qcoeffk <- rbinom(n_outliers, 1, kprob)
    qcoeffk[qcoeffk == 0] <- -1
    qcoeffk <- qcoeffk*q
    indicator <- sapply(runif(n_outliers, a, b),
                        function(x) (tt >= x) )
    #y[, true_outliers] <- (muu+ t(L)%*%e) + (indicator*rep(qcoeffk, rep(p, n_outliers)))
    
    y[, true_outliers] <- t(normalData[true_outliers,]) + (indicator*rep(qcoeffk, rep(p, n_outliers)))
  }
  y[,-true_outliers] <- t(normalData[-true_outliers,])
  if(plot){
    plot_dtt(y, tt, p, true_outliers, show_legend, plot_title,
             title_cex, ylabel, xlabel )
  }
  return(list(data = t(y), true_outliers = true_outliers))
}


add_outlier_type3 <- function(normalData, p = 50, n_outliers, outlier_rate = 0.1,
                              mu = 4, q = 1.8, kprob = 0.5,
                              a = 0.25, b = 0.75,
                              cov_alpha = 1, cov_beta = 1, cov_nu = 1,
                              pi_coeff = 0.02, exp_pow = 2,
                              exp_coeff = 50,
                              deterministic = TRUE, seed = NULL,
                              plot = F,
                              plot_title = "Outlier Type 3",
                              title_cex = 1.5,
                              show_legend = T,
                              ylabel = "",
                              xlabel = "gridpoints"){
  tt <- seq(0, 1, length.out = p)
  covfun <- covfunexp(tt, cov_alpha, cov_beta, cov_nu)
  muu <- mu*tt
  L <- chol(covfun)
  
  # set seed
  if(!is.null(seed)){
    set.seed(seed)
  }
  e <- matrix(rnorm(n*p), nrow = p, ncol = n)
  y <- t(normalData) #The original y is in p*n form.
  
  true_outliers <- sample(1:dim(normalData)[1],n_outliers,replace = F)
  
  if(n_outliers > 0){
    e <- matrix(rnorm(p*n_outliers), nrow = p)
    u <- rbinom(n_outliers, 1, kprob)
    qcoeffk <- (-1)^u
    qcoeffk <- qcoeffk*q
    qcoeffk2 <- (-1)^(1-u)
    unif_values <- rep(runif(n_outliers, a, b), each = p)
    indicator <- exp(-exp_coeff*((tt-unif_values)^exp_pow))/sqrt(pi_coeff*pi) *
      rep(qcoeffk2, rep(p, n_outliers))
    y[, true_outliers] <- t(normalData[true_outliers,]) + indicator + rep(qcoeffk, rep(p, n_outliers))
  }
  y[,-true_outliers] <- t(normalData[-true_outliers,])
  
  if(plot){
    plot_dtt(y, tt, p, true_outliers, show_legend, plot_title,
             title_cex, ylabel, xlabel, legend_pos = "bottomright" )
  }
  return(list(data = t(y), true_outliers = true_outliers))
}

add_outlier_type2 <- function(normalData, p = 50, n_outliers, 
                              mu = 4,  sin_coeff = 9, pi_coeff = 4,
                              a = .25, b = .75,
                              cov_alpha = 1, cov_beta = 1, cov_nu = 1,
                              deterministic = TRUE, seed = NULL,
                              plot = F,
                              plot_title = "Simulation Model 7",
                              title_cex = 1.5,
                              show_legend = T,
                              ylabel = "",
                              xlabel = "gridpoints"){
  # ana model3
  tt <- seq(0, 1, length.out = p)
  covfun <- covfunexp(tt, cov_alpha, cov_beta, cov_nu)
  muu <- mu*tt
  L <- chol(covfun)
  
  ### Generate Data
  if(!is.null(seed)){
    set.seed(seed)
  }
  e <- matrix(rnorm(n*p), nrow = p, ncol = n)
  y <- t(normalData)
  
  # set seed
  
  
  # check deterministic or probabilistic
  #dtt <- determine(deterministic, n, outlier_rate)
  #n_outliers <- dtt$n_outliers
  true_outliers <- sample(1:dim(normalData)[1],n_outliers,replace = F)
  
  ## Generate outlier
  if(n_outliers){
    e <- matrix(rnorm(p*n_outliers), nrow = p)
    theta <- rep(runif(n_outliers, a, b), each = p)
    indicator <- sin_coeff * sin(pi_coeff*pi*(tt+ theta))
    y[, true_outliers] <- t(normalData[true_outliers,]) + indicator
  }
  y[,-true_outliers] <- t(normalData[-true_outliers,])
  if(plot){
    plot_dtt(y, tt, p, true_outliers, show_legend, plot_title,
             title_cex, ylabel, xlabel, legend_pos = "bottomright" )
  }
  
  return(list(data = t(y), true_outliers = true_outliers) )
}



add_outlier_type1 <- function(normalData, p = 50, n_outliers, outlier_rate = 0.1,
                              mu = 4, zero_len = 5, zero_times = 3,
                              cov_alpha = 1, cov_beta = 1, cov_nu = 1,
                              pi_coeff = 0.02, exp_pow = 2,
                              exp_coeff = 50,
                              deterministic = TRUE, seed = NULL,
                              plot = F,
                              plot_title = "Simulation Model 6",
                              title_cex = 1.5,
                              show_legend = T,
                              ylabel = "",
                              xlabel = "gridpoints"){
  tt <- seq(0, 1, length.out = p)
  covfun <- covfunexp(tt, cov_alpha, cov_beta, cov_nu)
  muu <- mu*tt
  L <- chol(covfun)
  
  # set seed
  if(!is.null(seed)){
    set.seed(seed)
  }
  e <- matrix(rnorm(n*p), nrow = p, ncol = n)
  y <- t(normalData) #The original y is in p*n form.
  
  true_outliers <- sample(1:dim(normalData)[1],n_outliers,replace = F)
  
  if(n_outliers > 0){
    for (index_out in true_outliers) {
      zero_start = sample(1:(p-zero_len), zero_times)
      zero_point = c()
      for (len in 1:zero_len-1) {
        zero_point = c(zero_point, zero_start + len)
      }
      
      y[zero_point, index_out] <- 0
    }
    
  }
  y[,-true_outliers] <- t(normalData[-true_outliers,])
  
  if(plot){
    plot_dtt(y, tt, p, true_outliers, show_legend, plot_title,
             title_cex, ylabel, xlabel, legend_pos = "bottomright" )
  }
  return(list(data = t(y), true_outliers = true_outliers))
}



### 3. Generate the data (with outliers), and plot graph.
argval = seq(0, 1, length.out = 50)
sim = generate_samples("no outlier", eigenvalue_type, argval) #Note that the outliers here are generated separately.
normalData <- sim$observed_curve[[3]] #Corresponding to Model 1

#Correspond to Model 3 in fdaoutlier
dt4 <- add_outlier_type4(normalData, p=50, q=8, n_outliers=5,  plot = TRUE,plot_title = "Outlier Type 4",xlabel ="Time")

#Correspond to Model 6 in fdaoutlier
dt3 <- add_outlier_type3(normalData, p=50, q=8, n_outliers=5,  plot = TRUE,plot_title = "Outlier Type 3",xlabel ="Time")

#Correspond to Model 7 in fdaoutlier
dt2 <- add_outlier_type2(normalData, p=50,  n_outliers=5, plot = TRUE,plot_title = "Outlier Type 2",xlabel ="Time")

#Correspond to Model 2 in fdaoutlier
dt1 <- add_outlier_type1(normalData, p=50, zero_len = 5, zero_times = 3, n_outliers=1, plot = TRUE,plot_title = "Outlier Type 1",xlabel ="Time")

# Now we convert the data into sparse data.
aa <- list(dt1$data,dt2$data,dt3$data)
bb <- list(dt1,dt2,dt3)
missingness <-0.65
MISS_PER <- 1-missingness
p<-length(aa)
sparse_data = multi_sparse(aa[[1]], p_size = rep(1, p),
                           p_curve = rep(MISS_PER, p),
                           sparsity = c("point"))
observed_curve = aa
TYPE = 2 #Specify the selected Outlier type here!
out_index = bb[[TYPE]]$true_outliers
plot_sparse(TYPE)
#




