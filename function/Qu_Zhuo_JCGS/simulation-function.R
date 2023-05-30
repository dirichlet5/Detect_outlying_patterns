one_sample <- function(outlierst, randomnumber, eigenvalue_type, argval) {
  standard_mu <- c()
  standard_epsilon <- c()
  standard_u <- c()
  for (j in 1:p) {
    standard_mu <- c(standard_mu, mu(j, argval))
    standard_epsilon <- c(standard_epsilon, epsilon(j, argval))
    standard_u <- c(standard_u, u(j, argval, eigenvalue_type))
  }
  k <- rbinom(p, 1, 0.5)
  sign <- sapply(k, function(oo) {
    ifelse(oo == 1, 1, -1)})
  if (randomnumber == 1) {
    if (outlierst == "persistent magnitude outlier") {
      standard_u <- standard_u + 12 * rep(sign, each = length(argval))
    } else if (outlierst == "isolated magnitude outlier") {
      st <- sample(1:(0.9 * length(argval)), p, replace = TRUE)
      period <- rep(0, (p * length(argval)))
      for (i in 1:p) {
        start <- (i - 1) * length(argval) + st[i]
        period[start:(start + round(0.1 * length(argval)))] <- sign[i]
      }
      standard_u <- standard_u + 12 * period
    } else if (outlierst == "shape outlier I") {
      standard_mu <- c()
      shift_val <- c(-0.3, -0.2, -0.5)
      for (j in 1:p) {
        standard_mu <- c(standard_mu, mu(j, argval + shift_val[j]))
      } 
    } else if (outlierst == "mixed outlier") {
      standard_mu <- c()
      for (j in 1:p) {
        interruption <- rexp(1, 2)
        addvalue <- (2 + interruption) * mu(j, argval)
        if (j == p) {
          addvalue <- addvalue - 6
        }
        standard_mu <- c(standard_mu, addvalue)
      }
    } 
  } 
  if (outlierst == "shape outlier II") {
    if (randomnumber == 1) {
      shift_value <- c(2 * sin(4 * pi * argval), 
                       2 * cos(4 * pi * argval),
                       2 * cos(8 * pi * argval))
    } else {
      shift_value <- runif(3 * length(argval), -2.1, 2.1)
    }
    standard_u <- standard_u + shift_value
  } else if (outlierst == "joint outlier") {
    if (randomnumber == 1) {
      z <- runif(p, 2, 8)
    } else {
      rgen <- runif (1, 2, 8)
      #z <- c(rgen, (3 * log(rgen) + 2), (rgen^2 / 4))
      z <- c(rgen, 8 - rgen, rgen - 2)
    }
    addvalue <- c(z[1] * argval * sin(pi * argval), z[2] * argval * cos(pi * argval),
                  z[3] * argval * sin(2* pi * argval))
    standard_u <-  addvalue
    #standard_u <- standard_u + addvalue
  } else if (outlierst == "covariance outlier") {
    if (randomnumber == 1) {
      nu_candidate <- runif(p, 0.1, 0.2)
    } else {
      nu_candidate <- runif(p, 2, 3)
    }
    matern_covar <- c()
    for (i in 1:p) {
      row_mat <- c()
      for (j in 1:p) {
        row_mat <- cbind(row_mat,
                         whittle_block(i, j, nu_candidate, std, rho_cor))
      }
      matern_covar <- rbind(matern_covar, row_mat)
    }
    standard_epsilon <- mvrnorm(n = 1, mu = rep(0, p * length(argval)),
                                Sigma = matern_covar)
  } 
  true_func <- standard_mu + standard_u    
  observed_func <- standard_mu + standard_u + standard_epsilon
  
  return (list(true_func = true_func,
               observed_func = observed_func))
}




mu <- function(j, argval) {
  if (j == 1) {
    func <- 5 * sin(2 * pi * argval)+10
  } else if (j == 2) {
    func <- 5 * cos(2 * pi * argval)
  } else if (j == 3){
    func <- 5 * (argval - 1.4)^2+5
  } else {
    func <- 5 * argval ^ j * log (argval + 1)
  }
  return (func)
}




epsilon <- function(j, argval) {
  func <- rnorm(length(argval), 0, sd = std[j])
  return (func)}




eigenscore <- function(m, eigenvalue_type) {
  if (eigenvalue_type == "linear") {
    func <- rnorm(1, 0, sd = sqrt((M + m - 1) / M))
  } else {
    func <- rnorm(1, 0, sd = exp(-(m + 1)/2) )
  }
  return (func)
}




u <- function(j, argval, eigenvalue_type) {
  # if (method == "MCov") {
  #   func <- mvrnorm(1, rep(0, p * length(argval)),
  #                   Sigma = covar()) ## vector 
  # } else if (method == "MFPCA") {
  eig_fun <- eigenfun(j, argval) ### (50*3) * M matrix
  scorevector <- sapply(1:M, function(i) {
    eigenscore(i, eigenvalue_type)})
  func <- as.numeric(eig_fun %*% scorevector) ### vector
  # }
  return (func)
}

eigenfun <- function(j, argval) {
  whole_time <- unmatrix(sapply(1:p, function(l) {l - 1 + argval}))
  basis <- fourier(whole_time, M)
  k <- runif(1, -1, 1)
  if (k < 0) {
    fu <- basis[((j - 1) * length(argval) + 1):(j * length(argval)), ]
  } else { 
    fu <- (-1) * basis[((j - 1) * length(argval) + 1):(j * length(argval)), ]
  } 
  return (fu)
}


generate_samples <- function(outlierst, eigenvalue_type, argval) {
  rand_nb <- rbinom(n, 1, contamination_level)
  sim <- lapply(1:n, function (i) {
    sample_generate <- one_sample(outlierst, rand_nb[i], eigenvalue_type, argval)
    return(sample_generate)
  })
  mat_true <- matrix(unlist(lapply(sim, function(i) {i$true_func})),
                     nrow = length(argval) * p, ncol = n, byrow = FALSE)
  mat_true_list <- lapply(1:p, function(varr) {
    t(mat_true[((varr - 1) * length(argval) + 1):
                 (varr * length(argval)), ])
  })
  mat_observed <- matrix(unlist(lapply(sim, function(i) {i$observed_func})),
                         nrow = length(argval) * p, ncol = n, byrow = FALSE)
  mat_observed_list <- lapply(1:p, function(varr) {
    t(mat_observed[((varr - 1) * length(argval) + 1):
                     (varr * length(argval)), ])
  })
  out_index <- which(rand_nb == 1)
  if (outlierst == "no outlier") {
    out_index <- NULL
  }
  return (list(true_curve = mat_true_list, observed_curve = mat_observed_list, outlier = out_index))
}


multi_sparse <- function(mat, p_size, p_curve, sparsity) {
  ####### assumes p_size and  p_curve are vectors
  if (length(mat) == length(p_size) && 
      length(p_curve) == length(sparsity) && 
      length(p_size) == length(p_curve)) {
    func <- function(i) {
      y <- uni_sparse(mat[[i]], p_size[i], p_curve[i], sparsity[i])
      return (y)
    }
    myfunc <- cmpfun(func)
    res <- lapply(1:p, FUN = myfunc)
    return (res)
  } else {
    stop ("not valid input")
  }
  ## res is a list where each element contains sparse matrix n * length(sett)
}


uni_sparse <- function(dat, p_size, p_curve, sparsity) {
  ### dat is a matrix of n * length(sett) 
  sparsecase <- dat
  nc <- length(sett)
  n <- nrow(dat)
  if (sparsity != "dense" && p_size > 0 && p_curve > 0) {
    sparse_ind <- sample.int(n, ceiling(n * p_size))
  } else {
    sparse_ind <- NULL
  }
  if (length(sparse_ind) > 1) {
    if (sparsity != "partial") {
      vep <- apply(sparsecase[sparse_ind, ], 1, function(l) {
        sp_fun <- observed_measure(l, sparsity, p_curve)
        sp_index <- sp_fun$sparse_index
        sp_data <- sp_fun$value
        return(sp_data)
      })
      sparsecase[sparse_ind, ] <- t(vep)
    } else {
      sp_fun <- observed_measure(sparsecase[sparse_ind[1], ], sparsity, p_curve)
      sp_index <- sp_fun$sparse_index
      sparsecase[sparse_ind, sp_index] <- NA
    }
  }
  return (sparsecase)  ### n * length(sett) matrix
}



observed_measure <- function(vec, sparsity, p_curve) {
  nc <- length(sett)
  if (sparsity == "point") {
    sparse_s <- sort(sample.int(nc, ceiling(nc * p_curve))) ## point sparsity with 20% missing in curves
  } else if (sparsity == "partial" || sparsity == "peak") {
    sparse_st <- base::sample(1:ceiling(nc - nc * p_curve), 1)
    sparse_s <- sparse_st:(sparse_st + ceiling(p_curve * nc))
  }
  if (p_curve == 0 || sparsity == "dense") {
    sparse_s <- length(sett) + 1
  }
  vec[sparse_s] <- NA
  vec <- vec[1:nc]
  return (list(sparse_index = sparse_s, value = vec))
}




plotsim_sparse <- function(sim, sparse_data, argval, 
                           xlab, ylab, title) {
  true_curve <- sim[[1]]
  observed_curve <- sim[[2]]
  out_index <- sim[[3]]
  
  #complete_result <- bootstrap_mfpca(30, sparsify, argval)
  par(mai = c(0.7, 0.8, 0.4, 0.1), mar = c(4.5, 4.2, 2.5, 1),
      mgp = c(3, 1.5, 0), mfrow = c(1, ifelse(p <= 3, p, 3)))
  for (j in 1:p) {
    sam <- observed_curve[[j]]
    sp <- sparse_data[[j]]
    plot(argval, sam[1, ], type = "n", 
         xlab = xlab[j], ylab = ylab[j], 
         cex.axis = 1.2, ylim = range(observed_curve), 
         cex.lab = 1.8, cex.main = 2, 
         main = title[j])
    apply(sam, 1, function(i) {
      lines(argval, i, col = "grey", lty = 2)
    })
    
    apply(sp[setdiff(1:n, out_index), ], 1,
          function(i) {
            lines(argval, i)})
    
    apply(sp, 1, function(i) {
      points(argval, i, cex = 0.4, pch = 1)}) ## pure
    if (length(out_index) >= 1) {
      apply(sp[out_index, ], 1, function(i) {
        lines(argval, i, col = "red")}) ## pure
    }
  }
  # plot_3d <- scatterplot3d(observed_curve[[1]][1, ], observed_curve[[2]][1, ], observed_curve[[3]][1, ],
  #                          type = "n", xlim = range(observed_curve[[1]]), ylim = range(observed_curve[[2]]),
  #                          zlim = range(observed_curve[[3]]),
  #               xlab = "Variable 1", ylab = "Variable 2", zlab = "Variable 3", main = "Three-dimensional Data")
  # lapply(1:n, function(k) {
  #   plot_3d$points3d(observed_curve[[1]][k, ], observed_curve[[2]][k, ], observed_curve[[3]][k, ],
  #                    type = "l", col = "grey", lty = 2)
  # })
  # lapply(setdiff(1:n, out_index), function(k) {
  #   data_transform <- data.frame(sparse_data[[1]][k, ], sparse_data[[2]][k, ], sparse_data[[3]][k, ])
  #   plot_3d$points3d(data_transform, type = "l")
  # })
  # lapply(1:n, function(k) {
  #   plot_3d$points3d(sparse_data[[1]][k, ], sparse_data[[2]][k, ], sparse_data[[3]][k, ],
  #                    type = "p", cex = 0.4, pch = 1)
  # })
  # lapply(out_index, function(k) {
  #   data_transform <- data.frame(sparse_data[[1]][k, ], sparse_data[[2]][k, ], sparse_data[[3]][k, ])
  #   plot_3d$points3d(data_transform, type = "l", col = "red")
  # })
  # lapply(out_index, function(k) {
  #   data_transform <- data.frame(observed_curve[[1]][k, ], observed_curve[[2]][k, ], observed_curve[[3]][k, ])
  #   plot_3d$points3d(data_transform, type = "l", col = "red")
  # })
  invisible(sim)
}

