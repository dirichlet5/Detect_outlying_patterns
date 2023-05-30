library("gdata")
library("fda")
library("compiler")
library("fdapace")
library("fdaoutlier")

rm(list = ls())

source("./function/Qu_Zhuo_JCGS/simulation-function.R")
source("./function/Qu_Zhuo_JCGS/simulation_parameter.R")
source("./function/Qu_Zhuo_JCGS/one_sample.R")
source("./function/Qu_Zhuo_JCGS/multi_sparse.R")
source("./function/MM-FOD.R")
source("./function/evaluation.R")
source("./function/Qu_Zhuo_JCGS/outl_detection.R")
source("./function/WXK-simulation_models.R")
source("./function/2-PlotOutlierData_v4.R")

n <- 100
OUTLIER_LIST <- c(0.05, 0.1, 0.15)
OUTLIER_LIST <- c(0.1, 0.15)
# OUTLIER_LIST <- c(0.05)
# OUTLIER_LIST <- c(0.3, 0.4)
LENGTH_OUT <- 50 # p

LOOPNUM <- 10

start_time <- Sys.time()

## --------MMFOD---TYPE1--------
TYPE <- 1
zero_len_list <- seq(2, 10, 1)
for (MODEL in 1:2) {
  for (outlier_rate in OUTLIER_LIST) {
    dfMean <- data.frame()
    rowname_list <- c()
    for (zero_len in zero_len_list) {
      parameter_name <- paste0("zero_len=", zero_len)
      rowname_list <- c(rowname_list, parameter_name)
      result_tp_fp <- data.frame()
      for (i in 1:LOOPNUM) {
        argval <- seq(0, 1, length.out = LENGTH_OUT)
        sim <- generate_samples("no outlier", "linear", argval)
        normalData <- sim$observed_curve[[5 - 2 * MODEL]]

        dt <- add_outlier_type1(normalData,
          p = 50, zero_len = zero_len, zero_times = 3, n_outliers = n * outlier_rate,
          plot = TRUE, plot_title = "Outlier Type 1", xlabel = "Time"
        )

        true_outlier <- dt$true_outliers
        sparse_data <- uni_sparse(dt$data,
          p_size = 1,
          p_curve = 0.3, sparsity = "point"
        )

        sample_1 <- as.data.frame(sparse_data)
        Lt <- list()
        Ly <- list()
        for (i in 1:nrow(sample_1)) {
          Lt[[i]] <- as.numeric(which(!is.na(sample_1[i, ])))
          Ly[[i]] <- as.numeric(sample_1[i, which(!is.na(sample_1[i, ]))])
        }
        y_1 <- list(Lt = Lt, Ly = Ly)

        ## MMFOD
        try(
          {
            out_MM_FOD_1 <- MM_FOD(y_1,
              iter = 1, k_sub = 0.6, alpha0 = 0.05,
              n = n, nout = 5, p = LENGTH_OUT
            )

            tp_MM_FOD_1 <- evaluation(n, out_MM_FOD_1, true_outlier)$tp
            fp_MM_FOD_1 <- evaluation(n, out_MM_FOD_1, true_outlier)$fp
            bcr_MM_FOD_1 <- BCR(n, out_MM_FOD_1, true_outlier)
            result_line <- c(tp_MM_FOD_1, fp_MM_FOD_1, bcr_MM_FOD_1)
            result_tp_fp <- rbind(result_tp_fp, result_line)
          },
          silent = TRUE
        )
      }
      dfMean <- rbind(dfMean, colMeans(result_tp_fp))
    }
    names(dfMean) <- c("point_tp", "point_fp", "bcr")
    row.names(dfMean) <- rowname_list
    filepath <- paste0(
      "./result/ROC_result/MMFOD/mod_", MODEL, "_TYPE_", TYPE,
      "_outlierRate_", outlier_rate, ".csv"
    )
    write.csv(dfMean, filepath, row.names = TRUE)
  }
}

## --------MMFOD---TYPE2--------
TYPE <- 2
q_list <- seq(6, 10, 0.5)
for (MODEL in 1:2) {
  for (outlier_rate in OUTLIER_LIST) {
    dfMean <- data.frame()
    rowname_list <- c()
    for (q in q_list) {
      parameter_name <- paste0("q=", q)
      rowname_list <- c(rowname_list, parameter_name)
      result_tp_fp <- data.frame()
      for (i in 1:LOOPNUM) {
        argval <- seq(0, 1, length.out = LENGTH_OUT)
        sim <- generate_samples("no outlier", "linear", argval)
        normalData <- sim$observed_curve[[5 - 2 * MODEL]]
        
        dt <- add_outlier_type4(normalData,
                                p = 50, q = q, n_outliers = n * outlier_rate,
                                plot = TRUE, plot_title = "Outlier Type 2", xlabel = "Time"
        )
        
        true_outlier <- dt$true_outliers
        sparse_data <- uni_sparse(dt$data,
                                  p_size = 1,
                                  p_curve = 0.3, sparsity = "point"
        )
        
        sample_1 <- as.data.frame(sparse_data)
        Lt <- list()
        Ly <- list()
        for (i in 1:nrow(sample_1)) {
          Lt[[i]] <- as.numeric(which(!is.na(sample_1[i, ])))
          Ly[[i]] <- as.numeric(sample_1[i, which(!is.na(sample_1[i, ]))])
        }
        y_1 <- list(Lt = Lt, Ly = Ly)
        
        ## MMFOD
        try(
          {
            out_MM_FOD_1 <- MM_FOD(y_1,
                                   iter = 1, k_sub = 0.6, alpha0 = 0.05,
                                   n = n, nout = 5, p = LENGTH_OUT
            )
            
            tp_MM_FOD_1 <- evaluation(n, out_MM_FOD_1, true_outlier)$tp
            fp_MM_FOD_1 <- evaluation(n, out_MM_FOD_1, true_outlier)$fp
            bcr_MM_FOD_1 <- BCR(n, out_MM_FOD_1, true_outlier)
            result_line <- c(tp_MM_FOD_1, fp_MM_FOD_1, bcr_MM_FOD_1)
            result_tp_fp <- rbind(result_tp_fp, result_line)
          },
          silent = TRUE
        )
      }
      dfMean <- rbind(dfMean, colMeans(result_tp_fp))
    }
    names(dfMean) <- c("point_tp", "point_fp", "bcr")
    row.names(dfMean) <- rowname_list
    filepath <- paste0(
      "./result/ROC_result/MMFOD/mod_", MODEL, "_TYPE_", TYPE,
      "_outlierRate_", outlier_rate, ".csv"
    )
    write.csv(dfMean, filepath, row.names = TRUE)
  }
}


## --------MMFOD---TYPE3--------
TYPE <- 3
sin_coeff_list <- seq(6, 10, 0.5)
for (MODEL in 1:2) {
  for (outlier_rate in OUTLIER_LIST) {
    dfMean <- data.frame()
    rowname_list <- c()
    for (sin_coeff in sin_coeff_list) {
      parameter_name <- paste0("sin_coeff=", sin_coeff)
      rowname_list <- c(rowname_list, parameter_name)
      result_tp_fp <- data.frame()
      for (i in 1:LOOPNUM) {
        argval <- seq(0, 1, length.out = LENGTH_OUT)
        sim <- generate_samples("no outlier", "linear", argval)
        normalData <- sim$observed_curve[[5 - 2 * MODEL]]

        dt <- add_outlier_type2(normalData,
          p = 50, sin_coeff = sin_coeff, n_outliers = n * outlier_rate,
          plot = TRUE, plot_title = "Outlier Type 3", xlabel = "Time"
        )

        true_outlier <- dt$true_outliers
        sparse_data <- uni_sparse(dt$data,
          p_size = 1,
          p_curve = 0.3, sparsity = "point"
        )

        sample_1 <- as.data.frame(sparse_data)
        Lt <- list()
        Ly <- list()
        for (i in 1:nrow(sample_1)) {
          Lt[[i]] <- as.numeric(which(!is.na(sample_1[i, ])))
          Ly[[i]] <- as.numeric(sample_1[i, which(!is.na(sample_1[i, ]))])
        }
        y_1 <- list(Lt = Lt, Ly = Ly)

        ## MMFOD
        try(
          {
            out_MM_FOD_1 <- MM_FOD(y_1,
              iter = 1, k_sub = 0.6, alpha0 = 0.05,
              n = n, nout = 5, p = LENGTH_OUT
            )

            tp_MM_FOD_1 <- evaluation(n, out_MM_FOD_1, true_outlier)$tp
            fp_MM_FOD_1 <- evaluation(n, out_MM_FOD_1, true_outlier)$fp
            bcr_MM_FOD_1 <- BCR(n, out_MM_FOD_1, true_outlier)
            result_line <- c(tp_MM_FOD_1, fp_MM_FOD_1, bcr_MM_FOD_1)
            result_tp_fp <- rbind(result_tp_fp, result_line)
          },
          silent = TRUE
        )
      }
      dfMean <- rbind(dfMean, colMeans(result_tp_fp))
    }
    names(dfMean) <- c("point_tp", "point_fp", "bcr")
    row.names(dfMean) <- rowname_list
    filepath <- paste0(
      "./result/ROC_result/MMFOD/mod_", MODEL, "_TYPE_", TYPE,
      "_outlierRate_", outlier_rate, ".csv"
    )
    write.csv(dfMean, filepath, row.names = TRUE)
  }
}


## -------JCGS----TYPE1--------
depth_list <- c("MFHD", "MFPD", "MFDPD", "MFSPD")
TYPE <- 1
zero_len_list <- seq(2, 10, 1)
for (MODEL in 1:2) {
  for (depth in depth_list) {
    for (outlier_rate in OUTLIER_LIST) {
      dfMean <- data.frame()
      rowname_list <- c()

      for (zero_len in zero_len_list) {
        parameter_name <- paste0("zero_len=", zero_len)
        rowname_list <- c(rowname_list, parameter_name)
        result_tp_fp <- data.frame()

        for (i in 1:LOOPNUM) {
          argval <- seq(0, 1, length.out = LENGTH_OUT)
          sim <- generate_samples("no outlier", "linear", argval)
          normalData <- sim$observed_curve[[5 - 2 * MODEL]]

          dt <- add_outlier_type1(normalData,
            p = 50, zero_len = zero_len, zero_times = 3, n_outliers = n * outlier_rate,
            plot = TRUE, plot_title = "Outlier Type 1", xlabel = "Time"
          )

          true_outlier <- dt$true_outliers
          sparse_data <- uni_sparse(dt$data,
            p_size = 1,
            p_curve = 0.3, sparsity = "point"
          )

          try(
            {
              bootstrapTimes <- 20
              confid_alpha <- 0.05
              result_bs <- bootstrap_implementation(
                bootstrapTimes, argval,
                list(dt$data),
                list(sparse_data), confid_alpha
              )

              bmfpca_list <- lapply(1:length(argval), function(k) {
                sapply(result_bs$bootstrap_fit, function(l) {
                  l[, k]
                })
              })

              depth_bmfpc <- multi_depth(bmfpca_list, depth)

              fit <- as.data.frame(result_bs$bootstrap_fit)
              sp <- sparse_fbplot(
                fit = fit,
                sparse = sparse_data,
                time_index = NULL, depth = depth_bmfpc,
                two_stage = FALSE, sq_vo = sq_vo, plot = FALSE,
                xlab = NULL, ylab = NULL, title = NULL,
                yrange = NULL,
                cex.main = 1.3, cex.lab = 1.3, cex.axis = 1.3,
                medlabel = TRUE, outlabel = TRUE,
                prob = 0.5, factor = 1.5,
                color = 6, outliercolor.fb = 2, barcol = 4,
                outliercolor.dir = 3, fullout = FALSE
              )

              JCGS_outlier <- unique(unlist(sp$fb_outlier_index))

              tp_JCGS_1 <- evaluation(n, JCGS_outlier, true_outlier)$tp
              fp_JCGS_1 <- evaluation(n, JCGS_outlier, true_outlier)$fp
              bcr_JCGS_1 <- BCR(n, JCGS_outlier, true_outlier)
              result_line <- c(tp_JCGS_1, fp_JCGS_1, bcr_JCGS_1)
              result_tp_fp <- rbind(result_tp_fp, result_line)
            },
            silent = TRUE
          )
        }
        dfMean <- rbind(dfMean, colMeans(result_tp_fp))
      }
      names(dfMean) <- c("point_tp", "point_fp", "bcr")
      row.names(dfMean) <- rowname_list
      filepath <- paste0(
        "./result/ROC_result/JCGS_", depth, "/mod_", MODEL, "_TYPE_", TYPE,
        "_outlierRate_", outlier_rate, ".csv"
      )
      write.csv(dfMean, filepath, row.names = TRUE)
    }
  }
}


## -------JCGS----TYPE2--------
TYPE <- 2
q_list <- seq(6, 10, 0.5)
for (MODEL in 1:2) {
  for (depth in depth_list) {
    for (outlier_rate in OUTLIER_LIST) {
      dfMean <- data.frame()
      rowname_list <- c()
      
      for (q in q_list) {
        parameter_name <- paste0("q=", q)
        rowname_list <- c(rowname_list, parameter_name)
        result_tp_fp <- data.frame()
        
        for (i in 1:LOOPNUM) {
          argval <- seq(0, 1, length.out = LENGTH_OUT)
          sim <- generate_samples("no outlier", "linear", argval)
          normalData <- sim$observed_curve[[5 - 2 * MODEL]]
          
          dt <- add_outlier_type4(normalData,
                                  p = 50, q = q, n_outliers = n * outlier_rate,
                                  plot = TRUE, plot_title = "Outlier Type 2", xlabel = "Time"
          )
          
          true_outlier <- dt$true_outliers
          sparse_data <- uni_sparse(dt$data,
                                    p_size = 1,
                                    p_curve = 0.3, sparsity = "point"
          )
          
          try(
            {
              bootstrapTimes <- 20
              confid_alpha <- 0.05
              result_bs <- bootstrap_implementation(
                bootstrapTimes, argval,
                list(dt$data),
                list(sparse_data), confid_alpha
              )
              
              bmfpca_list <- lapply(1:length(argval), function(k) {
                sapply(result_bs$bootstrap_fit, function(l) {
                  l[, k]
                })
              })
              
              depth_bmfpc <- multi_depth(bmfpca_list, depth)
              
              fit <- as.data.frame(result_bs$bootstrap_fit)
              sp <- sparse_fbplot(
                fit = fit,
                sparse = sparse_data,
                time_index = NULL, depth = depth_bmfpc,
                two_stage = FALSE, sq_vo = sq_vo, plot = FALSE,
                xlab = NULL, ylab = NULL, title = NULL,
                yrange = NULL,
                cex.main = 1.3, cex.lab = 1.3, cex.axis = 1.3,
                medlabel = TRUE, outlabel = TRUE,
                prob = 0.5, factor = 1.5,
                color = 6, outliercolor.fb = 2, barcol = 4,
                outliercolor.dir = 3, fullout = FALSE
              )
              
              JCGS_outlier <- unique(unlist(sp$fb_outlier_index))
              
              tp_JCGS_1 <- evaluation(n, JCGS_outlier, true_outlier)$tp
              fp_JCGS_1 <- evaluation(n, JCGS_outlier, true_outlier)$fp
              bcr_JCGS_1 <- BCR(n, JCGS_outlier, true_outlier)
              result_line <- c(tp_JCGS_1, fp_JCGS_1, bcr_JCGS_1)
              result_tp_fp <- rbind(result_tp_fp, result_line)
            },
            silent = TRUE
          )
        }
        dfMean <- rbind(dfMean, colMeans(result_tp_fp))
      }
      names(dfMean) <- c("point_tp", "point_fp", "bcr")
      row.names(dfMean) <- rowname_list
      filepath <- paste0(
        "./result/ROC_result/JCGS_", depth, "/mod_", MODEL, "_TYPE_", TYPE,
        "_outlierRate_", outlier_rate, ".csv"
      )
      write.csv(dfMean, filepath, row.names = TRUE)
    }
  }
}



## -------JCGS----TYPE3--------
TYPE <- 3
sin_coeff_list <- seq(6, 10, 0.5)
for (MODEL in 1:2) {
  for (depth in depth_list) {
    for (outlier_rate in OUTLIER_LIST) {
      dfMean <- data.frame()
      rowname_list <- c()

      for (sin_coeff in sin_coeff_list) {
        parameter_name <- paste0("sin_coeff=", sin_coeff)
        rowname_list <- c(rowname_list, parameter_name)
        result_tp_fp <- data.frame()

        for (i in 1:LOOPNUM) {
          argval <- seq(0, 1, length.out = LENGTH_OUT)
          sim <- generate_samples("no outlier", "linear", argval)
          normalData <- sim$observed_curve[[5 - 2 * MODEL]]

          dt <- add_outlier_type2(normalData,
            p = 50, sin_coeff = sin_coeff, n_outliers = n * outlier_rate,
            plot = TRUE, plot_title = "Outlier Type 3", xlabel = "Time"
          )

          true_outlier <- dt$true_outliers
          sparse_data <- uni_sparse(dt$data,
            p_size = 1,
            p_curve = 0.3, sparsity = "point"
          )

          try(
            {
              bootstrapTimes <- 20
              confid_alpha <- 0.05
              result_bs <- bootstrap_implementation(
                bootstrapTimes, argval,
                list(dt$data),
                list(sparse_data), confid_alpha
              )

              bmfpca_list <- lapply(1:length(argval), function(k) {
                sapply(result_bs$bootstrap_fit, function(l) {
                  l[, k]
                })
              })

              depth_bmfpc <- multi_depth(bmfpca_list, depth)

              fit <- as.data.frame(result_bs$bootstrap_fit)
              sp <- sparse_fbplot(
                fit = fit,
                sparse = sparse_data,
                time_index = NULL, depth = depth_bmfpc,
                two_stage = FALSE, sq_vo = sq_vo, plot = FALSE,
                xlab = NULL, ylab = NULL, title = NULL,
                yrange = NULL,
                cex.main = 1.3, cex.lab = 1.3, cex.axis = 1.3,
                medlabel = TRUE, outlabel = TRUE,
                prob = 0.5, factor = 1.5,
                color = 6, outliercolor.fb = 2, barcol = 4,
                outliercolor.dir = 3, fullout = FALSE
              )

              JCGS_outlier <- unique(unlist(sp$fb_outlier_index))

              tp_JCGS_1 <- evaluation(n, JCGS_outlier, true_outlier)$tp
              fp_JCGS_1 <- evaluation(n, JCGS_outlier, true_outlier)$fp
              bcr_JCGS_1 <- BCR(n, JCGS_outlier, true_outlier)
              result_line <- c(tp_JCGS_1, fp_JCGS_1, bcr_JCGS_1)
              result_tp_fp <- rbind(result_tp_fp, result_line)
            },
            silent = TRUE
          )
        }
        dfMean <- rbind(dfMean, colMeans(result_tp_fp))
      }
      names(dfMean) <- c("point_tp", "point_fp", "bcr")
      row.names(dfMean) <- rowname_list
      filepath <- paste0(
        "./result/ROC_result/JCGS_", depth, "/mod_", MODEL, "_TYPE_", TYPE,
        "_outlierRate_", outlier_rate, ".csv"
      )
      write.csv(dfMean, filepath, row.names = TRUE)
    }
  }
}

end_time <- Sys.time()
print(end_time - start_time)

