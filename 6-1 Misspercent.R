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
source('./function/MM-FOD.R')
source('./function/evaluation.R')
source("./function/Qu_Zhuo_JCGS/outl_detection.R")
source("./function/WXK-simulation_models.R")
source("./function/2-PlotOutlierData_v4.R")

## ----------Hyper_parameter-------------------
n = 100
LENGTH_OUT = 50 # p

model_list = c("Cos", "Flat", "Spike")

outliertype_list = c("I", "II", "III")

#MISS_PER = seq(0.1, 0.9, length.out = 9)
MISS_PER = c(0.05,0.1,0.3,0.5,0.7,0.9)
#MISS_PER = 0.05
TYPEs <- c(1,2,3)
#OUTLIER_PER = 0.05

LOOPNUM = 30

## -------------Test-----------
MODEL = 1
for(TYPE in TYPEs){
  cat("Type: ",TYPE,"\n")
for (missingness in MISS_PER) {
  print(missingness)
  result_tp_fp_1 = data.frame()
  result_tp_fp_2 = data.frame()
  for (i in 1:LOOPNUM) {

    argval = seq(0, 1, length.out = LENGTH_OUT)
    # (\mu(t) = (5 sin(2\pi t), 5 cos(2\pi t), 5(t − 1)^2))
    p=3
    sim = generate_samples("no outlier", "linear", argval)

    
    normalData = sim$observed_curve[[MODEL]]
    

    
    if (TYPE == 1) {
      dt = add_outlier_type1(normalData, p=50, zero_len = 10, zero_times = 2, n_outliers=10,
                             plot = TRUE,plot_title = "Outlier Type 1",xlabel ="Time")
    }else if (TYPE == 2) {
      #After I made some changes, outlier type 4 becomes the new type 2.
      dt = add_outlier_type4(normalData, p=50,q=10,  n_outliers=10,
                             plot = TRUE,plot_title = "Outlier Type 2",xlabel ="Time")
    }else if (TYPE == 3) {
      dt = add_outlier_type2(normalData, p=50,n_outliers=10,
                             plot = TRUE,plot_title = "Outlier Type 3",xlabel ="Time")
    }else if (TYPE == 4) {
      dt = add_outlier_type4(normalData, p=50, q=8, n_outliers=5,
                             plot = TRUE,plot_title = "Outlier Type 4",xlabel ="Time")
    }

    true_outlier  = dt$true_outliers
    sparse_data = uni_sparse(dt$data, p_size = 1,
                             p_curve = 1 - missingness, sparsity = "point")
  
    sample_1 = as.data.frame(sparse_data)
    Lt = list()
    Ly = list()
    for (i in 1:nrow(sample_1)) {
      Lt[[i]] = as.numeric(which(!is.na(sample_1[i,])))
      Ly[[i]] = as.numeric(sample_1[i,which(!is.na(sample_1[i,]))])
    }
    y_1 = list(Lt = Lt, Ly = Ly)

    ## MMFOD
    try({
      out_MM_FOD_1 = MM_FOD(y_1, iter = 1, k_sub = 0.6, alpha0 = 0.05,
                            n = n, nout = 5, p = LENGTH_OUT)

      tp_MM_FOD_1 = evaluation(n, out_MM_FOD_1, true_outlier)$tp
      fp_MM_FOD_1 = evaluation(n, out_MM_FOD_1, true_outlier)$fp
      result_line = c(tp_MM_FOD_1, fp_MM_FOD_1)
      result_tp_fp_1 = rbind(result_tp_fp_1, result_line)

    }, silent = TRUE)
  
    
    ## JCGS
    try({
      bootstrapTimes = 20
      confid_alpha = 0.05
      result_bs <- bootstrap_implementation(bootstrapTimes, argval,
                                            list(dt$data),
                                            list(sparse_data), confid_alpha)

      bmfpca_list <- lapply(1:length(argval), function(k) {
        sapply(result_bs$bootstrap_fit, function(l) {l[, k]})
      })
      
      depth_bmfpc <- multi_depth(bmfpca_list, "MFHD")

      fit = as.data.frame(result_bs$bootstrap_fit)
      sp <- sparse_fbplot(fit = fit,
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

      JCGS_outlier = unique(unlist(sp$fb_outlier_index))

      tp_JCGS_1 = evaluation(n, JCGS_outlier, true_outlier)$tp
      fp_JCGS_1 = evaluation(n, JCGS_outlier, true_outlier)$fp
      result_line = c(tp_JCGS_1, fp_JCGS_1)
      result_tp_fp_2 = rbind(result_tp_fp_2, result_line)
    }, silent = TRUE)


    }
  names(result_tp_fp_1) = c("point_tp", "point_fp")
  filepath1 = paste0("E:\\ICU-patient0529\\result\\MMFOD_", MODEL,
                     "_", TYPE, "_", missingness, ".csv")
  write.csv(result_tp_fp_1, filepath1, row.names = TRUE)

  names(result_tp_fp_2) = c("point_tp", "point_fp")
  filepath2 = paste0("E:\\ICU-patient0529\\result\\JCGS_", MODEL,
                     "_", TYPE, "_", missingness, ".csv")
  write.csv(result_tp_fp_2, filepath2, row.names = TRUE)
}
}
  
# ## ----------------Simulation----------------
start_time = Sys.time()
for (MODEL in 1:2) {
  for (TYPE in 1:4) {
    for (missingness in MISS_PER) {
      print(missingness)
      result_tp_fp_1 = data.frame()
      result_tp_fp_2 = data.frame()
      for (i in 1:LOOPNUM) {
        
        argval = seq(0, 1, length.out = LENGTH_OUT)
        # (\mu(t) = (5 sin(2\pi t), 5 cos(2\pi t), 5(t − 1)^2))
        sim = generate_samples("no outlier", "linear", argval)
        
        if (MODEL != 3) {
          normalData = sim$observed_curve[[MODEL + 1]]
        }else if (MODLE == 3) {
          
          df = read.csv("./data/combine-data-set-a/all_AST.csv")
          
          normalData = sim$observed_curve[[MODEL]] # use AST true data
        }
        
        
        if (TYPE == 1) {
          dt = add_outlier_type1(normalData, p=50, zero_len = 6, zero_times = 3, n_outliers=5, 
                                 plot = TRUE,plot_title = "Outlier Type 1",xlabel ="Time")
        }else if (TYPE == 2) {
          dt = add_outlier_type2(normalData, p=50, sin_coeff = 7, n_outliers=5, 
                                 plot = TRUE,plot_title = "Outlier Type 2",xlabel ="Time")
        }else if (TYPE == 3) {
          dt = add_outlier_type3(normalData, p=50, q=6, n_outliers=5,  
                                 plot = TRUE,plot_title = "Outlier Type 3",xlabel ="Time")
        }else if (TYPE == 4) {
          dt = add_outlier_type4(normalData, p=50, q=7, n_outliers=5,  
                                 plot = TRUE,plot_title = "Outlier Type 4",xlabel ="Time")
        }
        
        true_outlier  = dt$true_outliers
        sparse_data = uni_sparse(dt$data, p_size = 1,
                                 p_curve = 1 - missingness, sparsity = "point")
        
        sample_1 = as.data.frame(sparse_data)
        Lt = list()
        Ly = list()
        for (i in 1:nrow(sample_1)) {
          Lt[[i]] = as.numeric(which(!is.na(sample_1[i,])))
          Ly[[i]] = as.numeric(sample_1[i,which(!is.na(sample_1[i,]))])
        }
        y_1 = list(Lt = Lt, Ly = Ly)
        
        ## MMFOD
        try({
          out_MM_FOD_1 = MM_FOD(y_1, iter = 1, k_sub = 0.6, alpha0 = 0.05,
                                n = n, nout = 5, p = LENGTH_OUT)
          
          tp_MM_FOD_1 = evaluation(n, out_MM_FOD_1, true_outlier)$tp
          fp_MM_FOD_1 = evaluation(n, out_MM_FOD_1, true_outlier)$fp
          result_line = c(tp_MM_FOD_1, fp_MM_FOD_1)
          result_tp_fp_1 = rbind(result_tp_fp_1, result_line)
          
        }, silent = TRUE)
        
        ## JCGS
        try({
          bootstrapTimes = 20
          confid_alpha = 0.05
          result_bs <- bootstrap_implementation(bootstrapTimes, argval,
                                                list(dt$data),
                                                list(sparse_data), confid_alpha)
          
          bmfpca_list <- lapply(1:length(argval), function(k) {
            sapply(result_bs$bootstrap_fit, function(l) {l[, k]})
          })
          
          depth_bmfpc <- multi_depth(bmfpca_list, "MFHD")
          
          fit = as.data.frame(result_bs$bootstrap_fit)
          sp <- sparse_fbplot(fit = fit,
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
          
          JCGS_outlier = unique(unlist(sp$fb_outlier_index))
          
          tp_JCGS_1 = evaluation(n, JCGS_outlier, true_outlier)$tp
          fp_JCGS_1 = evaluation(n, JCGS_outlier, true_outlier)$fp
          result_line = c(tp_JCGS_1, fp_JCGS_1)
          result_tp_fp_2 = rbind(result_tp_fp_2, result_line)
        }, silent = TRUE)
        
      }
      names(result_tp_fp_1) = c("point_tp", "point_fp")
      filepath1 = paste0("./result/chap5result/Missbox/MMFOD_", MODEL,
                         "_", TYPE, "_", missingness, ".csv")
      write.csv(result_tp_fp_1, filepath1, row.names = TRUE)
      
      names(result_tp_fp_2) = c("point_tp", "point_fp")
      filepath2 = paste0("./result/chap5result/Missbox/JCGS_", MODEL,
                         "_", TYPE, "_", missingness, ".csv")
      write.csv(result_tp_fp_2, filepath2, row.names = TRUE)
      
    }
  }
}
end_time = Sys.time()
run_time = end_time - start_time








