### 1. Load functions and packages.
source("./function/Qu_Zhuo_JCGS/simulation-function.R")
source("./function/Qu_Zhuo_JCGS/simulation_parameter.R")
source("./function/Qu_Zhuo_JCGS/one_sample.R")
source("./function/Qu_Zhuo_JCGS/multi_sparse.R")

### 2. Parameters
n = 10  # 30
MISS_PER = 0.65
MISS_PER1 = 0.35
argval = seq(0, 1, length.out = 50)
outlierst = "no outlier"
eigenvalue_type = "linear"

##--- 3. Defined Functions
plot_sparse = function(MODEL = 2){
  observed_curve = sim[[2]]
  out_index = sim[[3]]
  sp = sparse_data[[MODEL]]
  sam= observed_curve[[MODEL]]
  
  title = paste0("N:",n, " MISS_PER:", MISS_PER)
  par(mar=c(4,6,1,1))
  plot(argval, sam[1, ], type = "n", 
       xlab = "Time", ylab = "Value", 
       cex.axis = 2.5, ylim = range(sam)*1.2, 
       cex.lab = 3, cex.main = 2,mai=c(2,5,2,2),
       main = "")
  apply(sam, 1, function(i) {
    # lines(argval, i, col = "grey", lty = 3,lwd=6)
    points(argval, i, col = "grey", cex = 0.8, pch = 19)
  })
  
  apply(sp[setdiff(1:n, out_index), ], 1,
        function(i) {
          lines(argval, i,lwd=2)})
  
  apply(sp, 1, function(i) {
    points(argval, i, cex = 0.8, pch = 19)}) 
  if (length(out_index) >= 1) {
    apply(sp[out_index, ], 1, function(i) {
      lines(argval, i, col = "red", lty = 3,lwd=3)}) 
  }
}



##--- 4. Generate Data

# (\mu(t) = (5 sin(2\pi t)+10, 5 cos(2\pi t), 5(t − 1.4)^2+5))
sim = generate_samples(outlierst, eigenvalue_type, argval) #Note that the outliers here are generated separately.

# p_size: sparse percent in N samples
# p_curve: missing percent in one sample 
# sparsity: "point", "peak", "partial"
MISS_PER1 = MISS_PER[1] #输入进去的是个大于0.5的数
sparse_data = multi_sparse(sim$observed_curve, p_size = rep(1, p),
                           p_curve = rep((1-MISS_PER1), p), 
                           sparsity = c("point", "point", "point"))

#true_curve = sim[[1]]


##----- 5. Plot graphs

# mu2, 5 sin(2\pi t)+8
plot_sparse(MODEL = 1)

# mu3, 5(t − 1.4)^2+5
plot_sparse(MODEL = 3)


#####----- 6. Get the real data that is similar to the simulated data. 
{
  library(pracma)
  VAR_NAME = "HR"
  START_TIME = 27*60
  END_TIME = 42*60
  
  # VAR_NAME = "Urine"
  # START_TIME = 10*60
  # END_TIME = 30*60
  
  path = paste0("./data/combine-data-set-a/all_", VAR_NAME, ".csv")
  df = read.csv(path)
  
  df$Time = as.POSIXct(df$Time, format = '%Y-%m-%d %H:%M:%S')
  
  
  
  df_sec = df[START_TIME:END_TIME, -1]
  df_time = df[START_TIME:END_TIME, 1]
  
  sample_Num = 15
  samp = sample(1:ncol(df_sec), sample_Num, replace = FALSE)
  rand_samp = df_sec[,samp]
  
  argval = seq(0, 1, length.out = nrow(rand_samp))
  par(mar=c(4,6,1,1))
  plot(argval, rand_samp[, 1], type = "n",
       xlab = "Time", ylab = "Value",
       cex.axis = 2.5, ylim = c(range(rand_samp, na.rm = T)[1]-15,
                                range(rand_samp, na.rm = T)[2]+15), 
       cex.lab = 3, cex.main = 2,
       main = "")
  
  for (i in 1:sample_Num) {
    grey_data = interp1(argval, rand_samp[, i])
    lines(argval, grey_data, col = "grey", lty = 2,lwd=6)
    pos = as.numeric(which(!is.na(rand_samp[, i])))
    pos1 = pos
    for (len in 1:10) {
      pos1 = c(pos1, pos + len)
    }
    for (j in 1:length(grey_data)) {
      if (!(j %in% pos1)) {
        grey_data[j] = NA
      }
    }
    lines(argval, grey_data, type = "l",lwd=6)
    points(argval, rand_samp[, i], cex = 1.5, pch = 1) 
  }
  
}

#####--- 7. Get the real data for AST, which looks like the triangular waveform.
{
  library(dplyr) # 使用管道函数
  library(ggplot2)
  library(plotly)
  df_Type = read.csv('./data/Outcomes-a-Type.csv')
  var_name_list = read.table('./data/var_name2.txt',sep = '\t',header = FALSE)
  var_name = 'AST'
  path = paste0('./data/combine-data-set-a/all_', var_name, '.csv')
  df = read.csv(path)
  df$Time = as.POSIXct(df$Time, format = '%Y-%m-%d %H:%M:%S')
  N = 100
  dense_data = function(type, sample_num = 100){
    ID_type = df_Type[df_Type$Survival_type == type, 1]
    var_ID = as.numeric(gsub('X', '', colnames(df)[-1]))
    
    # 从该变量的生存类型为1的病例ID中随机抽取样本量N
    ID_sample = sample(intersect(ID_type, var_ID), size = N)
    
    # 读取这些ID所在的位置并根据位置选择数据
    ID_pos = c()
    for (i in ID_sample) {
      ID_pos = c(ID_pos, which(var_ID == i) + 1)
    }
    data_type = df[, ID_pos]
    
    # 计算这些样本量的均值并忽视其中的NA
    data_mean = data_type %>% apply(1, mean, na.rm = T)
    
    # 评估该操作过后的缺失值比例
    NA_rate = sum(is.na(data_mean))/length(data_mean)
    # sprintf('The NA rate in data is: %.2f%%',NA_rate * 100)
    return(data_mean)
  }
  
  res <- matrix(0,120,300) #Here we use a sampling method
  for(i in 1:120)
  {
    data_mean1 = dense_data(type = 1, sample_num = N)
    data_mean22<-na.omit(data_mean1)
    res[i,1:length(data_mean22)] <- data_mean22
  }
  
  
  #Here I generate the triangular waveform, plot the dense data.
  sam = res[sample(1:length(res[,1]),20,replace = F),1:50]
 
  sp <- matrix(NA,dim(sam)[1],dim(sam)[2])
  for(i in 1:length(sam[,1]))
  {
    index<- sort(sample(1:dim(sam)[2],floor(dim(sam)[2]*(1-MISS_PER)),replace = F))
    sp[i,index] <- sam[i,index]                                                                                                  
  } 
  
  par(mar=c(4,6,1,1))
  plot(argval, sam[1, ], type = "n", 
       xlab = "Time", ylab = "Value", 
       cex.axis = 2.5, ylim = range(sam), 
       cex.lab = 3, cex.main = 2,mai=c(2,5,2,2),
       main = "")
  apply(sam, 1, function(i) {
    lines(argval, i, col = "grey", lty = 3,lwd=5)
    })
  apply(sp, 1, function(i) {
    points(argval, i, cex = 1.5, pch = 16)})
  
}







