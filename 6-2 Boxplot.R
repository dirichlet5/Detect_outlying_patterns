

library(ggplot2)
library(patchwork)


rm(list = ls())

# demo
MODEL = 2
TYPE = 3
MISS_PER = c(0.05,0.1,0.3,0.5,0.7)
#missingness = 0.1

box_df = data.frame()
for (missingness in MISS_PER) {
  filepath = paste0("E:\\ICU-patient0529\\result\\MMFOD_", MODEL, "_", 
                    TYPE, "_", missingness, ".csv")
  dt = read.csv(filepath)
  
  temp = data.frame(dt$point_fp, missingness, method = "MMFOD")
  box_df = rbind(box_df, temp)
}

for (missingness in MISS_PER) {
  filepath = paste0("E:\\ICU-patient0529\\result\\JCGS_", MODEL, "_", 
                    TYPE, "_", missingness, ".csv")
  dt = read.csv(filepath)
  
  temp = data.frame(dt$point_fp, missingness, method = "JCGS")
  box_df = rbind(box_df, temp)
}

box_df$method = as.factor(box_df$method)
box_df$missingness = as.factor(box_df$missingness)

ggplot(box_df, aes(x = box_df$missingness, y = box_df$dt.point_fp, 
                   fill = box_df$method)) + 
  geom_boxplot()+scale_fill_brewer(palette="Set1")+
  theme(axis.title.x =element_text(size=30),axis.title.y =element_text(size=30),
        axis.text.x =element_text(size=30),axis.text.y =element_text(size=30)
        ,legend.position = "none")+
  labs(x="Missingness", y="FPR")


## -----------boxplot function--------------
twobox = function(MODEL = 1, TYPE = 1){
  box_df = data.frame()
  for (missingness in MISS_PER) {
    filepath = paste0("./result/chap5result/Missbox/MMFOD_", MODEL, "_", 
                      TYPE, "_", missingness, ".csv")
    dt = read.csv(filepath)
    
    temp = data.frame(dt$point_tp, missingness, method = "MMFOD")
    box_df = rbind(box_df, temp)
  }
  
  for (missingness in MISS_PER) {
    filepath = paste0("./result/chap5result/Missbox/JCGS_", MODEL, "_", 
                      TYPE, "_", missingness, ".csv")
    dt = read.csv(filepath)
    
    temp = data.frame(dt$point_tp, missingness, method = "JCGS")
    box_df = rbind(box_df, temp)
  }
  
  box_df$method = as.factor(box_df$method)
  box_df$missingness = as.factor(box_df$missingness)
  
  ggplot(box_df, aes(x = box_df$missingness, y = box_df$dt.point_tp, 
                     fill = box_df$method)) + 
    geom_boxplot() + 
    ggtitle(paste0("Model:", MODEL, " Type:", TYPE))
  
}

## --------loop-----------
pdf(file = "boxplot.pdf", width = 19, height = 15)

twobox(1, 1) + twobox(2, 1) +
twobox(1, 2) + twobox(2, 2) +
twobox(1, 3) + twobox(2, 3) +
twobox(1, 4) + twobox(2, 4) +
  plot_layout(ncol = 2)
dev.off()


