# 第二部分
# 对4000个病例的存活情况进行描述性分析
library(ggplot2)
# 图中的超参数设置

rm(list=ls())
TITLE_SIZE = 30
TEXT_SIZE = 30

plothist = function(data, xtitle, title){
  precent = quantile(data$Length_of_stay, 0.95, na.rm = TRUE)
  p = ggplot(data = data,aes(x = Length_of_stay, y = ..density.. * nrow(data))) + 
    geom_histogram(bins = 50, color = 'black', fill = '#4784b7') + 
    ylab('count') +
    xlab(xtitle) + 
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.65, size = TITLE_SIZE, face = "bold")) + 
    theme(axis.text.x = element_text(size = TEXT_SIZE, face = "bold")) + 
    theme(axis.text.y = element_text(size = TEXT_SIZE, face = "bold")) + 
    theme(axis.title.x = element_text(size = TITLE_SIZE, face = "bold")) + 
    theme(axis.title.y = element_text(size = TITLE_SIZE, face = "bold")) + 
    geom_density(size = 0.5) + 
    geom_vline(xintercept = precent, color = 'red', linetype = 5, size = 0.5) +
    theme(plot.margin=unit(rep(1,4),'lines'))
  return(p)
}



df_a = read.csv('./data/Outcomes-a.txt',sep = ',')
df_b = read.csv('./data/Outcomes-b.txt',sep = ',')
df = rbind(df_a)

N = 4000
# 将-1替换为缺失值
df[df == -1] = NA
table(df$In.hospital_death) / nrow(df)

In_death = df[df$In.hospital_death == 1,]
Out_death = df[df$Survival > df$Length_of_stay & !is.na(df$Survival),]
survival = df[is.na(df$Survival),]
death = rbind(In_death, Out_death)

# 计算各部分比例
rate1 = c(nrow(death),nrow(survival)) / N
rate2 = c(nrow(In_death),nrow(Out_death)) / nrow(death)

# 总体的在医院停留天数的直方图
plothist(df, 'day', 'Distribution of Hospital Stay Days')
mean(df$Length_of_stay, na.rm = TRUE)

# 院内死亡病例的在医院停留天数的直方图
plothist(In_death, 'day', 'Length of stay (In Death)')
mean(In_death$Length_of_stay, na.rm = TRUE)

# 存活病例的在医院停留天数的直方图
plothist(survival, 'day', 'Length of stay (Survival)')
mean(survival$Length_of_stay, na.rm = TRUE)

# 院外死亡病例的在医院停留天数的直方图
plothist(Out_death, 'day', 'Length of stay (Out Death)')
mean(Out_death$Length_of_stay, na.rm = TRUE)

# 死亡病例的在医院停留天数的直方图
plothist(death, 'day', 'Length of stay (Death)')
mean(death$Length_of_stay, na.rm = TRUE)












