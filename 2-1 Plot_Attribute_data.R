library(ggplot2)
library(plotrix)

rm(list=ls())
# ---------------Age------------------
data_Age = read.csv('./data/combine-data-set-a/all_Age.csv')
N = ncol(data_Age) - 1
data_Age = as.numeric(data_Age[-1]) # 第一列为时间，删除

PLOT_SIZE = 30

# 分区间计算频数
age_count = table(cut(data_Age, breaks = c(0, 40, 60, 80, max(data_Age))))

# 年龄分布的饼图
info = as.numeric(age_count)
info / N
names = c("under 40 10.86%", "40-60 27.15%", "60-80 42.20%", "over 80 19.78%")
cols = c("#2f83e4", "#00e5c1", "#23cbff", "#29527e")
pie3D(info,
      labels = names,
      explode = 0.1,
      col = cols,
      main = "Age distribution map")

# 年龄分布的直方图
Age = data.frame(data_Age)
# 统计年龄变量中的缺失值
table(Age == -1) 
table(Age < 15)  # 年龄中不存在15岁以下，故身高中的100以下认为是异常值

# 计算Age95%分位数点
age_95 = quantile(as.numeric(Age[,1]), 0.95)  # 89

p1 = ggplot(data = Age,aes(x = data_Age, y = ..density.. * N)) + 
  geom_histogram(bins = 25, color = 'black', fill = '#4784b7') + 
  ylab('count') +
  xlab('Age') + 
  ggtitle('Age distribution of ICU patients') + 
  theme(plot.title = element_text(hjust = 0.8, size = PLOT_SIZE, face = "bold")) + 
  theme(axis.text.x = element_text(size = PLOT_SIZE, face = "bold")) + 
  theme(axis.text.y = element_text(size = PLOT_SIZE, face = "bold")) + 
  theme(axis.title.x = element_text(size = PLOT_SIZE, face = "bold")) + 
  theme(axis.title.y = element_text(size = PLOT_SIZE, face = "bold")) +
  geom_density(size = 0.5)
p1


# -----------------Height---------------
# 身高分布的直方图
data_Gender = read.csv('./data/combine-data-set-a/all_Gender.csv')
data_Gender = as.numeric(data_Gender[-1])
data_Height = read.csv('./data/combine-data-set-a/all_Height.csv')
data_Height = as.numeric(data_Height[-1])
Height = data.frame(data_Gender, data_Height)

# 统计身高变量中的缺失值情况
table(Height$data_Height == -1)  # 1894个缺失值


Height[Height$data_Height == -1 | Height$data_Height > 300 | Height$data_Height < 100, 2] = NA
man_Height = data.frame(Height = Height[data_Gender == 1, 2])
woman_Height = data.frame(Height = Height[data_Gender == 0, 2])

# 男性身高
p2 = ggplot(data = man_Height,aes(x = Height, y = ..density.. * N)) + 
  geom_histogram(bins = 25, color = 'black', fill = '#4784b7') + 
  ylab('count') + 
  xlab('Height/cm') + 
  ggtitle('Male Height Distribution') + 
  geom_density(size = 0.5) +
  theme(plot.title = element_text(hjust = 0.5, size = PLOT_SIZE, face = "bold")) + 
  theme(axis.text.x = element_text(size = PLOT_SIZE, face = "bold")) + 
  theme(axis.text.y = element_text(size = PLOT_SIZE, face = "bold")) + 
  theme(axis.title.x = element_text(size = PLOT_SIZE, face = "bold")) + 
  theme(axis.title.y = element_text(size = PLOT_SIZE, face = "bold"))
p2

# 女性身高
p3 = ggplot(data = woman_Height,aes(x = Height, y = ..density.. * N)) + 
  geom_histogram(bins = 25, color = 'black', fill = '#4784b7') + 
  ylab('count') + 
  xlab('Height/cm') + 
  ggtitle('Female Height Distribution') + 
  geom_density(size = 0.5) + 
  theme(plot.title = element_text(hjust = 0.5, size = PLOT_SIZE, face = "bold")) + 
  theme(axis.text.x = element_text(size = PLOT_SIZE, face = "bold")) + 
  theme(axis.text.y = element_text(size = PLOT_SIZE, face = "bold")) + 
  theme(axis.title.x = element_text(size = PLOT_SIZE, face = "bold")) + 
  theme(axis.title.y = element_text(size = PLOT_SIZE, face = "bold"))
p3


# --------------Weight--------------
data_Weight = read.csv('./data/combine-data-set-a/all_Weight.csv')
data_Weight = as.numeric(data_Weight[1, -1])
Weight = data.frame(data_Gender, data_Weight)

# 统计体重中的缺失值
table(Weight$data_Weight == -1) # 326个缺失值

Weight[Weight$data_Weight == -1, 2] = NA
man_Weight = data.frame(Weight = Weight[data_Gender == 1, 2])
woman_Weight = data.frame(Weight = Weight[data_Gender == 0, 2])

# 男性体重
p4 = ggplot(data = man_Weight,aes(x = Weight, y = ..density.. * N)) + 
  geom_histogram(bins = 25, color = 'black', fill = '#4784b7') + 
  ylab('count') + 
  xlab('Weight/kg') + 
  ggtitle('Male Weight Distribution') + 
  geom_density(size = 0.5) + 
  theme(plot.title = element_text(hjust = 0.5, size = PLOT_SIZE, face = "bold")) + 
  theme(axis.text.x = element_text(size = PLOT_SIZE, face = "bold")) + 
  theme(axis.text.y = element_text(size = PLOT_SIZE, face = "bold")) + 
  theme(axis.title.x = element_text(size = PLOT_SIZE, face = "bold")) + 
  theme(axis.title.y = element_text(size = PLOT_SIZE, face = "bold"))
p4

# 女性体重
p5 = ggplot(data = woman_Weight,aes(x = Weight, y = ..density.. * N)) + 
  geom_histogram(bins = 25, color = 'black', fill = '#4784b7') + 
  ylab('count') + 
  xlab('Weight/kg') + 
  ggtitle('Female Weight Distribution') + 
  geom_density(size = 0.5) + 
  theme(plot.title = element_text(hjust = 0.5, size = PLOT_SIZE, face = "bold")) + 
  theme(axis.text.x = element_text(size = PLOT_SIZE, face = "bold")) + 
  theme(axis.text.y = element_text(size = PLOT_SIZE, face = "bold")) + 
  theme(axis.title.x = element_text(size = PLOT_SIZE, face = "bold")) + 
  theme(axis.title.y = element_text(size = PLOT_SIZE, face = "bold"))
p5



# ---------------ICU-Type-----------------
data_ICUType = read.csv('./data/combine-data-set-a/all_ICUType.csv')
data_ICUType = as.numeric(data_ICUType[-1]) # 第一列为时间，删除
table(data_ICUType)
# 1: Coronary Care Unit
# 2: Cardiac Surgery Recovery Unit
# 3: Medical ICU
# 4: Surgical ICU


# 棘状图
df_a = read.csv('./data/Outcomes-a-Type.csv')
Age_range = c()
for (i in 1:4000) {
  if (Age[i, 1] < 40) {
    Age_range = c(Age_range, 1)
  }else if (Age[i, 1] >= 40 & Age[i, 1] < 60) {
    Age_range = c(Age_range, 2)
  }else if (Age[i, 1] >= 60 & Age[i, 1] < 80) {
    Age_range = c(Age_range, 3)
  }else if (Age[i, 1] >= 80) {
    Age_range = c(Age_range, 4) 
  }
}
df_a = cbind(df_a, Age_range)
count_mat = table(df_a$In.hospital_death, df_a$Age_range)
p_view_sex = spineplot(t(count_mat), 
                       col =c("#4784b7" ,"#9fcdff"),
                       xaxlabels = c('under 40','40-60', '60-80', 'over 80'),
                       yaxlabels = c('survivor', 'died in-hospital'))

















