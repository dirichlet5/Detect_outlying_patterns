
library(ggplot2)
library(patchwork)

rm(list=ls())


var_name = 'HR'
path = paste0('./data/combine-data-set-a/all_', var_name, '.csv')
df = read.csv(path)
df$Time = as.POSIXct(df$Time, format = '%Y-%m-%d %H:%M:%S')


START_TIME = 500
# START_TIME = nrow(df)

data = df[1:START_TIME,-1]
data_time = df[1:START_TIME,1]

# Lt and Ly record the location and value of outliers
Lt = list()
Ly = list()
for (i in 1:ncol(data)) {
  
  Lt[[i]] = as.numeric(which(!is.na(data[i])))
  Ly[[i]] = as.numeric(data[which(!is.na(data[i])),i])
  
}


SAMPLE = 96
plot_data = cbind(as.data.frame(data_time), df[1:START_TIME, SAMPLE])
colnames(plot_data)[2] = colnames(df)[SAMPLE]

which(colnames(df) == 'X132610')
# Increase the y of spline interpolation
interpolation = spline(1:START_TIME,as.numeric(plot_data[,2]), n = START_TIME)
plot_data = cbind(plot_data, Y = interpolation$y) 

# Increase the time of spline interpolation

X = as.POSIXct(-8*60*60 + interpolation$x, origin = '1900-01-01 00:00:00')
plot_data = cbind(plot_data, X) 

difftime(data_time[100], 
         as.POSIXct(-8*60*60, origin = '1900-01-01 00:00:00'), 
         units = 'secs')

# No sampling
p = ggplot(plot_data) +
  geom_point(aes(x = data_time, y = plot_data[,2]),color = '#0000fd', size = 4) +
  geom_line(aes(x = X, y = Y), linetype = 2, size = 1) +
  scale_x_datetime(date_breaks = '2 min',date_labels = "%M:%S") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1) ) + 
  # ggtitle('Univariate Sampled') + 
  ylab('X(t)') +
  xlab('Time') +
  theme(plot.title = element_text(hjust = 0.5, size = 20)) + 
  theme(axis.text.x = element_text(size = 30, face = "bold")) + 
  theme(axis.text.y = element_text(size = 30, face = "bold")) + 
  theme(axis.title.x = element_text(size = 30, face = "bold")) + 
  theme(axis.title.y = element_text(size = 30, face = "bold")) + 
  theme(panel.background = element_blank()) +
  theme(panel.border = element_blank(), axis.line = element_line()) + 
  theme(plot.margin=unit(rep(1,4),'lines'))
p



# -----------multivariate irregulary sampled(unaligened)----------
for (i in 2:4) {
  SAMPLE = i
  plot_data = cbind(as.data.frame(data_time), df[1:START_TIME, SAMPLE])
  colnames(plot_data)[2] = colnames(df)[SAMPLE]

  # Increase the y of spline interpolation
  interpolation = spline(1:START_TIME,as.numeric(plot_data[,2]), n = START_TIME)
  plot_data = cbind(plot_data, Y = interpolation$y)

  # Increase the time of spline interpolation

  X = as.POSIXct(-8*60*60 + interpolation$x, origin = '1900-01-01 00:00:00')
  plot_data = cbind(plot_data, X)

  difftime(data_time[100],
           as.POSIXct(-8*60*60, origin = '1900-01-01 00:00:00'),
           units = 'secs')

  # No sampling
  p = ggplot(plot_data) +
    geom_point(aes(x = data_time, y = plot_data[,2]),color = '#0000fd', size = 4) +
    geom_line(aes(x = X, y = Y), linetype = 2, size = 1) +
    scale_x_datetime(date_breaks = '2 min',date_labels = "%M:%S") +
    theme(axis.text.x = element_text(angle = 45,hjust = 1) ) +
    # ggtitle('Univariate Sampled') +
    ylab('X(t)') +
    xlab('Time') +
    theme(plot.title = element_text(hjust = 0.5, size = 20)) +
    theme(axis.text.x = element_text(size = 30, face = "bold")) +
    theme(axis.text.y = element_text(size = 30, face = "bold")) +
    theme(axis.title.x = element_text(size = 30, face = "bold")) +
    theme(axis.title.y = element_text(size = 30, face = "bold")) +
    theme(panel.background = element_blank()) +
    theme(panel.border = element_blank(), axis.line = element_line())

  fig_name = paste0('./figure/Sparse_plot/',var_name,'_',colnames(plot_data)[2],'.png')
  ggsave(fig_name, plot=p)
}


## ------ggplot function-----------
sparse_plot = function(sample_id){
  SAMPLE = sample_id
  plot_data = cbind(as.data.frame(data_time), df[1:START_TIME, SAMPLE])
  colnames(plot_data)[2] = colnames(df)[SAMPLE]
  
  which(colnames(df) == 'X132610')
  # Increase the y of spline interpolation
  interpolation = spline(1:START_TIME,as.numeric(plot_data[,2]), n = START_TIME)
  plot_data = cbind(plot_data, Y = interpolation$y) 
  
  # Increase the time of spline interpolation
  
  X = as.POSIXct(-8*60*60 + interpolation$x, origin = '1900-01-01 00:00:00')
  plot_data = cbind(plot_data, X) 
  
  difftime(data_time[100], 
           as.POSIXct(-8*60*60, origin = '1900-01-01 00:00:00'), 
           units = 'secs')
  
  # No sampling
  p = ggplot(plot_data) +
    geom_point(aes(x = data_time, y = plot_data[,2]),color = '#0000fd', size = 4) +
    geom_line(aes(x = X, y = Y), linetype = 2, size = 1) +
    scale_x_datetime(date_breaks = '2 min',date_labels = "%M:%S") +
    theme(axis.text.x = element_text(angle = 45,hjust = 1) ) + 
    # ggtitle('Univariate Sampled') + 
    ylab('X(t)') +
    xlab('Time') +
    theme(plot.title = element_text(hjust = 0.5, size = 20)) + 
    theme(axis.text.x = element_text(size = 30, face = "bold")) + 
    theme(axis.text.y = element_text(size = 30, face = "bold")) + 
    theme(axis.title.x = element_text(size = 30, face = "bold")) + 
    theme(axis.title.y = element_text(size = 30, face = "bold")) + 
    theme(panel.background = element_blank()) +
    theme(panel.border = element_blank(), axis.line = element_line()) + 
    theme(plot.margin=unit(rep(1,4),'lines'))
  p
  
}



## -----output result to pdf-------------
id_list = c('X132610', 'X132682', 'X132780', 'X133125', 'X133086', 'X132835')
pos_list = c()
for (i in id_list) {
  pos = which(colnames(df) == i)
  pos_list = c(pos_list, pos)
}
print(pos_list)
# write in pdf
pdf(file = "Sparse_data.pdf", width = 19, height = 15)

sparse_plot(32) + sparse_plot(60) + 
sparse_plot(96) + sparse_plot(218) + 
sparse_plot(207) + sparse_plot(119) + 
plot_layout(ncol = 2)

dev.off()




