library("ggplot2")
library("rlist")
library("MASS")
library("tidyverse")
library("plotly")

rm(list = ls())

# Acquisition of variable name
var_name <- read.csv("./data/var_name2.txt", header = F)
print(var_name)

## ------------plot---------------------
variable <- "HR"
# variable <- "HCT"
# variable <- "Glucose"
# variable <- "Na"
# variable <- "NIDiasABP"



path <- paste0("./data/combine-data-set-a/all_", variable, ".csv")
df <- read.csv(path)
var_ID <- as.numeric(gsub("X", "", colnames(df)[-1]))
time <- as.data.frame(as.POSIXct(df[, 1]))

NUM <- 3

if (variable == "HR") {
  sample_pos <- c(245, 589, 267) # choose, HR outlier
  sample_pos <- c(24, 128, 358) # choose, HR outlier
}else if (variable == "Na") {
  sample_pos <- c(744, 795, 74) # choose, HCT outlier
}else if (variable == "Glucose") {
  sample_pos <- c(193, 399, 41) # choose, HCT outlier
}else if (variable == "Temp") {
  sample_pos <- c(62, 413, 177) # choose, HCT outlier
}else if (variable == "NIDiasABP") {
  sample_pos <- c(56, 21, 138)
  sample_pos <- c(20, 56, 156)
  sample_pos <- c(20, 56, 130)
}


sample_ID <- var_ID[sample_pos]
print(sample_ID)

sample_df <- cbind(time, df[, 1 + sample_pos])
colnames(sample_df)[1] <- "time"

df_mean <- rowMeans(df[, -1], na.rm = T)
df_plot_line <- cbind(time, df_mean)

MEAN_SPAR <- 35
spar_df <- df_plot_line[seq(1, length(df_mean), MEAN_SPAR), ]
names(spar_df) <- c("time", "mean")

df_plot <- left_join(sample_df, spar_df, by = "time")
for (i in 1:NUM) {
  interpolation <- spline(1:nrow(df_plot), as.numeric(df_plot[, i + 1]), n = nrow(df_plot))
  X <- as.POSIXct(-8 * 60 * 60 + interpolation$x, origin = "1900-01-01 00:00:00")
  df_plot <- cbind(df_plot, X)
  df_plot <- cbind(df_plot, Y = interpolation$y)
  colnames(df_plot)[2 * i + 4] <- paste0(sample_ID[i], "_time")
  colnames(df_plot)[2 * i + 5] <- sample_ID[i]
}
interpolation <- spline(1:nrow(df_plot), as.numeric(df_plot[, 5]), n = nrow(df_plot))
X <- as.POSIXct(-8 * 60 * 60 + interpolation$x, origin = "1900-01-01 00:00:00")
df_plot <- cbind(df_plot, mean_time = X)
df_plot <- cbind(df_plot, mean_smooth = interpolation$y)

p <- ggplot(df_plot) +
  geom_line(aes(x = df_plot[, 6], y = df_plot[, 7]),
            color = "grey", linetype = "longdash", size = 1
  ) +
  geom_line(aes(x = df_plot[, 8], y = df_plot[, 9]),
            color = "grey", linetype = "longdash", size = 1
  ) +
  geom_line(aes(x = df_plot[, 10], y = df_plot[, 11]),
            color = "grey", linetype = "longdash", size = 1
  ) +
  geom_line(aes(x = df_plot[, 12], y = df_plot[, 13]),
            color = "black", linetype = "solid", size = 1
  ) +
  geom_point(aes(x = time, y = df_plot[, 2]),
             color = "#EE052E", size = 2.5
  ) +
  geom_point(aes(x = time, y = df_plot[, 3]),
             color = "#009400", size = 2.5
  ) +
  geom_point(aes(x = time, y = df_plot[, 4]),
             color = "#0083FF", size = 2.5
  ) +
  scale_x_datetime(date_breaks = "8 min", date_labels = "%M:%S") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab(variable) +
  xlab("Time") +
  theme(plot.title = element_text(hjust = 0.5, size = 20)) +
  theme(axis.text.x = element_text(size = 30, face = "bold")) +
  theme(axis.text.y = element_text(size = 30, face = "bold")) +
  theme(axis.title.x = element_text(size = 30, face = "bold")) +
  theme(axis.title.y = element_text(size = 30, face = "bold")) +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, 
                                    linetype="solid"), axis.line = element_line()) +
  theme(plot.margin = unit(rep(1, 4), "lines"))
p





