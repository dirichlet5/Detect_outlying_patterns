rm(list=ls())

################################
#Load the packages
################################
list.of.packages <- c("gratis","feasts","dplyr","ggplot2","tsibble",
                      "ggfortify","tsfeatures")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(gratis)
library(feasts)
library(dplyr)
library(ggplot2)
library(tsibble)
library(ggfortify)
library(tsfeatures)


################################
#Read in a time series, simulate another one using specified features.
################################

#Read in the data
var_name = 'HR'
path = paste0('data/combine-data-set-a/all_', var_name, '.csv')
df = read.csv(path)
mpg <- na.omit(data_mean1)
mpg = mpg[1000:1300]  #<------put here the time series, in a vector form. This can not contain NA's.
mpg1 <- ts(mpg)
sum(is.na(mpg1))
plot(mpg1)


#This specifies what features we want to simulate. We used 4 features as an example.
# Here is the link for the available features:
# https://cran.r-project.org/web/packages/tsfeatures/vignettes/tsfeatures.html
my_features1 <- function(y) {
  c(stl_features(y)[c( "spike","e_acf1","e_acf10")], lumpiness(y)[c("lumpiness")])
}
# This is the simulation code, I found this to be rather slow for long time series.
y <- simulate_target(
  length = length(mpg1),
  seasonal_periods = frequency(mpg1),
  feature_function = my_features1, target = my_features1(mpg1)
)

# 

# Make new series same scale and frequency as the benchmark data.
y <- ts(scale(y) * sd(mpg1) + mean(mpg1))
tsp(y) <- tsp(mpg1)
cbind(as.ts(mpg1), as.ts(y)) %>% autoplot()
cbind(mpg1, y) %>% autoplot() # Here "y" is the final data that we generated. This is in "ts" format.





