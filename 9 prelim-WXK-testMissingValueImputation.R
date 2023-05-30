library(VIM)
library(magrittr)

dataset <- sleep[, c("Dream", "NonD", "BodyWgt", "Span")]
dataset$BodyWgt <- log(dataset$BodyWgt)
dataset$Span <- log(dataset$Span)
aggr(dataset)

imp_hotdeck <- hotdeck(dataset[,1:2], variable = "NonD")  # hotdeck imputation
imp_knn <- kNN(dataset[,1:2], variable = "NonD") # kNN imputation
imp_match <- matchImpute(dataset, variable = "NonD", match_var = c("BodyWgt","Span")) # match imputation
aggr(imp_knn, delimiter = "_imp")
aggr(imp_match, delimiter = "_imp")

library(Amelia)
data(freetrade)
freetrade1 <- freetrade[,1:4]
a.out <- amelia(freetrade1, m = 5, ts = "year", cs = "country", p2s = 0)
plot(a.out, which.vars = 3)





