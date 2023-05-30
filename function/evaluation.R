
evaluation = function(n, results, true){
  tp = length(intersect(results, true))/length(true)
  fp = length(setdiff(results, true))/(n - length(true))
  list(tp = tp, fp = fp)
}

BCR = function(n, results, true){
  TP = length(intersect(results, true))/length(true)
  
  
  bcr = 0.5 * ((TP / (TP + FN)) +(TN / (TN + FP)))
  return(bcr)
}














