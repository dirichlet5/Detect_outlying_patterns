library("refund")

month <- c(-18:42)
list("cd4")
CD4 <- cd4[sample(1:366,30,replace = F),] ## load the cd4 data which is a matrix of 366 * 61

#pdf("observed_cd4.pdf", height = 4, width = 4)
#par(
#  mfrow = c(1, 1), mai = c(0.5, 0.55, 0.3, 0.07),
#  mar = c(3.5, 3.5, 2, 1), mgp = c(2, 1, 0)
#)

par(mar=c(4,6,1,1))
plot(month, CD4[1, ],
  type = "n",
  ylim = range(CD4, na.rm = TRUE), ylab = "CD4 Cell Counts",
  xlab = "Months", main = "",cex.lab = 1.8,cex.axis = 1.8
)

# lty change the type of line
apply(CD4, 1, function(ij) {
  cd4_df <- data.frame(month, ij)
  lines(na.omit(cd4_df), col = "grey", lty = 4,lwd=4)
})

# pch change the shape of point
# cex change the size of point
apply(CD4, 1, function(ij) {
  cd4_df <- data.frame(month, ij)
  lines(month, ij, col = "black", lty = 1,lwd=3)
  points(month, ij, col = "black", cex = 2, pch = 19)
})

dev.off()

