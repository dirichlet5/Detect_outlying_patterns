rm(list = ls())


METHOD_LIST <- c("MMFOD", "JCGS_MFDPD", "JCGS_MFHD", "JCGS_MFPD", "JCGS_MFSPD")

### Plotting function
plot_dtt <- function(y, flag, grid_points, p, true_outliers, show_legend, plot_title,
                     title_cex, ylabel, xlabel, legend_pos = "topleft", labels) {
  ### -----###
  par(mai = c(1.5, 1, 1, 0.5))
  
  if (flag == 1) {
    lim <- c(0, 1.1)
  } else {
    lim <- range(y) + c(-.15 * sd(y[p, ]), .15 * sd(y[p, ]))
  }
  
  matplot(y,
          type = "o",
          # col = c("black", "blue", "red"),
          col = c("black", "blue", "red", "darkolivegreen3", "#af9a8b"),
          # col = c("black","blue","red","darkolivegreen3","#af9a8b","#594675"),
          lty = "solid",
          lwd = 2.5,
          # ylim = range(y) + c(-.15*sd(y[p,]), .15*sd(y[p,])),
          ylim = lim,
          ylab = ylabel,
          xlab = xlabel,
          axes = F, pch = 13:18,
          cex = 1.3, cex.lab = 2,
          col.lab = "gray20"
  )
  grid(col = "grey75", lwd = .3)
  axis(1,
       col = "white", at = c(2, 4, 6, 8), labels = labels,
       col.ticks = "grey61",
       lwd.ticks = .3, tck = -0.025,
       cex.axis = 1.5, cex.lab = 8, col.axis = "gray30"
  )
  # at=c(0,0.2,0.4,0.6,0.8,1),
  # labels=c(0,0.2,0.4,0.6,0.8,1),
  axis(2,
       col = "white",
       col.ticks = "grey61",
       lwd.ticks = .2, tck = -0.025,
       cex.axis = 1.5, cex.lab = 8,col.axis = "gray30"
  )
  
  box(col = "grey51")
  # ,"FB-MBD","ADA","OUG"
  # ,"darkolivegreen3","#af9a8b","#594675"
  if (show_legend) {
    legend(legend_pos,
           legend = c("SFOD-MM","Qu and Genton - MFDPD","- MFHD",
                      "- MFPD","- MFSPD"),
           lty = c("solid", "solid", "solid"),
           lwd = c(2, 2),
           col = c("black", "blue", "red", "darkolivegreen3", "#af9a8b"),
           text.col = "gray40", bty = "n",
           box.lwd = .1, pch = 13:18, xjust = -1, cex = 1, inset = .01
    )
  }
  mtext(plot_title, 3,
        adj = 0.5, line = 1, cex = title_cex,
        col = "gray20"
  )
}

# OUTLIER_LIST <- c(0.05, 0.1, 0.15)
outlier_rate <- 0.1

for (MODEL in 1:2) {
  for (TYPE in 1:3) {
    
    dev.off()
     if (TYPE == 1) {
      TPR_data <- data.frame(parameter = seq(2, 10, 1))
      FPR_data = TPR_data
      xlabel <- "zero_len"
      labels <- c(3, 5, 7, 9)
    } else if (TYPE == 2) {
      TPR_data <- data.frame(parameter = seq(6, 10, 0.5))
      FPR_data = TPR_data
      xlabel <- "k"
      labels <- c(6.5, 7.5, 8.5, 9.5)
    } else if (TYPE == 3) {
      TPR_data <- data.frame(parameter = seq(6, 10, 0.5))
      FPR_data = TPR_data
      xlabel <- "q"
      labels <- c(6.5, 7.5, 8.5, 9.5)
    }
    for (method in METHOD_LIST) {
      filepath <- paste0(
        "./ROC_result/", method, "/mod_", MODEL, "_TYPE_", TYPE,
        "_outlierRate_", outlier_rate, ".csv"
      )
      df <- read.csv(filepath)
      TPR_data <- cbind(TPR_data, df[2])
      FPR_data = cbind(FPR_data, df[3])
    }
    TPR_data <- TPR_data[, -1]
    FPR_data <- FPR_data[, -1]
    names(TPR_data) <- METHOD_LIST
    names(FPR_data) <- METHOD_LIST
    dev.new()
    par(mfrow=c(1,1))
    title <- ""
    
    plot_dtt(
      TPR_data, 1, tt1,
      p = 9,
      true_outliers = c(),
      show_legend = T, plot_title = title,
      title_cex = 1.5,
      ylabel = "TPR", xlabel = xlabel,
      labels = labels
    )
    
    
    plot_dtt(
      FPR_data, 2, tt1,
      p = 9,
      true_outliers = c(),
      show_legend = F, plot_title = title,
      title_cex = 1.5,
      ylabel = "FPR", xlabel = xlabel,
      labels = labels
    )
    
  }
}
