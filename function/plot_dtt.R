
plot_dtt <- function(y, grid_points, p, true_outliers, show_legend, plot_title,
                     title_cex, ylabel, xlabel, legend_pos = "bottomright") {
  if (length(true_outliers) > 0) {
    dttout <- y[, true_outliers, drop = F]
    dttnorm <- y[, -true_outliers, drop = F]
    plot(
      x = grid_points, type = "n", ylab = ylabel, xlab = xlabel,
      ylim = range(y) + c(-0.5 * sd(y[p, ]), 0.5 * sd(y[p, ])), col.lab = "gray20", axes = F
    )
    grid(col = "grey75", lwd = 0.3)
    matlines(dttnorm, col = "grey61", lty = "solid", lwd = 0.4)
    matlines(dttout, col = "#D55E00", lty = "solid", lwd = 1.3)
  } else {
    matplot(y,
      type = "l", col = "grey61", lty = "solid",
      lwd = 0.4, ylim = range(y) + c(
        -0.5 * sd(y[p, ]),
        0.5 * sd(y[p, ])
      ), ylab = ylabel, xlab = xlabel,
      axes = F, col.lab = "gray20"
    )
    grid(col = "grey75", lwd = 0.3)
  }
  axis(1,
    col = "white", col.ticks = "grey61", lwd.ticks = 0.5,
    tck = -0.025, cex.axis = 0.9, col.axis = "gray30"
  )
  axis(2,
    col = "white", col.ticks = "grey61", lwd.ticks = 0.5,
    tck = -0.025, cex.axis = 0.9, col.axis = "gray30"
  )
  box(col = "grey51")
  if (show_legend) {
    legend(legend_pos,
      legend = c("normal", "outlier"),
      lty = c("solid", "solid"), lwd = c(0.4, 1.3), col = c(
        "grey61",
        "#D55E00"
      ), text.col = "gray40", bty = "n",
      box.lwd = 0.1, xjust = 0, inset = 0.01
    )
  }
  mtext(plot_title, 3,
    adj = 0.5, line = 1, cex = title_cex,
    col = "gray20"
  )
}
