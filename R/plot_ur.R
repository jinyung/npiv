# --- plot U* against r ---
# result: (list) object from calc_r() function
# ... pass arguments into matplot()
# col: col for each time frame, default is rainbow color if not specified
plot_ur <- function(result, col, ...) {
  
  # prevent bug from log with zeros
  if (any(result$r == 0)) {
    result$r[result$r == 0] <- NA
    warning("zero values in radius are NA-ed")
  }
  
  # calculate b
  b <- apply(result$r, 2, calc_power, result$u)
  
  # range of power b
  b1 <- floor(min(b, na.rm = TRUE))
  b2 <- ceiling(min(b, na.rm = TRUE))
  
  # calculate intercept
  a1 <- mean(log10(result$u) - b1 * log10(result$r[, which.min(b)]))
  a2 <- mean(log10(result$u) - b2 * log10(result$r[, which.min(b)]))
  
  # determine range of plot
  ymin.power <- min(floor(log10(result$u)))
  ymax.power <- max(ceiling(log10(result$u)))
  xmin.power <- floor(log10(result$r))
  xmin.power <- min(xmin.power[is.finite(xmin.power)])
  xmax.power <- ceiling(log10(result$r))
  xmax.power <- max(xmax.power[is.finite(xmax.power)])
  
  if(missing(col))
    col = rainbow(dim(result$r)[2])
                  
  # plot
  matplot(result$r, result$u, log = "xy", pch = 16,
          col = col, xaxt = "n", yaxt = "n",
          ylim = c(10^(ymin.power), 10^(ymax.power)),
          xlim = c(10^(xmin.power), 10^(xmax.power)), ...)
  yticks.power <- seq(ymin.power, ymax.power, 1)
  axis(side = 2, at = 10^(yticks.power),
       sapply(yticks.power, function(x) parse(text = paste0("10^", x))),
       las = 2)
  xticks.power <- seq(xmin.power, xmax.power, 1)
  axis(side = 1, at = 10^(xticks.power),
       sapply(xticks.power, function(x) parse(text = paste0("10^", x))))
  
  # show slopes
  abline(a = a1, b = b1, col = 2)
  abline(a = a2, b = b2, col = 3)
  
  # legend
  legend("bottomleft", lty = 1, col = c(2, 3), bty = "n",
         legend = c(parse(text = paste0("r^", b1)),
                    parse(text = paste0("r^", b2))))
  
}
