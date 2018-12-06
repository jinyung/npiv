calc_power <- function(r, u) {
  fit <- lm(log10(u)~log10(r))  # power = slope of log-log model
  pval <- anova(fit)$P[1]
  ifelse(pval < 0.05, coef(fit)[2], NA)
}