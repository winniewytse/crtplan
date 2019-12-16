geticc <- function(model){
  # 1. Obtain posterior draws of tau and sigma
  sd_model <- VarCorr(model, summary = FALSE)
  draws_tau <- sd_model[[1]][["sd"]][ , "Intercept"]
  draws_sigma <- sd_model$residual__$sd[ , 1]
  # 2. Compute draws for ICC
  draws_icc <- draws_tau^2 / (draws_tau^2 + draws_sigma^2)
  # Summarize the ICC distribution
  icc <- mean(draws_icc)
  icc_sd <- sd(draws_icc)
  # Plot the ICC
  plot <- qplot(draws_icc_m1, geom = "density", xlab = "ICC") +
    labs(title = "Posterior Distribution of ICC", y = "Density")
  print(plot)
  return(data.frame(icc = icc, 
                    icc_sd = icc_sd))
}
