# Posterior distributions of ICC, squared SE, and confidence interval width for the treatment effect
# Compute the squared SE and C.I. width for the specified ICC mean, sd, cluster size (n), 
# and number of clusters

post_dists_brms <- function(model, n, k, post_q = .8, conf_int = .95, tvar = 1){  
  # 1. Obtain posterior draws of tau and sigma
  sd_model <- VarCorr(model, summary = FALSE)
  draws_tau <- sd_model[[1]][["sd"]][ , "Intercept"]
  draws_sigma <- sd_model$residual__$sd[ , 1]
  # 2. Compute draws for ICC
  draws_icc <- draws_tau^2 / (draws_tau^2 + draws_sigma^2)
  
  iccq <- as.double(quantile(draws_icc, post_q))
  sq_se <- ((n-1)*iccq + 1) * 4*tvar/(n*k)
  crit_q <- qnorm((1-conf_int)/2, lower.tail = FALSE)
  w <- 2 * crit_q * sqrt(sq_se)
  
  # plot the posterior distribution of ICC
  dens1 <- density(draws_icc)
  df1 <- data.frame(x = dens1$x, y = dens1$y)
  exp1 <- bquote(rho[.(post_q)] == .(round(iccq, 4)))
  plot1 <- ggplot(df1, aes(x, y)) + 
    geom_line() +
    geom_area(mapping = aes(ifelse(x < iccq, x, 0)), fill = "turquoise3", alpha = .5) +
    geom_vline(xintercept = iccq, linetype = "dashed", color = "turquoise3", alpha = .5, size = 1) +
    geom_text(aes(x = iccq, y = max(df1$y)), 
              label = deparse(exp1), parse = TRUE, color = "turquoise4", size = 4) +
    ylim(0, max(df1$y) + .5) +
    labs(title = "Posterior Distribution of ICC", x = expression(rho), y = "Density")
  
  # plot the posterior distribution of squared SE
  sq_se_dis <- ((n-1)*draws_icc + 1) * 4*tvar/(n*k)
  dens2 <- density(sq_se_dis)
  df2 <- data.frame(x = dens2$x, y = dens2$y)
  exp2 <- bquote(SE[.(post_q)]^2 == .(round(sq_se, 4)))
  plot2 <- ggplot(df2, aes(x, y)) + 
    geom_line() +
    geom_area(mapping = aes(ifelse(x < sq_se, x, 0)), fill = "coral2", alpha = .5) +
    geom_vline(xintercept = sq_se, linetype = "dashed", color = "coral2", alpha = .5, size = 1) +
    geom_text(aes(x = sq_se, y = max(df2$y)), 
              label = deparse(exp2), parse = TRUE, color = "coral3", size = 4) +
    ylim(0, max(df2$y) + .5) +
    labs(title = expression(Posterior~Distribution~of~SE^2), 
         x = expression(SE^2), y = "Density")
  
  # plot the posterior distirbution of width
  w_dis <- 2 * crit_q * sqrt(sq_se_dis)
  dens3 <- density(w_dis)
  df3 <- data.frame(x = dens3$x, y = dens3$y)
  plot3 <- ggplot(df3, aes(x, y)) +
    geom_line() +
    geom_area(mapping = aes(ifelse(x < w, x, 0)), fill = "purple2", alpha = .4) +
    geom_vline(xintercept = w, linetype = "dashed", color = "purple2", alpha = .4, size = 1) +
    geom_text(aes(x = w, y = max(df3$y)), label = deparse(exp3), parse = TRUE,  
              color = "purple4", size = 4) +
    ylim(0, max(df3$y) + .5) +
    labs(title = "Posterior Distribution of C.I. Width", x = "Width", y = "Density")
  
  allplots <- grid.arrange(plot1, plot2, plot3, nrow=2)
  allplots
  cat(paste("The squared standard error of the treatment effect will be", 
            round(sq_se, 4), "for the specified mean, sd, and credible interval of ICC.\n"))
  return(data.frame(icc_q = iccq, 
                    squared_se_q = sq_se))
}
