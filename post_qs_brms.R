post_qs_brms <- function(model, n, k, post_q = .8, conf_int = .95, tvar = 1){  
  # 1. Obtain posterior draws of tau and sigma
  sd_model <- VarCorr(model, summary = FALSE)
  draws_tau <- sd_model[[1]][["sd"]][ , "Intercept"]
  draws_sigma <- sd_model$residual__$sd[ , 1]
  # 2. Compute draws for ICC
  draws_icc <- draws_tau^2 / (draws_tau^2 + draws_sigma^2)
  
  iccq <- as.double(quantile(draws_icc, post_q))
  sq_se <- ((n-1)*iccq + 1) * 4*tvar/(n*k)
  crit_q <- qnorm((1-conf_int)/2, lower.tail = FALSE)
  width <- 2 * crit_q * sqrt(sq_se)
  return(data.frame(icc_q = iccq, 
                    squared_se_q = sq_se,
                    ci_width_q = width))
}
