min_k_brms <- function(model, width = .35, post_q = .8, 
                       conf_int = .95, tvar = 1){
  # 1. Obtain posterior draws of tau and sigma
  sd_model <- VarCorr(model, summary = FALSE)
  draws_tau <- sd_model[[1]][["sd"]][ , "Intercept"]
  draws_sigma <- sd_model$residual__$sd[ , 1]
  # 2. Compute draws for ICC
  draws_icc <- draws_tau^2 / (draws_tau^2 + draws_sigma^2)
  
  iccq <- as.double(quantile(draws_icc, post_q))
  crit_q <- qnorm((1-conf_int)/2, lower.tail = FALSE)
  sq_se <- ((width/2)/crit_q)^2
  min_k <- 4*tvar*iccq / sq_se
  return(data.frame(min_k = ceiling(min_k)))
}
