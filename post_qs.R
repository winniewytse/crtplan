# ICC, squared SE, and confidence interval widths for the treatment effect based on the
# posterior distribution of ICC, a beta distribution with parameters alpha and beta
# alpha = icc^2*(1-icc)/icc_sd^2 - icc
# beta = alpha * (1/icc -1)

post_qs <- function(icc, icc_sd, n, k, post_q = .8, conf_int = .95, tvar = 1){  
  alpha <- icc^2*(1-icc)/icc_sd^2 - icc
  beta <- alpha * (1/icc -1)
  iccq <- qbeta(post_q, alpha, beta)
  sq_se <- ((n-1)*iccq + 1) * 4*tvar/(n*k)
  crit_q <- qnorm((1-conf_int)/2, lower.tail = FALSE)
  width <- 2 * crit_q * sqrt(sq_se)
  return(data.frame(icc_q = iccq, 
                    squared_se_q = sq_se,
                    ci_width_q = width))
}
