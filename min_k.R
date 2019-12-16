min_k <- function(icc, icc_sd, width = .35, post_q = .8, 
                  conf_int = .95, tvar = 1){
  alpha <- icc^2*(1-icc)/icc_sd^2 - icc
  beta <- alpha * (1/icc -1)
  iccq <- qbeta(post_q, alpha, beta)
  crit_q <- qnorm((1-conf_int)/2, lower.tail = FALSE)
  sq_se <- ((width/2)/crit_q)^2
  min_k <- 4*tvar*iccq / sq_se
  return(data.frame(min_k = ceiling(min_k)))
}
