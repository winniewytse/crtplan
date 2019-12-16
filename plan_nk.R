plan_nk <- function(icc, icc_sd, n = NULL, k = NULL, width = .35, post_q = .8, 
                    conf_int = .95, tvar = 1){
  alpha <- icc^2*(1-icc)/icc_sd^2 - icc
  beta <- alpha * (1/icc -1)
  iccq <- qbeta(post_q, alpha, beta)
  crit_q <- qnorm((1-conf_int)/2, lower.tail = FALSE)
  sq_se <- ((width/2)/crit_q)^2
  min_k <- 4*tvar*iccq / sq_se
  min_k_int <- ceiling(min_k)
  
  if(hasArg(k)){
    if(min_k > k){
      warning("A minimum of ", ceiling(min_k), " clusters is required for the desired confidence interval width. ", 
              "Increase k or width. ")
    } else {
      # Plot CI width vs n for varying k
      k1 <- post_qs(icc, icc_sd, n = 1:500, k = min_k_int, post_q = post_q)$ci_width_q
      k2 <- post_qs(icc, icc_sd, n = 1:500, k = ceiling(min_k_int*1.5), post_q = post_q)$ci_width_q
      k3 <- post_qs(icc, icc_sd, n = 1:500, k = min_k_int*2, post_q = post_q)$ci_width_q
      df1 <- data.frame(n = rep(1:500, 3), 
                        se = c(k1, k2, k3), 
                        k = rep(c("k1", "k2", "k3"), each = 500))
      plot1 <- ggplot(df1, aes(x=n, y=se, col=k)) + 
        geom_line() +
        labs(title = "C.I. Width vs n", x = "n", y = "C.I. Width") +
        scale_color_discrete(name = "Number of Cluters (k)",
                             labels = c(min_k_int, ceiling(min_k_int*1.5), min_k_int*2)) +
        geom_hline(yintercept = width) + 
        geom_text(aes(x = 200, y = width), label = paste0("desired width = ", width), 
                  color = "black", size=3, vjust=-0.4, hjust=0)  +
        theme(legend.position = c(0.7, 0.8))
      
      # Plot CI width vs n for varying post_q
      credq50 <- post_qs(icc, icc_sd, n = 1:500, k = k, post_q = .5)$ci_width_q
      credq80 <- post_qs(icc, icc_sd, n = 1:500, k = k, post_q = .8)$ci_width_q
      credq95 <- post_qs(icc, icc_sd, n = 1:500, k = k, post_q = .95)$ci_width_q
      df2 <- data.frame(n = rep(1:500, 3), 
                        credq = c(credq50, credq80, credq95), 
                        conf = rep(c("credq50", "credq80", "credq95"), each = 500))
      plot2 <- ggplot(df2, aes(x=n, y=credq, col=conf)) + 
        geom_line() +
        labs(x = "n", y = "C.I. Width", 
             title = paste0("For k = ", ceiling(k), ", C.I. Width against n")) +
        scale_color_discrete(name = paste0("% of chance an observed \nC.I. width < ", width), 
                             labels = c("50%", "80%", "95%")) +
        geom_hline(yintercept = width) + 
        geom_text(aes(x = 200, y = width), label = paste0("desired width = ", width), color = "black",
                  size=3, vjust=-0.4, hjust=0)  +
        theme(legend.position = c(0.7, 0.8))
      allplots <- grid.arrange(plot1, plot2, ncol=2)
      
      n <- 4*tvar*(1-iccq) / (sq_se*k - 4*tvar*iccq)
      allplots
      return(data.frame(n = ceiling(n)))
    }
  } else {
    k <- ((n-1)*iccq + 1) * 4*tvar/(n*sq_se)
    
    # Plot CI width vs n for varying k
    k1 <- post_qs(icc, icc_sd, n = 1:500, k = ceiling(k/1.5), post_q = post_q)$ci_width_q
    k2 <- post_qs(icc, icc_sd, n = 1:500, k = k, post_q = post_q)$ci_width_q
    k3 <- post_qs(icc, icc_sd, n = 1:500, k = ceiling(k*1.5), post_q = post_q)$ci_width_q
    df1 <- data.frame(n = rep(1:500, 3), 
                      se = c(k1, k2, k3), 
                      k = rep(c("k1", "k2", "k3"), each = 500))
    plot1 <- ggplot(df1, aes(x=n, y=se, col=k)) + 
      geom_line() +
      labs(title = "C.I. Width against n", x = "n", y = "C.I. Width") +
      scale_color_discrete(name = "Number of Cluters (k)",
                           labels = c(min_k_int, ceiling(min_k_int*1.5), min_k_int*2)) +
      geom_hline(yintercept = width) + 
      geom_text(aes(x = 200, y = width), label = paste0("desired width = ", width), 
                color = "black", size=3, vjust=-0.4, hjust=0)  +
      theme(legend.position = c(0.7, 0.8))
    
    # Plot CI width vs n for varying post_q
    credq50 <- post_qs(icc, icc_sd, n = 1:500, k = k, post_q = .5)$ci_width_q
    credq80 <- post_qs(icc, icc_sd, n = 1:500, k = k, post_q = .8)$ci_width_q
    credq95 <- post_qs(icc, icc_sd, n = 1:500, k = k, post_q = .95)$ci_width_q
    df2 <- data.frame(n = rep(1:500, 3), 
                      credq = c(credq50, credq80, credq95), 
                      conf = rep(c("credq50", "credq80", "credq95"), each = 500))
    plot2 <- ggplot(df2, aes(x=n, y=credq, col=conf)) + 
      geom_line() +
      labs(x = "n", y = "C.I. Width", 
           title = paste0("For k = ", ceiling(k), ", C.I. Width against n")) +
      scale_color_discrete(name = paste0("Level of certainty that \nC.I. width < ", width), 
                           labels = c("50%", "80%", "95%")) +
      geom_hline(yintercept = width) + 
      geom_text(aes(x = 200, y = width), label = paste0("desired width = ", width), color = "black",
                size=3, vjust=-0.4, hjust=0)  +
      theme(legend.position = c(0.7, 0.8))
    allplots <- grid.arrange(plot1, plot2, ncol=2)
    
    allplots
    return(data.frame(k = ceiling(k)))
  }
}
