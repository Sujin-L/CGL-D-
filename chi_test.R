#### chi-squared test ####

chi.test.fdr <- function(fit.sgg, X, y, Theta.hat, ind.gg, gr.ind, M){
  
  n <- nrow(X); p <- ncol(X)
  r <- p/M
  
  sigma.hat <- fit.sgg$hsigma
  Sigma.hat <- t(X) %*% X / n
  debias.var <- (sigma.hat^2) * Theta.hat %*% Sigma.hat %*% Theta.hat
  beta.debias <- fit.sgg$coefficients + Theta.hat %*% t(X) %*% (y - X %*% fit.sgg$coefficients) / n
  
  #### group inference ------------------------------------------------------------------------------------------------------------
  
  chi.test <- NULL
  for(j in 1:M){
    
    gr.ind <- which(ind.gg==j)
    betahat.Gj <- beta.debias[gr.ind]
    V.Gj <- debias.var[gr.ind, gr.ind]
    
    chisq.j <- n * t(betahat.Gj)%*% inv(V.Gj) %*% betahat.Gj
    chi.test <- c(chi.test, chisq.j)
  }
  
  names(chi.test) <- 1:M
  chi.result <- sort(chi.test[chi.test > qchisq(0.95, df = r)], decreasing = TRUE)
  
  p_value <- 1-pchisq(chi.test, df = r) 
  p_value_fdr <- p.adjust(p_value, method="BH")
  chi.fdr.result <- round(p_value_fdr[p_value_fdr < 0.05], 3)
  
  chi.result <- list(beta.debias = beta.debias, debias.var = debias.var,
                     p_value_fdr = p_value_fdr, chi.fdr.result = chi.fdr.result)
  
  return(chi.result)
  
}

