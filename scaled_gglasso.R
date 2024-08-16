#### scaled group lasso ####
require(gglasso)
source("cv_gglasso.R")


slassoEst_gg <- function(X, y, M, weights, group, lam0 = "univ"){
  
  nX = dim(X)[1]; pX = dim(X)[2]
  r = pX / M
  

  Sigma.hat <- t(X) %*% X / nX
  trace <- spec <- NULL
  
  for(i in 1:M){
    ind.i <- which(group==i)
    Sigma.j <- Sigma.hat[ind.i, ind.i]
    trace <- c(trace, sum(diag(Sigma.j)))
    spec <- c(spec, norm(Sigma.j, "2"))
    
  }

  q=0; b=0
  while((abs(1-b)>1e-30) & (( 2*sqrt(max(trace) +2*max(spec)*(2*q*log(M) + sqrt(r*q*log(M)))) / sqrt(nX)) < 1)){
    q <- q+1
    b <- 1-2*M^(1-q)
  }
  
  if((( 2*sqrt(max(trace) +2*max(spec)*(2*q*log(M) + sqrt(r*q*log(M)))) / sqrt(nX)) >1)){
    q <- q-1
  }
  if(q < 0){
    q = 0
  }

  if (lam0 == "univ" | lam0 == "universal"){
    lam0 <- 2*sqrt(max(trace) +2*max(spec)*(2*q*log(M) + sqrt(r*q*log(M))))/sqrt(nX)
  } 
  
  
  
  #### change objective function to group lasso (gglasso) #### 
  objlasso <- cv.gglasso.mod(x = X, y = y, group = group, r = r, loss = "wls", weight = weights)
  
  sigmaint = 0.1
  sigmanew = 5
  flag = 0
  
  while (abs(sigmaint - sigmanew) > 1e-10 & flag <= 100) {
    flag = flag + 1
    sigmaint = sigmanew
    lam = lam0 * sigmaint
    hy = X %*% as.vector(gglasso(x = X, y = y, group = group, loss = "wls", weight = diag(weights, nX), intercept = FALSE, lambda = lam)$beta) 
    sigmanew = sqrt(mean((y - hy)^2))
  }


  hsigma = sigmanew
  hlam = lam
  hbeta = as.vector(gglasso(x = X, y = y, group = group, loss = "wls", weight = diag(weights, nX), intercept = FALSE, lambda = lam)$beta)
  hy = X %*% hbeta  
  
  return(list(hsigma = hsigma, coefficients = hbeta, residuals = y - hy, fitted.values = hy, lambda.seq = objlasso$lambda, lambda = hlam, lam0 = lam0))
  
}



#### scaled lasso inversion 
slassoInv_gg <- function(X, y, weights, group, lam0 = NULL, LSE = F) {
  
  nX = dim(X)[1]
  pX = dim(X)[2]
  
  hsigma = rep(0, pX)
  Beta = matrix(-1, pX, pX)
  res = matrix(0, nX, pX)
  
  if (LSE == T) {
    hsigma.lse = rep(0, pX)
    Beta.lse = matrix(-1, pX, pX)
    res.lse = matrix(0, nX, pX)
  }
  
  for (j in 1:pX) {
    scalefac = sqrt(colSums(X[, -j]^2)/nX)
    Xj = t(t(X[, -j])/scalefac)
    
    objtmp = slassoEst_gg(X = Xj, y = X[, j], weights, group, lam0)
    
    hsigma[j] = objtmp$hsigma
    Beta[-j, j] = objtmp$coefficients/scalefac
    res[, j] = objtmp$residuals
    
    if (LSE == T) {
      lsetmp = lse(Xj, X[, j], indexset = which(objtmp$coefficients != 0))
      hsigma.lse[j] = lsetmp$hsigma
      Beta.lse[-j, j] = lsetmp$coefficients/scalefac
      res.lse[, j] = lsetmp$residuals
    }
  }
  
  tTheta = diag(hsigma^(-2))
  tTheta = -Beta %*% tTheta
  hTheta = tTheta * (abs(tTheta) <= abs(t(tTheta))) + t(tTheta) * (abs(tTheta) > abs(t(tTheta)))
  est = list(precision = hTheta, hsigma = hsigma)
  
  if (LSE == T) {
    tTheta.lse = diag(hsigma.lse^(-2))
    tTheta.lse = -Beta.lse %*% tTheta.lse
    hTheta.lse = tTheta.lse * (abs(tTheta.lse) <= abs(t(tTheta.lse))) + 
      t(tTheta.lse) * (abs(tTheta.lse) > abs(t(tTheta.lse)))
    lse = list(precision = hTheta.lse, hsigma = hsigma.lse)
    est$lse = lse
  }
  
  est
  
}

scalreg_gg <- function(X, y, M, weights, group, lam0 = "univ", LSE = FALSE) {
  

  if (!is.null(y)) {
    
    est <- slassoEst_gg(X, y, M, weights, group, lam0 = lam0)
    
    est$fitted.values <- as.vector(X %*% est$coefficients)
    est$residuals <- y - est$fitted.values
    est$type = "regression"
    
    if (LSE == TRUE) {
      lse = lse(X, y, indexset = which(est$coefficients != 0))
      est$lse = lse
    }
  }
  if (is.null(y)) {
    est <- slassoInv_gg(X, y, C.weight, group, lam0, LSE)
    est$type = "precision matrix"
  }
  est$call <- match.call()
  class(est) <- "scalreg"
  return(est)
}

