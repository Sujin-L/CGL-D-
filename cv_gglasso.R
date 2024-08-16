library(gglasso)


cv.gglasso.mod <- function(x, y, r, group, pred.loss = "L1", 
                           nfolds = 5, foldid, delta=1, loss = loss, weights = weights) {
  if (missing(pred.loss)) 
    pred.loss <- "default" else pred.loss <- match.arg(pred.loss)
    
    N <- nrow(x)
    y <- drop(y)
    
    cv.x <- x[1:(N-r), ]; cv.y <- y[1:(N-r)]; cv.weights <- weights[1:(N-r)]
    weight.x <- tail(x, r); weight.y <- tail(y, r); w.weights <- tail(weights, r)
  
    ###Fit the model once to get dimensions etc of output
    gglasso.object <- gglasso(x, y, group, delta = delta, loss = loss, weight = diag(weights, N), intercept = FALSE)
    lambda <- gglasso.object$lambda
    
    # predict -> coef
    if (missing(foldid)) 
      foldid <- sample(rep(seq(nfolds), length = (N-r))) else nfolds <- max(foldid)
    if (nfolds < 3) 
      stop("nfolds must be bigger than 3; nfolds=10 recommended")
    outlist <- as.list(seq(nfolds))
    
    ###Now fit the nfold models and store them
    for (i in seq(nfolds)) {
      which <- foldid == i
      
      x_sub <- rbind(cv.x[!which, ], weight.x)
      y_sub <- c(cv.y[!which], weight.y)
      weight_sub <- c(cv.weights[!which], w.weights) 
        
      outlist[[i]] <- gglasso(x = x_sub, y = y_sub, group = group, 
                              lambda = lambda, delta = 1, loss = loss, weight = diag(weight_sub, nrow(x_sub)))
    }
    
    
    ###What to do depends on the pred.loss and the model fit
    fun <- paste("cv", class(gglasso.object)[[2]], sep = ".")
    cvstuff <- do.call(fun, list(outlist, lambda, x, y, foldid, pred.loss, delta=1))
    cvm <- cvstuff$cvm
    cvsd <- cvstuff$cvsd
    cvname <- cvstuff$name
    out <- list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvupper = cvm + cvsd, 
                cvlo = cvm - cvsd, name = cvname, gglasso.fit = gglasso.object)
    lamin <- getmin(lambda, cvm, cvsd)
    obj <- c(out, as.list(lamin))
    class(obj) <- "cv.gglasso"
    obj
}


getmin <- function(lambda, cvm, cvsd) {
  cvmin <- min(cvm)
  idmin <- cvm <= cvmin
  lambda.min <- max(lambda[idmin])
  idmin <- match(lambda.min, lambda)
  semin <- (cvm + cvsd)[idmin]
  idmin <- cvm <= semin
  lambda.1se <- max(lambda[idmin])
  list(lambda.min = lambda.min, lambda.1se = lambda.1se)
}



