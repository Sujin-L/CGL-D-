library(gglasso)
source("cv_gglasso_mod.R")

# compositional group lasso (CGL)

# y: centered response variable
# Z: centered compositional covariate

CGL <- function(y, Z, M, r){
  
ind.g <- rep(1:M, r)

# constraint matrix 
C <- matrix(0, r, (r*M))
ind <- rep(1:r, each=M)
for(i in 1:r){
  C[i,  which(ind==i)] <- 1
}

weights <- c(rep(1, nrow(Z)), rep(1000, nrow(C)))
ind.gg <- rep(1:M, each = r)

# rearrange the order of variables
x <- rbind(Z, C)
x.mod <- NULL
for(i in 1:M){
  x.mod <- cbind(x.mod, x[,which(ind.g == i)])
}


fit.cvgg <- cv.gglasso.mod(x = x.mod, y = c(y, rep(0, nrow(C))), r, 
                           group = ind.gg, loss = "wls", weight = weights)
fit.gg <- gglasso(x = x.mod, y = c(y, rep(0, nrow(C))), group = ind.gg, 
                  loss = "wls", lambda = fit.cvgg$lambda.1se, 
                  weight = diag(weights, nrow(x.mod)), intercept = FALSE)

return(fit.gg$beta)

}


