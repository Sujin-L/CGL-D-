
#### Estimation of precision matrix ####

library(CVXR)
est.precision <- function(X, P_C, lambda, group, d=0.01){
  
  n <- dim(X)[1]; p <- dim(X)[2]
  num.gr <- length(unique(group)); r <- p/num.gr
  
  Q <- (diag(1, nrow = p) - P_C)
  Sigma.hat <- t(X) %*% X / n
  
  # tuning parameter
  gamma = d*lambda
  b1 = rep(1, p)*gamma
  
  M <- matrix(NA, p, p)
  
  for(i in 1:p){
    
    b2 <- Q[, i]

      # find m with initial gamma
      m <- Variable(p)
      obj <- Minimize(quad_form(m, Sigma.hat)) ## m^T Sigma m
      
      constr <- list()
      for(j in 1:num.gr){
        lg <- (group == j)
        constr.g <- list(norm((Sigma.hat %*% m - b2)[lg,] , type="2") <= b1)
        constr <- append(constr, constr.g)
      }
      
      prob <- Problem(obj, constr)
      result <- solve(prob, solver="SCS", verbose = FALSE)
      
      M[i,] <- result$getValue(m)


  # if optimization problem is unsolved, update gamma (gamma2)
      if(result$status != "optimal"){

        gamma2 <- Variable(1)
        m <- Variable(p)

        obj2 <- Minimize(gamma2)
        constr2 <- list()
        for(j in 1:num.gr){
          lg <- (group == j)
          constr2.g <- list(norm((Sigma.hat %*% m - b2)[lg,] , type="2") <= gamma2)
          constr2 <- append(constr2, constr2.g)
        }

        prob2 <- Problem(obj2, constr2)
        result2 <- solve(prob2, solver="SCS", verbose = FALSE)

        gamma2 <- result2$getValue(gamma2)
        b3 <-  rep(1,p)*max(1.1*gamma2, 1.2*gamma)

          m <- Variable(p)
          obj <- Minimize(quad_form(m, Sigma.hat))
          constr3 <- list()
          for(j in 1:num.gr){
            lg <- (group == j)
            constr3.g <- list(norm((Sigma.hat %*% m - b2)[lg,] , type="2") <= b3)
            constr3 <- append(constr3, constr3.g)
          }
          prob <- Problem(obj, constr3)
          result <- solve(prob, solver="SCS", verbose = FALSE)


        if(result$status =="optimal"){
          M[i,] <- result$getValue(m)
        }
      }

    }
  M.tilde <- Q %*% M %*% Q # transform M matrix
  return(M.tilde)
}











