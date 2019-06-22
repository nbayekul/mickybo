# This function calculates the Frobenius norm of a given matrix

Frobenius <- function(argMatrix){
  result <- sqrt(sum(argMatrix^2))
  return(result)
}

# We want to approximate, w.r.t. the Frobenius distance, non-negative p*q matrix V
# by a product W.H: W (p*k) and H (k*q)    (1)
# subject to all entries in W and H being non-negative

# Popular idea: Alternating least squares 
# (alternately regard W or H as a given and optimise for the other one)
# The so-called 'multiplicative updating rules' were first promoted in
# https://papers.nips.cc/paper/1861-algorithms-for-non-negative-matrix-factorization.pdf
# See https://en.wikipedia.org/wiki/Non-negative_matrix_factorization (Algorithms section)

# V: matrix to be factorised
# rank: common dimensionality of the two factors (number of latent factors) ; k from (1)
# min.iter: minimum number of iterations
# max.iter: maximum number of iterations
# tol: distance threshold for convergence

solveNMF <- function(V, rank, min.iter, max.iter, tol){
  
  # Initialising the factors W and H
  
  p <- nrow(V)
  q <- ncol(V)
  vec_W <- runif(p*rank)
  vec_H <- runif(rank*q)
  
  W <- matrix(vec_W, nrow=p)
  H <- matrix(vec_H, ncol=q)
  
  old.W <- W
  old.H <- H
  
  # Using the updating rules
  
  for (k in 1:max.iter){
    # update W
    old.W <- W
    for (i in 1:nrow(W)){
      for (j in 1:ncol(W)){
        W[i,j] <- (W[i,j]*(V%*%t(H))[i,j])/(W%*%H%*%t(H))[i,j]
      }
    }
    # update H
    old.H <- H
    for (i in 1:nrow(H)){
      for (j in 1:ncol(H)){
        H[i,j] <- (H[i,j]*(t(W)%*%V)[i,j])/(t(W)%*%W%*%H)[i,j]
      }
    }
    if( (k>= min.iter) && (Frobenius(W-old.W)+Frobenius(H-old.H)) <= tol){
      print(paste("Convergence in ", k, " iterations"))
      print("Convergence occurred before the set limit for the number of iterations.")
      break
    }
  }
  if (k==max.iter){
     print("Algorithm stopped because maximum number of iterations was reached.")
  }
  return(list(W,H))
}


# Example 
set.seed(2000)
V <- matrix(runif(14*20),ncol=14)
# Low number of iterations allowed
result1 <- solveNMF(V=V, rank=5, min.iter=5, max.iter=7, tol=0.01)
# Higher number of iterations allowed
result2 <- solveNMF(V=V, rank=5, min.iter=100, max.iter=500, tol=0.01)


