rm(list = ls())

## Problem 1
logdet <- function(mat) {
  
  return(sum(log(eigen(mat, only.values = TRUE)$values)))
  
}

mat <- matrix(c(2, -1, 0, -1, 2, -1, 0, -1, 2), nrow = 3, ncol = 3)

print(paste("My function's result:", logdet(mat)))
print(paste("Result using R's log(det(mat)):", log(det(mat))))

logdet_bad <- function(mat) {
  
  return(log(prod(eigen(mat, only.values = TRUE)$values)))
  
}
print(logdet_bad(mat))

set.seed(503)
n <- 1000
p <- qr.Q(qr(matrix(rnorm(n^2), n)))
large_matrix <- crossprod(p, p*(n:1))

print(paste("My function's result:", logdet(large_matrix)))
print(paste("My bad function's result:", logdet_bad(large_matrix)))

print(paste("Result using R's log(det(large_matrix)):",
            log(det(large_matrix))))
print(paste("Result using R's determinant(large_matrix, logarithm = TRUE):",
            determinant(large_matrix, logarithm = TRUE)$modulus[1]))

## Problem 2
logmarglik <- function(data, A) {
  
  n <- nrow(data)
  
  D_1 <- as.matrix(data[, 1], nrow = n)
  D_A <- as.matrix(data[, A], nrow = n)
  M_A <- diag(length(A)) + t(D_A)%*%(D_A)
  
  lml <- lgamma((n + length(A) + 2)/2) - lgamma((length(A) + 2)/2) - 
    (1/2)*logdet(M_A) - ((n + length(A) + 2)/2) * 
    log((1 + t(D_1)%*%D_1 - t(D_1)%*%D_A%*%solve(M_A)%*%t(D_A)%*%D_1))
  
  return(as.numeric(lml))
  
}

erdata <- read.table("erdata.txt")
print(paste("logmarglik(data, c(2,5,10)) = ", logmarglik(erdata, c(2,5,10))))
