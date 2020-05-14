################################################################################
## Allen Roberts
## April 23, 2020
## Stat 534
## Homework 4: Helper functions
################################################################################

require(MASS)

## Calculate the log determinant of a matrix
logdet <- function(mat)
{
  return(sum(log(eigen(mat)$values)))	
}


#this function uses 'glm' to fit a logistic regression
#and obtain the MLEs of the two coefficients beta0 and beta1
getcoefglm <- function(response,explanatory,data)
{
  return(coef(glm(data[,response] ~ data[,explanatory],
                  family=binomial(link=logit))));
}

#the inverse of the logit function
inverseLogit <- function(x)
{
  return(exp(x)/(1+exp(x))); 
}

#function for the computation of the Hessian
inverseLogit2 <- function(x)
{
  return(exp(x)/(1+exp(x))^2); 
}


#computes pi_i = P(y_i = 1 | x_i)
getPi <- function(x,beta)
{
  x0 = cbind(rep(1,length(x)),x);
  return(inverseLogit(x0%*%beta));
}

#another function for the computation of the Hessian
getPi2 <- function(x,beta)
{
  x0 = cbind(rep(1,length(x)),x);
  return(inverseLogit2(x0%*%beta));
}

#logistic log-likelihood (formula (3) in your handout)
logisticLoglik <- function(y,x,beta)
{
  Pi = getPi(x,beta);
  return(sum(y*log(Pi))+sum((1-y)*log(1-Pi)));
}

## Logistic log-likelihood star (from Bayesian logistic regression eq 2.5)
logisticLoglikStar <- function(y, x, beta) {
  
  loglik <- logisticLoglik(y,x, beta)
  
  return(-log(2*pi) - 1/2*(sum(beta^2)) + loglik)
  
}

#obtain the gradient for Newton-Raphson
getGradient <- function(y,x,beta)
{
  gradient = matrix(0,2,1);
  Pi = getPi(x,beta);
  
  # gradient[1,1] = sum(y-Pi);
  # gradient[2,1] = sum((y-Pi)*x);
  
  ## Updated to work with Bayesian model
  gradient[1,1] = sum(y-Pi) - beta[1];
  gradient[2,1] = sum((y-Pi)*x) - beta[2];
  
  return(gradient);
}

#obtain the Hessian for Newton-Raphson
getHessian <- function(y,x,beta)
{
  hessian = matrix(0,2,2);
  Pi2 = getPi2(x,beta);
  
  # hessian[1,1] = sum(Pi2);
  # hessian[1,2] = sum(Pi2*x);
  # hessian[2,1] = hessian[1,2];
  # hessian[2,2] = sum(Pi2*x^2);
  
  ## Updated to work with Bayesian model
  hessian[1,1] = sum(Pi2) + 1;
  hessian[1,2] = sum(Pi2*x);
  hessian[2,1] = hessian[1,2];
  hessian[2,2] = sum(Pi2*x^2) + 1;
  
  return(-hessian);
}

#this function implements our own Newton-Raphson procedure
getcoefNR <- function(response,explanatory,data, maxIter = 1000)
{
  #2x1 matrix of coefficients`
  beta = matrix(0,2,1);
  y = data[,response];
  x = data[,explanatory];
  
  #current value of log-likelihood
  currentLoglik = logisticLoglikStar(y,x,beta);
  
  iter <- 0
  
  #infinite loop unless we stop it someplace inside
  while(iter < maxIter)
  {
    iter <- iter + 1
    
    newBeta = beta - solve(getHessian(y,x,beta))%*%getGradient(y,x,beta);
    newLoglik = logisticLoglikStar(y,x,newBeta);
    
    #at each iteration the log-likelihood must increase
    if(newLoglik<currentLoglik)
    {
      cat("CODING ERROR!!\n");
      break;
    }
    beta = newBeta;
    #stop if the log-likelihood does not improve by too much
    if(newLoglik-currentLoglik<1e-6)
    {
      break; 
    }
    currentLoglik = newLoglik;
  }
  
  return(beta);
}


## getLaplaceApprox uses the Laplace Approximation to calcuate an approximate 
## marginal likelihood for univariate Bayesian logistic regression. Note that 
## this function calculates and returns the log-likehood, as specified by the 
## TA on the discussion board
getLaplaceApprox <- function(response, explanatory, data, betaMode) {
  
  y <- data[,response]
  x <- data[,explanatory]
  
  ## Calculate l*(beta_0, beta_1)
  loglik_star <- logisticLoglikStar(y,x, 
                                    beta = betaMode)
  
  ## Calculate the log marginal likelihood
  log_p_d <- log(2*pi) + loglik_star - (1/2)*logdet(-getHessian(y,x,betaMode))
  
  return(log_p_d)
  
}

## sampleMH performs Metropolis-Hastings sampling from the posterior distribution
## P(beta0, beta1 | D) of Bayesian univariate logistic regression. It returns
## a matrix with two columns (beta0 and beta1) and niter rows (one for each)
## sample.
sampleMH <- function(response, explanatory, data, betaMode, niter) {
  
  y <- data[, response]
  x <- data[, explanatory]
  
  samples <- matrix(0, nrow = niter, ncol = 2)
  
  ## Proposal distribution covariance matrix
  covMat <- -solve(getHessian(y, x, beta = betaMode))
  
  ## Initial state
  currentBeta <- betaMode
  
  ## Start Markov chain
  for(k in 1:niter) {
    
    ## Sample candidate beta from proposal distribution
    candidateBeta <- mvrnorm(n = 1, mu = currentBeta, Sigma = covMat)
    
    ## Accept or reject candidate beta
    score <- logisticLoglikStar(y, x, candidateBeta) - logisticLoglikStar(y, x, currentBeta)
    
    if(score >= 1) {
      
      currentBeta <- candidateBeta
      
    } else {
      
      u <- runif(n = 1, min = 0, max = 1)
      
      if(log(u) <= score) {
        
        currentBeta <- candidateBeta
      }
    }
    
    ## Update chain
    samples[k, ] <- currentBeta
  }
  
  # hist(samples[, 1])
  # plot(density(samples[,1]))
  # abline(v = mean(samples[,1]))
  # hist(samples[, 2])
  # plot(density(samples[,2]))
  # abline(v = mean(samples[,2]))
  
  return(samples)
}

## Calculates the posterior means of niter samples from the joint distribution
## of the betas given the observed data. Sampling is done via Metropolis-
## Hastings. Returns a matrix of two beta values.
getPosteriorMeans <- function(response, explanatory, data, betaMode, niter) {
  
  ## Simulate 1000 samples using Metropolis Hastings
  samples <- sampleMH(response, explanatory, data, betaMode, niter)
  
  ## Calculate sample means
  sampleMeans <- apply(samples, 2, mean)
  
  ## Return vector of beta values
  return(sampleMeans)
}
