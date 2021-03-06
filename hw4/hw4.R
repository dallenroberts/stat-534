################################################################################
## Allen Roberts
## April 23, 2020
## Stat 534
## Homework 4: Main file
################################################################################

bayesLogistic <- function(apredictor, response, data, NumberOfIterations) {
  
  ## Load helper functions
  source("helperFunc.R")
  
  ## Calculate MLEs
  betaMLE <- getcoefglm(response, apredictor, data)
  beta0mle <- betaMLE[1]
  beta1mle <- betaMLE[2]
  
  ## Calculate beta modes
  betaMode <- getcoefNR(response, apredictor, data)
  
  ## Calculate posterior means for betas
  posteriorMeans <- getPosteriorMeans(response, apredictor, data, betaMode, niter = NumberOfIterations)
  beta0bayes <- posteriorMeans[1]
  beta1bayes <- posteriorMeans[2]
  
  ## Calculate log marginal likelihood
  logmarglik <- getLaplaceApprox(response, apredictor, data, betaMode)
  
  return(list("apredictor" = apredictor, 
              "logmarglik" = logmarglik,
              "beta0bayes" = beta0bayes,
              "beta1bayes" = beta1bayes,
              "beta0mle" = beta0mle,
              "beta1mle" = beta1mle))
  
}

main <- function(datafile, NumberOfIterations, clusterSize) {
  
  #read the data
  data <- read.table(datafile,header=FALSE);
  response <- ncol(data)
  lastPredictor <- ncol(data)-1
  
  #initialize a cluster for parallel computing
  cluster <- makeCluster(clusterSize, type = "SOCK")
  
  #run the MC3 algorithm from several times
  results <- clusterApply(cluster, 1:lastPredictor, bayesLogistic,
                         response,data,NumberOfIterations);
  
  ## Print results
  for(i in 1:lastPredictor) {
    cat('Regression of Y on explanatory variable ',results[[i]]$apredictor,
        ' has log marginal likelihood ',results[[i]]$logmarglik,
        ' with beta0 = ',results[[i]]$beta0bayes,' (',results[[i]]$beta0mle,')',
        ' and beta1 = ',results[[i]]$beta1bayes,' (',results[[i]]$beta1mle,')',
        '\n');    
  }
  
  #destroy the cluster
  stopCluster(cluster); 
}

require(snow)

main('534binarydata.txt', 10000, 10)

## Comments: In general, the posterior means of the regression coefficients are 
## similar but generally shrunk toward zero with respect to the MLEs, though 
## this is not always the case. I would expect the posterior means to be closer 
## to zero because the priors on each of the coefficients were centered at zero.
## It is also possible that some of the variability is due to sampling error for
## the calculation of the posterior means.
