#'
#' @title Computes the regression for microbiome analysis based on multivariate zero-inflated logistic normal model
#' @description This function is similar to the native R function from the Zoib package
#' @details The function calls the server-side function \code{microbiomeZoibDS} that computes the
#' model from a data.frame object.
#' @param df is a string character describing the data.frame object to be analyzed
#' @param model symbolic description of the model in the format of formula, such as y ~ x, y1 | y2 ~ x1+x2, or y1 ~ x | z.
#' @param n.response Number of response variables. Default is 1.
#' @param joint is a logical. Whether to jointly model multiple responses if n.response >= 2. Default is TRUE.
#' @param zero.inflation A vector that contains n.response values of TRUE or FALSE on whether each of the response variables has inflation at zero. Default is TRUE.
#' @param one.inflation A vector that contains n.response values of TRUE or FALSE on whether each of the response variables has inflation at one. Default is TRUE.
#' @param random Specifies which linear predictors or link functions contain a random component. Default is 0 (no random component).
#' @param EUID Listing of the experimental unit ID in each row of the data set.
#' @param link.mu Link function for the mean of the beta distribution piece of the zoib model. Choices are "logit", "probit" and "cloglog".
#' @param link.x0 Link function for Pr(y=0). Choices are "logit", "probit" and "cloglog".
#' @param link.x1 Link function for Pr(y=1 | y>0). Choices are "logit", "probit" and "cloglog".
#' @param prec.int Precision parameter of the prior distribution (diffuse normal) of the intercept in the linear predictor of each link function. Default is 0.0001 in all 4 link functions for all response variables.
#' @param prior.beta Prior choice for the regression coefficients other than the intercepts in each of the 4 link functions. Default is "DN" (diffuse normal) in all 4 link functions for all response variables. Other options are "L1" (L1-like shrinkage prior), "L2" (L2-like shrinkage prior), "ARD" (ARD-like shrinkage prior).
#' @param prec.DN Precision parameters of the normal distribution if the diffuse normal prior is chosen as the prior distribtions of the regression coefficients in all 4 linear predictors for all response variables. Default precision is 0.0001.
#' @param lambda.L2 Scale parameter of the prior distributions of the regression coefficients in the linear predictors for all response variables if the L2-like prior is chosen.
#' @param lambda.L1 Scale parameter of the prior distributions of the regression coefficients in the linear predictors for all response variables if the L1-like prior is chosen.
#' @param lambda.ARD Scale parameter of the prior distributions of the regression coefficients in the linear predictors for all response variables if the ARD-like prior is chosen.
#' @param prior.Sigma Prior choice for the variance or the variance-covariance in the case of a single random variable and multiple random variables, respectively. The default is "VC.unif". When there is a single random variable, choose from "VC.unif" and "VC.halfcauchy".
#' @param scale.unif Upper bound of the uniform prior for the standard deviation of each random variable if prior.Sigma="VC.unif" is specified. Default is 20.
#' @param scale.halfcauchy Scale parameter of the half-Cauchy prior for the standard deviation of each random variable if prior.Sigma = "VC.halfCauchy" is specified. Default is 20.
#' @param n.chain Number of Markov chains from which posterior samples will be drawn (>=1; default = 2).
#' @param n.iter Number of iterations per chain in the MCMC sampling (default = 5000) before burning-in and thinning.
#' @param n.burn Burning in period of the MCMC chains (default = 200).
#' @param n.thin Thinning period of the MCMC chains after the burn-in (default = 5).
#' @param inits optional specification of initial values for regression coefficients and variance/covariance parameters in the form of a list. If omitted, initial values will be generated automatically.
#' @param seeds a vector of dimension n.chain that contains seeds for the initial values and the random number generators of the MCMC chains, if users wish to make the output from the model reproducible.
#' @param datasources a list of \code{\link{DSConnection-class}} objects obtained after login
#' @return \code{ds.microbiomeZoib} returns the multivariate zero-inflated logistic normal model
#' @author Florian Schwarz for the German Institute of Human Nutrition
#' @import DSI
#' @import dsBaseClient
#' @import methods
#' @export
#'

ds.microbiomeZoib <- function(df = NULL, model = NULL, n.response = 1, joint = TRUE, zero.inflation = TRUE, one.inflation = TRUE, random = 0, EUID, link.mu = "logit",
                               link.x0 = "logit", link.x1 = "logit", prec.int = matrix(1e-3, n.response, 4), prior.beta = matrix("DN", n.response, 4), prec.DN = matrix(1e-3, n.response,4),
                               lambda.L2 = matrix(1e-3, n.response, 4), lambda.L1 = matrix(1e-3, n.response, 4), lambda.ARD = matrix(1e-3, n.response,4), prior.Sigma = "VC.unif",
                               scale.unif = 20, scale.halfcauchy = 20, n.chain = 2, n.iter = 5000, n.burn = 200, n.thin = 2, inits = NULL, seeds = NULL, datasources = NULL){


  # look for DS connections
  if(is.null(datasources)){
    datasources <- datashield.connections_find()
  }


  # ensure datasources is a list of DSConnection-class
  if(!(is.list(datasources) && all(unlist(lapply(datasources, function(d) {methods::is(d,"DSConnection")}))))){
    stop("The 'datasources' were expected to be a list of DSConnection-class objects", call.=FALSE)
  }


  if(is.null(df)){
    stop("Please provide the name of the data.frame object!", call.=FALSE)
  }


  if(is.null(model)){
    stop("Please provide the a symbolic description of the model!", call.=FALSE)
  }



  # call the internal function that checks the input object is of the same class in all studies.
  typ <- dsBaseClient::ds.class(df, datasources)



  # Check whether the input is either of type data frame or matrix
  if(!('data.frame' %in% typ)){
    stop("Only objects of type 'data.frame' are allowed.", call.=FALSE)
  }



  # call the server side function that does the operation
  cally <- call("microbiomeZoibDS", df, model, n.response, joint, zero.inflation, one.inflation, random, EUID, link.mu, link.x0, link.x1, prec.int, prior.beta, prec.DN, lambda.L1, lambda.L2, lambda.ARD, prior.Sigma, scale.unif, scale.halfcauchy, n.chain, n.iter, n.burn, n.thin, inits, seeds)
  DSI::datashield.aggregate(datasources, cally)



}

# aggregate function


