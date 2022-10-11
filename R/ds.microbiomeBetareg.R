#### Documentation needs to be adjusted to Betareg!





#'
#' @title Fits beta regression models for rates and proportions via maximum likelihood
#' @description This function is similar to the native R function betareg from the betareg package
#' @details The function calls the server-side function \code{microbiomeBetaregDS} that computes the beta regression model
#' from data.frame. The function uses parametrization with mean and precision parameter.
#' @param df is a string character for the data.frame object
#' @param formula symbolic description of the model.
#' @param na.action is a string character specifying how NA data shall be treated.
#' @param weights optional numeric vector of case weights.
#' @param offset optional numeric vector with an a priori known component to be included in the linear predictor for the mean. In betareg.fit, offset may also be a list of two offsets for the mean and precision equation, respectively.
#' @param link character specification of the link function in the mean model (mu). Currently, "logit", "probit", "cloglog", "cauchit", "log", "loglog" are supported. Alternatively, an object of class "link-glm" can be supplied.
#' @param link.phi character specification of the link function in the precision model (phi). Currently, "identity", "log", "sqrt" are supported. The default is "log" unless formula is of type y ~ x where the default is "identity" (for backward compatibility). Alternatively, an object of class "link-glm" can be supplied.
#' @param type character specification of the type of estimator. Currently, maximum likelihood ("ML"), ML with bias correction ("BC"), and ML with bias reduction ("BR") are supported.
#' @param control a list of control arguments specified via betareg.control.
#' @param model is a logical. If TRUE, the model frame from the model fit will be returned.
#' @param x is a logical. If TRUE, the response from the model fit will be returned.
#' @param y is a logical. If TRUE, the model matrix from the model fit will be returned.
#' @param datasources a list of \code{\link{DSConnection-class}} objects obtained after login. If the \code{datasources} argument is not specified the default set of connections will be used: see \code{\link{datashield.connections_default}}.
#' @return \code{ds.microbiomeBetareg} returns XXXXXXXXXXXXXXXXXXXXXX
#' @author Florian Schwarz for the German Institute of Human Nutrition
#' @import DSI
#' @import dsBaseClient
#' @import methods
#' @export
#'

ds.microbiomeBetareg <- function(df = NULL, formula = NULL, na.action, weights, offset, link = c("logit", "probit", "cloglog", "cauchit", "log", "loglog"),
                                 link.phi = NULL, type = c("ML", "BC", "BR"), control = betareg.control(...), model = TRUE, y = TRUE, x = FALSE, datasources = NULL){


  # look for DS connections
  if(is.null(datasources)){
    datasources <- datashield.connections_find()
  }


  # ensure datasources is a list of DSConnection-class
  if(!(is.list(datasources) && all(unlist(lapply(datasources, function(d) {methods::is(d,"DSConnection")}))))){
    stop("The 'datasources' were expected to be a list of DSConnection-class objects", call.=FALSE)
  }


  if(is.null(df)){
    stop("Please provide the name of the data.frame!", call.=FALSE)
  }


  if(is.null(formula)){
    stop("Please provide a correct formula !", call.=FALSE)
  }



  # call the internal function that checks the input object is of the same class in all studies.
  typ <- dsBaseClient::ds.class(df, datasources)



  # Check whether the input is either of type data frame or matrix
  if(!('data.frame' %in% typ)){
    stop("Only objects of type 'data.frame' are allowed.", call.=FALSE)
  }



  # call the server side function that does the operation
  cally <- call("microbiomeMZILNDS", SumExp, taxa, covariates, sampleIDname, adjust_method, fdrRate, paraJobs, bootB, taxDropThresh, standardize, sequentialRun, verbose, seed)
  DSI::datashield.aggregate(datasources, cally)



}

# aggregate function





















