#'
#' @title Computes the regression for microbiome analysis based on multivariate zero-inflated logistic normal model
#' @description This function is similar to the native R function from the IFAA package
#' @details The function calls the server-side function \code{microbiomeMZILNDS} that computes the
#' regression from a SummarizedExperiment object. SummarizedExperiment objects can be computed using the \code{ds.summarizedExperiment} function.
#' @param SumExp is a string character describing the SummarizedExperiment object
#' @param taxa is a string character for the microbiome variable denominator (can also be a vector of microbiome variables)
#' @param covariates is a string character of covariates to be examined along the microbiome variables (can also be a vector of covariates).
#' @param datasources a list of \code{\link{DSConnection-class}} objects obtained after login
#' @return XXXXX
#' @author Florian Schwarz for the German Institute of Human Nutrition
#' @import DSI
#' @import dsBaseClient
#' @import methods
#' @export
#'

ds.microbiomeMZILN <- function(SumExp = NULL, taxa = NULL, covariates = NULL, datasources = NULL){


  # look for DS connections
  if(is.null(datasources)){
    datasources <- datashield.connections_find()
  }


  # ensure datasources is a list of DSConnection-class
  if(!(is.list(datasources) && all(unlist(lapply(datasources, function(d) {methods::is(d,"DSConnection")}))))){
    stop("The 'datasources' were expected to be a list of DSConnection-class objects", call.=FALSE)
  }


  if(is.null(SumExp)){
    stop("Please provide the name of the SummarizedExperiment object!", call.=FALSE)
  }


  if(is.null(taxa)){
    stop("Please provide the name(s) of the microbiome denominator data!", call.=FALSE)
  }



  # call the internal function that checks the input object is of the same class in all studies.
  typ <- dsBaseClient::ds.class(SumExp, datasources)



  # Check whether the input is either of type data frame or matrix
  if(!('SummarizedExperiment' %in% typ)){
    stop("Only objects of type 'SummarizedExperiment' are allowed.", call.=FALSE)
  }



  # call the server side function that does the operation
  cally <- call("microbiomeMZILNDS", SumExp, taxa, covariates)
  DSI::datashield.aggregate(datasources, cally)



}

# aggregate function
















