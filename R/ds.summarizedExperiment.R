#'
#' @title Computes the SummarizedExperiment object
#' @description This function is similar to the native R function from Bioconductor
#' @details The function calls the server-side function \code{summarizedExperimentDS} that computes the
#' SummarizedExperiment from a data.frame and assigns the object to the server-side.
#' The new object is named by the user using the \code{newobj} argument, otherwise it is named \code{sumExp.newobj} by default.
#' @param df is a string character of the data.frame
#' @param microbiomeData is a string character of a microbiome variable of interest (can also be a vector of microbiome variables)
#' @param covariateData is a string character of a covariate variable of interest (can also be a vector of covariates)
#' @param newobj is the name of the new object which is created with this function
#' @param datasources a list of \code{\link{DSConnection-class}} objects obtained after login
#' @return the object specified by the \code{newobj} argument of \code{ds.summarizedExperiment} or default name \code{sumExp.newobj}
#' @author Florian Schwarz for the German Institute of Human Nutrition
#' @import DSI
#' @import dsBaseClient
#' @import methods
#' @export
#'

ds.summarizedExperiment <- function(df = NULL, microbiomeData = NULL, covariateData = NULL, newobj = NULL, datasources = NULL){


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


  if(is.null(microbiomeData)){
    stop("Please provide the name(s) of the microbiome data!", call.=FALSE)
  }


  if(is.null(covariateData)){
    stop("Please provide the name(s) of the covariate data!", call.=FALSE)
  }


  # call the internal function that checks the input object is of the same class in all studies.
  typ <- dsBaseClient::ds.class(df, datasources)


  # Check whether the input is either of type data frame or matrix
  if(!('data.frame' %in% typ)){
    stop("Only objects of type 'data frame' are allowed.", call.=FALSE)
  }


  if(is.null(newobj)){
    newobj <- "sumExp.newobj"
  }


  # call the server side function that does the operation
  cally <- call("summarizedExperimentDS", df, microbiomeData, covariateData)
  DSI::datashield.assign(datasources, newobj, cally)



}

# assign function






