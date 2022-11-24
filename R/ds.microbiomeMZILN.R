#'
#' @title Computes the regression for microbiome analysis based on multivariate zero-inflated logistic normal model
#' @description This function is similar to the native R function from the IFAA package
#' @details The function calls the server-side function \code{microbiomeMZILNDS} that computes the
#' regression from a SummarizedExperiment object. SummarizedExperiment objects can be computed using the \code{ds.summarizedExperiment} function.
#' @param SumExp is a string character describing the SummarizedExperiment object
#' @param taxa is a string character for the microbiome variable denominator (can also be a vector of microbiome variables)
#' @param covariates is a string character of covariates to be examined along the microbiome variables (can also be a vector of covariates).
#' @param sampleIDname is a string character for the sample ID variable.
#' @param adjust_method The adjusting method for p value adjustment. Default is "BY" for dependent FDR adjustment. It can take any adjustment method for the p.adjust function in R.
#' @param fdrRate The false discovery rate for identifying taxa/OTU/ASV associated with 'covariates'.
#' @param paraJobs If 'sequentialRun' is FALSE, this specifies the number of parallel jobs that will be registered to run the algoithm. If specified as NULL, it will automatically detect the cores to decide the number of parallel jobs.
#' @param bootB Number of bootstrap samples for obtaining confidence interval of estimates for the high dimensional regression.Default is 500.
#' @param taxDropThresh The threshold of number of non-zero sequencing reads for each taxon to be dropped from the analysis. Default is 0 which means that taxon without any sequencing reads will be dropped from the analysis.
#' @param standardize is a logical. If TRUE, the design matrix for X will be standardized in the analyses and the results. Default is FALSE.
#' @param sequentialRun is a logical. Default is TRUE. It can be set to be FALSE to increase speed if there are multiple taxa in the argument 'taxa'.
#' @param verbose Whether the process message is printed out to the console. Default is TRUE.
#' @param seed Random seed for reproducibility. Can be set to NULL to remove seeding.
#' @param datasources a list of \code{\link{DSConnection-class}} objects obtained after login
#' @return \code{ds.microbiomeMZILN} returns the outcome of the specified multivariate zero-inflated logistic normal model
#' @author Florian Schwarz for the German Institute of Human Nutrition
#' @import DSI
#' @import dsBaseClient
#' @import methods
#' @export
#'

ds.microbiomeMZILN <- function(SumExp = NULL, taxa = NULL, covariates = NULL, sampleIDname = NULL, adjust_method = "BY", fdrRate = 0.15, paraJobs = NULL,
                               bootB = 500, taxDropThresh = 0, standardize = FALSE, sequentialRun = TRUE, verbose = TRUE, seed = 1, datasources = NULL){


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
  cally <- call("microbiomeMZILNDS", SumExp, taxa, covariates, sampleIDname, adjust_method, fdrRate, paraJobs, bootB, taxDropThresh, standardize, sequentialRun, verbose, seed)
  DSI::datashield.aggregate(datasources, cally)



}

# aggregate function


