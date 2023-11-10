#'
#' @title Computes the SummarizedExperiment object
#' @description This function is similar to the native R function from Bioconductor
#' @details The function calls the server-side function \code{summarizedExperimentDS} that computes the
#' SummarizedExperiment from a data.frame and assigns the object to the server-side.
#' The new object is named by the user using the \code{newobj} argument, otherwise it is named \code{sumExp.newobj} by default.
#' @param microbiomeData is a string character for a data.frame containing the microbiome variables of interest
#' @param covariateData is a string character for a data.frame containing the covariates variables of interest
#' @param newobj is the name of the new object which is created with this function
#' @param datasources a list of \code{\link{DSConnection-class}} objects obtained after login
#' @return the object specified by the \code{newobj} argument of \code{ds.summarizedExperiment} or default name \code{sumExp.newobj}
#' @author Florian Schwarz for the German Institute of Human Nutrition
#' @import DSI
#' @import dsBaseClient
#' @import methods
#' @export
#' @examples
#' \dontrun{
#'
#'   # connecting to the Opal servers
#'
#'   require('DSI')
#'   require('DSOpal')
#'   require('dsBaseClient')
#'   require('dsIntestinalMicrobiomicsClient')
#'
#'   builder <- DSI::newDSLoginBuilder()
#'   builder$append(server = "study1",
#'                  url = "http://192.168.56.100:8080/",
#'                  user = "administrator", password = "datashield_test&",
#'                  table = "MicrobSIM.MicrobSIM1", driver = "OpalDriver")
#'   builder$append(server = "study2",
#'                  url = "http://192.168.56.100:8080/",
#'                  user = "administrator", password = "datashield_test&",
#'                  table = "MicrobSIM.MicrobSIM2", driver = "OpalDriver")
#'   builder$append(server = "study3",
#'                  url = "http://192.168.56.100:8080/",
#'                  user = "administrator", password = "datashield_test&",
#'                  table = "MicrobSIM.MicrobSIM3", driver = "OpalDriver")
#'   logindata <- builder$build()
#'
#'   connections <- DSI::datashield.login(logins = logindata, assign = TRUE, symbol = "D")
#'
#'
#'   # Create data.frames with microbiome and covariate data of interest
#'
#'   ds.dataFrame(x = c("D$P_ACTINOBACTERIA",
#'                      "D$P_BACTEROIDETES",
#'                      "D$P_FIRMICUTES",
#'                      "D$P_VERRUCOMICROBIA"),
#'                newobj = "microbdata",
#'                stringsAsFactors = FALSE)
#'
#'   ds.dataFrame(x = c("D$AGE",
#'                      "D$SEX",
#'                      "D$WEIGHT",
#'                      "D$HEIGHT"),
#'                newobj = "covdata",
#'                stringsAsFactors = FALSE)
#'
#'
#'   # Create the summarizedExperiment object on the server-side based on the microbiome and covariate data.frames
#'
#'   ds.summarizedExperiment(microbiomeData = "microbdata",
#'                           covariateData = "covdata",
#'                           newobj = "SumExpT")
#'
#'   # clear the Datashield R sessions and logout
#'   datashield.logout(connections)
#' }
#'

ds.summarizedExperiment <- function(microbiomeData = NULL, covariateData = NULL, newobj = NULL, datasources = NULL){


  # look for DS connections
  if(is.null(datasources)){
    datasources <- datashield.connections_find()
  }


  # ensure datasources is a list of DSConnection-class
  if(!(is.list(datasources) && all(unlist(lapply(datasources, function(d) {methods::is(d,"DSConnection")}))))){
    stop("The 'datasources' were expected to be a list of DSConnection-class objects", call.=FALSE)
  }




  if(is.null(microbiomeData)){
    stop("Please provide the name(s) of the microbiome data!", call.=FALSE)
  }


  if(is.null(covariateData)){
    stop("Please provide the name(s) of the covariate data!", call.=FALSE)
  }




  if(is.null(newobj)){
    newobj <- "sumExp.newobj"
  }


  # call the server side function that does the operation
  cally <- call("summarizedExperimentDS", microbiomeData, covariateData)
  DSI::datashield.assign(datasources, newobj, cally)



}

# assign function






