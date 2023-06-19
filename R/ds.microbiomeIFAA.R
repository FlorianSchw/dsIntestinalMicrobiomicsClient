#'
#' @title Computes the association of microbiome data with covariates
#' @description This function is similar to the native R function from the IFAA package
#' @details The function calls the server-side function \code{microbiomeIFAADS} that computes the
#' association analysis from a SummarizedExperiment object. SummarizedExperiment objects can be computed using the \code{ds.summarizedExperiment} function.
#' @param SumExp is a string character describing the SummarizedExperiment object
#' @param microbVar This takes a single or vector of microbiome variable names (e.g., taxa, OTU and ASV names) of interest. Default is "all" meaning all microbiome variables will be analyzed. If a subset of microbiome variables is specified, the output will only contain the specified variables, and p-value adjustment for multiple testing will only be applied to the subset.
#' @param testCov is a string character of covariates to be examined along the microbiome variables (can also be a vector of covariates).
#' @param ctrlCov is a string character for the covariates that will be adjusted in the model (can also be a vector of microbiome variables)
#' @param sampleIDname is a string character for the sample ID variable.
#' @param testMany is a logical. If 'TRUE' and 'covariates' are set to NULL, then all variables in the 'covariates' will be used.
#' @param ctrlMany is a logical. If 'TRUE' and 'confounders' are set to NULL, then all variables except the 'coviariates' will be used as confounders.
#' @param nRef number of randomly picked reference taxa used in phase 1.
#' @param nRefMaxForEsti maximum number of final reference taxa used in phase 2.
#' @param refTaxa vector of taxa or OTU or ASV names. Theses are reference taxa specified by the user to be used in phase 1.
#' @param adjust_method The adjusting method for p value adjustment. Default is "BY" for dependent FDR adjustment. It can take any adjustment method for the p.adjust function in R.
#' @param fdrRate The false discovery rate for identifying taxa/OTU/ASV associated with 'covariates'.
#' @param paraJobs If 'sequentialRun' is FALSE, this specifies the number of parallel jobs that will be registered to run the algoithm. If specified as NULL, it will automatically detect the cores to decide the number of parallel jobs.
#' @param bootB number of bootstrap samples for obtaining confidence interval of estimates in phase 2 for the high dimensional regression. Default is 500.
#' @param standardize is a logical. If 'TRUE', the design matrix for X will be standardized in the analyses and the results. Default is FALSE.
#' @param sequentialRun is a logical. Defines whether there are parallel jobs or not.
#' @param refReadsThresh The threshold of proportion of non-zero sequencing reads for choosing the reference taxon in phase 2. Default is 0.2 meaning that at least 20\% non-zero sequencing reads are necessary.
#' @param taxDropThresh The threshold of number of non-zero sequencing reads for each taxon to be dropped from the analysis. Default is 0 which means that taxon without any sequencing reads will be dropped from the analysis.
#' @param SDTresh The threshold of standard deviations of sequencing reads for being chosen as the reference taxon in phase 2. The default is 0.05 which means the standard deviation of sequencing reads should be at least 0.05 in order to be chosen as a reference taxon.
#' @param SDquantilThresh The threshold of the quantile of standard deviation of sequencing reads, above which could be selected as a reference taxon. Default is 0.
#' @param balanceCut The threshold of the proportion of non-zero sequencing reads in each group of a binary variable for choosing the final reference taxa in phase 2. The default is 0.2 meaning at least 20\% non-zero sequencing reads in each group are needed to be eligible for being chosen as a final reference taxon.
#' @param verbose Whether the process message is printed out to the console. Default is TRUE.
#' @param datasources a list of \code{\link{DSConnection-class}} objects obtained after login
#' @return \code{ds.microbiomeIFAA} returns the association of the microbiome data with the covariates
#' @author Florian Schwarz for the German Institute of Human Nutrition
#' @import DSI
#' @import dsBaseClient
#' @import methods
#' @export
#'

ds.microbiomeIFAA <- function(SumExp = NULL, microbVar = "all", testCov = NULL, ctrlCov = NULL, sampleIDname = NULL, testMany = TRUE, ctrlMany = FALSE,
                              nRef = 40, nRefMaxForEsti = 2, refTaxa = NULL, adjust_method = "BY", fdrRate = 0.15, paraJobs = NULL, bootB = 500,
                              standardize = FALSE, sequentialRun = FALSE, refReadsThresh = 0.2, taxDropThresh = 0, SDThresh = 0.05, SDquantilThresh = 0,
                              balanceCut = 0.2, verbose = TRUE, type = c("split", "pooled"), datasources = NULL){


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

  if(is.null(testCov)){
    stop("Please provide name(s) of the covariate(s)!", call.=FALSE)
  }

  if(is.null(ctrlCov)){
    stop("Please provide name(s) of the confounders!", call.=FALSE)
  }

  # call the internal function that checks the input object is of the same class in all studies.
  typ <- dsBaseClient::ds.class(SumExp, datasources)



  # Check whether the input is either of type data frame or matrix
  if(!('SummarizedExperiment' %in% typ)){
    stop("Only objects of type 'SummarizedExperiment' are allowed.", call.=FALSE)
  }

  if(type == "split"){

    # call the server side function that does the operation
    cally <- call("microbiomeIFAADS", SumExp, microbVar, testCov, ctrlCov, sampleIDname, testMany, ctrlMany, nRef, nRefMaxForEsti, refTaxa, adjust_method,
                  fdrRate, paraJobs, bootB, standardize, sequentialRun, refReadsThresh, taxDropThresh, SDThresh, SDquantilThresh, balanceCut, verbose)
    outcome <- DSI::datashield.aggregate(datasources, cally)

  }


  if(type == "pooled"){

    # call the server side function that does the operation
    cally <- call("microbiomeIFAAPooledDS", SumExp, microbVar, testCov, ctrlCov, sampleIDname, testMany, ctrlMany, nRef, nRefMaxForEsti, refTaxa, adjust_method,
                  fdrRate, paraJobs, bootB, standardize, sequentialRun, refReadsThresh, taxDropThresh, SDThresh, SDquantilThresh, balanceCut, verbose)
    outcome <- DSI::datashield.aggregate(datasources, cally)




    XTX_comb_list <- list()
    Xy_comb_list <- list()
    yy_comb_list <- list()
    dfr_comb_list <- list()

    results <- list()
    results_origD <- list()
    selectRegroup <- list()
    selecList <- list()
    estList <- list()

    #### hard coded nRuns for now

    for (k in 1:length(outcome[[1]]$analysisResults)){
      for (q in 1:length(datasources)){

        XTX_comb_list[[q]] <- outcome[[q]][["analysisResults"]][[k]][[1]][[1]]
        Xy_comb_list[[q]] <- outcome[[q]][["analysisResults"]][[k]][[1]][[2]]
        yy_comb_list[[q]] <- outcome[[q]][["analysisResults"]][[k]][[1]][[3]]
        dfr_comb_list[[q]] <- outcome[[q]][["analysisResults"]][[k]][[1]][[5]]

      }

      XTX_comb <- Reduce('+', XTX_comb_list)
      Xy_comb <- Reduce('+', Xy_comb_list)
      yy_comb <- Reduce('+', yy_comb_list)
      dfr_comb <- Reduce('+', dfr_comb_list)

      nvar <- outcome[[1]][["analysisResults"]][[k]][[1]][[4]]
      keep <- outcome[[1]][["analysisResults"]][[k]][[1]][[6]]
      xnames <- outcome[[1]][["analysisResults"]][[k]][[1]][[7]]
      intercept <- outcome[[1]][["analysisResults"]][[k]][[1]][[8]]
      nPredics <- outcome[[1]][["analysisResults"]][[k]][[1]][[9]]
      fwerRate <- outcome[[1]][["analysisResults"]][[k]][[1]][[10]]
      nRuns <- outcome[[1]][["analysisResults"]][[k]][[1]][[11]]
      nAlphaSelec <- outcome[[1]][["analysisResults"]][[k]][[1]][[12]]
      nAlphaNoInt <- outcome[[1]][["analysisResults"]][[k]][[1]][[13]]
      nTaxa <- outcome[[1]][["analysisResults"]][[k]][[1]][[14]]
      ii <- outcome[[1]][["analysisResults"]][[k]][[1]][[15]]
      testCovInd <- outcome[[1]][["analysisResults"]][[k]][[1]][[16]]
      taxaNames <- outcome[[1]][["analysisResults"]][[k]][[1]][[17]]
      goodRefTaxaCandi <- outcome[[1]][["analysisResults"]][[k]][[1]][[18]]
      refTaxa_DS <- outcome[[1]][["analysisResults"]][[k]][[1]][[19]]



      coef_ds   <- Matrix::solve(XTX_comb, Xy_comb, tol = 1e-7)

      coefficients_ds <- rep(NA, nvar)
      coefficients_ds <- round(coef_ds@x,6)

      RSS_ds <- yy_comb - 2 * MatrixExtra::crossprod(coef_ds, Xy_comb) + MatrixExtra::crossprod(coef_ds, MatrixExtra::crossprod(XTX_comb, coef_ds))

      var_res_ds <- as.numeric(RSS_ds)/(dfr_comb + 2*nvar)

      se_coef_ds <- rep(NA, nvar)
      inv_ds     <- Matrix::solve(XTX_comb, diag(nrow(XTX_comb)), tol = 1e-7)

      se_coef_ds[keep] <- sqrt(var_res_ds * Matrix::diag(inv_ds))
      t1            <- coefficients_ds/se_coef_ds
      p             <- 2 * pt(abs(t1), df = dfr_comb, lower.tail = FALSE)



      coefMat<-data.frame(estimate  = coefficients_ds,
                          std.error = se_coef_ds,
                          t = t1,
                          p.value   = p)
      if(length(xnames)>0){
        if(intercept) xnames[1] <- "intercept"
        row.names(coefMat) <- xnames
      }else
      {if(intercept) row.names(coefMat) <- c("intercept",seq(nvar-1))
      }



      bootResu <- coefMat


      p_value_est<-bootResu[,4]
      p_value_est_noint<-p_value_est[-seq(1,length(p_value_est),by=(nPredics+1))]
      p_value_est_noint_adj<-p.adjust(p_value_est_noint,adjust_method)
      p_value_est_noint_adj[is.na(p_value_est_noint_adj)]<-1

      coef_est<-abs(bootResu[,1])
      coef_est_noint<-coef_est[-seq(1,length(coef_est),by=(nPredics+1))]
      coef_est_noint[is.na(coef_est_noint)]<-max(coef_est_noint,na.rm = TRUE)


      # return
      results$betaNoInt=p_value_est_noint_adj<fwerRate
      results$betaInt=p_value_est
      results$coef_est_noint=coef_est_noint


      #### originDataScreen part

      BetaNoInt.k <- 0 + (results$betaNoInt != 0)
      EstNoInt.k <- abs(results$coef_est_noint)


      #### the nRuns check below is the last bit of the nRuns loop in originDataScreen
      #### nRuns is not properly coded yet to work with nRuns > 1
      if (nRuns == 1) {
        BetaNoInt.i <- BetaNoInt.k
        EstNoInt.i <- EstNoInt.k
      }
      if (nRuns > 1) {
        BetaNoInt.i <- BetaNoInt.i + BetaNoInt.k
        EstNoInt.i <- EstNoInt.i + EstNoInt.k
      }


      BetaNoInt.i <- BetaNoInt.i / nRuns
      EstNoInt.i <- EstNoInt.i / nRuns


      selection.i <- rep(0, nAlphaSelec)
      coef.i <- rep(0, nAlphaSelec)

      if (ii == 1) {
        selection.i[-seq(1, nPredics)] <- BetaNoInt.i
        coef.i[-seq(1, nPredics)] <- EstNoInt.i
      }

      if (ii == nTaxa) {
        selection.i[-seq((nAlphaSelec - nPredics + 1), nAlphaSelec)] <- BetaNoInt.i
        coef.i[-seq((nAlphaSelec - nPredics + 1), nAlphaSelec)] <- EstNoInt.i
      }

      if ((ii > 1) & (ii < nTaxa)) {
        selection.i[seq_len((nPredics * (ii - 1)))] <- BetaNoInt.i[seq_len((nPredics * (ii - 1)))]
        selection.i[(nPredics * ii + 1):nAlphaSelec] <- BetaNoInt.i[(nPredics * (ii - 1) + 1):nAlphaNoInt]
        coef.i[seq_len((nPredics * (ii - 1)))] <- EstNoInt.i[seq_len((nPredics * (ii - 1)))]
        coef.i[(nPredics * ii + 1):nAlphaSelec] <- EstNoInt.i[(nPredics * (ii - 1) + 1):nAlphaNoInt]
      }


      selecList[[k]] <- selection.i
      estList[[k]] <- coef.i

      #### loop for each taxa ends below

    }

    scr1ResuSelec <- DescTools::DoCall(cbind, selecList)
    scr1ResuEst <- DescTools::DoCall(cbind, estList)



    # create count of selection for individual testCov
    countOfSelecForAllPred <- matrix(rowSums(scr1ResuSelec), nrow = nPredics)
    EstOfAllPred <- matrix(rowMeans(scr1ResuEst), nrow = nPredics)

    testCovCountMat <- countOfSelecForAllPred[testCovInd, , drop = FALSE]
    testEstMat <- EstOfAllPred[testCovInd, , drop = FALSE]


    # create overall count of selection for all testCov as a whole
    countOfSelecForAPred <- matrix(colSums(testCovCountMat), nrow = 1)
    estOfSelectForAPred <- matrix(colSums(testEstMat), nrow = 1)

    colnames(countOfSelecForAPred) <- taxaNames
    colnames(estOfSelectForAPred) <- taxaNames


    #### runScrParal part

    results_origD$testCovCountMat <- testCovCountMat
    results_origD$testEstMat <- testEstMat
    results_origD$countOfSelecForAPred <- countOfSelecForAPred
    results_origD$estOfSelectForAPred <- estOfSelectForAPred

    nTestCov <- length(testCovInd)
    results_origD$nTestCov <- nTestCov
    results_origD$nTaxa <- nTaxa
    results_origD$nPredics <- nPredics

    results_origD$taxaNames <- taxaNames
    results_origD$refTaxa <- refTaxa_DS
    results_origD$goodRefTaxaCandi <- goodRefTaxaCandi




    #### getScrResu part
    #### can probably be shortened a lot



    selectRegroup$refTaxa <- results_origD$refTaxa


    if (nTestCov == 1) {
      selectRegroup$selecCountMatIndv <- results_origD$countOfSelecForAPred
      selectRegroup$selecEstMatIndv <- results_origD$estOfSelectForAPred
    }

    if (nTestCov > 1) {
      selectRegroup$selecCountMatIndv <- results_origD$testCovCountMat
      selectRegroup$selecEstMatIndv <- results_origD$testEstMat

    }

    goodIndpRefTaxWithCount <- results_origD$countOfSelecForAPred[1, (colnames(results_origD$countOfSelecForAPred) %in% goodRefTaxaCandi)]
    goodIndpRefTaxWithEst <- results_origD$estOfSelectForAPred[1, (colnames(results_origD$estOfSelectForAPred) %in% goodRefTaxaCandi)]

    restRefTaxWithCount <- results_origD$countOfSelecForAPred[1, !(colnames(results_origD$countOfSelecForAPred) %in% goodRefTaxaCandi)]
    restRefTaxWithEst <- results_origD$estOfSelectForAPred[1, !(colnames(results_origD$estOfSelectForAPred) %in% goodRefTaxaCandi)]

    sort_goodIndpRefTaxWithCount <- goodIndpRefTaxWithCount[order(goodIndpRefTaxWithCount, abs(goodIndpRefTaxWithEst))]
    sort_goodIndpRefTaxWithEst <- abs(goodIndpRefTaxWithEst[order(goodIndpRefTaxWithCount, abs(goodIndpRefTaxWithEst))])

    sort_restRefTaxWithCount <- restRefTaxWithCount[order(restRefTaxWithCount, abs(restRefTaxWithEst))]
    sort_restRefTaxWithEst <- abs(restRefTaxWithEst[order(restRefTaxWithCount, abs(restRefTaxWithEst))])


    selectRegroup$finalIndpRefTax <- names(c(sort_goodIndpRefTaxWithCount, sort_restRefTaxWithCount))[seq_len(2)]

    selectRegroup$selecCountOverall <- results_origD$countOfSelecForAPred

    selectRegroup$goodIndpRefTaxWithCount <- c(sort_goodIndpRefTaxWithCount, sort_restRefTaxWithCount)
    selectRegroup$goodIndpRefTaxWithEst <- c(sort_goodIndpRefTaxWithEst, sort_restRefTaxWithEst)

    selectRegroup$goodRefTaxaCandi <- goodRefTaxaCandi


    #### regulariz part before the while loop
    nRef_smaller <- max(2, ceiling(nRef / 2))
    while_loop_ind <- FALSE
    loop_num <- 0

    nRef_smaller_num <- nRef_smaller
    nRef_smaller <- paste0(as.character(nRef_smaller))

    #### potential for an internal function? it is roughly the same procedure as above

    while (while_loop_ind == FALSE) {
      if (loop_num >= 2) {
        break
      }
      loop_num <- loop_num + 1

      refTaxa_smaller <- head(names(selectRegroup$goodIndpRefTaxWithCount), n = nRef_smaller_num)
      refTaxa_smaller <- paste0(as.character(refTaxa_smaller), collapse = ",")

      fin_ref_1 <- selectRegroup$finalIndpRefTax
      ref_taxa_1 <- selectRegroup$refTaxa

      #### here starts the long section

      # call the server side function that does the operation
      callz <- call("microbiomeIFAAPooledDS2", SumExp, microbVar, testCov, ctrlCov, sampleIDname, testMany, ctrlMany, nRef, nRefMaxForEsti, refTaxa, adjust_method,
                    fdrRate, paraJobs, bootB, standardize, sequentialRun, refReadsThresh, taxDropThresh, SDThresh, SDquantilThresh, balanceCut, verbose, nRef_smaller, refTaxa_smaller)
      outcome_while <- DSI::datashield.aggregate(datasources, callz)


      XTX_comb_list <- list()
      Xy_comb_list <- list()
      yy_comb_list <- list()
      dfr_comb_list <- list()

      results <- list()
      results_origD <- list()
      selectRegroup <- list()
      selecList <- list()
      estList <- list()

      #### hard coded nRuns for now

      for (k in 1:length(outcome_while[[1]]$analysisResults)){
        for (q in 1:length(datasources)){

          XTX_comb_list[[q]] <- outcome_while[[q]][["analysisResults"]][[k]][[1]][[1]]
          Xy_comb_list[[q]] <- outcome_while[[q]][["analysisResults"]][[k]][[1]][[2]]
          yy_comb_list[[q]] <- outcome_while[[q]][["analysisResults"]][[k]][[1]][[3]]
          dfr_comb_list[[q]] <- outcome_while[[q]][["analysisResults"]][[k]][[1]][[5]]

        }

        XTX_comb <- Reduce('+', XTX_comb_list)
        Xy_comb <- Reduce('+', Xy_comb_list)
        yy_comb <- Reduce('+', yy_comb_list)
        dfr_comb <- Reduce('+', dfr_comb_list)

        nvar <- outcome_while[[1]][["analysisResults"]][[k]][[1]][[4]]
        keep <- outcome_while[[1]][["analysisResults"]][[k]][[1]][[6]]
        xnames <- outcome_while[[1]][["analysisResults"]][[k]][[1]][[7]]
        intercept <- outcome_while[[1]][["analysisResults"]][[k]][[1]][[8]]
        nPredics <- outcome_while[[1]][["analysisResults"]][[k]][[1]][[9]]
        fwerRate <- outcome_while[[1]][["analysisResults"]][[k]][[1]][[10]]
        nRuns <- outcome_while[[1]][["analysisResults"]][[k]][[1]][[11]]
        nAlphaSelec <- outcome_while[[1]][["analysisResults"]][[k]][[1]][[12]]
        nAlphaNoInt <- outcome_while[[1]][["analysisResults"]][[k]][[1]][[13]]
        nTaxa <- outcome_while[[1]][["analysisResults"]][[k]][[1]][[14]]
        ii <- outcome_while[[1]][["analysisResults"]][[k]][[1]][[15]]
        testCovInd <- outcome_while[[1]][["analysisResults"]][[k]][[1]][[16]]
        taxaNames <- outcome_while[[1]][["analysisResults"]][[k]][[1]][[17]]
        goodRefTaxaCandi <- outcome_while[[1]][["analysisResults"]][[k]][[1]][[18]]
        refTaxa_DS <- outcome_while[[1]][["analysisResults"]][[k]][[1]][[19]]



        coef_ds   <- Matrix::solve(XTX_comb, Xy_comb, tol = 1e-7)

        coefficients_ds <- rep(NA, nvar)
        coefficients_ds <- round(coef_ds@x,6)

        RSS_ds <- yy_comb - 2 * MatrixExtra::crossprod(coef_ds, Xy_comb) + MatrixExtra::crossprod(coef_ds, MatrixExtra::crossprod(XTX_comb, coef_ds))

        var_res_ds <- as.numeric(RSS_ds)/(dfr_comb + 2*nvar)

        se_coef_ds <- rep(NA, nvar)
        inv_ds     <- Matrix::solve(XTX_comb, diag(nrow(XTX_comb)), tol = 1e-7)

        se_coef_ds[keep] <- sqrt(var_res_ds * Matrix::diag(inv_ds))
        t1            <- coefficients_ds/se_coef_ds
        p             <- 2 * pt(abs(t1), df = dfr_comb, lower.tail = FALSE)



        coefMat<-data.frame(estimate  = coefficients_ds,
                            std.error = se_coef_ds,
                            t = t1,
                            p.value   = p)
        if(length(xnames)>0){
          if(intercept) xnames[1] <- "intercept"
          row.names(coefMat) <- xnames
        }else
        {if(intercept) row.names(coefMat) <- c("intercept",seq(nvar-1))
        }



        bootResu <- coefMat


        p_value_est<-bootResu[,4]
        p_value_est_noint<-p_value_est[-seq(1,length(p_value_est),by=(nPredics+1))]
        p_value_est_noint_adj<-p.adjust(p_value_est_noint,adjust_method)
        p_value_est_noint_adj[is.na(p_value_est_noint_adj)]<-1

        coef_est<-abs(bootResu[,1])
        coef_est_noint<-coef_est[-seq(1,length(coef_est),by=(nPredics+1))]
        coef_est_noint[is.na(coef_est_noint)]<-max(coef_est_noint,na.rm = TRUE)


        # return
        results$betaNoInt=p_value_est_noint_adj<fwerRate
        results$betaInt=p_value_est
        results$coef_est_noint=coef_est_noint


        #### originDataScreen part

        BetaNoInt.k <- 0 + (results$betaNoInt != 0)
        EstNoInt.k <- abs(results$coef_est_noint)


        #### the nRuns check below is the last bit of the nRuns loop in originDataScreen
        #### nRuns is not properly coded yet to work with nRuns > 1
        if (nRuns == 1) {
          BetaNoInt.i <- BetaNoInt.k
          EstNoInt.i <- EstNoInt.k
        }
        if (nRuns > 1) {
          BetaNoInt.i <- BetaNoInt.i + BetaNoInt.k
          EstNoInt.i <- EstNoInt.i + EstNoInt.k
        }


        BetaNoInt.i <- BetaNoInt.i / nRuns
        EstNoInt.i <- EstNoInt.i / nRuns


        selection.i <- rep(0, nAlphaSelec)
        coef.i <- rep(0, nAlphaSelec)

        if (ii == 1) {
          selection.i[-seq(1, nPredics)] <- BetaNoInt.i
          coef.i[-seq(1, nPredics)] <- EstNoInt.i
        }

        if (ii == nTaxa) {
          selection.i[-seq((nAlphaSelec - nPredics + 1), nAlphaSelec)] <- BetaNoInt.i
          coef.i[-seq((nAlphaSelec - nPredics + 1), nAlphaSelec)] <- EstNoInt.i
        }

        if ((ii > 1) & (ii < nTaxa)) {
          selection.i[seq_len((nPredics * (ii - 1)))] <- BetaNoInt.i[seq_len((nPredics * (ii - 1)))]
          selection.i[(nPredics * ii + 1):nAlphaSelec] <- BetaNoInt.i[(nPredics * (ii - 1) + 1):nAlphaNoInt]
          coef.i[seq_len((nPredics * (ii - 1)))] <- EstNoInt.i[seq_len((nPredics * (ii - 1)))]
          coef.i[(nPredics * ii + 1):nAlphaSelec] <- EstNoInt.i[(nPredics * (ii - 1) + 1):nAlphaNoInt]
        }


        selecList[[k]] <- selection.i
        estList[[k]] <- coef.i

        #### loop for each taxa ends below

      }

      scr1ResuSelec <- DescTools::DoCall(cbind, selecList)
      scr1ResuEst <- DescTools::DoCall(cbind, estList)



      # create count of selection for individual testCov
      countOfSelecForAllPred <- matrix(rowSums(scr1ResuSelec), nrow = nPredics)
      EstOfAllPred <- matrix(rowMeans(scr1ResuEst), nrow = nPredics)

      testCovCountMat <- countOfSelecForAllPred[testCovInd, , drop = FALSE]
      testEstMat <- EstOfAllPred[testCovInd, , drop = FALSE]


      # create overall count of selection for all testCov as a whole
      countOfSelecForAPred <- matrix(colSums(testCovCountMat), nrow = 1)
      estOfSelectForAPred <- matrix(colSums(testEstMat), nrow = 1)

      colnames(countOfSelecForAPred) <- taxaNames
      colnames(estOfSelectForAPred) <- taxaNames


      #### runScrParal part

      results_origD$testCovCountMat <- testCovCountMat
      results_origD$testEstMat <- testEstMat
      results_origD$countOfSelecForAPred <- countOfSelecForAPred
      results_origD$estOfSelectForAPred <- estOfSelectForAPred

      nTestCov <- length(testCovInd)
      results_origD$nTestCov <- nTestCov
      results_origD$nTaxa <- nTaxa
      results_origD$nPredics <- nPredics

      results_origD$taxaNames <- taxaNames
      results_origD$refTaxa <- refTaxa_DS
      results_origD$goodRefTaxaCandi <- goodRefTaxaCandi




      #### getScrResu part
      #### can probably be shortened a lot



      selectRegroup$refTaxa <- results_origD$refTaxa


      if (nTestCov == 1) {
        selectRegroup$selecCountMatIndv <- results_origD$countOfSelecForAPred
        selectRegroup$selecEstMatIndv <- results_origD$estOfSelectForAPred
      }

      if (nTestCov > 1) {
        selectRegroup$selecCountMatIndv <- results_origD$testCovCountMat
        selectRegroup$selecEstMatIndv <- results_origD$testEstMat

      }

      goodIndpRefTaxWithCount <- results_origD$countOfSelecForAPred[1, (colnames(results_origD$countOfSelecForAPred) %in% goodRefTaxaCandi)]
      goodIndpRefTaxWithEst <- results_origD$estOfSelectForAPred[1, (colnames(results_origD$estOfSelectForAPred) %in% goodRefTaxaCandi)]

      restRefTaxWithCount <- results_origD$countOfSelecForAPred[1, !(colnames(results_origD$countOfSelecForAPred) %in% goodRefTaxaCandi)]
      restRefTaxWithEst <- results_origD$estOfSelectForAPred[1, !(colnames(results_origD$estOfSelectForAPred) %in% goodRefTaxaCandi)]

      sort_goodIndpRefTaxWithCount <- goodIndpRefTaxWithCount[order(goodIndpRefTaxWithCount, abs(goodIndpRefTaxWithEst))]
      sort_goodIndpRefTaxWithEst <- abs(goodIndpRefTaxWithEst[order(goodIndpRefTaxWithCount, abs(goodIndpRefTaxWithEst))])

      sort_restRefTaxWithCount <- restRefTaxWithCount[order(restRefTaxWithCount, abs(restRefTaxWithEst))]
      sort_restRefTaxWithEst <- abs(restRefTaxWithEst[order(restRefTaxWithCount, abs(restRefTaxWithEst))])


      selectRegroup$finalIndpRefTax <- names(c(sort_goodIndpRefTaxWithCount, sort_restRefTaxWithCount))[seq_len(2)]

      selectRegroup$selecCountOverall <- results_origD$countOfSelecForAPred

      selectRegroup$goodIndpRefTaxWithCount <- c(sort_goodIndpRefTaxWithCount, sort_restRefTaxWithCount)
      selectRegroup$goodIndpRefTaxWithEst <- c(sort_goodIndpRefTaxWithEst, sort_restRefTaxWithEst)

      selectRegroup$goodRefTaxaCandi <- goodRefTaxaCandi

      #### here ends the long section


      fin_ref_2 <- selectRegroup$finalIndpRefTax
      ref_taxa_2 <- selectRegroup$refTaxa
      while_loop_ind <- identical(fin_ref_1, fin_ref_2) || identical(ref_taxa_1, ref_taxa_2)

      if (while_loop_ind == FALSE) {
        message("Looping for DS IFAA part1 still ongoing.")
      }
      if (while_loop_ind == TRUE) {
        message("Looping for DS IFAA part1 done.")
      }
    }








  #### type = pooled ends below
  }

  return(selectRegroup)

}

# aggregate function



