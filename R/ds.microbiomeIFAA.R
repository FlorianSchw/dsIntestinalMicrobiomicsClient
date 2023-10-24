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
    outcome_part3 <- DSI::datashield.aggregate(datasources, cally)



    #### k stands for selected microbiome


    XTX_comb_list <- list()
    Xy_comb_list <- list()
    yy_comb_list <- list()
    dfr_comb_list <- list()

    results <- list()
    results$estiList <- list()

    for (k in 1:length(outcome_part3[[1]]$analysisResults$estiList)){
      for (q in 1:length(datasources)){

        XTX_comb_list[[q]] <- outcome_part3[[q]][["analysisResults"]][["estiList"]][[k]][[1]][[1]][[1]]
        Xy_comb_list[[q]] <- outcome_part3[[q]][["analysisResults"]][["estiList"]][[k]][[1]][[1]][[2]]
        yy_comb_list[[q]] <- outcome_part3[[q]][["analysisResults"]][["estiList"]][[k]][[1]][[1]][[3]]
        dfr_comb_list[[q]] <- outcome_part3[[q]][["analysisResults"]][["estiList"]][[k]][[1]][[1]][[5]]
      }

      XTX_comb <- Reduce('+', XTX_comb_list)
      Xy_comb <- Reduce('+', Xy_comb_list)
      yy_comb <- Reduce('+', yy_comb_list)
      dfr_comb <- Reduce('+', dfr_comb_list)

      nvar <- outcome_part3[[1]][["analysisResults"]][["estiList"]][[k]][[1]][[1]][[4]]
      keep <- outcome_part3[[1]][["analysisResults"]][["estiList"]][[k]][[1]][[1]][[6]]
      xnames <- outcome_part3[[1]][["analysisResults"]][["estiList"]][[k]][[1]][[1]][[7]]
      intercept <- outcome_part3[[1]][["analysisResults"]][["estiList"]][[k]][[1]][[1]][[8]]
      originRefTaxNam <- outcome_part3[[1]][["analysisResults"]][["estiList"]][[k]][[1]][[1]][[9]]
      testCovInOrder <- outcome_part3[[1]][["analysisResults"]][["estiList"]][[k]][[1]][[1]][[10]]
      nRuns <- outcome_part3[[1]][["analysisResults"]][["estiList"]][[k]][[1]][[1]][[11]]
      taxa_sepname_list <- outcome_part3[[1]][["analysisResults"]][["estiList"]][[k]][[1]][[1]][[12]]
      nPredics <- outcome_part3[[1]][["analysisResults"]][["estiList"]][[k]][[1]][[1]][[13]]
      i <- outcome_part3[[1]][["analysisResults"]][["estiList"]][[k]][[1]][[1]][[14]]
      microbName <- outcome_part3[[1]][["analysisResults"]][["estiList"]][[k]][[1]][[1]][[15]]
      unbalanceTaxa_ori_name <- outcome_part3[[1]][["analysisResults"]][["estiList"]][[k]][[1]][[1]][[16]]
      unbalancePred_ori_name <- outcome_part3[[1]][["analysisResults"]][["estiList"]][[k]][[1]][[1]][[17]]
      fwerRate <- outcome_part3[[1]][["analysisResults"]][["estiList"]][[k]][[1]][[1]][[18]]

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


      bootResu_k <- coefMat

      fin_ref_taxon_name <- originRefTaxNam
      nTestcov <- length(testCovInOrder)


      boot_est <- bootResu_k[, 1] / nRuns
      se_est_all <- bootResu_k[, 2] / nRuns



      ref_taxon_name <- originRefTaxNam
      p_value_save_mat <- matrix(nrow = nTestcov, ncol = (length(taxa_sepname_list[[i]]) - 1))
      est_save_mat <- matrix(nrow = nTestcov, ncol = (length(taxa_sepname_list[[i]]) - 1))
      CI_up_mat <- matrix(nrow = nTestcov, ncol = (length(taxa_sepname_list[[i]]) - 1))
      CI_low_mat <- matrix(nrow = nTestcov, ncol = (length(taxa_sepname_list[[i]]) - 1))
      se_mat <- matrix(nrow = nTestcov, ncol = (length(taxa_sepname_list[[i]]) - 1))

      adjust_method = "BY"

      for (ii in seq_len(nTestcov)) {
        se_est <- se_est_all[seq(ii + 1, length(se_est_all), nPredics + 1)]
        boot_est_par <- boot_est[seq(ii + 1, length(boot_est), nPredics + 1)]

        p_value_unadj <- (1 - pnorm(abs(boot_est_par / se_est))) * 2

        boot_est_CI_low <- boot_est_par - 1.96 * se_est
        boot_est_CI_up <- boot_est_par + 1.96 * se_est

        p_value_adj <- p.adjust(p_value_unadj, adjust_method)

        p_value_save_mat[ii, ] <- p_value_unadj
        est_save_mat[ii, ] <- boot_est_par
        CI_low_mat[ii, ] <- boot_est_CI_low
        CI_up_mat[ii, ] <- boot_est_CI_up
        se_mat[ii, ] <- se_est

      }

      rownames(p_value_save_mat) <- testCovInOrder
      rownames(est_save_mat) <- testCovInOrder
      rownames(CI_low_mat) <- testCovInOrder
      rownames(CI_up_mat) <- testCovInOrder
      rownames(se_mat) <- testCovInOrder

      ref_taxon_name <- originRefTaxNam
      colname_use <- microbName[microbName != ref_taxon_name]
      colnames(p_value_save_mat) <- colname_use
      colnames(est_save_mat) <- colname_use
      colnames(CI_low_mat) <- colname_use
      colnames(CI_up_mat) <- colname_use
      colnames(se_mat) <- colname_use


      if (length(unbalanceTaxa_ori_name) > 0) {
        est_save_mat[cbind(unbalancePred_ori_name, unbalanceTaxa_ori_name)] <- NA
        CI_low_mat[cbind(unbalancePred_ori_name, unbalanceTaxa_ori_name)] <- NA
        CI_up_mat[cbind(unbalancePred_ori_name, unbalanceTaxa_ori_name)] <- NA
        se_mat[cbind(unbalancePred_ori_name, unbalanceTaxa_ori_name)] <- NA
        p_value_save_mat[cbind(unbalancePred_ori_name, unbalanceTaxa_ori_name)] <- NA
      }

      sig_ind <- which(p_value_save_mat < fwerRate, arr.ind = TRUE, useNames = FALSE)
      est_sig <- est_save_mat[sig_ind]
      CI_low_sig <- CI_low_mat[sig_ind]
      CI_up_sig <- CI_up_mat[sig_ind]
      p_adj_sig <- p_value_save_mat[sig_ind]
      se_sig <- se_mat[sig_ind]


      cov_sig_index <- sort(unique(sig_ind[, 1]))
      sig_list_each <- list()
      if (length(cov_sig_index) > 0) {
        for (iii in seq_len(length(cov_sig_index))) {
          sig_loc <- which(sig_ind[, 1] == cov_sig_index[iii])
          est_spe_cov <- est_sig[sig_loc]
          CI_low_spe_cov <- CI_low_sig[sig_loc]
          CI_up_spe_cov <- CI_up_sig[sig_loc]
          p_adj_spe_cov <- p_adj_sig[sig_loc]
          se_spe_cov <- se_sig[sig_loc]
          cov_sig_mat <- matrix(nrow = length(sig_loc), ncol = 5)
          colnames(cov_sig_mat) <-
            c("estimate", "SE est", "CI low", "CI up", "adj p-value")
          cov_sig_mat[, 1] <- est_spe_cov
          cov_sig_mat[, 2] <- se_spe_cov
          cov_sig_mat[, 3] <- CI_low_spe_cov
          cov_sig_mat[, 4] <- CI_up_spe_cov
          cov_sig_mat[, 5] <- p_adj_spe_cov
          rownames(cov_sig_mat) <- colname_use[sig_ind[sig_loc, 2]]
          sig_list_each[[testCovInOrder[cov_sig_index[iii]]]] <-
            cov_sig_mat
        }
      }
      all_cov_list <- list()


      all_cov_list$est_save_mat <- est_save_mat
      all_cov_list$p_value_save_mat <- p_value_save_mat
      all_cov_list$CI_low_mat <- CI_low_mat
      all_cov_list$CI_up_mat <- CI_up_mat
      all_cov_list$se_mat <- se_mat


      current_name <- names(outcome_part3[[1]]$analysisResults$estiList[k])



      results$estiList[[current_name]]$sig_list_each <- sig_list_each
      results$estiList[[current_name]]$all_cov_list <- all_cov_list


    }





    #### objects for regulariz + MZILN part of MZILN

    sub_taxa <- microbVar
    nSub <- outcome_part3[[1]]$analysisResults$nSub
    nTaxa <- outcome_part3[[1]]$analysisResults$nTaxa
    #covariatesData <- outcome_part3[[1]]$covariatesData
    linkIDname <- outcome_part3[[1]]$linkIDname


    #### regulariz starts here
    #### there need to be a few changes

    fin_ref_taxon_name <- names(results$estiList)
    all_cov_list_sep <- list()

    for (i in seq_len(length(fin_ref_taxon_name))) {
      save_list_temp <- results$estiList[[i]]$all_cov_list
      rearrage_res_list <- list()
      for (j in testCovInOrder) {
        est_res_save_all <-
          cbind(
            save_list_temp$est_save_mat[j,],
            save_list_temp$se_mat[j,],
            save_list_temp$CI_low_mat[j,],
            save_list_temp$CI_up_mat[j,],
            save_list_temp$p_value_save_mat[j,]
          )
        est_res_save_all <-
          data.frame(fin_ref_taxon_name[i],
                     rownames(est_res_save_all),
                     j,
                     est_res_save_all)

        colnames(est_res_save_all) <-
          c("ref_tax",
            "taxon",
            "cov",
            "estimate",
            "SE est",
            "CI low",
            "CI up",
            "adj p-value")

        rearrage_res_list[[j]] <- est_res_save_all

      }

      unorder_long <- DescTools::DoCall("rbind", rearrage_res_list)
      all_cov_list_sep[[fin_ref_taxon_name[i]]] <- data.frame(unorder_long[stringr::str_order(unorder_long[, c("taxon")], decreasing = FALSE, numeric = TRUE),], row.names = NULL)
    }


    all_cov_list_sep <- DescTools::DoCall("rbind", all_cov_list_sep)
    rownames(all_cov_list_sep) <- NULL
    all_cov_list_sep$sig_ind <- all_cov_list_sep$adj.p.value < fwerRate
    results$all_cov_list_sep <- S4Vectors::DataFrame(all_cov_list_sep)



    all_cov_list <- list()

    for (i in seq_len(length(fin_ref_taxon_name))) {
      all_cov_list[[fin_ref_taxon_name[i]]] <-
        results$estiList[[i]]$all_cov_list
    }
    results$all_cov_list <- all_cov_list


    ref_taxon_name <- names(all_cov_list)

    exclu_1 <- !colnames(all_cov_list[[1]]$est_save_mat) %in% ref_taxon_name[2]
    exclu_2 <- !colnames(all_cov_list[[2]]$est_save_mat) %in% ref_taxon_name[1]



    #### hard coded mean estimations from the original!
    est_save_mat_mean <- (all_cov_list[[1]]$est_save_mat[, exclu_1, drop = FALSE] + all_cov_list[[2]]$est_save_mat[, exclu_2, drop = FALSE]) / 2
    se_mat_mean <- (all_cov_list[[1]]$se_mat[, exclu_1, drop = FALSE] + all_cov_list[[2]]$se_mat[, exclu_2, drop = FALSE]) / 2
    CI_low_mat_mean <- (all_cov_list[[1]]$CI_low_mat[, exclu_1, drop = FALSE] + all_cov_list[[2]]$CI_low_mat[, exclu_2, drop = FALSE]) / 2
    CI_up_mat_mean <- (all_cov_list[[1]]$CI_up_mat[, exclu_1, drop = FALSE] + all_cov_list[[2]]$CI_up_mat[, exclu_2, drop = FALSE]) / 2


    p_value_unadj_mean <- t(apply(est_save_mat_mean / se_mat_mean, 1, function(x) {
        (1 - pnorm(abs(x))) * 2
      }))
    p_value_adj_mean <- t(apply(p_value_unadj_mean, 1, function(x) {
        p.adjust(x, method = adjust_method)
      }))
    colname_use <- colnames(est_save_mat_mean)


    sig_ind <- which(p_value_adj_mean <= fwerRate, arr.ind = TRUE, useNames = FALSE)


    full_results <- list()
    for (j in testCovInOrder) {
      est_res_save_all <-
        data.frame(
          colname_use,
          j,
          est_save_mat_mean[j, ],
          se_mat_mean[j, ],
          CI_low_mat_mean[j, ],
          CI_up_mat_mean[j, ],
          p_value_unadj_mean[j, ],
          p_value_adj_mean[j, ],
          row.names = NULL
        )
      colnames(est_res_save_all) <-
        c(
          "taxon",
          "cov",
          "estimate",
          "SE est",
          "CI low",
          "CI up",
          "unadj p-value",
          "adj p-value"
        )
      fin_ref_taxon_dat <-
        data.frame(
          taxon = fin_ref_taxon_name,
          cov = j,
          estimate = 0,
          SE.est = NA,
          CI.low = NA,
          CI.up = NA,
          unadj.p.value = 1,
          adj.p.value = 1
        )

      res <- rbind(data.frame(est_res_save_all), fin_ref_taxon_dat)

      if (any(sub_taxa!="all")) {
        res <- res[res$taxon %in% sub_taxa, , drop = FALSE]
        res$adj.p.value <- p.adjust(res$unadj.p.value, method = adjust_method)
      }

      full_results[[j]] <- res
    }



    full_results <- DescTools::DoCall("rbind", full_results)
    rownames(full_results) <- NULL
    full_results <- full_results[stringr::str_order(full_results$taxon,
                                                    decreasing = FALSE,
                                                    numeric = TRUE), ]
    full_results$sig_ind <- full_results$adj.p.value < fwerRate
    full_results$sig_ind[is.na(full_results$sig_ind)] <- FALSE
    results$full_results <- S4Vectors::DataFrame(full_results)

    results$nTaxa <- nTaxa
    results$nPredics <- nPredics
    results$fin_ref_taxon_name <- fin_ref_taxon_name

    rm(data)

    # return results

    results$nRef <- nRef



    output_obj <- as.data.frame(results[["full_results"]]@listData)
    output_obj$confounder <- list(c(ctrlCov))
    output_obj$RefTaxa <- list(c(refTaxa))
    output_obj$nTaxa <- results[["nTaxa"]]
    output_obj$nPredics <- results[["nPredics"]]














  #### type = pooled ends below
  }

  return(outcome_part3)

}

# aggregate function



