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
#' @param datasources a list of \code{\link{DSConnection-class}} objects obtained after login
#' @return \code{ds.microbiomeMZILN} returns the outcome of the specified multivariate zero-inflated logistic normal model
#' @author Florian Schwarz for the German Institute of Human Nutrition
#' @import DSI
#' @import dsBaseClient
#' @import methods
#' @import Matrix
#' @import S4Vectors
#' @import MatrixExtra
#' @import DescTools
#' @import stringr
#' @export
#'

ds.microbiomeMZILN <- function(SumExp = NULL,
                               microbVar = NULL,
                               refTaxa = NULL,
                               allCov = NULL,
                               sampleIDname = NULL,
                               adjust_method = "BY",
                               fdrRate = 0.05,
                               paraJobs = NULL,
                               bootB = 500,
                               taxDropThresh = 0,
                               standardize = FALSE,
                               verbose = TRUE,
                               type = c("split", "pooled"),
                               datasources = NULL){


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




  # call the internal function that checks the input object is of the same class in all studies.
  typ <- dsBaseClient::ds.class(SumExp, datasources)



  # Check whether the input is either of type data frame or matrix
  if(!('SummarizedExperiment' %in% typ)){
    stop("Only objects of type 'SummarizedExperiment' are allowed.", call.=FALSE)
  }


  if(!(is.null(microbVar))){
    microbVar_ds <- paste0(microbVar, collapse = ",")
  } else {
    microbVar_ds <- microbVar
  }

  if(!(is.null(refTaxa))){
    refTaxa_ds <- paste0(refTaxa, collapse = ",")
  } else {
    refTaxa_ds <- refTaxa
  }

  if(!(is.null(allCov))){
    allCov_ds <- paste0(allCov, collapse = ",")
  } else {
    allCov_ds <- allCov
  }







  # Check whether the input is either of type data frame or matrix
  if(type == "split"){

    start.time <- proc.time()[3]


    # call the server side function that does the operation
    cally <- call("microbiomeMZILNDS", SumExp, microbVar_ds, refTaxa_ds, allCov_ds, sampleIDname, adjust_method, fdrRate, paraJobs, taxDropThresh, standardize, verbose)
    output <- DSI::datashield.aggregate(datasources, cally)

    TotalTime_Split <- (proc.time()[3] - start.time) / 60

    for (t in 1:length(datasources)){

      output[[t]]$Analysis <- datasources[[t]]@name

    }

    output_obj <- dplyr::bind_rows(output)
    output_obj$Time_min <- TotalTime_Split

  }



  if(type == "pooled"){

    start.time <- proc.time()[3]


    # call the server side function that does the operation
      cally <- call("microbiomeMZILNPooledDS", SumExp, microbVar_ds, refTaxa_ds, allCov_ds, sampleIDname, adjust_method, fdrRate, paraJobs, taxDropThresh, standardize, verbose)
      outcome <- DSI::datashield.aggregate(datasources, cally)




      #### from different section
      intermediate_results <- list()

      #### k stands for selected microbiome


      XTX_comb_list <- list()
      Xy_comb_list <- list()
      yy_comb_list <- list()
      dfr_comb_list <- list()

      results <- list()
      results$estiList <- list()

      for (k in 1:length(outcome[[1]]$analysisResults$estiList)){
        for (q in 1:length(datasources)){

          XTX_comb_list[[q]] <- outcome[[q]][["analysisResults"]][["estiList"]][[k]][[1]][[1]][[1]]
          Xy_comb_list[[q]] <- outcome[[q]][["analysisResults"]][["estiList"]][[k]][[1]][[1]][[2]]
          yy_comb_list[[q]] <- outcome[[q]][["analysisResults"]][["estiList"]][[k]][[1]][[1]][[3]]
          dfr_comb_list[[q]] <- outcome[[q]][["analysisResults"]][["estiList"]][[k]][[1]][[1]][[5]]
        }

        XTX_comb <- Reduce('+', XTX_comb_list)
        Xy_comb <- Reduce('+', Xy_comb_list)
        yy_comb <- Reduce('+', yy_comb_list)
        dfr_comb <- Reduce('+', dfr_comb_list)

        nvar <- outcome[[1]][["analysisResults"]][["estiList"]][[k]][[1]][[1]][[4]]
        keep <- outcome[[1]][["analysisResults"]][["estiList"]][[k]][[1]][[1]][[6]]
        xnames <- outcome[[1]][["analysisResults"]][["estiList"]][[k]][[1]][[1]][[7]]
        intercept <- outcome[[1]][["analysisResults"]][["estiList"]][[k]][[1]][[1]][[8]]
        originRefTaxNam <- outcome[[1]][["analysisResults"]][["estiList"]][[k]][[1]][[1]][[9]]
        testCovInOrder <- outcome[[1]][["analysisResults"]][["estiList"]][[k]][[1]][[1]][[10]]
        nRuns <- outcome[[1]][["analysisResults"]][["estiList"]][[k]][[1]][[1]][[11]]
        taxa_sepname_list <- outcome[[1]][["analysisResults"]][["estiList"]][[k]][[1]][[1]][[12]]
        nPredics <- outcome[[1]][["analysisResults"]][["estiList"]][[k]][[1]][[1]][[13]]
        i <- outcome[[1]][["analysisResults"]][["estiList"]][[k]][[1]][[1]][[14]]
        microbName <- outcome[[1]][["analysisResults"]][["estiList"]][[k]][[1]][[1]][[15]]
        unbalanceTaxa_ori_name <- outcome[[1]][["analysisResults"]][["estiList"]][[k]][[1]][[1]][[16]]
        unbalancePred_ori_name <- outcome[[1]][["analysisResults"]][["estiList"]][[k]][[1]][[1]][[17]]
        fwerRate <- outcome[[1]][["analysisResults"]][["estiList"]][[k]][[1]][[1]][[18]]

        dfr_corr <- XTX_comb@Dim[1]


        coef_ds   <- Matrix::solve(XTX_comb, Xy_comb, tol = 1e-7)

        coefficients_ds <- rep(NA, nvar)
        coefficients_ds[as.vector(keep)] <- coef_ds

        RSS_ds <- yy_comb - 2 * MatrixExtra::crossprod(coef_ds, Xy_comb) + MatrixExtra::crossprod(coef_ds, MatrixExtra::crossprod(XTX_comb, coef_ds))

        var_res_ds <- as.numeric(RSS_ds)/(dfr_comb + (length(datasources)-1)*dfr_corr)

        se_coef_ds <- rep(NA, nvar)
        inv_ds     <- Matrix::solve(XTX_comb, diag(nrow(XTX_comb)), tol = 1e-7)

        se_coef_ds[keep] <- sqrt(var_res_ds * Matrix::diag(inv_ds))
        t1            <- coefficients_ds/se_coef_ds
        p             <- 2 * pt(abs(t1), df = dfr_comb + (length(datasources)-1)*dfr_corr, lower.tail = FALSE)



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
          est_save_mat[cbind(unbalancePred_ori_name, unbalanceTaxa_ori_name)] <-
            NA
          CI_low_mat[cbind(unbalancePred_ori_name, unbalanceTaxa_ori_name)] <-
            NA
          CI_up_mat[cbind(unbalancePred_ori_name, unbalanceTaxa_ori_name)] <-
            NA
          se_mat[cbind(unbalancePred_ori_name, unbalanceTaxa_ori_name)] <-
            NA
          p_value_save_mat[cbind(unbalancePred_ori_name, unbalanceTaxa_ori_name)] <-
            NA
        }

        sig_ind <-
          which(p_value_save_mat < fwerRate,
                arr.ind = TRUE,
                useNames = FALSE
          )
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


        current_name <- names(outcome[[1]]$analysisResults$estiList[k])



        results$estiList[[current_name]]$sig_list_each <- sig_list_each
        results$estiList[[current_name]]$all_cov_list <- all_cov_list


      }





        #### objects for regulariz + MZILN part of MZILN

        sub_taxa <- microbVar
        nSub <- outcome[[1]]$analysisResults$nSub
        nTaxa <- outcome[[1]]$analysisResults$nTaxa
        #covariatesData <- outcome[[1]]$covariatesData
        linkIDname <- outcome[[1]]$linkIDname


        #### regulariz_MZILN starts here

        fin_ref_taxon_name <- names(results$estiList)
        all_cov_sig_list <- list()
        all_cov_list <- list()

        for (i in seq_len(length(fin_ref_taxon_name))) {
          all_cov_sig_list[[fin_ref_taxon_name[i]]] <-
            results$estiList[[i]]$sig_list_each
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
                "unadj p-value")

            if (any(sub_taxa!="all")) {
              est_res_save_all <-
                est_res_save_all[est_res_save_all$taxon %in% sub_taxa, , drop = FALSE]
            }

            est_res_save_all$adj.p.value <-
              p.adjust(est_res_save_all$`unadj p-value`, adjust_method)
            rearrage_res_list[[j]] <- est_res_save_all
          }

          unorder_long <- do.call("rbind", rearrage_res_list)
          all_cov_list[[fin_ref_taxon_name[i]]] <-
            data.frame(unorder_long[stringr::str_order(unorder_long[, c("taxon")], decreasing = FALSE, numeric = TRUE),], row.names = NULL)
        }

        all_cov_list <- do.call("rbind", all_cov_list)
        rownames(all_cov_list) <- NULL
        all_cov_list$sig_ind <- all_cov_list$adj.p.value < fdrRate
        all_cov_list$sig_ind[is.na(all_cov_list$sig_ind)] <- FALSE




        results$full_results <- S4Vectors::DataFrame(all_cov_list)
        results$nSub <- nSub
        results$nTaxa <- nTaxa
        results$nPredics <- nPredics




        totalTimeMins <- (proc.time()[3] - start.time) / 60
        message("The entire analysis took ", round(totalTimeMins, 2), " minutes")


        output_obj_orig <- list(full_results = results$full_results,
                                metadata = list(totalTimeMins = totalTimeMins,
                                                fdrRate = fdrRate,
                                                adjust_method = adjust_method))



        #### Transforming the output object from original type to desired data.frame

        RefTaxa <- unlist(output_obj_orig[[1]][1])
        Taxa <- unlist(output_obj_orig[[1]][2])
        Covariate <- unlist(output_obj_orig[[1]][3])
        Estimate <- unlist(output_obj_orig[[1]][4])
        Std.Error <- unlist(output_obj_orig[[1]][5])
        CI_lower <- unlist(output_obj_orig[[1]][6])
        CI_upper <- unlist(output_obj_orig[[1]][7])
        p.value_adj <- unlist(output_obj_orig[[1]][8])
        p.value_unadj <- unlist(output_obj_orig[[1]][9])
        Significance <- unlist(output_obj_orig[[1]][10])
        Time_min <- unlist(output_obj_orig[[2]][1])
        FDR <- unlist(output_obj_orig[[2]][2])
        Adjustment_Method <- unlist(output_obj_orig[[2]][3])
        Analysis <- "pooled"


        output_obj <- data.frame(Taxa,
                                 RefTaxa,
                                 Covariate,
                                 Estimate,
                                 Std.Error,
                                 CI_lower,
                                 CI_upper,
                                 p.value_unadj,
                                 p.value_adj,
                                 FDR,
                                 Significance,
                                 Adjustment_Method,
                                 Analysis,
                                 Time_min)


  }

  return(output_obj)

}

# aggregate function


