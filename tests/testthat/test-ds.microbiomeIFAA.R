

test_that("ds.microbiomeIFAA Errors", {



  # Create data.frames with microbiome and covariate data of interest
    ds.dataFrame(x = c("D$P_ACTINOBACTERIA",
                       "D$P_BACTEROIDETES",
                       "D$P_FIRMICUTES",
                       "D$P_VERRUCOMICROBIA"),
                 newobj = "microbdata",
                 stringsAsFactors = FALSE)

    ds.dataFrame(x = c("D$Age",
                       "D$Sex",
                       "D$Weight",
                       "D$Height"),
                 newobj = "covdata",
                 stringsAsFactors = FALSE)


    # Create the summarizedExperiment object on the server-side based on the microbiome and covariate data.frames

    ds.summarizedExperiment(microbiomeData = "microbdata",
                            covariateData = "covdata",
                            newobj = "SumExpT")

    # Calculate the associations of covariates with the microbiome ratios

    results <- ds.microbiomeIFAA(SumExp = "SumExpT",
                                 microbVar = "P_BACTEROIDETES",
                                 refTaxa = "P_VERRUCOMICROBIA",
                                 testCov = "Weight",
                                 ctrlCov = c("Age", "Sex"),
                                 adjust_method = "BY",
                                 fdrRate = 0.05,
                                 type = "both")

  expect_silent(results)

})


