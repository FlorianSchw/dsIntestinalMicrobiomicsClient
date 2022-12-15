

#context("ds.dist - errors")


connect.studies.dataset.cnsim(list("LAB_TSC", "GENDER", "PM_BMI_CATEGORICAL", "DIS_DIAB", "LAB_TRIG"))

test_that("ds.dist errors", {

  #### subset of the microbiome and covariate data needs to be done prior to calculating the summarized experiment
  ds.dataFrame(x = c("D$LAB_TSC", "D$DIS_DIAB"), newobj = "microbiome_df1", datasources = ds.test_env$connections)
  ds.dataFrame(x = c("D$GENDER", "D$PM_BMI_CATEGORICAL"), newobj = "covariates_df1", datasources = ds.test_env$connections)
  ds.summarizedExperiment(df = "D", microbiomeData = "microbiome_df1", covariateData = "covariates_df1", newobj = "sumexp_df1")

  # Actual Test Start
  expect_error(ds.microbiomeMZILN(), "Please provide the name of the SummarizedExperiment object!", fixed = TRUE)
  expect_error(ds.microbiomeMZILN(SumExp = "sumexp_df1"), "Please provide the name(s) of the microbiome denominator data!", fixed = TRUE)
  expect_error(ds.microbiomeMZILN(SumExp = "microbiome_df1", taxa = "LAB_TSC"), "Only objects of type 'SummarizedExperiment' are allowed", fixed = TRUE)
  expect_no_error(ds.microbiomeMZILN(SumExp = "sumexp_df1", taxa = "LAB_TSC"))
  expect_no_error(ds.microbiomeMZILN(SumExp = "sumexp_df1", taxa = "LAB_TSC", allCov = "GENDER"))

})


disconnect.studies.dataset.cnsim()

