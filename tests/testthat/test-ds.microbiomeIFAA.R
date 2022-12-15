

#context("ds.microbiomeIFAA - errors")


connect.studies.dataset.cnsim(list("LAB_TSC", "GENDER", "PM_BMI_CATEGORICAL", "DIS_DIAB", "LAB_TRIG"))

test_that("ds.microbiomeIFAA errors", {

  #### subset of the microbiome and covariate data needs to be done prior to calculating the summarized experiment
  ds.dataFrame(x = c("D$LAB_TSC", "D$DIS_DIAB"), newobj = "microbiome_df1", datasources = ds.test_env$connections)
  ds.dataFrame(x = c("D$GENDER", "D$PM_BMI_CATEGORICAL"), newobj = "covariates_df1", datasources = ds.test_env$connections)
  ds.summarizedExperiment(df = "D", microbiomeData = "microbiome_df1", covariateData = "covariates_df1", newobj = "sumexp_df1")

  # Actual Test Start
  expect_error(ds.microbiomeIFAA(), "Please provide the name of the SummarizedExperiment object!", fixed = TRUE)
  expect_error(ds.microbiomeIFAA(SumExp = "sumexp_df1"), "Please provide name(s) of the covariate(s)!", fixed = TRUE)
  expect_error(ds.microbiomeIFAA(SumExp = "sumexp_df1", covariates = "GENDER"), "Please provide name(s) of the confounders!", fixed = TRUE)
  expect_error(ds.microbiomeIFAA(SumExp = "microbiome_df1", covariates = "GENDER", confounders = "PM_BMI_CATEGORICAL"), "Only objects of type 'SummarizedExperiment' are allowed", fixed = TRUE)
  expect_no_error(ds.microbiomeIFAA(SumExp = "sumexp_df1", covariates = "GENDER", confounders = "PM_BMI_CATEGORICAL"))

})


disconnect.studies.dataset.cnsim()

