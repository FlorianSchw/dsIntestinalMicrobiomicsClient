




connect.studies.dataset.cnsim(list("LAB_TSC", "GENDER", "PM_BMI_CATEGORICAL", "DIS_DIAB", "LAB_TRIG"))

test_that("ds.summarizedExperiment errors", {


  #### subset of the microbiome and covariate data needs to be done prior to calculating the summarized experiment
  ds.dataFrame(x = c("D$LAB_TSC", "D$DIS_DIAB"), newobj = "microbiome_df1", datasources = ds.test_env$connections)
  ds.dataFrame(x = c("D$GENDER", "D$PM_BMI_CATEGORICAL"), newobj = "covariates_df1", datasources = ds.test_env$connections)


  # Actual Test Start
  expect_error(ds.summarizedExperiment(), "Please provide the name of the data.frame!", fixed = TRUE)
  expect_error(ds.summarizedExperiment(df = "D"), "Please provide the name(s) of the microbiome data!", fixed = TRUE)
  expect_error(ds.summarizedExperiment(df = "D", microbiomeData = "microbiome_df1"), "Please provide the name(s) of the covariate data!", fixed = TRUE)
  expect_silent(ds.summarizedExperiment(df = "D", microbiomeData = "microbiome_df1", covariateData = "covariates_df1"))
  expect_error(ds.summarizedExperiment(df = "D$LAB_TSC", microbiomeData = "microbiome_df1", covariateData = "covariates_df1"), "Only objects of type 'data frame' are allowed.", fixed = TRUE)

})


disconnect.studies.dataset.cnsim()

