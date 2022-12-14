




connect.studies.dataset.cnsim(list("LAB_TSC", "GENDER", "PM_BMI_CATEGORICAL", "DIS_DIAB", "LAB_TRIG"))

test_that("ds.summarizedExperiment errors", {

  # Creating differing test data.frames in the respective connections for testing if columns check works correctly


  #### a potential microbiome data needs to be simulated

  #### the covariate data need to be created

  #### it could work by using the ds. data frame subset function

  #### subset of the microbiome and covariate data needs to be done prior to calculating the summarized experiment



  ds.dataFrame(x = c("D$LAB_TSC", "D$DIS_DIAB"), newobj = "microbiome_df1", datasources = ds.test_env$connections)
  ds.dataFrame(x = c("D$GENDER", "D$PM_BMI_CATEGORICAL"), newobj = "covariates_df1", datasources = ds.test_env$connections)






  ds.dataFrame(x = c("D$LAB_TSC", "D$LAB_TRIG", "D$PM_BMI_CATEGORICAL", "D$DIS_DIAB"), newobj = "test_df2", datasources = ds.test_env$connections[1])
  ds.dataFrame(x = c("D$LAB_TSC", "D$PM_BMI_CATEGORICAL", "D$LAB_TRIG", "D$DIS_DIAB"), newobj = "test_df2", datasources = ds.test_env$connections[2])
  ds.dataFrame(x = c("D$LAB_TSC", "D$PM_BMI_CATEGORICAL", "D$LAB_TRIG", "D$DIS_DIAB"), newobj = "test_df2", datasources = ds.test_env$connections[3])


  ds.dataFrame(x = c("D$LAB_TSC", "D$PM_BMI_CATEGORICAL", "D$LAB_TRIG", "D$DIS_DIAB"), newobj = "test_df3", datasources = ds.test_env$connections[1])
  ds.dataFrame(x = c("D$LAB_TSC", "D$PM_BMI_CATEGORICAL", "D$LAB_TRIG"), newobj = "test_df3", datasources = ds.test_env$connections[2])
  ds.dataFrame(x = c("D$LAB_TSC", "D$PM_BMI_CATEGORICAL", "D$LAB_TRIG"), newobj = "test_df3", datasources = ds.test_env$connections[3])


  ds.asCharacter(x.name = "D$LAB_TSC", newobj = "char.obj", datasources = ds.test_env$connections[1])
  ds.asNumeric(x.name = "D$LAB_TSC", newobj = "char.obj", datasources = ds.test_env$connections[2])
  ds.asNumeric(x.name = "D$LAB_TSC", newobj = "char.obj", datasources = ds.test_env$connections[3])
  ds.dataFrame(x = c("char.obj", "D$GENDER", "D$PM_BMI_CATEGORICAL", "D$DIS_DIAB"), newobj = "test_df4", datasources = ds.test_env$connections)


  ds.completeCases("D", newobj = "D_clean", datasources = ds.test_env$connections)

  ds.dataFrame(x = c("D$LAB_TSC", "D$DIS_DIAB"), newobj = "test_df6", datasources = ds.test_env$connections[1])

  # Actual Test Start
  expect_error(ds.summarizedExperiment(), "Please provide the name of the data.frame!", fixed = TRUE)
  expect_error(ds.summarizedExperiment(df = "D"), "Please provide the name(s) of the microbiome data!", fixed = TRUE)
  expect_error(ds.summarizedExperiment(df = "D", microbiomeData = "microbiome_df1"), "The data frames do not have the same columns. There are columns missing in some data frames!", fixed = TRUE)
  expect_silent(ds.summarizedExperiment("test_df2"))
  expect_error(ds.summarizedExperiment("test_df3"), "The data frames do not have the same columns. There are columns missing in some data frames!", fixed = TRUE)
  expect_error(ds.summarizedExperiment("test_df4"), "The data frames contain columns which are not of type 'numeric' or 'integer'.", fixed = TRUE)
  expect_error(ds.summarizedExperiment("D_clean", method = "euclidea"), "Method needs to be one of the following: 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary' or 'minkowski'.", fixed = TRUE)
  expect_error(ds.summarizedExperiment("D_clean", method = "x"), "Method needs to be one of the following: 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary' or 'minkowski'.", fixed = TRUE)
  expect_error(ds.summarizedExperiment("test_df6"))
})


disconnect.studies.dataset.cnsim()

