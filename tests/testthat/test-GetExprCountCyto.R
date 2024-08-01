test_that("GetExprCountCyto", {
  data(SimData)
  res <- GetCytoLocation(cyto_feature = "chr20(q11.1-q13.1)")
  expect <- GetExprCountCyto(cytoloc_output = res, Counts = as.matrix(SimData), normalization = TRUE, qt_cutoff = 0.99)
  
  expect_type(res, "list")
  expect_type(expect, "list")
})
