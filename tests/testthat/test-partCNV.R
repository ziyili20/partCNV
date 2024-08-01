test_that("partCNV", {
  data(SimData)
  cytoloc <- GetCytoLocation(cyto_feature = "chr20(q11.1-q13.1)")
  exprout <- GetExprCountCyto(cytoloc_output = cytoloc, Counts = as.matrix(SimData), normalization = TRUE, qt_cutoff = 0.99)
  status <- partCNV(int_counts = exprout$ProcessedCount, cyto_type = "del", cyto_p = 0.2)
  
  expect_type(exprout, "list")
  expect_type(status, "double")
})
