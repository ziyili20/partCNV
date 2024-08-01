test_that("NormalizeCounts", {
  data(SimDataSce)
  counts <- NormalizeCounts(SimDataSce)
  
  expect_type(counts, "double")
})
