test_that("GetCytoLocation", {
  expect1 <- GetCytoLocation(cyto_feature = "chr20(q11.1-q13.1)")
  expect2 <- GetCytoLocation(chr = "chr20", start = 25600000, end = 49800000)  
  expect_type(expect1, "list")
  expect_type(expect2, "list")
})
