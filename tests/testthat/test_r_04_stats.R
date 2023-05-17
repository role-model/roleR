test_that("getSumStats defaults return several metrics", {
    model <- quickModel()
    stats <- getSumStats(model)
    expect_true(ncol(stats) > 1)
})
