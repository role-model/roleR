
test_that("use cases vignette compiles without error"){
    expect_error(rmarkdown::render("vignettes/roleR_use_cases.Rmd"), NA)
}