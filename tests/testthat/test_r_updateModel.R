# # tests for updateLocalComm, updateMetaComm, updateRolePhylo ----
# 
# # helpers
# make_model <- function(sm = 3, J = 10, n = 100) {
#     roleModel(roleParams(species_meta = sm,
#                          individuals_local = J,
#                          individuals_meta = 100,
#                          niter = n,
#                          niterTimestep = n))
# }
# 
# 
# # assertPreRun blocks post-run models ----
# 
# test_that("update functions reject post-run models", {
#     m <- make_model()
#     mr <- runRole(m)
#     expect_error(updateLocalComm(mr, rep(1L, 10), slot = "indSpecies"),
#                  "pre-run")
#     expect_error(updateMetaComm(mr, rep(1/3, 3), slot = "spAbund"),
#                  "pre-run")
# })
# 
# 
# # updateLocalComm: slot = "indSpecies" ----
# 
# test_that("updateLocalComm indSpecies rejects wrong length", {
#     m <- make_model(sm = 3, J = 10)
#     expect_error(updateLocalComm(m, rep(1L, 5), slot = "indSpecies"),
#                  "individuals_local")
# })
# 
# test_that("updateLocalComm indSpecies rejects out-of-range IDs", {
#     m <- make_model(sm = 3, J = 10)
#     expect_error(updateLocalComm(m, rep(4L, 10), slot = "indSpecies"),
#                  "outside")
# })
# 
# test_that("updateLocalComm indSpecies syncs spAbund", {
#     m <- make_model(sm = 3, J = 10)
#     m <- updateLocalComm(m, c(rep(1L, 5), rep(2L, 5)), slot = "indSpecies")
#     expect_equal(m@modelSteps[[1]]@localComm@spAbund, c(5, 5, 0))
# })
# 
# test_that("updateLocalComm indSpecies syncs spAbundHarmMean", {
#     m <- make_model(sm = 3, J = 10)
#     m <- updateLocalComm(m, c(rep(1L, 7), rep(3L, 3)), slot = "indSpecies")
#     expect_equal(m@modelSteps[[1]]@localComm@spAbundHarmMean, c(7, 0, 3))
# })
# 
# test_that("updateLocalComm indSpecies syncs spPropHarmMean", {
#     m <- make_model(sm = 2, J = 10)
#     m <- updateLocalComm(m, c(rep(1L, 4), rep(2L, 6)), slot = "indSpecies")
#     expect_equal(m@modelSteps[[1]]@localComm@spPropHarmMean, c(0.4, 0.6))
# })
# 
# test_that("updateLocalComm indSpecies sets spLastOriginStep to zero", {
#     m <- make_model(sm = 3, J = 10)
#     m <- updateLocalComm(m, c(rep(1L, 5), rep(2L, 5)), slot = "indSpecies")
#     expect_equal(m@modelSteps[[1]]@localComm@spLastOriginStep, c(0, 0, 0))
# })
# 
# test_that("updateLocalComm indSpecies sets spLastOriginGen to zero", {
#     m <- make_model(sm = 3, J = 10)
#     m <- updateLocalComm(m, c(rep(1L, 5), rep(2L, 5)), slot = "indSpecies")
#     expect_equal(m@modelSteps[[1]]@localComm@spLastOriginGen, c(0, 0, 0))
# })
# 
# test_that("updateLocalComm indSpecies produces species-level vectors of length sm", {
#     sm <- 4
#     m <- make_model(sm = sm, J = 20)
#     # only two species present; all spp-level vectors must still be length sm
#     m <- updateLocalComm(m, c(rep(1L, 10), rep(2L, 10)), slot = "indSpecies")
#     lc <- m@modelSteps[[1]]@localComm
#     expect_length(lc@spAbund, sm)
#     expect_length(lc@spAbundHarmMean, sm)
#     expect_length(lc@spPropHarmMean, sm)
#     expect_length(lc@spLastOriginStep, sm)
#     expect_length(lc@spLastOriginGen, sm)
# })
# 
# 
# # updateLocalComm: slot = "indTrait" ----
# 
# test_that("updateLocalComm indTrait rejects non-matrix", {
#     m <- make_model(sm = 2, J = 10)
#     expect_error(updateLocalComm(m, rep(0, 10), slot = "indTrait"),
#                  "matrix")
# })
# 
# test_that("updateLocalComm indTrait rejects wrong nrow", {
#     m <- make_model(sm = 2, J = 10)
#     expect_error(updateLocalComm(m, matrix(0, nrow = 5, ncol = 1), slot = "indTrait"),
#                  "individuals_local")
# })
# 
# test_that("updateLocalComm indTrait updates successfully", {
#     m <- make_model(sm = 2, J = 10)
#     trt <- matrix(seq(-1, 1, length.out = 10), ncol = 1)
#     m <- updateLocalComm(m, trt, slot = "indTrait")
#     expect_equal(m@modelSteps[[1]]@localComm@indTrait, trt)
# })
# 
# 
# # updateMetaComm: slot = "spAbund" ----
# 
# test_that("updateMetaComm spAbund rejects wrong length", {
#     m <- make_model(sm = 3)
#     expect_error(updateMetaComm(m, c(0.5, 0.5), slot = "spAbund"),
#                  "species_meta")
# })
# 
# test_that("updateMetaComm spAbund rejects values that don't sum to 1", {
#     m <- make_model(sm = 3)
#     expect_error(updateMetaComm(m, c(0.5, 0.5, 0.5), slot = "spAbund"),
#                  "sum")
# })
# 
# test_that("updateMetaComm spAbund rejects negative values", {
#     m <- make_model(sm = 3)
#     expect_error(updateMetaComm(m, c(-0.1, 0.6, 0.5), slot = "spAbund"),
#                  "non-negative")
# })
# 
# test_that("updateMetaComm spAbund updates successfully", {
#     m <- make_model(sm = 3)
#     ab <- c(0.2, 0.3, 0.5)
#     m <- updateMetaComm(m, ab, slot = "spAbund")
#     expect_equal(m@modelSteps[[1]]@metaComm@spAbund, ab)
# })
# 
# 
# # updateMetaComm: slot = "spTrait" ----
# 
# test_that("updateMetaComm spTrait rejects non-matrix", {
#     m <- make_model(sm = 3)
#     expect_error(updateMetaComm(m, c(0, 1, 2), slot = "spTrait"),
#                  "matrix")
# })
# 
# test_that("updateMetaComm spTrait rejects wrong nrow", {
#     m <- make_model(sm = 3)
#     expect_error(updateMetaComm(m, matrix(0, nrow = 2, ncol = 1), slot = "spTrait"),
#                  "species_meta")
# })
# 
# test_that("updateMetaComm spTrait updates successfully", {
#     m <- make_model(sm = 3)
#     trt <- matrix(c(-5, 0, 5), ncol = 1)
#     m <- updateMetaComm(m, trt, slot = "spTrait")
#     expect_equal(m@modelSteps[[1]]@metaComm@spTrait, trt)
# })
# 
# 
# # consistency: after updateLocalComm indSpecies, model runs without error ----
# 
# test_that("model with updated indSpecies runs without error", {
#     sm <- 2
#     J  <- 20
#     m  <- make_model(sm = sm, J = J, n = 100)
#     m  <- updateLocalComm(m, c(rep(1L, J / 2), rep(2L, J / 2)), slot = "indSpecies")
#     expect_no_error(runRole(m))
# })
