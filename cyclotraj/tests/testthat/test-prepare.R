test_that("prepare_clone_data standardizes columns and levels", {
  d <- prepare_clone_data(read_example_data(), treatment_ref = "arm_A")
  expect_s3_class(d, "cyclotraj_data")
  expect_true(all(c("clone_id", "patient_id", "treatment", "timepoint_num",
                    "timepoint_cat", "count", "total", "perc", "trajectory")
                  %in% names(d)))
  expect_identical(levels(d$treatment)[1], "arm_A")
  expect_false("Other" %in% levels(d$trajectory))
  expect_true(all(d$count <= d$total))
})

test_that("missing required columns raise an informative error", {
  raw <- read_example_data()
  raw$n <- NULL
  expect_error(prepare_clone_data(raw), "Missing required column")
})

test_that("count > total is rejected", {
  raw <- read_example_data()
  raw$n[1] <- raw$n_tcells_tp[1] + 1
  expect_error(prepare_clone_data(raw), "cannot exceed")
})

test_that("custom column mapping works on a generic data set", {
  generic <- data.frame(
    unit = rep(c("u1", "u2"), each = 3),
    subj = "s1",
    arm = rep(c("A", "B"), each = 3),
    wk = rep(c(0, 1, 2), 2),
    y = c(1, 2, 3, 0, 1, 0),
    N = 100
  )
  d <- prepare_clone_data(generic, clone = "unit", patient = "subj",
                          treatment = "arm", timepoint = "wk",
                          count = "y", total = "N", trajectory = NULL)
  expect_s3_class(d, "cyclotraj_data")
  expect_false("trajectory" %in% names(d))
  expect_equal(nlevels(d$treatment), 2)
})

test_that("baseline is the reference level and treatment_ref is validated", {
  ## the difference-in-differences estimator is only interpretable if the
  ## earliest timepoint is the reference level of timepoint_cat, and all
  ## treatment effects are reported against treatment_ref, so both are
  ## worth pinning down
  d <- prepare_clone_data(read_example_data(), treatment_ref = "arm_C")
  expect_identical(levels(d$timepoint_cat), c("0", "3", "6", "9"))
  expect_identical(levels(d$treatment)[1], "arm_C")
  expect_equal(d$perc, 100 * d$count / d$total)
  expect_error(prepare_clone_data(read_example_data(), treatment_ref = "arm_Z"),
               "not a value of")
})
