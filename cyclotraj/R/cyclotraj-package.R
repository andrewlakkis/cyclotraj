#' cyclotraj: Mixed-model testing of treatment effects on clone trajectories
#'
#' A reusable framework for testing treatment effects in nested longitudinal
#' count data: units measured repeatedly within subjects, where subjects are
#' assigned to treatment groups and each observation is a count out of a
#' sample-specific total. Counts are modeled with binomial generalized linear
#' mixed-effects models (logit link) in two complementary parameterizations of
#' time, read together:
#'
#' \itemize{
#'   \item \strong{Continuous-time (linear-slope)} model: estimates each
#'         treatment's overall rate of expansion or contraction.
#'   \item \strong{Difference-in-differences (categorical-time)} model: tests, at
#'         each timepoint, whether a treatment's change from baseline differs
#'         from the reference arm's.
#' }
#'
#' The core workflow is: \code{\link{prepare_clone_data}} ->
#' \code{\link{run_clone_analysis}} -> (\code{\link{results_table}},
#' \code{\link{plot_trajectory}}). See the package vignette for a worked example.
#'
#' @keywords internal
#' @importFrom ggplot2 .data
"_PACKAGE"

## quiet R CMD check notes for non-standard evaluation used in ggplot2/dplyr-free code
utils::globalVariables(c(
  "clone_id", "patient_id", "treatment", "timepoint_num", "timepoint_cat",
  "count", "total", "perc", "trajectory", ".fitted", ".lower", ".upper"
))
