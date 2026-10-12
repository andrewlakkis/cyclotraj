#' Fit the continuous-time (linear-slope) binomial GLMM
#'
#' Treats time as a numeric covariate, so the treatment-by-time interaction
#' measures each treatment's overall rate (slope) of clonal expansion or
#' contraction relative to the reference arm.
#'
#' Model (logit link):
#' \code{cbind(count, total - count) ~ treatment * timepoint_num +
#' (1 | clone_id) + (1 | patient_id)}.
#'
#' @param data A \code{cyclotraj_data} frame (see \code{\link{prepare_clone_data}}),
#'   typically already subset to a single trajectory.
#' @param drop_patient If \code{TRUE}, omit the patient random intercept (useful
#'   when it is singular, e.g. very few patients per arm).
#' @param optimizer,maxfun Passed to \code{lme4::glmerControl}.
#'
#' @return A \code{glmerMod} object.
#' @seealso \code{\link{fit_did}}
#' @export
fit_continuous <- function(data, drop_patient = FALSE,
                           optimizer = "bobyqa", maxfun = 2e5) {
  .fit_glmm(data, time_term = "timepoint_num",
            drop_patient = drop_patient,
            optimizer = optimizer, maxfun = maxfun)
}

#' Fit the difference-in-differences (categorical-time) binomial GLMM
#'
#' Treats time as a factor with the baseline timepoint as reference, so the
#' treatment-by-timepoint interaction is the difference-in-differences estimator:
#' how much a treatment's change from baseline to each timepoint exceeds the
#' reference arm's change over the same interval. Makes no linearity assumption.
#'
#' Model (logit link):
#' \code{cbind(count, total - count) ~ treatment * timepoint_cat +
#' (1 | clone_id) + (1 | patient_id)}.
#'
#' @inheritParams fit_continuous
#' @return A \code{glmerMod} object.
#' @seealso \code{\link{fit_continuous}}
#' @export
fit_did <- function(data, drop_patient = FALSE,
                    optimizer = "bobyqa", maxfun = 2e5) {
  .fit_glmm(data, time_term = "timepoint_cat",
            drop_patient = drop_patient,
            optimizer = optimizer, maxfun = maxfun)
}

## internal fitter shared by both parameterizations
.fit_glmm <- function(data, time_term, drop_patient, optimizer, maxfun) {
  if (!inherits(data, "cyclotraj_data")) {
    warning("`data` is not a cyclotraj_data object; did you run prepare_clone_data()?")
  }
  needed <- c("count", "total", "treatment", "clone_id", time_term)
  if (!drop_patient) needed <- c(needed, "patient_id")
  miss <- setdiff(needed, names(data))
  if (length(miss)) {
    stop("`data` is missing required column(s): ", paste(miss, collapse = ", "),
         ". Run prepare_clone_data() first.", call. = FALSE)
  }

  trt_before <- levels(data$treatment)
  data <- droplevels(as.data.frame(data))
  if (!nrow(data)) stop("`data` has no rows to fit.", call. = FALSE)

  ## Subsetting to one trajectory can empty a whole arm. droplevels() then
  ## renumbers the factor, and if the dropped arm was the reference every
  ## coefficient is suddenly relative to a different baseline.
  trt_after <- levels(data$treatment)
  if (length(trt_before) && !identical(trt_before, trt_after)) {
    msg <- paste0("Treatment level(s) ",
                  paste(sQuote(setdiff(trt_before, trt_after)), collapse = ", "),
                  " have no rows in this subset and were dropped")
    if (!trt_before[1] %in% trt_after) {
      msg <- paste0(msg, "; the reference arm is now ", sQuote(trt_after[1]),
                    ", so all treatment effects are relative to it")
    }
    warning(msg, ".", call. = FALSE)
  }

  ## Guard the identifiability conditions for this design up front: without
  ## these, lme4 either aborts with a message about contrasts or silently drops
  ## rank-deficient interaction columns, which would quietly shrink the
  ## hypothesis tested by joint_wald().
  n_trt <- length(unique(data$treatment))
  if (n_trt < 2L) {
    stop("`treatment` has only one level in this subset; at least two groups ",
         "are needed to estimate a treatment effect.", call. = FALSE)
  }
  n_time <- length(unique(data[[time_term]]))
  if (n_time < 2L) {
    stop("`", time_term, "` has only one distinct value in this subset; at ",
         "least two timepoints are needed to estimate a time effect.",
         call. = FALSE)
  }
  tp_per_arm <- tapply(data[[time_term]], factor(data$treatment),
                       function(z) length(unique(z)))
  thin <- names(tp_per_arm)[is.na(tp_per_arm) | tp_per_arm < 2L]
  if (length(thin)) {
    stop("Treatment arm(s) ", paste(sQuote(thin), collapse = ", "),
         " are observed at a single timepoint, so their treatment-by-time ",
         "terms are not estimable.", call. = FALSE)
  }
  if (length(unique(data$clone_id)) < 2L) {
    stop("`clone_id` has only one level in this subset; a unit-level random ",
         "intercept is not identifiable.", call. = FALSE)
  }
  if (!drop_patient && length(unique(data$patient_id)) < 2L) {
    stop("`patient_id` has only one level in this subset; the subject random ",
         "intercept cannot be estimated. Re-fit with drop_patient = TRUE.",
         call. = FALSE)
  }

  re <- if (drop_patient) "(1 | clone_id)" else "(1 | clone_id) + (1 | patient_id)"
  form <- stats::as.formula(sprintf(
    "cbind(count, total - count) ~ treatment * %s + %s", time_term, re))
  ctrl <- lme4::glmerControl(optimizer = optimizer,
                             optCtrl = list(maxfun = maxfun))
  fit <- lme4::glmer(form, family = stats::binomial(link = "logit"),
                     data = data, control = ctrl)
  attr(fit, "cyclotraj_time") <- time_term
  fit
}

#' Convergence and singularity diagnostics for a fitted GLMM
#'
#' @param model A \code{glmerMod} from \code{\link{fit_continuous}} or
#'   \code{\link{fit_did}}.
#' @return A one-row data frame with the number of observations, whether the fit
#'   is singular, any convergence-warning messages, and the estimated
#'   random-effect standard deviations for clone and patient.
#' @export
glmm_diagnostics <- function(model) {
  vc <- as.data.frame(lme4::VarCorr(model))
  sd_of <- function(grp) {
    v <- vc$sdcor[vc$grp == grp & is.na(vc$var2)]
    if (length(v)) v[1] else NA_real_
  }
  warns <- model@optinfo$conv$lme4$messages
  data.frame(
    n_obs = stats::nobs(model),
    singular = lme4::isSingular(model),
    converged = is.null(warns),
    message = if (is.null(warns)) "" else paste(warns, collapse = "; "),
    sd_clone = sd_of("clone_id"),
    sd_patient = sd_of("patient_id"),
    stringsAsFactors = FALSE
  )
}
