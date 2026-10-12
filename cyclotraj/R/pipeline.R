#' Analyze a single trajectory with both GLMM parameterizations
#'
#' Fits the continuous-time and difference-in-differences models to one
#' trajectory's data and computes the associated inference (joint Wald tests,
#' slope contrasts, difference-in-differences contrasts, marginal means) and
#' diagnostics.
#'
#' @param data A single-trajectory \code{cyclotraj_data} frame.
#' @param drop_patient Passed to the fitters; drop the patient random intercept.
#' @param label Optional label stored with the result (e.g. the trajectory name).
#'
#' @return A list of class \code{cyclotraj_traj} with elements \code{label},
#'   \code{continuous}, \code{did} (the two \code{glmerMod} fits),
#'   \code{joint_wald}, \code{slopes}, \code{did_contrasts}, \code{emmeans},
#'   and \code{diagnostics}.
#' @export
analyze_trajectory <- function(data, drop_patient = FALSE, label = NULL) {
  m_cont <- fit_continuous(data, drop_patient = drop_patient)
  m_did  <- fit_did(data, drop_patient = drop_patient)

  diags <- rbind(
    cbind(model = "continuous", glmm_diagnostics(m_cont)),
    cbind(model = "did",        glmm_diagnostics(m_did))
  )

  out <- list(
    label = label,
    continuous = m_cont,
    did = m_did,
    joint_wald = list(
      continuous = joint_wald(m_cont),
      did = joint_wald(m_did)
    ),
    slopes = slope_contrasts(m_cont),
    did_contrasts = did_contrasts(m_did),
    emmeans = emm_means(m_did),
    diagnostics = diags
  )
  class(out) <- "cyclotraj_traj"
  out
}

#' Run the full trajectory analysis pipeline
#'
#' The single entry point. Optionally prepares the data, then fits and analyzes
#' every trajectory (or the whole data set as one group if no trajectory column
#' is present).
#'
#' Trajectories are independent analyses, so a group that cannot be fit (too few
#' arms, a single timepoint in an arm, a model that will not converge) is
#' reported as a warning and skipped rather than aborting the whole run; the
#' returned \code{trajectories} list contains only the groups that succeeded.
#'
#' @param data Either a raw data frame (in which case \code{prepare = TRUE}
#'   standardizes it) or an already-prepared \code{cyclotraj_data} frame.
#' @param prepare If \code{TRUE}, run \code{\link{prepare_clone_data}} on
#'   \code{data} first; extra arguments in \code{...} are forwarded to it.
#' @param drop_patient Drop the patient random intercept in all fits.
#' @param ... Passed to \code{\link{prepare_clone_data}} when \code{prepare = TRUE}.
#'
#' @return An object of class \code{cyclotraj_analysis}: a list with \code{data}
#'   (the prepared data) and \code{trajectories} (a named list of
#'   \code{cyclotraj_traj} results, one per successfully fitted group).
#'
#' @examples
#' \donttest{
#' res <- run_clone_analysis(read_example_data())
#' results_table(res)
#' }
#' @export
run_clone_analysis <- function(data, prepare = TRUE, drop_patient = FALSE, ...) {
  if (prepare) data <- prepare_clone_data(data, ...)
  if (!inherits(data, "cyclotraj_data")) {
    stop("`data` must be prepared; set prepare = TRUE or call prepare_clone_data().",
         call. = FALSE)
  }
  has_traj <- "trajectory" %in% names(data)
  traj_vals <- if (has_traj) levels(droplevels(data$trajectory)) else "all"

  results <- lapply(traj_vals, function(tv) {
    sub <- if (has_traj) data[data$trajectory == tv, , drop = FALSE] else data
    message("Fitting trajectory: ", tv,
            " (", length(unique(sub$clone_id)), " clones)")
    ## one group's failure (unestimable design, no convergence) must not take
    ## the other groups with it; warn, record nothing, and carry on
    tryCatch(
      analyze_trajectory(sub, drop_patient = drop_patient, label = tv),
      error = function(e) {
        warning("Trajectory ", sQuote(tv), " was skipped: ",
                conditionMessage(e), call. = FALSE)
        NULL
      })
  })
  names(results) <- traj_vals
  failed <- names(results)[vapply(results, is.null, logical(1))]
  results <- results[!names(results) %in% failed]
  if (!length(results)) {
    stop("No trajectory could be analyzed; see the warnings above.",
         call. = FALSE)
  }

  out <- list(data = data, trajectories = results)
  class(out) <- "cyclotraj_analysis"
  out
}

#' Combined tidy results table across trajectories
#'
#' Collates the \code{\link{joint_wald}} omnibus test (does each treatment arm's
#' trajectory differ from the reference arm's, in level and/or in shape?) for
#' both parameterizations into a single long table.
#'
#' @param x A \code{cyclotraj_analysis} object.
#' @return A data frame with columns \code{trajectory}, \code{model}
#'   (\code{"continuous"} or \code{"did"}), \code{treatment} (the non-reference
#'   arm being compared with the reference), \code{df} (degrees of freedom of the
#'   joint test), \code{chisq} and \code{p_value}. The p-values are the raw
#'   omnibus p-values: they are not adjusted across arms, parameterizations or
#'   trajectories, so apply a correction if the whole table is screened at once.
#' @export
results_table <- function(x) {
  stopifnot(inherits(x, "cyclotraj_analysis"))
  do.call(rbind, lapply(names(x$trajectories), function(tv) {
    tr <- x$trajectories[[tv]]
    rbind(
      cbind(trajectory = tv, model = "continuous", tr$joint_wald$continuous),
      cbind(trajectory = tv, model = "did",        tr$joint_wald$did)
    )
  }))
}

#' @export
print.cyclotraj_analysis <- function(x, ...) {
  cat("<cyclotraj_analysis>\n")
  cat("  trajectories:", length(x$trajectories),
      paste0("(", paste(names(x$trajectories), collapse = ", "), ")"), "\n")
  cat("  clones:", length(unique(x$data$clone_id)),
      " patients:", length(unique(x$data$patient_id)),
      " treatments:", nlevels(x$data$treatment), "\n")
  cat("  reference treatment:", levels(x$data$treatment)[1], "\n")
  cat("  use results_table(), plot_trajectory(), or inspect $trajectories\n")
  invisible(x)
}

#' @export
summary.cyclotraj_analysis <- function(object, ...) {
  rt <- results_table(object)
  cat("Joint Wald tests (treatment alters trajectory vs. reference):\n\n")
  rt$signif <- ifelse(is.na(rt$p_value), "",
                ifelse(rt$p_value < 0.001, "***",
                ifelse(rt$p_value < 0.01, "**",
                ifelse(rt$p_value < 0.05, "*", ""))))
  rt$p_value <- signif(rt$p_value, 3)
  print(rt, row.names = FALSE)
  invisible(rt)
}

#' @export
print.cyclotraj_traj <- function(x, ...) {
  cat("<cyclotraj_traj>", if (!is.null(x$label)) x$label else "", "\n")
  cat("  continuous joint Wald:\n"); print(x$joint_wald$continuous, row.names = FALSE)
  cat("  DiD contrasts with adjusted p < 0.05:",
      sum(x$did_contrasts$p.value < 0.05, na.rm = TRUE),
      "of", nrow(x$did_contrasts), "\n")
  invisible(x)
}
