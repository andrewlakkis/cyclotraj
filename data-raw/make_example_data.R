#!/usr/bin/env Rscript
# =============================================================================
# Generate the synthetic example dataset bundled with cyclotraj.
#
# The data are SIMULATED from the same binomial GLMM the package fits, so the
# example exercises every exported function on data whose generating process is
# known. No real subject data is used anywhere in this package.
#
# Parameters are chosen to reproduce the structural features that make this
# inference problem hard, because those are what the modeling choices are for:
#
#   * sample totals spanning ~23-fold (drawn log-uniform over a 27-fold design
#     range), so per-observation precision varies
#   * ~23% of counts exactly zero, so transforms would need ad hoc offsets
#   * small median count (2) at frequencies well under 1% (median 0.29%)
#   * unbalanced arms (4/4/7/3 subjects) with few subjects per arm
#   * unit-level variance >> subject-level variance, so some fits are singular
#     and glmm_diagnostics() has something real to report
#
# Six trajectory groups are simulated with deliberately different shapes so the
# continuous-time and difference-in-differences models disagree in interpretable
# ways (see the vignette):
#
#   Traj 1  contracting, arm_D steeper slope
#   Traj 2  contracting, no arm effects at all (a true null group)
#   Traj 3  expanding, arm_D gets a transient one-timepoint spike at week 3
#   Traj 4  expanding, arm_B gets a delayed step up at weeks 6 and 9
#   Traj 5  expanding, arm_D steeper slope, sustained
#   Traj 6  expanding, arm_C slope dampened
#
# These are the generating mechanisms, not predictions about which test wins:
# a non-monotone mechanism (Traj 3) is where the two encodings diverge most,
# since a single slope averages a one-timepoint spike away.
#
# Usage:  Rscript data-raw/make_example_data.R   (run from the repository root)
# Output: cyclotraj/inst/extdata/example_trajectories.csv
#
# Re-running this script must reproduce the committed CSV; the R-CMD-check
# workflow regenerates it and fails if it differs. The random-number
# generators are named explicitly rather than left to the defaults, so the
# output does not depend on which R release happens to run it (these are the
# defaults for R >= 3.6.0).
# =============================================================================

set.seed(20270601, kind = "Mersenne-Twister",
         normal.kind = "Inversion", sample.kind = "Rejection")

ARMS   <- c(arm_A = 4L, arm_B = 4L, arm_C = 7L, arm_D = 3L)   # ref arm first
WEEKS  <- c(0, 3, 6, 9)
MIN_UNITS_PER_ARM <- 6L   # floor per arm per group; smallest group has 29
NUNITS <- c("Traj 1" = 111L, "Traj 2" = 29L, "Traj 3" = 104L,
            "Traj 4" = 104L, "Traj 5" = 64L, "Traj 6" = 88L)

SD_UNIT    <- 0.95   # between-unit baseline spread (large), truncated below
SD_SUBJECT <- 0.06   # between-subject spread (near zero by design)

## ---- subjects -------------------------------------------------------------
subjects <- data.frame(
  patient_alias = sprintf("subj_%02d", seq_len(sum(ARMS))),
  treatment     = rep(names(ARMS), times = ARMS),
  stringsAsFactors = FALSE
)
subjects$u_subject <- stats::rnorm(nrow(subjects), 0, SD_SUBJECT)

## ---- sequencing depth: one total per subject x timepoint ------------------
## log-uniform over the observed 27-fold range, so depth is unrelated to arm.
depth <- expand.grid(patient_alias = subjects$patient_alias, timepoint = WEEKS,
                     stringsAsFactors = FALSE)
depth$n_tcells_tp <- round(exp(stats::runif(nrow(depth), log(140), log(3800))))

## ---- fixed effects, per trajectory ----------------------------------------
## b0 = baseline logit; b1 = reference-arm slope per week;
## slope = extra slope per week for an arm; bump = extra logit at given weeks.
spec <- list(
  "Traj 1" = list(b0 = -5.60, b1 = -0.075,
                  slope = c(arm_B = -0.030, arm_C = -0.015, arm_D = -0.145)),
  "Traj 2" = list(b0 = -6.25, b1 = -0.060, slope = c()),
  "Traj 3" = list(b0 = -6.10, b1 =  0.055,
                  slope = c(arm_C = 0.015),
                  bump  = list(arm_D = c(`3` = 0.95))),
  "Traj 4" = list(b0 = -6.15, b1 =  0.050,
                  bump  = list(arm_B = c(`6` = 0.80, `9` = 0.95))),
  "Traj 5" = list(b0 = -5.90, b1 =  0.070,
                  slope = c(arm_B = 0.020, arm_D = 0.135)),
  "Traj 6" = list(b0 = -6.00, b1 =  0.085,
                  slope = c(arm_C = -0.090, arm_D = 0.025))
)

## ---- simulate --------------------------------------------------------------
rows <- list()
uid  <- 0L

for (k in names(NUNITS)) {
  s <- spec[[k]]
  ## Assign units to subjects, so each unit sits in exactly one subject (and
  ## therefore inherits exactly one arm), as in real clone data. Units are
  ## allocated per arm -- a floor of MIN_UNITS_PER_ARM, then the remainder split
  ## in proportion to arm size -- rather than by sampling subjects uniformly.
  ## Uniform sampling left the smallest arm with 2 units in the smallest group,
  ## whose counts reached zero at the last two timepoints; the categorical-time
  ## model then has no events to estimate those interactions from and separates.
  base    <- stats::setNames(rep(MIN_UNITS_PER_ARM, length(ARMS)), names(ARMS))
  rest    <- NUNITS[[k]] - sum(base)
  extra   <- diff(c(0, round(cumsum(ARMS / sum(ARMS)) * rest)))
  per_arm <- base + extra
  owner <- unlist(lapply(names(ARMS), function(a)
    sample(subjects$patient_alias[subjects$treatment == a], per_arm[[a]],
           replace = TRUE)))

  for (i in seq_len(NUNITS[[k]])) {
    uid <- uid + 1L
    pid <- owner[i]
    arm <- subjects$treatment[match(pid, subjects$patient_alias)]
    ## truncated at +/-2.2 SD: an untruncated normal occasionally yields a unit
    ## an order of magnitude more abundant than anything in real repertoire data
    u_i <- max(-2.2 * SD_UNIT, min(2.2 * SD_UNIT, stats::rnorm(1, 0, SD_UNIT)))
    v_m <- subjects$u_subject[match(pid, subjects$patient_alias)]

    extra_slope <- if (length(s$slope) && arm %in% names(s$slope)) s$slope[[arm]] else 0
    bumps <- if (!is.null(s$bump) && arm %in% names(s$bump)) s$bump[[arm]] else NULL

    for (w in WEEKS) {
      bump <- if (!is.null(bumps) && as.character(w) %in% names(bumps))
                bumps[[as.character(w)]] else 0
      eta <- s$b0 + u_i + v_m + (s$b1 + extra_slope) * w + bump
      p   <- stats::plogis(eta)
      N   <- depth$n_tcells_tp[depth$patient_alias == pid & depth$timepoint == w]
      rows[[length(rows) + 1L]] <- data.frame(
        clone_id      = sprintf("unit_%04d", uid),
        patient_alias = pid,
        treatment     = arm,
        timepoint     = w,
        n             = stats::rbinom(1, N, p),
        n_tcells_tp   = N,
        k             = k,
        stringsAsFactors = FALSE
      )
    }
  }
}

out <- do.call(rbind, rows)
out$perc <- 100 * out$n / out$n_tcells_tp
out <- out[order(out$k, out$clone_id, out$timepoint), ]

## paths are relative, so refuse to write a stray CSV somewhere else
if (!file.exists(file.path("cyclotraj", "DESCRIPTION"))) {
  stop("Run this script from the repository root (the directory holding ",
       "cyclotraj/DESCRIPTION).", call. = FALSE)
}
dest <- file.path("cyclotraj", "inst", "extdata", "example_trajectories.csv")
dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(out, dest, row.names = FALSE)

## ---- report the structural features, so drift is visible -------------------
cat(sprintf("rows                %d\n", nrow(out)))
cat(sprintf("units               %d across %d groups\n",
            length(unique(out$clone_id)), length(unique(out$k))))
cat(sprintf("subjects            %d in %d arms (%s)\n",
            length(unique(out$patient_alias)), length(ARMS),
            paste(ARMS, collapse = "/")))
cat(sprintf("total N             %d to %d (%.1f-fold)\n",
            min(out$n_tcells_tp), max(out$n_tcells_tp),
            max(out$n_tcells_tp) / min(out$n_tcells_tp)))
cat(sprintf("count n             median %g, max %d\n",
            stats::median(out$n), max(out$n)))
cat(sprintf("zero counts         %.1f%%\n", 100 * mean(out$n == 0)))
cat(sprintf("frequency           median %.3f%%, max %.2f%%\n",
            stats::median(out$perc), max(out$perc)))
cat(sprintf("wrote               %s\n", dest))
