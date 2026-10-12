#' Prepare longitudinal clone-count data for modeling
#'
#' Validates and standardizes a long-format data set (one row per
#' clone x patient x timepoint) into the canonical columns used by the rest of
#' the package. The function is deliberately generic: any nested longitudinal
#' count/proportion data can be supplied by mapping its column names through the
#' arguments below.
#'
#' @param data A data frame in long format.
#' @param clone,patient,treatment,timepoint,count,total Column names (as strings)
#'   for, respectively, the tracked unit (clone), the grouping subject (patient),
#'   the treatment/condition, the observation time (numeric weeks or similar),
#'   the observed count of the unit, and the sample total the count is out of.
#' @param trajectory Optional column name giving a pre-computed trajectory/cluster
#'   label (e.g. Cyclone's output). If \code{NULL}, all rows are treated as a
#'   single group.
#' @param treatment_ref Reference level for \code{treatment}. Defaults to the
#'   first level encountered. All treatment effects are reported relative to this
#'   arm.
#' @param trajectory_exclude Character vector of trajectory labels to drop
#'   (e.g. \code{"Other"} for non-selected clones). Ignored if \code{trajectory}
#'   is \code{NULL}.
#'
#' @details Validation is strict by design, because the failure modes it guards
#'   against are silent: a missing identifier or treatment value, or a row with
#'   \code{total == 0}, would otherwise be dropped (or turned into an undefined
#'   proportion) somewhere downstream rather than reported here. All such rows
#'   raise an error. Numeric coercion of \code{count} and \code{total} goes
#'   through character for non-numeric columns, so a factor-valued count yields
#'   its labels rather than its integer level codes; non-integer values are
#'   rounded with a warning.
#'
#' @return A data frame of class \code{cyclotraj_data} with standardized columns:
#'   \code{clone_id}, \code{patient_id}, \code{treatment} (factor, reference
#'   first), \code{timepoint_num} (numeric), \code{timepoint_cat} (factor,
#'   baseline first), \code{count}, \code{total}, \code{perc}, and (if supplied)
#'   \code{trajectory} (factor).
#'
#' @examples
#' d <- prepare_clone_data(read_example_data())
#' levels(d$treatment)
#'
#' @export
prepare_clone_data <- function(data,
                               clone = "clone_id",
                               patient = "patient_alias",
                               treatment = "treatment",
                               timepoint = "timepoint",
                               count = "n",
                               total = "n_tcells_tp",
                               trajectory = "k",
                               treatment_ref = NULL,
                               trajectory_exclude = "Other") {
  stopifnot(is.data.frame(data))
  req <- c(clone, patient, treatment, timepoint, count, total)
  missing_cols <- setdiff(req, names(data))
  if (length(missing_cols)) {
    stop("Missing required column(s): ", paste(missing_cols, collapse = ", "),
         call. = FALSE)
  }
  has_traj <- !is.null(trajectory) && trajectory %in% names(data)

  ## Coerce via character when the column is not already numeric: as.numeric()
  ## applied to a factor returns the integer level codes, which would silently
  ## substitute ranks for counts.
  as_num <- function(x) {
    if (is.numeric(x)) as.numeric(x) else suppressWarnings(as.numeric(as.character(x)))
  }

  out <- data.frame(
    clone_id = as.factor(data[[clone]]),
    patient_id = as.factor(data[[patient]]),
    treatment = as.character(data[[treatment]]),
    timepoint_num = as_num(data[[timepoint]]),
    count = as_num(data[[count]]),
    total = as_num(data[[total]]),
    stringsAsFactors = FALSE
  )
  if (has_traj) out$trajectory <- as.character(data[[trajectory]])

  ## ---- validation ----------------------------------------------------------
  if (any(is.na(out$timepoint_num))) {
    stop("`timepoint` must be numeric (or coercible to numeric).", call. = FALSE)
  }
  ## A missing identifier or treatment label is not benign: lme4 would drop
  ## those rows silently, changing the analysis set without saying so.
  id_na <- c(clone = anyNA(out$clone_id), patient = anyNA(out$patient_id),
             treatment = anyNA(out$treatment))
  if (any(id_na)) {
    stop("Missing values in: ", paste(names(id_na)[id_na], collapse = ", "),
         ". Rows with a missing identifier or treatment label would be dropped ",
         "silently when fitting; remove or impute them first.", call. = FALSE)
  }
  n_bad <- sum(is.na(out$count) | is.na(out$total))
  if (n_bad) {
    stop("`count` and `total` must be numeric and non-missing (", n_bad,
         " row(s) are missing or not coercible to a number).", call. = FALSE)
  }
  if (any(out$count < 0) || any(out$total < 0)) {
    stop("`count` and `total` must be non-negative.", call. = FALSE)
  }
  if (any(out$count > out$total)) {
    stop("`count` cannot exceed `total` in any row.", call. = FALSE)
  }
  ## total == 0 passes every check above but is a zero-trial binomial
  ## observation: it carries no information and makes `perc` undefined.
  n_zero <- sum(out$total == 0)
  if (n_zero) {
    stop(n_zero, " row(s) have `total` == 0. A sample with no observed units ",
         "contributes nothing to a binomial likelihood and gives an undefined ",
         "proportion; drop those rows before preparing the data.", call. = FALSE)
  }
  if (any(abs(out$count - round(out$count)) > 1e-8) ||
      any(abs(out$total - round(out$total)) > 1e-8)) {
    warning("Non-integer values in `count`/`total` were rounded; check that ",
            "these columns are counts and not proportions.", call. = FALSE)
  }
  out$count <- as.integer(round(out$count))
  out$total <- as.integer(round(out$total))
  if (anyNA(out$count) || anyNA(out$total)) {
    stop("`count`/`total` exceed the largest representable integer (",
         .Machine$integer.max, ").", call. = FALSE)
  }

  ## ---- trajectory filtering ------------------------------------------------
  if (has_traj) {
    ## NA labels survive `trajectory_exclude` but are silently skipped later by
    ## run_clone_analysis(), so drop them here and say so.
    na_traj <- is.na(out$trajectory)
    if (any(na_traj)) {
      warning(sum(na_traj), " row(s) have a missing trajectory label and were ",
              "excluded.", call. = FALSE)
      out <- out[!na_traj, , drop = FALSE]
    }
    if (length(trajectory_exclude)) {
      out <- out[!out$trajectory %in% trajectory_exclude, , drop = FALSE]
    }
    out$trajectory <- factor(out$trajectory)
  }
  if (!nrow(out)) {
    stop("No rows remain after trajectory filtering; check `trajectory_exclude`.",
         call. = FALSE)
  }

  ## ---- factor levels -------------------------------------------------------
  trt_levels <- unique(out$treatment)
  if (!is.null(treatment_ref)) {
    if (!treatment_ref %in% trt_levels) {
      stop("`treatment_ref` (\"", treatment_ref,
           "\") is not a value of `treatment`.", call. = FALSE)
    }
    trt_levels <- c(treatment_ref, setdiff(trt_levels, treatment_ref))
  }
  out$treatment <- factor(out$treatment, levels = trt_levels)
  out$clone_id <- droplevels(out$clone_id)
  out$patient_id <- droplevels(out$patient_id)

  tp_levels <- sort(unique(out$timepoint_num))
  out$timepoint_cat <- factor(out$timepoint_num, levels = tp_levels)
  out$perc <- 100 * out$count / out$total

  class(out) <- c("cyclotraj_data", "data.frame")
  out
}

#' Load the bundled example trajectory data
#'
#' Reads the synthetic example data distributed with the package: 500 units in
#' 6 trajectory groups, nested in 18 subjects across 4 unbalanced treatment arms,
#' each observed at weeks 0, 3, 6 and 9.
#'
#' The data are simulated from the binomial mixed model this package fits, with
#' parameters chosen to reproduce the structural features that motivate it: a
#' ~23-fold spread in sample totals, ~23\% zero counts, median frequency well
#' under 1\%, and unit-level variance far exceeding subject-level variance. The
#' generating script is \code{data-raw/make_example_data.R} in the source
#' repository. No real subject data is distributed with this package.
#'
#' @return A raw data frame with columns \code{clone_id}, \code{patient_alias},
#'   \code{treatment}, \code{timepoint}, \code{n}, \code{n_tcells_tp},
#'   \code{perc}, \code{k}. Pass it to
#'   \code{\link{prepare_clone_data}} before modeling. Treatment arms are
#'   \code{arm_A} (the natural reference) through \code{arm_D}; trajectory
#'   labels are \code{"Traj 1"} to \code{"Traj 6"}.
#' @export
read_example_data <- function() {
  path <- system.file("extdata", "example_trajectories.csv",
                      package = "cyclotraj")
  if (!nzchar(path)) {
    stop("Example data not found; is the package installed correctly?",
         call. = FALSE)
  }
  utils::read.csv(path, stringsAsFactors = FALSE)
}
