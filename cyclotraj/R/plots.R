#' Population-level predictions with confidence bands
#'
#' Computes fitted proportions and 95% confidence limits for a fitted GLMM.
#' Intervals are built on the linear predictor, where the standard error is
#' \code{sqrt(t(x) \%*\% V \%*\% x)} for the fixed-effect covariance \code{V},
#' and the limits
#' are then back-transformed through the inverse logit; this keeps them inside
#' (0, 1) and avoids a delta-method approximation on the response scale. The
#' random effects are held at zero, so the curve describes a typical unit in a
#' typical subject and the band reflects fixed-effect uncertainty only. This
#' backs the confidence ribbons drawn by \code{\link{plot_trajectory}}.
#'
#' @param model A \code{glmerMod} from \code{\link{fit_continuous}} or
#'   \code{\link{fit_did}}.
#' @param data The \code{cyclotraj_data} the model was fit on (used to build the
#'   prediction grid and align factor levels).
#' @param level Confidence level (default 0.95).
#' @param n_grid Number of time points in the grid for the continuous model.
#'
#' @return A data frame with \code{treatment}, \code{timepoint_num},
#'   \code{fitted}, \code{lower}, \code{upper} (all on the percentage scale,
#'   0-100).
#' @export
predict_with_ci <- function(model, data, level = 0.95, n_grid = 100) {
  time_term <- attr(model, "cyclotraj_time")
  if (is.null(time_term)) time_term <- "timepoint_num"
  ## Subsetting to one trajectory leaves the factor levels of the full data set
  ## behind, while the fitter drops the unused ones. Without this droplevels()
  ## the grid would contain arms (or timepoints) the model never saw; their
  ## design-matrix columns are not in `beta` and would be silently discarded
  ## below, giving those rows the reference arm's / baseline's prediction.
  data <- droplevels(as.data.frame(data))
  trt_levels <- levels(data$treatment)
  tp_levels <- levels(data$timepoint_cat)
  tp_num <- as.numeric(as.character(tp_levels))

  if (identical(time_term, "timepoint_num")) {
    grid <- expand.grid(
      treatment = factor(trt_levels, levels = trt_levels),
      timepoint_num = seq(min(tp_num), max(tp_num), length.out = n_grid)
    )
    fixed <- stats::formula(~ treatment * timepoint_num)
  } else {
    grid <- expand.grid(
      treatment = factor(trt_levels, levels = trt_levels),
      timepoint_cat = factor(tp_levels, levels = tp_levels)
    )
    grid$timepoint_num <- as.numeric(as.character(grid$timepoint_cat))
    fixed <- stats::formula(~ treatment * timepoint_cat)
  }

  X <- stats::model.matrix(fixed, grid)
  beta <- lme4::fixef(model)
  ## the grid must reproduce every column the model estimated, otherwise the
  ## prediction would silently correspond to a different model than was fit
  absent <- setdiff(names(beta), colnames(X))
  if (length(absent)) {
    stop("`data` does not reproduce the model's design: no grid column for ",
         paste(absent, collapse = ", "),
         ". Pass the same data the model was fit on.", call. = FALSE)
  }
  X <- X[, names(beta), drop = FALSE]
  V <- as.matrix(stats::vcov(model))

  eta <- as.vector(X %*% beta)
  se <- sqrt(pmax(0, rowSums((X %*% V) * X)))
  z <- stats::qnorm(1 - (1 - level) / 2)

  inv <- stats::plogis
  data.frame(
    treatment = grid$treatment,
    timepoint_num = grid$timepoint_num,
    fitted = 100 * inv(eta),
    lower = 100 * inv(eta - z * se),
    upper = 100 * inv(eta + z * se)
  )
}

#' Plot a single trajectory: observed clones with a fitted GLMM overlay
#'
#' Draws individual clones (grey spaghetti) and the population-level GLMM fit
#' with a 95% confidence band, faceted by treatment. Works for either model; the
#' continuous fit is drawn as a smooth line, the difference-in-differences fit as
#' points connected across the discrete timepoints.
#'
#' @param data A single-trajectory \code{cyclotraj_data} frame.
#' @param model A fitted GLMM (\code{\link{fit_continuous}} or
#'   \code{\link{fit_did}}).
#' @param title,subtitle Plot annotations.
#' @param fit_color Colour of the fitted line/band.
#' @return A \code{ggplot} object.
#' @export
plot_trajectory <- function(data, model,
                            title = NULL, subtitle = NULL,
                            fit_color = "#2166AC") {
  time_term <- attr(model, "cyclotraj_time")
  is_cont <- identical(time_term, "timepoint_num")
  ## drop arms with no rows in this subset, so the facets match the arms the
  ## model was fit on instead of adding an empty NA-labelled panel
  data <- droplevels(as.data.frame(data))
  pred <- predict_with_ci(model, data)

  counts <- stats::aggregate(clone_id ~ treatment, data = data,
                             FUN = function(x) length(unique(x)))
  names(counts)[2] <- "n_clones"
  lab <- stats::setNames(
    sprintf("%s (n=%d)", counts$treatment, counts$n_clones),
    as.character(counts$treatment))
  data$treatment_label <- factor(lab[as.character(data$treatment)],
                                 levels = lab[levels(data$treatment)])
  pred$treatment_label <- factor(lab[as.character(pred$treatment)],
                                 levels = lab[levels(data$treatment)])

  p <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = data,
      ggplot2::aes(.data$timepoint_num, .data$perc, group = .data$clone_id),
      alpha = 0.25, colour = "grey55") +
    ggplot2::geom_ribbon(
      data = pred,
      ggplot2::aes(.data$timepoint_num, ymin = .data$lower, ymax = .data$upper),
      fill = fit_color, alpha = 0.20) +
    ggplot2::geom_line(
      data = pred,
      ggplot2::aes(.data$timepoint_num, .data$fitted),
      colour = fit_color, linewidth = 1.1) +
    ggplot2::facet_wrap(~ treatment_label) +
    ggplot2::labs(
      title = title, subtitle = subtitle,
      x = "Time (weeks)", y = "Clone frequency (%)") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold"))

  if (!is_cont) {
    p <- p + ggplot2::geom_point(
      data = pred,
      ggplot2::aes(.data$timepoint_num, .data$fitted),
      colour = fit_color, size = 2)
  }
  p
}
