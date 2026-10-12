#' Joint Wald test that a treatment alters the trajectory at all
#'
#' For each non-reference treatment, jointly tests all of its coefficients (its
#' main effect plus every treatment-by-time interaction term) against zero with a
#' multi-degree-of-freedom Wald chi-square test. This asks "does this treatment
#' change the trajectory relative to the reference arm in any way?", which is more
#' powerful and less multiplicity-prone than inspecting individual coefficients.
#' Works for either parameterization (the interaction terms differ, but they are
#' detected automatically).
#'
#' Including the main effect is deliberate: the omnibus question is whether an
#' arm's whole trajectory differs from the reference arm's, and the main effect
#' carries the difference in level (the vertical offset) while the interactions
#' carry the difference in shape. A consequence is that a group difference
#' already present at the reference timepoint contributes to the test, so a
#' significant result means "differs in level and/or in shape". Note also that
#' the main effect is the group difference at \code{timepoint_num == 0}, which is
#' interpretable only because time is coded with baseline at zero. For tests that
#' condition on baseline and so isolate the shape, use \code{\link{slope_contrasts}}
#' (continuous model) or \code{\link{did_contrasts}} (categorical model).
#'
#' @param model A \code{glmerMod} from \code{\link{fit_continuous}} or
#'   \code{\link{fit_did}}.
#' @return A data frame with one row per non-reference treatment: \code{treatment},
#'   \code{df}, \code{chisq}, \code{p_value}. P-values are per treatment and are
#'   not adjusted for the number of treatments compared.
#' @export
joint_wald <- function(model) {
  beta <- lme4::fixef(model)
  cf <- names(beta)
  ## split each coefficient name into its interaction components so we can match
  ## a treatment term exactly. Component-wise equality (rather than a regex or
  ## prefix match on the whole name) keeps levels that are prefixes of one
  ## another apart, e.g. "drug-A" vs "drug-A (naive)", and tolerates spaces,
  ## hyphens and parentheses in level names. A ":" inside a level name is the
  ## one case it cannot resolve; that is caught and reported below.
  parts <- strsplit(cf, ":", fixed = TRUE)
  trt <- model@frame$treatment
  if (is.null(trt)) {
    stop("`model` has no `treatment` variable; joint_wald() expects a fit from ",
         "fit_continuous() or fit_did().", call. = FALSE)
  }
  ## a character/logical treatment is converted to a factor with sorted levels
  ## when the model matrix is built, so mirror that to recover the arm names
  trt_levels <- if (is.factor(trt)) levels(trt) else sort(unique(as.character(trt)))
  others <- trt_levels[-1]
  if (!length(others)) {
    return(data.frame(treatment = character(0), df = integer(0),
                      chisq = numeric(0), p_value = numeric(0),
                      stringsAsFactors = FALSE))
  }

  do.call(rbind, lapply(others, function(tt) {
    tag <- paste0("treatment", tt)
    idx <- which(vapply(parts, function(p) tag %in% p, logical(1)))
    if (!length(idx)) {
      warning("No fixed-effect coefficients matched treatment level ",
              sQuote(tt), "; returning NA. If the level name contains a \":\" ",
              "it cannot be matched against coefficient names: rename the ",
              "level and re-fit.", call. = FALSE)
      return(data.frame(treatment = tt, df = NA_integer_,
                        chisq = NA_real_, p_value = NA_real_,
                        stringsAsFactors = FALSE))
    }
    ## build a contrast matrix selecting exactly this treatment's coefficients
    L <- matrix(0, nrow = length(idx), ncol = length(beta),
                dimnames = list(NULL, cf))
    L[cbind(seq_along(idx), idx)] <- 1
    lh <- car::linearHypothesis(model, hypothesis.matrix = L,
                                rhs = rep(0, nrow(L)), test = "Chisq")
    data.frame(
      treatment = tt,
      df = lh$Df[2],
      chisq = lh$Chisq[2],
      p_value = lh$`Pr(>Chisq)`[2],
      stringsAsFactors = FALSE
    )
  }))
}

#' Pairwise comparison of treatment slopes (continuous-time model)
#'
#' Uses \code{emmeans::emtrends} to compare the estimated time slopes (rate of
#' change per unit time on the logit scale) between every pair of treatments,
#' with Tukey adjustment for multiplicity.
#'
#' @param model A \code{glmerMod} from \code{\link{fit_continuous}}.
#' @return A tidy data frame of pairwise slope contrasts. Estimates are
#'   differences of log-odds slopes, so \code{exp(estimate)} is the ratio of the
#'   two arms' per-unit-time odds ratios. The Tukey adjustment covers the
#'   pairwise comparisons within this one model; it does not account for fitting
#'   several trajectories.
#' @export
slope_contrasts <- function(model) {
  emt <- emmeans::emtrends(model, ~ treatment, var = "timepoint_num")
  as.data.frame(emmeans::contrast(emt, method = "pairwise", adjust = "tukey"))
}

#' Difference-in-differences contrasts (categorical-time model)
#'
#' Computes the contrast-of-contrasts: for each pair of treatments and each pair
#' of timepoints, how much the change between those timepoints differs between
#' the two treatments. Tukey-adjusted within each pair of timepoints, i.e. the
#' adjustment covers the treatment comparisons made at a given pair of
#' timepoints, not the full set of timepoint pairs.
#'
#' @param model A \code{glmerMod} from \code{\link{fit_did}}.
#' @param alpha Optional significance threshold applied to the adjusted
#'   p-values; if not \code{NULL}, only rows with \code{p.value < alpha} are
#'   returned.
#' @return A tidy data frame of difference-in-differences estimates with columns
#'   \code{treatments} (the pair of arms), \code{timepoints} (the pair of
#'   timepoints) and the estimate on the log odds ratio scale, so
#'   \code{exp(estimate)} is a ratio of odds ratios.
#' @export
did_contrasts <- function(model, alpha = NULL) {
  emm <- emmeans::emmeans(model, ~ timepoint_cat | treatment)
  ## first contrast: change between each pair of timepoints, within an arm
  time_diffs <- emmeans::contrast(emm, "pairwise", by = "treatment")
  ## second contrast: compare those changes between arms. emmeans keeps the
  ## inner contrast labels in `contrast` and names the new (treatment) ones
  ## `contrast1`, hence the renaming below.
  did <- emmeans::contrast(time_diffs, "pairwise", by = "contrast")
  ## drop to a plain data frame before renaming: as.data.frame() on an emmeans
  ## object keeps class "summary_emm" and a "by.vars" attribute pointing at the
  ## old "contrast" column, so renaming it leaves print() looking for a column
  ## that no longer exists.
  out <- as.data.frame(did)
  class(out) <- "data.frame"
  attr(out, "by.vars") <- NULL
  names(out)[names(out) == "contrast1"] <- "treatments"
  names(out)[names(out) == "contrast"] <- "timepoints"
  if (!is.null(alpha)) out <- out[out$p.value < alpha, , drop = FALSE]
  out
}

#' Estimated marginal means by treatment and timepoint (categorical-time model)
#'
#' Means are formed on the logit scale and back-transformed, with the random
#' effects held at zero, so they describe a typical unit in a typical subject
#' rather than the average proportion across units.
#'
#' @param model A \code{glmerMod} from \code{\link{fit_did}}.
#' @return A data frame of estimated marginal proportions (response scale) with
#'   confidence limits, by treatment and timepoint.
#' @export
emm_means <- function(model) {
  emm <- emmeans::emmeans(model, ~ treatment | timepoint_cat, type = "response")
  as.data.frame(emm)
}
