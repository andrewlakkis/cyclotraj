# =============================================================================
# Helper functions for the cyclotraj Shiny explorer.
# Plotting + table utilities. Namespaces are fully qualified so these load
# regardless of package attach order (shiny autoloads this R/ directory).
# =============================================================================

# ---- data upload ------------------------------------------------------------
# Guess the field separator of a delimited text file from its header line.
sniff_delim <- function(path, default = ",") {
  line <- tryCatch(readLines(path, n = 1L, warn = FALSE),
                   error = function(e) character(0))
  if (!length(line) || !nzchar(line)) return(default)
  cands <- c(",", "\t", ";", "|")
  hits <- vapply(cands, function(d) {
    m <- gregexpr(d, line[1], fixed = TRUE)[[1]]
    if (m[1] == -1L) 0L else length(m)
  }, integer(1))
  if (max(hits) == 0L) return(default)
  cands[[which.max(hits)]]
}

# Read an uploaded data file, dispatching on the ORIGINAL filename's extension
# (fileInput saves the upload under a temporary name that has none). Supports
# CSV, TSV/TXT (delimiter sniffed from the header), .rds holding a data.frame,
# and Excel if the optional readxl package happens to be installed. Column names
# are kept verbatim because the sidebar offers them to the user for mapping.
# Every stop() here surfaces as a readable message in the UI (see raw()).
read_upload <- function(path, filename = path) {
  if (!file.exists(path) || file.size(path) == 0) {
    stop("The uploaded file is empty.")
  }
  ext <- tolower(tools::file_ext(filename))
  df <- switch(ext,
    csv = , tsv = , txt = tryCatch(
      utils::read.table(path, header = TRUE, sep = sniff_delim(path),
                        quote = "\"", comment.char = "", check.names = FALSE,
                        stringsAsFactors = FALSE),
      error = function(e) stop(conditionMessage(e),
        ". Check that the first row holds the column names and that every row ",
        "has the same number of fields.", call. = FALSE)),
    xlsx = , xls = {
      if (!requireNamespace("readxl", quietly = TRUE)) {
        stop("Reading Excel files needs the optional 'readxl' package; ",
             "install it or save the sheet as CSV.")
      }
      as.data.frame(readxl::read_excel(path))
    },
    rds = {
      obj <- readRDS(path)
      if (!is.data.frame(obj)) stop("The .rds file must contain a data.frame.")
      as.data.frame(obj)
    },
    stop(sprintf(
      "Unsupported file type '.%s'. Use csv, tsv, txt, xlsx, xls, or rds.", ext))
  )
  if (!ncol(df)) stop("No columns were read from the file.")
  if (!nrow(df)) stop("The file has a header row but no data rows.")
  ## drop a UTF-8 byte-order mark so the first column name stays matchable
  names(df)[1] <- sub("^﻿", "", names(df)[1])
  df
}

# ---- column helpers ---------------------------------------------------------
# Candidate reference arms in a raw (unprepared) data frame, in the order the
# sidebar offers them. Shared by the selector and by the server's resolved
# reference level so the two always agree on the default.
treatment_levels <- function(df, col) {
  if (is.null(col) || !length(col) || !col %in% names(df)) return(character(0))
  sort(unique(as.character(df[[col]])))
}

# First of `candidates` that is actually a column of `df`. emmeans names its
# estimate and confidence-limit columns differently depending on the response
# scale and how degrees of freedom are computed, so look the names up instead of
# assuming one spelling.
first_col <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (!length(hit)) {
    stop("Expected one of the columns: ", paste(candidates, collapse = ", "),
         call. = FALSE)
  }
  hit[[1]]
}

# ---- palette / theme --------------------------------------------------------
.cyclo <- list(
  teal   = "#16736B",
  blue   = "#2166AC",
  red    = "#B2182B",
  amber  = "#D9A441",
  grey   = "#7F7F7F",
  ink    = "#1c2b2b"
)

# Colours for treatment arms, assigned by POSITION rather than by name so that
# any data set gets the same reading: the reference arm (always the first level)
# is grey, the arms compared against it take a qualitative scale.
arm_palette <- function(levels) {
  n <- length(levels)
  if (!n) return(stats::setNames(character(0), character(0)))
  base <- c("#E08214", "#4DAF4A", "#6A51A3", "#2166AC", "#B2182B")
  rest <- if (n - 1L <= length(base)) {
    base[seq_len(max(n - 1L, 0L))]
  } else {
    grDevices::hcl.colors(n - 1L, palette = "Dark 3")
  }
  stats::setNames(c(.cyclo$grey, rest), levels)
}

signif_stars <- function(p) {
  ifelse(is.na(p), "",
  ifelse(p < 0.001, "***",
  ifelse(p < 0.01,  "**",
  ifelse(p < 0.05,  "*", ""))))
}

fmt_p <- function(p) {
  ifelse(is.na(p), "—",
  ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3)))
}

theme_cyclo <- function(base = 13) {
  ggplot2::theme_minimal(base_size = base) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title    = ggplot2::element_text(face = "bold"),
      strip.text    = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )
}

# ---- cross-trajectory joint-Wald heatmap -----------------------------------
# long table of the omnibus test (does each arm alter each trajectory?)
jointwald_long <- function(res, model = c("continuous", "did")) {
  model <- match.arg(model)
  rt <- results_table(res)
  ## a fit can fail for every group (or there may be no non-reference arm left),
  ## in which case results_table() has nothing to collate
  if (is.null(rt) || !nrow(rt)) {
    return(data.frame(trajectory = character(0), model = character(0),
                      treatment = character(0), p_value = numeric(0),
                      neglogp = numeric(0), stars = character(0),
                      stringsAsFactors = FALSE))
  }
  rt <- rt[rt$model == model, , drop = FALSE]
  rt$neglogp <- -log10(pmax(rt$p_value, 1e-16))
  rt$stars   <- signif_stars(rt$p_value)
  rt$trajectory <- factor(rt$trajectory,
                          levels = rev(sort(unique(rt$trajectory))))
  rt$treatment  <- factor(rt$treatment, levels = unique(rt$treatment))
  rt
}

plot_jointwald_heatmap <- function(res, model = c("continuous", "did")) {
  d <- jointwald_long(res, model)
  if (!nrow(d)) {
    return(ggplot2::ggplot() +
             ggplot2::annotate("text", x = 0, y = 0,
                               label = "No joint Wald results to display.") +
             ggplot2::theme_void())
  }
  d$lab <- sprintf("%s\n%s", d$stars, fmt_p(d$p_value))
  ggplot2::ggplot(d, ggplot2::aes(.data$treatment, .data$trajectory,
                                  fill = .data$neglogp)) +
    ggplot2::geom_tile(colour = "white", linewidth = 1.1) +
    ggplot2::geom_text(ggplot2::aes(label = .data$lab),
                       size = 3.2, lineheight = 0.9,
                       colour = ifelse(d$neglogp > 3, "white", .cyclo$ink)) +
    ggplot2::scale_fill_gradient(low = "#EAF2F1", high = .cyclo$teal,
                                 name = expression(-log[10](p))) +
    ggplot2::labs(
      x = NULL, y = NULL,
      subtitle = "Does each treatment alter the trajectory? (joint Wald test vs. reference)") +
    theme_cyclo() +
    ggplot2::theme(legend.position = "right",
                   axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))
}

# ---- fitted-trajectory plot (wraps the package plot) ------------------------
plot_fit <- function(sub_data, model_obj, title, color) {
  plot_trajectory(sub_data, model_obj, title = title, fit_color = color) +
    theme_cyclo(12)
}

# ---- forest plot of pairwise treatment slope contrasts (continuous model) ---
# Note: the error bars are a normal-approximation interval built here for
# display from the estimate and SE that emmeans returns, whereas `p.value` is
# Tukey-adjusted by the package. The two therefore answer slightly different
# questions and can disagree at the margin; the plot caption says so.
forest_slopes_df <- function(tr) {
  s <- tr$slopes
  s$lower <- s$estimate - 1.96 * s$SE
  s$upper <- s$estimate + 1.96 * s$SE
  s$sig   <- s$p.value < 0.05
  s$contrast <- gsub("[()]", "", s$contrast)
  s[order(s$estimate), , drop = FALSE]
}

plot_forest_slopes <- function(tr) {
  d <- forest_slopes_df(tr)
  d$contrast <- factor(d$contrast, levels = d$contrast)
  ggplot2::ggplot(d, ggplot2::aes(.data$estimate, .data$contrast,
                                  colour = .data$sig)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        colour = .cyclo$grey) +
    ggplot2::geom_linerange(ggplot2::aes(xmin = .data$lower, xmax = .data$upper),
                            linewidth = 0.9) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = .cyclo$red,
                                            `FALSE` = .cyclo$grey),
                                 guide = "none") +
    ggplot2::labs(
      x = "Difference in slope (log-odds per week)", y = NULL,
      subtitle = "Pairwise treatment slope contrasts (continuous model)",
      caption = "Bars: unadjusted 95% Wald intervals. Colour: Tukey-adjusted p < 0.05.") +
    theme_cyclo(12)
}

# ---- estimated marginal means over time (DiD model) -------------------------
plot_emm <- function(tr) {
  e <- as.data.frame(tr$emmeans)
  e$week <- as.numeric(as.character(e$timepoint_cat))
  ## emmeans spells these columns differently depending on the response scale
  ## and the degrees-of-freedom method, so resolve them rather than hardcoding.
  est <- first_col(e, c("prob", "response", "emmean"))
  lo  <- first_col(e, c("asymp.LCL", "lower.CL"))
  hi  <- first_col(e, c("asymp.UCL", "upper.CL"))
  pal <- arm_palette(levels(droplevels(factor(e$treatment))))
  ggplot2::ggplot(e, ggplot2::aes(.data$week, 100 * .data[[est]],
                                  colour = .data$treatment, fill = .data$treatment)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = 100 * .data[[lo]],
                                      ymax = 100 * .data[[hi]]),
                         alpha = 0.15, colour = NA) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_colour_manual(values = pal, name = NULL) +
    ggplot2::scale_fill_manual(values = pal, name = NULL) +
    ggplot2::labs(x = "Week", y = "Model-estimated frequency (%)",
                  subtitle = "Estimated marginal means by treatment and timepoint (95% CI)") +
    theme_cyclo(12)
}

# ---- one arm's slope effect across all trajectories (forest) ---------------
# For each trajectory, the continuous-model slope difference of `arm` vs the
# reference arm (positive = arm expands faster than reference).
arm_across_traj_df <- function(res, arm) {
  ref <- levels(res$data$treatment)[1]
  # emmeans wraps contrast levels in parentheses, and we strip them from the
  # contrast strings below; strip the keys the SAME way so arms whose names
  # contain parentheses (e.g. "arm B (naive)") still match.
  strip <- function(x) gsub("[()]", "", x)
  key1 <- paste(strip(ref), "-", strip(arm))
  key2 <- paste(strip(arm), "-", strip(ref))
  rows <- lapply(names(res$trajectories), function(tv) {
    s <- res$trajectories[[tv]]$slopes
    if (is.null(s) || !nrow(s)) return(NULL)
    s$contrast <- strip(s$contrast)
    r <- s[s$contrast %in% c(key1, key2), , drop = FALSE]
    if (!nrow(r)) return(NULL)
    est <- r$estimate[1]; se <- r$SE[1]
    if (identical(r$contrast[1], key1)) est <- -est  # orient as arm - ref
    data.frame(trajectory = tv, estimate = est,
               lower = est - 1.96 * se, upper = est + 1.96 * se,
               p.value = r$p.value[1], stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (is.null(out) || !nrow(out))
    return(data.frame(trajectory = factor(character(0)), estimate = numeric(0),
                      lower = numeric(0), upper = numeric(0),
                      p.value = numeric(0), sig = logical(0)))
  out$sig <- out$p.value < 0.05
  out$trajectory <- factor(out$trajectory, levels = rev(sort(unique(out$trajectory))))
  out
}

plot_arm_across_traj <- function(res, arm) {
  d <- arm_across_traj_df(res, arm)
  ref <- levels(res$data$treatment)[1]
  if (!nrow(d))
    return(ggplot2::ggplot() +
             ggplot2::annotate("text", x = 0, y = 0,
               label = sprintf("No slope contrasts available for %s.", arm)) +
             ggplot2::theme_void())
  ggplot2::ggplot(d, ggplot2::aes(.data$estimate, .data$trajectory,
                                  colour = .data$sig)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = .cyclo$grey) +
    ggplot2::geom_linerange(ggplot2::aes(xmin = .data$lower, xmax = .data$upper),
                            linewidth = 1) +
    ggplot2::geom_point(size = 3.4) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = .cyclo$teal, `FALSE` = .cyclo$grey),
                                 guide = "none") +
    ggplot2::labs(
      x = "Slope difference vs. reference (log-odds per week)", y = NULL,
      subtitle = sprintf("%s vs. %s across trajectories: right of the line = faster expansion",
                         arm, ref),
      caption = "Bars: unadjusted 95% Wald intervals. Colour: Tukey-adjusted p < 0.05.") +
    theme_cyclo(12)
}

# shared reactable theme: a light "sheet" with dark text, readable in light AND
# dark mode (otherwise dark-mode makes the text light-on-light / foggy).
cyclo_reactable_theme <- function() {
  reactable::reactableTheme(
    color           = "#1c2b2b",
    backgroundColor = "#ffffff",
    borderColor     = "#e6ebee",
    stripedColor    = "#f6f8f9",
    highlightColor  = "#eaf2f1",
    cellPadding     = "6px 9px",
    style           = list(fontSize = "0.9rem"),
    headerStyle     = list(background = "#EEF3F5", color = "#12333a",
                           fontWeight = 700, borderColor = "#dbe4e8"),
    searchInputStyle = list(color = "#1c2b2b", background = "#fff",
                            borderColor = "#cdd8dc"),
    paginationStyle = list(color = "#1c2b2b")
  )
}

# ---- reactable of a results table with significance highlighting ------------
# `p_col` names the p-value column to format and highlight; pass NULL for tables
# that have none (e.g. marginal means, diagnostics).
signif_reactable <- function(df, p_col = "p_value", digits = 3) {
  df <- as.data.frame(df)
  ## round numeric columns for display, but leave p-values alone: fmt_p() below
  ## renders them, and rounding would turn small p-values into 0
  num <- setdiff(names(df)[vapply(df, is.numeric, logical(1))],
                 c("p_value", "p.value"))
  for (nm in num) df[[nm]] <- round(df[[nm]], digits)
  cols <- list()
  if (!is.null(p_col) && length(p_col) && p_col %in% names(df)) {
    cols[[p_col]] <- reactable::colDef(
      name = "p",
      cell = function(v) fmt_p(v),
      style = function(v) {
        if (is.na(v)) return(list())
        if (v < 0.001) list(background = "#16736B22", fontWeight = 600)
        else if (v < 0.05) list(background = "#D9A44133", fontWeight = 600)
        else list(color = "#888")
      })
  }
  reactable::reactable(
    df, columns = cols, striped = TRUE, highlight = TRUE, compact = TRUE,
    defaultPageSize = 12, searchable = nrow(df) > 12, resizable = TRUE,
    theme = cyclo_reactable_theme())
}
