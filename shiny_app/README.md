# cyclotraj explorer (Shiny app)

An interactive dashboard over the `cyclotraj` package. It fits both binomial
GLMMs within each trajectory group and lets you explore the resulting treatment
effects visually, on the bundled synthetic example or on your own uploaded
long-format data.

All statistics come from the package's exported functions. This app contains no
modeling logic of its own, which is why the package is unit-tested and the app
is not a second implementation to keep in sync.

## Run

```r
install.packages(c("shiny", "bslib", "reactable", "fontawesome", "patchwork"))
shiny::runApp("shiny_app")      # from the repository root
```

If `cyclotraj` is not installed, the app loads it from `../cyclotraj` with
`devtools::load_all()`, so it runs from a fresh clone without installing the
package first (that path needs `devtools` as well).

## What's in it

- **Overview.** Four value boxes counting patients (subjects), clones (units),
  trajectories (groups) and treatment arms, and a cross-group heatmap: for every
  arm and trajectory group, the joint Wald test of whether that arm alters the
  trajectory relative to the reference. Toggle between the continuous-time and
  difference-in-differences views; download as PNG.
- **Across trajectories.** Pick an arm and see a forest plot of its slope
  difference versus the reference, in log-odds per week, across all groups
  (right of the line is faster expansion, left is faster contraction), showing
  how one treatment acts differently depending on the underlying trajectory
  shape.
- **Explore a trajectory.** For a selected group, in tabs: fitted curves from
  both models with 95% confidence bands (downloadable as PNG), a forest plot of
  pairwise slope contrasts alongside estimated marginal means over time by arm,
  and sortable, filterable tables with significance highlighting for the joint
  Wald tests, slope contrasts, difference-in-differences contrasts, marginal
  means (each downloadable as CSV) and the per-fit diagnostics.
- **Data.** A filterable preview of the prepared data, one row per clone and
  timepoint.
- **About.** A short description of the two models and the inference.

A light/dark mode toggle sits in the navbar.

## Using your own data

Choose "Upload file" in the sidebar and provide a long-format table, one row
per unit and timepoint, as CSV, TSV/TXT, an `.rds` holding a data frame, or
Excel (`.xlsx`/`.xls`, which needs the optional `readxl` package). If your column names differ from the defaults, open the
**Column mapping** panel below it and map each field (unit, subject, treatment,
timepoint, count, total, trajectory), then set the reference arm. The "Drop
patient random effect" switch refits without the subject-level random intercept,
which is useful when that variance component is estimated at or near zero and
the fit is singular.

## Files

- `app.R` — UI and server (`shiny`, `bslib`, `reactable`, `fontawesome`).
- `R/helpers.R` — plotting and table helpers, auto-loaded by Shiny.
