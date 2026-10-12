# cyclotraj

**Mixed-model testing of treatment effects on longitudinal count trajectories.**

`cyclotraj` is an R package for a specific, recurring inference problem: you have
many *units* tracked repeatedly over time, nested inside *subjects*, where each
subject is assigned to one *treatment group*, and each observation is a **count
out of a sample-specific total**. You want to test whether the treatment groups
differ in how the units change over time.

The motivating application is immune monitoring during cancer therapy, where
the units are T-cell clones tracked in a patient's blood. Nothing in the method
is specific to that domain: the same structure arises in flow-cytometry subset
frequencies, microbiome abundances, and any repeated "successes out of trials"
measurement compared across arms. The API keeps the vocabulary of the
motivating application, so throughout it a `clone` is the unit and a `patient`
is the subject it sits in.

---

## The statistical problem

Suppose unit *i*, in subject *m*, at timepoint *j*, yields a count
$Y_{ijm}$ out of a total $N_{jm}$. Three features of this data break the
obvious analyses:

1. **The denominator varies.** Sample totals in the bundled example span a
   23-fold range (145 to 3,387). Comparing raw proportions treats a 1/145
   observation and a 24/3387 observation as equally precise. They are not.
2. **Observations are not independent.** Repeated measurements on the same unit
   are correlated, and units within a subject are correlated. Treatment is
   assigned at the *subject* level, so the effective sample size for a treatment
   effect is the number of subjects, not the number of rows.
3. **The outcome is a bounded, sparse count.** 23% of observations in the
   example are exactly zero. Log-ratio or fold-change transforms require ad hoc
   offsets; a *t*-test assumes a scale the data do not live on.

The conventional workaround is a marginal rank test (Wilcoxon) on per-unit
frequencies, one timepoint at a time. That discards the count structure, the
nesting, and the information shared across timepoints, and it gives you a
*P*-value without an effect size.

## The modeling approach

`cyclotraj` fits a binomial generalized linear mixed model with a logit link:

$$Y_{ijm}\sim\mathrm{Binomial}(N_{jm},\,p_{ijm}),\qquad
\mathrm{logit}(p_{ijm}) = x_{ijm}^{\top}\beta + u_i + v_m,$$

with random intercepts $u_i\sim N(0,\sigma_u^2)$ for unit and
$v_m\sim N(0,\sigma_v^2)$ for subject. The binomial likelihood weights each
observation by its own denominator, the logit keeps fitted values in $(0,1)$ and
handles zeros without transformation, and the random intercepts absorb the
nesting (the unit intercept also soaks up extra-binomial dispersion, so no
separate quasi-binomial scale parameter is needed).

The design matrix $x_{ijm}$ is built two ways, and **both are fit and read
together**:

| | Time encoding | Treatment effect parameter | Answers |
|---|---|---|---|
| **Continuous** | numeric | $\beta_{3k}$, a `treatment:time` interaction | Does arm *k* change at a different *rate* than the reference? |
| **Difference-in-differences** | factor, baseline = reference | $\beta_{3kt}$, one per post-baseline timepoint | Does arm *k*'s change *from baseline to time t* exceed the reference's? |

The continuous model spends a single degree of freedom on the treatment-by-time
interaction, which gives it power against sustained monotone trends but averages
non-monotone shapes away. The DiD model makes no linearity assumption and
localizes effects in time, at the cost of more parameters. Agreement between
the two corroborates a monotone trend; disagreement is itself the finding, and
localizes a non-linear one.

Inference, per group:

- **`joint_wald()`** — a multi-degree-of-freedom omnibus Wald test of all of an
  arm's coefficients at once, its baseline offset together with every
  treatment-by-time interaction (via `car::linearHypothesis`), as a screen
  before looking at individual coefficients. That is 2 df under the continuous
  encoding and 4 df under the difference-in-differences encoding here.
- **`slope_contrasts()`** — pairwise comparisons of arms' slopes
  (`emmeans::emtrends`), Tukey-adjusted.
- **`did_contrasts()`** — a contrast of contrasts: the difference between two
  arms in their baseline-to-*t* change, Tukey-adjusted.
- Every contrast carries an **effect size**, not just a *P*-value. Estimates are
  on the log-odds scale, so `exp(estimate)` is an odds ratio for the slope
  contrasts and a ratio of odds ratios for the DiD contrasts.
  `emm_means()` returns marginal proportions with 95% confidence limits and
  `predict_with_ci()` returns fitted curves with 95% bands.
- **`glmm_diagnostics()`** — per-fit convergence, singularity, and variance
  components, because with few subjects per arm the subject-level variance is
  often estimated at or near zero and you need to know that.

## Install

```r
# install.packages("remotes")
remotes::install_github("andrewlakkis/cyclotraj", subdir = "cyclotraj")
```

Requires R >= 4.1. Imports `lme4`, `emmeans`, `car`, `ggplot2`.

## Runnable example

```r
library(cyclotraj)

raw <- read_example_data()        # bundled synthetic example, 500 units x 4 timepoints
res <- run_clone_analysis(raw, treatment_ref = "arm_A")

results_table(res)                # one row per group x arm x model: df, chisq, P
tr <- res$trajectories[["Traj 3"]]
tr$joint_wald                     # omnibus test, both parameterizations
tr$slopes                         # Tukey-adjusted pairwise slope contrasts
```

The example is simulated so the two time encodings disagree somewhere, because
that disagreement is the point. In group `Traj 3`, `arm_D` is given a transient
one-timepoint spike at week 3. The difference-in-differences model, which
estimates each timepoint separately, localizes it decisively. The continuous-time
model has to express the same effect as a slope, which averages a one-timepoint
spike over the whole follow-up, so it registers the effect far more weakly on the
same data:

```
Traj 3   continuous   arm_D    chisq = 12.0    P = 0.0025
Traj 3   did          arm_D    chisq = 62.9    P = 7.3e-13
```

Group `Traj 2` is simulated with no treatment effect at all, and nothing in it
reaches significance in either model, which is the behaviour you want to be able
to confirm on data where you know the truth.

Lower-level control, one group at a time:

```r
d  <- prepare_clone_data(raw, treatment_ref = "arm_A")   # validate + standardize
d1 <- d[d$trajectory == "Traj 1", ]

m_cont <- fit_continuous(d1)
m_did  <- fit_did(d1)

joint_wald(m_cont)                                       # all non-reference arms
slope_contrasts(m_cont)
did_contrasts(m_did, alpha = 0.05)                       # only the ones that survive
glmm_diagnostics(m_cont)
```

`prepare_clone_data()` is the generic entry point. It takes one column-name
argument per field (`clone=`, `patient=`, `treatment=`, `timepoint=`, `count=`,
`total=`, `trajectory=`), so a dataset with entirely different column names maps
in without renaming anything upstream. It validates the schema and fails loudly
on the mistakes that otherwise surface as a confusing `glmer` error or, worse,
as rows silently dropped from the analysis set: non-numeric time, a missing
identifier or treatment label, missing or non-numeric counts, negative counts,
`count > total`, and a `total` of zero.

## Shiny app

`shiny_app/` is a `bslib` dashboard over the same package functions, for
exploring fits without writing R:

```r
install.packages(c("shiny", "bslib", "reactable", "fontawesome", "patchwork"))
shiny::runApp("shiny_app")
```

Tabs: an overview with a cross-group omnibus-test heatmap (continuous/DiD
toggle); a forest plot of one arm's effect across all groups; a per-group view
with both model fits, confidence bands, a slope-contrast forest plot, marginal
means, and significance-highlighted tables with CSV/PNG export; and a data tab
previewing the prepared rows. The sidebar takes a **file upload with
interactive column mapping** (CSV, TSV, Excel or `.rds`), so the app runs on a
user's own data, not just the example. The app is a thin UI layer: all
statistics come from the exported package functions, which is why the package is
tested and the app is not duplicated logic. See `shiny_app/README.md` for a
tab-by-tab description.

## Tests and checks

```sh
R CMD build cyclotraj && R CMD check cyclotraj_*.tar.gz
Rscript -e 'devtools::test("cyclotraj")'
```

Current status on R 4.5.0: **55 assertions in 13 `test_that` blocks, all
passing** (~3 s), and `R CMD check --as-cran` reports **0 errors, 0 warnings,
1 NOTE**. The NOTE is the CRAN incoming-feasibility check flagging a new
submission and a `URL:` field that does not resolve until this repository is
public.

`testthat` (edition 3) covers, for `prepare_clone_data()`, the standardized
output schema, which levels become the treatment and baseline references, and
the errors raised by a missing required column, a count exceeding its total and
an unknown `treatment_ref`, plus a custom column mapping on a data frame with
unrelated column names. The model tests go past "the code runs". Because the
example data are simulated from known parameters, they assert *recovery*: that
`slope_contrasts()` returns the simulated slope gap of 0.145 logits per week in
`Traj 1` while leaving a planted 0.015 gap non-significant, and that the DiD
coefficients in `Traj 4` put the step at the weeks it was simulated at (0 at
week 3, then 0.80 and 0.95) rather than spreading it. They also assert
*calibration*: in `Traj 2`, the group generated with no treatment effect at all,
every omnibus test and every difference-in-differences contrast is
non-significant. The joint Wald degrees of freedom (2 continuous, 4 DiD), the
variance-component ordering, and the shape and column names of the returned
tables are pinned down too.

## Repository layout

```
cyclotraj/            the R package
  R/                  prep, fit, inference, plots, pipeline
  man/                15 roxygen-generated help pages
  tests/testthat/     schema tests + model-recovery tests
  vignettes/          worked end-to-end example
  inst/extdata/       the synthetic example dataset
data-raw/             make_example_data.R -- the seeded generator for that data
shiny_app/            bslib dashboard (app.R + R/helpers.R)
.github/workflows/    R CMD check on push and pull request
```

Fourteen exported functions. `prepare_clone_data()` and `run_clone_analysis()`
are the entry points; everything else is callable on its own so the pieces can
be used without the wrapper.

## Example data

The bundled dataset is **synthetic**, generated by `data-raw/make_example_data.R`
under a fixed seed. It is simulated from the same binomial GLMM the package
fits, with parameters chosen to reproduce the structural features that make the
problem hard rather than to look tidy:

| | |
|---|---|
| 500 units, 4 timepoints | 2,000 rows, balanced panel |
| 18 subjects, 4 arms | deliberately unbalanced: 4 / 4 / 7 / 3 |
| sample totals 145 to 3,387 | 23-fold spread, unrelated to arm |
| 23% of counts exactly zero | median count 2, median frequency 0.29% |
| unit variance >> subject variance | some fits are singular, by design |

Because treatment is assigned at the subject level and there are 3 to 7 subjects
per arm, the subject-level variance is often estimated at or near zero: it is
numerically zero in 8 of the 12 fits the pipeline produces on this data, and 9
of the 12 are flagged singular. That is a real property of this design, not a simulation
artifact, and it is why `glmm_diagnostics()` exists and why the vignette reads
the two models together instead of trusting one.

No real subject data is distributed with this package.

## License

MIT. See `LICENSE`.
