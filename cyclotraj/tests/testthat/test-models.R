## The bundled example data are simulated from known parameters by
## data-raw/make_example_data.R, so the tests below assert that the pipeline
## recovers the effects that were actually put in, and finds nothing in the
## group that was simulated with no treatment effects at all.
##
## Fits are memoised because each test_that() block runs in its own frame: a
## naive helper would refit the same GLMM several times.
.cache <- new.env(parent = emptyenv())
cached <- function(key, expr) {
  if (!exists(key, envir = .cache, inherits = FALSE)) {
    assign(key, suppressWarnings(suppressMessages(expr)), envir = .cache)
  }
  get(key, envir = .cache, inherits = FALSE)
}

## one trajectory of the prepared data
prep_traj <- function(label = "Traj 1") {
  d <- cached("prepared",
              prepare_clone_data(read_example_data(), treatment_ref = "arm_A"))
  cached(paste0("data_", label), d[d$trajectory == label, , drop = FALSE])
}

fit_cached <- function(label, type = c("continuous", "did")) {
  type <- match.arg(type)
  cached(paste0(type, "_", label),
         if (type == "continuous") fit_continuous(prep_traj(label))
         else fit_did(prep_traj(label)))
}

## kept for the pre-existing tests below
prep_one <- function() prep_traj("Traj 1")


test_that("both GLMMs fit and carry a time attribute", {
  mc <- fit_cached("Traj 1", "continuous")
  md <- fit_cached("Traj 1", "did")
  expect_s4_class(mc, "glmerMod")
  expect_s4_class(md, "glmerMod")
  expect_identical(attr(mc, "cyclotraj_time"), "timepoint_num")
  expect_identical(attr(md, "cyclotraj_time"), "timepoint_cat")
})

test_that("joint_wald returns one row per non-reference treatment with valid p", {
  sub <- prep_one()
  mc <- fit_cached("Traj 1", "continuous")
  jw <- joint_wald(mc)
  expect_equal(nrow(jw), nlevels(sub$treatment) - 1L)
  expect_true(all(jw$p_value >= 0 & jw$p_value <= 1))
  ## Traj 1: arm_D is simulated with a strong slope effect, so the
  ## omnibus test should detect it (see data-raw/make_example_data.R)
  expect_lt(jw$p_value[jw$treatment == "arm_D"], 0.05)
})

test_that("predict_with_ci returns ordered, in-range bands", {
  sub <- prep_one()
  mc <- fit_cached("Traj 1", "continuous")
  pc <- predict_with_ci(mc, sub)
  expect_true(all(pc$lower <= pc$fitted + 1e-8))
  expect_true(all(pc$fitted <= pc$upper + 1e-8))
  expect_true(all(pc$fitted >= 0 & pc$fitted <= 100))
})

test_that("plot_trajectory returns a ggplot", {
  sub <- prep_one()
  mc <- fit_cached("Traj 1", "continuous")
  p <- plot_trajectory(sub, mc)
  expect_s3_class(p, "ggplot")
})

test_that("slope_contrasts recovers the simulated slope gap (Traj 1)", {
  ## Traj 1 is simulated with a reference slope of -0.075 per week and an
  ## extra -0.145 per week for arm_D, so the arm_A - arm_D slope contrast
  ## should land near +0.145 and be clearly significant. arm_B and arm_C get
  ## extra slopes of -0.030 and -0.015, a gap of only 0.015, which should not
  ## be detectable with this many units.
  sc <- suppressWarnings(slope_contrasts(fit_cached("Traj 1", "continuous")))
  expect_equal(nrow(sc), choose(4, 2))
  expect_true(all(c("contrast", "estimate", "p.value") %in% names(sc)))

  ad <- sc[sc$contrast == "arm_A - arm_D", ]
  expect_equal(nrow(ad), 1L)
  ## tolerance is relative here: 0.2 allows +/- 0.029, about one standard
  ## error of the estimate
  expect_equal(ad$estimate, 0.145, tolerance = 0.2)
  expect_lt(ad$p.value, 0.01)

  bc <- sc[sc$contrast == "arm_B - arm_C", ]
  expect_equal(bc$estimate, -0.015, tolerance = 0.05)
  expect_gt(bc$p.value, 0.05)
})

test_that("the DiD model places the simulated effect at the right weeks", {
  ## Traj 4 gives arm_B a delayed step of +0.80 logits at week 6 and +0.95 at
  ## week 9, with nothing at week 3, and gives arm_C and arm_D no effect at
  ## all. The DiD parameterization should place the effect at the right weeks.
  md <- fit_cached("Traj 4", "did")
  b <- lme4::fixef(md)
  step <- function(w) unname(b[paste0("treatmentarm_B:timepoint_cat", w)])
  expect_equal(step(3), 0.00, tolerance = 0.3)
  expect_equal(step(6), 0.80, tolerance = 0.3)
  expect_equal(step(9), 0.95, tolerance = 0.3)

  jw_did <- joint_wald(md)
  expect_equal(unique(jw_did$df), 4)
  expect_lt(jw_did$p_value[jw_did$treatment == "arm_B"], 0.001)

  ## the continuous model also flags arm_B (a step has a non-zero average
  ## slope), while the two arms with no simulated effect stay non-significant
  jw_cont <- joint_wald(fit_cached("Traj 4", "continuous"))
  expect_lt(jw_cont$p_value[jw_cont$treatment == "arm_B"], 0.001)
  expect_true(all(jw_cont$p_value[jw_cont$treatment != "arm_B"] > 0.05))

  ## categorical branch of predict_with_ci: one row per arm x observed week
  pc <- predict_with_ci(md, prep_traj("Traj 4"))
  expect_equal(nrow(pc), 4L * 4L)
  expect_identical(sort(unique(pc$timepoint_num)), c(0, 3, 6, 9))
})

test_that("the full pipeline runs end to end and returns a tidy table", {
  res <- cached("analysis_traj2", {
    raw <- read_example_data()
    run_clone_analysis(raw[raw$k == "Traj 2", ], treatment_ref = "arm_A")
  })
  expect_s3_class(res, "cyclotraj_analysis")
  expect_identical(names(res$trajectories), "Traj 2")

  rt <- results_table(res)
  expect_identical(names(rt), c("trajectory", "model", "treatment",
                                "df", "chisq", "p_value"))
  ## both parameterizations x each non-reference arm
  expect_equal(nrow(rt), 2L * 3L)
  expect_setequal(unique(rt$model), c("continuous", "did"))

  tr <- res$trajectories[["Traj 2"]]
  expect_equal(nrow(tr$emmeans), 4L * 4L)
  expect_equal(nrow(tr$diagnostics), 2L)
  ## the generator sets between-unit SD to 0.95 and between-subject SD to
  ## 0.06, so the unit intercept should carry nearly all the variance
  sdc <- tr$diagnostics$sd_clone
  expect_true(all(sdc > 0.6 & sdc < 1.3))
  expect_true(all(tr$diagnostics$sd_clone > tr$diagnostics$sd_patient))
})

test_that("joint_wald finds nothing in the group simulated as a true null", {
  ## Traj 2 is the one group given no arm-specific slopes and no bumps, so
  ## every omnibus test is being evaluated under the null hypothesis.
  res <- cached("analysis_traj2", {
    raw <- read_example_data()
    run_clone_analysis(raw[raw$k == "Traj 2", ], treatment_ref = "arm_A")
  })
  rt <- results_table(res)
  expect_equal(rt$df[rt$model == "continuous"], rep(2, 3))
  expect_equal(rt$df[rt$model == "did"], rep(4, 3))
  expect_true(all(rt$p_value > 0.1))

  ## and no post-hoc difference-in-differences contrast survives either
  tr <- res$trajectories[["Traj 2"]]
  expect_equal(sum(tr$did_contrasts$p.value < 0.05), 0L)
  expect_equal(nrow(suppressWarnings(did_contrasts(tr$did, alpha = 0.05))), 0L)
})
