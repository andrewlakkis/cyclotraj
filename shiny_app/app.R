# =============================================================================
#  cyclotraj explorer  --  an interactive dashboard for the cyclotraj pipeline.
#
#  Fits the two binomial GLMMs within each trajectory group and lets you explore
#  treatment effects: a cross-trajectory overview (which arms modulate which
#  trajectories), a per-trajectory deep dive (fitted curves, forest plots,
#  marginal means, and every results table with downloads), on the bundled
#  synthetic example or any uploaded long-format data.
#
#  Every model fit, test and contrast comes from cyclotraj's exported functions;
#  this file and R/helpers.R only arrange and draw their output.
#
#  Run:  shiny::runApp("shiny_app")   from the repository root.
#  Helper plot/table functions live in R/helpers.R (auto-loaded by Shiny).
# =============================================================================

suppressMessages({
  library(shiny)
  library(bslib)
  library(ggplot2)
  library(reactable)
  library(fontawesome)
  library(patchwork)
})

## the default 5 MB cap is smaller than the clone tables this app is meant for
options(shiny.maxRequestSize = 100 * 1024^2)

## ---- make the cyclotraj package available ----------------------------------
if (!requireNamespace("cyclotraj", quietly = TRUE)) {
  pkg_dir <- file.path("..", "cyclotraj")
  if (dir.exists(pkg_dir) && requireNamespace("devtools", quietly = TRUE)) {
    suppressMessages(devtools::load_all(pkg_dir, quiet = TRUE))
  } else stop("Install 'cyclotraj', or run the app from the repository root ",
              "so that ../cyclotraj can be loaded with devtools.")
} else library(cyclotraj)

DEFAULT_COLS <- list(clone = "clone_id", patient = "patient_alias",
                     treatment = "treatment", timepoint = "timepoint",
                     count = "n", total = "n_tcells_tp", trajectory = "k")

theme_app <- bs_theme(
  version = 5, primary = "#16736B", secondary = "#2166AC",
  "navbar-bg" = "#16736B"
)

# =============================================================================
# UI
# =============================================================================
ui <- page_navbar(
  title = tagList(
    tags$img(src = "favicon.svg", height = "24px",
             style = "margin-right:7px;vertical-align:text-bottom;"),
    "cyclotraj"),
  window_title = "cyclotraj · clonal trajectory explorer",
  theme = theme_app,
  fillable = TRUE,
  header = tags$head(
    tags$link(rel = "icon", type = "image/svg+xml", href = "favicon.svg")),

  sidebar = sidebar(
    width = 300, title = "Data & settings",
    radioButtons("source", NULL,
      c("Example data" = "example", "Upload file" = "upload")),
    conditionalPanel("input.source == 'upload'",
      fileInput("file", "Long-format data file",
                accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls", ".rds")),
      helpText(class = "text-muted small",
               "CSV, TSV/TXT, or .rds (a data.frame).",
               "Excel (.xlsx/.xls) also works if the readxl package is installed."),
      accordion(open = FALSE,
        accordion_panel("Column mapping", icon = fa("table-columns"),
          helpText("Pick which of your file's columns maps to each field cyclotraj expects."),
          uiOutput("colmap_ui")))),
    uiOutput("ref_ui"),
    input_switch("drop_patient", "Drop patient random effect", FALSE),
    hr(),
    uiOutput("status_ui"),
    hr(),
    div(class = "text-muted small",
        "Two binomial GLMMs are fit per trajectory: a continuous-time",
        "(slope) model and a difference-in-differences (timepoint) model.")
  ),

  # ---- Overview -------------------------------------------------------------
  nav_panel("Overview", icon = fa("gauge-high"),
    layout_columns(
      fill = FALSE, gap = "0.6rem",
      value_box("Patients", textOutput("vb_patients"), showcase = fa("users"),
                theme = "primary"),
      value_box("Clones", textOutput("vb_clones"), showcase = fa("dna"),
                theme = "secondary"),
      value_box("Trajectories", textOutput("vb_traj"), showcase = fa("chart-line"),
                theme = value_box_theme(bg = "#D9A441", fg = "white")),
      value_box("Treatment arms", textOutput("vb_arms"), showcase = fa("flask-vial"),
                theme = value_box_theme(bg = "#6A51A3", fg = "white"))
    ),
    layout_columns(
      col_widths = c(7, 5),
      card(full_screen = TRUE,
        card_header(
          div(class = "d-flex justify-content-between align-items-center",
              span(fa("border-all"), " Which treatments modulate which trajectories?"),
              radioButtons("hm_model", NULL, inline = TRUE,
                c("Continuous" = "continuous", "Diff-in-diff" = "did")))),
        plotOutput("heatmap", height = "440px"),
        downloadButton("dl_heatmap", "PNG", class = "btn-sm")),
      card(
        card_header(fa("lightbulb"), " Reading this"),
        uiOutput("overview_note"))
    )
  ),

  # ---- Across trajectories --------------------------------------------------
  nav_panel("Across trajectories", icon = fa("layer-group"),
    card(full_screen = TRUE,
      card_header(
        div(class = "d-flex justify-content-between align-items-center",
            span(fa("layer-group"), " How does one treatment act across trajectories?"),
            uiOutput("arm_select_ui"))),
      layout_columns(col_widths = c(8, 4),
        plotOutput("across_forest", height = "440px"),
        card(card_header(fa("lightbulb"), " Reading this"),
             uiOutput("across_note"))))
  ),

  # ---- Explore a trajectory -------------------------------------------------
  nav_panel("Explore a trajectory", icon = fa("magnifying-glass-chart"),
    layout_columns(
      fill = FALSE, col_widths = c(4, 8),
      card(card_body(uiOutput("traj_ui"), class = "py-2")),
      card(card_body(uiOutput("traj_verdict"), class = "py-2"))
    ),
    navset_card_tab(
      full_screen = TRUE,
      nav_panel("Fitted trajectories", icon = fa("chart-area"),
        layout_columns(
          plotOutput("fit_cont", height = "380px"),
          plotOutput("fit_did",  height = "380px")),
        downloadButton("dl_fits", "Download figure (PNG)", class = "btn-sm")),
      nav_panel("Effect sizes", icon = fa("arrows-left-right-to-line"),
        layout_columns(
          plotOutput("forest", height = "330px"),
          plotOutput("emm", height = "330px"))),
      nav_panel("Joint Wald", icon = fa("filter"),
        card_body(reactableOutput("tbl_wald"),
                  downloadButton("dl_wald", "CSV", class = "btn-sm"))),
      nav_panel("Slope contrasts", icon = fa("angle-right"),
        card_body(reactableOutput("tbl_slopes"),
                  downloadButton("dl_slopes", "CSV", class = "btn-sm"))),
      nav_panel("Diff-in-diff", icon = fa("timeline"),
        card_body(reactableOutput("tbl_did"),
                  downloadButton("dl_did", "CSV", class = "btn-sm"))),
      nav_panel("Marginal means", icon = fa("chart-simple"),
        card_body(reactableOutput("tbl_emm"),
                  downloadButton("dl_emm", "CSV", class = "btn-sm"))),
      nav_panel("Diagnostics", icon = fa("stethoscope"),
        card_body(reactableOutput("tbl_diag")))
    )
  ),

  # ---- Data -----------------------------------------------------------------
  nav_panel("Data", icon = fa("table"),
    card(full_screen = TRUE,
      card_header(fa("table"), " Prepared data (one row per clone-timepoint)"),
      card_body(reactableOutput("data_preview")))),

  # ---- About ----------------------------------------------------------------
  nav_panel("About", icon = fa("circle-info"),
    card(card_body(uiOutput("about")))),

  nav_spacer(),
  nav_item(input_dark_mode(id = "dark_mode")),
  nav_item(tags$a(href = "https://github.com/andrewlakkis/cyclotraj", target = "_blank",
                  fa("github"), " code"))
)

# =============================================================================
# SERVER
# =============================================================================
server <- function(input, output, session) {

  # Which column of the input data plays each role. For an upload these come
  # from the mapping dropdowns, which do not exist until the file has been read,
  # so fall back to the canonical names until they do (that is also what the
  # dropdowns preselect) rather than briefly reporting "column not found".
  cols <- reactive({
    if (input$source == "example") return(DEFAULT_COLS)
    g <- function(x, default) if (is.null(x)) default else x
    list(clone      = g(input$c_clone,     DEFAULT_COLS$clone),
         patient    = g(input$c_patient,   DEFAULT_COLS$patient),
         treatment  = g(input$c_treatment, DEFAULT_COLS$treatment),
         timepoint  = g(input$c_timepoint, DEFAULT_COLS$timepoint),
         count      = g(input$c_count,     DEFAULT_COLS$count),
         total      = g(input$c_total,     DEFAULT_COLS$total),
         # "" is the "(none)" choice: no trajectory column, one group overall
         trajectory = g(input$c_traj,      DEFAULT_COLS$trajectory))
  })

  raw <- reactive({
    if (input$source == "example") read_example_data()
    else {
      req(input$file)
      tryCatch(read_upload(input$file$datapath, input$file$name),
               error = function(e)
                 validate(need(FALSE, paste("Could not read file:", conditionMessage(e)))))
    }
  })

  # The reference arm every effect is reported against. input$ref arrives one
  # flush after the selector is rendered, and is stale for one flush after a new
  # upload, so resolve it to the same default the selector shows. Keying the
  # cached fits on this resolved value (not on input$ref) is what stops the
  # whole GLMM run from happening twice at start-up.
  ref_level <- reactive({
    lv <- treatment_levels(raw(), cols()$treatment)
    if (!length(lv)) return(NULL)
    if (!is.null(input$ref) && input$ref %in% lv) input$ref else lv[[1]]
  })

  # fast: prepared data (drives value boxes, trajectory list, data tab)
  prepared <- reactive({
    df <- raw(); cc <- cols()
    validate(need(cc$treatment %in% names(df),
                  "Treatment column not found; check the column mapping."))
    d <- tryCatch(
      prepare_clone_data(df, clone = cc$clone, patient = cc$patient,
        treatment = cc$treatment, timepoint = cc$timepoint, count = cc$count,
        total = cc$total, trajectory = cc$trajectory,
        treatment_ref = ref_level()),
      error = function(e) validate(need(FALSE, paste("Data problem:", conditionMessage(e)))))
    # the two conditions the design cannot do without; checking them here gives
    # a plain message instead of an lme4 error from deep inside the fit
    validate(
      need(nlevels(d$treatment) > 1,
           "Only one treatment arm in these data; at least two are needed to compare arms."),
      need(nlevels(d$timepoint_cat) > 1,
           "Only one timepoint in these data; at least two are needed to model change over time."))
    d
  })

  # slow: fit all trajectories once, cached by data + settings
  analysis <- reactive({
    d <- prepared()
    withProgress(message = "Fitting GLMMs across trajectories...", value = 0.3, {
      tryCatch(
        suppressMessages(suppressWarnings(
          run_clone_analysis(d, prepare = FALSE, drop_patient = input$drop_patient))),
        error = function(e) validate(need(FALSE, paste0(
          "Model fitting failed: ", conditionMessage(e),
          " (a group needs units in at least two arms, each observed at two or more timepoints)"))))
    })
  }) |> bindCache(input$source,
                  # datapath is a fresh temp file per upload, so re-uploading an
                  # edited file with the same name is not served from the cache
                  if (isTRUE(input$source == "upload")) input$file$datapath else "ex",
                  ref_level(), input$drop_patient, cols())

  # column-mapping dropdowns, auto-populated from the uploaded CSV's columns
  output$colmap_ui <- renderUI({
    req(input$source == "upload", input$file)
    cn <- names(raw())
    validate(need(length(cn) > 0, "Could not read column names from the file."))
    pick <- function(default) if (default %in% cn) default else cn[[1]]
    tagList(
      selectInput("c_clone",     "Clone ID",   cn, selected = pick(DEFAULT_COLS$clone)),
      selectInput("c_patient",   "Patient ID", cn, selected = pick(DEFAULT_COLS$patient)),
      selectInput("c_treatment", "Treatment",  cn, selected = pick(DEFAULT_COLS$treatment)),
      selectInput("c_timepoint", "Timepoint",  cn, selected = pick(DEFAULT_COLS$timepoint)),
      selectInput("c_count",     "Count",      cn, selected = pick(DEFAULT_COLS$count)),
      selectInput("c_total",     "Total",      cn, selected = pick(DEFAULT_COLS$total)),
      # a trajectory/cluster label is optional: "(none)" analyzes the whole
      # data set as a single group
      selectInput("c_traj", "Trajectory", c("(none)" = "", cn),
                  selected = if (DEFAULT_COLS$trajectory %in% cn)
                    DEFAULT_COLS$trajectory else ""),
      helpText(class = "text-muted small",
               "Rows whose trajectory label is \"Other\" are excluded as",
               "unselected units, which is cyclotraj's default.")
    )
  })

  output$ref_ui <- renderUI({
    ## default to the first level alphabetically; the user can change it, and
    ## every effect in the app is reported relative to whatever is chosen here
    lv <- treatment_levels(raw(), cols()$treatment)
    if (!length(lv)) return(NULL)
    selectInput("ref", "Reference treatment", lv, selected = lv[[1]])
  })

  output$status_ui <- renderUI({
    d <- prepared()
    tags$div(class = "small",
      tags$b("Loaded: "), sprintf("%d clones, %d patients, %d arms",
        nlevels(d$clone_id), nlevels(d$patient_id), nlevels(d$treatment)),
      br(), tags$b("Reference: "), levels(d$treatment)[1])
  })

  # ---- value boxes ----
  output$vb_patients <- renderText(nlevels(prepared()$patient_id))
  output$vb_clones   <- renderText(format(nlevels(prepared()$clone_id), big.mark = ","))
  output$vb_traj     <- renderText(
    if ("trajectory" %in% names(prepared())) nlevels(droplevels(prepared()$trajectory)) else 1)
  output$vb_arms     <- renderText(nlevels(prepared()$treatment))

  # ---- overview heatmap ----
  # one reactive shared by the plot and its download so the two cannot drift
  heatmap_plot <- reactive(plot_jointwald_heatmap(analysis(), input$hm_model))
  output$heatmap <- renderPlot(heatmap_plot(), res = 100)

  output$overview_note <- renderUI({
    res <- analysis()
    rt <- results_table(res)
    validate(need(!is.null(rt) && nrow(rt) > 0,
                  "No joint Wald results available for these data."))
    ## which() so that a treatment whose test returned NA is not counted
    sig <- rt[which(rt$model == "continuous" & rt$p_value < 0.05), , drop = FALSE]
    top <- sort(table(sig$treatment), decreasing = TRUE)
    dom <- if (length(top)) names(top)[1] else "no arm"
    ntr <- if (length(top)) top[1] else 0
    HTML(sprintf(
      "<p>Each cell asks whether a treatment arm's trajectory differs from the
       <b>reference arm</b>'s (%s) in level and/or in shape, via a joint Wald
       test over that arm's coefficients. Darker and starred cells are stronger
       effects.</p>
       <p>In this dataset, <b>%s</b> is the most pervasive modulator, significantly
       altering <b>%d</b> of %d trajectories under the continuous-time model.</p>
       <p>Switch the toggle to compare the continuous-time (overall rate) and
       difference-in-differences (timepoint-specific) views. Open a trajectory in
       <i>Explore</i> to see the fitted curves and effect sizes behind each cell.</p>",
      levels(res$data$treatment)[1], dom, as.integer(ntr),
      length(unique(rt$trajectory))))
  })

  output$dl_heatmap <- downloadHandler(
    filename = function() sprintf("cyclotraj_heatmap_%s.png", input$hm_model),
    content = function(f) ggsave(f, heatmap_plot(), width = 7, height = 4.6,
                                 dpi = 150))

  # ---- across trajectories ----
  output$arm_select_ui <- renderUI(
    selectInput("across_arm", NULL, choices = levels(prepared()$treatment)[-1],
                width = "220px"))
  output$across_forest <- renderPlot({
    req(input$across_arm); plot_arm_across_traj(analysis(), input$across_arm)
  }, res = 100)
  output$across_note <- renderUI({
    req(input$across_arm)
    d <- arm_across_traj_df(analysis(), input$across_arm)
    up <- sum(d$sig & d$estimate > 0); dn <- sum(d$sig & d$estimate < 0)
    HTML(sprintf(
      "<p>Each row is the continuous-time slope difference of <b>%s</b> versus the
       reference arm in one trajectory. Points to the right of the dashed line mean
       %s drives faster expansion there; points to the left mean faster
       contraction.</p><p><b>%s</b> significantly accelerates <b>%d</b> and dampens
       <b>%d</b> trajectories relative to the reference.</p>",
      input$across_arm, input$across_arm, input$across_arm, up, dn))
  })

  # ---- per-trajectory ----
  output$traj_ui <- renderUI({
    d <- prepared()
    if (!"trajectory" %in% names(d))
      return(helpText("No trajectory column; whole dataset treated as one group."))
    tags$div(
      tags$label(class = "form-label fw-bold", fa("chart-line"), " Trajectory"),
      selectInput("traj", NULL, levels(droplevels(d$trajectory)), width = "100%",
                  selectize = FALSE))
  })

  cur_traj <- reactive({
    d <- prepared()
    if (!"trajectory" %in% names(d)) return("all")
    ## ignore a selection left over from the previous data set (input$traj lags
    ## the picker by a flush) instead of looking up a trajectory that is gone
    req(input$traj %in% levels(droplevels(d$trajectory)))
    input$traj
  })
  tr <- reactive({
    out <- analysis()$trajectories[[cur_traj()]]
    validate(need(!is.null(out), "This group could not be fitted."))
    out
  })
  sub_data <- reactive({
    d <- analysis()$data
    if ("trajectory" %in% names(d)) d[d$trajectory == cur_traj(), , drop = FALSE] else d
  })

  output$traj_verdict <- renderUI({
    jw <- tr()$joint_wald$continuous
    sig <- jw[which(jw$p_value < 0.05), , drop = FALSE]
    if (!nrow(sig)) return(HTML(paste0("<span class='text-muted'>",
      fa("circle-check"), " No treatment significantly alters this trajectory.</span>")))
    chips <- paste(sprintf(
      "<span class='badge rounded-pill' style='background:#16736B'>%s (p=%s)</span>",
      sig$treatment, fmt_p(sig$p_value)), collapse = " ")
    HTML(paste0("<b>Significant vs. reference (continuous model):</b><br>", chips))
  })

  # built once each so the on-screen panels and the downloaded figure are the
  # same plots
  fit_plots <- reactive(list(
    continuous = plot_fit(sub_data(), tr()$continuous,
                          "Continuous-time fit", "#2166AC"),
    did = plot_fit(sub_data(), tr()$did,
                   "Difference-in-differences fit", "#B2182B")))
  output$fit_cont <- renderPlot(fit_plots()$continuous, res = 100)
  output$fit_did  <- renderPlot(fit_plots()$did, res = 100)
  output$forest <- renderPlot(plot_forest_slopes(tr()), res = 100)
  output$emm    <- renderPlot(plot_emm(tr()), res = 100)

  # ---- tables + downloads ----
  wald_df <- reactive(rbind(
    cbind(model = "continuous", tr()$joint_wald$continuous),
    cbind(model = "did",        tr()$joint_wald$did)))
  output$tbl_wald   <- renderReactable(signif_reactable(wald_df(), "p_value"))
  output$tbl_slopes <- renderReactable(signif_reactable(tr()$slopes, "p.value"))
  output$tbl_did    <- renderReactable({
    d <- tr()$did_contrasts; signif_reactable(d[order(d$p.value), ], "p.value") })
  output$tbl_emm    <- renderReactable(signif_reactable(tr()$emmeans, NULL))
  output$tbl_diag   <- renderReactable(signif_reactable(tr()$diagnostics, NULL))

  dl <- function(name, data_fun)
    downloadHandler(
      filename = function() sprintf("%s_%s.csv", cur_traj(), name),
      content = function(f) utils::write.csv(data_fun(), f, row.names = FALSE))
  output$dl_wald   <- dl("joint_wald", wald_df)
  output$dl_slopes <- dl("slope_contrasts", reactive(tr()$slopes))
  output$dl_did    <- dl("did_contrasts", reactive(tr()$did_contrasts))
  output$dl_emm    <- dl("marginal_means", reactive(tr()$emmeans))
  output$dl_fits   <- downloadHandler(
    filename = function() sprintf("%s_fits.png", cur_traj()),
    content = function(f) {
      p <- patchwork::wrap_plots(fit_plots()$continuous, fit_plots()$did, ncol = 2)
      ggsave(f, p, width = 13, height = 5, dpi = 150)
    })

  # ---- data preview ----
  output$data_preview <- renderReactable({
    d <- as.data.frame(prepared())
    d$perc <- round(d$perc, 3)
    reactable::reactable(d, striped = TRUE, compact = TRUE, filterable = TRUE,
      defaultPageSize = 15, showPageSizeOptions = TRUE,
      theme = cyclo_reactable_theme())
  })

  # ---- about ----
  output$about <- renderUI(HTML(
    "<h4>About cyclotraj</h4>
     <p><b>cyclotraj</b> tests whether treatment groups differ in the longitudinal
     dynamics of tracked clones (e.g. the trajectories produced by the
     <i>Cyclone</i> algorithm). Within each trajectory it fits two binomial
     generalized linear mixed-effects models (logit link; random intercepts for
     clone and patient):</p>
     <ul>
       <li><b>Continuous-time model</b> &mdash; time as a numeric covariate; the
           treatment&times;time interaction is each arm's difference in slope
           (overall rate of expansion or contraction) versus the reference arm.</li>
       <li><b>Difference-in-differences model</b> &mdash; time as a factor; the
           treatment&times;timepoint interaction is how much more an arm changes
           from baseline to a given timepoint than the reference arm.</li>
     </ul>
     <p>Reading the two together diagnoses the <i>shape</i> of a response:
     agreement indicates a monotone trend, disagreement reveals non-linear
     dynamics. Inference uses a joint Wald test per arm plus Tukey-adjusted
     post-hoc contrasts. Contrast estimates are on the log-odds scale, so
     exponentiating one gives an odds ratio.</p>
     <p>Built on <code>lme4</code> and <code>emmeans</code>. This explorer is a
     thin interface over the cyclotraj R package: every fit, test and contrast
     shown here comes from the package's exported functions.</p>"))
}

shinyApp(ui, server)
