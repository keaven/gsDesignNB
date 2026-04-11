parse_numeric_csv <- function(x) {
  parts <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  parts <- parts[nzchar(parts)]
  vals <- suppressWarnings(as.numeric(parts))
  if (length(vals) == 0L || anyNA(vals)) {
    stop("Analysis times must be a comma-separated list of numbers.", call. = FALSE)
  }
  vals
}

event_gap_in_unit <- function(event_gap_days, time_unit) {
  if (identical(time_unit, "years")) {
    event_gap_days / 365.25
  } else {
    event_gap_days / 30.4375
  }
}

round_numeric_df <- function(x, digits = 3) {
  out <- x
  is_num <- vapply(out, is.numeric, logical(1))
  out[is_num] <- lapply(out[is_num], round, digits = digits)
  out
}

metric_choices <- c(
  "Expected sample size" = "mean_n",
  "Expected participants with events" = "expected_participants_with_events",
  "Expected events observed" = "expected_events_observed",
  "Rejection rate" = "rejection_rate",
  "Futility rate" = "futility_rate"
)

ui <- shiny::fluidPage(
  shiny::titlePanel("gsDesignNB SSR Explorer"),
  shiny::sidebarLayout(
    shiny::sidebarPanel(
      shiny::helpText(
        "Prototype app: the interface is exploratory, but all calculations use",
        "package functions such as sample_size_nbinom(), gsNBCalendar(),",
        "sim_ssr_nbinom(), and summarize_ssr_sim()."
      ),
      shiny::helpText(
        "Keep rates, accrual durations, trial duration, follow-up, and analysis",
        "times on the same time scale. Event-gap days are converted to the",
        "selected time unit."
      ),
      shiny::selectInput(
        "time_unit",
        "Primary time unit",
        choices = c("Months" = "months", "Years" = "years"),
        selected = "months"
      ),
      shiny::tags$hr(),
      shiny::h4("Planned design"),
      shiny::numericInput("lambda1_plan", "Planned control rate", value = 0.5, min = 0.001, step = 0.05),
      shiny::numericInput("rr_plan", "Planned rate ratio", value = 0.7, min = 0.05, max = 2, step = 0.05),
      shiny::numericInput("k_plan", "Planned dispersion", value = 0.5, min = 0, step = 0.1),
      shiny::numericInput("power_plan", "Target power", value = 0.9, min = 0.5, max = 0.999, step = 0.01),
      shiny::numericInput("alpha_plan", "One-sided alpha", value = 0.025, min = 0.0001, max = 0.25, step = 0.005),
      shiny::numericInput("accrual_rate_plan", "Planned accrual rate", value = 30, min = 0.1, step = 1),
      shiny::numericInput("accrual_duration", "Planned accrual duration", value = 12, min = 0.1, step = 1),
      shiny::numericInput("max_followup", "Maximum follow-up", value = 12, min = 0.1, step = 1),
      shiny::textInput("analysis_times", "Analysis times (comma-separated)", value = "9, 14, 24"),
      shiny::numericInput("upper_spending", "Upper spending parameter (HSD)", value = -2, step = 0.5),
      shiny::numericInput("lower_spending", "Lower spending parameter (HSD)", value = 1, step = 0.5),
      shiny::tags$hr(),
      shiny::h4("True scenario"),
      shiny::numericInput("lambda1_true", "True control rate", value = 0.3, min = 0.001, step = 0.05),
      shiny::numericInput("rr_true", "True rate ratio", value = 0.7, min = 0.05, max = 2, step = 0.05),
      shiny::numericInput("k_true", "True dispersion", value = 1.0, min = 0, step = 0.1),
      shiny::numericInput("accrual_true", "True accrual rate", value = 18, min = 0.1, step = 1),
      shiny::tags$hr(),
      shiny::h4("Simulation"),
      shiny::numericInput("n_sims", "Simulation replicates", value = 200, min = 10, step = 10),
      shiny::checkboxGroupInput(
        "strategies",
        "Strategies",
        choices = c("No adaptation", "Blinded SSR", "Unblinded SSR"),
        selected = c("No adaptation", "Blinded SSR", "Unblinded SSR")
      ),
      shiny::selectInput(
        "bound_info",
        "Information for bounds",
        choices = c(
          "Unblinded ML" = "unblinded_ml",
          "Blinded ML" = "blinded_ml",
          "Unblinded MoM" = "unblinded_mom",
          "Blinded MoM" = "blinded_mom"
        ),
        selected = "unblinded_ml"
      ),
      shiny::selectInput("adapt_analysis", "SSR adaptation look", choices = c("2" = 2), selected = 2),
      shiny::numericInput("max_n_multiplier", "Maximum N multiplier", value = 1.5, min = 1, step = 0.1),
      shiny::numericInput("min_if_futility", "Minimum IF before futility", value = 0.3, min = 0, max = 1, step = 0.05),
      shiny::numericInput("max_enrollment_frac_for_adapt", "Max enrollment fraction for SSR", value = 1, min = 0, max = 1, step = 0.05),
      shiny::numericInput("min_months_to_close_for_adapt", "Min time to close for SSR", value = 2, min = 0, step = 0.5),
      shiny::numericInput("analysis_lag_months", "Enrollment lag after futility", value = 2, min = 0, step = 0.5),
      shiny::numericInput("event_gap_days", "Event gap (days)", value = 20, min = 0, step = 1),
      shiny::numericInput("seed", "Random seed", value = 1234, min = 1, step = 1),
      shiny::actionButton("run", "Run simulation", class = "btn-primary")
    ),
    shiny::mainPanel(
      shiny::tabsetPanel(
        shiny::tabPanel(
          "Overview",
          shiny::br(),
          shiny::verbatimTextOutput("status"),
          shiny::tableOutput("design_table"),
          shiny::br(),
          shiny::selectInput("metric", "Plot metric", choices = metric_choices, selected = "mean_n"),
          shiny::plotOutput("summary_plot", height = "320px")
        ),
        shiny::tabPanel(
          "Trial summary",
          shiny::br(),
          shiny::tableOutput("trial_summary")
        ),
        shiny::tabPanel(
          "Analysis summary",
          shiny::br(),
          shiny::tableOutput("analysis_summary")
        ),
        shiny::tabPanel(
          "Example trials",
          shiny::br(),
          shiny::tableOutput("example_rows")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  analysis_times_input <- shiny::reactive({
    tryCatch(parse_numeric_csv(input$analysis_times), error = function(e) numeric(0))
  })

  shiny::observeEvent(analysis_times_input(), {
    times <- analysis_times_input()
    if (length(times) >= 2L) {
      interim_choices <- seq_len(length(times) - 1L)
      selected <- min(2L, max(interim_choices))
      shiny::updateSelectInput(
        session,
        "adapt_analysis",
        choices = stats::setNames(interim_choices, paste("Analysis", interim_choices)),
        selected = selected
      )
    }
  }, ignoreInit = FALSE)

  run_results <- shiny::eventReactive(input$run, {
    analysis_times <- parse_numeric_csv(input$analysis_times)
    if (length(analysis_times) < 2L) {
      stop("Provide at least one interim and one final analysis time.", call. = FALSE)
    }
    if (is.unsorted(analysis_times, strictly = TRUE)) {
      stop("Analysis times must be strictly increasing.", call. = FALSE)
    }

    trial_duration <- max(analysis_times)
    if (abs(tail(analysis_times, 1) - trial_duration) > sqrt(.Machine$double.eps)) {
      stop("The final analysis time could not be determined.", call. = FALSE)
    }

    if (length(input$strategies) == 0L) {
      stop("Select at least one SSR strategy.", call. = FALSE)
    }

    event_gap <- event_gap_in_unit(input$event_gap_days, input$time_unit)

    shiny::withProgress(message = "Running SSR prototype", value = 0, {
      shiny::incProgress(0.15, detail = "Building fixed design")
      fixed_design <- gsDesignNB::sample_size_nbinom(
        lambda1 = input$lambda1_plan,
        lambda2 = input$lambda1_plan * input$rr_plan,
        dispersion = input$k_plan,
        power = input$power_plan,
        alpha = input$alpha_plan,
        accrual_rate = input$accrual_rate_plan,
        accrual_duration = input$accrual_duration,
        trial_duration = trial_duration,
        max_followup = input$max_followup,
        event_gap = event_gap
      )

      shiny::incProgress(0.3, detail = "Building group sequential design")
      gs_design <- gsDesignNB::gsNBCalendar(
        fixed_design,
        k = length(analysis_times),
        test.type = 4,
        alpha = input$alpha_plan,
        sfu = gsDesignNB::sfHSD,
        sfupar = input$upper_spending,
        sfl = gsDesignNB::sfHSD,
        sflpar = input$lower_spending,
        analysis_times = analysis_times
      ) |>
        gsDesignNB::toInteger()

      shiny::incProgress(0.6, detail = "Running simulation")
      n_max <- ceiling(input$max_n_multiplier * gs_design$n_total[gs_design$k])
      sim_res <- gsDesignNB::sim_ssr_nbinom(
        n_sims = as.integer(input$n_sims),
        enroll_rate = data.frame(rate = input$accrual_true, duration = n_max / input$accrual_true),
        fail_rate = data.frame(
          treatment = c("Control", "Experimental"),
          rate = c(input$lambda1_true, input$lambda1_true * input$rr_true),
          dispersion = input$k_true
        ),
        dropout_rate = data.frame(
          treatment = c("Control", "Experimental"),
          rate = c(0, 0),
          duration = c(100, 100)
        ),
        max_followup = input$max_followup,
        design = gs_design,
        n_max = n_max,
        strategies = input$strategies,
        adapt_analysis = as.integer(input$adapt_analysis),
        min_if_futility = input$min_if_futility,
        max_enrollment_frac_for_adapt = input$max_enrollment_frac_for_adapt,
        min_months_to_close_for_adapt = input$min_months_to_close_for_adapt,
        analysis_lag_months = input$analysis_lag_months,
        event_gap = event_gap,
        bound_info = input$bound_info,
        seed = as.integer(input$seed)
      )

      shiny::incProgress(0.85, detail = "Summarizing")
      summary_res <- gsDesignNB::summarize_ssr_sim(sim_res, by = "strategy")

      list(
        fixed_design = fixed_design,
        gs_design = gs_design,
        sim_res = sim_res,
        summary = summary_res,
        analysis_times = analysis_times,
        event_gap = event_gap,
        n_max = n_max,
        trial_duration = trial_duration
      )
    })
  })

  output$status <- shiny::renderText({
    res <- run_results()
    shiny::req(res)
    paste(
      "Planned fixed N:", round(res$fixed_design$n_total, 1),
      "| Planned GS N:", ceiling(res$gs_design$n_total[res$gs_design$k]),
      "| Target information:", round(res$gs_design$n.I[res$gs_design$k], 2),
      "| Max simulated N:", res$n_max,
      "| Event gap:", round(res$event_gap, 4), input$time_unit
    )
  })

  output$design_table <- shiny::renderTable({
    res <- run_results()
    shiny::req(res)
    round_numeric_df(data.frame(
      Metric = c(
        "Planned control rate",
        "Planned rate ratio",
        "Planned dispersion",
        "Accrual duration",
        "Maximum follow-up",
        "Analysis times",
        "Fixed design N",
        "Group sequential N"
      ),
      Value = c(
        input$lambda1_plan,
        input$rr_plan,
        input$k_plan,
        input$accrual_duration,
        input$max_followup,
        paste(res$analysis_times, collapse = ", "),
        res$fixed_design$n_total,
        ceiling(res$gs_design$n_total[res$gs_design$k])
      ),
      stringsAsFactors = FALSE
    ), digits = 3)
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  output$summary_plot <- shiny::renderPlot({
    res <- run_results()
    shiny::req(res)
    sm <- res$summary$trial_summary
    metric <- input$metric
    shiny::req(metric %in% names(sm))

    y <- sm[[metric]]
    names(y) <- sm$strategy
    ylab <- names(metric_choices)[match(metric, metric_choices)]

    if (metric %in% c("rejection_rate", "futility_rate")) {
      y <- 100 * y
      ylab <- paste0(ylab, " (%)")
    }

    graphics::barplot(
      height = y,
      names.arg = names(y),
      col = c("#6baed6", "#3182bd", "#08519c")[seq_along(y)],
      ylab = ylab,
      las = 2,
      main = "Strategy comparison"
    )
  })

  output$trial_summary <- shiny::renderTable({
    res <- run_results()
    shiny::req(res)
    sm <- res$summary$trial_summary
    keep <- intersect(
      c(
        "strategy", "rejection_rate", "futility_rate", "mean_n", "pct_adapted",
        "expected_participants_with_events", "expected_events_observed",
        "mean_if_ia1", "mean_if_ia2", "mean_if_final",
        "mean_ia1_time", "mean_ia2_time", "mean_final_time"
      ),
      names(sm)
    )
    round_numeric_df(sm[, keep, drop = FALSE], digits = 3)
  }, striped = TRUE, bordered = TRUE, spacing = "xs")

  output$analysis_summary <- shiny::renderTable({
    res <- run_results()
    shiny::req(res)
    round_numeric_df(res$summary$analysis_summary, digits = 3)
  }, striped = TRUE, bordered = TRUE, spacing = "xs")

  output$example_rows <- shiny::renderTable({
    res <- run_results()
    shiny::req(res)
    trial_rows <- res$sim_res$trial_results
    keep <- intersect(
      c(
        "sim", "strategy", "reject_stage", "futility_stage", "n_adapted",
        "participants_with_events_stop", "events_observed_stop",
        "if_ia1", "if_ia2", "if_final"
      ),
      names(trial_rows)
    )
    round_numeric_df(utils::head(trial_rows[, keep, drop = FALSE], 8), digits = 3)
  }, striped = TRUE, bordered = TRUE, spacing = "xs")
}

app <- shiny::shinyApp(ui = ui, server = server)
app
