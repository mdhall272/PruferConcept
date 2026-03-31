library(shiny)
library(tidyverse)
library(Rcpp)

sourceCpp("forest_rcpp.cpp")
source("forest_v2.5.R", local = TRUE)

# Override R functions with faster C++ versions
log.forest.count <- function(N, n, k1, k2, roots.known = TRUE) {
  log_forest_count_cpp(N, n, k1, k2, roots.known)
}

log.likelihood <- function(data, sequence, components, seq.type, cluster = FALSE) {
  log_likelihood_cpp(as.integer(data), sequence, as.integer(components), seq.type, cluster)
}

# UI
ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      body {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        min-height: 100vh;
      }
      .container-fluid {
        background: transparent;
      }
      h2 {
        color: white;
        text-shadow: 2px 2px 4px rgba(0,0,0,0.3);
      }
      .well {
        background: rgba(224, 242, 255, 0.95);
        border: none;
        border-radius: 15px;
        box-shadow: 0 8px 32px rgba(0,0,0,0.15);
        color: #333;
      }
      .well h4 {
        color: #667eea;
        border-bottom: 2px solid #e0e0e0;
        padding-bottom: 8px;
        margin-bottom: 15px;
      }
      .well hr {
        border-color: #e0e0e0;
      }
      .form-control {
        background-color: #f0f8ff;
        border: 1px solid #b8d4e8;
        border-radius: 8px;
      }
      .form-control:focus {
        border-color: #667eea;
        box-shadow: 0 0 8px rgba(102,126,234,0.3);
      }
      .selectize-input {
        background-color: #f0f8ff !important;
        border: 1px solid #b8d4e8 !important;
        border-radius: 8px !important;
      }
      .btn-primary {
        background: linear-gradient(45deg, #667eea, #764ba2);
        border: none;
        border-radius: 25px;
        font-weight: bold;
        color: white;
        padding: 12px 30px;
        box-shadow: 0 4px 15px rgba(102,126,234,0.4);
        transition: all 0.3s ease;
      }
      .btn-primary:hover {
        transform: translateY(-2px);
        box-shadow: 0 6px 20px rgba(102,126,234,0.6);
        background: linear-gradient(45deg, #764ba2, #667eea);
      }
      .shiny-output-error { color: #e74c3c; }
      pre {
        background-color: #e8f4fc;
        border: 1px solid #b8d4e8;
        border-radius: 10px;
        color: #333;
      }
      .shiny-plot-output {
        background: #e0f2ff;
        border-radius: 15px;
        padding: 25px;
        box-shadow: 0 4px 20px rgba(0,0,0,0.1);
        margin-bottom: 20px;
        overflow: hidden;
        display: flex;
        justify-content: center;
        align-items: center;
      }
      .shiny-plot-output img {
        border-radius: 10px;
        max-width: 100%;
        height: auto;
        display: block;
        margin: 0 auto;
      }
    "))
  ),

  titlePanel("Proof of Concept: Transmission Forest Inference Dashboard"),

  tags$p(style = "color: white; font-size: 16px; margin-bottom: 10px;",
         "For an infectious disease epidemic, infer the total size of a samplable infected population (N) and the number of independent lineage introductions (k) from observed cluster data."),

  tags$p(style = "color: white; margin-bottom: 20px;",
         tags$a(href = "https://github.com/mdhall272/PruferConcept/blob/main/Lighting_Up_Random_Trees.pdf",
                style = "color: #a0d8ef;",
                target = "_blank",
                "Read the accompanying paper")),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("Data"),
      tags$p(tags$small("Enter the clusters and singletons found, separated by commas. For example, if you observed two clusters of size 3, one of size 2, and three singletons, enter '3,3,2,1,1,1'. The order does not matter.")),
      textInput("cluster_sizes", "Cluster sizes (comma-separated)", value = "3, 3, 2, 1, 1, 1"),
      checkboxInput("cluster_mode", "Unresolved clusters (marginalize over subtree structures within clusters)", value = FALSE),

      hr(),
      actionButton("run", "Run Inference", class = "btn-primary btn-lg"),

      hr(),
      h4("Proposal Distribution"),
      tags$p(tags$small("Controls importance sampling proposal; tune to match expected posterior")),
      tags$p(tags$small("Number of unsampled individuals N-n ~ NegBinom(mean, size)")),
      numericInput("nb_mean", "Mean", value = 18, min = 1, step = 1),
      numericInput("nb_size", "Size (dispersion)", value = 0.5, min = 0.01, step = 0.1),
      tags$p(tags$small("Number of introductions k ~ 1 + Binomial(N-1, (mean-1)/(N-1))")),
      numericInput("comp_mean", "Mean", value = 3, min = 1, step = 0.5),

      conditionalPanel(
        condition = "!input.fix_N",
        hr(),
        h4("Prior for N - n (Unsampled Individuals)"),
        tags$p(tags$small("Number of infected individuals not sampled. N - n = 0 means everyone was sampled.")),
        selectInput("prior_N_dist", "Distribution",
                    choices = c("Poisson" = "poisson", "Negative Binomial" = "nbinom"),
                    selected = "poisson"),
        numericInput("prior_N_mean", "Mean of N - n", value = 18, min = 0, step = 1),
        conditionalPanel(
          condition = "input.prior_N_dist == 'nbinom'",
          numericInput("prior_N_size", "Size (dispersion)", value = 2, min = 0.01, step = 0.1)
        )
      ),

      conditionalPanel(
        condition = "!input.fix_k",
        hr(),
        h4("Prior for k - 1 (Extra Introductions)"),
        tags$p(tags$small("Number of introductions minus one. k - 1 = 0 means a single introduction.")),
        selectInput("prior_k_dist", "Distribution",
                    choices = c("Poisson" = "poisson", "Negative Binomial" = "nbinom"),
                    selected = "poisson"),
        numericInput("prior_k_mean", "Mean of k - 1", value = 2, min = 0, step = 0.5),
        conditionalPanel(
          condition = "input.prior_k_dist == 'nbinom'",
          numericInput("prior_k_size", "Size (dispersion)", value = 2, min = 0.01, step = 0.1)
        )
      ),

      hr(),
      h4("Inference Settings"),
      numericInput("n_samples", "Number of IS samples", value = 50000, min = 100, max = 1000000, step = 1000),

      checkboxInput("fix_N", "Fix infected population size (N)", value = FALSE),
      conditionalPanel(
        condition = "input.fix_N",
        numericInput("fixed_N", "Fixed N", value = 20, min = 1, step = 1)
      ),

      checkboxInput("fix_k", "Fix introductions (k)", value = FALSE),
      conditionalPanel(
        condition = "input.fix_k",
        numericInput("fixed_k", "Fixed k", value = 1, min = 1, step = 1)
      )
    ),

    mainPanel(
      width = 9,

      fluidRow(
        column(6,
               h4("Data Summary"),
               verbatimTextOutput("data_summary")
        ),
        column(6,
               h4("Diagnostics"),
               verbatimTextOutput("diagnostics")
        )
      ),

      hr(),

      fluidRow(
        conditionalPanel(
          condition = "!input.fix_N",
          column(6, plotOutput("posterior_size", height = "400px"))
        ),
        conditionalPanel(
          condition = "!input.fix_k",
          column(6, plotOutput("posterior_components", height = "400px"))
        )
      ),
      conditionalPanel(
        condition = "!input.fix_N && !input.fix_k",
        fluidRow(
          column(12, plotOutput("joint_posterior", height = "500px"))
        )
      )
    )
  )
)

# Server
server <- function(input, output, session) {

  parsed_data <- reactive({
    sizes_str <- input$cluster_sizes
    sizes <- as.integer(strsplit(gsub(" ", "", sizes_str), ",")[[1]])
    sizes <- sizes[!is.na(sizes) & sizes > 0]
    if (length(sizes) == 0) sizes <- c(1)
    sizes
  })

  results <- eventReactive(input$run, {
    data <- parsed_data()
    n <- length(data)

    withProgress(message = 'Running inference...', value = 0, {

      nb_mean <- input$nb_mean
      comp_mean <- input$comp_mean
      nb_size <- input$nb_size

      fix_N_val <- if (input$fix_N) as.integer(input$fixed_N - sum(data) + n) else 0L
      fix_k_val <- if (input$fix_k) as.integer(input$fixed_k) else 0L

      incProgress(0.1, detail = "Generating samples...")

      samples <- generate_samples_batch_cpp(
        Nreps = as.integer(input$n_samples),
        n_sampled = as.integer(n),
        nb_mean = nb_mean,
        nb_size = nb_size,
        components_mean = comp_mean,
        fix_N = fix_N_val,
        fix_k = fix_k_val
      )

      incProgress(0.3, detail = "Computing posterior weights...")

      params <- list(
        prior_N_dist = input$prior_N_dist,
        prior_N_mean = input$prior_N_mean,
        prior_N_size = if (input$prior_N_dist == "nbinom") input$prior_N_size else 1,
        prior_k_dist = input$prior_k_dist,
        prior_k_mean = input$prior_k_mean,
        prior_k_size = if (input$prior_k_dist == "nbinom") input$prior_k_size else 1,
        nb_mean = nb_mean,
        nb_size = nb_size,
        components_mean = comp_mean,
        fix_N = input$fix_N,
        fix_k = input$fix_k,
        cluster = input$cluster_mode
      )

      result <- compute_posterior_batch_cpp(as.integer(data), samples, params)

      incProgress(0.4, detail = "Normalizing weights...")

      df <- tibble(
        actual.size = as.integer(result$actual.size),
        components = samples$components,
        k1 = samples$k1,
        k2 = samples$k2,
        type = as.character(samples$type),
        likelihood = result$likelihood,
        prior = result$prior,
        samp.prob = result$samp.prob,
        log.weight = result$log.weight
      )

      log_weights <- df$log.weight - max(df$log.weight)
      weights <- exp(log_weights)
      weights <- weights / sum(weights)

      ESS <- 1 / sum(weights^2)
      post_mean_size <- sum(weights * df$actual.size)
      post_mean_comp <- sum(weights * df$components)

      size_dist <- df %>%
        mutate(weight = weights) %>%
        group_by(actual.size) %>%
        summarise(prob = sum(weight), .groups = "drop")

      comp_dist <- df %>%
        mutate(weight = weights) %>%
        group_by(components) %>%
        summarise(prob = sum(weight), .groups = "drop")

      # KL divergence
      min_size <- sum(data)
      kl_N <- NA
      kl_k <- NA

      if (!input$fix_N) {
        size_prior <- size_dist %>%
          mutate(prior_prob = if (input$prior_N_dist == "poisson") {
            dpois(actual.size - min_size, input$prior_N_mean)
          } else {
            dnbinom(actual.size - min_size, mu = input$prior_N_mean, size = input$prior_N_size)
          })
        valid <- size_prior %>% filter(prob > 0, prior_prob > 0)
        kl_N <- sum(valid$prob * log(valid$prob / valid$prior_prob))
      }

      if (!input$fix_k) {
        comp_prior <- comp_dist %>%
          mutate(prior_prob = if (input$prior_k_dist == "poisson") {
            dpois(components - 1, input$prior_k_mean)
          } else {
            dnbinom(components - 1, mu = input$prior_k_mean, size = input$prior_k_size)
          })
        valid <- comp_prior %>% filter(prob > 0, prior_prob > 0)
        kl_k <- sum(valid$prob * log(valid$prob / valid$prior_prob))
      }

      incProgress(0.2, detail = "Done!")

      list(
        kl_N = kl_N,
        kl_k = kl_k,
        result = df,
        weights = weights,
        ESS = ESS,
        post_mean_size = post_mean_size,
        post_mean_comp = post_mean_comp,
        size_dist = size_dist,
        comp_dist = comp_dist,
        data = data,
        fix_N = input$fix_N,
        fixed_N = input$fixed_N,
        fix_k = input$fix_k,
        fixed_k = input$fixed_k,
        n_samples_used = nrow(df)
      )
    })
  })

  output$data_summary <- renderPrint({
    data <- parsed_data()
    cat("Number of clusters (n):", length(data), "\n")
    cat("Cluster sizes:", paste(data, collapse = ", "), "\n")
    cat("Total observed individuals:", sum(data), "\n")
    cat("Cluster mode:", ifelse(input$cluster_mode, "Unresolved", "Fully resolved"), "\n")
  })

  output$diagnostics <- renderPrint({
    req(results())
    res <- results()

    if (res$fix_N || res$fix_k) {
      cat("Constraints:\n")
      if (res$fix_N) cat("  Infected population size (N) fixed at", res$fixed_N, "\n")
      if (res$fix_k) cat("  Introductions (k) fixed at", res$fixed_k, "\n")
      cat("\n")
    }

    cat("=== Importance Sampling Results ===\n\n")

    n_used <- res$n_samples_used
    cat(sprintf("Total samples drawn: %d\n", n_used))
    cat("Effective Sample Size (ESS):", round(res$ESS, 1),
        sprintf("(%.1f%%)\n", 100 * res$ESS / n_used))
    cat("\nPosterior means:\n")
    if (res$fix_N) {
      cat("  Infected population size (N):", res$fixed_N, "(fixed)\n")
    } else {
      cat("  Infected population size (N):", round(res$post_mean_size, 2), "\n")
    }
    if (res$fix_k) {
      cat("  Introductions (k):", res$fixed_k, "(fixed)\n")
    } else {
      cat("  Introductions (k):", round(res$post_mean_comp, 2), "\n")
    }

    cat("\nKL divergence (posterior || prior):\n")
    if (!res$fix_N) {
      cat(sprintf("  N: %.4f nats (%.4f bits)\n", res$kl_N, res$kl_N / log(2)))
    }
    if (!res$fix_k) {
      cat(sprintf("  k: %.4f nats (%.4f bits)\n", res$kl_k, res$kl_k / log(2)))
    }

    if (res$ESS < 100) {
      cat("\n Warning: ESS is low. Consider:\n")
      cat("  - Decreasing nb.size for heavier tails\n")
      cat("  - Adjusting nb.mean closer to posterior\n")
    }
  })

  output$posterior_size <- renderPlot({
    req(results())
    res <- results()

    ggplot(res$size_dist, aes(x = actual.size, y = prob)) +
      geom_col(fill = "#6A5ACD", alpha = 0.8) +
      geom_vline(xintercept = res$post_mean_size, color = "#FF6B6B", linetype = "dashed", linewidth = 1.2) +
      labs(
        title = "Posterior: Infected Population Size (N)",
        subtitle = sprintf("Posterior mean: %.1f", res$post_mean_size),
        x = "Infected population size (N)",
        y = "Probability"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(face = "bold", color = "#333"),
        plot.subtitle = element_text(color = "#666"),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_line(color = "#e0e0e0"),
        panel.grid.minor = element_line(color = "#f0f0f0"),
        axis.text = element_text(color = "#666"),
        axis.title = element_text(color = "#333")
      )
  })

  output$posterior_components <- renderPlot({
    req(results())
    res <- results()

    ggplot(res$comp_dist, aes(x = components, y = prob)) +
      geom_col(fill = "#20B2AA", alpha = 0.8) +
      geom_vline(xintercept = res$post_mean_comp, color = "#FF6B6B", linetype = "dashed", linewidth = 1.2) +
      labs(
        title = "Posterior: Lineage Introductions (k)",
        subtitle = sprintf("Posterior mean: %.1f", res$post_mean_comp),
        x = "Number of Introductions (k)",
        y = "Probability"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(face = "bold", color = "#333"),
        plot.subtitle = element_text(color = "#666"),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_line(color = "#e0e0e0"),
        panel.grid.minor = element_line(color = "#f0f0f0"),
        axis.text = element_text(color = "#666"),
        axis.title = element_text(color = "#333")
      )
  })

  output$joint_posterior <- renderPlot({
    req(results())
    res <- results()

    joint_dist <- res$result %>%
      mutate(weight = res$weights) %>%
      group_by(actual.size, components) %>%
      summarise(prob = sum(weight), .groups = "drop")

    ggplot(joint_dist, aes(x = actual.size, y = components, fill = prob)) +
      geom_tile() +
      scale_fill_viridis_c(name = "Probability", trans = "log10") +
      geom_point(aes(x = res$post_mean_size, y = res$post_mean_comp),
                 color = "white", size = 5, shape = 4, stroke = 2.5) +
      labs(
        title = "Joint Posterior Distribution",
        subtitle = sprintf("ESS: %.1f (%.1f%%)", res$ESS, 100 * res$ESS / res$n_samples_used),
        x = "Infected population size (N)",
        y = "Introductions (k)"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(face = "bold", color = "#333"),
        plot.subtitle = element_text(color = "#666"),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_line(color = "#e0e0e0"),
        panel.grid.minor = element_line(color = "#f0f0f0"),
        axis.text = element_text(color = "#666"),
        axis.title = element_text(color = "#333"),
        legend.background = element_rect(fill = "white", color = NA),
        legend.text = element_text(color = "#666"),
        legend.title = element_text(color = "#333")
      )
  })
}

shinyApp(ui = ui, server = server)
