library(shiny)
library(cmdstanr)
library(ggplot2)
library(parallel)
library(posterior)
library(polished)
library(httr)


# --- Global Code (runs once at startup) ---

# Define the Stan model code
stan_code <- "
data {
  int<lower=0> N_drug;
  array[N_drug] real drug_x_data;

  int<lower=0> N_placebo;
  array[N_placebo] real placebo_data;

  real prior_mu_drug_mean;
  real<lower=0> prior_mu_drug_sd;
  real prior_mu_placebo_mean;
  real<lower=0> prior_mu_placebo_sd;
  real<lower=0> prior_sigma_scale;
}

parameters {
  real mu_drug;
  real mu_placebo;
  real<lower=0> sigma;
}

model {
  // Priors
  mu_drug ~ normal(prior_mu_drug_mean, prior_mu_drug_sd);
  mu_placebo ~ normal(prior_mu_placebo_mean, prior_mu_placebo_sd);
  sigma ~ cauchy(0, prior_sigma_scale);

  // Likelihood
  drug_x_data ~ normal(mu_drug, sigma);
  placebo_data ~ normal(mu_placebo, sigma);
}
"


#cmdstanr::install_cmdstan()
# cmdstanr::set_cmdstan_path("/Users/chris/.cmdstan/cmdstan-2.36.0")
# cmdstanr::cmdstan_version()

# Write the Stan model to a file and compile it
stan_file <- write_stan_file(stan_code)
cmdstan_model_obj <- cmdstan_model(stan_file)

# --- End Global Code ---

# UI definition
ui <- fluidPage(
  titlePanel("Bayesian Hierarchical Modeling: (Drug X vs Placebo)"),
  sidebarLayout(
    sidebarPanel(
      h4("1. Input Priors (e.g., percent reduction)"),
      numericInput("prior_mu_drug_mean", "Prior Mean for Drug X (%)", value = 30, step = 1),
      numericInput("prior_mu_drug_sd", "Prior SD for Drug X (%)", value = 10, step = 1),
      numericInput("prior_mu_placebo_mean", "Prior Mean for Placebo (%)", value = 0, step = 1),
      numericInput("prior_mu_placebo_sd", "Prior SD for Placebo (%)", value = 10, step = 1),
      numericInput("prior_sigma_scale", "Prior Scale for Sigma (Cauchy)", value = 10, step = 1),
      
      h4("2. Enter Observed Data"),
      helpText("Enter percent reduction values (in %), comma-separated. Example: 35, 28, 42, ..."),
      textAreaInput("drug_data", "Drug X Data (%)", 
                    value = "35, 30, 45, 25, 40, 39, 32, 50, 42, 27, 34, 36, 29, 37, 41, 38, 33, 46, 48, 31",
                    rows = 5),
      textAreaInput("placebo_data", "Placebo Data (%)",
                    value = "5, 15, -2, 10, 0, 8, 2, 4, 3, 7, -5, 1, 6, 9, 0, 11, -1, 12, 3, 1",
                    rows = 5),
      
      h4("3. Effectiveness Threshold"),
      numericInput("threshold", "Threshold for (mu_drug - mu_placebo) (%)", value = 10, step = 1),
      
      actionButton("runAnalysis", "Run Bayesian Analysis")
    ),
    
    mainPanel(
      h3("Results"),
      verbatimTextOutput("summaryOutput"),
      plotOutput("posteriorPlot"),
      br(), br(), br(), br(), br(), br(), br(), br(),
      
      h5(
        span("All questions can be sent to Christian Dide-Agossou, PhD:", style = "color:black; display: inline-block;"),
        span(htmlOutput("uicmt14"), style = "display: inline-block;")
      )
    )
  )
)

# SERVER definition
server <- function(input, output, session) {
  
  runBayes <- eventReactive(input$runAnalysis, {
    # Convert text input into numeric vectors
    drug_vals <- suppressWarnings(as.numeric(unlist(strsplit(input$drug_data, "[,\\s]+"))))
    placebo_vals <- suppressWarnings(as.numeric(unlist(strsplit(input$placebo_data, "[,\\s]+"))))
    
    # Remove NA values
    drug_vals <- na.omit(drug_vals)
    placebo_vals <- na.omit(placebo_vals)
    
    # Validate that we have at least 2 data points per group
    if (length(drug_vals) < 2 || length(placebo_vals) < 2) {
      showNotification("Error: Please enter at least two values for each group.", type = "error")
      return(NULL)
    }
    
    # Build the data list for Stan – ensure names match those in the Stan model
    stan_data_list <- list(
      N_drug = length(drug_vals),
      drug_x_data = drug_vals,
      N_placebo = length(placebo_vals),
      placebo_data = placebo_vals,
      prior_mu_drug_mean = input$prior_mu_drug_mean,
      prior_mu_drug_sd = input$prior_mu_drug_sd,
      prior_mu_placebo_mean = input$prior_mu_placebo_mean,
      prior_mu_placebo_sd = input$prior_mu_placebo_sd,
      prior_sigma_scale = input$prior_sigma_scale
    )
    
    # Run sampling using cmdstanr's model object from the global environment.
    # Suppress progress output by setting refresh = 0.
    fit <- cmdstan_model_obj$sample(
      data = stan_data_list,
      chains = 4,
      iter_warmup = 1000,
      iter_sampling = 2000,  # Total iterations = 3000 (1000 warmup + 2000 sampling)
      seed = 42,
      adapt_delta = 0.95,
      max_treedepth = 15,
      refresh = 0      # Disables progress output
    )
    
    return(fit)
  })
  
  # Output summary of the posterior samples
  output$summaryOutput <- renderPrint({
    fit <- runBayes()
    req(fit)
    
    # Convert draws to a tibble for easier manipulation using the posterior package
    library(posterior)
    draws <- as_draws_df(fit$draws())
    
    # Check that the required parameters are present
    if (!("mu_drug" %in% names(draws)) || !("mu_placebo" %in% names(draws))) {
      stop("Columns 'mu_drug' or 'mu_placebo' not found in the draws.")
    }
    
    # Compute the difference between mu_drug and mu_placebo
    drug_minus_placebo <- draws$mu_drug - draws$mu_placebo
    threshold <- input$threshold
    prob_effective <- mean(drug_minus_placebo > threshold)
    
    # Print summary for key parameters using cmdstanr's built-in summary
    cat("Posterior Summary:\n")
    print(fit$summary(c("mu_drug", "mu_placebo", "sigma")))
    
    cat("\n-----------------------------------------------------------\n")
    cat(sprintf("Probability(Drug X - Placebo > %0.2f) = %1.2f%%\n",
                threshold, 100 * prob_effective))
  })
  
  # Posterior Plot
  output$posteriorPlot <- renderPlot({
    fit <- runBayes()
    req(fit)
    
    # Convert draws to a tibble
    library(posterior)
    draws <- as_draws_df(fit$draws())
    
    # Check that required columns exist
    if (!("mu_drug" %in% names(draws)) || !("mu_placebo" %in% names(draws))) {
      stop("Columns 'mu_drug' or 'mu_placebo' not found in the draws.")
    }
    
    # Compute the difference between mu_drug and mu_placebo
    draws$mu_diff <- draws$mu_drug - draws$mu_placebo
    
    # Use the number of draws (rows) for plotting
    n_draws <- nrow(draws)
    
    # Build a combined data frame for plotting densities
    df_all <- data.frame(
      value = c(draws$mu_drug, draws$mu_placebo, draws$mu_diff),
      group = rep(c("mu_drug", "mu_placebo", "mu_drug - mu_placebo"), each = n_draws)
    )
    
    # Create the density plot using ggplot2
    ggplot(df_all, aes(x = value, fill = group)) +
      geom_density(alpha = 0.4) +
      facet_wrap(~ group, scales = "free") +
      theme_minimal() +
      xlab("Posterior Estimate") +
      ylab("Density") +
      ggtitle("Posterior Distributions")
  })
  
  output$uicmt14 <- renderUI({
    email <- "christian.dideagossou@gmail.com"
    link <- paste0("mailto:", email)
    tags$a(href = link, email)  # Use tags$a to create the hyperlink
  })
  
}

# Launch the app
shinyApp(ui = ui, server = server)

