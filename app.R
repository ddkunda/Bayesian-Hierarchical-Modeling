# app.R
library(shiny)
library(StanHeaders)
library(ggplot2)
library(parallel)
library(rstan)  

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


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

# Compile the Stan model once globally
stan_model_obj <- stan_model(model_code = stan_code)

# my_model <- stan_model(file = "../Bayesian_Hierarchical_Modeling/MyStanPkg1/inst/stan/my_stan_model1.stan")
# 
# saveRDS(my_model, file = "my_model_compiled.rds")
# 
# # Adjust the filename as needed for your OS (".so" vs ".dll")
# dyn.load("my_model.so")  
# my_model <- readRDS("my_model_compiled.rds")

# --------------------------------------------------------------------------
# UI
# --------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Bayesian Hierarchical Model (Drug X vs Placebo)"),
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
      br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(),
      span(h5(textOutput("uicmt14")), style="color:red")
    )
  )
)

# --------------------------------------------------------------------------
# SERVER
# --------------------------------------------------------------------------
server <- function(input, output, session) {
  
  # Load your pre-compiled model from RDS
  
  runBayes <- eventReactive(input$runAnalysis, {
    
    # Convert text input into numeric vectors
    drug_vals <- suppressWarnings(as.numeric(unlist(strsplit(input$drug_data, "[,\\s]+"))))
    placebo_vals <- suppressWarnings(as.numeric(unlist(strsplit(input$placebo_data, "[,\\s]+"))))
    
    # Remove NA values (caused by bad formatting)
    drug_vals <- na.omit(drug_vals)
    placebo_vals <- na.omit(placebo_vals)
    
    # Ensure at least 2 data points exist in each group
    if (length(drug_vals) < 2 || length(placebo_vals) < 2) {
      showNotification("Error: Please enter at least two values for each group.", type = "error")
      return(NULL)
    }
    
    # Build the data list for Stan – make sure these names match your Stan model data block!
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
    
    # Run Stan sampling using the precompiled model from MyStanPkg1.
    # This calls your package function that internally uses the compiled Stan model.
    # 3. Fit with Stan
    fit <- sampling(
      object = stan_model_obj,
      data   = stan_data_list,
      iter   = 1000,
      warmup = 500,
      chains = 4,
      seed   = 42,
      control = list(adapt_delta = 0.95, max_treedepth = 15)
    )
    
    return(fit)
  })
  
  # Output summary
  output$summaryOutput <- renderPrint({
    req(runBayes())  # Ensure results exist before running
    
    fit <- runBayes()
    posterior_samples <- extract(fit)
    drug_minus_placebo <- posterior_samples$mu_drug - posterior_samples$mu_placebo
    threshold <- input$threshold
    
    prob_effective <- mean(drug_minus_placebo > threshold)
    
    cat("Posterior Summary:\n")
    print(summary(fit, probs = c(0.025, 0.5, 0.975))$summary)
    
    cat("\n-----------------------------------------------------------\n")
    cat(sprintf("Probability(Drug X - Placebo > %0.2f) = %1.2f%%\n",
                threshold, 100 * prob_effective))
  })
  
  # Posterior Plot
  output$posteriorPlot <- renderPlot({
    req(runBayes())
    
    fit <- runBayes()
    posterior_samples <- extract(fit)
    
    df_all <- data.frame(
      value = c(
        posterior_samples$mu_drug, 
        posterior_samples$mu_placebo, 
        posterior_samples$mu_drug - posterior_samples$mu_placebo
      ),
      group = rep(
        c("mu_drug", "mu_placebo", "mu_drug - mu_placebo"), 
        each = length(posterior_samples$mu_drug)
      )
    )
    
    ggplot(df_all, aes(x = value, fill = group)) +
      geom_density(alpha = 0.4) +
      facet_wrap(~ group, scales = "free") +
      theme_minimal() +
      xlab("Posterior Estimate") + ylab("Density") +
      ggtitle("Posterior Distributions")
  })
  
  output$uicmt14 <- renderText({
    "Any question can be sent to Christian Dide-Agossou, PhD: christian.dideagossou@gmail.com"
  })
  
}

shinyApp(ui = ui, server = server)
