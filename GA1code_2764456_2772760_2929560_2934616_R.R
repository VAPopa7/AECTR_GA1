library(tidyverse)  
library(dplyr) 
library(lubridate) 
library(moments) 
library(ggplot2)
library(gridExtra)


#Question 2 
alpha <- 0.4
deltas <- c(0.3, 0.1, 0, -0.3)
gammas <- c(0.01, 0.1, 1)
x <- seq(-3, 3, length.out = 500)

# Create a data frame to hold all combinations
plot_data <- expand.grid(x = x, delta = deltas, gamma = gammas)

# Calculate the News-Impact Curve values
plot_data <- plot_data %>%
  mutate(
    nic = (alpha + delta * tanh(-gamma * x)) * x^2,
    delta_label = paste0("δ = ", delta),
    gamma_label = paste0("γ = ", gamma)
  )

# Generate the 2x2 grid plot
nic <- ggplot(plot_data, aes(x = x, y = nic, color = as.factor(gamma))) +
  geom_line(linewidth = 1) +
  facet_wrap(~ delta_label, scales = "free_y", ncol = 2) +
  labs(
    title = "News-Impact Curves for GARCH-M-L Model",
    subtitle = expression(paste("Holding ", sigma[t-1]^2, " = 1, ", mu, " = 0, ", lambda, " = 0, ", alpha, " = 0.4")),
    x = expression(x[t-1] ~ "(Shock)"),
    y = expression(sigma[t]^2 ~ "(Conditional Volatility)"),
    color = expression(gamma)
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave("nic.png", nic, width = 6.5, height = 5.2, dpi =300)

#Question 3 
tic_data<-read_csv("ticker_data.csv", show_col_types = FALSE) 
tic_data <- tic_data |> 
  mutate(RET = RET*100)|>
  arrange(date) |> 
  filter(!is.na(RET)) 

summary = tic_data |> 
  group_by(TICKER) |> 
  summarise( 
    Length = n(),
    Mean = mean(RET), 
    Med = median(RET),  
    SD = sd(RET), 
    Skew = skewness(RET), 
    Excess_Kurt = kurtosis(RET) - 3,  
    Min = min(RET), 
    Max = max(RET)
  )

plot <-ggplot(data = tic_data, aes(x=date, y=RET)) + 
  geom_line(linewidth = 0.5) + 
  facet_wrap(~ TICKER, ncol = 2, scales = "free")+ 
  scale_x_date(date_breaks = "1 years", date_labels = "%Y") +
  labs( 
    
    x = NULL, 
    y = "% per day"
  ) +
  theme_bw() + 
  theme( 
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold") 
  )
ggsave("plot_ret.png", plot, width = 10, height = 8, dpi =300)

#Question 4

pfe_ret  <- tic_data[tic_data$TICKER=="PFE", ][1:2500, "RET"]
jnj_ret  <- tic_data[tic_data$TICKER=="JNJ", ][1:2500, "RET"]
mrk_ret  <- tic_data[tic_data$TICKER=="MRK", ][1:2500, "RET"]

pfe_full <- tic_data[tic_data$TICKER=="PFE", ]["RET"]
jnj_full <- tic_data[tic_data$TICKER=="JNJ", ]["RET"]
mrk_full <- tic_data[tic_data$TICKER=="MRK", ]["RET"]

asset_rets <- list(PFE = pfe_ret, JNJ = jnj_ret, MRK = mrk_ret)
asset_full <- list(PFE = pfe_full, JNJ = jnj_full, MRK = mrk_full)

vol_garch <- function(par, data){
  # Define parameters
  mu <- par[1]
  omega <- par[2]
  alpha <- par[3]
  beta <- par[4]
  
  # Implementing statistical model
  data <- data[["RET"]]
  n <- length(data)
  sig2 <- rep(0,n)
  x_sub <- data[1:50]
  xbar <- mean(x_sub)
  sig2[1] <- mean((x_sub - xbar)^2)
  
  # Recursive variance loop
  for (t in 2:n){
    sig2[t] <- omega + alpha * ((data[t-1] - mu) / sqrt(sig2[t-1]))^2 + beta * sig2[t-1]
  }
  
  return(sqrt(sig2))
}


# Function for computing the loglikelihood for given GARCH model
ll_garch <- function(par, data){
  # Define parameters
  mu <- par[1]
  omega <- par[2]
  alpha <- par[3]
  beta <- par[4]
  v <- par[5]
  
  # Check constraints
  if (v <= 2) {
    return(1e10)
  }
  
  # Implementing statistical model
  data <- data[["RET"]]
  n <- length(data)
  sig2 <- rep(0,n)
  x_sub <- data[1:50]
  xbar <- mean(x_sub)
  sig2[1] <- mean((x_sub - xbar)^2)
  
  # Recursive variance loop
  for (t in 2:n){
    sig2[t] <- omega + alpha * ((data[t-1] - mu) / sqrt(sig2[t-1]))^2 + beta * sig2[t-1]
  }
  if (any(!is.finite(sig2)) || any(sig2 <= 0)) return(1e10)
  
  # Compute loglikelihood
  ll <- n * (lgamma((v + 1) / 2) - lgamma(v / 2) - 0.5 * log(v * pi)) - 0.5 * sum(log(sig2)) - ((v + 1) / 2) * sum(log(1 + (data - mu)^2 / (v * sig2)))
  
  return(-ll)
}

# GARCH-M Model

# Function for computing the volatilities for given GARCH-M model
vol_garch_m <- function(par, data){
  # Define parameters
  mu <- par[1]
  lambda <- par[2]
  omega <- par[3]
  alpha <- par[4]
  beta <- par[5]
  
  # Implementing statistical model
  data <- data[["RET"]]
  n <- length(data)
  sig2 <- rep(0,n)
  x_sub <- data[1:50]
  xbar <- mean(x_sub)
  sig2[1] <- mean((x_sub - xbar)^2)
  
  # Recursive variance loop
  for (t in 2:n){
    sig2[t] <- omega + alpha * ((data[t-1] - mu - lambda * sig2[t-1]) / sqrt(sig2[t-1]))^2 + beta * sig2[t-1]
  }
  
  return(sqrt(sig2))
}


# Function for computing the loglikelihood for given GARCH-M model
ll_garch_m <- function(par, data){
  # Define parameters
  mu <- par[1]
  lambda <- par[2]
  omega <- par[3]
  alpha <- par[4]
  beta <- par[5]
  v <- par[6]
  
  # Check constraints
  if (v <= 2) {
    return(1e10)
  }
  
  # Implementing statistical model
  data <- data[["RET"]]
  n <- length(data)
  sig2 <- rep(0,n)
  x_sub <- data[1:50]
  xbar <- mean(x_sub)
  sig2[1] <- mean((x_sub - xbar)^2)
  
  # Recursive variance loop
  for (t in 2:n){
    sig2[t] <- omega + alpha * ((data[t-1] - mu - lambda * sig2[t-1]) / sqrt(sig2[t-1]))^2 + beta * sig2[t-1]
  }
  if (any(!is.finite(sig2)) || any(sig2 <= 0)) return(1e10)
  
  # Compute loglikelihood
  ll <- n * (lgamma((v + 1) / 2) - lgamma(v / 2) - 0.5 * log(v * pi)) - 0.5 * sum(log(sig2)) - ((v + 1) / 2) * sum(log(1 + (data - mu - lambda * sig2)^2 / (v * sig2)))
  return(-ll)
}

# GARCH-M-L Model

# Function for computing the loglikelihood for given GARCH-M-L model
vol_garch_ml <- function(par, data){
  # Define parameters
  mu <- par[1]
  lambda <- par[2]
  omega <- par[3]
  alpha <- par[4]
  beta <- par[5]
  delta <- par[6]
  gamma <- par[7]
  
  # Implementing statistical model
  data <- data[["RET"]]
  n <- length(data)
  sig2 <- rep(0,n)
  x_sub <- data[1:50]
  xbar <- mean(x_sub)
  sig2[1] <- mean((x_sub - xbar)^2)
  
  # Recursive variance loop
  for (t in 2:n){
    sig2[t] <- omega + (alpha + delta * tanh(-gamma * data[t-1])) * ((data[t-1] - mu - lambda * sig2[t-1]) / sqrt(sig2[t-1]))^2 + beta * sig2[t-1]
  }
  
  return(sqrt(sig2))
}


# Function for computing the loglikelihood for given GARCH-M-L model
ll_garch_ml <- function(par, data){
  # Define parameters
  mu <- par[1]
  lambda <- par[2]
  omega <- par[3]
  alpha <- par[4]
  beta <- par[5]
  delta <- par[6]
  gamma <- par[7]
  v <- par[8]
  
  # Check constraints
  if (v <= 2) {
    return(1e10)
  }
  
  # Implementing statistical model
  data <- data[["RET"]]
  n <- length(data)
  sig2 <- rep(0,n)
  x_sub <- data[1:50]
  xbar <- mean(x_sub)
  sig2[1] <- mean((x_sub - xbar)^2)
  
  # Recursive variance loop
  for (t in 2:n){
    sig2[t] <- omega + (alpha + delta * tanh(-gamma * data[t-1])) * ((data[t-1] - mu - lambda * sig2[t-1]) / sqrt(sig2[t-1]))^2 + beta * sig2[t-1]
  }
  if (any(!is.finite(sig2)) || any(sig2 <= 0)) return(1e10)
  
  # Compute penalized loglikelihood
  ll <- n * (lgamma((v + 1) / 2) - lgamma(v / 2) - 0.5 * log(v * pi)) - 0.5 * sum(log(sig2)) - ((v + 1) / 2) * sum(log(1 + (data - mu - lambda * sig2)^2 / (v * sig2)))
  ll_pen <- ll - 0.001 * gamma^2
  
  return(-ll_pen)
}

# Function for computing AIC and BIC
compute_AIC_BIC <- function(logL, n, k){
  aic = 2 * k - 2 * logL 
  bic =  k * log(n) - 2 * logL
  return(c(aic, bic))
}

# AIC, BIC

# Function for computing AIC and BIC
compute_AIC_BIC <- function(logL, n, k){
  aic = 2 * k - 2 * logL 
  bic =  k * log(n) - 2 * logL
  return(c(aic, bic))
}

# NIC

# Set up shock grid
x <- seq(-3, 3, length.out = 500)

# Function for computing NIC for GARCH
nic_garch <- function(est_par, data, x) {
  est_omega <- est_par[[2]]
  est_alpha <- est_par[[3]]
  est_beta <- est_par[[4]]
  sig2_bar <- var(data[["RET"]])
  
  nic <- est_omega + est_alpha * (x^2) + est_beta * sig2_bar
  return(nic)
}

# Function for computing NIC for GARCH-M
nic_garch_m <- function(est_par, data, x) {
  est_omega <- est_par[[3]]
  est_alpha <- est_par[[4]]
  est_beta <- est_par[[5]]
  sig2_bar <- var(data[["RET"]])
  
  nic <- est_omega + est_alpha * (x^2) + est_beta * sig2_bar
  return(nic)
}

# Function for computing NIC for GARCH-M-L
nic_garch_ml <- function(est_par, data, x) {
  est_mu <- est_par[[1]]
  est_lambda <- est_par[[2]]
  est_omega <- est_par[[3]]
  est_alpha <- est_par[[4]]
  est_beta <- est_par[[5]]
  est_delta <- est_par[[6]]
  est_gamma <- est_par[[7]]
  sig2_bar <- var(data[["RET"]])
  
  x_lag <- est_mu + est_lambda * sig2_bar + sqrt(sig2_bar) * x
  
  nic <- est_omega + (est_alpha + est_delta * tanh(-est_gamma * x_lag)) * (x^2) + est_beta * sig2_bar
  return(nic)
}

# Results

# Estimate parameters, AIC and BIC for each ticker using initial parameters
ticker_vols <- c()
ticker_nics <- c()
for(ticker in names(asset_rets)){
  ini_s2 <- var(asset_rets[[ticker]][["RET"]])
  
  # GARCH model
  print("GARCH | Ticker, Estimated Parameters: mu, omega, alpha, beta, v, LogL, convergence, AIC, BIC")
  garch_par <- c(mu = 0, omega = ini_s2/50, alpha = 0.05, beta = 0.9, v = 10) # Define initial parameters
  est_garch <- optim(par = garch_par, fn = function(p) ll_garch(p, asset_rets[[ticker]]), method = "BFGS") # Estimate model 
  aic_bic_garch <- compute_AIC_BIC(-est_garch$value, nrow(asset_rets[[ticker]]), length(garch_par)) # Compute AIC and BIC
  cat(ticker, round(est_garch$par, 3), round(-est_garch$value, 0), est_garch$convergence, round(aic_bic_garch[1], 0), round(aic_bic_garch[2], 0), "\n")
  ticker_nics[[ticker]][["GARCH"]] <- nic_garch(est_garch$par, asset_full[[ticker]], x) # Compute NIC
  
  # GARCH-M model
  print("GARCH-M | Ticker, Estimated Parameters: mu, lambda, omega, alpha, beta, v, LogL, convergence, AIC, BIC")
  garch_m_par <- c(mu = 0, lambda = 0, omega = ini_s2/50, alpha = 0.05, beta = 0.9, v = 10)
  est_garch_m <- optim(par = garch_m_par, fn = function(p) ll_garch_m(p, asset_rets[[ticker]]), method = "BFGS")
  aic_bic_garch_m <- compute_AIC_BIC(-est_garch_m$value, nrow(asset_rets[[ticker]]), length(garch_m_par))
  cat(ticker, round(est_garch_m$par, 3), round(-est_garch_m$value, 0), est_garch_m$convergence, round(aic_bic_garch_m[1], 0), round(aic_bic_garch_m[2], 0), "\n") 
  ticker_nics[[ticker]][["GARCH-M"]] <- nic_garch_m(est_garch_m$par, asset_full[[ticker]], x)
  
  # GARCH-M-L model
  print("GARCH-M-L | Ticker, Estimated Parameters: mu, lambda, omega, alpha, beta, delta, gamma, v, LogL, convergence, AIC, BIC")
  garch_ml_par <- c(mu = 0, lambda = 0, omega = ini_s2/50, alpha = 0.05, beta = 0.9, delta = 0.01, gamma = 0.01, v = 10)
  sig_garch_ml <- vol_garch_ml(garch_ml_par, asset_full[[ticker]])
  ticker_vols[[ticker]] <- sig_garch_ml
  est_garch_ml <- optim(par = garch_ml_par, fn = function(p) ll_garch_ml(p, asset_rets[[ticker]]), method = "BFGS")
  aic_bic_garch_ml <- compute_AIC_BIC(-est_garch_ml$value, nrow(asset_rets[[ticker]]), length(garch_ml_par))
  cat(ticker, round(est_garch_ml$par, 3), round(-est_garch_ml$value, 0), est_garch_ml$convergence, round(aic_bic_garch_ml[1], 0), round(aic_bic_garch_ml[2], 0), "\n")
  ticker_nics[[ticker]][["GARCH-M-L"]] <- nic_garch_ml(est_garch_ml$par, asset_full[[ticker]], x)
  
}

#Question 5

plot_dates <- tic_data[tic_data$TICKER=="AAPL", ][["date"]]
oos_line <- tic_data[tic_data$TICKER=="AAPL", ][[2500, "date"]]

combined_plots <- list()
for(ticker in names(ticker_vols)){
  # Plot the NIC for each ticker, for each model
  df_nic <- data.frame(
    x = x,
    GARCH = ticker_nics[[ticker]][["GARCH"]],
    GARCH_M = ticker_nics[[ticker]][["GARCH-M"]],
    GARCH_ML = ticker_nics[[ticker]][["GARCH-M-L"]]
  )
  
  plot_nic <- ggplot(df_nic, aes(x = x)) +
    geom_line(aes(y = GARCH, color = "GARCH"), linetype = "dashed") +
    geom_line(aes(y = GARCH_M, color = "GARCH-M"), linetype = "dotdash") +
    geom_line(aes(y = GARCH_ML, color = "GARCH-M-L"), linetype = "solid") +
    scale_color_manual(values = c("GARCH" = "orange", "GARCH-M" = "green", "GARCH-M-L" = "blue")) +
    labs(title = paste("News-Impact Curves for", ticker), x = expression(x[t-1] ~ "(Shock)"), y = expression(sigma[t]^2 ~ "(TV Volatility)"), color = "Model") +
    theme_bw() +
    theme(legend.position = "bottom")
  
  # Plot the volatility for each ticker for the GARCH-M-L model
  plot_vol <- ggplot() +
    geom_line(aes(x = plot_dates, y = ticker_vols[[ticker]], color = "GARCH-M-L")) +
    geom_line(aes(x = plot_dates, y = asset_full[[ticker]][["RET"]], color = "Returns"), alpha = 0.4) +
    geom_vline(xintercept = oos_line, linetype = "dashed", color = "black", linewidth = 1) +
    scale_color_manual(values = c("Returns" = "gray", "GARCH-M-L" = "red")) +
    labs(title = paste("Filtered volatilities for", ticker), x = "", y = "", color = "Series") +
    theme_bw() +
    theme(legend.position = "bottom")
  
  # Store both plots in the list sequentially
  combined_plots[[paste0(ticker, "_vol")]] <- plot_nic
  combined_plots[[paste0(ticker, "_nic")]] <- plot_vol
}

grid <- grid.arrange(grobs = combined_plots, ncol = 2, nrow = 3)
dir.create("output", showWarnings = FALSE)
ggsave("output/q5.png", grid, width = 10, height = 8, dpi =300)

#Question 6

appl_ret <- tic_data$RET[tic_data$TICKER == "AAPL"]

# Function for computing the filtered variances for GARCH-M-L model
filter_sig2 <- function(par, x){
  # Define parameters
  mu <- par[1]
  lambda <- par[2]
  omega <- par[3]
  alpha <- par[4]
  beta <- par[5]
  delta <- par[6]
  gamma <- par[7]
  
  # Implementing statistical model
  n <- length(x)
  sig2 <- rep(0,n)
  x_sub <- x[1:50]
  xbar <- mean(x_sub)
  sig2[1] <- mean((x_sub - xbar)^2)
  
  # Recursive variance loop
  for (t in 2:n){
    sig2[t] <- omega + (alpha + delta * tanh(-gamma * x[t-1])) * ((x[t-1] - mu - lambda * sig2[t-1]) / sqrt(sig2[t-1]))^2 + beta * sig2[t-1]
  }
  
  return(sig2)
}

# Function for computing simulation based VaR at t0 for horizons 1, 5, 20 and levels 1%, 5%, 10%
simulate_var <- function(par, x, t0 = 2518, S = 10000, H = 20, horizons = c(1, 5, 20), levels = c(0.01, 0.05, 0.10)){
  # Define parameters
  mu <- par[1]
  lambda <- par[2]
  omega <- par[3]
  alpha <- par[4]
  beta <- par[5]
  delta <- par[6]
  gamma <- par[7]
  v <- par[8]
  
  # Filtered variance at t0 (Jan 4, 2021), reproduces Table 2 for AAPL
  sig2 <- filter_sig2(par[1:7], x[1:t0])
  sig2_now <- rep(sig2[t0], S)
  
  # Simulate S paths of H daily returns
  x_sim <- matrix(NA, S, H)
  for (j in 1:H){
    eps <- rt(S, df = v)
    x_sim[, j] <- mu + lambda * sig2_now + sqrt(sig2_now) * eps
    sig2_now <- omega + (alpha + delta * tanh(-gamma * x_sim[, j])) * eps^2 + beta * sig2_now
  }
  
  # Compound returns per horizon and their quantiles
  var_mat <- matrix(NA, length(horizons), length(levels), dimnames = list(paste(horizons, "step"), paste0(levels * 100, "%")))
  for (i in seq_along(horizons)){
    h <- horizons[i]
    comp_ret <- 100 * (apply(1 + x_sim[, 1:h, drop = FALSE] / 100, 1, prod) - 1)
    var_mat[i, ] <- quantile(comp_ret, levels)
  }
  
  return(var_mat)
}

# Test using Table 1 parameters for APPL, compare with Table 2
set.seed(42)
test_par <- list(
  M1 = c(mu = 0.154, lambda = 0, omega = 0.038, alpha = 0.090, beta = 0.873, delta = 0, gamma = 0, v = 4.146),
  M2 = c(mu = 0.072, lambda = 0.061, omega = 0.037, alpha = 0.089, beta = 0.875, delta = 0, gamma = 0, v = 4.138),
  M3 = c(mu = 0.108, lambda = 0.022, omega = 0.012, alpha = 0.073, beta = 0.915, delta = 0.071, gamma = 0.439, v = 4.402)
)
for (m in names(test_par)){
  cat("TEST RUN APPL", m, "\n")
  print(round(simulate_var(test_par[[m]], appl_ret), 2))
}

# Estimated parameters from Q4 (first 2500 obs): mu, lambda, omega, alpha, beta, delta, gamma, v
est_par <- list(
  PFE = list(
    M1 = c(0.037, 0, -0.027, 0.037, 0.961, 0, 0, 4.654),
    M2 = c(-0.006, 0.087, -0.025, 0.036, 0.959, 0, 0, 4.682),
    M3 = c(-0.036, 0.110, -0.024, 0.042, 0.952, 0.029, 0.408, 4.885)),
  JNJ = list(
    M1 = c(0.072, 0, -0.002, 0.026, 0.924, 0, 0, 4.430),
    M2 = c(0.058, 0.031, -0.002, 0.025, 0.924, 0, 0, 4.397),
    M3 = c(0.033, 0.063, 0.002, 0.026, 0.916, 0.017, 1.101, 4.550)),
  MRK = list(
    M1 = c(0.066, 0, 0.014, 0.039, 0.906, 0, 0, 4.499),
    M2 = c(0.078, -0.015, 0.015, 0.040, 0.904, 0, 0, 4.503),
    M3 = c(-0.045, 0.116, 0.020, 0.041, 0.905, 0.039, 2.465, 4.723))
)

# Compute VaR for each ticker and model
set.seed(42)
print("Ticker, Model, VaR: 1 step (1%, 5%, 10%), 5 steps (1%, 5%, 10%), 20 steps (1%, 5%, 10%)")
var_table <- c()
for (ticker in names(est_par)){
  for (m in names(est_par[[ticker]])){
    var_mat <- simulate_var(est_par[[ticker]][[m]], asset_full[[ticker]][["RET"]])
    var_row <- round(as.vector(t(var_mat)), 2)
    cat(ticker, m, var_row, "\n")
    var_table <- rbind(var_table, c(ticker, m, var_row))
  }
}
colnames(var_table) <- c("Ticker", "Model", paste(rep(c("1 step", "5 step", "20 step"), each = 3), rep(c("1%", "5%", "10%"), 3)))
dir.create("output", showWarnings = FALSE)
write.csv(var_table, "output/q6_var.csv", row.names = FALSE)
