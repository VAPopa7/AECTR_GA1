library(readr)

#load data
tic_data <- read_csv("input/ticker_data.csv", show_col_types = FALSE)
appl_ret <- tic_data$RET[tic_data$TICKER == "AAPL"] * 100
pfe_ret  <- tic_data$RET[tic_data$TICKER == "PFE"] * 100
jnj_ret  <- tic_data$RET[tic_data$TICKER == "JNJ"] * 100
mrk_ret  <- tic_data$RET[tic_data$TICKER == "MRK"] * 100

asset_rets <- list(PFE = pfe_ret, JNJ = jnj_ret, MRK = mrk_ret)

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
        var_mat <- simulate_var(est_par[[ticker]][[m]], asset_rets[[ticker]])
        var_row <- round(as.vector(t(var_mat)), 2)
        cat(ticker, m, var_row, "\n")
        var_table <- rbind(var_table, c(ticker, m, var_row))
    }
}
colnames(var_table) <- c("Ticker", "Model", paste(rep(c("1 step", "5 step", "20 step"), each = 3), rep(c("1%", "5%", "10%"), 3)))
write.csv(var_table, "output/q6_var.csv", row.names = FALSE)

