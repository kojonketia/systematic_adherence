
rm(list = ls())

pacman::p_load(Rcpp, parallel, tidyverse, abc, dplyr, Hmisc)
sourceCpp(file.path("code", "Transfil.cpp"))

set.seed(123)
start.time <- Sys.time()

# observed baseline prevalence
prev_df <- read.csv(file.path("data", "cleaned", "prevalence.csv"))
prevalence_obs <- prev_df$worm_prevalence[prev_df$year == 2001]

# TRANSFIL model
simulate_baseline_prevalence <- function(k, vth) {
  pop_size <- 1000
  model <- new(Data, pop_size)
  model$set_parameter("k", k)
  model$set_parameter("vth", vth)
  model$set_parameter("importation_rate", 0.0003)
  
  model$burn_in()
  
  return(model$worm_prev)
}

# generate prior parameter values and simulations
simulate_results <- function(i) {
  k_prior <- runif(1, min = 0.001, max = 0.3)
  vth_prior <- runif(1, min = 30, max = 150)
  
  sim_data <- simulate_baseline_prevalence(k = k_prior, vth = vth_prior)
  data.frame(simulation = i, k = k_prior, vth = vth_prior, prev = sim_data)
}

# run simulations in parallel
N <- 1000000
n_cores <- max(1, parallel::detectCores() - 1)

output_file <- file.path("output", "simulated_baseline_prevalence.csv")
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

if (file.exists(output_file)) {
  message("Reading existing simulation results...")
  results <- read.csv(output_file)
} else {
  message(sprintf("Running %d simulations using %d cores...", N, n_cores) )
  
  cl <- makeCluster(n_cores)
  clusterEvalQ(cl, {library(Rcpp)
    Rcpp::sourceCpp(file.path("code", "Transfil.cpp"))
  })
  
  clusterExport(cl, 
                varlist = c("simulate_results", "simulate_baseline_prevalence"),
                envir = environment())
  
  results_list <- tryCatch(parLapply(cl, seq_len(N), simulate_results),
                           finally = stopCluster(cl))
  
  results <- bind_rows(results_list)
  
  # retain valid simulation results
  write.csv(results, output_file, row.names = FALSE)
}

# prepare inputs for ABC

parameter_matrix <- as.matrix(results[, c("k", "vth")])
summary_matrix <- matrix(results$prev, ncol = 1, dimnames = list(NULL, "worm_prevalence"))
target_stat <- c(as.numeric(prevalence_obs)) # observed summary statistic

# fit parameters using ABC rejection
# Retain the 1% of simulations closest to the
# observed baseline prevalence.
abc_tol <- 0.01
abc_fit <- abc::abc(target = target_stat,
                    param = parameter_matrix,
                    sumstat = summary_matrix,
                    tol = abc_tol,
                    method = "rejection",
                    transf = "none")

# extract posterior parameter samples
posterior_df <- as.data.frame(abc_fit$unadj.values)
colnames(posterior_df) <- c("k", "vth")
posterior_path <- file.path("output", "posterior_samples.csv")
write.csv(posterior_df, posterior_path, row.names = FALSE)

# --------------------------------------------------
# 8. Summarise posterior distributions
# --------------------------------------------------
# 
posterior_summary <- posterior_df %>%
  summarise(k_mean = mean(k),
            k_median = median(k),
            k_lower = quantile(k, 0.025),
            k_upper = quantile(k, 0.975),
    
            vth_mean = mean(vth),
            vth_median = median(vth),
            vth_lower = quantile(vth, 0.025),
            vth_upper = quantile(vth, 0.975))

print(posterior_summary)

write.csv(posterior_summary, file.path("output", "posterior_summary.csv"), 
          row.names = FALSE)


# --------------------------------------------------
# 10. Posterior predictive simulations
# --------------------------------------------------

out <- file.path("output")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

posterior <- posterior_df
pp <- t(vapply(seq_len(nrow(posterior)), function(i) {
  pred <- simulate_baseline_prevalence(k = posterior$k[i], vth = posterior$vth[i])
  c(k = posterior$k[i], vth = posterior$vth[i], worm_prev = pred)
}, numeric(3)))

pp <- as.data.frame(pp)
write.csv(pp, file.path(out, "posterior_predictive.csv"), row.names = FALSE)

cat("Observed baseline prevalence:", target_stat, "\n")
cat("Posterior predictive mean:", mean(pp$worm_prev), "\n")
cat("Posterior predictive median:", median(pp$worm_prev), "\n")
cat("Posterior predictive 95% interval:", quantile(pp$worm_prev, 0.025),
    "to", quantile(pp$worm_prev, 0.975), "\n")

# Save posterior predictive summary
pp_summary <- data.frame(
  observed_prevalence = target_stat,
  predictive_mean = mean(pp$worm_prev),
  predictive_median = median(pp$worm_prev),
  predictive_lower = quantile(pp$worm_prev, 0.025),
  predictive_upper = quantile(pp$worm_prev, 0.975)
)

write.csv(pp_summary, file.path(out, "posterior_predictive_summary.csv"), 
          row.names = FALSE)


posterior_prev <- results %>%
  mutate(distance = abs(prev - target_stat)) %>%
  arrange(distance) %>%
  slice_head(n = nrow(posterior_df))

cat("Observed prevalence:", target_stat, "\n")
cat("Valid simulations:", nrow(results), "\n")
cat("Posterior samples:", nrow(posterior_df), "\n")
cat("ABC tolerance:", abc_tol, "\n")

cat("Runtime:", difftime(Sys.time(), start.time, units = "mins"), "minutes\n")
