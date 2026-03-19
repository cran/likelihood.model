## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>",
    fig.width = 6,
    fig.height = 4
)
old_opts <- options(digits = 4)

## ----install, eval=FALSE------------------------------------------------------
# # Install from GitHub
# if (!require(devtools)) install.packages("devtools")
# devtools::install_github("queelius/likelihood.model")

## ----load, message=FALSE, warning=FALSE---------------------------------------
library(likelihood.model)

## ----exponential-example------------------------------------------------------
# Generate exponential survival data
set.seed(99)
df_exp <- data.frame(t = rexp(200, rate = 3.0))

# Create model and fit -- no initial guess needed!
model_exp <- exponential_lifetime("t")

# View model assumptions
assumptions(model_exp)

# Fit the model
mle_exp <- fit(model_exp)(df_exp)

# View results
summary(mle_exp)

## ----exponential-results------------------------------------------------------
# Parameter estimates
cat("Estimated parameters:\n")
print(coef(mle_exp))

# Confidence intervals (Wald-based)
cat("\n95% Confidence Intervals:\n")
print(confint(mle_exp))

# Standard errors
cat("\nStandard Errors:\n")
print(se(mle_exp))

# Log-likelihood value
cat("\nLog-likelihood:", as.numeric(logLik(mle_exp)), "\n")

# AIC for model selection
cat("AIC:", AIC(mle_exp), "\n")

# The score at the MLE is exactly zero (by construction)
cat("Score at MLE:", score_val(mle_exp), "\n")

## ----exponential-censored-----------------------------------------------------
# Generate data with right-censoring at time 0.5
set.seed(99)
true_lambda <- 3.0
x <- rexp(200, rate = true_lambda)
censored <- x > 0.5
df_cens <- data.frame(
  t = pmin(x, 0.5),
  status = ifelse(censored, "right", "exact")
)

cat("Censoring rate:", mean(censored) * 100, "%\n")

model_cens <- exponential_lifetime("t", censor_col = "status")
mle_cens <- fit(model_cens)(df_cens)

cat("MLE (censored):", coef(mle_cens), "(true:", true_lambda, ")\n")
cat("95% CI:", confint(mle_cens)[1, ], "\n")

## ----score-check--------------------------------------------------------------
model <- exponential_lifetime("t")
df <- data.frame(t = rexp(100, rate = 2.0))

# Fit MLE
mle <- fit(model)(df)

# Evaluate score at MLE
score_func <- score(model)
score_at_mle <- score_func(df, coef(mle))

cat("Score at MLE (should be near zero):\n")
print(score_at_mle)

# The score is also stored in the MLE object
cat("\nScore stored in MLE object:\n")
print(score_val(mle))

## ----fisherian-example--------------------------------------------------------
# Fit a model
model <- exponential_lifetime("t")
df <- data.frame(t = rexp(200, rate = 2.0))
mle <- fit(model)(df)

# Support function: log relative likelihood
# S(theta) = logL(theta) - logL(theta_hat)
# Support at MLE is always 0
s_at_mle <- support(mle, coef(mle), df, model)
cat("Support at MLE:", s_at_mle, "\n")

# Support at alternative parameter values (negative = less supported)
s_alt <- support(mle, c(lambda = 1.0), df, model)
cat("Support at lambda=1.0:", s_alt, "\n")

# Relative likelihood = exp(support)
rl <- relative_likelihood(mle, c(lambda = 1.0), df, model)
cat("Relative likelihood at lambda=1.0:", rl, "\n")

## ----likelihood-intervals-----------------------------------------------------
# Compute 1/8 likelihood interval (roughly equivalent to 95% CI)
# This is the set of theta where R(theta) >= 1/8
li <- likelihood_interval(mle, df, model, k = 8)
print(li)

# Compare with Wald confidence interval
cat("\nWald 95% CI:\n")
print(confint(mle))

## ----evidence-example---------------------------------------------------------
model <- exponential_lifetime("t")
set.seed(42)
df <- data.frame(t = rexp(100, rate = 2.0))

# Evidence for lambda=2 vs lambda=1
ev <- evidence(model, df, c(lambda = 2), c(lambda = 1))
cat("Evidence for lambda=2 vs lambda=1:", ev, "\n")

# Positive evidence supports the first hypothesis
if (ev > log(8)) {
  cat("Strong evidence favoring lambda=2\n")
}

## ----bootstrap-example, cache=TRUE--------------------------------------------
model <- exponential_lifetime("t")
set.seed(42)
df <- data.frame(t = rexp(100, rate = 2.0))

# Create bootstrap sampler
boot_sampler <- sampler(model, df = df, par = c(lambda = 2))

# Generate bootstrap samples (100 replicates for speed)
boot_result <- boot_sampler(n = 100)

# Compare asymptotic vs bootstrap standard errors
mle <- fit(model)(df)
asymp_se <- se(mle)
boot_se <- se(boot_result)

cat("Standard Error Comparison:\n")
cat("           Asymptotic  Bootstrap\n")
cat(sprintf("Lambda:    %10.4f  %9.4f\n", asymp_se[1], boot_se[1]))

# Bootstrap bias estimate
cat("\nBootstrap Bias Estimate:\n")
print(bias(boot_result))

## ----session------------------------------------------------------------------
sessionInfo()

## ----cleanup, include=FALSE---------------------------------------------------
options(old_opts)

