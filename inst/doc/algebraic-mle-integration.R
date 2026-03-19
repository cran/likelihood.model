## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>",
    fig.width = 6,
    fig.height = 4
)
old_opts <- options(digits = 4)

## ----load, message=FALSE, warning=FALSE---------------------------------------
library(likelihood.model)
library(algebraic.mle)

## ----basic-fit----------------------------------------------------------------
set.seed(42)
true_lambda <- 2.0
n <- 200
df <- data.frame(t = rexp(n, rate = true_lambda))

model <- exponential_lifetime("t")
mle_result <- fit(model)(df)

# algebraic.mle generics -- these work because fisher_mle inherits from mle
params(mle_result)
nparams(mle_result)

## ----observed-fim-------------------------------------------------------------
observed_fim(mle_result)

## ----fim-pd-------------------------------------------------------------------
cat("Observed FIM:", observed_fim(mle_result)[1,1], "\n")
cat("Positive:", observed_fim(mle_result)[1,1] > 0, "\n")

## ----vcov-vs-fim--------------------------------------------------------------
vcov(mle_result)
1 / observed_fim(mle_result)

## ----rmap---------------------------------------------------------------------
# Transform rate -> mean lifetime
mean_life_mle <- rmap(mle_result, function(p) {
  c(mean_lifetime = 1 / p[1])
})

params(mean_life_mle)
se(mean_life_mle)

## ----rmap-compare-------------------------------------------------------------
true_mean <- 1 / true_lambda
cat("True mean lifetime:", true_mean, "\n")
cat("Estimated mean lifetime:", params(mean_life_mle), "\n")
cat("95% CI:", confint(mean_life_mle), "\n")

## ----rmap-multi---------------------------------------------------------------
# Derive mean, variance, and median of the exponential distribution
derived_mle <- rmap(mle_result, function(p) {
  lam <- p[1]
  c(
    mean   = 1 / lam,
    var    = 1 / lam^2,
    median = log(2) / lam
  )
})

params(derived_mle)
se(derived_mle)

## ----expectation--------------------------------------------------------------
set.seed(123)
# E[lambda^2] under the asymptotic distribution
e_lam_sq <- expectation(mle_result, function(p) p[1]^2,
                         control = list(n = 10000L))
cat("E[lambda^2]:", e_lam_sq, "\n")
cat("lambda^2 at MLE:", params(mle_result)[1]^2, "\n")

# Probability that lambda > 1.5 (under asymptotic distribution)
pr_lam_gt <- expectation(mle_result, function(p) as.numeric(p[1] > 1.5),
                          control = list(n = 10000L))
cat("P(lambda > 1.5):", pr_lam_gt, "\n")

## ----mse----------------------------------------------------------------------
mse(mle_result)
all.equal(mse(mle_result), vcov(mle_result))

## ----bootstrap, cache=TRUE----------------------------------------------------
set.seed(42)
boot_sampler <- sampler(model, df = df, par = c(lambda = 2))
boot_result <- boot_sampler(n = 200)

# Same algebraic.mle generics work
params(boot_result)
nparams(boot_result)
se(boot_result)
bias(boot_result)

## ----boot-compare-------------------------------------------------------------
cat("Asymptotic 95% CI:\n")
confint(mle_result)

cat("\nBootstrap percentile 95% CI:\n")
confint(boot_result, type = "perc")

## ----dist-check, include=FALSE------------------------------------------------
has_dist <- requireNamespace("algebraic.dist", quietly = TRUE)

## ----dist-section, eval=has_dist----------------------------------------------
library(algebraic.dist)

## ----dist-compare, eval=has_dist----------------------------------------------
set.seed(42)

cat("MLE:", params(mle_result), "\n")
cat("SE:", se(mle_result), "\n")

# Theoretical asymptotic distribution of the MLE
asymp_var <- true_lambda^2 / n
asymp_dist <- normal(mu = true_lambda, var = asymp_var)
cat("\nTheoretical asymptotic distribution:\n")
cat("  Mean:", params(asymp_dist)[1], "\n")
cat("  Variance:", params(asymp_dist)[2], "\n")

# Compare: sample from the MLE's estimated distribution
mle_sampler <- sampler(mle_result)
set.seed(1)
mle_samples <- mle_sampler(5000)

# vs. sample from the theoretical distribution
dist_sampler <- sampler(asymp_dist)
set.seed(1)
dist_samples <- dist_sampler(5000)

cat("\nMLE sampler mean:", mean(mle_samples), "\n")
cat("Theoretical sampler mean:", mean(dist_samples), "\n")
cat("MLE sampler sd:", sd(mle_samples), "\n")
cat("Theoretical sampler sd:", sd(dist_samples), "\n")

## ----cleanup, include=FALSE---------------------------------------------------
options(old_opts)

