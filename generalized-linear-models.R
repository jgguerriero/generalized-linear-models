
# general linear model ------

# setting seed and simulating data
set.seed(973)

lm_df = data.frame(y = rnorm(100),x1 = rnorm(100), x2 = rnorm(100))

## coefs -----
# running model with lm function
lm_mod = lm(y~ x1 + x2, lm_df)

# deriving coefficients manually
X_matrix = matrix(c(rep(1, 100),lm_df$x1,lm_df$x2),
                  nrow=100, ncol=3)

betas = solve(t(X_matrix) %*% X_matrix) %*% t(X_matrix) %*% lm_df$y

# compare manual betas to lm betas 
betas;coef(lm_mod)

## SEs -----
# deriving standard errors manually
predicted_ys = X_matrix %*% betas

se = sqrt(diag((sum((lm_df$y - predicted_ys)^2)/(100-2-1)) *
            solve(t(X_matrix) %*% X_matrix)))

# compare manual SEs to lm SEs 
se;summary(lm_mod)$coefficients[, "Std. Error"]

## translating lm to 3x2 ANOVA stats -----
## type 2/3 sums of squares -----

lm_df$group1 = sample(c("A","B"), 100, replace=T)
lm_df$group2 = sample(c("X","Y", "Z"), 100, replace=T)

# this is the model we want to translate to ANOVA
lm_4_anova_mod = lm(y ~ group1 + group2, lm_df)
summary(lm_4_anova_mod)

# partition sums of squares

t_sse = sum((lm_df$y - predict(lm(y~1, lm_df)))^2);t_sse

sse_both_groups = t_sse - sum((lm_df$y - predict(lm(y~1 + group1 + group2, lm_df)))^2);sse_both_groups

sse_group1 = t_sse - sum((lm_df$y - predict(lm(y~1 + group1, lm_df)))^2); sse_group1

sse_group2 = t_sse - sum((lm_df$y - predict(lm(y~1 + group2, lm_df)))^2);sse_group2

common_sse = (sse_group1+sse_group2) - sse_both_groups;common_sse

sse_group1 - common_sse
sse_group2 - common_sse

# compare sses to anova mod
car::Anova(lm_4_anova_mod, type = 2) # type 3 gives same results

# how one decides to partition sums of squares makes a difference, but sometimes it doesn't
# https://onlinestatbook.com/2/analysis_of_variance/unequal.html


## weighted coefs ----
# if there were weights to apply, the formula below acccounts for this in the solution to the coefficients

lm_df$weights = rnorm(100)^2

weighted_betas = solve(t(X_matrix) %*% diag(lm_df$weights) %*% X_matrix) %*% t(X_matrix) %*% diag(lm_df$weights) %*% lm_df$y

weighted_lm_mod = lm(y~ x1 + x2, lm_df, weights = weights)

# compare manual and lm weighed coefficients 
coef(weighted_lm_mod);weighted_betas

## weighted SEs ----
# deriving standard errors of weighed coefs manually

predicted_weighted_ys = X_matrix %*% weighted_betas

weighted_se = sqrt(diag((sum(lm_df$weights*(lm_df$y - predicted_weighted_ys)^2)/(100-2-1)) *
                 solve(t(X_matrix) %*% diag(lm_df$weights) %*% X_matrix)))

# compare manual weighted SEs to lm weighted SEs 
weighted_se;summary(weighted_lm_mod)$coefficients[, "Std. Error"]


# logistic regression -----
##  coefs -----
lr_df = data.frame(y = rep(c(0,1), times = 50),x1 = rnorm(100), x2 = rnorm(100))

# running model with glm function
lr_mod = glm(y ~ x1 + x2, family = "binomial", lr_df)

# get analytic solution for coefficients
X_matrix = matrix(c(rep(1, 100),lr_df$x1,lr_df$x2),
                  nrow=100, ncol=3)

logistic_regressions_coefs = function(params){
  # params = c(.1,.1,.1)
  
  betas = matrix(c(params[1],params[2],params[3]),
                      ncol = 1, nrow = 3)
  
  predicted_ys_log_odds = X_matrix %*% betas
  
  predicted_probabilities = exp(predicted_ys_log_odds) / (1 + exp(predicted_ys_log_odds))
  
  log_likelihoods = ifelse(lr_df$y == 1, log(c(predicted_probabilities)), log(1-predicted_probabilities))
  
  log_likelihood = -1*sum(log_likelihoods)
  
  return(log_likelihood)
  
}

optimized_betas = optim(par = c(.1,.1,.1), fn = logistic_regressions_coefs)$par

coef(lr_mod);optimized_betas

##  SEs -----
# deriving standard errors manually

predicted_ys_log_odds = X_matrix %*% optimized_betas
predicted_probs = c(exp(predicted_ys_log_odds) / (1 + exp(predicted_ys_log_odds)))

se = sqrt(diag(solve(t(X_matrix) %*% diag(predicted_probs * (1-predicted_probs)) %*% X_matrix)))

# compare manual SEs to lm SEs 
se;summary(lr_mod)$coefficients[, "Std. Error"]
