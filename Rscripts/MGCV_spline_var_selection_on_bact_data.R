# This script fits GAM models to the data using different spline
# basis types and iteratively refines the models by removing
# non-significant smooth terms until all remaining terms are significant.

# Load necessary libraries for generalized additive models and plotting
library(mgcv)
library(ggplot2)
library(tidyr)

# Load the dataset from a CSV file
Bact_data <- read.csv("data/Bacteremia_public_S2.csv", header = TRUE)

# Convert the 'BloodCulture' column into a binary response variable 'Y'
Bact_data$Y <- ifelse(Bact_data$BloodCulture == "no", 0, 1)

# Function to fit GAM models with automatic term selection and comparison
fit_gam_auto <- function(data, response_var, terms, smooth_terms_ps, smooth_terms_ts, smooth_terms_tp) {
  
  # Construct formulas for each model with different spline basis
  formula_ps <- as.formula(paste(response_var, "~", paste(paste(terms, collapse = " + "), paste(smooth_terms_ps, collapse = " + "), sep = "+")))
  formula_ts <- as.formula(paste(response_var, "~", paste(paste(terms, collapse = " + "), paste(smooth_terms_ts, collapse = " + "), sep = "+")))
  formula_tp <- as.formula(paste(response_var, "~", paste(paste(terms, collapse = " + "), paste(smooth_terms_tp, collapse = " + "), sep = "+")))
  
  # Fit the GAM models with each basis
  fit_ps <- gam(formula_ps, family = "binomial", data = data, select = TRUE, method = "REML") # p-splines
  fit_ts <- gam(formula_ts, family = "binomial", data = data, method = "REML") # thin plate regression with var selection 1
  fit_tp <- gam(formula_tp, family = "binomial", data = data, select = TRUE, method = "REML") # thin plate regression with var selection 1
  
  models <- list(ps = fit_ps, ts = fit_ts, tp = fit_tp)
  
  best_aic <- Inf # Initialize AIC for comparison
  best_model <- NULL
  
  # Iterate over each model to optimize and check term significance
  for (basis in c("ps", "ts", "tp")) {
    fit <- models[[basis]]
    step <- 1
    
    repeat {
      summary_fit <- summary(fit)
      
      # Extract p-values to determine significant terms
      p_values <- summary_fit$s.table[, "p-value"]
      significant_terms <- names(p_values)[p_values <= 0.05 & !is.na(p_values)]
      
      # Break if all terms are significant
      if (length(significant_terms) == length(fit$smooth)) {
        message(paste("All terms are significant for basis", basis, ". Model fitting complete."))
        break
      }
      
      # Update the formula with significant terms and refit
      smooth_terms <- significant_terms
      formula_string <- paste(response_var, "~", paste(paste(terms, collapse = " + "), paste(smooth_terms, collapse = " + "), sep = "+"))
      formula <- as.formula(formula_string)
      fit <- gam(formula, family = "binomial", data = data, method = "REML")
      step <- step + 1
    }
    
    models[[basis]] <- fit
    # Record the best model based on AIC
    if (fit$aic < best_aic) {
      best_aic <- fit$aic
      best_model <- fit
    }
  }
  
  return(list(ps_model = models[["ps"]], ts_model = models[["ts"]], tp_model = models[["tp"]], best_model = best_model))
}

# Specify linear terms for the fixed part of the model
terms <- "SEX"

# Specify smooth terms for each model with different basis functions
smooth_terms_ps <- c("s(AGE, bs = 'ps')", "s(MCV, bs = 'ps')", "s(MCHC, bs = 'ps')", 
                     "s(MONO, bs = 'ps')", "s(FIB, bs = 'ps')", "s(SODIUM, bs = 'ps')", 
                     "s(CA, bs = 'ps')", "s(PHOS, bs = 'ps')", "s(MG, bs = 'ps')", 
                     "s(CREA, bs = 'ps')", "s(BUN, bs = 'ps')", "s(HS, bs = 'ps')", 
                     "s(GBIL, bs = 'ps')", "s(ALB, bs = 'ps')", "s(AMY, bs = 'ps')", 
                     "s(AP, bs = 'ps')", "s(GGT, bs = 'ps')", "s(LDH, bs = 'ps')", 
                     "s(CHOL, bs = 'ps')", "s(CRP, bs = 'ps')", "s(EOSR, bs = 'ps')", 
                     "s(LYMR, bs = 'ps')")

# Change basis types for other models
smooth_terms_ts <- gsub("bs = 'ps'", "bs = 'ts'", smooth_terms_ps)
smooth_terms_tp <- gsub("bs = 'ps'", "bs = 'tp'", smooth_terms_ps)

# Measure the time taken to fit models and fit them ##
## This takes about 30' on my pc, not sure how slow it might be on other machines
system.time(models <- fit_gam_auto(Bact_data, "Y", terms, smooth_terms_ps, smooth_terms_ts, smooth_terms_tp))

# Print summaries for the best model and all models
summary(models$best_model)
summary(models$tp_model)
summary(models$ts_model)
summary(models$ps_model)

# Extract the names of significant smooth terms from each model for comparison
smooth_terms_tp <- rownames(summary(models$tp_model)$s.table)
smooth_terms_ts <- rownames(summary(models$ts_model)$s.table)
smooth_terms_ps <- rownames(summary(models$ps_model)$s.table)

# Generate a table to compare which terms are significant in each model
all_terms <- unique(c(smooth_terms_tp, smooth_terms_ts, smooth_terms_ps))
comparison_table <- data.frame(
  Term = all_terms,
  TP = all_terms %in% smooth_terms_tp,
  TS = all_terms %in% smooth_terms_ts,
  PS = all_terms %in% smooth_terms_ps
)

# Print or view the comparison table and sort it alphabetically
print(comparison_table)
comparison_table <- comparison_table[order(comparison_table$Term), ]
print(comparison_table)


### selected models
## 
#      Term    TP    TS    PS
#       SEX   TRUE  TRUE  TRUE  
#     s(AGE)  TRUE  TRUE  TRUE
#     s(AMY) FALSE FALSE  TRUE
#      s(AP)  TRUE  TRUE  TRUE
#     s(BUN)  TRUE  TRUE  TRUE
#      s(CA) FALSE  TRUE FALSE
#    s(CHOL) FALSE  TRUE FALSE
#     s(CRP)  TRUE  TRUE  TRUE
#    s(EOSR)  TRUE  TRUE  TRUE
#     s(FIB)  TRUE  TRUE FALSE
#    s(GBIL)  TRUE  TRUE  TRUE
#     s(GGT) FALSE FALSE  TRUE
#     s(LDH)  TRUE  TRUE  TRUE
#    s(LYMR)  TRUE  TRUE  TRUE
#    s(MCHC) FALSE FALSE  TRUE
#     s(MCV)  TRUE  TRUE  TRUE
#      s(MG)  TRUE  TRUE  TRUE
#    s(MONO)  TRUE  TRUE  TRUE
#    s(PHOS)  TRUE  TRUE FALSE
#  s(SODIUM)  TRUE  TRUE  TRUE

# Generate predictions for each model
preds_ps <- predict(models$ps_model, type = "terms", se = TRUE)
preds_tp <- predict(models$tp_model, type = "terms", se = TRUE)
preds_ts <- predict(models$ts_model, type = "terms", se = TRUE)

# Prepare model data frames by removing response variable and intercept
ps_m <- model.frame(models$ps_model)[, -c(1, 2)]
tp_m <- model.frame(models$tp_model)[, -c(1, 2)]
ts_m <- model.frame(models$ts_model)[, -c(1, 2)]

# Add prediction and standard error columns to each model's data frame
pred_se_ps <- preds_ps$se.fit
colnames(pred_se_ps) <- paste("se.", colnames(pred_se_ps), sep = "")
pred_se_tp <- preds_tp$se.fit
colnames(pred_se_tp) <- paste("se.", colnames(pred_se_tp), sep = "")
pred_se_ts <- preds_ts$se.fit
colnames(pred_se_ts) <- paste("se.", colnames(pred_se_ts), sep = "")

# Combine model frame data with predictions
ps_mf <- cbind(ps_m, preds_ps$fit, pred_se_ps)
tp_mf <- cbind(tp_m, preds_tp$fit, pred_se_tp)
ts_mf <- cbind(ts_m, preds_ts$fit, pred_se_ts)

# Function to standardize data for plotting smooth terms
prepare_data <- function(df, model_name, var_name) {
  # Construct column names using the variable name
  s_col_name <- paste0("s(", var_name, ")")
  se_col_name <- paste0("se.", s_col_name)
  
  df %>%
    select(
      !!var_name := all_of(var_name),  # Select the specified variable
      contains(s_col_name),            # Select the s(var_name) column
      contains(se_col_name)            # Select the se.s(var_name) column
    ) %>%
    rename(
      fit = all_of(s_col_name),        # Rename s(var_name) to fit
      se = all_of(se_col_name)         # Rename se.s(var_name) to se
    ) %>%
    mutate(
      lower = fit - 1.96 * se,
      upper = fit + 1.96 * se,
      model = model_name
    )
}

# Prepare data for plotting smooth terms for "MCV"
combined_data <- bind_rows(
  prepare_data(ps_mf, "PS Model", "MCV"),
  prepare_data(tp_mf, "TP Model", "MCV"),
  prepare_data(ts_mf, "TS Model", "MCV")
)

# Plot model predictions for "MCV"
ggplot(combined_data, aes(x = MCV, y = fit, color = model, fill = model)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, linetype = 2) +
  facet_wrap(~model) +  # Separate panels for each model
  labs(title = "Model Predictions for s(MCV) by Model",
       x = "MCV",
       y = "s(MCV) Fit") +
  theme_bw() # Use a clean theme

# Plot with ribbons only for "TP Model"
ggplot(combined_data, aes(x = MCV, y = fit, color = model)) +
  geom_line() +
  geom_ribbon(data = subset(combined_data, model == "TP Model"),  # Filter for TP Model
              aes(ymin = lower, ymax = upper, fill = model), 
              alpha = 0.1, linetype = 2) +
  labs(title = "Model Predictions for s(MCV)",
       x = "MCV",
       y = "s(MCV) Fit") +
  theme_bw()

# Repeat the preparation and plotting for a different variable "MG"
combined_data <- bind_rows(
  prepare_data(ps_mf, "PS Model", "MG"),
  prepare_data(tp_mf, "TP Model", "MG"),
  prepare_data(ts_mf, "TS Model", "MG")
)

# Plot model predictions for "MG" with a ribbon only for "TP Model"
ggplot(combined_data, aes(x = MG, y = fit, color = model)) +
  geom_line() +
  geom_ribbon(data = subset(combined_data, model == "TP Model"),  
              aes(ymin = lower, ymax = upper, fill = model), 
              alpha = 0.1, linetype = 2) +
  labs(title = "Model Predictions for s(MG)",
       x = "MG",
       y = "s(MG) Fit") +
  theme_bw()




# Plot model predictions for "MG"
ggplot(combined_data, aes(x = MG, y = fit, color = model, fill = model)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, linetype = 2) +
  facet_wrap(~model) +  # Separate panels for each model
  labs(title = "Model Predictions for s(MCV) by Model",
       x = "MG",
       y = "s(MG) Fit") +
  theme_bw() # Use a clean theme






# Extract the effective degrees of freedom (edf) and term names for comparison
get_terms_and_edf <- function(model_summary) {
  s_table <- model_summary$s.table
  terms <- rownames(s_table)
  edf <- s_table[, "edf"]
  return(data.frame(Term = terms, EDF = edf))
}

# Get terms and edf for each model
terms_tp <- get_terms_and_edf(summary(models$tp_model))
terms_ts <- get_terms_and_edf(summary(models$ts_model))
terms_ps <- get_terms_and_edf(summary(models$ps_model))

# Combine all unique terms
all_terms <- unique(c(terms_tp$Term, terms_ts$Term, terms_ps$Term))

# Initialize columns for term types (Line/Smooth) and presence (TRUE/FALSE)
comparison_info <- data.frame(
  Term = all_terms,
  TP = all_terms %in% terms_tp$Term,
  TP_Type = NA,
  TS = all_terms %in% terms_ts$Term,
  TS_Type = NA,
  PS = all_terms %in% terms_ps$Term,
  PS_Type = NA
)

# Function to determine if a term is line or smooth
get_term_type <- function(edf) {
  if (edf <= 1.5) "Line" else "Smooth"
}

# Populate the types based on edf values
for (i in seq_along(all_terms)) {
  term <- all_terms[i]
  
  if (comparison_info$TP[i]) {
    tp_edf <- terms_tp$EDF[terms_tp$Term == term]
    comparison_info$TP_Type[i] <- get_term_type(tp_edf)
  }
  
  if (comparison_info$TS[i]) {
    ts_edf <- terms_ts$EDF[terms_ts$Term == term]
    comparison_info$TS_Type[i] <- get_term_type(ts_edf)
  }
  
  if (comparison_info$PS[i]) {
    ps_edf <- terms_ps$EDF[terms_ps$Term == term]
    comparison_info$PS_Type[i] <- get_term_type(ps_edf)
  }
}

# Sort table alphabetically by term name
comparison_info <- comparison_info[order(comparison_info$Term), ]

# Print the comparison table with type information
print(comparison_info)
#> print(comparison_info)
#       Term    TP TP_Type    TS TS_Type    PS PS_Type
#     s(AGE)  TRUE  Smooth  TRUE    Line  TRUE    Line
#     s(AMY) FALSE    <NA> FALSE    <NA>  TRUE  Smooth
#      s(AP)  TRUE  Smooth  TRUE  Smooth  TRUE  Smooth
#      s(BUN)  TRUE    Line  TRUE  Smooth  TRUE  Smooth
#       s(CA) FALSE    <NA>  TRUE  Smooth FALSE    <NA>
#     s(CHOL) FALSE    <NA>  TRUE    Line FALSE    <NA>
#     s(CRP)  TRUE  Smooth  TRUE  Smooth  TRUE  Smooth
#    s(EOSR)  TRUE  Smooth  TRUE  Smooth  TRUE  Smooth
#     s(FIB)  TRUE  Smooth  TRUE  Smooth FALSE    <NA>
#    s(GBIL)  TRUE  Smooth  TRUE  Smooth  TRUE  Smooth
#     s(GGT) FALSE    <NA> FALSE    <NA>  TRUE  Smooth
#     s(LDH)  TRUE    Line  TRUE    Line  TRUE    Line
#    s(LYMR)  TRUE  Smooth  TRUE  Smooth  TRUE  Smooth
#    s(MCHC) FALSE    <NA> FALSE    <NA>  TRUE  Smooth
#     s(MCV)  TRUE    Line  TRUE    Line  TRUE    Line
#      s(MG)  TRUE    Line  TRUE    Line  TRUE  Smooth
#    s(MONO)  TRUE  Smooth  TRUE  Smooth  TRUE  Smooth
#    s(PHOS)  TRUE  Smooth  TRUE  Smooth FALSE    <NA>
#  s(SODIUM)  TRUE  Smooth  TRUE  Smooth  TRUE  Smooth
