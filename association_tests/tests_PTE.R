library (tidyverse)
library(purrr)
library(broom)
library (datasets)
library(rcompanion)


setwd ("~/clinical_variables/")

vars <- read.csv("cat_PTE.csv", row.names = 1)


df <- vars

# Specify the column to be tested
target_col <- "MiR_9_5p"

# Initialize the results data frame
results <- data.frame(variable = character(),
                      chisq_stat = numeric(),
                      p_value = numeric(),
                      cramer_v = numeric(),
                      stringsAsFactors = FALSE)

# Loop through each column in the data frame
for (col in names(df)[-1]) {
  
  # Skip the target column
  if (col == target_col) next
  
  # Create a contingency table
  tab <- table(df[[target_col]], df[[col]])
  
  # Perform the chi-squared test and calculate Cramer's V
  chisq <- chisq.test(tab)
  cramer_v <- cramerV(tab)
  
  # Store the results in the data frame
  new_row <- data.frame(variable = col,
                        chisq_stat = chisq$statistic,
                        p_value = chisq$p.value,
                        cramer_v = cramer_v,
                        stringsAsFactors = FALSE)
  
  results <- rbind(results, new_row)
}
results
dfR1 <- results

write.csv(dfR1, "chi_square_Results.csv")


##################################
# logit regression - categorical #
##################################

vars <- read.csv("cat_PTE.csv", row.names = 1)
df <- vars

# Create an empty dataframe to store results
results <- data.frame(column = character(),
                      coefficient = numeric(),
                      std_error = numeric(),
                      z_value = numeric(),
                      p_value = numeric(),
                      stringsAsFactors = FALSE)

model <- glm(PTE ~ Alcohol, data = df, family = binomial, na.action = "na.exclude")

coef <- summary(model)$coefficients[2,]
coef[4]

# Loop through columns B to E and perform logistic regression
for (col in names(df)[-1]) {
  model <- glm(PTE ~ df[[col]], data = df, family = binomial, na.action = "na.exclude")
  coef <- summary(model)$coefficients[2,]
  new_row <- data.frame(column = col,
                        coefficient = coef[1],
                        std_error = coef[2],
                        z_value = coef[3],
                        p_value = coef[4],
                        stringsAsFactors = FALSE)
  results <- rbind(results, new_row)
}

# View the results
results
dfR2 <- results
write.csv(dfR2, "logit_Results_PTE.csv")

results <- lapply(colnames(df)[-1], function(col) {
  model <- glm(PTE ~ df[[col]], data = df,family = binomial, na.action = "na.exclude")
  odds_ratio <- exp(coef(model)[2])
  ci <- confint(model)[2,]
  ci_low <- exp(ci[1])
  ci_high <- exp(ci[2])
  or_ci <- paste(round(odds_ratio, 3), "(", round(ci_low, 3), "-", round(ci_high, 3), ")")
  return(data.frame(column = col, odds_ratio = or_ci))
})
results_df <- do.call(rbind, results)

write.csv(results_df, "OddsRatio_logit_PTE1.csv")



###############################
# logit regression - cat-cont #
###############################

vars <- read.csv("cont_PTE.csv", row.names = 1)
df <- vars


# Create an empty dataframe to store results
results <- data.frame(column = character(),
                      coefficient = numeric(),
                      std_error = numeric(),
                      z_value = numeric(),
                      p_value = numeric(),
                      stringsAsFactors = FALSE)

# Loop through columns B to E and perform logistic regression
for (col in names(df)[-1]) {
  model <- glm(PTE ~ df[[col]], data = df, family = binomial, na.action = "na.exclude")
  coef <- summary(model)$coefficients[2,]
  new_row <- data.frame(column = col,
                        coefficient = coef[1],
                        std_error = coef[2],
                        z_value = coef[3],
                        p_value = coef[4],
                        stringsAsFactors = FALSE)
  results <- rbind(results, new_row)
}

# View the results
results
dfR3 <- results
write.csv(dfR3, "glm_Results_PTE.csv")


results <- lapply(colnames(df)[-1], function(col) {
  model <- glm(PTE ~ df[[col]], family = binomial, data = df, na.action = "na.exclude")
  odds_ratio <- exp(coef(model)[2])
  ci <- confint(model)[2,]
  ci_low <- exp(ci[1])
  ci_high <- exp(ci[2])
  or_ci <- paste(round(odds_ratio, 3), "(", round(ci_low, 3), "-", round(ci_high, 3), ")")
  return(data.frame(column = col, odds_ratio = or_ci))
})
results_df <- do.call(rbind, results)
write.csv(results_df, "OddsRatio_glm_PTE.csv")

warnings()


##############################
# OPT #######################
#############################


###############################
# logit regression - cat-cont #
###############################

vars <- read.csv("OTI.csv", row.names = 1)
df <- vars


# Create an empty dataframe to store results
results <- data.frame(column = character(),
                      coefficient = numeric(),
                      std_error = numeric(),
                      z_value = numeric(),
                      p_value = numeric(),
                      stringsAsFactors = FALSE)

# Loop through columns B to E and perform logistic regression
for (col in names(df)[-1]) {
  model <- glm(PTE ~ df[[col]], data = df, family = binomial, na.action = "na.exclude")
  coef <- summary(model)$coefficients[2,]
  new_row <- data.frame(column = col,
                        coefficient = coef[1],
                        std_error = coef[2],
                        z_value = coef[3],
                        p_value = coef[4],
                        stringsAsFactors = FALSE)
  results <- rbind(results, new_row)
}

# View the results
results
dfR3 <- results
write.csv(dfR3, "glm_Results_PTE_OTI.csv")


results <- lapply(colnames(df)[-1], function(col) {
  model <- glm(PTE ~ df[[col]], family = binomial, data = df, na.action = "na.exclude")
  odds_ratio <- exp(coef(model)[2])
  ci <- confint(model)[2,]
  ci_low <- exp(ci[1])
  ci_high <- exp(ci[2])
  or_ci <- paste(round(odds_ratio, 3), "(", round(ci_low, 3), "-", round(ci_high, 3), ")")
  return(data.frame(column = col, odds_ratio = or_ci))
})
results_df <- do.call(rbind, results)
write.csv(results_df, "OddsRatio_glm_PTE_OTI.csv")



####################################


model <- glm(PTE ~ XXXX, data = df,family = binomial, na.action = "na.exclude")
odds_ratio <- exp(coef(model)[2])
ci <- confint(model)[2,]
ci_low <- exp(ci[1])
ci_high <- exp(ci[2])
or_ci <- paste(round(odds_ratio, 3), "(", round(ci_low, 3), "-", round(ci_high, 3), ")")
