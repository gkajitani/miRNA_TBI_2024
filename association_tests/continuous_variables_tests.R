library (tidyverse)
library(dplyr)   # For data manipulation
library(tidyr)   # For data tidying

setwd ("~/clinical_variables/")


####################################
# lm regression - cat-cont- miRge3 #
####################################

vars <- read.csv("miR_9_CONT_EPT_CONT.csv", row.names = 1)
df <- vars

results <- NULL
for (col in names(df)[-1]) {
  model <- lm(MiR_9_5p ~ df[[col]], data = df, na.action = "na.exclude")
  coef <- summary(model)$coefficients[2,]
  predicted <- predict(model, df)
  conf_int <- confint(model)[2,]
  ci_string <- paste0(round(coef[1], 3), " (", round(conf_int[1], 3), ", ", round(conf_int[2], 3), ")")
  new_row <- data.frame(column = col,
                        coef_ci = ci_string,
                        std_error = coef[2],
                        z_value = coef[3],
                        p_value = coef[4],
                        stringsAsFactors = FALSE)
  
  results <- rbind(results, new_row)
}

results
dfR <- results



# View the results
results
dfR <- results
write.csv(dfR, "cont_lm_Results.csv")

#
#############################################
# lm regression without outliers - cat-cont #
#############################################

vars <- read.csv("miR_9_CONT_EPT_CONT_2.csv", row.names = 1)
df <- vars


# Create an empty dataframe to store results
results <- data.frame(column = character(),
                      coefficient = numeric(),
                      std_error = numeric(),
                      z_value = numeric(),
                      p_value = numeric(),
                      stringsAsFactors = FALSE)

results <- NULL
for (col in names(df)[-1]) {
  model <- lm(MiR_9_5p ~ df[[col]], data = df, na.action = "na.exclude")
  coef <- summary(model)$coefficients[2,]
  predicted <- predict(model, df)
  conf_int <- confint(model)[2,]
  ci_string <- paste0(round(coef[1], 3), " (", round(conf_int[1], 3), ", ", round(conf_int[2], 3), ")")
  new_row <- data.frame(column = col,
                        coef_ci = ci_string,
                        std_error = coef[2],
                        z_value = coef[3],
                        p_value = coef[4],
                        stringsAsFactors = FALSE)
  
  results <- rbind(results, new_row)
}


# View the results
results
dfR <- results
write.csv(dfR, "cont_lm_Results2.csv")


####################
# lm - categorical #
####################


vars <- read.csv("miR_9_CONT_cat_EPT2.csv", row.names = 1)
df <- vars

results <- NULL
for (col in names(df)[-1]) {
  model <- lm(MiR_9_5p ~ df[[col]], data = df, na.action = "na.exclude")
  coef <- summary(model)$coefficients[2,]
  predicted <- predict(model, df)
  conf_int <- confint(model)[2,]
  ci_string <- paste0(round(coef[1], 3), " (", round(conf_int[1], 3), ", ", round(conf_int[2], 3), ")")
  new_row <- data.frame(column = col,
                        coef_ci = ci_string,
                        std_error = coef[2],
                        z_value = coef[3],
                        p_value = coef[4],
                        stringsAsFactors = FALSE)
  
  results <- rbind(results, new_row)
}

results
dfR <- results



# View the results
results
dfR <- results


#####################################
# lm regression - cat-cont- exceRpt #
#####################################

vars <- read.csv("miR_9_CONT_EPT_CONT_EXCERPT.csv", row.names = 1)
df <- vars

results <- NULL
for (col in names(df)[-1]) {
  model <- lm(MiR_9_5p ~ df[[col]], data = df, na.action = "na.exclude")
  coef <- summary(model)$coefficients[2,]
  predicted <- predict(model, df)
  conf_int <- confint(model)[2,]
  ci_string <- paste0(round(coef[1], 3), " (", round(conf_int[1], 3), ", ", round(conf_int[2], 3), ")")
  new_row <- data.frame(column = col,
                        coef_ci = ci_string,
                        std_error = coef[2],
                        z_value = coef[3],
                        p_value = coef[4],
                        stringsAsFactors = FALSE)
  
  results <- rbind(results, new_row)
}

results
dfR <- results



# View the results
results
dfR <- results
write.csv(dfR, "cont_lm_Results_exceRpt.csv")


a <- lm(MiR_9_5p ~ Tratamento, data = df, na.action = "na.exclude")
coef <- summary(a)$coefficients[2,]
p_value = coef[4]

b <- as.data.frame(a)

####################
# lm - categorical #
####################


vars <- read.csv("miR_9_CONT_cat_EPT2_EXCERPT.csv", row.names = 1)
df <- vars

results <- NULL
for (col in names(df)[-1]) {
  model <- lm(MiR_9_5p ~ df[[col]], data = df, na.action = "na.exclude")
  coef <- summary(model)$coefficients[2,]
  predicted <- predict(model, df)
  conf_int <- confint(model)[2,]
  ci_string <- paste0(round(coef[1], 3), " (", round(conf_int[1], 3), ", ", round(conf_int[2], 3), ")")
  new_row <- data.frame(column = col,
                        coef_ci = ci_string,
                        std_error = coef[2],
                        z_value = coef[3],
                        p_value = coef[4],
                        stringsAsFactors = FALSE)
  
  results <- rbind(results, new_row)
}


for (col in names(df)[-1]) {
  tryCatch({
    model <- lm(MiR_9_5p ~ df[[col]], data = df, na.action = "na.exclude")
    coef <- summary(model)$coefficients[2,]
    predicted <- predict(model, df)
    conf_int <- confint(model)[2,]
    ci_string <- paste0(round(coef[1], 3), " (", round(conf_int[1], 3), ", ", round(conf_int[2], 3), ")")
    new_row <- data.frame(column = col,
                          coef_ci = ci_string,
                          std_error = coef[2],
                          z_value = coef[3],
                          p_value = coef[4],
                          stringsAsFactors = FALSE)
    
    results <- rbind(results, new_row)
  }, error = function(e) {
    cat("Error in model for column", col, ":", conditionMessage(e), "\n")
  })
}

results
dfR2 <- results

write.csv(dfR2, "excerpt_miR9_cat.csv")

# View the results
results
dfR2 <- results


a <- lm(MiR_9_5p ~ EPT, data = df, na.action = "na.exclude")
coef <- summary(a)$coefficients[2,]
p_value = coef[4]


##################
# Kruskal-Wallis #
##################


vars <- read.csv("miR_9_CONT_cat_EPT2.csv", row.names = 1)
df <- vars

# Create an empty dataframe to store the results
results <- data.frame()

# Loop through the columns of the dataframe, excluding the first column (the outcome variable)
for (col in names(df)[-1]) {
  # Perform a Kruskal-Wallis test comparing the outcome variable to the current column
  test_result <- kruskal.test(df[[col]], df$MiR_9_5p)
  
  # Extract the p-value from the test result
  p_value <- test_result$p.value
  
  # Create a new row with the column name and the p-value
  new_row <- data.frame(column = col,
                        p_value = p_value,
                        stringsAsFactors = FALSE)
  
  # Append the new row to the results dataframe
  results <- rbind(results, new_row)
}

dfR2 <-results
write.csv(dfR2, "Kruskal-Wallis.csv")
