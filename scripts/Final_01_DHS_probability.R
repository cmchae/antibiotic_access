# main analysis 
rm(list = ls())
gc()

requiredPackages = c("rdhs", "dplyr", "tidyr", "janitor", "epitools", "tibble","sandwich","lmtest", "gt")
lib.path <- .libPaths()
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p, lib = lib.path)
  library(p,character.only = TRUE)
}

set.seed(123)
options(scipen = 999, digits = 3, print.max = 1000)

analysis_data <- readRDS("outputs/analysis_data.rds")
str(analysis_data)

## =============crude model with robust SE: glm (by country)===============

str(analysis_data)
analysis_data_crude <- analysis_data %>% dplyr::select(c("Antibiotics", "ARI", "Cough",  "Fever", "CountryYear"))
analysis_data_crude <- na.omit(analysis_data_crude)
analysis_data_crude <- analysis_data_crude %>% mutate(PersonTime = 1) 
analysis_data_crude <- analysis_data_crude %>% mutate(ARI = factor(ARI, levels = c(3, 0, 1, 2)))
analysis_data_crude <- analysis_data_crude %>% mutate(Cough = factor(Cough, levels = c(1,0)))
analysis_data_crude_2 <- analysis_data_crude
analysis_data_crude <- analysis_data_crude %>% mutate(Fever = factor(Fever, levels = c(1, 0)))
analysis_data_crude_nonFever <- analysis_data_crude[analysis_data_crude$Fever==0,]
analysis_data_crude_Fever <- analysis_data_crude[analysis_data_crude$Fever==1,]
analysis_data_crude_nonARI <- analysis_data_crude[analysis_data_crude$ARI!=3,]

country_list <- unique(analysis_data_crude$CountryYear)  # Get the unique list of CountryYear combinations

table_country_Ab_ARI_glm <- # Initialize an empty tibble to store results
  tibble::tibble(CountryYear = character(), Rate = numeric(), Lower = numeric(), Upper = numeric(), SE = numeric(), numerator = numeric(), denominator = numeric())

table_country_Fever_Cough_glm <- 
  table_country_Fever_RapidBreath_glm <- 
  table_country_Fever_ARI_glm <- 
  table_country_Ab_Both_glm <- 
  table_country_Ab_nonFever_glm <- 
  table_country_Ab_nonARI_glm <- 
  table_country_Ab_Fever_glm <- 
  table_country_Ab_ARI_glm



# -----------------P(Antibiotics | Syndrome)------------------

# P(Antibiotics | ARI with or without Fever)
for (country in country_list) {  # Loop through each CountryYear
  country_data <- analysis_data_crude %>% filter(CountryYear == country)   # Filter data for the current country-year

  model <- tryCatch(   # Fit the Poisson regression model
    glm(Antibiotics ~ ARI, 
        family = poisson(link = "log"), offset = log(PersonTime), data = country_data),
    error = function(e) {
      cat("Error fitting model for CountryYear:", country, "-", e$message, "\n")
      return(NULL)
    }
  )
  robust_se <- vcovHC(model, type = "HC0")
  
  n_LRI <- nrow(country_data[country_data$ARI == 3,])
  n_Ab <- nrow(country_data[country_data$ARI == 3 & country_data$Antibiotics == 1,])
  
  if (is.null(model)) next  # Skip if the model failed
  model_robust <- coeftest(model, vcov. = robust_se)
  effect <- model_robust[1, "Estimate"]
  se <- model_robust[1, "Std. Error"]
  
  z_value <- 1.96  # 95% confidence level
  lower <- effect - z_value * se
  upper <- effect + z_value * se
  
  rate <- exp(effect)  # Convert to incidence rate
  lower <- exp(lower)  # Lower bound
  upper <- exp(upper)  # Upper bound
  
  
  table_country_Ab_ARI_glm <- table_country_Ab_ARI_glm %>%   # Add results to the summary table
    add_row(CountryYear = country, Rate = rate, Lower = lower, Upper = upper, SE = se, numerator = n_Ab, denominator = n_LRI)
}
AbLRI <- table_country_Ab_ARI_glm
colnames(AbLRI) <- c("CountryYear", "AbLRI", "AbLRI_L", "AbLRI_U", "AbLRI_SE", "n_AbLRI", "n_LRI")



# P(Antibiotics | Fever or ARI)
for (country in country_list) {  # Loop through each CountryYear
  country_data <- analysis_data_crude %>% filter(CountryYear == country)   # Filter data for the current country-year
  
  model <- tryCatch(   # Fit the Poisson regression model
    glm(Antibiotics ~ Fever, 
        family = poisson(link = "log"), offset = log(PersonTime), data = country_data),
    error = function(e) {
      cat("Error fitting model for CountryYear:", country, "-", e$message, "\n")
      return(NULL)
    }
  )
  robust_se <- vcovHC(model, type = "HC0")

  n_Fever_or_LRI <- nrow(country_data[country_data$Fever == 1,])
  n_Ab <- nrow(country_data[country_data$Fever == 1 & country_data$Antibiotics == 1,])
    
  if (is.null(model)) next  # Skip if the model failed
  model_robust <- coeftest(model, vcov. = robust_se)
  effect <- model_robust[1, "Estimate"]
  se <- model_robust[1, "Std. Error"]
  
  z_value <- 1.96  # 95% confidence level
  lower <- effect - z_value * se
  upper <- effect + z_value * se
  
  rate <- exp(effect)  # Convert to incidence rate
  lower <- exp(lower)  # Lower bound
  upper <- exp(upper)  # Upper bound
  
  table_country_Ab_Fever_glm <- table_country_Ab_Fever_glm %>%   # Add results to the summary table
    add_row(CountryYear = country, Rate = rate, Lower = lower, Upper = upper, SE = se, numerator = n_Ab, denominator = n_Fever_or_LRI)
}
AbFever<- table_country_Ab_Fever_glm
colnames(AbFever) <- c("CountryYear", "AbFever", "AbFever_L", "AbFever_U", "AbFever_SE", "n_Ab_Fever_or_ARI", "n_Fever_or_ARI")



# P(Antibiotics | Fever without ARI)
for (country in country_list) {  # Loop through each CountryYear
  country_data <- analysis_data_crude_nonARI %>% filter(CountryYear == country)  # Filter data for the current country-year
  
  model <- tryCatch(   # Fit the Poisson regression model
    glm(Antibiotics ~ Fever, 
        family = poisson(link = "log"), offset = log(PersonTime), data = country_data),
    error = function(e) {
      cat("Error fitting model for CountryYear:", country, "-", e$message, "\n")
      return(NULL)
    }
  )
  robust_se <- vcovHC(model, type = "HC0")

  n_Fever_or_LRI <- nrow(country_data[country_data$Fever == 1,])
  n_Ab <- nrow(country_data[country_data$Fever == 1 & country_data$Antibiotics == 1,])
  
  if (is.null(model)) next  # Skip if the model failed
  model_robust <- coeftest(model, vcov. = robust_se)
  effect <- model_robust[1, "Estimate"]
  se <- model_robust[1, "Std. Error"]
  
  z_value <- 1.96  # 95% confidence level
  lower <- effect - z_value * se
  upper <- effect + z_value * se
  
  rate <- exp(effect)  # Convert to incidence rate
  lower <- exp(lower)  # Lower bound
  upper <- exp(upper)  # Upper bound
  
  table_country_Ab_nonARI_glm <- table_country_Ab_nonARI_glm %>%   # Add results to the summary table
    add_row(CountryYear = country, Rate = rate, Lower = lower, Upper = upper, SE = se, numerator = n_Ab, denominator = n_Fever_or_LRI)
}
AbFever.nonARI <- table_country_Ab_nonARI_glm
colnames(AbFever.nonARI) <- c("CountryYear", "AbFever.nonARI", "AbFever.nonARI_L", "AbFever.nonARI_U", "AbFever.nonARI_SE", "n_Ab_Fever.nonARI", "n_Fever.nonARI")



# P(Antibiotics | ARI without Fever)
for (country in country_list) {  # Loop through each CountryYear
  country_data <- analysis_data_crude_nonFever %>% filter(CountryYear == country)  # Filter data for the current country-year
  
  model <- tryCatch(   # Fit the Poisson regression model
    glm(Antibiotics ~ ARI, 
        family = poisson(link = "log"), offset = log(PersonTime), data = country_data),
    error = function(e) {
      cat("Error fitting model for CountryYear:", country, "-", e$message, "\n")
      return(NULL)
    }
  )
  robust_se <- vcovHC(model, type = "HC0")
  
  n_LRI <- nrow(country_data[country_data$ARI == 3,])
  n_Ab <- nrow(country_data[country_data$ARI == 3 & country_data$Antibiotics == 1,])
  
  if (is.null(model)) next  # Skip if the model failed
  model_robust <- coeftest(model, vcov. = robust_se)
  effect <- model_robust[1, "Estimate"]
  se <- model_robust[1, "Std. Error"]
  
  z_value <- 1.96  # 95% confidence level
  lower <- effect - z_value * se
  upper <- effect + z_value * se
  
  rate <- exp(effect)  # Convert to incidence rate
  lower <- exp(lower)  # Lower bound
  upper <- exp(upper)  # Upper bound
  
  table_country_Ab_nonFever_glm <- table_country_Ab_nonFever_glm %>%   # Add results to the summary table
    add_row(CountryYear = country, Rate = rate, Lower = lower, Upper = upper, SE = se, numerator = n_Ab, denominator = n_LRI)
}
AbLRI.nonFever <- table_country_Ab_nonFever_glm
colnames(AbLRI.nonFever) <- c("CountryYear", "AbLRI.nonFever", "AbLRI.nonFever_L", "AbLRI.nonFever_U", "AbLRI.nonFever_SE", "n_Ab_ARI.nonFever", "n_ARI.nonFever")



# P(Antibiotics | Both ARI and Fever)
for (country in country_list) {  # Loop through each CountryYear
  country_data <- analysis_data_crude_Fever %>% filter(CountryYear == country)   # Filter data for the current country-year
  
  model <- tryCatch(   # Fit the Poisson regression model
    glm(Antibiotics ~ ARI, 
        family = poisson(link = "log"), offset = log(PersonTime), data = country_data),
    error = function(e) {
      cat("Error fitting model for CountryYear:", country, "-", e$message, "\n")
      return(NULL)
    }
  )
  robust_se <- vcovHC(model, type = "HC0")
  
  n_LRI <- nrow(country_data[country_data$ARI == 3,])
  n_Ab <- nrow(country_data[country_data$ARI == 3 & country_data$Antibiotics == 1,])
  
  if (is.null(model)) next   # Skip if the model failed
  model_robust <- coeftest(model, vcov. = robust_se)
  effect <- model_robust[1, "Estimate"]
  se <- model_robust[1, "Std. Error"]
  
  z_value <- 1.96  # 95% confidence level
  lower <- effect - z_value * se
  upper <- effect + z_value * se
  
  rate <- exp(effect)  # Convert to incidence rate
  lower <- exp(lower)  # Lower bound
  upper <- exp(upper)  # Upper bound
  
  table_country_Ab_Both_glm <- table_country_Ab_Both_glm %>%   # Add results to the summary table
    add_row(CountryYear = country, Rate = rate, Lower = lower, Upper = upper,SE = se, numerator = n_Ab, denominator = n_LRI)
}

AbBoth<- table_country_Ab_Both_glm
colnames(AbBoth) <- c("CountryYear", "AbBoth", "AbBoth_L", "AbBoth_U", "AbBoth_SE", "n_Ab_ARI_and_Fever", "n_ARI_and_Fever")

# Check the probability of Ab. by ARI w/o Fever, Fever w/o ARI, and Both : assign same P(Antibiotics) to ARI with and without Fever
Ab_prob <- AbLRI %>% 
  inner_join(AbFever, by = "CountryYear")  %>% 
  inner_join(AbFever.nonARI, by = "CountryYear") %>% 
  inner_join(AbLRI.nonFever, by = "CountryYear") %>% 
  inner_join(AbBoth, by = "CountryYear")

par(mfrow = c(1,2))
plot(Ab_prob$AbLRI, Ab_prob$AbBoth, xlim = c(0,1), ylim = c(0,1)); abline(a = 0, b = 1, col = "red", lty = 1)
plot(Ab_prob$AbFever, Ab_prob$AbBoth, xlim = c(0,1), ylim = c(0,1)); abline(a = 0, b = 1, col = "red", lty = 1)
cor(Ab_prob$AbLRI, Ab_prob$AbBoth, use = "complete.obs")  #0.982
cor(Ab_prob$AbFever, Ab_prob$AbBoth, use = "complete.obs") #0.951

plot(Ab_prob$AbLRI.nonFever, Ab_prob$AbBoth, xlim = c(0,1), ylim = c(0,1)); abline(a = 0, b = 1, col = "red", lty = 1)
plot(Ab_prob$AbFever.nonARI, Ab_prob$AbBoth, xlim = c(0,1), ylim = c(0,1)); abline(a = 0, b = 1, col = "red", lty = 1)
cor(Ab_prob$AbLRI.nonFever, Ab_prob$AbBoth, use = "complete.obs") 
cor(Ab_prob$AbFever.nonARI, Ab_prob$AbBoth, use = "complete.obs") 
par(mfrow = c(1,1))

analysis_data <- analysis_data %>% left_join(AbLRI[,c(1,2,5,6,7)], by = "CountryYear") # include P(Antibiotics | ARI with or without Fever)
analysis_data <- analysis_data %>% left_join(AbFever.nonARI[,c(1,2,5,6,7)], by = "CountryYear")  # include P(Antibiotics | Fever without ARI)

par(mfrow = c(4,1))
hist(unique(AbFever$AbFever), main ='P(Ab treatment|LRI or Fever)', xlab = '', ylab='', xlim=c(0,1))
hist(unique(AbFever.nonARI$AbFever.nonARI), main ='P(Ab treatment|Fever without LRI)', xlab = '', ylab='', xlim=c(0,1))
hist(unique(AbLRI.nonFever$AbLRI.nonFever), main ='P(Ab treatment|LRI without Fever)', xlab = '', ylab='', xlim=c(0,1))
hist(unique(AbBoth$AbBoth), main ='P(Ab treatment|Fever and LRI)', xlab = '', ylab='', xlim=c(0,1))


# -----------------Distribution of ARI and Fever------------------

# P(Fever | ARI)
analysis_data_crude_2 <- analysis_data_crude_2 %>% mutate(Fever = as.numeric(Fever))

for (country in country_list) {  # Loop through each CountryYear
  country_data <- analysis_data_crude_2 %>% filter(CountryYear == country)    # Filter data for the current country-year
  model <- tryCatch(    # Fit the Poisson regression model
    glm(Fever ~ ARI, 
        family = poisson(link = "log"), offset = log(PersonTime), data = country_data),
    error = function(e) {
      cat("Error fitting model for CountryYear:", country, "-", e$message, "\n")
      return(NULL)
    }
  )
  robust_se <- vcovHC(model, type = "HC0")

  n_LRI <- nrow(country_data[country_data$ARI == 3, ])
  n_Fever_LRI <- nrow(country_data[country_data$ARI == 3 & country_data$Fever == 1, ]) 
  
  if (is.null(model)) next    # Skip if the model failed
  model_robust <- coeftest(model, vcov. = robust_se)
  effect <- model_robust[1, "Estimate"]
  se <- model_robust[1, "Std. Error"]
  
  z_value <- 1.96  # 95% confidence level
  lower <- effect - z_value * se
  upper <- effect + z_value * se
  
  rate <- exp(effect)  # Convert to incidence rate
  lower <- exp(lower)  # Lower bound
  upper <- exp(upper)  # Upper bound
  
  table_country_Fever_ARI_glm <- table_country_Fever_ARI_glm %>%   # Add results to the summary table
    add_row(CountryYear = country, Rate = rate, Lower = lower, Upper = upper, SE = se, numerator = n_Fever_LRI, denominator = n_LRI)
}

FeverARI <- table_country_Fever_ARI_glm
colnames(FeverARI) <- c("CountryYear", "Fever_LRI", "Fever_LRI_L", "Fever_LRI_U", "Fever_LRI_SE", "n_Fever_LRI", "n_LRI")
analysis_data <- analysis_data %>% left_join(FeverARI[,c(1,2,5,6)], by = "CountryYear")



# P(Fever | RapidBreath)
analysis_data_crude_2 <- analysis_data_crude_2 %>% mutate(ARI = factor(ARI, levels = c(2, 3, 0, 1)))

for (country in country_list) {  # Loop through each CountryYear
  country_data <- analysis_data_crude_2 %>% filter(CountryYear == country)    # Filter data for the current country-year

  model <- tryCatch(   # Fit the Poisson regression model
    glm(Fever ~ ARI, 
        family = poisson(link = "log"), offset = log(PersonTime), data = country_data),
    error = function(e) {
      cat("Error fitting model for CountryYear:", country, "-", e$message, "\n")
      return(NULL)
    }
  )
  robust_se <- vcovHC(model, type = "HC0")
  
  n_Fever_Rapid <- nrow(country_data[country_data$ARI == 2 & country_data$Fever == 1, ]) 
  n_Rapid <- nrow(country_data[country_data$ARI == 2, ])
  
  if (is.null(model)) next    # Skip if the model failed
  model_robust <- coeftest(model, vcov. = robust_se)
  effect <- model_robust[1, "Estimate"]
  se <- model_robust[1, "Std. Error"]
  
  z_value <- 1.96  # 95% confidence level
  lower <- effect - z_value * se
  upper <- effect + z_value * se
  
  rate <- exp(effect)  # Convert to incidence rate
  lower <- exp(lower)  # Lower bound
  upper <- exp(upper)  # Upper bound
  
  table_country_Fever_RapidBreath_glm <- table_country_Fever_RapidBreath_glm %>%  # Add results to the summary table
    add_row(CountryYear = country, Rate = rate, Lower = lower, Upper = upper, SE = se, numerator = n_Fever_Rapid, denominator = n_Rapid)
}

Fever_RapidBreath <- table_country_Fever_RapidBreath_glm
colnames(Fever_RapidBreath) <- c("CountryYear", "Fever_RapidBreath", "Fever_RapidBreath_L", "Fever_RapidBreath_U", "Fever_RapidBreath_SE", "n_Fever_RapidBreath", "n_RapidBreath")

analysis_data <- analysis_data %>% left_join(Fever_RapidBreath[,c(1,2,5,6,7)], by = "CountryYear")



# P(Fever | Cough)
analysis_data_crude_2_nonARI <- analysis_data_crude_2[analysis_data_crude_2$ARI!=3,]

for (country in country_list) {   # Loop through each CountryYear
  country_data <- analysis_data_crude_2_nonARI %>% filter(CountryYear == country)    # Filter data for the current country-year

  model <- tryCatch(   # Fit the Poisson regression model
    glm(Fever ~ Cough, 
        family = poisson(link = "log"), offset = log(PersonTime), data = country_data),
    error = function(e) {
      cat("Error fitting model for CountryYear:", country, "-", e$message, "\n")
      return(NULL)
    }
  )
  robust_se <- vcovHC(model, type = "HC0")
  
  n_Fever_Cough <- nrow(country_data[country_data$Cough == 1 & country_data$Fever == 1, ]) 
  n_Cough <- nrow(country_data[country_data$Cough == 1, ]) 

  if (is.null(model)) next   # Skip if the model failed

  model_robust <- coeftest(model, vcov. = robust_se)
  effect <- model_robust[1, "Estimate"]
  se <- model_robust[1, "Std. Error"]
  
  z_value <- 1.96  # 95% confidence level
  lower <- effect - z_value * se
  upper <- effect + z_value * se
  
  rate <- exp(effect)  # Convert to incidence rate
  lower <- exp(lower)  # Lower bound
  upper <- exp(upper)  # Upper bound
  
  table_country_Fever_Cough_glm <- table_country_Fever_Cough_glm %>%   # Add results to the summary table
    add_row(CountryYear = country, Rate = rate, Lower = lower, Upper = upper, SE = se, numerator = n_Fever_Cough, denominator = n_Cough)
}

Fever_Cough <- table_country_Fever_Cough_glm
colnames(Fever_Cough) <- c("CountryYear", "Fever_Cough", "Fever_Cough_L", "Fever_Cough_U", "Fever_Cough_SE", "n_Fever_Cough", "n_Cough")
analysis_data <- analysis_data %>% left_join(Fever_Cough[,c(1,2,5,6,7)], by = "CountryYear")



# assign lower value between P(Fever|Rapid breath) or P(Fever|Cough)
FeverARI$'Fever_NonARI'<- NA
for (i in 1:nrow(FeverARI)){
  FeverARI$'Fever_NonARI'[i] <-  min(Fever_RapidBreath$`Fever_RapidBreath`[i], Fever_Cough$`Fever_Cough`[i])
}

plot(Fever_RapidBreath$Fever_RapidBreath,Fever_Cough$Fever_Cough, xlim = c(0,1), ylim = c(0,1)); abline(a = 0, b = 1, col = "red", lty = 1)
par(mfrow = c(3,1))
plot(FeverARI$Fever_LRI,Fever_Cough$Fever_Cough, xlim = c(0,1), ylim = c(0,1)); abline(a = 0, b = 1, col = "red", lty = 1)
plot(FeverARI$Fever_LRI,Fever_RapidBreath$Fever_RapidBreath, xlim = c(0,1), ylim = c(0,1)); abline(a = 0, b = 1, col = "red", lty = 1)
plot(FeverARI$Fever_LRI,FeverARI$Fever_NonARI, xlim = c(0,1), ylim = c(0,1)); abline(a = 0, b = 1, col = "red", lty = 1)
par(mfrow = c(1,1))

analysis_data <- analysis_data %>% mutate(
  Fever_URI = pmin(`Fever_RapidBreath`, `Fever_Cough`, na.rm = TRUE),
  Fever_URI_SE = if_else(`Fever_RapidBreath` <= `Fever_Cough` | is.na(`Fever_Cough`), 
                         `Fever_RapidBreath_SE`, `Fever_Cough_SE`),
  n_Fever_URI =  if_else(`Fever_RapidBreath` <= `Fever_Cough` | is.na(`Fever_Cough`), 
                         `n_Fever_RapidBreath`, `n_Fever_Cough`),
  n_URI =        if_else(`Fever_RapidBreath` <= `Fever_Cough` | is.na(`Fever_Cough`), 
                         `n_RapidBreath`, `n_Cough`)
)

plot_bar_hist <- function(x, main_title){
  h <- hist(x, breaks = breaks_seq, plot = FALSE)
  
  # Create barplot without axes
  bound <- barplot(h$counts, space = 0,
                  col = rgb(0/255,114/255,178/255,0.7), border = rgb(0/255,114/255,178/255),
                  xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = "")
  
  # Custom x-axis (0 to 1)
  axis(1, at = seq(0, 10, by = 1), labels = round(breaks_seq,2), cex.axis = 0.8)
  
  # Custom y-axis at x = 0
  y_ticks <- pretty(c(0, max(h$counts)))
  axis(2, at = y_ticks, labels = y_ticks, las = 2, cex.axis = 0.8, pos = 0)
  
  # Draw x = 0 line for y-axis
  abline(v = 0, col = "black", lwd = 1)
  
  # Title
  title(main = main_title, cex.main = 1.1, line = 1.5)
  mtext("Count", side = 2, line = 1.0, cex = 0.7)
}



breaks_seq <- seq(0, 1, by = 0.1)


# ============================================================
# Export source data for Figure S3
# ============================================================

dir.create(
  "outputs/plot_objects",
  recursive = TRUE,
  showWarnings = FALSE
)

make_hist_source <- function(x, panel, measure) {
  
  x <- unique(x)
  x <- x[is.finite(x)]
  
  h <- hist( x, breaks = breaks_seq, plot = FALSE)
  
  data.frame(
    Panel = panel,
    Measure = measure,
    Bin_lower = head(h$breaks, -1),
    Bin_upper = tail(h$breaks, -1),
    Bin_midpoint = h$mids,
    Count = h$counts
  )
}


figS3_A <- make_hist_source(
  analysis_data$Fever_LRI,
  panel = "A",
  measure = "P(Fever|LRTI)"
)

figS3_B <- make_hist_source(
  Fever_RapidBreath$Fever_RapidBreath,
  panel = "B",
  measure = "P(Fever|Rapid Breath)"
)

figS3_C <- make_hist_source(
  Fever_Cough$Fever_Cough,
  panel = "C",
  measure = "P(Fever|Cough)"
)

figS3_D <- make_hist_source(
  analysis_data$Fever_URI,
  panel = "D",
  measure = "Minimum of P(Fever|Rapid Breath) and P(Fever|Cough)"
)

Figure_S3_source_data <- bind_rows(
  figS3_A,
  figS3_B,
  figS3_C,
  figS3_D
)

write.csv(Figure_S3_source_data, "outputs/plot_objects/Figure_S2_source_data.csv", row.names = FALSE)

dir.create(
  "figures",
  recursive = TRUE,
  showWarnings = FALSE
)

pdf("figures/Figure_S2.pdf", width = 7.5, height = 5.5)
par(mfrow = c(4,1))
par(mar = c(4, 2.5, 2.5, 1))  # Bottom, Left, Top, Right
par(mgp = c(3, 0.3, 0))     # Axis label and tick distance
par(tck = -0.025)
par(lwd = 1)

# Plot all 4 histograms
plot_bar_hist(unique(analysis_data$Fever_LRI), "P(Fever|LRTI)")
plot_bar_hist(unique(Fever_RapidBreath$Fever_RapidBreath), "P(Fever|Rapid Breath)")
plot_bar_hist(unique(Fever_Cough$Fever_Cough), "P(Fever|Cough)")
plot_bar_hist(unique(analysis_data$Fever_URI), "Minimum of P(Fever|Rapid Breath) and P(Fever|Cough)")
mtext("Proportion", side = 1, line = 2.5, cex = 0.7)

dev.off()

saveRDS(analysis_data, file = "outputs/analysis_data_prop.rds")
rm(list = ls())
gc()