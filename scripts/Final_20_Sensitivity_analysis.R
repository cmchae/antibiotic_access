########## 1. Probability ##########

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

## =============crude model with robust SE: glm (by country)===============

analysis_data <- readRDS("outputs/analysis_data.rds")

# sensitivity analysis: exclude children treated with antimalarial agents
n_before <- nrow(analysis_data)
n_antimalarial_excluded <- sum(analysis_data$Antimalaria == 1, na.rm = TRUE)

analysis_data <- analysis_data %>%
  filter(is.na(Antimalaria) | Antimalaria != 1)

c(
  n_before = n_before,
  n_antimalarial_excluded = n_antimalarial_excluded,
  n_after = nrow(analysis_data)
)

stopifnot(!any(analysis_data$Antimalaria == 1, na.rm = TRUE))

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

country_list <- unique(analysis_data_crude$CountryYear)  

table_country_Ab_ARI_glm <- 
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
for (country in country_list) {  
  country_data <- analysis_data_crude %>% filter(CountryYear == country)   
  
  model <- tryCatch(   
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
  
  if (is.null(model)) next  
  model_robust <- coeftest(model, vcov. = robust_se)
  effect <- model_robust[1, "Estimate"]
  se <- model_robust[1, "Std. Error"]
  
  z_value <- 1.96  
  lower <- effect - z_value * se
  upper <- effect + z_value * se
  
  rate <- exp(effect)  
  lower <- exp(lower)  
  upper <- exp(upper)  
  
  
  table_country_Ab_ARI_glm <- table_country_Ab_ARI_glm %>%   
    add_row(CountryYear = country, Rate = rate, Lower = lower, Upper = upper, SE = se, numerator = n_Ab, denominator = n_LRI)
}
AbLRI <- table_country_Ab_ARI_glm
colnames(AbLRI) <- c("CountryYear", "AbLRI", "AbLRI_L", "AbLRI_U", "AbLRI_SE", "n_AbLRI", "n_LRI")



# P(Antibiotics | Fever or ARI)
for (country in country_list) {  
  country_data <- analysis_data_crude %>% filter(CountryYear == country)   
  
  model <- tryCatch(   
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
  
  if (is.null(model)) next  
  model_robust <- coeftest(model, vcov. = robust_se)
  effect <- model_robust[1, "Estimate"]
  se <- model_robust[1, "Std. Error"]
  
  z_value <- 1.96  
  lower <- effect - z_value * se
  upper <- effect + z_value * se
  
  rate <- exp(effect)  
  lower <- exp(lower)  
  upper <- exp(upper)  
  
  table_country_Ab_Fever_glm <- table_country_Ab_Fever_glm %>%   
    add_row(CountryYear = country, Rate = rate, Lower = lower, Upper = upper, SE = se, numerator = n_Ab, denominator = n_Fever_or_LRI)
}
AbFever<- table_country_Ab_Fever_glm
colnames(AbFever) <- c("CountryYear", "AbFever", "AbFever_L", "AbFever_U", "AbFever_SE", "n_Ab_Fever_or_ARI", "n_Fever_or_ARI")



# P(Antibiotics | Fever without ARI)
for (country in country_list) {  
  country_data <- analysis_data_crude_nonARI %>% filter(CountryYear == country)  
  
  model <- tryCatch(   
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
  
  if (is.null(model)) next  
  model_robust <- coeftest(model, vcov. = robust_se)
  effect <- model_robust[1, "Estimate"]
  se <- model_robust[1, "Std. Error"]
  
  z_value <- 1.96  
  lower <- effect - z_value * se
  upper <- effect + z_value * se
  
  rate <- exp(effect)  
  lower <- exp(lower)  
  upper <- exp(upper)  
  
  table_country_Ab_nonARI_glm <- table_country_Ab_nonARI_glm %>%   
    add_row(CountryYear = country, Rate = rate, Lower = lower, Upper = upper, SE = se, numerator = n_Ab, denominator = n_Fever_or_LRI)
}
AbFever.nonARI <- table_country_Ab_nonARI_glm
colnames(AbFever.nonARI) <- c("CountryYear", "AbFever.nonARI", "AbFever.nonARI_L", "AbFever.nonARI_U", "AbFever.nonARI_SE", "n_Ab_Fever.nonARI", "n_Fever.nonARI")



# P(Antibiotics | ARI without Fever)
for (country in country_list) {  
  country_data <- analysis_data_crude_nonFever %>% filter(CountryYear == country)  
  
  model <- tryCatch(   
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
  
  if (is.null(model)) next  
  model_robust <- coeftest(model, vcov. = robust_se)
  effect <- model_robust[1, "Estimate"]
  se <- model_robust[1, "Std. Error"]
  
  z_value <- 1.96  
  lower <- effect - z_value * se
  upper <- effect + z_value * se
  
  rate <- exp(effect)  
  lower <- exp(lower)  
  upper <- exp(upper)  
  
  table_country_Ab_nonFever_glm <- table_country_Ab_nonFever_glm %>%   
    add_row(CountryYear = country, Rate = rate, Lower = lower, Upper = upper, SE = se, numerator = n_Ab, denominator = n_LRI)
}
AbLRI.nonFever <- table_country_Ab_nonFever_glm
colnames(AbLRI.nonFever) <- c("CountryYear", "AbLRI.nonFever", "AbLRI.nonFever_L", "AbLRI.nonFever_U", "AbLRI.nonFever_SE", "n_Ab_ARI.nonFever", "n_ARI.nonFever")



# P(Antibiotics | Both ARI and Fever)
for (country in country_list) {  
  country_data <- analysis_data_crude_Fever %>% filter(CountryYear == country)   
  
  model <- tryCatch(   
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
  
  if (is.null(model)) next   
  model_robust <- coeftest(model, vcov. = robust_se)
  effect <- model_robust[1, "Estimate"]
  se <- model_robust[1, "Std. Error"]
  
  z_value <- 1.96  
  lower <- effect - z_value * se
  upper <- effect + z_value * se
  
  rate <- exp(effect)  
  lower <- exp(lower)  
  upper <- exp(upper)  
  
  table_country_Ab_Both_glm <- table_country_Ab_Both_glm %>%   
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
cor(Ab_prob$AbLRI, Ab_prob$AbBoth, use = "complete.obs")
cor(Ab_prob$AbFever, Ab_prob$AbBoth, use = "complete.obs")

plot(Ab_prob$AbLRI.nonFever, Ab_prob$AbBoth, xlim = c(0,1), ylim = c(0,1)); abline(a = 0, b = 1, col = "red", lty = 1)
plot(Ab_prob$AbFever.nonARI, Ab_prob$AbBoth, xlim = c(0,1), ylim = c(0,1)); abline(a = 0, b = 1, col = "red", lty = 1)
cor(Ab_prob$AbLRI.nonFever, Ab_prob$AbBoth, use = "complete.obs")
cor(Ab_prob$AbFever.nonARI, Ab_prob$AbBoth, use = "complete.obs")
par(mfrow = c(1,1))

analysis_data <- analysis_data %>% left_join(AbLRI[,c(1,2,5,6,7)], by = "CountryYear") 
analysis_data <- analysis_data %>% left_join(AbFever.nonARI[,c(1,2,5,6,7)], by = "CountryYear")  

par(mfrow = c(4,1))
hist(unique(AbFever$AbFever), main ='P(Ab treatment|LRI or Fever)', xlab = '', ylab='', xlim=c(0,1))
hist(unique(AbFever.nonARI$AbFever.nonARI), main ='P(Ab treatment|Fever without LRI)', xlab = '', ylab='', xlim=c(0,1))
hist(unique(AbLRI.nonFever$AbLRI.nonFever), main ='P(Ab treatment|LRI without Fever)', xlab = '', ylab='', xlim=c(0,1))
hist(unique(AbBoth$AbBoth), main ='P(Ab treatment|Fever and LRI)', xlab = '', ylab='', xlim=c(0,1))



# -----------------Distribution of ARI and Fever------------------

# P(Fever | ARI)
analysis_data_crude_2 <- analysis_data_crude_2 %>% mutate(Fever = as.numeric(Fever))

for (country in country_list) {  
  country_data <- analysis_data_crude_2 %>% filter(CountryYear == country)    
  model <- tryCatch(    
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
  
  if (is.null(model)) next    
  model_robust <- coeftest(model, vcov. = robust_se)
  effect <- model_robust[1, "Estimate"]
  se <- model_robust[1, "Std. Error"]
  
  z_value <- 1.96  
  lower <- effect - z_value * se
  upper <- effect + z_value * se
  
  rate <- exp(effect)  
  lower <- exp(lower)  
  upper <- exp(upper)  
  
  table_country_Fever_ARI_glm <- table_country_Fever_ARI_glm %>%   
    add_row(CountryYear = country, Rate = rate, Lower = lower, Upper = upper, SE = se, numerator = n_Fever_LRI, denominator = n_LRI)
}

FeverARI <- table_country_Fever_ARI_glm
colnames(FeverARI) <- c("CountryYear", "Fever_LRI", "Fever_LRI_L", "Fever_LRI_U", "Fever_LRI_SE", "n_Fever_LRI", "n_LRI")
analysis_data <- analysis_data %>% left_join(FeverARI[,c(1,2,5,6)], by = "CountryYear")



# P(Fever | RapidBreath)
analysis_data_crude_2 <- analysis_data_crude_2 %>% mutate(ARI = factor(ARI, levels = c(2, 3, 0, 1)))

for (country in country_list) {  
  country_data <- analysis_data_crude_2 %>% filter(CountryYear == country)    
  
  model <- tryCatch(   
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
  
  if (is.null(model)) next    
  model_robust <- coeftest(model, vcov. = robust_se)
  effect <- model_robust[1, "Estimate"]
  se <- model_robust[1, "Std. Error"]
  
  z_value <- 1.96  
  lower <- effect - z_value * se
  upper <- effect + z_value * se
  
  rate <- exp(effect)  
  lower <- exp(lower)  
  upper <- exp(upper)  
  
  table_country_Fever_RapidBreath_glm <- table_country_Fever_RapidBreath_glm %>%  
    add_row(CountryYear = country, Rate = rate, Lower = lower, Upper = upper, SE = se, numerator = n_Fever_Rapid, denominator = n_Rapid)
}

Fever_RapidBreath <- table_country_Fever_RapidBreath_glm
colnames(Fever_RapidBreath) <- c("CountryYear", "Fever_RapidBreath", "Fever_RapidBreath_L", "Fever_RapidBreath_U", "Fever_RapidBreath_SE", "n_Fever_RapidBreath", "n_RapidBreath")

analysis_data <- analysis_data %>% left_join(Fever_RapidBreath[,c(1,2,5,6,7)], by = "CountryYear")



# P(Fever | Cough)
analysis_data_crude_2_nonARI <- analysis_data_crude_2[analysis_data_crude_2$ARI!=3,]

for (country in country_list) {   
  country_data <- analysis_data_crude_2_nonARI %>% filter(CountryYear == country)    
  
  model <- tryCatch(   
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
  
  if (is.null(model)) next   
  
  model_robust <- coeftest(model, vcov. = robust_se)
  effect <- model_robust[1, "Estimate"]
  se <- model_robust[1, "Std. Error"]
  
  z_value <- 1.96  
  lower <- effect - z_value * se
  upper <- effect + z_value * se
  
  rate <- exp(effect)  
  lower <- exp(lower)  
  upper <- exp(upper)  
  
  table_country_Fever_Cough_glm <- table_country_Fever_Cough_glm %>%   
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

length(unique(analysis_data$Country))
unique(analysis_data$SurveyYear)

dir.create(
  "outputs/sensitivity_analysis",
  recursive = TRUE,
  showWarnings = FALSE
)

### save the cleaned version of data frame
saveRDS(analysis_data, file = "outputs/sensitivity_analysis/2-2.SA_analysis_data_prop.rds")
rm(list = ls())
gc()
























########## 2. Fever_given_RI ##########
rm(list = ls())
gc()

options(scipen = 999, digits = 4)

requiredPackages = c("dplyr","tidyr","readr","betareg", "MASS")
lib.path <- .libPaths()
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p, lib = lib.path)
  library(p,character.only = TRUE)
}
set.seed(123)
n_sim <- 10000

analysis_data <- readRDS("outputs/sensitivity_analysis/2-2.SA_analysis_data_prop.rds")


UNI <- readRDS("outputs/pc_gdp_ten.rds")

prop.dat <- analysis_data %>% dplyr::select(c(Country, Year, Fever_LRI, Fever_URI, Fever_LRI_SE, Fever_URI_SE, n_Fever_LRI, n_Fever_URI, n_LRI, n_URI))
a <- unique(prop.dat$Country); b <- unique(UNI$Country); 
setdiff(a,b); setdiff(b,a)

prop.data <- prop.dat %>% inner_join(UNI, by = c("Country","Year"))
prop.data <- prop.data %>% distinct() %>% na.omit()

x <- unique(distinct(analysis_data[,c("Country","Year")]))
y <- unique(distinct(prop.data[,c("Country","Year")]))
setdiff(x,y) # Yemen 2013 was excluded as there isn't PC aligned to this country-year (only in 2020)

UNI.data <- UNI[,1:3]



# Malaria proportion :IHME data (GBD2021, adjust the proportion of malaria in 2021 to 2022-2024) - https://vizhub.healthdata.org/gbd-results/
mal.prop <- read.csv("data/malaria_proportion.csv")
mal.prop <- mal.prop[,-1]
colnames(mal.prop) <- c("Country", "Age", "Year", "mal.prop", "RI.prop")

a <- unique(UNI.data$Country); b <- unique(prop.data$Country); c <- unique(mal.prop$Country); 
setdiff(a,c); setdiff(c,a)
setdiff(b,c); setdiff(c,b)

country_mapping_mal <- c("Taiwan (Province of China)" = "Taiwan",   "Lao People's Democratic Republic" = "Lao PDR", "Kyrgyzstan" = "Kyrgyz Republic", "Democratic People's Republic of Korea" = "Korea (Democratic People's Rep. of)", "Micronesia (Federated States of)" = "Micronesia, Fed. Sts.", "Republic of Moldova" = "Moldova", "Republic of Korea" = "Korea, Rep.", "Slovakia" = "Slovak Republic", "Bahamas" = "Bahamas, The", "Venezuela (Bolivarian Republic of)" = "Venezuela, RB", "United States of America" = "United States", "Saint Vincent and the Grenadines" = "St. Vincent and the Grenadines", "Bolivia (Plurinational State of)" = "Bolivia", "Saint Lucia" = "St. Lucia", "Iran (Islamic Republic of)" = "Iran, Islamic Rep.", "Palestine" = "West Bank and Gaza", "Egypt" = "Egypt, Arab Rep.", "Türkiye" = "Turkiye", "Yemen" = "Yemen, Rep.", "Congo" = "Congo, Rep.", "Democratic Republic of the Congo" = "Congo, Dem. Rep.", "United Republic of Tanzania" = "Tanzania", "Côte d'Ivoire" = "Cote d'Ivoire", "Gambia" = "Gambia, The", "Saint Kitts and Nevis" = "St. Kitts and Nevis", "United States Virgin Islands" = "Virgin Islands (U.S.)"
)

mal.prop <- mal.prop %>% mutate(Country = ifelse(Country %in% names(country_mapping_mal), country_mapping_mal[Country], Country))

mal.prop_under5 <- mal.prop[mal.prop$Age == "<5 years", ]
mal.prop_all.ages <- mal.prop[mal.prop$Age == "All ages", ]

prop.data <- prop.data %>% inner_join(mal.prop_all.ages, by = c("Country" = "Country", "Year" = "Year"))
UNI.data <- UNI.data %>% inner_join(mal.prop_all.ages, by = c("Country" = "Country", "Year" = "Year"))



#Temperature: World bank - Climate Change Knowledge Portal, https://climateknowledgeportal.worldbank.org/
Temperature <- read.csv("data/Temperature_2024.csv")
Temperature <- Temperature %>% pivot_longer(cols = starts_with("X"), names_to = "year_month", values_to = "temperature")  
Temperature <- Temperature %>%mutate(year_month = sub("^X", "", year_month)) %>%  separate(year_month, into = c("Year", "Month"), sep = "\\.") 
Temperature <- Temperature %>% mutate(Year = as.integer(Year), Month = as.integer(Month)) 
Temperature <- Temperature %>% group_by(name, Year, Month) %>% summarise(temperature = mean(temperature), .groups = "drop") 
Temperature <- Temperature %>% pivot_wider(names_from = Month, names_prefix = "Month", values_from = temperature) 
Temperature <- Temperature %>%  arrange(name, Year) %>% rename(Country = name) 
Temperature$Temperature <- apply(Temperature[, 3:14], 1, mean) 
Temperature$Temp_maximum <- apply(Temperature[, 3:14], 1, max) 
plot(Temperature$Temperature, Temperature$Temp_maximum)
Temperature <- Temperature[, c("Country", "Year", "Temperature")]

a <- unique(UNI.data$Country); b <- unique(prop.data$Country); c <- unique(Temperature$Country); 
setdiff(a,c); setdiff(c,a)
setdiff(b,c); setdiff(c,b)

country_mapping_Temp <- c("Aruba (Neth.)" = "Aruba", "The Bahamas" = "Bahamas, The", "Bermuda (U.K.)" = "Bermuda", "Cayman Islands (U.K.)" = "Cayman Islands", "Democratic Republic of Congo" = "Congo, Dem. Rep.", "Congo" = "Congo, Rep.", "C\x99te d'Ivoire" = "Cote d'Ivoire", "Cura\x8dao (Neth.)" = "Curacao", "Czech Republic" = "Czechia", "Arab Republic of Egypt" = "Egypt, Arab Rep.", "Faroe Islands (Den.)" = "Faroe Islands", "The Gambia" = "Gambia, The", "Greenland (Den.)" = "Greenland", "Hong Kong, SAR" = "Hong Kong SAR, China", "Islamic Republic of Iran" = "Iran, Islamic Rep.", "Republic of Korea" = "Korea, Rep.", "Lao People's Democratic Republic" = "Lao PDR", "Macau, SAR" = "Macao SAR, China", "Federated States of Micronesia" = "Micronesia, Fed. Sts.", "Puerto Rico (U.S.)" = "Puerto Rico", "S\x8bo Tom\x8e and Pr\x92ncipe" = "Sao Tome and Principe", "Sint Maarten (Neth.)" = "Sint Maarten (Dutch part)", "Saint Kitts and Nevis" = "St. Kitts and Nevis", "Saint Lucia" = "St. Lucia", "Saint Vincent and the Grenadines" = "St. Vincent and the Grenadines", "T\x9frkiye" = "Turkiye", "Turks and Caicos Islands (U.K.)" = "Turks and Caicos Islands", "United States of America" = "United States", "R. B. de Venezuela" = "Venezuela, RB", "Vietnam" = "Viet Nam", "United States Virgin Islands (U.S.)" = "Virgin Islands (U.S.)", "Republic of Yemen" = "Yemen, Rep."
)

Temperature <- Temperature %>% mutate(Country = ifelse(Country %in% names(country_mapping_Temp), country_mapping_Temp[Country], Country))
Temperature_filtered <- Temperature %>% filter(Country %in% b)
prop.data <- prop.data %>% inner_join(Temperature, by = c("Country" = "Country", "Year" = "Year"))
UNI.data <- UNI.data %>% inner_join(Temperature, by = c("Country" = "Country", "Year" = "Year"))



#Precipitation: Our world in data - Annual precipitation, 1940 to 2024, https://ourworldindata.org/grapher/average-precipitation-per-year?tab=table
Precipitation <- read.csv("data/Precipitation.csv")
Precipitation <- as.data.frame(Precipitation)
Precipitation$Year <- as.numeric(Precipitation$Year) 

a <- unique(UNI.data$Country); b <- unique(prop.data$Country); c <- unique(Precipitation$Country); 
setdiff(a,c); setdiff(c,a)
setdiff(b,c); setdiff(c,b)

country_mapping_Prec <- c("Bahamas" = "Bahamas, The", "Brunei" = "Brunei Darussalam", "Cape Verde" = "Cabo Verde", "Democratic Republic of Congo" = "Congo, Dem. Rep.", "Congo" = "Congo, Rep.", "Curacao" = "Curacao", "Dominica" = "Dominica", "Egypt" = "Egypt, Arab Rep.", "Gambia" = "Gambia, The", "Grenada" = "Grenada", "Hong Kong" = "Hong Kong SAR, China", "Iran" = "Iran, Islamic Rep.", "South Korea" = "Korea, Rep.", "Kyrgyzstan" = "Kyrgyz Republic", "Laos" = "Lao PDR", "Macau" = "Macao SAR, China", "Micronesia" = "Micronesia, Fed. Sts.", "Russia" = "Russian Federation", "Slovakia" = "Slovak Republic", "Saint Vincent and the Grenadines" = "St. Vincent and the Grenadines", "Saint Kitts and Nevis" = "St. Kitts and Nevis", "Saint Lucia" = "St. Lucia", "Syria" = "Syrian Arab Republic", "East Timor" = "Timor-Leste", "Turkey" = "Turkiye", "Turks and Caicos Islands" = "Turks and Caicos Islands", "United States Virgin Islands" = "Virgin Islands (U.S.)", "Venezuela" = "Venezuela, RB", "Vietnam" = "Viet Nam", "Yemen" = "Yemen, Rep.", "Palestine" = "West Bank and Gaza"
)

Precipitation <- Precipitation %>% mutate(Country = ifelse(Country %in% names(country_mapping_Prec), country_mapping_Prec[Country], Country))
prop.data <- prop.data %>% inner_join(Precipitation, by = c("Country" = "Country", "Year" = "Year"))
UNI.data <- UNI.data %>% inner_join(Precipitation, by = c("Country" = "Country", "Year" = "Year"))


### Fever/URI: Logistic binomial --------------
mod_sim_URI_logit <- glm(
  cbind(n_Fever_URI,n_URI-n_Fever_URI) ~ Year + GDP + Temperature + Precipitation + mal.prop,
  family = binomial(link="logit"), data = prop.data)

coef_URI <- coef(mod_sim_URI_logit); vcov_URI <- vcov(mod_sim_URI_logit)

new.data.UNI <- UNI.data %>% dplyr::select(c("Year", "GDP", "Temperature", "Precipitation", "mal.prop"))
X <- cbind(1, as.matrix(new.data.UNI))  
n_obs <- nrow(X)
beta_sim <- MASS::mvrnorm(n = n_sim, mu = coef_URI, Sigma = vcov_URI)  
logit_p_sim <- X %*% t(beta_sim)
pred.prop.URI <- plogis(logit_p_sim)

pred.prop.med_URI <- apply(pred.prop.URI, 1, median)
pred.prop.lower_URI <- apply(pred.prop.URI, 1, quantile, probs = 0.025)
pred.prop.upper_URI <- apply(pred.prop.URI, 1, quantile, probs = 0.975)

prop.med_URI <- round(pred.prop.med_URI, 3)
prop.lower_URI <- round(pred.prop.lower_URI, 3)
prop.upper_URI <- round(pred.prop.upper_URI, 3)


### Fever/LRI
mod_sim_LRI_logit <- glm(
  cbind(n_Fever_LRI,n_LRI-n_Fever_LRI) ~ Year + GDP + Temperature + Precipitation + mal.prop,
  family = binomial(link="logit"), data = prop.data)

coef_LRI <- coef(mod_sim_LRI_logit); vcov_LRI <- vcov(mod_sim_LRI_logit)

new.data.UNI <- UNI.data %>% dplyr::select(c("Year", "GDP", "Temperature", "Precipitation", "mal.prop"))
X <- cbind(1, as.matrix(new.data.UNI))  
n_obs <- nrow(X)
beta_sim <- MASS::mvrnorm(n = n_sim, mu = coef_LRI, Sigma = vcov_LRI)  
logit_p_sim <- X %*% t(beta_sim)
pred.prop.LRI <- plogis(logit_p_sim)

pred.prop.med_LRI <- apply(pred.prop.LRI, 1, median)
pred.prop.lower_LRI <- apply(pred.prop.LRI, 1, quantile, probs = 0.025)
pred.prop.upper_LRI <- apply(pred.prop.LRI, 1, quantile, probs = 0.975)

prop.med_LRI <- round(pred.prop.med_LRI, 3)
prop.lower_LRI <- round(pred.prop.lower_LRI, 3)
prop.upper_LRI <- round(pred.prop.upper_LRI, 3)



sum(residuals(mod_sim_LRI_logit, type = "pearson")^2) /
  df.residual(mod_sim_LRI_logit)
sum(residuals(mod_sim_URI_logit, type = "pearson")^2) /
  df.residual(mod_sim_URI_logit)



### Result
pred.prop_fever_RI <- data.frame(
  Country = UNI.data$Country,
  Year = UNI.data$Year,
  Median_LRI = prop.med_LRI,
  Lower_LRI = prop.lower_LRI,
  Upper_LRI = prop.upper_LRI,
  Median_URI = prop.med_URI,
  Lower_URI = prop.lower_URI,
  Upper_URI = prop.upper_URI
)

unique(pred.prop_fever_RI$Country)
unique(pred.prop_fever_RI$Year)

par(mfrow = c(2, 1))

hist(unique(analysis_data$Fever_LRI), xlim = c(0, 1), ylim = c(0, 15), breaks = seq(0,1,by=0.05), main = "P(Fever|LRI)", xlab = "")
hist(unique(analysis_data$Fever_URI), xlim = c(0, 1), ylim = c(0, 15), breaks = seq(0,1,by=0.05), main = "P(Fever|URI)", xlab = "")

hist(prop.med_URI, xlim = c(0, 1), breaks = seq(0,1,by=0.05), ylim = c(0, 1000), main = "P(Fever|URI)", xlab = "")
hist(prop.med_LRI, xlim = c(0, 1), breaks = seq(0,1,by=0.05), ylim = c(0, 1000), main = "P(Fever|LRI)", xlab = "")

par(mfrow = c(1, 1))

length(unique(pred.prop_fever_RI$Country))
unique(pred.prop_fever_RI$Year)


prop.Fever.RI <- data.frame(
  Country = UNI.data$Country,
  Year = UNI.data$Year,
  Fever_LRI = pred.prop.LRI,
  Fever_URI = pred.prop.URI)



saveRDS(pred.prop_fever_RI, "outputs/sensitivity_analysis/SA_Fever_RI_proportion_20yrs.rds")
saveRDS(prop.Fever.RI, "outputs/sensitivity_analysis/SA_prop.Fever.RI.rds")









########## 3. IHME ##########

rm(list = ls())
gc()

requiredPackages = c("dplyr", "tidyr", 'stringr', "pbapply", "parallel")
lib.path <- .libPaths()
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p, lib = lib.path)
  library(p,character.only = TRUE)
}

options(scipen = 100, digits = 4)

set.seed(123)
n_sim <- 10000

gbd.data <- read.csv("data/All_inf2.csv")
colnames(gbd.data) 
str(gbd.data)
unique(gbd.data$year) 
unique(gbd.data$location) 
unique(gbd.data$age) 
unique(gbd.data$cause)  

gbd.data <- gbd.data[gbd.data$cause == "Lower respiratory infections" | gbd.data$cause == "Upper respiratory infections",] 
gbd.data <- gbd.data[, c("location","year","cause","age","val","lower","upper")] 



# ------- fitting the LRI and URI point estimate and upper/lower limit to gamma distribution --------
estimate <- function(val, lower, upper) { 
  obj_fun <- function(par) {
    shape <- par[1]
    scale <- par[2]
    q_lower <- qgamma(0.025, shape = shape, scale = scale)
    q_upper <- qgamma(0.975, shape = shape, scale = scale)
    mean_val <- shape * scale
    (q_lower - lower)^2 + (q_upper - upper)^2 + (mean_val - val)^2
  }
  
  sd_est <- (upper - lower) / 4
  shape_init <- (val / sd_est)^2
  scale_init <- val / shape_init
  
  fit <- optim(
    par = c(shape_init, scale_init),
    fn = obj_fun,
    method = "L-BFGS-B",
    lower = c(1e-3, 1e-3)
  )
  
  shape <- fit$par[1]
  scale <- fit$par[2]
  
  fitted_mean <- shape * scale   
  fitted_lower <- qgamma(0.025, shape = shape, scale = scale)
  fitted_upper <- qgamma(0.975, shape = shape, scale = scale)
  
  error_val <- abs(fitted_mean - val)/val  
  error_lower <- abs(fitted_lower - lower)/lower
  error_upper <- abs(fitted_upper - upper)/upper
  
  return(c(shape = shape, scale = scale,
           error_val = error_val,
           error_lower = error_lower,
           error_upper = error_upper))
}

gbd.data <- gbd.data %>%
  rowwise() %>%
  mutate(
    gamma_est = list(estimate(val, lower, upper)),
    shape = gamma_est[[1]],
    scale = gamma_est[[2]],
    error_val = gamma_est[[3]],
    error_lower = gamma_est[[4]],
    error_upper = gamma_est[[5]]
  ) %>%
  dplyr::select(-gamma_est) %>%
  ungroup()



# Number of Respiratory infection (URI and LRI) among children ≤5 
par(mfrow = c(3,1))
hist(gbd.data$error_val)
hist(gbd.data$error_lower)
hist(gbd.data$error_upper)
par(mfrow = c(1,1))

gbd.data[which.max(gbd.data$error_val),]
gbd.data[which.max(gbd.data$error_lower),]
gbd.data[which.max(gbd.data$error_upper),]
gbd.data <- gbd.data %>% dplyr::select(location, year, cause, age, val, shape, scale)

gbd.data <- gbd.data %>% pivot_wider(names_from = age, values_from = c(shape, scale,val))
colnames(gbd.data) <- c("Country", "Year", "cause", "AllAge_shape", "U5_shape", "AllAge_scale", "U5_scale","AllAge_val", "U5_val")
gbd.data$U5_prop <- gbd.data$U5_val/gbd.data$AllAge_val
gbd.data <- gbd.data %>% dplyr::select(-c(AllAge_val,U5_val))
gbd.data <- gbd.data %>% pivot_wider(names_from = cause, values_from = c(AllAge_shape, U5_shape, AllAge_scale, U5_scale, U5_prop))
colnames(gbd.data) <- c("Country", "Year", "AllAge_shape_LRI", "AllAge_shape_URI", "U5_shape_LRI", "U5_shape_URI", "AllAge_scale_LRI", "AllAge_scale_URI", "U5_scale_LRI", "U5_scale_URI", "U5_prop_LRI", "U5_prop_URI")        
gbd.data <- gbd.data %>% dplyr::select(-c(U5_prop_URI, AllAge_shape_LRI, AllAge_shape_URI, AllAge_scale_LRI, AllAge_scale_URI))
gbd.data <- as.data.frame(gbd.data)

prop_fever_RI <- readRDS("outputs/sensitivity_analysis/SA_prop.Fever.RI.rds")
prop_RI_fever <- readRDS("outputs/prop.RI.Fever.RE.rds")

country_mapping <- c("Bahamas, The" = "Bahamas", "Bolivia" = "Bolivia (Plurinational State of)", "Congo, Dem. Rep." = "Democratic Republic of the Congo", "Congo, Rep." = "Congo", "Cote d'Ivoire" = "Côte d'Ivoire", "Egypt, Arab Rep." = "Egypt", "Gambia, The" = "Gambia", "Iran, Islamic Rep." = "Iran (Islamic Republic of)", "Korea, Rep." = "Republic of Korea", "Kyrgyz Republic" = "Kyrgyzstan", "Lao PDR" = "Lao People's Democratic Republic", "Moldova" = "Republic of Moldova", "Slovak Republic" = "Slovakia", "St. Vincent and the Grenadines" = "Saint Vincent and the Grenadines", "Tanzania" = "United Republic of Tanzania", "Turkiye" = "Türkiye", "United States" = "United States of America", "Venezuela, RB" = "Venezuela (Bolivarian Republic of)", "Virgin Islands (U.S.)" = "United States Virgin Islands", "West Bank and Gaza" = "Palestine", "Yemen, Rep." = "Yemen")

prop_fever_RI <- prop_fever_RI %>% mutate(Country = ifelse(Country %in% names(country_mapping), country_mapping[Country], Country))
prop_RI_fever <- prop_RI_fever %>% mutate(Country = ifelse(Country %in% names(country_mapping), country_mapping[Country], Country))

merged.data <- prop_fever_RI %>% inner_join(gbd.data, by = c("Country","Year"))
merged.data <- merged.data %>% inner_join(prop_RI_fever, by = c("Country","Year"))
dim(merged.data)

gbd.data <- merged.data[,colnames(gbd.data)]
prop_fever_RI <- merged.data[,colnames(prop_fever_RI)]
prop_RI_fever <- merged.data[,colnames(prop_RI_fever)]
prop_fever_RI_sim <- prop_fever_RI %>% dplyr::select(-c(Country, Year))
prop_fever_LRI_matrix <- prop_fever_RI_sim[,1:n_sim]
prop_fever_URI_matrix <- prop_fever_RI_sim[,(n_sim+1):(2*n_sim)]
prop_RI_fever_matrix <- prop_RI_fever %>% dplyr::select(-c(Country, Year))


# making distribution of LRI and URI incidence based on gamma distribution parameters
gbd.data.meta <- gbd.data %>% dplyr::select(Country, Year, U5_prop_LRI)
n_row <- nrow(gbd.data)
U5_URI_matrix <- U5_LRI_matrix <- as.data.frame(array(NA, dim = c(n_row, n_sim)))
last_progress <- -1 

start.time <- Sys.time()
n_cores <- detectCores() - 1

simulate_row <- function(i) {
  U5_LRI <- rgamma(n_sim, shape = gbd.data$U5_shape_LRI[i], scale = gbd.data$U5_scale_LRI[i])
  U5_URI <- rgamma(n_sim, shape = gbd.data$U5_shape_URI[i], scale = gbd.data$U5_scale_URI[i])
  list(LRI = U5_LRI, URI = U5_URI)
}

start.time <- Sys.time()
results <- pblapply(1:n_row, simulate_row, cl = n_cores)
end.time <- Sys.time()
end.time - start.time

U5_LRI_matrix <- do.call(rbind, lapply(results, function(x) x$LRI))
U5_URI_matrix <- do.call(rbind, lapply(results, function(x) x$URI))

# Compute the incidence of LRI with fever and URI with fever: RI * proportion of fever among children with RI
U5_fever_n_LRI_matrix <- U5_LRI_matrix*prop_fever_LRI_matrix
U5_fever_n_URI_matrix <- U5_URI_matrix*prop_fever_URI_matrix

# Compute the incidence of overall fever: (LRI with fever + URI with fever)/(proportion of RI among children with fever)
U5_fever_only_matrix <- ((U5_fever_n_LRI_matrix + U5_fever_n_URI_matrix) / prop_RI_fever_matrix ) - U5_fever_n_LRI_matrix
U5_LRI_only_matrix <- U5_LRI_matrix - U5_fever_n_LRI_matrix

# Incidence of LRI without Fever, Fever without LRI, and LRI with Fever
U5_matrix <- data.frame(
  Country = merged.data$Country,
  Year = merged.data$Year,
  U5_LRI_only = U5_LRI_only_matrix,
  U5_fever_only = U5_fever_only_matrix,
  U5_fever_LRI = U5_fever_n_LRI_matrix)

dim(U5_matrix)
saveRDS(U5_matrix, "outputs/sensitivity_analysis/SA_matrix.result.ver.1.rds")

U5_LRI_only_med <- apply(U5_LRI_only_matrix, 1, median)
U5_fever_only_med <- apply(U5_fever_only_matrix, 1, median)
U5_fever_LRI_med <- apply(U5_fever_n_LRI_matrix, 1, median)

U5_LRI_only_lower <- apply(U5_LRI_only_matrix, 1, quantile, probs = 0.025)
U5_fever_only_lower <- apply(U5_fever_only_matrix, 1, quantile, probs = 0.025)
U5_fever_LRI_lower <- apply(U5_fever_n_LRI_matrix, 1, quantile, probs = 0.025)

U5_LRI_only_upper <- apply(U5_LRI_only_matrix, 1, quantile, probs = 0.975)
U5_fever_only_upper <- apply(U5_fever_only_matrix, 1, quantile, probs = 0.975)
U5_fever_LRI_upper <- apply(U5_fever_n_LRI_matrix, 1, quantile, probs = 0.975)

U5_prop <- data.frame(
  Country = merged.data$Country,
  Year = merged.data$Year,
  U5_prop = merged.data$U5_prop
)

merged.data <- data.frame(
  Country = merged.data$Country,
  Year = merged.data$Year,
  U5_LRI_only = U5_LRI_only_med,
  U5_LRI_only_lower = U5_LRI_only_lower,
  U5_LRI_only_upper = U5_LRI_only_upper,
  U5_Fever_only = U5_fever_only_med,
  U5_Fever_only_lower = U5_fever_only_lower,
  U5_Fever_only_upper = U5_fever_only_upper,
  U5_Fever_n_LRI = U5_fever_LRI_med,
  U5_Fever_n_LRI_lower = U5_fever_LRI_lower,
  U5_Fever_n_LRI_upper = U5_fever_LRI_upper,
  U5_prop = merged.data$U5_prop
)

saveRDS(U5_prop, "outputs/sensitivity_analysis/SA_U5_prop.rds")
saveRDS(merged.data, "outputs/sensitivity_analysis/SA_Incidence_LRI_Fever_Both.ver.1.rds")



# Incidence of LRI (w/ and w/o fever) and Fever
U5_matrix <- data.frame(
  Country = merged.data$Country,
  Year = merged.data$Year,
  U5_LRI = U5_LRI_matrix,
  U5_fever_only = U5_fever_only_matrix)

saveRDS(U5_matrix, "outputs/sensitivity_analysis/SA_matrix.result.ver.2.rds")

dim(U5_matrix)

U5_LRI_med <- apply(U5_LRI_matrix, 1, median)
U5_fever_only_med <- apply(U5_fever_only_matrix, 1, median)

U5_LRI_lower <- apply(U5_LRI_matrix, 1, quantile, probs = 0.025)
U5_fever_only_lower <- apply(U5_fever_only_matrix, 1, quantile, probs = 0.025)

U5_LRI_upper <- apply(U5_LRI_matrix, 1, quantile, probs = 0.975)
U5_fever_only_upper <- apply(U5_fever_only_matrix, 1, quantile, probs = 0.975)

U5_prop <- data.frame(
  Country = merged.data$Country,
  Year = merged.data$Year,
  U5_prop = merged.data$U5_prop
)

merged.data <- data.frame(
  Country = merged.data$Country,
  Year = merged.data$Year,
  U5_LRI = U5_LRI_med,
  U5_LRI_lower = U5_LRI_lower,
  U5_LRI_upper = U5_LRI_upper,
  U5_Fever_only = U5_fever_only_med,
  U5_Fever_only_lower = U5_fever_only_lower,
  U5_Fever_only_upper = U5_fever_only_upper,
  U5_prop = merged.data$U5_prop
)

saveRDS(merged.data, "outputs/sensitivity_analysis/SA_Incidence_LRI_Fever_Both.ver.2.rds")









########## 4. Ab Consumption ##########



rm(list = ls())
gc()

requiredPackages = c("dplyr", "tidyr", 'stringr', "pbapply", "parallel")
lib.path <- .libPaths()
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p, lib = lib.path)
  library(p,character.only = TRUE)
}

options(scipen = 100, digits = 4) 

n_sim <- 10000
set.seed(123)

ab.data <- read.csv("data/DDD_country_year_class.csv")
colnames(ab.data)  
ab.data <- ab.data[,-3]
unique(ab.data$Year) 
unique(ab.data$Country) 
unique(ab.data$Sector) 
MainAb <- unique(ab.data$Antimicrobial)

merged.data <- readRDS("outputs/sensitivity_analysis/SA_Incidence_LRI_Fever_Both.ver.2.rds")

# matching country name
a <- unique(ab.data$Country)
b <- unique(merged.data$Country)
setdiff(a,b)
setdiff(b,a)
ab.data$Country <- trimws(ab.data$Country, which = "right")
country_mapping = c("Korea, Rep." = "Republic of Korea", "Czech Republic" = "Czechia", "Russia" = "Russian Federation", "United States" = "United States of America", "Venezuela" = "Venezuela (Bolivarian Republic of)", "Vietnam" = "Viet Nam")
ab.data <- ab.data %>% mutate(Country = ifelse(Country %in% names(country_mapping), country_mapping[Country], Country))



# Select included antibiotics
ab.data <- ab.data %>% group_by(Country, Year, Antimicrobial) %>% mutate(DDD_sum = sum(DDD, na.rm = TRUE))
ab.data <- ab.data[, c("Country", "Year", "Antimicrobial", "DDD_sum")]
ab.data <- ab.data %>% pivot_wider(names_from = Antimicrobial, values_from = DDD_sum)

colnames(ab.data) <- c("Country", "Year", "Aminoglycosides", "BSP", "Carbapenems",     
                       "Cephalosporins" , "Chloramphenicols" , "Fluoroquinolones" ,  "Glycopeptides"    , "Glycylcyclines"   ,          
                       "Lipopeptides"   , "Macrolides"       , "Monobactams"      ,  "NSP"              ,            "Others",          
                       "Oxazolidinones" , "Phosphonics"      , "Polymyxins"       ,  "Tetracyclines"    , "Trimethoprim"     ,          
                       "Rifampicin")

ab.data <- ab.data %>% dplyr::select(c(Country, Year, BSP, NSP, Cephalosporins, Macrolides, Aminoglycosides, Carbapenems, Glycopeptides, Trimethoprim))

# assume that the amount of antibiotics used for children aged 5 or younger is proportional to the incidence among children aged 5 or younger over overall incidence
U5_prop <- readRDS("outputs/sensitivity_analysis/SA_U5_prop.rds")
ab.data <- ab.data %>% inner_join(U5_prop, by = c("Country", "Year"))
ab.data[,3:10] <- ab.data[,3:10]*ab.data$U5_prop
ab.data

ab.data <- na.omit(ab.data)

for(i in 3:10){
  print(min(ab.data[,i]))
}

ab.data <- ab.data %>% filter(Carbapenems != 0) %>% filter(Glycopeptides != 0) %>% filter(NSP != 0)

U5_matrix <- readRDS("outputs/sensitivity_analysis/SA_matrix.result.ver.1.rds")
U5_matrix <- ab.data[,1:2] %>% inner_join(U5_matrix, by = c("Country", "Year"))
dim(U5_matrix)
U5_LRI_matrix <- U5_matrix[,3:(n_sim+2)] + U5_matrix[,(2*n_sim+3):(3*n_sim++2)]
U5_fever_only_matrix <- U5_matrix[,(n_sim+3):(2*n_sim+2)]
U5_all_matrix <- U5_matrix[,3:(n_sim+2)] + U5_matrix[,(n_sim+3):(2*n_sim+2)] + U5_matrix[,(2*n_sim+3):(3*n_sim++2)]

colnames(ab.data)

dim(U5_LRI_matrix)
dim(ab.data)

n_abx <- 8
U5_DDD_LRI <- U5_DDD_fever_only <- U5_DDD_all <- matrix(NA, nrow = nrow(ab.data), ncol = n_sim * n_abx)
U5_LRI_matrix <- apply(U5_LRI_matrix, 2, as.numeric)
U5_fever_only_matrix <- apply(U5_fever_only_matrix, 2, as.numeric)
U5_all_matrix <- apply(U5_all_matrix, 2, as.numeric)



# Loop to compute
for (i in 1:n_abx) {
  abx_col <- as.numeric(ab.data[[i + 2]]) 
  for (j in 1:n_sim) {
    col_index <- (i - 1) * n_sim + j
    U5_DDD_LRI[, col_index] <- abx_col / U5_LRI_matrix[, j]
    U5_DDD_fever_only[, col_index] <- abx_col / U5_fever_only_matrix[, j]
    U5_DDD_all[, col_index] <- abx_col / U5_all_matrix[, j]
  }
}

abx_names <- colnames(ab.data)[3:10]
col_names <- unlist(lapply(abx_names, function(name) paste0(name, "_sim", 1:n_sim)))
colnames(U5_DDD_LRI) <- col_names
U5_DDD_LRI_df <- cbind(ab.data[, 1:2], as.data.frame(U5_DDD_LRI))
U5_DDD_fever_only_df <- cbind(ab.data[, 1:2], as.data.frame(U5_DDD_fever_only))
U5_DDD_all_df <- cbind(ab.data[, 1:2], as.data.frame(U5_DDD_all))


chunk_size <- 10000  
total_cols <- ncol(U5_DDD_LRI) 
n_chunks <- total_cols / chunk_size  

LRI_median_list <- fever_only_median_list <- all_median_list <- 
  LRI_lower_list  <- fever_only_lower_list  <- all_lower_list  <- 
  LRI_upper_list  <- fever_only_upper_list  <- all_upper_list  <- list()

chunk_size <- n_sim; total_cols <- ncol(U5_DDD_LRI); n_chunks <- total_cols / chunk_size

# Loop through each chunk
for (i in 1:n_chunks) {
  cols <- ((i - 1) * chunk_size + 1):(i * chunk_size)
  
  LRI_median_list[[i]] <- apply(U5_DDD_LRI[, cols], 1, median)
  fever_only_median_list[[i]] <- apply(U5_DDD_fever_only[, cols], 1, median)
  all_median_list[[i]] <- apply(U5_DDD_all[, cols], 1, median)
  
  LRI_lower_list[[i]] <- apply(U5_DDD_LRI[, cols], 1, quantile, probs = 0.025)
  fever_only_lower_list[[i]] <- apply(U5_DDD_fever_only[, cols], 1, quantile, probs = 0.025)
  all_lower_list[[i]] <- apply(U5_DDD_all[, cols], 1, quantile, probs = 0.025)
  
  LRI_upper_list[[i]] <- apply(U5_DDD_LRI[, cols], 1, quantile, probs = 0.975)
  fever_only_upper_list[[i]] <- apply(U5_DDD_fever_only[, cols], 1, quantile, probs = 0.975)
  all_upper_list[[i]] <- apply(U5_DDD_all[, cols], 1, quantile, probs = 0.975)
}

#assigning antibiotics name
ab_names <- c("BSP", "NSP", "Cephalosporins", "Macrolides", 
              "Aminoglycosides", "Carbapenems", "Glycopeptides", "Trimethoprim")
names(LRI_median_list) <- names(LRI_lower_list) <- names(LRI_upper_list) <- 
  names(fever_only_median_list) <- names(fever_only_lower_list) <- names(fever_only_upper_list) <- 
  names(all_median_list) <- names(all_lower_list) <- names(all_upper_list) <- ab_names



ab.data_all <- ab.data_Fever_only <- ab.data_LRI <- ab.data_CountryYear <- data.frame(
  Country = ab.data$Country,
  Year = ab.data$Year)

for (ab in ab_names) {
  ab.data_LRI[[paste0("U5_DDD_LRI_", ab)]]         <- LRI_median_list[[ab]]
  ab.data_LRI[[paste0("U5_DDD_LRI_", ab, "_lower")]] <- LRI_lower_list[[ab]]
  ab.data_LRI[[paste0("U5_DDD_LRI_", ab, "_upper")]] <- LRI_upper_list[[ab]]
  
  ab.data_Fever_only[[paste0("U5_DDD_fever_", ab)]]         <- fever_only_median_list[[ab]]
  ab.data_Fever_only[[paste0("U5_DDD_fever_", ab, "_lower")]] <- fever_only_lower_list[[ab]]
  ab.data_Fever_only[[paste0("U5_DDD_fever_", ab, "_upper")]] <- fever_only_upper_list[[ab]]
  
  ab.data_all[[paste0("U5_DDD_fever_", ab)]]         <- all_median_list[[ab]]
  ab.data_all[[paste0("U5_DDD_fever_", ab, "_lower")]] <- all_lower_list[[ab]]
  ab.data_all[[paste0("U5_DDD_fever_", ab, "_upper")]] <- all_upper_list[[ab]]
  
  print(i)
}

GDP_PPP <- read.csv("data/GDP_PPP.csv")
GDP_PPP <- GDP_PPP[,1:22]
colnames(GDP_PPP) <- c("Country","2001",  "2002",  "2003",  "2004",  "2005",  "2006",  "2007",  "2008",  "2009",  "2010",  "2011",  "2012","2013", "2014",  "2015",  "2016",  "2017",  "2018",  "2019",  "2020",  "2021")
GDP_PPP <- GDP_PPP %>% pivot_longer(cols = -Country, names_to = "Year", values_to = "GDP")
GDP_PPP <- na.omit(GDP_PPP)

a <- unique(GDP_PPP$Country)
b <- unique(ab.data$Country)
setdiff(b,a)

country_mapping2 = c("Egypt" = "Egypt, Arab Rep.", "Republic of Korea" = "Korea, Rep.", "Slovakia" = "Slovak Republic", "Türkiye" ="Turkiye", "United States of America" = "United States")
GDP_PPP <- GDP_PPP %>% mutate(Country = ifelse(Country %in% names(country_mapping), country_mapping[Country], Country))
GDP_PPP$Year <- as.numeric(GDP_PPP$Year)
ab.data_LRI <- ab.data_LRI %>%left_join(GDP_PPP, by = c("Country","Year"))
ab.data_Fever_only <- ab.data_Fever_only %>%left_join(GDP_PPP, by = c("Country","Year"))
ab.data_all <- ab.data_all %>%left_join(GDP_PPP, by = c("Country","Year"))


U5_DDD_LRI       <- cbind(ab.data_CountryYear,U5_DDD_LRI)      
U5_DDD_Fever_only     <- cbind(ab.data_CountryYear,U5_DDD_fever_only) 
U5_DDD_all     <- cbind(ab.data_CountryYear,U5_DDD_all) 


saveRDS(U5_DDD_LRI      , "outputs/sensitivity_analysis/SA_U5_DDD_LRI_matrix.rds")
saveRDS(U5_DDD_Fever_only    , "outputs/sensitivity_analysis/SA_U5_DDD_Fever_only_matrix.rds")
saveRDS(U5_DDD_all    , "outputs/sensitivity_analysis/SA_U5_DDD_all_matrix.rds")

saveRDS(ab.data_LRI, "outputs/sensitivity_analysis/SA_ab.data_LRI.rds")
saveRDS(ab.data_Fever_only, "outputs/sensitivity_analysis/SA_ab.data_Fever_only.rds")
saveRDS(ab.data_all, "outputs/sensitivity_analysis/SA_ab.data_all.rds")








########## 5. adding Principal component ##########

rm(list = ls())
gc()

requiredPackages = c("dplyr", "readr", "data.table")
lib.path <- .libPaths()
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p, lib = lib.path)
  library(p,character.only = TRUE)
}

set.seed(123)
options(scipen = 999, digits = 4)
UNI <- readRDS("outputs/pc_gdp_ten.rds")

ab.data_LRI      <- readRDS("outputs/sensitivity_analysis/SA_ab.data_LRI.rds")
ab.data_Fever    <- readRDS("outputs/sensitivity_analysis/SA_ab.data_Fever_only.rds")
ab.data_all    <- readRDS("outputs/sensitivity_analysis/SA_ab.data_all.rds")
analysis_data <- readRDS("outputs/sensitivity_analysis/2-2.SA_analysis_data_prop.rds")

ab.data_LRI <- ab.data_LRI %>% dplyr::select(-GDP)
ab.data_Fever <- ab.data_Fever %>% dplyr::select(-GDP)
ab.data_all <- ab.data_all %>% dplyr::select(-GDP)

analysis_data <- analysis_data %>%
  dplyr::select(c("Antibiotics", "ARI", 
                  "ChildSex","ChildAgeMonth", 
                  "CountryYear", "Country", "Year", "AbLRI", "AbLRI_SE", "AbFever.nonARI", "AbFever.nonARI_SE", 
                  "n_AbLRI", "n_LRI", "n_Ab_Fever.nonARI", "n_Fever.nonARI"))

analysis_data <- na.omit(analysis_data) #188904
analysis_data <- analysis_data %>% mutate(PersonTime = 1/26) 
analysis_data$Year <- as.integer(analysis_data$Year)

head(UNI)
UNI_subset <- UNI
a <- unique(UNI_subset$Country)
b <- unique(ab.data_LRI$Country)
c <- unique(ab.data_Fever$Country)
d <- unique(analysis_data$Country)

analysis_data$Country[analysis_data$Country == "C\x99te d'Ivoire"] <- "Cote d'Ivoire"


country_mapping <- c("Korea, Rep."="Republic of Korea", "United States"="United States of America", "Egypt, Arab Rep."="Egypt", "Slovak Republic"="Slovakia", "Turkiye"="Türkiye")
UNI_subset <- UNI_subset %>% mutate(Country = ifelse(Country %in% names(country_mapping), country_mapping[Country], Country))
analysis_data$Country[analysis_data$Country == "Egypt, Arab Rep."] <- "Egypt"

analysis_data$Year <- as.integer(analysis_data$Year)

# ab.data mainly does not have data in african countries
ab.data_LRI          <- ab.data_LRI           %>% inner_join(UNI_subset, by = c("Country","Year"))
ab.data_Fever        <- ab.data_Fever         %>% inner_join(UNI_subset, by = c("Country","Year"))
ab.data_all          <- ab.data_all           %>% inner_join(UNI_subset, by = c("Country","Year"))
AbLRI.data           <- analysis_data         %>% inner_join(UNI_subset, by = c("Country","Year"))

lengths(ab.data_LRI) 
ab.data_LRI <- na.omit(ab.data_LRI) 
lengths(ab.data_Fever) 
ab.data_Fever <- na.omit(ab.data_Fever) 
lengths(ab.data_all) 
ab.data_all <- na.omit(ab.data_all) 




nrow(AbLRI.data) 
AbLRI.data <- na.omit(AbLRI.data)
nrow(AbLRI.data) 

years <- sort(unique(UNI$Year))
colors <- rainbow(length(years))  
year_colors <- colors[as.numeric(factor(UNI$Year, levels = years))]
plot(UNI$PC1,UNI$PC2, col = UNI$Year)

AbLRI.data_PC <- AbLRI.data %>% dplyr::select(PC1, PC2) %>% dplyr::filter(!duplicated(.))
points(AbLRI.data_PC$PC1, AbLRI.data_PC$PC2, col = 'red', pch = 19)
points(ab.data_LRI$PC1, ab.data_LRI$PC2, col = 'blue', pch = 19)

legend("topright", legend = years, col = colors, pch = 19, cex = 0.6, title = "Year")

saveRDS(AbLRI.data, "outputs/sensitivity_analysis/SA_AbLRI.data.UNI.rds")
saveRDS(ab.data_LRI, "outputs/sensitivity_analysis/SA_ab.data_LRI.UNI.rds")
saveRDS(ab.data_Fever, "outputs/sensitivity_analysis/SA_ab.data_Fever.UNI.rds")
saveRDS(ab.data_all, "outputs/sensitivity_analysis/SA_ab.data_all.UNI.rds")

dim(UNI)

write.csv(UNI, "outputs/sensitivity_analysis/SA_UNI.data.csv")










########## 6. regression ##########

#### 0. prep -------------
rm(list = ls())
gc()

requiredPackages = c("dplyr", "tidyr", "MASS")
lib.path <- .libPaths()
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p, lib = lib.path)
  library(p,character.only = TRUE)
}

n_sim <- 10000
set.seed(123)

#### 2. AbLRI ~ DDD per LRTI 
analysis_data <- readRDS("outputs/sensitivity_analysis/SA_AbLRI.data.UNI.rds")
AbAccess <- analysis_data %>%
  dplyr::select(c("ChildSex","ChildAgeMonth","Year","GDP","Income","Country",
                  "AbLRI","AbFever.nonARI","AbLRI_SE","AbFever.nonARI_SE","PC1","PC2","Endemic","gbd_region","n_AbLRI","n_LRI","n_Ab_Fever.nonARI","n_Fever.nonARI"))

AbAccess <- na.omit(AbAccess)
AbAccess <- AbAccess %>% mutate(PersonTime = 1/26) 
AbAccess$Year <- as.integer(AbAccess$Year)
unique(AbAccess$Country)

ab.data_LRI <- readRDS("outputs/sensitivity_analysis/SA_ab.data_LRI.UNI.rds")
ab.data_Fever <- readRDS("outputs/sensitivity_analysis/SA_ab.data_Fever.UNI.rds")
ab.data_all <- readRDS("outputs/sensitivity_analysis/SA_ab.data_all.UNI.rds")
matrix.LRI <- readRDS("outputs/sensitivity_analysis/SA_U5_DDD_LRI_matrix.rds")
matrix.Fever <- readRDS("outputs/sensitivity_analysis/SA_U5_DDD_Fever_only_matrix.rds")
matrix.all <- readRDS("outputs/sensitivity_analysis/SA_U5_DDD_all_matrix.rds")

matrix.meta <- ab.data_LRI %>% dplyr::select(Country,Year)
matrix.LRI <- matrix.meta %>% inner_join(matrix.LRI, by = c("Country","Year")) 
matrix.Fever <- matrix.meta %>% inner_join(matrix.Fever, by = c("Country","Year")) 
matrix.all <- matrix.meta %>% inner_join(matrix.all, by = c("Country","Year")) 

matrix.LRI.sum <- data.frame(matrix(NA_real_, nrow = nrow(matrix.LRI), ncol = (n_sim + 2)))
matrix.Fever.sum <- data.frame(matrix(NA_real_, nrow = nrow(matrix.Fever), ncol = (n_sim + 2)))
matrix.all.sum <- data.frame(matrix(NA_real_, nrow = nrow(matrix.all), ncol = (n_sim + 2)))
matrix.all.sum[, 1:2] <- matrix.Fever.sum[, 1:2] <- matrix.LRI.sum[, 1:2] <- as.matrix(matrix.LRI[, 1:2])

# sum antibiotic classes
for (i in 1:n_sim) {
  cols_to_sum <- seq(from = i + 2, to = (8*n_sim+2), by = n_sim)
  matrix.LRI.sum[, i + 2] <- as.numeric(rowSums(matrix.LRI[, cols_to_sum]))
  matrix.Fever.sum[, i + 2] <- as.numeric(rowSums(matrix.Fever[, cols_to_sum]))
  matrix.all.sum[, i + 2] <- as.numeric(rowSums(matrix.all[, cols_to_sum]))
}

dim(matrix.LRI.sum); dim(matrix.Fever.sum); dim(matrix.all.sum)
str(matrix.LRI.sum); str(matrix.Fever.sum); str(matrix.all.sum)
colnames(matrix.LRI.sum)[1:2] <- c("Country", "Year")
colnames(matrix.LRI.sum)[3:ncol(matrix.LRI.sum)] <- paste0("sim_", 1:n_sim)
colnames(matrix.Fever.sum)[1:2] <- c("Country", "Year")
colnames(matrix.Fever.sum)[3:ncol(matrix.Fever.sum)] <- paste0("sim_", 1:n_sim)
colnames(matrix.all.sum)[1:2] <- c("Country", "Year")
colnames(matrix.all.sum)[3:ncol(matrix.all.sum)] <- paste0("sim_", 1:n_sim)

ab.data_LRI$DDD_sum <- ab.data_LRI %>% dplyr::select(c("U5_DDD_LRI_BSP","U5_DDD_LRI_NSP","U5_DDD_LRI_Cephalosporins","U5_DDD_LRI_Macrolides",
                                                       "U5_DDD_LRI_Aminoglycosides","U5_DDD_LRI_Carbapenems","U5_DDD_LRI_Glycopeptides","U5_DDD_LRI_Trimethoprim")) %>% rowSums()

ab.data_Fever$DDD_sum <- ab.data_Fever %>% dplyr::select(c("U5_DDD_fever_BSP","U5_DDD_fever_NSP","U5_DDD_fever_Cephalosporins","U5_DDD_fever_Macrolides",
                                                           "U5_DDD_fever_Aminoglycosides","U5_DDD_fever_Carbapenems","U5_DDD_fever_Glycopeptides","U5_DDD_fever_Trimethoprim")) %>% rowSums()

ab.data_all$DDD_sum <- ab.data_all %>% dplyr::select(c("U5_DDD_fever_BSP","U5_DDD_fever_NSP","U5_DDD_fever_Cephalosporins","U5_DDD_fever_Macrolides",
                                                       "U5_DDD_fever_Aminoglycosides","U5_DDD_fever_Carbapenems","U5_DDD_fever_Glycopeptides","U5_DDD_fever_Trimethoprim")) %>% rowSums()





### 1-1. Ab_LRI ~ PC: logistic regression ---------------------------

AbAccess_Country_level <- AbAccess %>%
  dplyr::select(
    Country, Year, PC1, PC2,
    AbLRI, AbLRI_SE,
    n_AbLRI, n_LRI,
    n_Ab_Fever.nonARI, n_Fever.nonARI
  ) %>%
  distinct()

write.csv(
  AbAccess_Country_level,
  "outputs/draft_df/SA_Suppl_1_DHS_model_data.csv",
  row.names = FALSE
)

mod_AbLRI_logit <- glm(
  cbind(n_AbLRI, n_LRI - n_AbLRI) ~ PC1 + PC2,
  family = binomial(link = "logit"),
  data = AbAccess_Country_level
)

coef_AbLRI <- coef(mod_AbLRI_logit)
vcov_AbLRI <- vcov(mod_AbLRI_logit)

coef_matrix_LRI <- MASS::mvrnorm(
  n = n_sim,
  mu = coef_AbLRI,
  Sigma = vcov_AbLRI
)

coef_df_AbLRI <- as.data.frame(coef_matrix_LRI)



### 1-2. Ab_Fever ~ PC: logistic regression -----------------------

AbAccess_Country_level <- AbAccess %>%
  dplyr::select(
    Country, Year, PC1, PC2,
    AbFever.nonARI, AbFever.nonARI_SE,
    n_Ab_Fever.nonARI, n_Fever.nonARI
  ) %>%
  distinct()

mod_AbFever_logit <- glm(
  cbind(n_Ab_Fever.nonARI, n_Fever.nonARI - n_Ab_Fever.nonARI) ~ PC1 + PC2,
  family = binomial(link = "logit"),
  data = AbAccess_Country_level
)

coef_AbFever <- coef(mod_AbFever_logit)
vcov_AbFever <- vcov(mod_AbFever_logit)

coef_matrix_Fever <- MASS::mvrnorm(
  n = n_sim,
  mu = coef_AbFever,
  Sigma = vcov_AbFever
)

coef_df_AbFever <- as.data.frame(coef_matrix_Fever)



# Poisson for comparison ------------------------------
AbAccess_Country_level_pois <- AbAccess %>% dplyr::select(Country, Year, PC1, PC2, AbLRI, AbLRI_SE, n_AbLRI, n_LRI, n_Ab_Fever.nonARI, n_Fever.nonARI) %>% 
  distinct()

### 1-1. Ab_LRI ~ PC ###
mod_AbLRI_pois <- glm(n_AbLRI ~ PC1 + PC2 + offset(log(n_LRI)),
                      family = poisson, data = AbAccess_Country_level_pois)

coef_AbLRI_pois <- coef(mod_AbLRI_pois); vcov_AbLRI_pois <- vcov(mod_AbLRI_pois)
coef_matrix_LRI_pois <- MASS::mvrnorm(n = n_sim, mu = coef_AbLRI_pois, Sigma = vcov_AbLRI_pois)  
coef_df_AbLRI_pois <- as.data.frame(coef_matrix_LRI_pois)


### 1-2. Ab_Fever ~ PC ###
mod_AbFever_pois <- glm(n_Ab_Fever.nonARI ~ PC1 + PC2 + offset(log(n_Fever.nonARI)),
                        family = poisson, data = AbAccess_Country_level_pois)

coef_AbFever_pois <- coef(mod_AbFever_pois); vcov_AbFever_pois <- vcov(mod_AbFever_pois)
coef_matrix_Fever_pois <- MASS::mvrnorm(n = n_sim, mu = coef_AbFever_pois, Sigma = vcov_AbFever_pois)  
coef_df_AbFever_pois <- as.data.frame(coef_matrix_Fever_pois)





### comparing plots --------------------------------------
pdf("outputs/sensitivity_analysis/SA_(B) Suppl_Figure_logit_pois_Ab_LRTI_and_Fever.pdf", width = 6, height = 2.5)

par(mfrow = c(1, 2))
par(mar = c(3, 3.5, 1.5, 2.5))  
par(mgp = c(2.2, 0.43, 0))
par(tck = -0.025)
par(lwd = 1)

### LRTI

plot(
  mod_AbLRI_logit$fitted.values,
  mod_AbLRI_pois$fitted.values / AbAccess_Country_level_pois$n_LRI,
  xlim = c(0.5, 0.9),
  ylim = c(0.5, 0.9),
  pch = 16,
  cex = 0.65,
  axes = FALSE,
  xlab = "",
  ylab = ""
)

axis(
  1,
  at = seq(0.5, 0.9, by = 0.1),
  cex.axis = 0.8,
  mgp = c(3, 0.20, 0)
)

axis(
  2,
  at = seq(0.5, 0.9, by = 0.1),
  las = 1,
  cex.axis = 0.8
)

box(bty = "l")
abline(0, 1, lty = 2)

title(
  main = "P(Ab | LRTI)",
  cex.main = 0.9,
  font.main = 1
)

usr <- par("usr")

text(
  x = mean(usr[1:2]),
  y = usr[3] - 0.3 * diff(usr[3:4]),
  labels = "Predicted proportion: logistic regression",
  cex = 0.8,
  xpd = NA
)

text(
  x = usr[1] - 0.3 * diff(usr[1:2]),
  y = mean(usr[3:4]),
  labels = "Predicted proportion: Poisson",
  srt = 90,
  cex = 0.8,
  xpd = NA
)


### Fever

plot(
  mod_AbFever_logit$fitted.values,
  mod_AbFever_pois$fitted.values / AbAccess_Country_level_pois$n_Fever.nonARI,
  xlim = c(0.4, 0.7),
  ylim = c(0.4, 0.7),
  pch = 16,
  cex = 0.65,
  axes = FALSE,
  xlab = "",
  ylab = ""
)

axis(
  1,
  at = seq(0.4, 0.7, by = 0.1),
  cex.axis = 0.8,
  mgp = c(3, 0.20, 0)
)

axis(
  2,
  at = seq(0.4, 0.7, by = 0.1),
  las = 1,
  cex.axis = 0.8
)

box(bty = "l")
abline(0, 1, lty = 2)

title(
  main = "P(Ab | Fever)",
  cex.main = 0.9,
  font.main = 1
)

usr <- par("usr")

text(
  x = mean(usr[1:2]),
  y = usr[3] - 0.3 * diff(usr[3:4]),
  labels = "Predicted proportion: logistic regression",
  cex = 0.8,
  xpd = NA
)

dev.off()





# DDD analysis ---------------------------------------

meta_LRI <- ab.data_LRI %>% dplyr::select(c(Country, Year))
meta_Fever <- ab.data_Fever %>% dplyr::select(c(Country, Year))

meta_LRI$Year <- as.integer(meta_LRI$Year)
meta_Fever$Year <- as.integer(meta_Fever$Year)
matrix.LRI.sum$Year <- as.integer(matrix.LRI.sum$Year)
matrix.Fever.sum$Year <- as.integer(matrix.Fever.sum$Year)
matrix.LRI.sum <- meta_LRI %>% inner_join(matrix.LRI.sum, by = c("Country", "Year"))
matrix.Fever.sum <- meta_Fever %>% inner_join(matrix.Fever.sum, by = c("Country", "Year"))

tab_df <- as.data.frame(with(ab.data_LRI, table(Country, Year)))
tab_wide <- tab_df %>% pivot_wider(names_from = Year, values_from = Freq, values_fill = 0)
write.csv(tab_wide, "outputs/draft_df/SA_Suppl_2_IQVIA_IHME_data.csv", row.names = FALSE)


### 1-3. DDD(LRI) ~ PC -------------------------

LRI.DDD_sim <- matrix.LRI.sum %>% dplyr::select(-Country, -Year) %>% as.data.frame
mod <- glm(log(DDD_sum) ~ PC1 + PC2, data = ab.data_LRI)
coef_matrix_LRI_DDD <- coef_matrix_Fever_DDD <- coef_matrix_all_DDD <- matrix(NA, nrow = n_sim, ncol = length(coef(mod)))
colnames(coef_matrix_LRI_DDD) <- names(coef(mod))

for (i in 1:n_sim) {
  ab.data_LRI$LRI.DDD_sim <- LRI.DDD_sim[, i]
  mod <- glm(log(LRI.DDD_sim) ~ PC1 + PC2, data = ab.data_LRI)
  
  coef_mean <- coef(mod); coef_vcov <- vcov(mod)
  
  coef_matrix_LRI_DDD[i, ] <- MASS::mvrnorm(1, mu = coef_mean, Sigma = coef_vcov)
  
  if (i %% 100 == 0) cat("Completed simulation:", i, "\n")
}

coef_df_LRI_DDD <- as.data.frame(coef_matrix_LRI_DDD)



### 1-4. DDD(Fever) ~ PC --------------------------------

Fever.DDD_sim <- matrix.Fever.sum %>% dplyr::select(-Country, -Year) %>% as.data.frame
mod <- glm(log(DDD_sum) ~ PC1 + PC2, data = ab.data_Fever)
colnames(coef_matrix_Fever_DDD) <- names(coef(mod))

for (i in 1:n_sim) {
  ab.data_Fever$Fever.DDD_sim <- Fever.DDD_sim[, i]
  mod <- glm(log(Fever.DDD_sim) ~ PC1 + PC2, data = ab.data_Fever)
  
  coef_mean <- coef(mod); coef_vcov <- vcov(mod)
  
  coef_matrix_Fever_DDD[i, ] <- MASS::mvrnorm(1, mu = coef_mean, Sigma = coef_vcov)
  
  if (i %% 100 == 0) cat("Completed simulation:", i, "\n")
}

coef_df_Fever_DDD <- as.data.frame(coef_matrix_Fever_DDD)



### 1-5. DDD(all) ~ PC

all.DDD_sim <- matrix.all.sum %>% dplyr::select(-Country, -Year) %>% as.data.frame
mod <- glm(log(DDD_sum) ~ PC1 + PC2, data = ab.data_all)
colnames(coef_matrix_all_DDD) <- names(coef(mod))

for (i in 1:n_sim) {
  ab.data_all$all.DDD_sim <- all.DDD_sim[, i]
  mod <- glm(log(all.DDD_sim) ~ PC1 + PC2, data = ab.data_all)
  
  coef_mean <- coef(mod); coef_vcov <- vcov(mod)
  
  coef_matrix_all_DDD[i, ] <- MASS::mvrnorm(1, mu = coef_mean, Sigma = coef_vcov)
  
  if (i %% 100 == 0) cat("Completed simulation:", i, "\n")
}

coef_df_all_DDD <- as.data.frame(coef_matrix_all_DDD)



### 1.6 AbLRI/Fever ~ DDD

Abaccess_PC <- AbAccess %>% dplyr::select(Country, Year, PC1, PC2,GDP, Income, Endemic, gbd_region, AbLRI, AbFever.nonARI)
Abaccess_PC <- distinct(Abaccess_PC)
colnames(Abaccess_PC) <- c("Country", "Year","PC1","PC2","GDP","Income","Endemic","gbd_region","AbLRI","AbFever.nonARI")

ab.data_LRI_PC <- ab.data_LRI %>% dplyr::select(Country, Year, PC1, PC2, GDP, Income, Endemic, gbd_region, DDD_sum)
ab.data_Fever_PC <- ab.data_Fever %>% dplyr::select(Country, Year, PC1, PC2, GDP, Income, Endemic, gbd_region, DDD_sum)
ab.data_all_PC <- ab.data_all %>% dplyr::select(Country, Year, PC1, PC2, GDP, Income, Endemic, gbd_region, DDD_sum)

colnames(ab.data_LRI_PC) <- c("Country", "Year","PC1","PC2","GDP","Income","Endemic","gbd_region","DDD_sum_LRI")
colnames(ab.data_Fever_PC) <- c("Country", "Year","PC1","PC2","GDP","Income","Endemic","gbd_region","DDD_sum_Fever")
colnames(ab.data_all_PC) <- c("Country", "Year","PC1","PC2","GDP","Income","Endemic","gbd_region","DDD_sum_all")

UNI <- readRDS("outputs/pc_gdp_ten.rds")
colnames(Abaccess_PC); colnames(ab.data_LRI_PC); colnames(UNI)

UNI <- UNI %>% dplyr::select(c("Country","Year","PC1","PC2","GDP","Income","Endemic","gbd_region","Population"))

Abaccess_PC$Abaccess <- 1; ab.data_LRI_PC$ab.data <- 1

UNI_indexed <- UNI %>% left_join(Abaccess_PC %>% dplyr::select(Country, Year, Abaccess, AbLRI, AbFever.nonARI), by = c("Country", "Year"))
UNI_indexed <- UNI_indexed %>% mutate(Abaccess = ifelse(is.na(Abaccess), 0, Abaccess), AbLRI = ifelse(is.na(AbLRI), 0, AbLRI), AbFever.nonARI = ifelse(is.na(AbFever.nonARI), 0, AbFever.nonARI))

UNI_indexed <- UNI_indexed %>% left_join(ab.data_LRI_PC %>% dplyr::select(Country, Year, ab.data, DDD_sum_LRI), by = c("Country", "Year"))
UNI_indexed <- UNI_indexed %>% left_join(ab.data_Fever_PC %>% dplyr::select(Country, Year, DDD_sum_Fever), by = c("Country", "Year"))

UNI_indexed <- UNI_indexed %>% mutate(ab.data = ifelse(is.na(ab.data), 0, ab.data), DDD_sum_LRI = ifelse(is.na(DDD_sum_LRI), 0, DDD_sum_LRI), DDD_sum_Fever = ifelse(is.na(DDD_sum_Fever), 0, DDD_sum_Fever))

with(UNI_indexed, table(Abaccess, ab.data, useNA = 'always'))
UNI_indexed %>% filter(Abaccess == 1, ab.data == 1)

head(UNI_indexed); str(UNI_indexed)
head(coef_df_AbLRI); str(coef_df_AbLRI)
head(coef_df_AbFever); str(coef_df_AbFever)
head(coef_df_LRI_DDD); str(coef_df_LRI_DDD)
head(coef_df_Fever_DDD); str(coef_df_Fever_DDD)
head(coef_df_all_DDD); str(coef_df_all_DDD)

# confirm all final coefficient objects use intercept + PC1 + PC2
stopifnot(
  identical(names(coef_df_AbLRI), c("(Intercept)", "PC1", "PC2")),
  identical(names(coef_df_AbFever), c("(Intercept)", "PC1", "PC2")),
  identical(names(coef_df_LRI_DDD), c("(Intercept)", "PC1", "PC2")),
  identical(names(coef_df_Fever_DDD), c("(Intercept)", "PC1", "PC2")),
  identical(names(coef_df_all_DDD), c("(Intercept)", "PC1", "PC2"))
)

# List of PC variables

pc_vars <- paste0("PC", 1:2)
par(mfrow = c(2, 1))

for (pc in pc_vars) {
  pc_values <- UNI_indexed[[pc]]
  xlim_all <- range(pc_values, na.rm = TRUE)
  breaks_seq <- seq(floor(xlim_all[1]), ceiling(xlim_all[2]), by = 0.5)
  
  # Extract groups
  
  group_A <- pc_values[UNI_indexed$Abaccess == 1]
  group_B <- pc_values[UNI_indexed$ab.data == 1]
  group_C <- pc_values
  
  # Plot histograms with same breaks and xlim
  
  hist(group_A, breaks = breaks_seq, col = "gray", main = paste(pc, ": Abaccess = 1"), xlab = pc, xlim = xlim_all)
  hist(group_B, breaks = breaks_seq, col = "gray", main = paste(pc, ": ab.data_LRI = 1"), xlab = pc, xlim = xlim_all)
  hist(group_C, breaks = breaks_seq, col = "gray", main = paste(pc, ": All"), xlab = pc, xlim = xlim_all)
}

par(mfrow = c(1, 1))

saveRDS(UNI_indexed, "outputs/sensitivity_analysis/SA_UNI_indexed.rds")
saveRDS(coef_df_AbLRI, "outputs/sensitivity_analysis/SA_coef_df_LRI.rds")
saveRDS(coef_df_AbFever, "outputs/sensitivity_analysis/SA_coef_df_Fever.rds")
saveRDS(coef_df_LRI_DDD, "outputs/sensitivity_analysis/SA_coef_df_LRI_DDD.rds")
saveRDS(coef_df_Fever_DDD, "outputs/sensitivity_analysis/SA_coef_df_Fever_DDD.rds")
saveRDS(coef_df_all_DDD, "outputs/sensitivity_analysis/SA_coef_df_all_DDD.rds")

par(mfrow = c(1,1))



########## 7. Absolute number ##########

#### 0. prep

rm(list = ls())
gc()

requiredPackages = c("dplyr", "tidyr", "data.table", "stringr", "MASS", "naniar")
lib.path <- .libPaths()
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p, lib = lib.path)
  library(p,character.only = TRUE)
}

n_sim <- 10000
set.seed(123)

#---------------------------prepare data for figure------------------------
UNI_indexed <- readRDS("outputs/sensitivity_analysis/SA_UNI_indexed.rds")
coef_df_LRI <- readRDS("outputs/sensitivity_analysis/SA_coef_df_LRI.rds")
coef_df_Fever <- readRDS("outputs/sensitivity_analysis/SA_coef_df_Fever.rds")
coef_df_LRI_DDD <- readRDS("outputs/sensitivity_analysis/SA_coef_df_LRI_DDD.rds")
coef_df_Fever_DDD <- readRDS("outputs/sensitivity_analysis/SA_coef_df_Fever_DDD.rds")
coef_df_all_DDD <- readRDS("outputs/sensitivity_analysis/SA_coef_df_all_DDD.rds")

# confirm all models use intercept + PC1 + PC2
stopifnot(
  identical(names(coef_df_LRI), c("(Intercept)", "PC1", "PC2")),
  identical(names(coef_df_Fever), c("(Intercept)", "PC1", "PC2")),
  identical(names(coef_df_LRI_DDD), c("(Intercept)", "PC1", "PC2")),
  identical(names(coef_df_Fever_DDD), c("(Intercept)", "PC1", "PC2")),
  identical(names(coef_df_all_DDD), c("(Intercept)", "PC1", "PC2"))
)

UNI_indexed_Abaccess <- UNI_indexed %>% filter(Abaccess == 1)
UNI_indexed_Abdata <- UNI_indexed %>% filter(ab.data == 1)
UNI_indexed_both <- UNI_indexed %>% filter(Abaccess == 1, ab.data == 1)

with(UNI_indexed, table(Country, Year))

# Antibiotic treatment: PC1-PC2
X_pred <- model.matrix(~ PC1 + PC2, data = UNI_indexed)
X_pred_LRI <- X_pred[, names(coef_df_LRI), drop = FALSE]
X_pred_Fever <- X_pred[, names(coef_df_Fever), drop = FALSE]

# DDD: PC1-PC2
X_pred_DDD <- model.matrix(~ PC1 + PC2, data = UNI_indexed)
X_pred_LRI_DDD <- X_pred_DDD[, names(coef_df_LRI_DDD), drop = FALSE]
X_pred_Fever_DDD <- X_pred_DDD[, names(coef_df_Fever_DDD), drop = FALSE]
X_pred_all_DDD <- X_pred_DDD[, names(coef_df_all_DDD), drop = FALSE]

LRI_matrix <- as.matrix(X_pred_LRI) %*% t(as.matrix(coef_df_LRI))
Fever_matrix <- as.matrix(X_pred_Fever) %*% t(as.matrix(coef_df_Fever))
LRI_DDD_matrix <- as.matrix(X_pred_LRI_DDD) %*% t(as.matrix(coef_df_LRI_DDD))
Fever_DDD_matrix <- as.matrix(X_pred_Fever_DDD) %*% t(as.matrix(coef_df_Fever_DDD))
all_DDD_matrix <- as.matrix(X_pred_all_DDD) %*% t(as.matrix(coef_df_all_DDD))

# Convert logit-scale predictions to antibiotic treatment probabilities
AbLRI_sim <- plogis(LRI_matrix)  
AbFever_sim <- plogis(Fever_matrix)  

# Convert log-scale predictions to expected DDD
LRI_DDD_sim <- exp(LRI_DDD_matrix)  
Fever_DDD_sim <- exp(Fever_DDD_matrix)  
all_DDD_sim <- exp(all_DDD_matrix)  

epsilon_min <- 1e-6
epsilon_max <- 1e-6 + 1e-9

LRI_DDD_sim[LRI_DDD_sim <= 0] <- runif(sum(LRI_DDD_sim <= 0), min = epsilon_min, max = epsilon_max)
Fever_DDD_sim[Fever_DDD_sim <= 0] <- runif(sum(Fever_DDD_sim <= 0), min = epsilon_min, max = epsilon_max)
all_DDD_sim[all_DDD_sim <= 0] <- runif(sum(all_DDD_sim <= 0), min = epsilon_min, max = epsilon_max)

# Flatten the matrices into vectors (column-wise)

AbLRI_long <- as.vector(AbLRI_sim)
LRI_DDD_long <- as.vector(LRI_DDD_sim)
AbFever_long <- as.vector(AbFever_sim)
Fever_DDD_long <- as.vector(Fever_DDD_sim)
all_DDD_long <- as.vector(all_DDD_sim)

# Create matching indices

Year_index <- rep(UNI_indexed$Year, times = ncol(AbLRI_sim))
Country_index <- rep(UNI_indexed$Country, times = ncol(AbLRI_sim))
GDP_index <- rep(UNI_indexed$GDP, times = ncol(AbLRI_sim))
Income_index <- rep(UNI_indexed$Income, times = ncol(AbLRI_sim))
Malaria_index <- rep(UNI_indexed$Endemic, times = ncol(AbLRI_sim))
GBD_index <- rep(UNI_indexed$gbd_region, times = ncol(AbLRI_sim))
Population <- rep(UNI_indexed$Population, times = ncol(AbLRI_sim))

# Combine into dataframe

sim_df_LRI <- data.frame(
  Country = Country_index,
  Year = Year_index,
  GDP = GDP_index,
  Income = Income_index,
  Malaria = Malaria_index,
  GBD_region = GBD_index,
  AbLRI = AbLRI_long,
  Population = Population,
  sim = rep(1:n_sim, each = nrow(AbLRI_sim))
)

# Combine into dataframe

sim_df_Fever <- data.frame(
  Country = Country_index,
  Year = Year_index,
  GDP = GDP_index,
  Income = Income_index,
  Malaria = Malaria_index,
  GBD_region = GBD_index,
  AbFever = AbFever_long,
  Population = Population,
  sim = rep(1:n_sim, each = nrow(AbFever_sim))
)

AbLRI_matrix <- sim_df_LRI %>% pivot_wider(names_from = sim, values_from = AbLRI)
AbFever_matrix <- sim_df_Fever %>% pivot_wider(names_from = sim, values_from = AbFever)

Country_Year_Antibiotics <- data.frame(
  Country = AbLRI_matrix$Country,
  Year = AbLRI_matrix$Year)

U5_matrix <- readRDS("outputs/sensitivity_analysis/SA_matrix.result.ver.1.rds")

Country_Year_Incidence <- data.frame(
  Country = U5_matrix$Country,
  Year = U5_matrix$Year)

a <- unique(Country_Year_Antibiotics$Country)
b <- unique(Country_Year_Incidence$Country)
setdiff(a,b)
setdiff(b,a)

country_mapping <- c(
  "Bolivia" = "Bolivia (Plurinational State of)",
  "Congo, Dem. Rep." = "Democratic Republic of the Congo",
  "Congo, Rep." = "Congo",
  "Cote d'Ivoire" = "Côte d'Ivoire",
  "Egypt, Arab Rep." = "Egypt",
  "Gambia, The" = "Gambia",
  "Iran, Islamic Rep." = "Iran (Islamic Republic of)",
  "Kyrgyz Republic" = "Kyrgyzstan",
  "Lao PDR" = "Lao People's Democratic Republic",
  "Moldova" = "Republic of Moldova",
  "St. Vincent and the Grenadines" = "Saint Vincent and the Grenadines",
  "Tanzania" = "United Republic of Tanzania",
  "Turkiye" = "Türkiye",
  "Venezuela, RB" = "Venezuela (Bolivarian Republic of)",
  "West Bank and Gaza" = "Palestine",
  "Yemen, Rep." = "Yemen" 
)

Country_Year_Antibiotics <- Country_Year_Antibiotics %>% mutate(Country = ifelse(Country %in% names(country_mapping), country_mapping[Country], Country))
AbLRI_matrix <- AbLRI_matrix %>% mutate(Country = ifelse(Country %in% names(country_mapping), country_mapping[Country], Country))
AbFever_matrix <- AbFever_matrix %>% mutate(Country = ifelse(Country %in% names(country_mapping), country_mapping[Country], Country))

Country_Year_meta <- inner_join(Country_Year_Antibiotics, Country_Year_Incidence, by = c("Country", "Year"))

Incidence <- Country_Year_meta %>% inner_join(U5_matrix, by = c("Country","Year"))
Incidence_LRI_matrix <- Incidence[,3:(n_sim+2)] + Incidence[,(2*n_sim+3):(3*n_sim+2)]
Incidence_Fever_matrix <- Incidence[,(n_sim+3):(2*n_sim+2)]
Incidence_all_matrix <- Incidence[,3:(n_sim+2)] + Incidence[,(n_sim+3):(2*n_sim+2)] + Incidence[,(2*n_sim+3):(3*n_sim+2)]

Antibiotics_LRI <- Country_Year_meta %>% inner_join(AbLRI_matrix, by = c("Country","Year"))
Antibiotics_LRI_matrix <- Antibiotics_LRI[,8:(n_sim+7)]
Antibiotics_Fever <- Country_Year_meta %>% inner_join(AbFever_matrix, by = c("Country","Year"))
Antibiotics_Fever_matrix <- Antibiotics_Fever[,8:(n_sim+7)]

dim(Incidence_LRI_matrix)
dim(Antibiotics_LRI_matrix)
treated_LRI_matrix <- Incidence_LRI_matrix * Antibiotics_LRI_matrix
untreated_LRI_matrix <- Incidence_LRI_matrix - treated_LRI_matrix

dim(Incidence_Fever_matrix)
dim(Antibiotics_Fever_matrix)
treated_Fever_matrix <- Incidence_Fever_matrix * Antibiotics_Fever_matrix
untreated_Fever_matrix <- Incidence_Fever_matrix - treated_Fever_matrix


### for gini index

saveRDS(Country_Year_meta,             "outputs/global_map_df/gini/SA_Country_Year_meta.rds")
saveRDS(Antibiotics_LRI[,c(1,2,3,5)], "outputs/global_map_df/gini/SA_Antibiotics_LRI.rds")
saveRDS(treated_LRI_matrix,            "outputs/global_map_df/gini/SA_treated_LRI_matrix.rds")
saveRDS(Incidence_LRI_matrix,          "outputs/global_map_df/gini/SA_Incidence_LRI_matrix.rds")
saveRDS(Antibiotics_Fever[,c(1,2,3,5)], "outputs/global_map_df/gini/SA_Antibiotics_Fever.rds")
saveRDS(treated_Fever_matrix,          "outputs/global_map_df/gini/SA_treated_Fever_matrix.rds")
saveRDS(Incidence_Fever_matrix,        "outputs/global_map_df/gini/SA_Incidence_Fever_matrix.rds")


#-----------------------1.1 by country and year (treatment)------------------------

make_summary <- function(meta, income, treated, untreated, total) {
  
  # Compute proportion treated
  
  prop_treated <- treated / (treated + untreated)
  inc_rate_total <- (treated + untreated) / income$Population
  
  data.frame(
    Country = meta$Country,
    Year = meta$Year,
    Income = income$Income,
    Malaria = income$Malaria,
    GBD_region = income$GBD_region,
    Population = income$Population,
    
    Treated = sprintf("%s (%s  %s)",
                      formatC(apply(treated, 1, median, na.rm = TRUE), format = "f", big.mark = ",", digits = 0),
                      formatC(apply(treated, 1, quantile, probs = 0.025, na.rm = TRUE), format = "f", big.mark = ",", digits = 0),
                      formatC(apply(treated, 1, quantile, probs = 0.975, na.rm = TRUE), format = "f", big.mark = ",", digits = 0)),
    
    Untreated = sprintf("%s (%s  %s)",
                        formatC(apply(untreated, 1, median, na.rm = TRUE), format = "f", big.mark = ",", digits = 0),
                        formatC(apply(untreated, 1, quantile, probs = 0.025, na.rm = TRUE), format = "f", big.mark = ",", digits = 0),
                        formatC(apply(untreated, 1, quantile, probs = 0.975, na.rm = TRUE), format = "f", big.mark = ",", digits = 0)),
    
    Total = sprintf("%s (%s  %s)",
                    formatC(apply(total, 1, median, na.rm = TRUE), format = "f", big.mark = ",", digits = 0),
                    formatC(apply(total, 1, quantile, probs = 0.025, na.rm = TRUE), format = "f", big.mark = ",", digits = 0),
                    formatC(apply(total, 1, quantile, probs = 0.975, na.rm = TRUE), format = "f", big.mark = ",", digits = 0)),
    
    Incidence_rate = sprintf("%s (%s  %s)",
                             formatC(apply(inc_rate_total, 1, median, na.rm = TRUE), format = "f", big.mark = ",", digits = 3),
                             formatC(apply(inc_rate_total, 1, quantile, probs = 0.025, na.rm = TRUE), format = "f", big.mark = ",", digits = 3),
                             formatC(apply(inc_rate_total, 1, quantile, probs = 0.975, na.rm = TRUE), format = "f", big.mark = ",", digits = 3)),
    
    Treated_Proportion = sprintf("%.3f (%.3f, %.3f)",
                                 apply(prop_treated, 1, median, na.rm = TRUE),
                                 apply(prop_treated, 1, quantile, probs = 0.025, na.rm = TRUE),
                                 apply(prop_treated, 1, quantile, probs = 0.975, na.rm = TRUE))
  )
}

make_df <- function(meta, income, treated, untreated, total) {
  
  # Compute proportion treated
  
  prop_treated <- treated / (treated + untreated)
  inc_rate_total <- (treated + untreated) / income$Population
  
  data.frame(
    Country = meta$Country,
    Year = meta$Year,
    GDP = income$GDP,
    Malaria = income$Malaria,
    Income = income$Income,
    
    Treated_median   = apply(treated, 1, median, na.rm = TRUE),
    Treated_lower    = apply(treated, 1, quantile, probs = 0.025, na.rm = TRUE),
    Treated_upper    = apply(treated, 1, quantile, probs = 0.975, na.rm = TRUE),
    
    Untreated_median = apply(untreated, 1, median, na.rm = TRUE),
    Untreated_lower  = apply(untreated, 1, quantile, probs = 0.025, na.rm = TRUE),
    Untreated_upper  = apply(untreated, 1, quantile, probs = 0.975, na.rm = TRUE),
    
    Total_median     = apply(total, 1, median, na.rm = TRUE),
    Total_lower      = apply(total, 1, quantile, probs = 0.025, na.rm = TRUE),
    Total_upper      = apply(total, 1, quantile, probs = 0.975, na.rm = TRUE),
    
    Incidence_rate_median = apply(inc_rate_total, 1, median, na.rm = TRUE),
    Incidence_rate_lower  = apply(inc_rate_total, 1, quantile, probs = 0.025, na.rm = TRUE),
    Incidence_rate_upper  = apply(inc_rate_total, 1, quantile, probs = 0.975, na.rm = TRUE),
    
    Proportion_median = apply(prop_treated, 1, median, na.rm = TRUE),
    Proportion_lower  = apply(prop_treated, 1, quantile, probs = 0.025, na.rm = TRUE),
    Proportion_upper  = apply(prop_treated, 1, quantile, probs = 0.975, na.rm = TRUE)
  )
}

LRI_summary <- make_summary(Country_Year_meta, Antibiotics_LRI, treated_LRI_matrix, untreated_LRI_matrix, Incidence_LRI_matrix)
Fever_summary <- make_summary(Country_Year_meta, Antibiotics_Fever, treated_Fever_matrix, untreated_Fever_matrix, Incidence_Fever_matrix)

LRI_df <- make_df(Country_Year_meta, Antibiotics_LRI, treated_LRI_matrix, untreated_LRI_matrix, Incidence_LRI_matrix)
Fever_df <- make_df(Country_Year_meta, Antibiotics_Fever, treated_Fever_matrix, untreated_Fever_matrix, Incidence_Fever_matrix)

Fever_df_0 <- Fever_df[Fever_df$Malaria == 0,]
Fever_df_1 <- Fever_df[Fever_df$Malaria == 1,]

years_to_save <- 2003:2019

#### 1.1 CSV file

for (yr in years_to_save) {
  write.csv(LRI_summary[LRI_summary$Year == yr, ], paste0("outputs/table/by_country_year/probability/SA_LRI_", yr, ".csv"), row.names = FALSE)
  write.csv(Fever_summary[Fever_summary$Year == yr, ], paste0("outputs/table/by_country_year/probability/SA_Fever_", yr, ".csv"), row.names = FALSE)
}

for (yr in years_to_save) {
  write.csv(LRI_df[LRI_df$Year == yr, ], paste0("outputs/global_map_df/SA_1.1_Untreated_number_", yr, "_LRI.csv"), row.names = FALSE)
  write.csv(Fever_df_0[Fever_df_0$Year == yr, ], paste0("outputs/global_map_df/SA_1.1_Untreated_number_0_", yr, "_Fever.csv"), row.names = FALSE)
  write.csv(Fever_df_1[Fever_df_1$Year == yr, ], paste0("outputs/global_map_df/SA_1.1_Untreated_number_1_", yr, "_Fever.csv"), row.names = FALSE)
}







#----------------------------1.2 by year (treatment) ------------------------------

years <- 2003:2019
summary_list_LRI <- summary_list_LRI <- list()


#---LRI---
for (yr in years) {
  # Filter matrices by year
  treated_LRI <- cbind(Country_Year_meta, treated_LRI_matrix) %>%
    filter(Year == yr) %>%
    dplyr::select(-Country, -Year)
  
  untreated_LRI <- cbind(Country_Year_meta, untreated_LRI_matrix) %>%
    filter(Year == yr) %>%
    dplyr::select(-Country, -Year)
  
  total_LRI <- cbind(Country_Year_meta, Incidence_LRI_matrix) %>%
    filter(Year == yr) %>%
    dplyr::select(-Country, -Year)
  
  # Sum over countries for each simulation
  sum_treated_LRI   <- colSums(treated_LRI)
  sum_untreated_LRI <- colSums(untreated_LRI)
  sum_total_LRI     <- colSums(total_LRI)
  # --- population summed over countries (fixed, no sim dimension) ---
  sum_population_LRI <- sum(
    Antibiotics_LRI$Population[Country_Year_meta$Year == yr],
    na.rm = TRUE)
  incidence_LRI <- sum_total_LRI / sum_population_LRI
  # --- proportion treated ---
  proportion_treated_LRI <- sum_treated_LRI / 
    (sum_treated_LRI + sum_untreated_LRI)
  
  # Store results as long-format table
  summary_list_LRI[[as.character(yr)]] <- data.frame(
    Year = rep(yr, 15),
    Condition = rep(c("Treated LRI", "Untreated LRI", "Total LRI", "Incidence LRI", "Proportion LRI"), each = 3),
    Quantile = rep(c("Median", "2.5%", "97.5%"), times = 5),
    Value = c(
      quantile(sum_treated_LRI, probs = 0.50, na.rm = TRUE),
      quantile(sum_treated_LRI, probs = 0.025, na.rm = TRUE),
      quantile(sum_treated_LRI, probs = 0.975, na.rm = TRUE),
      
      quantile(sum_untreated_LRI, probs = 0.50, na.rm = TRUE),
      quantile(sum_untreated_LRI, probs = 0.025, na.rm = TRUE),
      quantile(sum_untreated_LRI, probs = 0.975, na.rm = TRUE),
      
      quantile(sum_total_LRI, probs = 0.50, na.rm = TRUE),
      quantile(sum_total_LRI, probs = 0.025, na.rm = TRUE),
      quantile(sum_total_LRI, probs = 0.975, na.rm = TRUE),
      
      quantile(incidence_LRI, probs = 0.50, na.rm = TRUE),
      quantile(incidence_LRI, probs = 0.025, na.rm = TRUE),
      quantile(incidence_LRI, probs = 0.975, na.rm = TRUE),
      
      quantile(proportion_treated_LRI, probs = 0.50, na.rm = TRUE),
      quantile(proportion_treated_LRI, probs = 0.025, na.rm = TRUE),
      quantile(proportion_treated_LRI, probs = 0.975, na.rm = TRUE)
    )
  )
}

# Combine all sublists
summary_all_years_LRI <- do.call(rbind, summary_list_LRI)

# Combine all years into one data frame
summary_wide_LRI <- summary_all_years_LRI %>%
  mutate(Column = paste0(gsub(" ", "_", Condition), "_", Quantile)) %>%
  dplyr::select(Year, Column, Value) %>%
  pivot_wider(names_from = Column, values_from = Value)

summary_wide_LRI_combined <- summary_wide_LRI %>%
  mutate(
    Treated_LRI = sprintf("%s (%s  %s)", 
                          formatC(Treated_LRI_Median, format = "f", digits = 0, big.mark = ","),
                          formatC(`Treated_LRI_2.5%`, format = "f", digits = 0, big.mark = ","),
                          formatC(`Treated_LRI_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Untreated_LRI = sprintf("%s (%s  %s)", 
                            formatC(Untreated_LRI_Median, format = "f", digits = 0, big.mark = ","),
                            formatC(`Untreated_LRI_2.5%`, format = "f", digits = 0, big.mark = ","),
                            formatC(`Untreated_LRI_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Total_LRI = sprintf("%s (%s  %s)", 
                        formatC(Total_LRI_Median, format = "f", digits = 0, big.mark = ","),
                        formatC(`Total_LRI_2.5%`, format = "f", digits = 0, big.mark = ","),
                        formatC(`Total_LRI_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Incidence_LRI = sprintf("%.2f (%.2f  %.2f)",
                            Incidence_LRI_Median,
                            `Incidence_LRI_2.5%`,
                            `Incidence_LRI_97.5%`
    ),
    
    Proportion_LRI = sprintf("%.3f (%.3f, %.3f)", 
                             Proportion_LRI_Median, 
                             `Proportion_LRI_2.5%`, 
                             `Proportion_LRI_97.5%`)
  ) %>%
  dplyr::select(Year, Treated_LRI, Untreated_LRI, Total_LRI, Incidence_LRI, Proportion_LRI)

LRI_year_proportion <- summary_wide_LRI %>%
  dplyr::select(Year, Proportion_LRI_Median, `Proportion_LRI_2.5%`, `Proportion_LRI_97.5%`) %>%
  dplyr::rename(Proportion_Median = Proportion_LRI_Median, Proportion_Lower = `Proportion_LRI_2.5%`, Proportion_Upper = `Proportion_LRI_97.5%`)




#---Fever---
years <- 2003:2019
summary_list_Fever <- list()

for (yr in years) {
  # Filter matrices by year
  treated_Fever <- cbind(Country_Year_meta, treated_Fever_matrix) %>%
    filter(Year == yr) %>%
    dplyr::select(-Country, -Year)
  
  untreated_Fever <- cbind(Country_Year_meta, untreated_Fever_matrix) %>%
    filter(Year == yr) %>%
    dplyr::select(-Country, -Year)
  
  total_Fever <- cbind(Country_Year_meta, Incidence_Fever_matrix) %>%
    filter(Year == yr) %>%
    dplyr::select(-Country, -Year)
  
  # Sum over countries for each simulation
  sum_treated_Fever <- colSums(treated_Fever)
  sum_untreated_Fever <- colSums(untreated_Fever)
  sum_total_Fever <- colSums(total_Fever)
  
  sum_population_Fever <- sum(
    Antibiotics_Fever$Population[Country_Year_meta$Year == yr],
    na.rm = TRUE)
  incidence_Fever <- sum_total_Fever / sum_population_Fever
  
  proportion_treated_Fever <- sum_treated_Fever / (sum_treated_Fever + sum_untreated_Fever)
  
  # Store results as long-format table
  summary_list_Fever[[as.character(yr)]] <- data.frame(
    Year = rep(yr, 15),
    Condition = rep(c("Treated Fever", "Untreated Fever", "Total Fever", "Incidence Fever","Proportion Fever"), each = 3),
    Quantile = rep(c("Median", "2.5%", "97.5%"), times = 5),
    Value = c(
      quantile(sum_treated_Fever, probs = 0.50, na.rm = TRUE),
      quantile(sum_treated_Fever, probs = 0.025, na.rm = TRUE),
      quantile(sum_treated_Fever, probs = 0.975, na.rm = TRUE),
      
      quantile(sum_untreated_Fever, probs = 0.50, na.rm = TRUE),
      quantile(sum_untreated_Fever, probs = 0.025, na.rm = TRUE),
      quantile(sum_untreated_Fever, probs = 0.975, na.rm = TRUE),
      
      quantile(sum_total_Fever, probs = 0.50, na.rm = TRUE),
      quantile(sum_total_Fever, probs = 0.025, na.rm = TRUE),
      quantile(sum_total_Fever, probs = 0.975, na.rm = TRUE),
      
      quantile(incidence_Fever, probs = 0.50, na.rm = TRUE),
      quantile(incidence_Fever, probs = 0.025, na.rm = TRUE),
      quantile(incidence_Fever, probs = 0.975, na.rm = TRUE),
      
      quantile(proportion_treated_Fever, probs = 0.50, na.rm = TRUE),
      quantile(proportion_treated_Fever, probs = 0.025, na.rm = TRUE),
      quantile(proportion_treated_Fever, probs = 0.975, na.rm = TRUE)
    )
  )
}

# Combine all sublists
summary_all_years_Fever <- do.call(rbind, summary_list_Fever)

# Combine all years into one data frame
summary_wide_Fever <- summary_all_years_Fever %>%
  mutate(Column = paste0(gsub(" ", "_", Condition), "_", Quantile)) %>%
  dplyr::select(Year, Column, Value) %>%
  pivot_wider(names_from = Column, values_from = Value)

summary_wide_Fever_combined <- summary_wide_Fever %>%
  mutate(
    Treated_Fever = sprintf("%s (%s  %s)", 
                            formatC(Treated_Fever_Median, format = "f", digits = 0, big.mark = ","),
                            formatC(`Treated_Fever_2.5%`, format = "f", digits = 0, big.mark = ","),
                            formatC(`Treated_Fever_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Untreated_Fever = sprintf("%s (%s  %s)", 
                              formatC(Untreated_Fever_Median, format = "f", digits = 0, big.mark = ","),
                              formatC(`Untreated_Fever_2.5%`, format = "f", digits = 0, big.mark = ","),
                              formatC(`Untreated_Fever_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Total_Fever = sprintf("%s (%s  %s)", 
                          formatC(Total_Fever_Median, format = "f", digits = 0, big.mark = ","),
                          formatC(`Total_Fever_2.5%`, format = "f", digits = 0, big.mark = ","),
                          formatC(`Total_Fever_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Incidence_Fever = sprintf("%.2f (%.2f  %.2f)",
                              Incidence_Fever_Median,
                              `Incidence_Fever_2.5%`,
                              `Incidence_Fever_97.5%`
    ),
    Proportion_Fever = sprintf("%.3f (%.3f, %.3f)", 
                               Proportion_Fever_Median, 
                               `Proportion_Fever_2.5%`, 
                               `Proportion_Fever_97.5%`)
  ) %>%
  dplyr::select(Year, Treated_Fever, Untreated_Fever, Total_Fever, Incidence_Fever, Proportion_Fever)


Fever_year_proportion <- summary_wide_Fever %>%
  dplyr::select(Year, Proportion_Fever_Median, `Proportion_Fever_2.5%`, `Proportion_Fever_97.5%`) %>%
  dplyr::rename(Proportion_Median = Proportion_Fever_Median, Proportion_Lower = `Proportion_Fever_2.5%`, Proportion_Upper = `Proportion_Fever_97.5%`)


colnames(summary_wide_LRI_combined) <- c("Year", "Treated_LRI", "Untreated_LRI", "Total_LRI", "Incidence_LRI", "Proportion_treated_LRI")
colnames(summary_wide_Fever_combined) <- c("Year", "Treated_Fever_without_LRI", "Untreated_Fever_without_LRI", "Total_Fever_without_LRI", "Incidence_Fever_without_LRI", "Proportion_treated_Fever_without_LRI")

write.csv(LRI_year_proportion   , "outputs/plot_df/SA_1.3_plot_Ab_LRI_year_proportion.csv")
write.csv(Fever_year_proportion   , "outputs/plot_df/SA_1.3_plot_Ab_Fever_year_proportion.csv")


#### 1.2 CSV file ####
write.csv(summary_wide_LRI_combined, "outputs/table/by_year/probability/SA_LRI.csv")
write.csv(summary_wide_Fever_combined, "outputs/table/by_year/probability/SA_Fever_without_LRI.csv")












#----------------------------1.3 by Malaria and year (treatment) ----------------------------

#------LRI-------
Country_Year_Malaria_meta <- Country_Year_meta %>% inner_join(AbLRI_matrix[,c(1,2,5)], by = c("Country", "Year"))
Malaria_level <- c(0,1)
years <- 2003:2019
summary_list_malaria_LRI <- list()
for (mal in Malaria_level) {
  for (yr in years) {
    
    treated_LRI <- cbind(Country_Year_Malaria_meta, treated_LRI_matrix) %>%
      filter(Malaria == mal, Year == yr) %>%
      dplyr::select(-Country, -Year, -Malaria)
    
    untreated_LRI <- cbind(Country_Year_Malaria_meta, untreated_LRI_matrix) %>%
      filter(Malaria == mal, Year == yr) %>%
      dplyr::select(-Country, -Year, -Malaria)
    
    total_LRI <- cbind(Country_Year_Malaria_meta, Incidence_LRI_matrix) %>%
      filter(Malaria == mal, Year == yr) %>%
      dplyr::select(-Country, -Year, -Malaria)
    
    # sums per simulation
    sum_treated   <- colSums(treated_LRI, na.rm = TRUE)
    sum_untreated <- colSums(untreated_LRI, na.rm = TRUE)
    sum_total     <- colSums(total_LRI, na.rm = TRUE)
    
    # population (fixed, stratified)
    sum_population_LRI <- sum(
      Antibiotics_LRI$Population[
        Antibiotics_LRI$Year == yr &
          Antibiotics_LRI$Malaria == mal
      ],
      na.rm = TRUE
    )
    
    # incidence per simulation
    incidence_LRI <- sum_total / sum_population_LRI
    
    # proportion treated
    proportion_treated <- sum_treated / (sum_treated + sum_untreated)
    
    summary_list_malaria_LRI[[paste(mal, yr, sep = "_")]] <- data.frame(
      Malaria = rep(mal, 15),
      Year = rep(yr, 15),
      Condition = rep(
        c("Treated LRI",
          "Untreated LRI",
          "Total LRI",
          "Incidence LRI",
          "Proportion LRI"),
        each = 3
      ),
      Quantile = rep(c("Median", "2.5%", "97.5%"), times = 5),
      Value = c(
        quantile(sum_treated, 0.50),
        quantile(sum_treated, 0.025),
        quantile(sum_treated, 0.975),
        
        quantile(sum_untreated, 0.50),
        quantile(sum_untreated, 0.025),
        quantile(sum_untreated, 0.975),
        
        quantile(sum_total, 0.50),
        quantile(sum_total, 0.025),
        quantile(sum_total, 0.975),
        
        quantile(incidence_LRI, 0.50),
        quantile(incidence_LRI, 0.025),
        quantile(incidence_LRI, 0.975),
        
        quantile(proportion_treated, 0.50),
        quantile(proportion_treated, 0.025),
        quantile(proportion_treated, 0.975)
      )
    )
  }
}

# Combine all sublists
summary_all_years_LRI <- do.call(rbind, summary_list_malaria_LRI)

# Convert to wide format
summary_wide_LRI <- summary_all_years_LRI %>%
  mutate(Column = paste0(gsub(" ", "_", Condition), "_", Quantile)) %>%
  dplyr::select(Year, Malaria, Column, Value) %>%
  pivot_wider(names_from = Column, values_from = Value)


summary_wide_LRI_by_malaria <- summary_wide_LRI %>%
  mutate(
    Treated_LRI = sprintf("%s (%s  %s)", 
                          formatC(Treated_LRI_Median, format = "f", digits = 0, big.mark = ","),
                          formatC(`Treated_LRI_2.5%`, format = "f", digits = 0, big.mark = ","),
                          formatC(`Treated_LRI_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Untreated_LRI = sprintf("%s (%s  %s)", 
                            formatC(Untreated_LRI_Median, format = "f", digits = 0, big.mark = ","),
                            formatC(`Untreated_LRI_2.5%`, format = "f", digits = 0, big.mark = ","),
                            formatC(`Untreated_LRI_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Total_LRI = sprintf("%s (%s  %s)", 
                        formatC(Total_LRI_Median, format = "f", digits = 0, big.mark = ","),
                        formatC(`Total_LRI_2.5%`, format = "f", digits = 0, big.mark = ","),
                        formatC(`Total_LRI_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Incidence_LRI = sprintf("%.2f (%.2f  %.2f)",
                            Incidence_LRI_Median,
                            `Incidence_LRI_2.5%`,
                            `Incidence_LRI_97.5%`
    ),
    
    Proportion_LRI = sprintf("%.3f (%.3f, %.3f)", 
                             Proportion_LRI_Median, 
                             `Proportion_LRI_2.5%`, 
                             `Proportion_LRI_97.5%`)
  ) %>%
  dplyr::select(Year, Malaria, Treated_LRI, Untreated_LRI, Total_LRI, Incidence_LRI, Proportion_LRI)


LRI_year_proportion_malaria_0 <- summary_wide_LRI %>%
  filter(Malaria == 0) %>%
  dplyr::select(Year, Malaria, Proportion_LRI_Median, `Proportion_LRI_2.5%`, `Proportion_LRI_97.5%`) %>%
  dplyr::rename(Proportion_Median = Proportion_LRI_Median, Proportion_Lower = `Proportion_LRI_2.5%`, Proportion_Upper = `Proportion_LRI_97.5%`)

LRI_year_proportion_malaria_1 <- summary_wide_LRI %>%
  filter(Malaria == 1) %>%
  dplyr::select(Year, Malaria, Proportion_LRI_Median, `Proportion_LRI_2.5%`, `Proportion_LRI_97.5%`) %>%
  dplyr::rename(Proportion_Median = Proportion_LRI_Median, Proportion_Lower = `Proportion_LRI_2.5%`, Proportion_Upper = `Proportion_LRI_97.5%`)








#------Fever-------
Country_Year_Malaria_meta <- Country_Year_meta %>% inner_join(AbFever_matrix[,c(1,2,5)], by = c("Country", "Year"))
Malaria_level <- c(0,1)
years <- 2003:2019
summary_list_malaria_Fever <- list()
for (mal in Malaria_level) {
  for (yr in years) {
    
    treated_Fever <- cbind(Country_Year_Malaria_meta, treated_Fever_matrix) %>%
      filter(Malaria == mal, Year == yr) %>%
      dplyr::select(-Country, -Year, -Malaria)
    
    untreated_Fever <- cbind(Country_Year_Malaria_meta, untreated_Fever_matrix) %>%
      filter(Malaria == mal, Year == yr) %>%
      dplyr::select(-Country, -Year, -Malaria)
    
    total_Fever <- cbind(Country_Year_Malaria_meta, Incidence_Fever_matrix) %>%
      filter(Malaria == mal, Year == yr) %>%
      dplyr::select(-Country, -Year, -Malaria)
    
    # sums per simulation
    sum_treated   <- colSums(treated_Fever, na.rm = TRUE)
    sum_untreated <- colSums(untreated_Fever, na.rm = TRUE)
    sum_total     <- colSums(total_Fever, na.rm = TRUE)
    
    # population (fixed, stratified)
    sum_population_Fever <- sum(
      Antibiotics_Fever$Population[
        Antibiotics_Fever$Year == yr &
          Antibiotics_Fever$Malaria == mal
      ],
      na.rm = TRUE
    )
    
    # incidence per simulation
    incidence_Fever <- sum_total / sum_population_Fever
    
    # proportion treated
    proportion_treated <- sum_treated / (sum_treated + sum_untreated)
    
    summary_list_malaria_Fever[[paste(mal, yr, sep = "_")]] <- data.frame(
      Malaria = rep(mal, 15),
      Year = rep(yr, 15),
      Condition = rep(
        c("Treated Fever",
          "Untreated Fever",
          "Total Fever",
          "Incidence Fever",
          "Proportion Fever"),
        each = 3
      ),
      Quantile = rep(c("Median", "2.5%", "97.5%"), times = 5),
      Value = c(
        quantile(sum_treated, 0.50),
        quantile(sum_treated, 0.025),
        quantile(sum_treated, 0.975),
        
        quantile(sum_untreated, 0.50),
        quantile(sum_untreated, 0.025),
        quantile(sum_untreated, 0.975),
        
        quantile(sum_total, 0.50),
        quantile(sum_total, 0.025),
        quantile(sum_total, 0.975),
        
        quantile(incidence_Fever, 0.50),
        quantile(incidence_Fever, 0.025),
        quantile(incidence_Fever, 0.975),
        
        quantile(proportion_treated, 0.50),
        quantile(proportion_treated, 0.025),
        quantile(proportion_treated, 0.975)
      )
    )
  }
}

# Combine all sublists
summary_all_years_Fever <- do.call(rbind, summary_list_malaria_Fever)

# Convert to wide format
summary_wide_Fever <- summary_all_years_Fever %>%
  mutate(Column = paste0(gsub(" ", "_", Condition), "_", Quantile)) %>%
  dplyr::select(Year, Malaria, Column, Value) %>%
  pivot_wider(names_from = Column, values_from = Value)


summary_wide_Fever_by_malaria <- summary_wide_Fever %>%
  mutate(
    Treated_Fever = sprintf("%s (%s  %s)", 
                            formatC(Treated_Fever_Median, format = "f", digits = 0, big.mark = ","),
                            formatC(`Treated_Fever_2.5%`, format = "f", digits = 0, big.mark = ","),
                            formatC(`Treated_Fever_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Untreated_Fever = sprintf("%s (%s  %s)", 
                              formatC(Untreated_Fever_Median, format = "f", digits = 0, big.mark = ","),
                              formatC(`Untreated_Fever_2.5%`, format = "f", digits = 0, big.mark = ","),
                              formatC(`Untreated_Fever_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Total_Fever = sprintf("%s (%s  %s)", 
                          formatC(Total_Fever_Median, format = "f", digits = 0, big.mark = ","),
                          formatC(`Total_Fever_2.5%`, format = "f", digits = 0, big.mark = ","),
                          formatC(`Total_Fever_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Incidence_Fever = sprintf("%.2f (%.2f  %.2f)",
                              Incidence_Fever_Median,
                              `Incidence_Fever_2.5%`,
                              `Incidence_Fever_97.5%`
    ),
    
    Proportion_Fever = sprintf("%.3f (%.3f, %.3f)", 
                               Proportion_Fever_Median, 
                               `Proportion_Fever_2.5%`, 
                               `Proportion_Fever_97.5%`)
  ) %>%
  dplyr::select(Year, Malaria, Treated_Fever, Untreated_Fever, Total_Fever, Incidence_Fever, Proportion_Fever)


Fever_year_proportion_malaria_0 <- summary_wide_Fever %>%
  filter(Malaria == 0) %>%
  dplyr::select(Year, Malaria, Proportion_Fever_Median, `Proportion_Fever_2.5%`, `Proportion_Fever_97.5%`) %>%
  dplyr::rename(Proportion_Median = Proportion_Fever_Median, Proportion_Lower = `Proportion_Fever_2.5%`, Proportion_Upper = `Proportion_Fever_97.5%`)

Fever_year_proportion_malaria_1 <- summary_wide_Fever %>%
  filter(Malaria == 1) %>%
  dplyr::select(Year, Malaria, Proportion_Fever_Median, `Proportion_Fever_2.5%`, `Proportion_Fever_97.5%`) %>%
  dplyr::rename(Proportion_Median = Proportion_Fever_Median, Proportion_Lower = `Proportion_Fever_2.5%`, Proportion_Upper = `Proportion_Fever_97.5%`)


colnames(summary_wide_LRI_by_malaria) <- c("Year", "Malaria", "Treated_LRI", "Untreated_LRI", "Total_LRI", "Incidence_LRI", "Proportion_treated_LRI")
colnames(summary_wide_Fever_by_malaria) <- c("Year", "Malaria", "Treated_Fever_without_LRI", "Untreated_Fever_without_LRI", "Total_Fever_without_LRI", "Incidence_Fever", "Proportion_treated_Fever_without_LRI")


# #### 1.4 CSV file ####
write.csv(summary_wide_LRI_by_malaria, "outputs/table/by_year_malaria/probability/SA_LRI.csv")
write.csv(summary_wide_Fever_by_malaria, "outputs/table/by_year_malaria/probability/SA_Fever_wihtout_LRI.csv")
Malaria_distribution <- as.matrix(with(Country_Year_Malaria_meta, addmargins(table(Year, Malaria), margin = c(1,2))))
write.csv(Malaria_distribution, "outputs/table/by_year_malaria/SA_Country_number_by_Malaria_Year.csv")


# save the plot data
write.csv(LRI_year_proportion_malaria_0  , "outputs/plot_df/SA_1.4_plot_Ab_LRI_year_proportion_malaria_0.csv")
write.csv(LRI_year_proportion_malaria_1  , "outputs/plot_df/SA_1.4_plot_Ab_LRI_year_proportion_malaria_1.csv")
write.csv(Fever_year_proportion_malaria_0, "outputs/plot_df/SA_1.4_plot_Ab_Fever_year_proportion_malaria_0.csv")
write.csv(Fever_year_proportion_malaria_1, "outputs/plot_df/SA_1.4_plot_Ab_Fever_year_proportion_malaria_1.csv")































#----------------------------1.4 by Income and year (treatment) ----------------------------



first_sim <- AbLRI_matrix[,c(1,2,4)]
income_dist <- pivot_wider(first_sim, names_from = Year, values_from = Income)


first5_years <- as.character(2003:2007)
last5_years  <- as.character(2015:2019)

income_dist_summary <- income_dist %>%
  rowwise() %>%
  mutate(
    First5 = {
      first_vals <- c_across(all_of(first5_years))
      tab <- table(first_vals)
      if(length(tab) == 0) NA_character_ else names(sort(tab, decreasing = TRUE))[1]
    },
    Last5 = {
      last_vals <- c_across(all_of(last5_years))
      tab <- table(last_vals)
      if(length(tab) == 0) NA_character_ else names(sort(tab, decreasing = TRUE))[1]
    }
  ) %>%
  ungroup()

income_dist_summary[which_na(income_dist_summary$First5),]
income_dist_summary[income_dist_summary$Country == "Djibouti",]$First5 <- "LM" 
AbLRI_matrix_2 <- AbLRI_matrix %>% inner_join(income_dist_summary[,c("Country","First5","Last5")], by = "Country") 

AbFever_matrix_2 <- AbFever_matrix %>% inner_join(income_dist_summary[,c("Country","First5","Last5")], by = "Country") 
rm(income_dist, income_dist_summary, first_sim, first5_years, last5_years)

Country_Year_info <- Country_Year_meta %>% inner_join(AbLRI_matrix_2[,c(1:7,10008:10009)], by = c("Country", "Year"))
Country_Year_info <- Country_Year_info[,c("Country","Year","GDP","Last5","Malaria","GBD_region","Population")]
Country_Year_info_2 <- Country_Year_info %>% dplyr::select(Country, Last5, Malaria, GBD_region, Population) %>% filter(!duplicated(Country))
write.csv(Country_Year_info, "outputs/draft_df/SA_Suppl_3_108LMICs.csv")
write.csv(Country_Year_info_2, "outputs/draft_df/SA_Suppl_3_2_108LMICs.csv")


#------LRI-------
Country_Year_Income_meta <- Country_Year_meta %>% inner_join(AbLRI_matrix_2[,c("Country","Year","Last5")], by = c("Country", "Year"))
Antibiotics_LRI_2 <- Antibiotics_LRI %>% inner_join(AbLRI_matrix_2[,c("Country","Year","Last5")], by = c("Country", "Year"))
Income_level <- c('L','LM','UM')
years <- 2003:2019
summary_list_LRI <- list()

for (inc in Income_level) {
  for (yr in years) {
    # Filter and select simulation matrices
    treated_LRI <- cbind(Country_Year_Income_meta, treated_LRI_matrix) %>%
      filter(Last5 == inc, Year == yr) %>%
      dplyr::select(-Country, -Year, -Last5)
    
    untreated_LRI <- cbind(Country_Year_Income_meta, untreated_LRI_matrix) %>%
      filter(Last5 == inc, Year == yr) %>%
      dplyr::select(-Country, -Year, -Last5)
    
    total_LRI <- cbind(Country_Year_Income_meta, Incidence_LRI_matrix) %>%
      filter(Last5 == inc, Year == yr) %>%
      dplyr::select(-Country, -Year, -Last5)
    
    # Sum over countries for each simulation
    sum_treated <- colSums(treated_LRI, na.rm = TRUE)
    sum_untreated <- colSums(untreated_LRI, na.rm = TRUE)
    sum_total <- colSums(total_LRI, na.rm = TRUE)
    # population (fixed, stratified)
    sum_population_LRI <- sum(
      Antibiotics_LRI_2$Population[
        Antibiotics_LRI_2$Year == yr &
          Antibiotics_LRI_2$Last5 == inc
      ],
      na.rm = TRUE
    )
    
    # incidence per simulation
    incidence_LRI <- sum_total / sum_population_LRI
    
    # proportion treated
    proportion_treated <- sum_treated / (sum_treated + sum_untreated)
    
    # Store results as long-format table
    summary_list_LRI[[paste(inc, yr, sep = "_")]] <- data.frame(
      Last5 = inc,
      Year = rep(yr, 15),
      Condition = rep(c("Treated LRI", "Untreated LRI", "Total LRI", "Incidence LRI", "Proportion LRI"), each = 3),
      Quantile = rep(c("Median", "2.5%", "97.5%"), times = 5),
      Value = c(
        quantile(sum_treated, probs = 0.50, na.rm = TRUE),
        quantile(sum_treated, probs = 0.025, na.rm = TRUE),
        quantile(sum_treated, probs = 0.975, na.rm = TRUE),
        
        quantile(sum_untreated, probs = 0.50, na.rm = TRUE),
        quantile(sum_untreated, probs = 0.025, na.rm = TRUE),
        quantile(sum_untreated, probs = 0.975, na.rm = TRUE),
        
        quantile(sum_total, probs = 0.50, na.rm = TRUE),
        quantile(sum_total, probs = 0.025, na.rm = TRUE),
        quantile(sum_total, probs = 0.975, na.rm = TRUE),
        
        quantile(incidence_LRI, 0.50),
        quantile(incidence_LRI, 0.025),
        quantile(incidence_LRI, 0.975),
        
        quantile(proportion_treated, probs = 0.50, na.rm = TRUE),
        quantile(proportion_treated, probs = 0.025, na.rm = TRUE),
        quantile(proportion_treated, probs = 0.975, na.rm = TRUE)
      )
    )
  }
}

# Combine all sublists
summary_all_years_LRI <- do.call(rbind, summary_list_LRI)

# Convert to wide format
summary_wide_LRI <- summary_all_years_LRI %>%
  mutate(Column = paste0(gsub(" ", "_", Condition), "_", Quantile)) %>%
  dplyr::select(Year, Last5, Column, Value) %>%
  pivot_wider(names_from = Column, values_from = Value)


summary_wide_LRI_by_income <- summary_wide_LRI %>%
  mutate(
    Treated_LRI = sprintf("%s (%s  %s)", 
                          formatC(Treated_LRI_Median, format = "f", digits = 0, big.mark = ","),
                          formatC(`Treated_LRI_2.5%`, format = "f", digits = 0, big.mark = ","),
                          formatC(`Treated_LRI_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Untreated_LRI = sprintf("%s (%s  %s)", 
                            formatC(Untreated_LRI_Median, format = "f", digits = 0, big.mark = ","),
                            formatC(`Untreated_LRI_2.5%`, format = "f", digits = 0, big.mark = ","),
                            formatC(`Untreated_LRI_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Total_LRI = sprintf("%s (%s  %s)", 
                        formatC(Total_LRI_Median, format = "f", digits = 0, big.mark = ","),
                        formatC(`Total_LRI_2.5%`, format = "f", digits = 0, big.mark = ","),
                        formatC(`Total_LRI_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Incidence_LRI = sprintf("%.2f (%.2f  %.2f)",
                            Incidence_LRI_Median,
                            `Incidence_LRI_2.5%`,
                            `Incidence_LRI_97.5%`),
    
    Proportion_LRI = sprintf("%.3f (%.3f, %.3f)", 
                             Proportion_LRI_Median, 
                             `Proportion_LRI_2.5%`, 
                             `Proportion_LRI_97.5%`)
  ) %>%
  dplyr::select(Year, Last5, Treated_LRI, Untreated_LRI, Total_LRI, Incidence_LRI, Proportion_LRI)


LRI_year_proportion_L <- summary_wide_LRI %>%
  filter(Last5 == 'L') %>%
  dplyr::select(Year, Last5, Proportion_LRI_Median, `Proportion_LRI_2.5%`, `Proportion_LRI_97.5%`) %>%
  dplyr::rename(Proportion_Median = Proportion_LRI_Median, Proportion_Lower = `Proportion_LRI_2.5%`, Proportion_Upper = `Proportion_LRI_97.5%`)

LRI_year_proportion_LM <- summary_wide_LRI %>%
  filter(Last5 == 'LM') %>%
  dplyr::select(Year, Last5, Proportion_LRI_Median, `Proportion_LRI_2.5%`, `Proportion_LRI_97.5%`) %>%
  dplyr::rename(Proportion_Median = Proportion_LRI_Median, Proportion_Lower = `Proportion_LRI_2.5%`, Proportion_Upper = `Proportion_LRI_97.5%`)

LRI_year_proportion_UM <- summary_wide_LRI %>%
  filter(Last5 == 'UM') %>%
  dplyr::select(Year, Last5, Proportion_LRI_Median, `Proportion_LRI_2.5%`, `Proportion_LRI_97.5%`) %>%
  dplyr::rename(Proportion_Median = Proportion_LRI_Median, Proportion_Lower = `Proportion_LRI_2.5%`, Proportion_Upper = `Proportion_LRI_97.5%`)









#------Fever-------
Country_Year_Income_meta <- Country_Year_meta %>% inner_join(AbFever_matrix_2[,c("Country","Year","Malaria","Last5")], by = c("Country", "Year"))
Antibiotics_Fever_2 <- Antibiotics_Fever %>% inner_join(AbFever_matrix_2[,c("Country","Year","Last5")], by = c("Country", "Year"))
Income_level <- c('L','LM','UM')
Malaria_level <- c(0,1)
years <- 2003:2019
summary_list_Fever <- list()

for (mal in Malaria_level) {
  for (inc in Income_level) {
    for (yr in years) {
      # Filter and select simulation matrices
      treated_Fever <- cbind(Country_Year_Income_meta, treated_Fever_matrix) %>%
        filter(Malaria == mal, Last5 == inc, Year == yr) %>%
        dplyr::select(-Country, -Year, -Malaria, -Last5)
      
      untreated_Fever <- cbind(Country_Year_Income_meta, untreated_Fever_matrix) %>%
        filter(Malaria == mal, Last5 == inc, Year == yr) %>%
        dplyr::select(-Country, -Year, -Malaria, -Last5)
      
      total_Fever <- cbind(Country_Year_Income_meta, Incidence_Fever_matrix) %>%
        filter(Malaria == mal, Last5 == inc, Year == yr) %>%
        dplyr::select(-Country, -Year, -Malaria, -Last5)
      
      # Sum over countries for each simulation
      sum_treated <- colSums(treated_Fever, na.rm = TRUE)
      sum_untreated <- colSums(untreated_Fever, na.rm = TRUE)
      sum_total <- colSums(total_Fever, na.rm = TRUE)
      # population (fixed, stratified)
      sum_population_Fever <- sum(
        Antibiotics_Fever_2$Population[
          Antibiotics_Fever_2$Year == yr &
            Antibiotics_Fever_2$Malaria == mal &
            Antibiotics_Fever_2$Last5 == inc
        ],
        na.rm = TRUE
      )
      
      # incidence per simulation
      incidence_Fever <- sum_total / sum_population_Fever
      
      # proportion treated
      proportion_treated <- sum_treated / (sum_treated + sum_untreated)
      
      # Store results as long-format table
      summary_list_Fever[[paste(mal, inc, yr, sep = "_")]] <- data.frame(
        Malaria = mal,
        Last5 = inc,
        Year = rep(yr, 15),
        Condition = rep(c("Treated Fever", "Untreated Fever", "Total Fever", "Incidence Fever", "Proportion Fever"), each = 3),
        Quantile = rep(c("Median", "2.5%", "97.5%"), times = 5),
        Value = c(
          quantile(sum_treated, probs = 0.50, na.rm = TRUE),
          quantile(sum_treated, probs = 0.025, na.rm = TRUE),
          quantile(sum_treated, probs = 0.975, na.rm = TRUE),
          
          quantile(sum_untreated, probs = 0.50, na.rm = TRUE),
          quantile(sum_untreated, probs = 0.025, na.rm = TRUE),
          quantile(sum_untreated, probs = 0.975, na.rm = TRUE),
          
          quantile(sum_total, probs = 0.50, na.rm = TRUE),
          quantile(sum_total, probs = 0.025, na.rm = TRUE),
          quantile(sum_total, probs = 0.975, na.rm = TRUE),
          
          quantile(incidence_Fever, 0.50),
          quantile(incidence_Fever, 0.025),
          quantile(incidence_Fever, 0.975),
          
          quantile(proportion_treated, probs = 0.50, na.rm = TRUE),
          quantile(proportion_treated, probs = 0.025, na.rm = TRUE),
          quantile(proportion_treated, probs = 0.975, na.rm = TRUE)
        )
      )
    }
  }
}  

# Combine all sublists
summary_all_years_Fever <- do.call(rbind, summary_list_Fever)

# Convert to wide format
summary_wide_Fever <- summary_all_years_Fever %>%
  mutate(Column = paste0(gsub(" ", "_", Condition), "_", Quantile)) %>%
  dplyr::select(Year, Malaria, Last5, Column, Value) %>%
  pivot_wider(names_from = Column, values_from = Value)


summary_wide_Fever_by_income <- summary_wide_Fever %>%
  mutate(
    Treated_Fever = sprintf("%s (%s  %s)", 
                            formatC(Treated_Fever_Median, format = "f", digits = 0, big.mark = ","),
                            formatC(`Treated_Fever_2.5%`, format = "f", digits = 0, big.mark = ","),
                            formatC(`Treated_Fever_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Untreated_Fever = sprintf("%s (%s  %s)", 
                              formatC(Untreated_Fever_Median, format = "f", digits = 0, big.mark = ","),
                              formatC(`Untreated_Fever_2.5%`, format = "f", digits = 0, big.mark = ","),
                              formatC(`Untreated_Fever_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Total_Fever = sprintf("%s (%s  %s)", 
                          formatC(Total_Fever_Median, format = "f", digits = 0, big.mark = ","),
                          formatC(`Total_Fever_2.5%`, format = "f", digits = 0, big.mark = ","),
                          formatC(`Total_Fever_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Incidence_Fever = sprintf("%.2f (%.2f  %.2f)",
                              Incidence_Fever_Median,
                              `Incidence_Fever_2.5%`,
                              `Incidence_Fever_97.5%`),
    
    Proportion_Fever = sprintf("%.3f (%.3f, %.3f)", 
                               Proportion_Fever_Median, 
                               `Proportion_Fever_2.5%`, 
                               `Proportion_Fever_97.5%`)
  ) %>%
  dplyr::select(Year, Malaria, Last5, Treated_Fever, Untreated_Fever, Total_Fever, Incidence_Fever, Proportion_Fever)


Fever_year_proportion_mal0_L <- summary_wide_Fever %>%
  filter(Malaria == 0, Last5 == 'L') %>%  
  
  dplyr::select(Year, Last5, Proportion_Fever_Median, `Proportion_Fever_2.5%`, `Proportion_Fever_97.5%`) %>%
  dplyr::rename(Proportion_Median = Proportion_Fever_Median, Proportion_Lower = `Proportion_Fever_2.5%`, Proportion_Upper = `Proportion_Fever_97.5%`)

Fever_year_proportion_mal0_LM <- summary_wide_Fever %>%
  filter(Malaria == 0, Last5 == 'LM') %>%
  dplyr::select(Year, Last5, Proportion_Fever_Median, `Proportion_Fever_2.5%`, `Proportion_Fever_97.5%`) %>%
  dplyr::rename(Proportion_Median = Proportion_Fever_Median, Proportion_Lower = `Proportion_Fever_2.5%`, Proportion_Upper = `Proportion_Fever_97.5%`)

Fever_year_proportion_mal0_UM <- summary_wide_Fever %>%
  filter(Malaria == 0, Last5 == 'UM') %>%
  dplyr::select(Year, Last5, Proportion_Fever_Median, `Proportion_Fever_2.5%`, `Proportion_Fever_97.5%`) %>%
  dplyr::rename(Proportion_Median = Proportion_Fever_Median, Proportion_Lower = `Proportion_Fever_2.5%`, Proportion_Upper = `Proportion_Fever_97.5%`)


Fever_year_proportion_mal1_L <- summary_wide_Fever %>%
  filter(Malaria == 1, Last5 == 'L') %>%
  dplyr::select(Year, Last5, Proportion_Fever_Median, `Proportion_Fever_2.5%`, `Proportion_Fever_97.5%`) %>%
  dplyr::rename(Proportion_Median = Proportion_Fever_Median, Proportion_Lower = `Proportion_Fever_2.5%`, Proportion_Upper = `Proportion_Fever_97.5%`)

Fever_year_proportion_mal1_LM <- summary_wide_Fever %>%
  filter(Malaria == 1, Last5 == 'LM') %>%
  dplyr::select(Year, Last5, Proportion_Fever_Median, `Proportion_Fever_2.5%`, `Proportion_Fever_97.5%`) %>%
  dplyr::rename(Proportion_Median = Proportion_Fever_Median, Proportion_Lower = `Proportion_Fever_2.5%`, Proportion_Upper = `Proportion_Fever_97.5%`)

Fever_year_proportion_mal1_UM <- summary_wide_Fever %>%
  filter(Malaria == 1, Last5 == 'UM') %>%
  dplyr::select(Year, Last5, Proportion_Fever_Median, `Proportion_Fever_2.5%`, `Proportion_Fever_97.5%`) %>%
  dplyr::rename(Proportion_Median = Proportion_Fever_Median, Proportion_Lower = `Proportion_Fever_2.5%`, Proportion_Upper = `Proportion_Fever_97.5%`)



# #### 1.3 CSV file ####
write.csv(summary_wide_LRI_by_income, "outputs/table/by_year_income/probability/SA_LRI.csv")
write.csv(summary_wide_Fever_by_income, "outputs/table/by_year_income/probability/SA_Fever_wihtout_LRI.csv")
Income_distribution <- as.matrix(with(Country_Year_Income_meta, addmargins(table(Year, Last5), margin = c(1,2))))
write.csv(Income_distribution, "outputs/table/by_year_income/SA_Country_number_by_Income_Year.csv")



# save the plot data
write.csv(LRI_year_proportion_L   , "outputs/plot_df/SA_1.3_plot_Ab_LRI_year_proportion_L.csv")
write.csv(LRI_year_proportion_LM  , "outputs/plot_df/SA_1.3_plot_Ab_LRI_year_proportion_LM.csv")
write.csv(LRI_year_proportion_UM  , "outputs/plot_df/SA_1.3_plot_Ab_LRI_year_proportion_UM.csv")
write.csv(Fever_year_proportion_mal0_L , "outputs/plot_df/SA_1.3_plot_Ab_Fever_year_proportion_mal0_L.csv")
write.csv(Fever_year_proportion_mal0_LM, "outputs/plot_df/SA_1.3_plot_Ab_Fever_year_proportion_mal0_LM.csv")
write.csv(Fever_year_proportion_mal0_UM, "outputs/plot_df/SA_1.3_plot_Ab_Fever_year_proportion_mal0_UM.csv")
write.csv(Fever_year_proportion_mal1_L , "outputs/plot_df/SA_1.3_plot_Ab_Fever_year_proportion_mal1_L.csv")
write.csv(Fever_year_proportion_mal1_LM, "outputs/plot_df/SA_1.3_plot_Ab_Fever_year_proportion_mal1_LM.csv")
write.csv(Fever_year_proportion_mal1_UM, "outputs/plot_df/SA_1.3_plot_Ab_Fever_year_proportion_mal1_UM.csv")














#----------------------------1.5 GBD_region (treatment) ----------------------------

#------LRI-------
Country_Year_GBD_region_meta <- Country_Year_meta %>% inner_join(AbLRI_matrix[,c(1,2,6)], by = c("Country", "Year"))
GBD_region_level <- c("North Africa and Middle East", "Central Europe, Eastern Europe, and Central Asia", "Sub-Saharan Africa", "High-income", "South Asia", "Latin America and Caribbean", "Southeast Asia, East Asia, and Oceania")

years <- 2003:2019
summary_list_GBD_region_LRI <- list()

for (gbd in GBD_region_level) {
  for (yr in years) {
    # Filter and select simulation matrices
    treated_LRI <- cbind(Country_Year_GBD_region_meta, treated_LRI_matrix) %>%
      filter(GBD_region == gbd, Year == yr) %>%
      dplyr::select(-Country, -Year, -GBD_region)
    
    untreated_LRI <- cbind(Country_Year_GBD_region_meta, untreated_LRI_matrix) %>%
      filter(GBD_region == gbd, Year == yr) %>%
      dplyr::select(-Country, -Year, -GBD_region)
    
    total_LRI <- cbind(Country_Year_GBD_region_meta, Incidence_LRI_matrix) %>%
      filter(GBD_region == gbd, Year == yr) %>%
      dplyr::select(-Country, -Year, -GBD_region)
    
    # Sum over countries for each simulation
    sum_treated <- colSums(treated_LRI, na.rm = TRUE)
    sum_untreated <- colSums(untreated_LRI, na.rm = TRUE)
    sum_total <- colSums(total_LRI, na.rm = TRUE)
    # population (fixed, stratified)
    sum_population_LRI <- sum(
      Antibiotics_LRI$Population[
        Antibiotics_LRI$Year == yr &
          Antibiotics_LRI$GBD_region == gbd
      ],
      na.rm = TRUE
    )
    
    # incidence per simulation
    incidence_LRI <- sum_total / sum_population_LRI
    
    proportion_treated <- sum_treated / (sum_treated + sum_untreated)
    
    # Store results as long-format table
    summary_list_GBD_region_LRI[[paste(gbd, yr, sep = "_")]] <- data.frame(
      GBD_region = gbd,
      Year = rep(yr, 15),
      Condition = rep(c("Treated LRI", "Untreated LRI", "Total LRI", "Incidence LRI", "Proportion LRI"), each = 3),
      Quantile = rep(c("Median", "2.5%", "97.5%"), times = 5),
      Value = c(
        quantile(sum_treated, probs = 0.50, na.rm = TRUE),
        quantile(sum_treated, probs = 0.025, na.rm = TRUE),
        quantile(sum_treated, probs = 0.975, na.rm = TRUE),
        
        quantile(sum_untreated, probs = 0.50, na.rm = TRUE),
        quantile(sum_untreated, probs = 0.025, na.rm = TRUE),
        quantile(sum_untreated, probs = 0.975, na.rm = TRUE),
        
        quantile(sum_total, probs = 0.50, na.rm = TRUE),
        quantile(sum_total, probs = 0.025, na.rm = TRUE),
        quantile(sum_total, probs = 0.975, na.rm = TRUE),
        
        quantile(incidence_LRI, 0.50),
        quantile(incidence_LRI, 0.025),
        quantile(incidence_LRI, 0.975),
        
        quantile(proportion_treated, probs = 0.50, na.rm = TRUE),
        quantile(proportion_treated, probs = 0.025, na.rm = TRUE),
        quantile(proportion_treated, probs = 0.975, na.rm = TRUE)
      )
    )
  }
}

# Combine all sublists
summary_all_years_LRI <- do.call(rbind, summary_list_GBD_region_LRI)

# Convert to wide format
summary_wide_LRI <- summary_all_years_LRI %>%
  mutate(Column = paste0(gsub(" ", "_", Condition), "_", Quantile)) %>%
  dplyr::select(Year, GBD_region, Column, Value) %>%
  pivot_wider(names_from = Column, values_from = Value)


summary_wide_LRI_by_GBD_region <- summary_wide_LRI %>%
  mutate(
    Treated_LRI = sprintf("%s (%s  %s)", 
                          formatC(Treated_LRI_Median, format = "f", digits = 0, big.mark = ","),
                          formatC(`Treated_LRI_2.5%`, format = "f", digits = 0, big.mark = ","),
                          formatC(`Treated_LRI_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Untreated_LRI = sprintf("%s (%s  %s)", 
                            formatC(Untreated_LRI_Median, format = "f", digits = 0, big.mark = ","),
                            formatC(`Untreated_LRI_2.5%`, format = "f", digits = 0, big.mark = ","),
                            formatC(`Untreated_LRI_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Total_LRI = sprintf("%s (%s  %s)", 
                        formatC(Total_LRI_Median, format = "f", digits = 0, big.mark = ","),
                        formatC(`Total_LRI_2.5%`, format = "f", digits = 0, big.mark = ","),
                        formatC(`Total_LRI_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Incidence_LRI = sprintf("%.2f (%.2f  %.2f)",
                            Incidence_LRI_Median,
                            `Incidence_LRI_2.5%`,
                            `Incidence_LRI_97.5%`),
    
    Proportion_LRI = sprintf("%.3f (%.3f, %.3f)", 
                             Proportion_LRI_Median, 
                             `Proportion_LRI_2.5%`, 
                             `Proportion_LRI_97.5%`)
  ) %>%
  dplyr::select(Year, GBD_region, Treated_LRI, Untreated_LRI, Total_LRI, Incidence_LRI, Proportion_LRI)


for (i in 1:length(GBD_region_level)) {
  temp <- summary_wide_LRI %>%
    filter(GBD_region == GBD_region_level[i]) %>%
    dplyr::select(Year, GBD_region, Proportion_LRI_Median, `Proportion_LRI_2.5%`, `Proportion_LRI_97.5%`) %>%
    rename(
      Proportion_Median = Proportion_LRI_Median,
      Proportion_Lower = `Proportion_LRI_2.5%`,
      Proportion_Upper = `Proportion_LRI_97.5%`
    )
  obj_name <- paste0("LRI_year_proportion_GBD_region_", i)
  assign(obj_name, temp)
}

for (i in 1:7) {
  obj_name <- paste0("LRI_year_proportion_GBD_region_", i)
  path_LRI <- paste0("outputs/plot_df/SA_1.5_plot_Ab_LRI_year_proportion_GBD_region_", i,".csv")
  
  write.csv(get(obj_name), path_LRI, row.names = FALSE)
}






#------Fever-------
Country_Year_GBD_region_meta <- Country_Year_meta %>% inner_join(AbFever_matrix[,c("Country","Year","Malaria","GBD_region")], by = c("Country", "Year"))
GBD_region_level <- c("North Africa and Middle East", "Central Europe, Eastern Europe, and Central Asia", "Sub-Saharan Africa", "High-income", "South Asia", "Latin America and Caribbean", "Southeast Asia, East Asia, and Oceania")
Malaria_level <- c(0,1)
years <- 2003:2019
summary_list_GBD_region_Fever <- list()

for (mal in Malaria_level) {
  for (gbd in GBD_region_level) {
    for (yr in years) {
      # Filter and select simulation matrices
      treated_Fever <- cbind(Country_Year_GBD_region_meta, treated_Fever_matrix) %>%
        filter(Malaria == mal, GBD_region == gbd, Year == yr) %>%
        dplyr::select(-Country, -Year, -Malaria, -GBD_region)
      
      untreated_Fever <- cbind(Country_Year_GBD_region_meta, untreated_Fever_matrix) %>%
        filter(Malaria == mal, GBD_region == gbd, Year == yr) %>%
        dplyr::select(-Country, -Year, -Malaria, -GBD_region)
      
      total_Fever <- cbind(Country_Year_GBD_region_meta, Incidence_Fever_matrix) %>%
        filter(Malaria == mal, GBD_region == gbd, Year == yr) %>%
        dplyr::select(-Country, -Year, -Malaria, -GBD_region)
      
      # Sum over countries for each simulation
      sum_treated <- colSums(treated_Fever, na.rm = TRUE)
      sum_untreated <- colSums(untreated_Fever, na.rm = TRUE)
      sum_total <- colSums(total_Fever, na.rm = TRUE)
      # population (fixed, stratified)
      sum_population_Fever <- sum(
        Antibiotics_Fever$Population[
          Antibiotics_Fever$Malaria == mal &
            Antibiotics_Fever$Year == yr &
            Antibiotics_Fever$GBD_region == gbd
        ],
        na.rm = TRUE
      )
      
      # incidence per simulation
      incidence_Fever <- sum_total / sum_population_Fever
      
      proportion_treated <- sum_treated / (sum_treated + sum_untreated)
      
      # Store results as long-format table
      summary_list_GBD_region_Fever[[paste(mal, gbd, yr, sep = "_")]] <- data.frame(
        Malaria = mal,
        GBD_region = gbd,
        Year = rep(yr, 15),
        Condition = rep(c("Treated Fever", "Untreated Fever", "Total Fever", "Incidence Fever", "Proportion Fever"), each = 3),
        Quantile = rep(c("Median", "2.5%", "97.5%"), times = 5),
        Value = c(
          quantile(sum_treated, probs = 0.50, na.rm = TRUE),
          quantile(sum_treated, probs = 0.025, na.rm = TRUE),
          quantile(sum_treated, probs = 0.975, na.rm = TRUE),
          
          quantile(sum_untreated, probs = 0.50, na.rm = TRUE),
          quantile(sum_untreated, probs = 0.025, na.rm = TRUE),
          quantile(sum_untreated, probs = 0.975, na.rm = TRUE),
          
          quantile(sum_total, probs = 0.50, na.rm = TRUE),
          quantile(sum_total, probs = 0.025, na.rm = TRUE),
          quantile(sum_total, probs = 0.975, na.rm = TRUE),
          
          quantile(incidence_Fever, 0.50, na.rm = TRUE),
          quantile(incidence_Fever, 0.025, na.rm = TRUE),
          quantile(incidence_Fever, 0.975, na.rm = TRUE),
          
          quantile(proportion_treated, probs = 0.50, na.rm = TRUE),
          quantile(proportion_treated, probs = 0.025, na.rm = TRUE),
          quantile(proportion_treated, probs = 0.975, na.rm = TRUE)
        )
      )
    }
  }
}  

# Combine all sublists
summary_all_years_Fever <- do.call(rbind, summary_list_GBD_region_Fever)

# Convert to wide format
summary_wide_Fever <- summary_all_years_Fever %>%
  mutate(Column = paste0(gsub(" ", "_", Condition), "_", Quantile)) %>%
  dplyr::select(Year, Malaria, GBD_region, Column, Value) %>%
  pivot_wider(names_from = Column, values_from = Value)


summary_wide_Fever_by_GBD_region <- summary_wide_Fever %>%
  mutate(
    Treated_Fever = sprintf("%s (%s  %s)", 
                            formatC(Treated_Fever_Median, format = "f", digits = 0, big.mark = ","),
                            formatC(`Treated_Fever_2.5%`, format = "f", digits = 0, big.mark = ","),
                            formatC(`Treated_Fever_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Untreated_Fever = sprintf("%s (%s  %s)", 
                              formatC(Untreated_Fever_Median, format = "f", digits = 0, big.mark = ","),
                              formatC(`Untreated_Fever_2.5%`, format = "f", digits = 0, big.mark = ","),
                              formatC(`Untreated_Fever_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Total_Fever = sprintf("%s (%s  %s)", 
                          formatC(Total_Fever_Median, format = "f", digits = 0, big.mark = ","),
                          formatC(`Total_Fever_2.5%`, format = "f", digits = 0, big.mark = ","),
                          formatC(`Total_Fever_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Incidence_Fever = sprintf("%.2f (%.2f  %.2f)",
                              Incidence_Fever_Median,
                              `Incidence_Fever_2.5%`,
                              `Incidence_Fever_97.5%`),
    
    Proportion_Fever = sprintf("%.3f (%.3f, %.3f)", 
                               Proportion_Fever_Median, 
                               `Proportion_Fever_2.5%`, 
                               `Proportion_Fever_97.5%`)
  ) %>%
  dplyr::select(Year, Malaria, GBD_region, Treated_Fever, Untreated_Fever, Total_Fever, Incidence_Fever, Proportion_Fever)




write.csv(summary_wide_LRI_by_GBD_region, "outputs/table/by_year_GBD_region/probability/SA_LRI.csv")
write.csv(summary_wide_Fever_by_GBD_region, "outputs/table/by_year_GBD_region/probability/SA_Fever_wihtout_LRI.csv")


for (mal in Malaria_level){
  for (i in 1:length(GBD_region_level)) {
    temp <- summary_wide_Fever %>%
      filter(Malaria == mal, GBD_region == GBD_region_level[i]) %>%
      dplyr::select(Year, GBD_region, Proportion_Fever_Median, `Proportion_Fever_2.5%`, `Proportion_Fever_97.5%`) %>%
      rename(
        Proportion_Median = Proportion_Fever_Median,
        Proportion_Lower = `Proportion_Fever_2.5%`,
        Proportion_Upper = `Proportion_Fever_97.5%`
      )
    obj_name <- paste0("Fever_year_proportion_GBD_region_", mal, "_", i)
    assign(obj_name, temp)
  }
}



for (mal in Malaria_level) {
  for (i in 1:7) {
    obj_Fever <- paste0("Fever_year_proportion_GBD_region_", mal,"_", i)
    path_Fever <- paste0("outputs/plot_df/SA_1.5_plot_Ab_Fever_year_proportion_GBD_region_", mal,"_", i,".csv")
    
    write.csv(get(obj_Fever), path_Fever, row.names = FALSE)
  }
}  



# #### 1.5 CSV file ####
GBD_region_distribution <- as.matrix(with(Country_Year_GBD_region_meta, addmargins(table(Year, GBD_region), margin = c(1,2))))
write.csv(GBD_region_distribution, "outputs/table/by_year_GBD_region/SA_Country_number_by_GBD_region_Year.csv")



































#___________________________________________________-----------------------

#-----------------------2.1 by country and year (DDD)------------------------

# Create matching indices
Year_index <- rep(UNI_indexed$Year, times = ncol(AbLRI_sim))
Country_index <- rep(UNI_indexed$Country, times = ncol(AbLRI_sim))
GDP_index <- rep(UNI_indexed$GDP, times = ncol(AbLRI_sim))
Income_index <- rep(UNI_indexed$Income, times = ncol(AbLRI_sim))


# Combine into dataframe
sim_df_LRI_DDD <- data.frame(
  Country = UNI_indexed$Country,
  Year = UNI_indexed$Year,
  GDP = GDP_index,
  Income = Income_index,
  LRI_DDD = LRI_DDD_long,
  sim = rep(1:10000, each = nrow(LRI_DDD_sim))
)

# Combine into dataframe
sim_df_Fever_DDD <- data.frame(
  Country = UNI_indexed$Country,
  Year = UNI_indexed$Year,
  GDP = GDP_index,
  Income = Income_index,
  Fever_DDD = Fever_DDD_long,
  sim = rep(1:10000, each = nrow(Fever_DDD_sim))
)

# Combine into dataframe
sim_df_all_DDD <- data.frame(
  Country = UNI_indexed$Country,
  Year = UNI_indexed$Year,
  GDP = GDP_index,
  Income = Income_index,
  all_DDD = all_DDD_long,
  sim = rep(1:10000, each = nrow(all_DDD_sim))
)

LRI_DDD_matrix <- sim_df_LRI_DDD %>% pivot_wider(names_from = sim, values_from = LRI_DDD)
Fever_DDD_matrix <- sim_df_Fever_DDD %>% pivot_wider(names_from = sim, values_from = Fever_DDD)
all_DDD_matrix <- sim_df_all_DDD %>% pivot_wider(names_from = sim, values_from = all_DDD)

LRI_DDD_matrix <- LRI_DDD_matrix %>% mutate(Country = ifelse(Country %in% names(country_mapping), country_mapping[Country], Country))
Fever_DDD_matrix <- Fever_DDD_matrix %>% mutate(Country = ifelse(Country %in% names(country_mapping), country_mapping[Country], Country))
all_DDD_matrix <- all_DDD_matrix %>% mutate(Country = ifelse(Country %in% names(country_mapping), country_mapping[Country], Country))

DDD_LRI <- Country_Year_meta %>% inner_join(LRI_DDD_matrix, by = c("Country","Year"))
DDD_LRI_matrix <- DDD_LRI[,5:(n_sim+4)]
DDD_Fever <- Country_Year_meta %>% inner_join(Fever_DDD_matrix, by = c("Country","Year"))
DDD_Fever_matrix <- DDD_Fever[,5:(n_sim+4)]
DDD_all <- Country_Year_meta %>% inner_join(all_DDD_matrix, by = c("Country","Year"))
DDD_all_matrix <- DDD_all[,5:(n_sim+4)]

total_DDD_LRI_matrix              <- Incidence_LRI_matrix * DDD_LRI_matrix
total_DDD_Fever_matrix            <- Incidence_Fever_matrix * DDD_Fever_matrix
total_DDD_all_matrix              <- Incidence_all_matrix * DDD_all_matrix
treated_LRI_matrix                <- Incidence_LRI_matrix*Antibiotics_LRI_matrix
treated_Fever_matrix              <- Incidence_Fever_matrix*Antibiotics_Fever_matrix
treated_all_matrix                <- Incidence_LRI_matrix*Antibiotics_LRI_matrix + Incidence_Fever_matrix*Antibiotics_Fever_matrix

DDD_per_treated_LRI               <- DDD_LRI_matrix/Antibiotics_LRI_matrix
DDD_per_treated_Fever             <- DDD_Fever_matrix/Antibiotics_Fever_matrix
DDD_per_treated_all               <- total_DDD_all_matrix /treated_all_matrix 



make_summary_2 <- function(meta, income, total, per_treated) {
  
  data.frame(
    Year = meta$Year,
    Country = meta$Country,
    Income = income$Income,
    
    Total = sprintf("%s (%s  %s)",
                    formatC(apply(total, 1, median, na.rm = TRUE), format = "f", big.mark = ",", digits = 0),
                    formatC(apply(total, 1, quantile, probs = 0.025, na.rm = TRUE), format = "f", big.mark = ",", digits = 0),
                    formatC(apply(total, 1, quantile, probs = 0.975, na.rm = TRUE), format = "f", big.mark = ",", digits = 0)),
    
    Per_treated = sprintf("%s (%s  %s)",
                          formatC(apply(per_treated, 1, median, na.rm = TRUE), format = "f", big.mark = ",", digits = 2),
                          formatC(apply(per_treated, 1, quantile, probs = 0.025, na.rm = TRUE), format = "f", big.mark = ",", digits = 2),
                          formatC(apply(per_treated, 1, quantile, probs = 0.975, na.rm = TRUE), format = "f", big.mark = ",", digits = 2))
  )
}

make_df_2 <- function(meta, income, total, per_treated) {
  data.frame(
    Year = meta$Year,
    Country = meta$Country,
    Income = income$Income,
    
    Total_median = apply(total, 1, median, na.rm = TRUE),
    Total_lower  = apply(total, 1, quantile, probs = 0.025, na.rm = TRUE),
    Total_upper  = apply(total, 1, quantile, probs = 0.975, na.rm = TRUE),
    
    Per_treated_median = apply(per_treated, 1, median, na.rm = TRUE),
    Per_treated_lower  = apply(per_treated, 1, quantile, probs = 0.025, na.rm = TRUE),
    Per_treated_upper  = apply(per_treated, 1, quantile, probs = 0.975, na.rm = TRUE)
  )
}


LRI_DDD_summary          <- make_summary_2(Country_Year_meta, DDD_LRI, total_DDD_LRI_matrix              , DDD_per_treated_LRI         )
Fever_DDD_summary        <- make_summary_2(Country_Year_meta, DDD_LRI, total_DDD_Fever_matrix            , DDD_per_treated_Fever       )
all_DDD_summary          <- make_summary_2(Country_Year_meta, DDD_LRI, total_DDD_all_matrix              , DDD_per_treated_all)

LRI_DDD_df          <- make_df_2(Country_Year_meta, DDD_LRI, total_DDD_LRI_matrix              , DDD_per_treated_LRI         )
Fever_DDD_df        <- make_df_2(Country_Year_meta, DDD_LRI, total_DDD_Fever_matrix            , DDD_per_treated_Fever       )
all_DDD_df          <- make_df_2(Country_Year_meta, DDD_LRI, total_DDD_all_matrix              , DDD_per_treated_all)

years_to_save <- 2003:2019

DDD_summary <- data.frame(
  'Year' = LRI_DDD_summary$Year,
  'Country' = LRI_DDD_summary$Country,
  'Income' = LRI_DDD_summary$Income,
  'DDD_per_treated_LRI' = LRI_DDD_summary$Per_treated,
  'DDD_per_treated_Fever_without_LRI' = Fever_DDD_summary$Per_treated,
  'DDD_per_treated_LRI_or_Fever' = all_DDD_summary$Per_treated 
)

#### 2.1 CSV file ####
for (yr in years_to_save) {
  write.csv(DDD_summary[DDD_summary$Year == yr, ], paste0("outputs/table/by_country_year/SA_DDD_per_treated_", yr, ".csv"), row.names = FALSE)
}

#### 2.1 CSV file ####
for (yr in years_to_save) {
  write.csv(LRI_DDD_df[LRI_DDD_df$Year == yr, ], paste0("outputs/global_map_df/SA_2.1_DDD_LRI_", yr, ".csv"), row.names = FALSE)
  write.csv(Fever_DDD_df[Fever_DDD_df$Year == yr, ], paste0("outputs/global_map_df/SA_2.1_DDD_Fever_", yr, ".csv"), row.names = FALSE)
  write.csv(all_DDD_df[all_DDD_df$Year == yr, ], paste0("outputs/global_map_df/SA_2.1_DDD_all_", yr, ".csv"), row.names = FALSE)
}



























#----------------------------2.2 by year (DDD) ------------------------------

years <- 2003:2019
summary_list_DDD_all <- summary_list_DDD_Fever <- summary_list_DDD_LRI <- list()

# ----- LRI------
for (yr in years) {
  # Filter matrices by year
  total_DDD_LRI_bind <- cbind(Country_Year_meta, total_DDD_LRI_matrix) %>%
    filter(Year == yr) %>%
    dplyr::select(-Country, -Year)
  
  treated_LRI_bind <- cbind(Country_Year_meta, treated_LRI_matrix) %>%
    filter(Year == yr) %>%
    dplyr::select(-Country, -Year)
  
  # Sum over countries for each simulation
  total_DDD_LRI <- colSums(total_DDD_LRI_bind)
  DDD_per_treated_LRI <- colSums(total_DDD_LRI_bind)/colSums(treated_LRI_bind)
  
  # Store results as long-format table
  summary_list_DDD_LRI[[as.character(yr)]] <- data.frame(
    Year = rep(yr, 6),
    Condition = rep(c("Total_DDD_LRI", "DDD_per_treated_LRI"), each = 3),
    Quantile = rep(c("Median", "2.5%", "97.5%"), times = 2),
    Value = c(
      quantile(total_DDD_LRI, probs = 0.50, na.rm = TRUE),
      quantile(total_DDD_LRI, probs = 0.025, na.rm = TRUE),
      quantile(total_DDD_LRI, probs = 0.975, na.rm = TRUE),
      
      quantile(DDD_per_treated_LRI, probs = 0.50, na.rm = TRUE),
      quantile(DDD_per_treated_LRI, probs = 0.025, na.rm = TRUE),
      quantile(DDD_per_treated_LRI, probs = 0.975, na.rm = TRUE)
    )
  )
}

# Combine all sublists
summary_all_years_DDD_LRI <- do.call(rbind, summary_list_DDD_LRI)

# Combine all years into one data frame
summary_wide_DDD_LRI <- summary_all_years_DDD_LRI %>%
  mutate(Column = paste0(gsub(" ", "_", Condition), "_", Quantile)) %>%
  dplyr::select(Year, Column, Value) %>%
  pivot_wider(names_from = Column, values_from = Value)

summary_wide_DDD_LRI_combined <- summary_wide_DDD_LRI %>%
  mutate(
    Total_DDD_LRI = sprintf("%s (%s  %s)", 
                            formatC(Total_DDD_LRI_Median, format = "f", digits = 0, big.mark = ","),
                            formatC(`Total_DDD_LRI_2.5%`, format = "f", digits = 0, big.mark = ","),
                            formatC(`Total_DDD_LRI_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    DDD_per_treated_LRI = sprintf("%s (%s  %s)", 
                                  formatC(DDD_per_treated_LRI_Median, format = "f", digits = 2, big.mark = ","),
                                  formatC(`DDD_per_treated_LRI_2.5%`, format = "f", digits = 2, big.mark = ","),
                                  formatC(`DDD_per_treated_LRI_97.5%`, format = "f", digits = 2, big.mark = ","))
  ) %>%
  dplyr::select(Year, Total_DDD_LRI, DDD_per_treated_LRI)

LRI_year_DDD_per_treated <- summary_wide_DDD_LRI %>%
  dplyr::select(Year, DDD_per_treated_LRI_Median, `DDD_per_treated_LRI_2.5%`, `DDD_per_treated_LRI_97.5%`) %>%
  dplyr::rename(DDD_per_treated_LRI = DDD_per_treated_LRI_Median, DDD_per_treated_LRI_lower = `DDD_per_treated_LRI_2.5%`, DDD_per_treated_LRI_upper = `DDD_per_treated_LRI_97.5%`)






# ----- Fever------
for (yr in years) {
  # Filter matrices by year
  total_DDD_Fever_bind <- cbind(Country_Year_meta, total_DDD_Fever_matrix) %>%
    filter(Year == yr) %>%
    dplyr::select(-Country, -Year)
  
  treated_Fever_bind <- cbind(Country_Year_meta, treated_Fever_matrix) %>%
    filter(Year == yr) %>%
    dplyr::select(-Country, -Year)
  
  # Sum over countries for each simulation
  total_DDD_Fever <- colSums(total_DDD_Fever_bind)
  DDD_per_treated_Fever <- colSums(total_DDD_Fever_bind)/colSums(treated_Fever_bind)
  
  # Store results as long-format table
  summary_list_DDD_Fever[[as.character(yr)]] <- data.frame(
    Year = rep(yr, 6),
    Condition = rep(c("Total_DDD_Fever", "DDD_per_treated_Fever"), each = 3),
    Quantile = rep(c("Median", "2.5%", "97.5%"), times = 2),
    Value = c(
      quantile(total_DDD_Fever, probs = 0.50, na.rm = TRUE),
      quantile(total_DDD_Fever, probs = 0.025, na.rm = TRUE),
      quantile(total_DDD_Fever, probs = 0.975, na.rm = TRUE),
      
      quantile(DDD_per_treated_Fever, probs = 0.50, na.rm = TRUE),
      quantile(DDD_per_treated_Fever, probs = 0.025, na.rm = TRUE),
      quantile(DDD_per_treated_Fever, probs = 0.975, na.rm = TRUE)
    )
  )
}

# Combine all sublists
summary_all_years_DDD_Fever <- do.call(rbind, summary_list_DDD_Fever)

# Combine all years into one data frame
summary_wide_DDD_Fever <- summary_all_years_DDD_Fever %>%
  mutate(Column = paste0(gsub(" ", "_", Condition), "_", Quantile)) %>%
  dplyr::select(Year, Column, Value) %>%
  pivot_wider(names_from = Column, values_from = Value)

summary_wide_DDD_Fever_combined <- summary_wide_DDD_Fever %>%
  mutate(
    Total_DDD_Fever = sprintf("%s (%s  %s)", 
                              formatC(Total_DDD_Fever_Median, format = "f", digits = 0, big.mark = ","),
                              formatC(`Total_DDD_Fever_2.5%`, format = "f", digits = 0, big.mark = ","),
                              formatC(`Total_DDD_Fever_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    DDD_per_treated_Fever = sprintf("%s (%s  %s)", 
                                    formatC(DDD_per_treated_Fever_Median, format = "f", digits = 2, big.mark = ","),
                                    formatC(`DDD_per_treated_Fever_2.5%`, format = "f", digits = 2, big.mark = ","),
                                    formatC(`DDD_per_treated_Fever_97.5%`, format = "f", digits = 2, big.mark = ","))
  ) %>%
  dplyr::select(Year, Total_DDD_Fever, DDD_per_treated_Fever)

Fever_year_DDD_per_treated <- summary_wide_DDD_Fever %>%
  dplyr::select(Year, DDD_per_treated_Fever_Median, `DDD_per_treated_Fever_2.5%`, `DDD_per_treated_Fever_97.5%`) %>%
  dplyr::rename(DDD_per_treated_Fever = DDD_per_treated_Fever_Median, DDD_per_treated_Fever_lower = `DDD_per_treated_Fever_2.5%`, DDD_per_treated_Fever_upper = `DDD_per_treated_Fever_97.5%`)









# ----- all------
for (yr in years) {
  # Filter matrices by year
  total_DDD_all_bind <- cbind(Country_Year_meta, total_DDD_all_matrix) %>%
    filter(Year == yr) %>%
    dplyr::select(-Country, -Year)
  
  treated_all_bind <- cbind(Country_Year_meta, treated_all_matrix) %>%
    filter(Year == yr) %>%
    dplyr::select(-Country, -Year)
  
  # Sum over countries for each simulation
  total_DDD_all <- colSums(total_DDD_all_bind)
  DDD_per_treated_all <- colSums(total_DDD_all_bind) / colSums(treated_all_bind)
  
  # Store results as long-format table
  summary_list_DDD_all[[as.character(yr)]] <- data.frame(
    Year = rep(yr, 6),
    Condition = rep(c("Total_DDD_all", "DDD_per_treated_all"), each = 3),
    Quantile = rep(c("Median", "2.5%", "97.5%"), times = 2),
    Value = c(
      quantile(total_DDD_all, probs = 0.50, na.rm = TRUE),
      quantile(total_DDD_all, probs = 0.025, na.rm = TRUE),
      quantile(total_DDD_all, probs = 0.975, na.rm = TRUE),
      
      quantile(DDD_per_treated_all, probs = 0.50, na.rm = TRUE),
      quantile(DDD_per_treated_all, probs = 0.025, na.rm = TRUE),
      quantile(DDD_per_treated_all, probs = 0.975, na.rm = TRUE)
    )
  )
}

# Combine all sublists
summary_all_years_DDD_all <- do.call(rbind, summary_list_DDD_all)

# Combine all years into one data frame
summary_wide_DDD_all <- summary_all_years_DDD_all %>%
  mutate(Column = paste0(gsub(" ", "_", Condition), "_", Quantile)) %>%
  dplyr::select(Year, Column, Value) %>%
  pivot_wider(names_from = Column, values_from = Value)

summary_wide_DDD_all_combined <- summary_wide_DDD_all %>%
  mutate(
    Total_DDD_all = sprintf("%s (%s  %s)", 
                            formatC(Total_DDD_all_Median, format = "f", digits = 0, big.mark = ","),
                            formatC(`Total_DDD_all_2.5%`, format = "f", digits = 0, big.mark = ","),
                            formatC(`Total_DDD_all_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    DDD_per_treated_all = sprintf("%s (%s  %s)", 
                                  formatC(DDD_per_treated_all_Median, format = "f", digits = 2, big.mark = ","),
                                  formatC(`DDD_per_treated_all_2.5%`, format = "f", digits = 2, big.mark = ","),
                                  formatC(`DDD_per_treated_all_97.5%`, format = "f", digits = 2, big.mark = ","))
  ) %>%
  dplyr::select(Year, Total_DDD_all, DDD_per_treated_all)

all_year_DDD_per_treated <- summary_wide_DDD_all %>%
  dplyr::select(Year, DDD_per_treated_all_Median, `DDD_per_treated_all_2.5%`, `DDD_per_treated_all_97.5%`) %>%
  dplyr::rename(DDD_per_treated_all = DDD_per_treated_all_Median, DDD_per_treated_all_lower = `DDD_per_treated_all_2.5%`, DDD_per_treated_all_upper = `DDD_per_treated_all_97.5%`)

DDD_by_year <- data.frame(
  'Year' = summary_wide_DDD_LRI_combined$Year,
  'DDD_per_treated_LRI' = summary_wide_DDD_LRI_combined$DDD_per_treated_LRI, 
  'DDD_per_treated_Fever_without_LRI' = summary_wide_DDD_Fever_combined$DDD_per_treated_Fever,
  'DDD_per_treated_LRI_or_Fever' = summary_wide_DDD_all_combined$DDD_per_treated_all 
)

#### 2.2 CSV file ####
write.csv(DDD_by_year, "outputs/table/by_year/SA_DDD_per_treated_by_year.csv", row.names = FALSE)







#----------------------------2.2 by year (DDD increase) ------------------------------



total_DDD_LRI_2003 <- cbind(Country_Year_meta, total_DDD_LRI_matrix) %>%
  filter(Year == 2003) %>%
  dplyr::select(-Country, -Year)

incidence_LRI_2003 <- cbind(Country_Year_meta, Incidence_LRI_matrix) %>%
  filter(Year == 2003) %>%
  dplyr::select(-Country, -Year)

DDD_per_LRI_2003 <- colSums(total_DDD_LRI_2003)/colSums(incidence_LRI_2003)

total_DDD_LRI_2019 <- cbind(Country_Year_meta, total_DDD_LRI_matrix) %>%
  filter(Year == 2019) %>%
  dplyr::select(-Country, -Year)

incidence_LRI_2019 <- cbind(Country_Year_meta, Incidence_LRI_matrix) %>%
  filter(Year == 2019) %>%
  dplyr::select(-Country, -Year)

DDD_per_LRI_2019 <- colSums(total_DDD_LRI_2019)/colSums(incidence_LRI_2019)

ratio_LRI <- DDD_per_LRI_2019/DDD_per_LRI_2003

round(quantile(ratio_LRI, probs = 0.50, na.rm = TRUE), digits = 3)
round(quantile(ratio_LRI, probs = 0.025, na.rm = TRUE), digits = 3)
round(quantile(ratio_LRI, probs = 0.975, na.rm = TRUE), digits = 3)




total_DDD_Fever_2003 <- cbind(Country_Year_meta, total_DDD_Fever_matrix) %>%
  filter(Year == 2003) %>%
  dplyr::select(-Country, -Year)

incidence_Fever_2003 <- cbind(Country_Year_meta, Incidence_Fever_matrix) %>%
  filter(Year == 2003) %>%
  dplyr::select(-Country, -Year)

DDD_per_Fever_2003 <- colSums(total_DDD_Fever_2003)/colSums(incidence_Fever_2003)

total_DDD_Fever_2019 <- cbind(Country_Year_meta, total_DDD_Fever_matrix) %>%
  filter(Year == 2019) %>%
  dplyr::select(-Country, -Year)

incidence_Fever_2019 <- cbind(Country_Year_meta, Incidence_Fever_matrix) %>%
  filter(Year == 2019) %>%
  dplyr::select(-Country, -Year)

DDD_per_Fever_2019 <- colSums(total_DDD_Fever_2019)/colSums(incidence_Fever_2019)

ratio_Fever <- DDD_per_Fever_2019/DDD_per_Fever_2003

round(quantile(ratio_Fever, probs = 0.50, na.rm = TRUE), digits = 3)
round(quantile(ratio_Fever, probs = 0.025, na.rm = TRUE), digits = 3)
round(quantile(ratio_Fever, probs = 0.975, na.rm = TRUE), digits = 3)




























#----------------------------2.3 by Income and year (treatment) ----------------------------

# -------LRI---------
Country_Year_Income_meta <- Country_Year_meta %>% inner_join(LRI_DDD_matrix[,c(1,2,4)], by = c("Country", "Year"))
Income_level <- c('L','LM','UM')
years <- 2003:2019
summary_list_DDD_LRI <- list()



for (inc in Income_level) {
  for (yr in years) {
    # Filter and select simulation matrices
    total_DDD_LRI_bind <- cbind(Country_Year_Income_meta, total_DDD_LRI_matrix) %>%
      filter(Income == inc, Year == yr) %>%
      dplyr::select(-Country, -Year, -Income)
    
    treated_LRI_bind <- cbind(Country_Year_Income_meta, treated_LRI_matrix) %>%
      filter(Income == inc, Year == yr) %>%
      dplyr::select(-Country, -Year, -Income)
    
    # Sum over countries for each simulation
    total_DDD_LRI <- colSums(total_DDD_LRI_bind)
    DDD_per_treated_LRI <- colSums(total_DDD_LRI_bind)/colSums(treated_LRI_bind)
    
    # Store results as long-format table
    summary_list_DDD_LRI[[paste(inc, yr, sep = "_")]] <- data.frame(
      Income = inc,
      Year = rep(yr, 6),
      Condition = rep(c("Total_DDD_LRI", "DDD_per_treated_LRI"), each = 3),
      Quantile = rep(c("Median", "2.5%", "97.5%"), times = 2),
      Value = c(
        quantile(total_DDD_LRI, probs = 0.50, na.rm = TRUE),
        quantile(total_DDD_LRI, probs = 0.025, na.rm = TRUE),
        quantile(total_DDD_LRI, probs = 0.975, na.rm = TRUE),
        
        quantile(DDD_per_treated_LRI, probs = 0.50, na.rm = TRUE),
        quantile(DDD_per_treated_LRI, probs = 0.025, na.rm = TRUE),
        quantile(DDD_per_treated_LRI, probs = 0.975, na.rm = TRUE)
      )
    )
  }
}

# Combine all sublists
summary_all_years_DDD_LRI <- do.call(rbind, summary_list_DDD_LRI)

# Combine all years into one data frame
summary_wide_DDD_LRI <- summary_all_years_DDD_LRI %>%
  mutate(Column = paste0(gsub(" ", "_", Condition), "_", Quantile)) %>%
  dplyr::select(Year, Income, Column, Value) %>%
  pivot_wider(names_from = Column, values_from = Value)


summary_wide_DDD_LRI_combined <- summary_wide_DDD_LRI %>%
  mutate(
    Total_DDD_LRI = sprintf("%s (%s  %s)", 
                            formatC(Total_DDD_LRI_Median, format = "f", digits = 0, big.mark = ","),
                            formatC(`Total_DDD_LRI_2.5%`, format = "f", digits = 0, big.mark = ","),
                            formatC(`Total_DDD_LRI_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    DDD_per_treated_LRI = sprintf("%s (%s  %s)", 
                                  formatC(DDD_per_treated_LRI_Median, format = "f", digits = 2, big.mark = ","),
                                  formatC(`DDD_per_treated_LRI_2.5%`, format = "f", digits = 2, big.mark = ","),
                                  formatC(`DDD_per_treated_LRI_97.5%`, format = "f", digits = 2, big.mark = ","))
  ) %>%
  dplyr::select(Year, Income, Total_DDD_LRI, DDD_per_treated_LRI)



LRI_year_DDD_per_treated_L <- summary_wide_DDD_LRI %>%
  filter(Income == 'L') %>%
  dplyr::select(Year, Income, DDD_per_treated_LRI_Median, `DDD_per_treated_LRI_2.5%`, `DDD_per_treated_LRI_97.5%`) %>%
  dplyr::rename(DDD_per_treated_Median = DDD_per_treated_LRI_Median, DDD_per_treated_LRI_Lower = `DDD_per_treated_LRI_2.5%`, DDD_per_treated_LRI_Upper = `DDD_per_treated_LRI_97.5%`)

LRI_year_DDD_per_treated_LM <- summary_wide_DDD_LRI %>%
  filter(Income == 'LM') %>%
  dplyr::select(Year, Income, DDD_per_treated_LRI_Median, `DDD_per_treated_LRI_2.5%`, `DDD_per_treated_LRI_97.5%`) %>%
  dplyr::rename(DDD_per_treated_Median = DDD_per_treated_LRI_Median, DDD_per_treated_LRI_Lower = `DDD_per_treated_LRI_2.5%`, DDD_per_treated_LRI_Upper = `DDD_per_treated_LRI_97.5%`)

LRI_year_DDD_per_treated_UM <- summary_wide_DDD_LRI %>%
  filter(Income == 'UM') %>%
  dplyr::select(Year, Income, DDD_per_treated_LRI_Median, `DDD_per_treated_LRI_2.5%`, `DDD_per_treated_LRI_97.5%`) %>%
  dplyr::rename(DDD_per_treated_Median = DDD_per_treated_LRI_Median, DDD_per_treated_LRI_Lower = `DDD_per_treated_LRI_2.5%`, DDD_per_treated_LRI_Upper = `DDD_per_treated_LRI_97.5%`)









# -------Fever---------
Country_Year_Income_meta <- Country_Year_meta %>% inner_join(Fever_DDD_matrix[,c(1,2,4)], by = c("Country", "Year"))
Income_level <- c('L','LM','UM')
years <- 2003:2019
summary_list_DDD_Fever <- list()



for (inc in Income_level) {
  for (yr in years) {
    # Filter and select simulation matrices
    total_DDD_Fever_bind <- cbind(Country_Year_Income_meta, total_DDD_Fever_matrix) %>%
      filter(Income == inc, Year == yr) %>%
      dplyr::select(-Country, -Year, -Income)
    
    treated_Fever_bind <- cbind(Country_Year_Income_meta, treated_Fever_matrix) %>%
      filter(Income == inc, Year == yr) %>%
      dplyr::select(-Country, -Year, -Income)
    
    # Sum over countries for each simulation
    total_DDD_Fever <- colSums(total_DDD_Fever_bind)
    DDD_per_treated_Fever <- colSums(total_DDD_Fever_bind)/colSums(treated_Fever_bind)
    
    # Store results as long-format table
    summary_list_DDD_Fever[[paste(inc, yr, sep = "_")]] <- data.frame(
      Income = inc,
      Year = rep(yr, 6),
      Condition = rep(c("Total_DDD_Fever", "DDD_per_treated_Fever"), each = 3),
      Quantile = rep(c("Median", "2.5%", "97.5%"), times = 2),
      Value = c(
        quantile(total_DDD_Fever, probs = 0.50, na.rm = TRUE),
        quantile(total_DDD_Fever, probs = 0.025, na.rm = TRUE),
        quantile(total_DDD_Fever, probs = 0.975, na.rm = TRUE),
        
        quantile(DDD_per_treated_Fever, probs = 0.50, na.rm = TRUE),
        quantile(DDD_per_treated_Fever, probs = 0.025, na.rm = TRUE),
        quantile(DDD_per_treated_Fever, probs = 0.975, na.rm = TRUE)
      )
    )
  }
}

# Combine all sublists
summary_all_years_DDD_Fever <- do.call(rbind, summary_list_DDD_Fever)

# Combine all years into one data frame
summary_wide_DDD_Fever <- summary_all_years_DDD_Fever %>%
  mutate(Column = paste0(gsub(" ", "_", Condition), "_", Quantile)) %>%
  dplyr::select(Year, Income, Column, Value) %>%
  pivot_wider(names_from = Column, values_from = Value)


summary_wide_DDD_Fever_combined <- summary_wide_DDD_Fever %>%
  mutate(
    Total_DDD_Fever = sprintf("%s (%s  %s)", 
                              formatC(Total_DDD_Fever_Median, format = "f", digits = 0, big.mark = ","),
                              formatC(`Total_DDD_Fever_2.5%`, format = "f", digits = 0, big.mark = ","),
                              formatC(`Total_DDD_Fever_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    DDD_per_treated_Fever = sprintf("%s (%s  %s)", 
                                    formatC(DDD_per_treated_Fever_Median, format = "f", digits = 2, big.mark = ","),
                                    formatC(`DDD_per_treated_Fever_2.5%`, format = "f", digits = 2, big.mark = ","),
                                    formatC(`DDD_per_treated_Fever_97.5%`, format = "f", digits = 2, big.mark = ","))
  ) %>%
  dplyr::select(Year, Income, Total_DDD_Fever, DDD_per_treated_Fever)



Fever_year_DDD_per_treated_L <- summary_wide_DDD_Fever %>%
  filter(Income == 'L') %>%
  dplyr::select(Year, Income, DDD_per_treated_Fever_Median, `DDD_per_treated_Fever_2.5%`, `DDD_per_treated_Fever_97.5%`) %>%
  dplyr::rename(DDD_per_treated_Median = DDD_per_treated_Fever_Median, DDD_per_treated_Fever_Lower = `DDD_per_treated_Fever_2.5%`, DDD_per_treated_Fever_Upper = `DDD_per_treated_Fever_97.5%`)

Fever_year_DDD_per_treated_LM <- summary_wide_DDD_Fever %>%
  filter(Income == 'LM') %>%
  dplyr::select(Year, Income, DDD_per_treated_Fever_Median, `DDD_per_treated_Fever_2.5%`, `DDD_per_treated_Fever_97.5%`) %>%
  dplyr::rename(DDD_per_treated_Median = DDD_per_treated_Fever_Median, DDD_per_treated_Fever_Lower = `DDD_per_treated_Fever_2.5%`, DDD_per_treated_Fever_Upper = `DDD_per_treated_Fever_97.5%`)

Fever_year_DDD_per_treated_UM <- summary_wide_DDD_Fever %>%
  filter(Income == 'UM') %>%
  dplyr::select(Year, Income, DDD_per_treated_Fever_Median, `DDD_per_treated_Fever_2.5%`, `DDD_per_treated_Fever_97.5%`) %>%
  dplyr::rename(DDD_per_treated_Median = DDD_per_treated_Fever_Median, DDD_per_treated_Fever_Lower = `DDD_per_treated_Fever_2.5%`, DDD_per_treated_Fever_Upper = `DDD_per_treated_Fever_97.5%`)









# -------all---------
Country_Year_Income_meta <- Country_Year_meta %>% inner_join(all_DDD_matrix[,c(1,2,4)], by = c("Country", "Year"))
Income_level <- c('L','LM','UM')
years <- 2003:2019
summary_list_DDD_all <- list()



for (inc in Income_level) {
  for (yr in years) {
    # Filter and select simulation matrices
    total_DDD_all_bind <- cbind(Country_Year_Income_meta, total_DDD_all_matrix) %>%
      filter(Income == inc, Year == yr) %>%
      dplyr::select(-Country, -Year, -Income)
    
    treated_all_bind <- cbind(Country_Year_Income_meta, treated_all_matrix) %>%
      filter(Income == inc, Year == yr) %>%
      dplyr::select(-Country, -Year, -Income)
    
    # Sum over countries for each simulation
    total_DDD_all <- colSums(total_DDD_all_bind)
    DDD_per_treated_all <- colSums(total_DDD_all_bind)/colSums(treated_all_bind)
    
    # Store results as long-format table
    summary_list_DDD_all[[paste(inc, yr, sep = "_")]] <- data.frame(
      Income = inc,
      Year = rep(yr, 6),
      Condition = rep(c("Total_DDD_all", "DDD_per_treated_all"), each = 3),
      Quantile = rep(c("Median", "2.5%", "97.5%"), times = 2),
      Value = c(
        quantile(total_DDD_all, probs = 0.50, na.rm = TRUE),
        quantile(total_DDD_all, probs = 0.025, na.rm = TRUE),
        quantile(total_DDD_all, probs = 0.975, na.rm = TRUE),
        
        quantile(DDD_per_treated_all, probs = 0.50, na.rm = TRUE),
        quantile(DDD_per_treated_all, probs = 0.025, na.rm = TRUE),
        quantile(DDD_per_treated_all, probs = 0.975, na.rm = TRUE)
      )
    )
  }
}

# Combine all sublists
summary_all_years_DDD_all <- do.call(rbind, summary_list_DDD_all)

# Combine all years into one data frame
summary_wide_DDD_all <- summary_all_years_DDD_all %>%
  mutate(Column = paste0(gsub(" ", "_", Condition), "_", Quantile)) %>%
  dplyr::select(Year, Income, Column, Value) %>%
  pivot_wider(names_from = Column, values_from = Value)


summary_wide_DDD_all_combined <- summary_wide_DDD_all %>%
  mutate(
    Total_DDD_all = sprintf("%s (%s  %s)", 
                            formatC(Total_DDD_all_Median, format = "f", digits = 0, big.mark = ","),
                            formatC(`Total_DDD_all_2.5%`, format = "f", digits = 0, big.mark = ","),
                            formatC(`Total_DDD_all_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    DDD_per_treated_all = sprintf("%s (%s  %s)", 
                                  formatC(DDD_per_treated_all_Median, format = "f", digits = 2, big.mark = ","),
                                  formatC(`DDD_per_treated_all_2.5%`, format = "f", digits = 2, big.mark = ","),
                                  formatC(`DDD_per_treated_all_97.5%`, format = "f", digits = 2, big.mark = ","))
  ) %>%
  dplyr::select(Year, Income, Total_DDD_all, DDD_per_treated_all)



all_year_DDD_per_treated_L <- summary_wide_DDD_all %>%
  filter(Income == 'L') %>%
  dplyr::select(Year, Income, DDD_per_treated_all_Median, `DDD_per_treated_all_2.5%`, `DDD_per_treated_all_97.5%`) %>%
  dplyr::rename(DDD_per_treated_Median = DDD_per_treated_all_Median, DDD_per_treated_all_Lower = `DDD_per_treated_all_2.5%`, DDD_per_treated_all_Upper = `DDD_per_treated_all_97.5%`)

all_year_DDD_per_treated_LM <- summary_wide_DDD_all %>%
  filter(Income == 'LM') %>%
  dplyr::select(Year, Income, DDD_per_treated_all_Median, `DDD_per_treated_all_2.5%`, `DDD_per_treated_all_97.5%`) %>%
  dplyr::rename(DDD_per_treated_Median = DDD_per_treated_all_Median, DDD_per_treated_all_Lower = `DDD_per_treated_all_2.5%`, DDD_per_treated_all_Upper = `DDD_per_treated_all_97.5%`)

all_year_DDD_per_treated_UM <- summary_wide_DDD_all %>%
  filter(Income == 'UM') %>%
  dplyr::select(Year, Income, DDD_per_treated_all_Median, `DDD_per_treated_all_2.5%`, `DDD_per_treated_all_97.5%`) %>%
  dplyr::rename(DDD_per_treated_Median = DDD_per_treated_all_Median, DDD_per_treated_all_Lower = `DDD_per_treated_all_2.5%`, DDD_per_treated_all_Upper = `DDD_per_treated_all_97.5%`)



DDD_by_year_income <- data.frame(
  'Year' = summary_wide_DDD_LRI_combined$Year,
  'Income' = summary_wide_DDD_LRI_combined$Income,
  'DDD_per_treated_LRI' = summary_wide_DDD_LRI_combined$DDD_per_treated_LRI, 
  'DDD_per_treated_Fever_without_LRI' = summary_wide_DDD_Fever_combined$DDD_per_treated_Fever,
  'DDD_per_treated_LRI_or_Fever' = summary_wide_DDD_all_combined$DDD_per_treated_all 
)

#### 2.2 CSV file ####
write.csv(DDD_by_year_income, "outputs/table/by_year_income/SA_DDD_per_treated_by_year_income.csv", row.names = FALSE)


# save the plot data
write.csv(LRI_year_DDD_per_treated     , "outputs/plot_df/SA_2.3_plot_Ab_LRI_year_DDD_per_treated.csv")
write.csv(LRI_year_DDD_per_treated_L   , "outputs/plot_df/SA_2.3_plot_Ab_LRI_year_DDD_per_treated_L.csv")
write.csv(LRI_year_DDD_per_treated_LM  , "outputs/plot_df/SA_2.3_plot_Ab_LRI_year_DDD_per_treated_LM.csv")
write.csv(LRI_year_DDD_per_treated_UM  , "outputs/plot_df/SA_2.3_plot_Ab_LRI_year_DDD_per_treated_UM.csv")
write.csv(Fever_year_DDD_per_treated   , "outputs/plot_df/SA_2.3_plot_Ab_Fever_year_DDD_per_treated.csv")
write.csv(Fever_year_DDD_per_treated_L , "outputs/plot_df/SA_2.3_plot_Ab_Fever_year_DDD_per_treated_L.csv")
write.csv(Fever_year_DDD_per_treated_LM, "outputs/plot_df/SA_2.3_plot_Ab_Fever_year_DDD_per_treated_LM.csv")
write.csv(Fever_year_DDD_per_treated_UM, "outputs/plot_df/SA_2.3_plot_Ab_Fever_year_DDD_per_treated_UM.csv")
write.csv(all_year_DDD_per_treated     , "outputs/plot_df/SA_2.3_plot_Ab_all_year_DDD_per_treated.csv")
write.csv(all_year_DDD_per_treated_L   , "outputs/plot_df/SA_2.3_plot_Ab_all_year_DDD_per_treated_L.csv")
write.csv(all_year_DDD_per_treated_LM  , "outputs/plot_df/SA_2.3_plot_Ab_all_year_DDD_per_treated_LM.csv")
write.csv(all_year_DDD_per_treated_UM  , "outputs/plot_df/SA_2.3_plot_Ab_all_year_DDD_per_treated_UM.csv")

#


######## 8. Figure.1 ########

#### 0. prep ####
rm(list = ls())
gc()

requiredPackages = c("dplyr", "tidyr", "rnaturalearth", "countrycode", "ggplot2", "patchwork","cowplot")
lib.path <- .libPaths()
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p, lib = lib.path)
  library(p,character.only = TRUE)
}

options(max.print = 1e6)

n_sim <- 10000
set.seed(123)


# confirm upstream extrapolation models use PC1-PC2
coef_df_LRI <- readRDS("outputs/sensitivity_analysis/SA_coef_df_LRI.rds")
coef_df_Fever <- readRDS("outputs/sensitivity_analysis/SA_coef_df_Fever.rds")
coef_df_LRI_DDD <- readRDS("outputs/sensitivity_analysis/SA_coef_df_LRI_DDD.rds")
coef_df_Fever_DDD <- readRDS("outputs/sensitivity_analysis/SA_coef_df_Fever_DDD.rds")

expected_coef_names <- c("(Intercept)", "PC1", "PC2")

stopifnot(
  identical(names(coef_df_LRI), expected_coef_names),
  identical(names(coef_df_Fever), expected_coef_names),
  identical(names(coef_df_LRI_DDD), expected_coef_names),
  identical(names(coef_df_Fever_DDD), expected_coef_names)
)

rm(
  coef_df_LRI,
  coef_df_Fever,
  coef_df_LRI_DDD,
  coef_df_Fever_DDD,
  expected_coef_names
)


# CSVs below should be regenerated from the current PC1-PC2 models

years_to_read <- c(2019, 2016, 2015, 2014, 2011, 2010, 2009, 2006, 2005, 2004)
Probability_LRI <- Probability_Fever_0 <- Probability_Fever_1 <- DDD_per_treated_LRI <- DDD_per_treated_Fever <- list()

for (yr in years_to_read) {
  Probability_LRI[[as.character(yr)]] <- read.csv(paste0("outputs/global_map_df/SA_1.1_Untreated_number_", yr, "_LRI.csv"))
  Probability_Fever_0[[as.character(yr)]] <- read.csv(paste0("outputs/global_map_df/SA_1.1_Untreated_number_0_", yr, "_Fever.csv"))
  Probability_Fever_1[[as.character(yr)]] <- read.csv(paste0("outputs/global_map_df/SA_1.1_Untreated_number_1_", yr, "_Fever.csv"))
  DDD_per_treated_LRI[[as.character(yr)]] <- read.csv(paste0("outputs/global_map_df/SA_2.1_DDD_LRI_", yr, ".csv"))
  DDD_per_treated_Fever[[as.character(yr)]] <- read.csv(paste0("outputs/global_map_df/SA_2.1_DDD_Fever_", yr, ".csv"))
}

Probability_LRI_raw <- do.call(rbind, Probability_LRI)
Probability_Fever_raw <- do.call(rbind, c(Probability_Fever_0, Probability_Fever_1))
DDD_per_treated_LRI_raw <- do.call(rbind, DDD_per_treated_LRI)
DDD_per_treated_Fever_raw <- do.call(rbind, DDD_per_treated_Fever)

Probability_LRI_sub <- Probability_LRI_raw %>% dplyr::select(Year, Country, Treated_incidence = Treated_median, Total_incidence = Total_median)
Probability_Fever_sub <- Probability_Fever_raw %>% dplyr::select(Year, Country, Treated_incidence = Treated_median, Total_incidence = Total_median)
DDD_per_treated_LRI_sub <- DDD_per_treated_LRI_raw %>% dplyr::select(Year, Country, Total_DDD = Total_median)
DDD_per_treated_Fever_sub <- DDD_per_treated_Fever_raw %>% dplyr::select(Year, Country, Total_DDD = Total_median)

LRI_joint <- Probability_LRI_sub %>% left_join(DDD_per_treated_LRI_sub, by = c("Country", "Year"))
Fever_joint <- Probability_Fever_sub %>% left_join(DDD_per_treated_Fever_sub, by = c("Country", "Year"))

LRI_joint <- LRI_joint %>%
  mutate(
    Proportion_median = Treated_incidence/Total_incidence,
    Per_treated_median = Total_DDD/Treated_incidence
  )

Fever_joint <- Fever_joint %>%
  mutate(
    Proportion_median = Treated_incidence/Total_incidence,
    Per_treated_median = Total_DDD/Treated_incidence
  )


world <- ne_countries()

a <- unique(world$admin)
b <- unique(LRI_joint$Country)
setdiff(a,b)
setdiff(b,a)

country_mapping <- c(
  "Bolivia (Plurinational State of)" = "Bolivia",
  "Congo" = "Republic of the Congo",
  "Côte d'Ivoire" = "Ivory Coast",
  "Eswatini" = "eSwatini",
  "Iran (Islamic Republic of)" = "Iran",
  "Lao People's Democratic Republic" = "Laos",
  "Republic of Moldova" = "Moldova",
  "Russian Federation" = "Russia",
  "Syrian Arab Republic" = "Syria",
  "Timor-Leste" = "East Timor",
  "Türkiye" = "Turkey",
  "Serbia" = "Republic of Serbia",
  "Viet Nam" = "Vietnam"
)

LRI_joint           <- LRI_joint       %>% mutate(Country = ifelse(Country %in% names(country_mapping), country_mapping[Country], Country))
Fever_joint         <- Fever_joint   %>% mutate(Country = ifelse(Country %in% names(country_mapping), country_mapping[Country], Country))




# ============================================================
# Export source data for Supplemental Figure SA1
# ============================================================

dir.create(
  "outputs/plot_objects",
  recursive = TRUE,
  showWarnings = FALSE
)


# Panel A: LRTI, 2019
fig_SA1_A <- LRI_joint %>%
  filter(Year == 2019) %>%
  group_by(Country) %>%
  summarise(
    Treated_incidence = sum(Treated_incidence, na.rm = TRUE),
    Total_incidence = sum(Total_incidence, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Panel = "A",
    Syndrome = "LRTI",
    Year = 2019,
    Treatment_probability = Treated_incidence / Total_incidence
  ) %>%
  dplyr::select(
    Panel,
    Syndrome,
    Country,
    Year,
    Treated_incidence,
    Total_incidence,
    Treatment_probability
  )


# Panel B: AFI, 2019
fig_SA1_B <- Fever_joint %>%
  filter(Year == 2019) %>%
  group_by(Country) %>%
  summarise(
    Treated_incidence = sum(Treated_incidence, na.rm = TRUE),
    Total_incidence = sum(Total_incidence, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Panel = "B",
    Syndrome = "AFI",
    Year = 2019,
    Treatment_probability = Treated_incidence / Total_incidence
  ) %>%
  dplyr::select(
    Panel,
    Syndrome,
    Country,
    Year,
    Treated_incidence,
    Total_incidence,
    Treatment_probability
  )


# Combine panels
Figure_SA1_source_data <- bind_rows(
  fig_SA1_A,
  fig_SA1_B
) %>%
  arrange(Panel, Country)


# Save
write.csv(
  Figure_SA1_source_data,
  "outputs/plot_objects/Figure_SA1_source_data.csv",
  row.names = FALSE
)




#----LRI-----
periods <- list(
  "2004_2006" = 2004:2006, "2009_2011" = 2009:2011,
  "2014_2016" = 2014:2016, "2019" = 2019)

plot_list_probability_LRI <- list()

for (period_name in names(periods)) {
  years <- periods[[period_name]]
  
  filtered_data <- LRI_joint %>% filter(Year %in% years) %>% group_by(Country) %>% 
    summarize(
      Proportion_median = sum(Treated_incidence)/sum(Total_incidence), 
      .groups = "drop")
  
  world_data <- world %>% left_join(filtered_data, by = c("admin" = "Country"))
  time_label <- if (period_name == "2019") {
    " in 2019"} else {paste0(" from ", gsub("_", " to ", period_name))
    }
  
  p <- ggplot(world_data) +
    geom_sf(aes(fill = Proportion_median)) +
    scale_fill_viridis_c(
      option = "viridis",
      limits = c(0.2, 1),
      breaks = c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
      labels = c("<0.3", "0.4", "0.5", "0.6", "0.7", "0.8", ">0.9"),
      oob = scales::squish,
      name = paste0("Antibiotic treatment probability for LRTI", time_label),
      na.value = "gray"
    ) +
    theme_minimal() +
    labs(
      title = paste("Antibiotic treatment probability for LRTI", time_label)
    ) +
    theme(
      legend.position = c(0.03, 0.15),        
      legend.justification = c(0, 0),         
      legend.title = element_text(size = 6), 
      legend.text = element_text(size = 6),
      legend.key.height = unit(0.5, "cm"),
      legend.key.width = unit(0.3, "cm"),
      plot.title = element_text(size = 8, face = "bold"),
      plot.subtitle = element_text(size = 8)
    )
  
  plot_list_probability_LRI[[period_name]] <- p
}

tiff("outputs/table/by_year_malaria/SA_probability_LRTI_2004_2006.tiff", width = 7.5, height = 5, units = "in", res = 300)
plot_list_probability_LRI[["2004_2006"]]
dev.off()
tiff("outputs/table/by_year_malaria/SA_probability_LRTI_2009_2011.tiff", width = 7.5, height = 5, units = "in", res = 300)
plot_list_probability_LRI[["2009_2011"]]
dev.off()
tiff("outputs/table/by_year_malaria/SA_probability_LRTI_2014_2016.tiff", width = 7.5, height = 5, units = "in", res = 300)
plot_list_probability_LRI[["2014_2016"]]
dev.off()
tiff("outputs/table/by_year_malaria/SA_probability_LRTI_2019.tiff", width = 7.5, height = 5, units = "in", res = 300)
plot_list_probability_LRI[["2019"]]
dev.off()






#----Fever-----
periods <- list(
  "2004_2006" = 2004:2006, "2009_2011" = 2009:2011,
  "2014_2016" = 2014:2016, "2019" = 2019)

plot_list_probability_Fever <- list()

for (period_name in names(periods)) {
  years <- periods[[period_name]]
  
  filtered_data <- Fever_joint %>% filter(Year %in% years) %>% group_by(Country) %>% 
    summarize(
      Proportion_median = sum(Treated_incidence)/sum(Total_incidence), 
      .groups = "drop")
  
  world_data <- world %>% left_join(filtered_data, by = c("admin" = "Country"))
  time_label <- if (period_name == "2019") {
    " in 2019"} else {paste0(" from ", gsub("_", " to ", period_name))
    }
  
  p <- ggplot(world_data) +
    geom_sf(aes(fill = Proportion_median)) +
    scale_fill_viridis_c(
      option = "viridis",
      limits = c(0.2,1),
      breaks = c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
      labels = c("<0.3", "0.4", "0.5", "0.6", "0.7", "0.8", ">0.9"),
      oob = scales::squish,
      name = paste0("Antibiotic treatment probability for AFI", time_label),
      na.value = "gray"
    ) +
    theme_minimal() +
    labs(
      title = paste("Antibiotic treatment probability for AFI", time_label),
    ) +
    theme(
      legend.position = c(0.03, 0.15),        
      legend.justification = c(0, 0),         
      legend.title = element_text(size = 6), 
      legend.text = element_text(size = 6),
      legend.key.height = unit(0.5, "cm"),
      legend.key.width = unit(0.3, "cm"),
      plot.title = element_text(size = 8, face = "bold"),
      plot.subtitle = element_text(size = 8)
    )
  
  plot_list_probability_Fever[[period_name]] <- p
}

tiff("outputs/table/by_year_malaria/SA_probability_Fever_2004_2006.tiff", width = 7.5, height = 5, units = "in", res = 300)
plot_list_probability_Fever[["2004_2006"]]
dev.off()
tiff("outputs/table/by_year_malaria/SA_probability_Fever_2009_2011.tiff", width = 7.5, height = 5, units = "in", res = 300)
plot_list_probability_Fever[["2009_2011"]]
dev.off()
tiff("outputs/table/by_year_malaria/SA_probability_Fever_2014_2016.tiff", width = 7.5, height = 5, units = "in", res = 300)
plot_list_probability_Fever[["2014_2016"]]
dev.off()
tiff("outputs/table/by_year_malaria/SA_probability_Fever_2019.tiff", width = 7.5, height = 5, units = "in", res = 300)
plot_list_probability_Fever[["2019"]]
dev.off()





combined_plot <- (
  plot_list_probability_LRI[["2019"]] + 
    ggtitle("A. Global antibiotic access for LRTI in 2019") + 
    theme(plot.margin = unit(c(0, 0, 0, 0), "pt"),
          plot.title = element_text(hjust = +0.03),
          panel.grid = element_blank()) 
) / plot_spacer() / (
  plot_list_probability_Fever[["2019"]] +
    ggtitle("B. Global antibiotic access for AFI in 2019") + 
    theme(plot.margin = unit(c(0, 0, 0, 0), "pt"),
          plot.title = element_text(hjust = +0.03),
          panel.grid = element_blank()) 
) +
  plot_layout(heights = c(1, 0.04, 1))

# 2. Save the combined plot to PDF
ggsave("figures/Figure_S4.pdf", plot = combined_plot, width = 7.5, height = 7.5)




#----LRTI-----
periods <- list(
  "2004_2006" = 2004:2006, "2009_2011" = 2009:2011,
  "2014_2016" = 2014:2016, "2019" = 2019)

plot_list_DDD_LRI <- list()

for (period_name in names(periods)) {
  years <- periods[[period_name]]
  
  filtered_data <- LRI_joint  %>% filter(Year %in% years) %>% group_by(Country) %>% 
    summarize(
      Per_treated_median = sum(Total_DDD)/sum(Treated_incidence), 
      .groups = "drop")
  world_data <- world %>% left_join(filtered_data, by = c("admin" = "Country"))
  time_label <- if (period_name == "2019") {
    " in 2019"} else {paste0(" from ", gsub("_", " to ", period_name))
    }
  
  p <- ggplot(world_data) +
    geom_sf(aes(fill = Per_treated_median)) +
    scale_fill_viridis_c(
      option = "inferno",
      limits = c(0, 600.1),
      breaks = c(0, 100, 200, 300, 400, 500, 600),
      labels = c("0","100","200","300","400","500", ">600"),
      oob = scales::squish,
      name = paste0("DDD per treated LRTI \n", time_label),
      na.value = "gray"
    ) +
    theme_minimal() +
    labs(
      title = paste("DDD per treated LRTI \n", time_label)
    ) +
    theme(
      legend.position = c(0.03, 0.15),        
      legend.justification = c(0, 0),         
      legend.title = element_text(size = 8, face = "bold"), 
      legend.text = element_text(size = 8),
      legend.key.height = unit(0.5, "cm"),
      legend.key.width = unit(0.3, "cm"),
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  
  plot_list_DDD_LRI[[period_name]] <- p
}

tiff("outputs/table/by_year_malaria/SA_DDD_per_treated_LRTI_2004_2006.tiff", width = 7.5, height = 5, units = "in", res = 300)
plot_list_DDD_LRI[["2004_2006"]]
dev.off()
tiff("outputs/table/by_year_malaria/SA_DDD_per_treated_LRTI_2009_2011.tiff", width = 7.5, height = 5, units = "in", res = 300)
plot_list_DDD_LRI[["2009_2011"]]
dev.off()
tiff("outputs/table/by_year_malaria/SA_DDD_per_treated_LRTI_2014_2016.tiff", width = 7.5, height = 5, units = "in", res = 300)
plot_list_DDD_LRI[["2014_2016"]]
dev.off()
tiff("outputs/table/by_year_malaria/SA_DDD_per_treated_LRTI_2019.tiff", width = 7.5, height = 5, units = "in", res = 300)
plot_list_DDD_LRI[["2019"]]
dev.off()



#----Fever-----
periods <- list(
  "2004_2006" = 2004:2006, "2009_2011" = 2009:2011,
  "2014_2016" = 2014:2016, "2019" = 2019)

plot_list_DDD_Fever <- list()

for (period_name in names(periods)) {
  years <- periods[[period_name]]
  
  filtered_data <- Fever_joint  %>% filter(Year %in% years) %>% group_by(Country) %>% 
    summarize(
      Per_treated_median = round(sum(Total_DDD)/sum(Treated_incidence),digits = 2), 
      .groups = "drop")
  
  world_data <- world %>% left_join(filtered_data, by = c("admin" = "Country"))
  time_label <- if (period_name == "2019") {
    " in 2019"} else {paste0(" from ", gsub("_", " to ", period_name))
    }
  
  p <- ggplot(world_data) +
    geom_sf(aes(fill = Per_treated_median)) +
    scale_fill_viridis_c(
      option = "inferno",
      limits = c(min(Fever_joint$Per_treated_median), max(Fever_joint$Per_treated_median)+0.01),
      breaks = c(0, 2, 4, 6, 8, 10, 12, 14),
      labels = c("0","2","4","6","8","10","12",">14"),
      oob = scales::squish,
      name = paste0("DDD per treated Fever without LRTI \n", time_label),
      na.value = "gray"
    ) +
    theme_minimal() +
    labs(
      title = paste("DDD per treated Fever without LRTI \n", time_label)
    ) +
    theme(
      legend.position = c(0.03, 0.15),        
      legend.justification = c(0, 0),         
      legend.title = element_text(size = 8, face = "bold"), 
      legend.text = element_text(size = 8),
      legend.key.height = unit(0.5, "cm"),
      legend.key.width = unit(0.3, "cm"),
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  
  plot_list_DDD_Fever[[period_name]] <- p
}

tiff("outputs/table/by_year_malaria/SA_DDD_per_treated_Fever_2004_2006.tiff", width = 7.5, height = 5, units = "in", res = 300)
plot_list_DDD_Fever[["2004_2006"]]
dev.off()
tiff("outputs/table/by_year_malaria/SA_DDD_per_treated_Fever_2009_2011.tiff", width = 7.5, height = 5, units = "in", res = 300)
plot_list_DDD_Fever[["2009_2011"]]
dev.off()
tiff("outputs/table/by_year_malaria/SA_DDD_per_treated_Fever_2014_2016.tiff", width = 7.5, height = 5, units = "in", res = 300)
plot_list_DDD_Fever[["2014_2016"]]
dev.off()
tiff("outputs/table/by_year_malaria/SA_DDD_per_treated_Fever_2019.tiff", width = 7.5, height = 5, units = "in", res = 300)
plot_list_DDD_Fever[["2019"]]
dev.off()

#






######## 9. Figure.2 ########

#### 0. prep ####
rm(list = ls())
gc()

requiredPackages = c("dplyr", "tidyr", "data.table", "stringr")
lib.path <- .libPaths()
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p, lib = lib.path)
  library(p,character.only = TRUE)
}

n_sim <- 10000
set.seed(123)


# confirm upstream extrapolation models use PC1-PC2

coef_df_LRI <- readRDS("outputs/sensitivity_analysis/SA_coef_df_LRI.rds")
coef_df_Fever <- readRDS("outputs/sensitivity_analysis/SA_coef_df_Fever.rds")
coef_df_LRI_DDD <- readRDS("outputs/sensitivity_analysis/SA_coef_df_LRI_DDD.rds")
coef_df_Fever_DDD <- readRDS("outputs/sensitivity_analysis/SA_coef_df_Fever_DDD.rds")
coef_df_all_DDD <- readRDS("outputs/sensitivity_analysis/SA_coef_df_all_DDD.rds")

expected_coef_names <- c("(Intercept)", "PC1", "PC2")

stopifnot(
  identical(names(coef_df_LRI), expected_coef_names),
  identical(names(coef_df_Fever), expected_coef_names),
  identical(names(coef_df_LRI_DDD), expected_coef_names),
  identical(names(coef_df_Fever_DDD), expected_coef_names),
  identical(names(coef_df_all_DDD), expected_coef_names)
)

rm(
  coef_df_LRI,
  coef_df_Fever,
  coef_df_LRI_DDD,
  coef_df_Fever_DDD,
  coef_df_all_DDD,
  expected_coef_names
)


# read the plot data generated from the current PC1-PC2 models: probability

LRI_year_proportion      <- read.csv("outputs/plot_df/SA_1.3_plot_Ab_LRI_year_proportion.csv")
LRI_year_proportion_L    <- read.csv("outputs/plot_df/SA_1.3_plot_Ab_LRI_year_proportion_L.csv")
LRI_year_proportion_LM   <- read.csv("outputs/plot_df/SA_1.3_plot_Ab_LRI_year_proportion_LM.csv")
LRI_year_proportion_UM   <- read.csv("outputs/plot_df/SA_1.3_plot_Ab_LRI_year_proportion_UM.csv")
LRI_year_proportion_malaria_0      <- read.csv("outputs/plot_df/SA_1.4_plot_Ab_LRI_year_proportion_malaria_0.csv")
LRI_year_proportion_malaria_1      <- read.csv("outputs/plot_df/SA_1.4_plot_Ab_LRI_year_proportion_malaria_1.csv")
Fever_year_proportion    <- read.csv("outputs/plot_df/SA_1.3_plot_Ab_Fever_year_proportion.csv")
Fever_year_proportion_mal0_L  <- read.csv("outputs/plot_df/SA_1.3_plot_Ab_Fever_year_proportion_mal0_L.csv")
Fever_year_proportion_mal0_LM <- read.csv("outputs/plot_df/SA_1.3_plot_Ab_Fever_year_proportion_mal0_LM.csv")
Fever_year_proportion_mal0_UM <- read.csv("outputs/plot_df/SA_1.3_plot_Ab_Fever_year_proportion_mal0_UM.csv")
Fever_year_proportion_mal1_L  <- read.csv("outputs/plot_df/SA_1.3_plot_Ab_Fever_year_proportion_mal1_L.csv")
Fever_year_proportion_mal1_LM <- read.csv("outputs/plot_df/SA_1.3_plot_Ab_Fever_year_proportion_mal1_LM.csv")
Fever_year_proportion_mal1_UM <- read.csv("outputs/plot_df/SA_1.3_plot_Ab_Fever_year_proportion_mal1_UM.csv")
Fever_year_proportion_malaria_0      <- read.csv("outputs/plot_df/SA_1.4_plot_Ab_Fever_year_proportion_malaria_0.csv")
Fever_year_proportion_malaria_1      <- read.csv("outputs/plot_df/SA_1.4_plot_Ab_Fever_year_proportion_malaria_1.csv")

for (i in 1:7) { # GBD_region
  path_LRI <- paste0("outputs/plot_df/SA_1.5_plot_Ab_LRI_year_proportion_GBD_region_", i, ".csv")
  obj_LRI <- paste0("LRI_year_proportion_GBD_region_", i)
  assign(obj_LRI, read.csv(path_LRI, stringsAsFactors = FALSE))
}

Malaria_level <- c(0,1)

for (mal in Malaria_level) {
  for (i in 1:7) {
    obj_Fever <- paste0("Fever_year_proportion_GBD_region_", mal,"_", i)
    path_Fever <- paste0("outputs/plot_df/SA_1.5_plot_Ab_Fever_year_proportion_GBD_region_", mal,"_", i,".csv")
    assign(obj_Fever, read.csv(path_Fever, stringsAsFactors = FALSE))
  }
}

LRI_year_DDD_per_treated      <- read.csv("outputs/plot_df/SA_2.3_plot_Ab_LRI_year_DDD_per_treated.csv")
LRI_year_DDD_per_treated_L    <- read.csv("outputs/plot_df/SA_2.3_plot_Ab_LRI_year_DDD_per_treated_L.csv")
LRI_year_DDD_per_treated_LM   <- read.csv("outputs/plot_df/SA_2.3_plot_Ab_LRI_year_DDD_per_treated_LM.csv")
LRI_year_DDD_per_treated_UM   <- read.csv("outputs/plot_df/SA_2.3_plot_Ab_LRI_year_DDD_per_treated_UM.csv")
Fever_year_DDD_per_treated    <- read.csv("outputs/plot_df/SA_2.3_plot_Ab_Fever_year_DDD_per_treated.csv")
Fever_year_DDD_per_treated_L  <- read.csv("outputs/plot_df/SA_2.3_plot_Ab_Fever_year_DDD_per_treated_L.csv")
Fever_year_DDD_per_treated_LM <- read.csv("outputs/plot_df/SA_2.3_plot_Ab_Fever_year_DDD_per_treated_LM.csv")
Fever_year_DDD_per_treated_UM <- read.csv("outputs/plot_df/SA_2.3_plot_Ab_Fever_year_DDD_per_treated_UM.csv")
all_year_DDD_per_treated      <- read.csv("outputs/plot_df/SA_2.3_plot_Ab_all_year_DDD_per_treated.csv")
all_year_DDD_per_treated_L    <- read.csv("outputs/plot_df/SA_2.3_plot_Ab_all_year_DDD_per_treated_L.csv")
all_year_DDD_per_treated_LM   <- read.csv("outputs/plot_df/SA_2.3_plot_Ab_all_year_DDD_per_treated_LM.csv")
all_year_DDD_per_treated_UM   <- read.csv("outputs/plot_df/SA_2.3_plot_Ab_all_year_DDD_per_treated_UM.csv")




# ============================================================
# Export source data for sensitivity-analysis Figure 2
# ============================================================

dir.create(
  "outputs/plot_objects",
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# Helper function
# ------------------------------------------------------------

format_SA_fig2 <- function(
    df,
    panel,
    syndrome,
    malaria_status = "Not applicable",
    stratification = "Overall",
    group = "Overall"
) {
  
  df %>%
    transmute(
      Panel = panel,
      Syndrome = syndrome,
      Malaria_status = malaria_status,
      Stratification = stratification,
      Group = group,
      Year = Year,
      Median = Proportion_Median,
      Lower = Proportion_Lower,
      Upper = Proportion_Upper
    )
}


# ------------------------------------------------------------
# Panels A-C: overall antibiotic access
# ------------------------------------------------------------

fig_A <- format_SA_fig2(
  LRI_year_proportion,
  panel = "A",
  syndrome = "LRTI"
)

fig_B <- format_SA_fig2(
  Fever_year_proportion_malaria_0,
  panel = "B",
  syndrome = "AFI",
  malaria_status = "Malaria non-endemic"
)

fig_C <- format_SA_fig2(
  Fever_year_proportion_malaria_1,
  panel = "C",
  syndrome = "AFI",
  malaria_status = "Malaria endemic"
)


# ------------------------------------------------------------
# Panels D-F: antibiotic access by income
# ------------------------------------------------------------

income_names <- c(
  "L"  = "Low income",
  "LM" = "Lower middle income",
  "UM" = "Upper middle income"
)


fig_D <- bind_rows(
  format_SA_fig2(
    LRI_year_proportion_L,
    "D", "LRTI",
    stratification = "Income",
    group = income_names["L"]
  ),
  
  format_SA_fig2(
    LRI_year_proportion_LM,
    "D", "LRTI",
    stratification = "Income",
    group = income_names["LM"]
  ),
  
  format_SA_fig2(
    LRI_year_proportion_UM,
    "D", "LRTI",
    stratification = "Income",
    group = income_names["UM"]
  )
)


fig_E <- bind_rows(
  format_SA_fig2(
    Fever_year_proportion_mal0_L,
    "E", "AFI",
    malaria_status = "Malaria non-endemic",
    stratification = "Income",
    group = income_names["L"]
  ),
  
  format_SA_fig2(
    Fever_year_proportion_mal0_LM,
    "E", "AFI",
    malaria_status = "Malaria non-endemic",
    stratification = "Income",
    group = income_names["LM"]
  ),
  
  format_SA_fig2(
    Fever_year_proportion_mal0_UM,
    "E", "AFI",
    malaria_status = "Malaria non-endemic",
    stratification = "Income",
    group = income_names["UM"]
  )
)


fig_F <- bind_rows(
  format_SA_fig2(
    Fever_year_proportion_mal1_L,
    "F", "AFI",
    malaria_status = "Malaria endemic",
    stratification = "Income",
    group = income_names["L"]
  ),
  
  format_SA_fig2(
    Fever_year_proportion_mal1_LM,
    "F", "AFI",
    malaria_status = "Malaria endemic",
    stratification = "Income",
    group = income_names["LM"]
  ),
  
  format_SA_fig2(
    Fever_year_proportion_mal1_UM,
    "F", "AFI",
    malaria_status = "Malaria endemic",
    stratification = "Income",
    group = income_names["UM"]
  )
)


# ------------------------------------------------------------
# Panels G-I: antibiotic access by GBD region
# ------------------------------------------------------------

GBD_region_names <- c(
  "1" = "North Africa and Middle East",
  "2" = "Central Europe, Eastern Europe, and Central Asia",
  "3" = "Sub-Saharan Africa",
  "5" = "South Asia",
  "6" = "Latin America and Caribbean",
  "7" = "Southeast Asia, East Asia, and Oceania"
)


fig_G <- bind_rows(
  lapply(c(1, 2, 3, 5, 6, 7), function(i) {
    
    df <- get(
      paste0("LRI_year_proportion_GBD_region_", i)
    )
    
    format_SA_fig2(
      df,
      panel = "G",
      syndrome = "LRTI",
      stratification = "GBD region",
      group = unname(
        GBD_region_names[as.character(i)]
      )
    )
  })
)


fig_H <- bind_rows(
  lapply(c(1, 2, 3, 5, 6, 7), function(i) {
    
    df <- get(
      paste0(
        "Fever_year_proportion_GBD_region_0_",
        i
      )
    )
    
    format_SA_fig2(
      df,
      panel = "H",
      syndrome = "AFI",
      malaria_status = "Malaria non-endemic",
      stratification = "GBD region",
      group = unname(
        GBD_region_names[as.character(i)]
      )
    )
  })
)


fig_I <- bind_rows(
  lapply(c(1, 2, 3, 5, 6, 7), function(i) {
    
    df <- get(
      paste0(
        "Fever_year_proportion_GBD_region_1_",
        i
      )
    )
    
    format_SA_fig2(
      df,
      panel = "I",
      syndrome = "AFI",
      malaria_status = "Malaria endemic",
      stratification = "GBD region",
      group = unname(
        GBD_region_names[as.character(i)]
      )
    )
  })
)


# ------------------------------------------------------------
# Combine all panels
# ------------------------------------------------------------

SA_Figure_2_source_data <- bind_rows(
  fig_A,
  fig_B,
  fig_C,
  fig_D,
  fig_E,
  fig_F,
  fig_G,
  fig_H,
  fig_I
) %>%
  arrange(
    Panel,
    Group,
    Year
  )


# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

write.csv(
  SA_Figure_2_source_data,
  "outputs/plot_objects/SA_Figure_5_source_data.csv",
  row.names = FALSE
)





pdf("figures/Figure_S5.pdf", width = 7.5, height = 7)

layout(matrix(c(1,1,1,1,1,
                0,4,5,6,0,
                2,2,2,2,2,
                0,7,8,9,10,
                3,3,3,3,3,
                0,11,12,13,14), nrow = 6, byrow = TRUE), widths = c(0.1 ,1, 1, 1, 0.8), heights = c(0.1,1,0.1,1,0.1,1))




#### 1. (title 1) ####

par(mar = c(0, 0, 0, 0))
plot.new()
text(x = -0.035, y = 0.3, labels = "Antibiotic access ", 
     adj = 0, cex = 0.7, font = 2, xpd = TRUE)

#### 2. (title 2)  ####
par(mar = c(0, 0, 0, 0))
plot.new()
text(x = -0.035, y = 0.3, labels = "Antibiotic access by Income ", 
     adj = 0, cex = 0.7, font = 2, xpd = TRUE)

#### 3. (title 3) ####
par(mar = c(0, 0, 0, 0))
plot.new()
text(x = -0.035, y = 0.3, labels = "Antibiotic access by GBD region ", 
     adj = 0, cex = 0.7, font = 2, xpd = TRUE)


#### 4. (Fig.A) ####

xmin <- 2003; xmax <- 2019
ymin <- 0; ymax <- 1

par(mar = c(2.5, 2.5, 2, 1))  
par(mgp = c(3, 0.3, 0))     
par(tck = -0.025)
par(lwd = 1)

plot(Proportion_Median ~ Year, 
     data =  LRI_year_proportion,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")



axis(1, at = seq(xmin, xmax, by = 4), labels = FALSE, las = 2, pos = ymin)
text(x = seq(xmin, xmax, by = 4), y = ymin-(ymax-ymin)*0.065,  
     labels = format(round(seq(xmin, xmax, by = (xmax-xmin)/6*0+4), digits = 0), nsmall = 0), xpd = TRUE, cex = 0.65)

axis(2,  at = seq(ymin, ymax, by = (ymax - ymin)/4), labels = FALSE, las = 2, pos = xmin)
text(y = seq(ymin, ymax, by = (ymax-ymin)/4), x = xmin-0.04*(xmax-xmin),  
     labels = format(round(seq(ymin, ymax, by = (ymax-ymin)/4), digits = 2), nsmall = 1), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines(Proportion_Median ~ Year, data = LRI_year_proportion, col = "#0072B2", lwd = 1.5)

polygon(c(LRI_year_proportion$Year, rev(LRI_year_proportion$Year)),
        c(LRI_year_proportion$Proportion_Lower, rev(LRI_year_proportion$Proportion_Upper)),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)
mtext(expression("Antibiotic treatment probability"), side = 2, line = 2, cex = 0.5)  
mtext("Year", side = 1, line = 1, cex = 0.5)
mtext("A. LRTI", at = xmin - 0.3*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.5, font = 1)






#### 5. (Fig.B) ####


plot(Proportion_Median ~ Year, 
     data =  Fever_year_proportion,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")



axis(1, at = seq(xmin, xmax, by = 4), labels = FALSE, las = 2, pos = ymin)
text(x = seq(xmin, xmax, by = 4), y = ymin-(ymax-ymin)*0.065,  
     labels = format(round(seq(xmin, xmax, by = (xmax-xmin)/6*0+4), digits = 0), nsmall = 0), xpd = TRUE, cex = 0.65)

axis(2,  at = seq(ymin, ymax, by = (ymax - ymin)/4), labels = FALSE, las = 2, pos = xmin)
text(y = seq(ymin, ymax, by = (ymax-ymin)/4), x = xmin-0.04*(xmax-xmin),  
     labels = format(round(seq(ymin, ymax, by = (ymax-ymin)/4), digits = 2), nsmall = 1), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)


lines(Proportion_Median ~ Year, data = Fever_year_proportion_malaria_0, col = "#0072B2", lwd = 1.5)

polygon(c(Fever_year_proportion_malaria_0$Year, rev(Fever_year_proportion_malaria_0$Year)),
        c(Fever_year_proportion_malaria_0$Proportion_Lower, rev(Fever_year_proportion_malaria_0$Proportion_Upper)),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)

mtext("Year", side = 1, line = 1, cex = 0.5)
mtext("B. AFI (malaria non-endemic countries)", at = xmin - 0.3*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.5, font = 1)





#### 6. (Fig.C) ####


plot(Proportion_Median ~ Year, 
     data =  Fever_year_proportion,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")

axis(1, at = seq(xmin, xmax, by = 4), labels = FALSE, las = 2, pos = ymin)
text(x = seq(xmin, xmax, by = 4), y = ymin-(ymax-ymin)*0.065,  
     labels = format(round(seq(xmin, xmax, by = (xmax-xmin)/6*0+4), digits = 0), nsmall = 0), xpd = TRUE, cex = 0.65)

axis(2,  at = seq(ymin, ymax, by = (ymax - ymin)/4), labels = FALSE, las = 2, pos = xmin)
text(y = seq(ymin, ymax, by = (ymax-ymin)/4), x = xmin-0.04*(xmax-xmin),  
     labels = format(round(seq(ymin, ymax, by = (ymax-ymin)/4), digits = 2), nsmall = 1), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)


lines(Proportion_Median ~ Year, data = Fever_year_proportion_malaria_1,  col = "#0072B2", lwd = 1.5)

polygon(c(Fever_year_proportion_malaria_1$Year, rev(Fever_year_proportion_malaria_1$Year)),
        c(Fever_year_proportion_malaria_1$Proportion_Lower, rev(Fever_year_proportion_malaria_1$Proportion_Upper)),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)

mtext("Year", side = 1, line = 1, cex = 0.5)
mtext("C. AFI (malaria-endemic countries)", at = xmin - 0.3*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.5, font = 1)





#### 7. (Fig.D) #####


plot(Proportion_Median ~ Year, 
     data =  LRI_year_proportion_L,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")



axis(1, at = seq(xmin, xmax, by = 4), labels = FALSE, las = 2, pos = ymin)
text(x = seq(xmin, xmax, by = 4), y = ymin-(ymax-ymin)*0.065,  
     labels = format(round(seq(xmin, xmax, by = (xmax-xmin)/6*0+4), digits = 0), nsmall = 0), xpd = TRUE, cex = 0.65)

axis(2,  at = seq(ymin, ymax, by = (ymax - ymin)/4), labels = FALSE, las = 2, pos = xmin)
text(y = seq(ymin, ymax, by = (ymax-ymin)/4), x = xmin-0.04*(xmax-xmin),  
     labels = format(round(seq(ymin, ymax, by = (ymax-ymin)/4), digits = 2), nsmall = 1), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines(Proportion_Median ~ Year, data = LRI_year_proportion_L,  col = "#CC79A7", lwd = 1.5)
lines(Proportion_Median ~ Year, data = LRI_year_proportion_LM, col = "#009E73", lwd = 1.5)
lines(Proportion_Median ~ Year, data = LRI_year_proportion_UM, col = "#0072B2", lwd = 1.5)

# Add shaded area for confidence interval
polygon(c(LRI_year_proportion_L$Year, rev(LRI_year_proportion_L$Year)),
        c(LRI_year_proportion_L$Proportion_Lower, rev(LRI_year_proportion_L$Proportion_Upper)),
        col = rgb(204/255, 121/255, 167/255, 0.25), border = NA)

polygon(c(LRI_year_proportion_LM$Year, rev(LRI_year_proportion_LM$Year)),
        c(LRI_year_proportion_LM$Proportion_Lower, rev(LRI_year_proportion_LM$Proportion_Upper)),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

polygon(c(LRI_year_proportion_UM$Year, rev(LRI_year_proportion_UM$Year)),
        c(LRI_year_proportion_UM$Proportion_Lower, rev(LRI_year_proportion_UM$Proportion_Upper)),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)
mtext(expression("Antibiotic treatment probability"), side = 2, line = 2, cex = 0.5)  
mtext("Year", side = 1, line = 1, cex = 0.5)
mtext("D. LRTI", at = xmin - 0.3*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.5, font = 1)


#### 8. (Fig.E) #####


plot(Proportion_Median ~ Year, 
     data =  Fever_year_proportion_mal0_L,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")



axis(1, at = seq(xmin, xmax, by = 4), labels = FALSE, las = 2, pos = ymin)
text(x = seq(xmin, xmax, by = 4), y = ymin-(ymax-ymin)*0.065,  
     labels = format(round(seq(xmin, xmax, by = (xmax-xmin)/6*0+4), digits = 0), nsmall = 0), xpd = TRUE, cex = 0.65)

axis(2,  at = seq(ymin, ymax, by = (ymax - ymin)/4), labels = FALSE, las = 2, pos = xmin)
text(y = seq(ymin, ymax, by = (ymax-ymin)/4), x = xmin-0.04*(xmax-xmin),  
     labels = format(round(seq(ymin, ymax, by = (ymax-ymin)/4), digits = 2), nsmall = 1), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines(Proportion_Median ~ Year, data = Fever_year_proportion_mal0_L,  col = "#CC79A7", lwd = 1.5)
lines(Proportion_Median ~ Year, data = Fever_year_proportion_mal0_LM, col = "#009E73", lwd = 1.5)
lines(Proportion_Median ~ Year, data = Fever_year_proportion_mal0_UM, col = "#0072B2", lwd = 1.5)

# Add shaded area for confidence interval
polygon(c(Fever_year_proportion_mal0_L$Year, rev(Fever_year_proportion_mal0_L$Year)),
        c(Fever_year_proportion_mal0_L$Proportion_Lower, rev(Fever_year_proportion_mal0_L$Proportion_Upper)),
        col = rgb(204/255, 121/255, 167/255, 0.25), border = NA)

polygon(c(Fever_year_proportion_mal0_LM$Year, rev(Fever_year_proportion_mal0_LM$Year)),
        c(Fever_year_proportion_mal0_LM$Proportion_Lower, rev(Fever_year_proportion_mal0_LM$Proportion_Upper)),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

polygon(c(Fever_year_proportion_mal0_UM$Year, rev(Fever_year_proportion_mal0_UM$Year)),
        c(Fever_year_proportion_mal0_UM$Proportion_Lower, rev(Fever_year_proportion_mal0_UM$Proportion_Upper)),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)
mtext("Year", side = 1, line = 1, cex = 0.5)
mtext("E. AFI (malaria non-endemic countries)", at = xmin - 0.3*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.5, font = 1)




#### 9. (Fig.F) #####


plot(Proportion_Median ~ Year, 
     data =  Fever_year_proportion_mal1_L,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")



axis(1, at = seq(xmin, xmax, by = 4), labels = FALSE, las = 2, pos = ymin)
text(x = seq(xmin, xmax, by = 4), y = ymin-(ymax-ymin)*0.065,  
     labels = format(round(seq(xmin, xmax, by = (xmax-xmin)/6*0+4), digits = 0), nsmall = 0), xpd = TRUE, cex = 0.65)

axis(2,  at = seq(ymin, ymax, by = (ymax - ymin)/4), labels = FALSE, las = 2, pos = xmin)
text(y = seq(ymin, ymax, by = (ymax-ymin)/4), x = xmin-0.04*(xmax-xmin),  
     labels = format(round(seq(ymin, ymax, by = (ymax-ymin)/4), digits = 2), nsmall = 1), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines(Proportion_Median ~ Year, data = Fever_year_proportion_mal1_L,  col = "#CC79A7", lwd = 1.5)
lines(Proportion_Median ~ Year, data = Fever_year_proportion_mal1_LM, col = "#009E73", lwd = 1.5)
lines(Proportion_Median ~ Year, data = Fever_year_proportion_mal1_UM, col = "#0072B2", lwd = 1.5)

# Add shaded area for confidence interval
polygon(c(Fever_year_proportion_mal1_L$Year, rev(Fever_year_proportion_mal1_L$Year)),
        c(Fever_year_proportion_mal1_L$Proportion_Lower, rev(Fever_year_proportion_mal1_L$Proportion_Upper)),
        col = rgb(204/255, 121/255, 167/255, 0.25), border = NA)

polygon(c(Fever_year_proportion_mal1_LM$Year, rev(Fever_year_proportion_mal1_LM$Year)),
        c(Fever_year_proportion_mal1_LM$Proportion_Lower, rev(Fever_year_proportion_mal1_LM$Proportion_Upper)),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

polygon(c(Fever_year_proportion_mal1_UM$Year, rev(Fever_year_proportion_mal1_UM$Year)),
        c(Fever_year_proportion_mal1_UM$Proportion_Lower, rev(Fever_year_proportion_mal1_UM$Proportion_Upper)),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)
mtext("Year", side = 1, line = 1, cex = 0.5)
mtext("F. AFI (malaria-endemic countries)", at = xmin - 0.3*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.5, font = 1)




#### 10. (legend) ####
par(mar = c(0, 0, 0, 0))
plot.new()

legend_labels <- c("Upper Middle Income Country", "Lower Middle Income Country", "Low Income Country")
legend_colors <- c("#0072B2", "#009E73", "#CC79A7")
n <- length(legend_labels)

# Coordinates for the custom legend
x0 <- 0.0  
x1 <- 0.07  
y_start <- 0.6
line_spacing <- 0.07

# Draw each line, label, and shaded polygon
for (i in seq_len(n)) {
  y <- y_start - (i - 1) * line_spacing
  

  polygon(
    x = c(x0-0.01, x1+0.01, x1+0.01, x0-0.01),
    y = c(y - 0.02, y - 0.02, y + 0.02, y + 0.02),
    col = adjustcolor(legend_colors[i], alpha.f = 0.25),
    border = NA
  )
  
 
  segments(x0, y, x1, y, col = legend_colors[i], lwd = 1.5)
  

  text(x = x1 + 0.04, y = y, labels = legend_labels[i], adj = 0, cex = 0.56)
}





#### 11. (Fig. G) ####

par(mar = c(2.5, 2.5, 2, 1))  
par(mgp = c(3, 0.3, 0))     
par(tck = -0.025)
par(lwd = 1)

plot(Proportion_Median ~ Year, 
     data =  LRI_year_proportion_GBD_region_2, 
     xlab = "", 
     ylab = "", 
     pch = 16, 
     col = rgb(0.8, 0.8, 0.8, 0), 
     xlim = c(xmin,xmax), 
     ylim = c(ymin,ymax), 
     lwd = "1", 
     axes = F,
     ann = F, 
     main = "")

axis(1, at = seq(xmin, xmax, by = 4), labels = FALSE, las = 2, pos = ymin)
text(x = seq(xmin, xmax, by = 4), y = ymin-(ymax-ymin)*0.065,  
     labels = format(round(seq(xmin, xmax, by = (xmax-xmin)/6*0+4), digits = 0), nsmall = 0), xpd = TRUE, cex = 0.65)

axis(2,  at = seq(ymin, ymax, by = (ymax - ymin)/4), labels = FALSE, las = 2, pos = xmin)
text(y = seq(ymin, ymax, by = (ymax-ymin)/4), x = xmin-0.04*(xmax-xmin),  
     labels = format(round(seq(ymin, ymax, by = (ymax-ymin)/4), digits = 2), nsmall = 1), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines(Proportion_Median ~ Year, data = LRI_year_proportion_GBD_region_2, col = "#CC79A7", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_2$Year, rev(LRI_year_proportion_GBD_region_2$Year)),
        c(LRI_year_proportion_GBD_region_2$Proportion_Lower, rev(LRI_year_proportion_GBD_region_2$Proportion_Upper)),
        col = rgb(204/255, 121/255, 167/255, 0.25), border = NA)

lines(Proportion_Median ~ Year, data = LRI_year_proportion_GBD_region_7, col = "#D55E00", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_7$Year, rev(LRI_year_proportion_GBD_region_7$Year)),
        c(LRI_year_proportion_GBD_region_7$Proportion_Lower, rev(LRI_year_proportion_GBD_region_7$Proportion_Upper)),
        col = rgb(213/255, 94/255, 0/255, 0.25), border = NA)

lines(Proportion_Median ~ Year, data = LRI_year_proportion_GBD_region_6, col = "#E69F00", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_6$Year, rev(LRI_year_proportion_GBD_region_6$Year)),
        c(LRI_year_proportion_GBD_region_6$Proportion_Lower, rev(LRI_year_proportion_GBD_region_6$Proportion_Upper)),
        col = rgb(230/255, 159/255, 0/255, 0.25), border = NA)

lines(Proportion_Median ~ Year, data = LRI_year_proportion_GBD_region_5, col = "#A7B97F", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_5$Year, rev(LRI_year_proportion_GBD_region_5$Year)),
        c(LRI_year_proportion_GBD_region_5$Proportion_Lower, rev(LRI_year_proportion_GBD_region_5$Proportion_Upper)),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

lines(Proportion_Median ~ Year, data = LRI_year_proportion_GBD_region_1, col = "#009E73", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_1$Year, rev(LRI_year_proportion_GBD_region_1$Year)),
        c(LRI_year_proportion_GBD_region_1$Proportion_Lower, rev(LRI_year_proportion_GBD_region_1$Proportion_Upper)),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

lines(Proportion_Median ~ Year, data = LRI_year_proportion_GBD_region_3, col = "#0072B2", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_3$Year, rev(LRI_year_proportion_GBD_region_3$Year)),
        c(LRI_year_proportion_GBD_region_3$Proportion_Lower, rev(LRI_year_proportion_GBD_region_3$Proportion_Upper)),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)

mtext("G. LRTI", at = xmin - 0.3*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.5, font = 1)
mtext(expression("Antibiotic treatment probability"), side = 2, line = 2, cex = 0.5)  
mtext("Year", side = 1, line = 1, cex = 0.5)











plot(Proportion_Median ~ Year, 
     data = Fever_year_proportion_GBD_region_0_2, 
     xlab = "", 
     ylab = "", 
     pch = 16, 
     col = rgb(0.8, 0.8, 0.8, 0), 
     xlim = c(xmin,xmax), 
     ylim = c(ymin,ymax), 
     lwd = "1", 
     axes = F,
     ann = F, 
     main = "")

axis(1, at = seq(xmin, xmax, by = 4), labels = FALSE, las = 2, pos = ymin)
text(x = seq(xmin, xmax, by = 4), y = ymin-(ymax-ymin)*0.065,  
     labels = format(round(seq(xmin, xmax, by = (xmax-xmin)/6*0+4), digits = 0), nsmall = 0), xpd = TRUE, cex = 0.65)

axis(2,  at = seq(ymin, ymax, by = (ymax - ymin)/4), labels = FALSE, las = 2, pos = xmin)
text(y = seq(ymin, ymax, by = (ymax-ymin)/4), x = xmin-0.04*(xmax-xmin),  
     labels = format(round(seq(ymin, ymax, by = (ymax-ymin)/4), digits = 2), nsmall = 1), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)


lines(Proportion_Median ~ Year, data = Fever_year_proportion_GBD_region_0_2, col = "#CC79A7", lwd = 1.5)
polygon(c(Fever_year_proportion_GBD_region_0_2$Year, rev(Fever_year_proportion_GBD_region_0_2$Year)),
        c(Fever_year_proportion_GBD_region_0_2$Proportion_Lower, rev(Fever_year_proportion_GBD_region_0_2$Proportion_Upper)),
        col = rgb(204/255, 121/255, 167/255, 0.25), border = NA)

lines(Proportion_Median ~ Year, data = Fever_year_proportion_GBD_region_0_7, col = "#D55E00", lwd = 1.5)
polygon(c(Fever_year_proportion_GBD_region_0_7$Year, rev(Fever_year_proportion_GBD_region_0_7$Year)),
        c(Fever_year_proportion_GBD_region_0_7$Proportion_Lower, rev(Fever_year_proportion_GBD_region_0_7$Proportion_Upper)),
        col = rgb(213/255, 94/255, 0/255, 0.25), border = NA)

lines(Proportion_Median ~ Year, data = Fever_year_proportion_GBD_region_0_6, col = "#E69F00", lwd = 1.5)
polygon(c(Fever_year_proportion_GBD_region_0_6$Year, rev(Fever_year_proportion_GBD_region_0_6$Year)),
        c(Fever_year_proportion_GBD_region_0_6$Proportion_Lower, rev(Fever_year_proportion_GBD_region_0_6$Proportion_Upper)),
        col = rgb(230/255, 159/255, 0/255, 0.25), border = NA)

lines(Proportion_Median ~ Year, data = Fever_year_proportion_GBD_region_0_5, col = "#A7B97F", lwd = 1.5)
polygon(c(Fever_year_proportion_GBD_region_0_5$Year, rev(Fever_year_proportion_GBD_region_0_5$Year)),
        c(Fever_year_proportion_GBD_region_0_5$Proportion_Lower, rev(Fever_year_proportion_GBD_region_0_5$Proportion_Upper)),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

lines(Proportion_Median ~ Year, data = Fever_year_proportion_GBD_region_0_1, col = "#009E73", lwd = 1.5)
polygon(c(Fever_year_proportion_GBD_region_0_1$Year, rev(Fever_year_proportion_GBD_region_0_1$Year)),
        c(Fever_year_proportion_GBD_region_0_1$Proportion_Lower, rev(Fever_year_proportion_GBD_region_0_1$Proportion_Upper)),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

lines(Proportion_Median ~ Year, data = Fever_year_proportion_GBD_region_0_3, col = "#0072B2", lwd = 1.5)
polygon(c(Fever_year_proportion_GBD_region_0_3$Year, rev(Fever_year_proportion_GBD_region_0_3$Year)),
        c(Fever_year_proportion_GBD_region_0_3$Proportion_Lower, rev(Fever_year_proportion_GBD_region_0_3$Proportion_Upper)),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)


mtext("Year", side = 1, line = 1, cex = 0.5)
mtext("H. AFI (malaria non-endemic countries)", at = xmin - 0.3*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.5, font = 1)





plot(Proportion_Median ~ Year, 
     data = Fever_year_proportion_GBD_region_1_2, 
     xlab = "", 
     ylab = "", 
     pch = 16, 
     col = rgb(0.8, 0.8, 0.8, 0), 
     xlim = c(xmin,xmax), 
     ylim = c(ymin,ymax), 
     lwd = "1", 
     axes = F,
     ann = F, 
     main = "")

axis(1, at = seq(xmin, xmax, by = 4), labels = FALSE, las = 2, pos = ymin)
text(x = seq(xmin, xmax, by = 4), y = ymin-(ymax-ymin)*0.065,  
     labels = format(round(seq(xmin, xmax, by = (xmax-xmin)/6*0+4), digits = 0), nsmall = 0), xpd = TRUE, cex = 0.65)

axis(2,  at = seq(ymin, ymax, by = (ymax - ymin)/4), labels = FALSE, las = 2, pos = xmin)
text(y = seq(ymin, ymax, by = (ymax-ymin)/4), x = xmin-0.04*(xmax-xmin),  
     labels = format(round(seq(ymin, ymax, by = (ymax-ymin)/4), digits = 2), nsmall = 1), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)


lines(Proportion_Median ~ Year, data = Fever_year_proportion_GBD_region_1_2, col = "#CC79A7", lwd = 1.5)
polygon(c(Fever_year_proportion_GBD_region_1_2$Year, rev(Fever_year_proportion_GBD_region_1_2$Year)),
        c(Fever_year_proportion_GBD_region_1_2$Proportion_Lower, rev(Fever_year_proportion_GBD_region_1_2$Proportion_Upper)),
        col = rgb(204/255, 121/255, 167/255, 0.25), border = NA)

lines(Proportion_Median ~ Year, data = Fever_year_proportion_GBD_region_1_7, col = "#D55E00", lwd = 1.5)
polygon(c(Fever_year_proportion_GBD_region_1_7$Year, rev(Fever_year_proportion_GBD_region_1_7$Year)),
        c(Fever_year_proportion_GBD_region_1_7$Proportion_Lower, rev(Fever_year_proportion_GBD_region_1_7$Proportion_Upper)),
        col = rgb(213/255, 94/255, 0/255, 0.25), border = NA)

lines(Proportion_Median ~ Year, data = Fever_year_proportion_GBD_region_1_6, col = "#E69F00", lwd = 1.5)
polygon(c(Fever_year_proportion_GBD_region_1_6$Year, rev(Fever_year_proportion_GBD_region_1_6$Year)),
        c(Fever_year_proportion_GBD_region_1_6$Proportion_Lower, rev(Fever_year_proportion_GBD_region_1_6$Proportion_Upper)),
        col = rgb(230/255, 159/255, 0/255, 0.25), border = NA)

lines(Proportion_Median ~ Year, data = Fever_year_proportion_GBD_region_1_5, col = "#A7B97F", lwd = 1.5)
polygon(c(Fever_year_proportion_GBD_region_1_5$Year, rev(Fever_year_proportion_GBD_region_1_5$Year)),
        c(Fever_year_proportion_GBD_region_1_5$Proportion_Lower, rev(Fever_year_proportion_GBD_region_1_5$Proportion_Upper)),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

lines(Proportion_Median ~ Year, data = Fever_year_proportion_GBD_region_1_1, col = "#009E73", lwd = 1.5)
polygon(c(Fever_year_proportion_GBD_region_1_1$Year, rev(Fever_year_proportion_GBD_region_1_1$Year)),
        c(Fever_year_proportion_GBD_region_1_1$Proportion_Lower, rev(Fever_year_proportion_GBD_region_1_1$Proportion_Upper)),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

lines(Proportion_Median ~ Year, data = Fever_year_proportion_GBD_region_1_3, col = "#0072B2", lwd = 1.5)
polygon(c(Fever_year_proportion_GBD_region_1_3$Year, rev(Fever_year_proportion_GBD_region_1_3$Year)),
        c(Fever_year_proportion_GBD_region_1_3$Proportion_Lower, rev(Fever_year_proportion_GBD_region_1_3$Proportion_Upper)),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)


mtext("Year", side = 1, line = 1, cex = 0.5)
mtext("I. AFI (malaria-endemic countries)", at = xmin - 0.3*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.5, font = 1)




#### 10. (legend) ####
par(mar = c(0, 0, 0, 0))
plot.new()

legend_labels <- c("Southeast Asia, East Asia, and Oceania",
                   "Central Europe, Eastern Europe, and Central Asia",
                   "Latin America and Caribbean",
                   "South Asia",
                   "North Africa and Middle East",
                   "Sub-Saharan Africa")
legend_colors <- c("#D55E00",
                   "#CC79A7",
                   "#E69F00",
                   "#A7B97F",
                   "#009E73",
                   "#0072B2")
n <- length(legend_labels)



# Coordinates for the custom legend
x0 <- 0.0  
x1 <- 0.07  
y_start <- 0.7
line_spacing <- 0.07

# Draw each line, label, and shaded polygon
for (i in seq_len(n)) {
  y <- y_start - (i - 1) * line_spacing
  

  polygon(
    x = c(x0-0.01, x1+0.01, x1+0.01, x0-0.01),
    y = c(y - 0.02, y - 0.02, y + 0.02, y + 0.02),
    col = adjustcolor(legend_colors[i], alpha.f = 0.25),
    border = NA
  )
  
 
  segments(x0, y, x1, y, col = legend_colors[i], lwd = 1.5)
  

  text(x = x1 + 0.04, y = y, labels = legend_labels[i], adj = 0, cex = 0.56)
}

dev.off()


######## 10. Figure.3 ########
# for health expenditure and GDP dataset

#### 0. prep ####
rm(list = ls())
gc()
par(mfrow = c(1,1))

requiredPackages = c("dplyr", "tidyr", "data.table", "stringr", "MASS")
lib.path <- .libPaths()
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p, lib = lib.path)
  library(p,character.only = TRUE)
}

n_sim <- 10000
set.seed(123)

#---------------------------prepare data for figure------------------------
UNI_indexed <- readRDS("outputs/sensitivity_analysis/SA_UNI_indexed.rds")
coef_df_LRI <- readRDS("outputs/sensitivity_analysis/SA_coef_df_LRI.rds")
coef_df_Fever <- readRDS("outputs/sensitivity_analysis/SA_coef_df_Fever.rds")
coef_df_LRI_DDD <- readRDS("outputs/sensitivity_analysis/SA_coef_df_LRI_DDD.rds")
coef_df_Fever_DDD <- readRDS("outputs/sensitivity_analysis/SA_coef_df_Fever_DDD.rds")

health_exp <- readRDS("outputs/health_expenditure.rds")
UNI_indexed <- UNI_indexed %>% inner_join(health_exp, by = c("Country", "Year"))

UNI_indexed_Abaccess <- UNI_indexed %>% filter(Abaccess == 1)
UNI_indexed_Abdata <- UNI_indexed %>% filter(ab.data == 1)
UNI_indexed_both <- UNI_indexed %>% filter(Abaccess == 1, ab.data == 1)

expected_coef_names <- c("(Intercept)", "PC1", "PC2")

stopifnot(
  identical(names(coef_df_LRI), expected_coef_names),
  identical(names(coef_df_Fever), expected_coef_names),
  identical(names(coef_df_LRI_DDD), expected_coef_names),
  identical(names(coef_df_Fever_DDD), expected_coef_names)
)

X_pred <- model.matrix(
  ~ PC1 + PC2,
  data = UNI_indexed
)

# Match design matrix columns to each coefficient matrix
X_pred_LRI <- X_pred[, names(coef_df_LRI), drop = FALSE]
X_pred_Fever <- X_pred[, names(coef_df_Fever), drop = FALSE]

X_pred_LRI_DDD <- X_pred[, names(coef_df_LRI_DDD), drop = FALSE]
X_pred_Fever_DDD <- X_pred[, names(coef_df_Fever_DDD), drop = FALSE]


# Linear predictors

LRI_matrix <- as.matrix(X_pred_LRI) %*%
  t(as.matrix(coef_df_LRI))

Fever_matrix <- as.matrix(X_pred_Fever) %*%
  t(as.matrix(coef_df_Fever))

LRI_DDD_matrix <- as.matrix(X_pred_LRI_DDD) %*%
  t(as.matrix(coef_df_LRI_DDD))

Fever_DDD_matrix <- as.matrix(X_pred_Fever_DDD) %*%
  t(as.matrix(coef_df_Fever_DDD))


# Antibiotic-treatment models: logistic link

AbLRI_sim <- plogis(LRI_matrix)
AbFever_sim <- plogis(Fever_matrix)


# DDD models: log-linear models

LRI_DDD_sim <- exp(LRI_DDD_matrix)
Fever_DDD_sim <- exp(Fever_DDD_matrix)

# Flatten the matrices into vectors (column-wise)
AbLRI_long <- as.vector(AbLRI_sim)
LRI_DDD_long <- as.vector(LRI_DDD_sim)
AbFever_long <- as.vector(AbFever_sim)
Fever_DDD_long <- as.vector(Fever_DDD_sim)

# Create matching indices
Year_index <- rep(UNI_indexed$Year, times = ncol(AbLRI_sim))
Country_index <- rep(UNI_indexed$Country, times = ncol(AbLRI_sim))
GDP_index <- rep(UNI_indexed$GDP, times = ncol(AbLRI_sim))
Malaria_index <- rep(UNI_indexed$Endemic, times = ncol(AbLRI_sim))
Health_exp_index <- rep(UNI_indexed$Current_health_expenditure_per_capita_PPP, times = ncol(AbLRI_sim))

# Combine into dataframe
sim_df_LRI <- data.frame(
  Country = Country_index,
  Year = Year_index,
  GDP = GDP_index,
  Malaria = Malaria_index,
  AbLRI = AbLRI_long,
  LRI_DDD = LRI_DDD_long,
  Health_exp = Health_exp_index
)

sim_df_Fever <- data.frame(
  Country = Country_index,
  Year = Year_index,
  GDP = GDP_index,
  Malaria = Malaria_index,
  AbFever = AbFever_long,
  Fever_DDD = Fever_DDD_long,
  Health_exp = Health_exp_index
)

sim_df_LRI_summary <- sim_df_LRI %>%
  group_by(Country, Year) %>%
  summarise(
    GDP = median(GDP, na.rm = TRUE),
    AbLRI_median = median(AbLRI, na.rm = TRUE),
    LRI_DDD_median = median(LRI_DDD, na.rm = TRUE),
    .groups = "drop")

sim_df_Fever_summary <- sim_df_Fever %>%
  group_by(Country, Year) %>%
  summarise(
    GDP = median(GDP, na.rm = TRUE),
    AbFever_median = median(AbFever, na.rm = TRUE),
    Fever_DDD_median = median(Fever_DDD, na.rm = TRUE),
    .groups = "drop")

sampled_df_ordered_LRI <- sim_df_LRI %>% arrange(LRI_DDD) %>%
  mutate(LRI_DDD_rank = cut(LRI_DDD,breaks = quantile(LRI_DDD, probs = seq(0, 1, 0.1), na.rm = TRUE),include.lowest = TRUE))

sampled_df_ordered_Fever <- sim_df_Fever %>% arrange(Fever_DDD) %>%
  mutate(Fever_DDD_rank = cut(Fever_DDD,breaks = quantile(Fever_DDD, probs = seq(0, 1, 0.1), na.rm = TRUE),include.lowest = TRUE))


group_df_LRI <- sampled_df_ordered_LRI %>%
  group_by(LRI_DDD_rank) %>%
  summarise(
    med_AbLRI   = median(AbLRI, na.rm = TRUE),
    lower_AbLRI = quantile(AbLRI, 0.025, na.rm = TRUE),
    upper_AbLRI = quantile(AbLRI, 0.975, na.rm = TRUE),
    seq_DDD     = median(LRI_DDD, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  na.omit()


group_df_Fever <- sampled_df_ordered_Fever %>%
  group_by(Fever_DDD_rank) %>%
  summarise(
    med_AbFever   = median(AbFever, na.rm = TRUE),
    lower_AbFever = quantile(AbFever, 0.025, na.rm = TRUE),
    upper_AbFever = quantile(AbFever, 0.975, na.rm = TRUE),
    seq_DDD       = median(Fever_DDD, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  na.omit()





# ------------------------1. P(Antibiotic treatment|infection) ~ Antibiotic consumption/incidence-----------------
sim_df_LRI_every_1000 <- sim_df_LRI[seq(1, nrow(sim_df_LRI), by = 10000), ]
year_range <- range(sim_df_LRI_every_1000$Year)
sim_df_Fever_every_1000 <- sim_df_Fever[seq(1, nrow(sim_df_Fever), by = 10000), ]

layout(matrix(c(1,2), nrow = 1, byrow = TRUE))
xmin <- min(group_df_LRI$seq_DDD); xmax <- max(group_df_LRI$seq_DDD)

ymin <- 0 ; ymax <- 1

par(oma = c(0, 0, 0, 0))
par(mar = c(2.5, 2.5, 2, 1))  
par(mgp = c(3, 0.3, 0))     
par(tck = -0.025)
par(lwd = 1)


plot(med_AbLRI ~ seq_DDD, 
     data = group_df_LRI,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.75, 0.75, 0.75, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")



axis(1, at = seq(xmin, xmax, by = (xmax - xmin)/4), labels = FALSE, las = 2, pos = ymin)
text(x = seq(xmin, xmax, by = (xmax - xmin)/4), y = ymin-0.08*(ymax-ymin),  
     labels = format(round(seq(xmin, xmax, by = ((xmax - xmin)/4)), digits = 0), nsmall = 0), xpd = TRUE, cex = 0.65)

axis(2,  at = seq(ymin, ymax, by = (ymax - ymin)/4), labels = FALSE, las = 2, pos = xmin)
text(y = seq(ymin, ymax, by = (ymax-ymin)/4), x = xmin-0.04*(xmax-xmin),  
     labels = format(round(seq(ymin, ymax, by = (ymax-ymin)/4), digits = 2), nsmall = 1), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)



lines(med_AbLRI ~ seq_DDD, data = group_df_LRI, col = "#0072B2", lwd = 1.5)

# Add shaded area for confidence interval
polygon(c(group_df_LRI$seq_DDD, rev(group_df_LRI$seq_DDD)),
        c(group_df_LRI$lower_AbLRI, rev(group_df_LRI$upper_AbLRI)),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)

# Probability of antibotics treatment by antibiotic consumption (DDD) per incidence
mtext("A. Lower Respiratory Infection", at = xmin - 0.2*(xmax-xmin), side = 3, line = 1, adj = 0, cex = 0.65, font = 2)
mtext("Antibiotic consumption (DDD) per incidence", side = 1, line = 1.5, cex = 0.65)
mtext(expression("Probability of antibiotic treatment"), side = 2, line = 1.5, cex = 0.65)  


xmin <-  min(group_df_Fever$seq_DDD); xmax <- max(group_df_Fever$seq_DDD)



plot(med_AbFever ~ seq_DDD, 
     data = group_df_Fever,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.75, 0.75, 0.75, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")



axis(1, at = seq(xmin, xmax, by = (xmax - xmin)/4), labels = FALSE, las = 2, pos = ymin)
text(x = seq(xmin, xmax, by = (xmax - xmin)/4), y = ymin-0.08*(ymax-ymin),  
     labels = format(round(seq(xmin, xmax, by = ((xmax - xmin)/4)), digits = 1), nsmall = 0), xpd = TRUE, cex = 0.65)

axis(2, at = seq(ymin, ymax, by = (ymax - ymin)/4), labels = FALSE, las = 2, pos = xmin)
text(y = seq(ymin, ymax, by = (ymax-ymin)/4), x = xmin-0.04*(xmax-xmin),  
     labels = format(round(seq(ymin, ymax, by = (ymax-ymin)/4), digits = 2), nsmall = 1), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines(med_AbFever ~ seq_DDD, data = group_df_Fever, col = "#0072B2", lwd = 1.5)

# Add shaded area for confidence interval
polygon(c(group_df_Fever$seq_DDD, rev(group_df_Fever$seq_DDD)),
        c(group_df_Fever$lower_AbFever, rev(group_df_Fever$upper_AbFever)),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)

mtext("B. Febrile Illness without Lower Respiratory Infection", at = xmin - 0.2*(xmax-xmin), side = 3, line = 1, adj = 0, cex = 0.65, font = 2)
mtext("Antibiotic consumption (DDD) per incidence", side = 1, line = 1.5, cex = 0.65)














# ------------------------2.Antibiotic consumption/incidence ~ GDP --------------------
LRI_DDD_sim_long <- as.vector(LRI_DDD_sim)
LRI_DDD_gdp <- as.data.frame(UNI_indexed$GDP)

for (i in 1:nrow(LRI_DDD_gdp)){
  LRI_DDD_gdp[i,2] <- quantile(LRI_DDD_sim[i,1:n_sim], probs = 0.5)
  LRI_DDD_gdp[i,3] <- quantile(LRI_DDD_sim[i,1:n_sim], probs = 0.025)
  LRI_DDD_gdp[i,4] <- quantile(LRI_DDD_sim[i,1:n_sim], probs = 0.975)
}

colnames(LRI_DDD_gdp) = c("GDP","med","lower","upper")



# Create and process the merged dataset
LRI_DDD_gdp_long <- rep(LRI_DDD_gdp$GDP, n_sim)
LRI_DDD_gdp_merged <- data.frame(LRI_DDD = LRI_DDD_sim_long, GDP = LRI_DDD_gdp_long) %>%
  arrange(GDP) %>% mutate(GDP_rank = cut(GDP,breaks = quantile(GDP, probs = seq(0, 1, 0.1), na.rm = TRUE),include.lowest = TRUE))

LRI_DDD_gdp_plot <- LRI_DDD_gdp_merged %>% group_by(GDP_rank) %>%
  summarize(
    med_DDD = median(LRI_DDD),
    lower_DDD = quantile(LRI_DDD, 0.025),
    upper_DDD = quantile(LRI_DDD, 0.975),
    GDP = mean(GDP),
    .groups = "drop")


Fever_DDD_gdp <- as.data.frame(UNI_indexed$GDP)
for (i in 1:nrow(Fever_DDD_gdp)){
  Fever_DDD_gdp[i,2] <- quantile(Fever_DDD_sim[i,1:n_sim], probs = 0.5)
  Fever_DDD_gdp[i,3] <- quantile(Fever_DDD_sim[i,1:n_sim], probs = 0.025)
  Fever_DDD_gdp[i,4] <- quantile(Fever_DDD_sim[i,1:n_sim], probs = 0.975)
}

colnames(Fever_DDD_gdp) = c("GDP","med","lower","upper")

Fever_DDD_gdp <- Fever_DDD_gdp %>% arrange(Fever_DDD_gdp$GDP)

Fever_DDD_sim_long <- as.vector(Fever_DDD_sim)
Fever_DDD_gdp_long <- rep(Fever_DDD_gdp$GDP, n_sim)

# Create and process the merged dataset
Fever_DDD_gdp_merged <- data.frame(Fever_DDD = Fever_DDD_sim_long, GDP = Fever_DDD_gdp_long) %>%
  arrange(GDP) %>% mutate(GDP_rank = cut(GDP,breaks = quantile(GDP, probs = seq(0, 1, 0.1), na.rm = TRUE),include.lowest = TRUE))

Fever_DDD_gdp_plot <- Fever_DDD_gdp_merged %>% group_by(GDP_rank) %>%
  summarize(
    med_DDD = median(Fever_DDD),
    lower_DDD = quantile(Fever_DDD, 0.025),
    upper_DDD = quantile(Fever_DDD, 0.975),
    GDP = mean(GDP),
    .groups = "drop")








#-----------------------P(Antibiotic treatment|infection) ~ GDP ----------------------
# 3.1 AbLRI ~ GDP -------------------------
AbLRI_gdp <- as.data.frame(UNI_indexed$GDP)
AbLRI_health_exp <- as.data.frame(UNI_indexed$Current_health_expenditure_per_capita_PPP)

# LRI_DDD_gdp <- LRI_DDD_gdp %>% arrange(LRI_DDD_gdp$GDP)
AbLRI_sim_long <- as.vector(AbLRI_sim)
AbLRI_gdp_long <- as.numeric(rep(AbLRI_gdp$'UNI_indexed$GDP', n_sim))
AbLRI_health_exp_long <- as.numeric(rep(AbLRI_health_exp$'UNI_indexed$Current_health_expenditure_per_capita_PPP', n_sim))

# Create and process the merged dataset
AbLRI_gdp_merged <- data.frame(AbLRI = AbLRI_sim_long, GDP = AbLRI_gdp_long) %>% 
  arrange(GDP) %>%
  mutate(GDP_rank = cut(GDP,breaks = quantile(GDP, probs = seq(0, 1, 0.1), na.rm = TRUE),include.lowest = TRUE))

AbLRI_health_exp_merged <- data.frame(AbLRI = AbLRI_sim_long, health_exp = AbLRI_health_exp_long) %>% 
  arrange(health_exp) %>%
  mutate(health_exp_rank = cut(health_exp, breaks = quantile(health_exp, probs = seq(0, 1, 0.1), na.rm = TRUE),include.lowest = TRUE)) 

AbLRI_gdp_plot <- AbLRI_gdp_merged %>% group_by(GDP_rank) %>% 
  summarize(
    med_AbLRI = median(AbLRI),
    lower_AbLRI = quantile(AbLRI, 0.025),
    upper_AbLRI = quantile(AbLRI, 0.975),
    GDP = mean(GDP),
    .groups = "drop")

AbLRI_health_exp_plot <- AbLRI_health_exp_merged %>% group_by(health_exp_rank) %>% 
  summarize(
    med_AbLRI = median(AbLRI),
    lower_AbLRI = quantile(AbLRI, 0.025),
    upper_AbLRI = quantile(AbLRI, 0.975),
    health_exp = mean(health_exp),
    .groups = "drop")



AbFever_gdp <- as.data.frame(UNI_indexed$GDP)
AbFever_health_exp <- as.data.frame(UNI_indexed$Current_health_expenditure_per_capita_PPP)

AbFever_sim_long <- as.vector(AbFever_sim)
AbFever_gdp_long <- as.numeric(rep(AbFever_gdp$'UNI_indexed$GDP', n_sim))
AbFever_health_exp_long <- as.numeric(rep(AbFever_health_exp$'UNI_indexed$Current_health_expenditure_per_capita_PPP', n_sim))

# Create and process the merged dataset
AbFever_gdp_merged <- data.frame(AbFever = AbFever_sim_long, GDP = AbFever_gdp_long) %>% 
  arrange(GDP) %>%
  mutate(GDP_rank = cut(GDP,breaks = quantile(GDP, probs = seq(0, 1, 0.1), na.rm = TRUE),include.lowest = TRUE))

AbFever_health_exp_merged <- data.frame(AbFever = AbFever_sim_long, health_exp = AbFever_health_exp_long) %>% 
  arrange(health_exp) %>%
  mutate(health_exp_rank = cut(health_exp, breaks = quantile(health_exp, probs = seq(0, 1, 0.1), na.rm = TRUE),include.lowest = TRUE)) 

AbFever_gdp_plot <- AbFever_gdp_merged %>% group_by(GDP_rank) %>% 
  summarize(
    med_AbFever = median(AbFever),
    lower_AbFever = quantile(AbFever, 0.025),
    upper_AbFever = quantile(AbFever, 0.975),
    GDP = mean(GDP),
    .groups = "drop")

AbFever_health_exp_plot <- AbFever_health_exp_merged %>% group_by(health_exp_rank) %>% 
  summarize(
    med_AbFever = median(AbFever),
    lower_AbFever = quantile(AbFever, 0.025),
    upper_AbFever = quantile(AbFever, 0.975),
    health_exp = mean(health_exp),
    .groups = "drop")

saveRDS(AbLRI_gdp_plot, "outputs/global_map_df/economic/SA_LRI_GDP.rds")
saveRDS(AbFever_gdp_plot, "outputs/global_map_df/economic/SA_Fever_GDP.rds")
saveRDS(AbLRI_health_exp_plot, "outputs/global_map_df/economic/SA_LRI_health_exp.rds")
saveRDS(AbFever_health_exp_plot, "outputs/global_map_df/economic/SA_Fever_health_exp.rds")



tiff("outputs/plot_df/SA_2_Probability_Consumption~GDP.tiff", width = 7.5, height = 3.5, units = "in", res = 300)
par(mfrow = c(1,2))

par(oma = c(0.5, 5, 0.5, 1))
par(mar = c(2.5, 1.5, 2, 1))  
par(mgp = c(3, 0.3, 0))     
par(tck = -0.025)
par(lwd = 1)

xmin <- min(AbLRI_gdp_plot$GDP)         ;     xmax <- max(UNI_indexed_Abaccess$GDP)       
ymin <- min(AbLRI_gdp_plot$lower_AbLRI) * 0 + 0 ; ymax <- max(AbLRI_gdp_plot$upper_AbLRI) * 0 + 1


plot(med_AbLRI ~ GDP, 
     data = AbLRI_gdp_plot,
     xlab = "Average AbLRI (LRI)", 
     ylab = "Median AbLRI",
     pch = 16, col = rgb(0.75, 0.75, 0.75, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")



axis(1, at = seq(xmin, xmax, length.out = 8), labels = FALSE, las = 2, pos = ymin)
text(x = seq(xmin, xmax, by = (xmax - xmin)/7), y = ymin-0.08*(ymax-ymin),  
     labels = format(round(seq(xmin, xmax, by = ((xmax - xmin)/7)), digits = 0), nsmall = 0), xpd = TRUE, cex = 0.65, srt = 30)

axis(2,  at = seq(ymin, ymax, length.out = 11), labels = FALSE, las = 2, pos = xmin)
text(y = seq(ymin, ymax, by = (ymax-ymin)*0.1), x = xmin-0.04*(xmax-xmin),  
     labels = format(round(seq(ymin, ymax, by = (ymax-ymin)*0.1), digits = 2), nsmall = 1), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines(med_AbLRI ~ GDP, data = AbLRI_gdp_plot, col = "#0072B2", lwd = 1.5)

# Add shaded area for confidence interval
polygon(c(AbLRI_gdp_plot$GDP, rev(AbLRI_gdp_plot$GDP)),
        c(AbLRI_gdp_plot$lower_AbLRI, rev(AbLRI_gdp_plot$upper_AbLRI)),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)

mtext("Probability of antibotics treatment by GDP per capita PPP", side = 3, line = 1, adj = 0, at = xmin - (xmax - xmin) * 0.53, cex = 0.65, font = 2)
text(x = mean(range(xmin, xmax)), y = ymax, labels = "LRI", cex = 0.65, font = 2)
mtext(expression("       Probability of \n antibiotic treatment"), side = 2, line = 2, cex = 0.65)  
mtext(expression("GDP per capita PPP"), side = 1, line = 1.5, cex = 0.65)  




xmin <- min(AbFever_gdp_plot$GDP)         ;     xmax <- max(UNI_indexed_Abaccess$GDP)       
ymin <- min(AbFever_gdp_plot$lower_AbFever) * 0 + 0; ymax <- max(AbFever_gdp_plot$upper_AbFever) * 0 + 1


plot(med_AbFever ~ GDP, 
     data = AbFever_gdp_plot,
     xlab = "Average AbFever (Fever)", 
     ylab = "Median AbFever",
     pch = 16, col = rgb(0.75, 0.75, 0.75, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")



axis(1, at = seq(xmin, xmax, length.out = 8), labels = FALSE, las = 2, pos = ymin)
text(x = seq(xmin, xmax, by = (xmax - xmin)/7), y = ymin-0.08*(ymax-ymin),  
     labels = format(round(seq(xmin, xmax, by = ((xmax - xmin)/7)), digits = 0), nsmall = 0), xpd = TRUE, cex = 0.65, srt = 30)

axis(2,  at = seq(ymin, ymax, length.out = 11), labels = FALSE, las = 2, pos = xmin)
text(y = seq(ymin, ymax, by = (ymax-ymin)*0.1), x = xmin-0.04*(xmax-xmin),  
     labels = format(round(seq(ymin, ymax, by = (ymax-ymin)*0.1), digits = 2), nsmall = 1), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines(med_AbFever ~ GDP, data = AbFever_gdp_plot, col = "#CC79A7", lwd = 1.5)

# Add shaded area for confidence interval
polygon(c(AbFever_gdp_plot$GDP, rev(AbFever_gdp_plot$GDP)),
        c(AbFever_gdp_plot$lower_AbFever, rev(AbFever_gdp_plot$upper_AbFever)),
        col = rgb(204/255, 121/255, 167/255, 0.25), border = NA)

text(x = mean(range(xmin, xmax)), y = ymax, labels = "Fever without LRI", cex = 0.65, font = 2)

mtext(expression("GDP per capita PPP"), side = 1, line = 1.5, cex = 0.65)  


dev.off()













tiff("outputs/plot_df/SA_S1. Coefficient.tiff", width = 10, height = 8, units = "in", res = 300)

par(mar = c(4.5, 4.5, 3, 1))  
par(mgp = c(3, 0.3, 0)) 
par(tck = -0.025)       
par(lwd = 1)            
par(mfrow = c(1,2))

plot(coef_df_LRI_DDD$PC1, coef_df_LRI_DDD$PC2, 
     xlim = c(0.03, 0.16), ylim = c(-0.01, 0.09),
     xlab = "Antibiotic Consumption (PC1)", ylab = "Antibiotic Consumption (PC2)",
     col = "#0072B2")
points(coef_df_Fever_DDD$PC1, coef_df_Fever_DDD$PC2, col = "#CC79A7")

legend("topright", legend = c(expression("LRI"), expression("Fever")), 
       col = c("#0072B2", "#CC79A7"), pch = 21)

plot(coef_df_LRI$PC1, coef_df_LRI$PC2, 
     xlim = c(0, 0.06), ylim = c(-0.03, 0.1),
     xlab = "Probability of antibiotic treatment (PC1)", ylab = "Probability of antibiotic treatment (PC2)",
     col = "#0072B2")
points(coef_df_Fever$PC1, coef_df_Fever$PC2, col = "#CC79A7")

legend("topright", legend = c(expression("LRI"), expression("Fever")), 
       col = c("#0072B2", "#CC79A7"), pch = 21)

par(mfrow = c(1,1))

dev.off()


#### 0. prep ####
rm(list = ls())
gc()

par(mfrow = c(1,1))

requiredPackages = c("dplyr", "tidyr", "data.table", "stringr")
lib.path <- .libPaths()
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p, lib = lib.path)
  library(p,character.only = TRUE)
}

# setwd() removed for GitHub repository compatibility

n_sim <- 10000
set.seed(123)

#---------------------------prepare data for figure------------------------
UNI_indexed     <- readRDS("outputs/sensitivity_analysis/SA_UNI_indexed.rds")
coef_df_LRI     <- readRDS("outputs/sensitivity_analysis/SA_coef_df_LRI.rds")
coef_df_Fever   <- readRDS("outputs/sensitivity_analysis/SA_coef_df_Fever.rds")
coef_df_LRI_DDD <- readRDS("outputs/sensitivity_analysis/SA_coef_df_LRI_DDD.rds")
coef_df_Fever_DDD <- readRDS("outputs/sensitivity_analysis/SA_coef_df_Fever_DDD.rds")

UNI_indexed_Abaccess <- UNI_indexed %>% filter(Abaccess == 1)
UNI_indexed_Abdata <- UNI_indexed %>% filter(ab.data == 1)
UNI_indexed_both <- UNI_indexed %>% filter(Abaccess == 1, ab.data == 1)

UNI_indexed_0 <- UNI_indexed[UNI_indexed$Endemic == 0,]
UNI_indexed_1 <- UNI_indexed[UNI_indexed$Endemic == 1,]


# all current models: PC1-PC2
stopifnot(
  identical(names(coef_df_LRI), c("(Intercept)", "PC1", "PC2")),
  identical(names(coef_df_Fever), c("(Intercept)", "PC1", "PC2")),
  identical(names(coef_df_LRI_DDD), c("(Intercept)", "PC1", "PC2")),
  identical(names(coef_df_Fever_DDD), c("(Intercept)", "PC1", "PC2"))
)

X_pred <- model.matrix(~ PC1 + PC2, data = UNI_indexed)
X_pred_0 <- model.matrix(~ PC1 + PC2, data = UNI_indexed_0)
X_pred_1 <- model.matrix(~ PC1 + PC2, data = UNI_indexed_1)

X_pred_LRI <- X_pred[, names(coef_df_LRI), drop = FALSE]
X_pred_Fever_0 <- X_pred_0[, names(coef_df_Fever), drop = FALSE]
X_pred_Fever_1 <- X_pred_1[, names(coef_df_Fever), drop = FALSE]

X_pred_LRI_DDD <- X_pred[, names(coef_df_LRI_DDD), drop = FALSE]
X_pred_Fever_DDD_0 <- X_pred_0[, names(coef_df_Fever_DDD), drop = FALSE]
X_pred_Fever_DDD_1 <- X_pred_1[, names(coef_df_Fever_DDD), drop = FALSE]

LRI_matrix <- as.matrix(X_pred_LRI) %*% t(as.matrix(coef_df_LRI))
Fever_matrix_0 <- as.matrix(X_pred_Fever_0) %*% t(as.matrix(coef_df_Fever))
Fever_matrix_1 <- as.matrix(X_pred_Fever_1) %*% t(as.matrix(coef_df_Fever))

LRI_DDD_matrix <- as.matrix(X_pred_LRI_DDD) %*% t(as.matrix(coef_df_LRI_DDD))
Fever_DDD_matrix_0 <- as.matrix(X_pred_Fever_DDD_0) %*% t(as.matrix(coef_df_Fever_DDD))
Fever_DDD_matrix_1 <- as.matrix(X_pred_Fever_DDD_1) %*% t(as.matrix(coef_df_Fever_DDD))

# convert to antibiotic treatment probability using logit link
AbLRI_sim <- plogis(LRI_matrix)
AbFever_sim_0 <- plogis(Fever_matrix_0)
AbFever_sim_1 <- plogis(Fever_matrix_1)

# convert to expected DDD using exp link function
LRI_DDD_sim <- exp(LRI_DDD_matrix)
Fever_DDD_sim_0 <- exp(Fever_DDD_matrix_0)
Fever_DDD_sim_1 <- exp(Fever_DDD_matrix_1)

# Flatten the matrices into vectors (column-wise)
AbLRI_long <- as.vector(AbLRI_sim)
LRI_DDD_long <- as.vector(LRI_DDD_sim)
AbFever_long_0 <- as.vector(AbFever_sim_0)
AbFever_long_1 <- as.vector(AbFever_sim_1)
Fever_DDD_long_0 <- as.vector(Fever_DDD_sim_0)
Fever_DDD_long_1 <- as.vector(Fever_DDD_sim_1)

# Create matching indices
Year_index <- rep(UNI_indexed$Year, times = ncol(AbLRI_sim))
Country_index <- rep(UNI_indexed$Country, times = ncol(AbLRI_sim))
GDP_index <- rep(UNI_indexed$GDP, times = ncol(AbLRI_sim))
Malaria_index <- rep(UNI_indexed$Endemic, times = ncol(AbLRI_sim))

Year_index_0 <- rep(UNI_indexed_0$Year, times = ncol(AbFever_sim_0))
Country_index_0 <- rep(UNI_indexed_0$Country, times = ncol(AbFever_sim_0))
GDP_index_0 <- rep(UNI_indexed_0$GDP, times = ncol(AbFever_sim_0))
Malaria_index_0 <- rep(UNI_indexed_0$Endemic, times = ncol(AbFever_sim_0))

Year_index_1 <- rep(UNI_indexed_1$Year, times = ncol(AbFever_sim_1))
Country_index_1 <- rep(UNI_indexed_1$Country, times = ncol(AbFever_sim_1))
GDP_index_1 <- rep(UNI_indexed_1$GDP, times = ncol(AbFever_sim_1))
Malaria_index_1 <- rep(UNI_indexed_1$Endemic, times = ncol(AbFever_sim_1))

# Combine into dataframe
sim_df_LRI <- data.frame(
  Country = Country_index,
  Year = Year_index,
  GDP = GDP_index,
  Malaria = Malaria_index,
  AbLRI = AbLRI_long,
  lack_AbLRI = 1-AbLRI_long,
  LRI_DDD = LRI_DDD_long
)

sim_df_Fever_0 <- data.frame(
  Country = Country_index_0,
  Year = Year_index_0,
  GDP = GDP_index_0,
  Malaria = Malaria_index_0,
  AbFever = AbFever_long_0,
  lack_AbFever = 1-AbFever_long_0,
  Fever_DDD = Fever_DDD_long_0
)

sim_df_Fever_1 <- data.frame(
  Country = Country_index_1,
  Year = Year_index_1,
  GDP = GDP_index_1,
  Malaria = Malaria_index_1,
  AbFever = AbFever_long_1,
  lack_AbFever = 1-AbFever_long_1,
  Fever_DDD = Fever_DDD_long_1
)

sim_df_LRI_summary <- sim_df_LRI %>%
  group_by(Country, Year) %>%
  summarise(
    GDP = median(GDP, na.rm = TRUE),
    AbLRI_median = median(AbLRI, na.rm = TRUE),
    AbLRI_median_lack = median(lack_AbLRI, na.rm = TRUE),
    LRI_DDD_median = median(LRI_DDD, na.rm = TRUE),
    .groups = "drop")

sim_df_Fever_summary_0 <- sim_df_Fever_0 %>%
  group_by(Country, Year) %>%
  summarise(
    GDP = median(GDP, na.rm = TRUE),
    AbFever_median = median(AbFever, na.rm = TRUE),
    AbFever_median_lack = median(lack_AbFever, na.rm = TRUE),
    Fever_DDD_median = median(Fever_DDD, na.rm = TRUE),
    .groups = "drop")

sim_df_Fever_summary_1 <- sim_df_Fever_1 %>%
  group_by(Country, Year) %>%
  summarise(
    GDP = median(GDP, na.rm = TRUE),
    AbFever_median = median(AbFever, na.rm = TRUE),
    AbFever_median_lack = median(lack_AbFever, na.rm = TRUE),
    Fever_DDD_median = median(Fever_DDD, na.rm = TRUE),
    .groups = "drop")

sampled_df_ordered_LRI <- sim_df_LRI %>% arrange(LRI_DDD) %>%
  mutate(LRI_DDD_rank = cut(LRI_DDD, breaks = quantile(LRI_DDD, probs = seq(0, 1, 0.1), na.rm = TRUE), include.lowest = TRUE))

sampled_df_ordered_Fever_0 <- sim_df_Fever_0 %>% arrange(Fever_DDD) %>%
  mutate(Fever_DDD_rank = cut(Fever_DDD, breaks = quantile(Fever_DDD, probs = seq(0, 1, 0.1), na.rm = TRUE), include.lowest = TRUE))

sampled_df_ordered_Fever_1 <- sim_df_Fever_1 %>% arrange(Fever_DDD) %>%
  mutate(Fever_DDD_rank = cut(Fever_DDD, breaks = quantile(Fever_DDD, probs = seq(0, 1, 0.1), na.rm = TRUE), include.lowest = TRUE))

group_df_LRI <- sampled_df_ordered_LRI %>% group_by(LRI_DDD_rank) %>% 
  summarize(
    med_AbLRI = median(AbLRI),
    lower_AbLRI = quantile(AbLRI, 0.025),
    upper_AbLRI = quantile(AbLRI, 0.975),
    lack_med_AbLRI = median(lack_AbLRI),
    lack_lower_AbLRI = quantile(lack_AbLRI, 0.025),
    lack_upper_AbLRI = quantile(lack_AbLRI, 0.975),
    seq_DDD = median(LRI_DDD),
    .groups = "drop") %>%  
  distinct() %>% na.omit()

group_df_Fever_0 <- sampled_df_ordered_Fever_0 %>% group_by(Fever_DDD_rank) %>% 
  summarize(
    med_AbFever = median(AbFever),
    lower_AbFever = quantile(AbFever, 0.025),
    upper_AbFever = quantile(AbFever, 0.975),
    lack_med_AbFever = median(lack_AbFever),
    lack_lower_AbFever = quantile(lack_AbFever, 0.025),
    lack_upper_AbFever = quantile(lack_AbFever, 0.975),
    seq_DDD = median(Fever_DDD),
    .groups = "drop") %>%  
  distinct() %>% na.omit()

group_df_Fever_1 <- sampled_df_ordered_Fever_1 %>% group_by(Fever_DDD_rank) %>% 
  summarize(
    med_AbFever = median(AbFever),
    lower_AbFever = quantile(AbFever, 0.025),
    upper_AbFever = quantile(AbFever, 0.975),
    lack_med_AbFever = median(lack_AbFever),
    lack_lower_AbFever = quantile(lack_AbFever, 0.025),
    lack_upper_AbFever = quantile(lack_AbFever, 0.975),
    seq_DDD = median(Fever_DDD),
    .groups = "drop") %>%  
  distinct() %>% na.omit()

group_df_LRI$ideal <- NA
for (i in 1:10){
  intercept <- group_df_LRI$med_AbLRI[1]
  slope <- (group_df_LRI$med_AbLRI[2] - group_df_LRI$med_AbLRI[1]) / (group_df_LRI$seq_DDD[2] - group_df_LRI$seq_DDD[1])
  group_df_LRI$ideal[i] <- intercept + slope * (group_df_LRI$seq_DDD[i] - group_df_LRI$seq_DDD[1])
}
group_df_LRI$ideal <- 1-group_df_LRI$ideal

group_df_Fever_0$ideal <- NA
for (i in 1:10){
  intercept <- group_df_Fever_0$med_AbFever[1]
  slope <- (group_df_Fever_0$med_AbFever[2] - group_df_Fever_0$med_AbFever[1]) / (group_df_Fever_0$seq_DDD[2] - group_df_Fever_0$seq_DDD[1])
  group_df_Fever_0$ideal[i] <- intercept + slope * (group_df_Fever_0$seq_DDD[i] - group_df_Fever_0$seq_DDD[1])
}
group_df_Fever_0$ideal <- 1-group_df_Fever_0$ideal

group_df_Fever_1$ideal <- NA
for (i in 1:10){
  intercept <- group_df_Fever_1$med_AbFever[1]
  slope <- (group_df_Fever_1$med_AbFever[2] - group_df_Fever_1$med_AbFever[1]) / (group_df_Fever_1$seq_DDD[2] - group_df_Fever_1$seq_DDD[1])
  group_df_Fever_1$ideal[i] <- intercept + slope * (group_df_Fever_1$seq_DDD[i] - group_df_Fever_1$seq_DDD[1])
}
group_df_Fever_1$ideal <- 1-group_df_Fever_1$ideal

LRI_gdp_plot   <- readRDS("outputs/global_map_df/economic/SA_LRI_GDP.rds")
Fever_gdp_plot <- readRDS("outputs/global_map_df/economic/SA_Fever_GDP.rds")

UNI_indexed <- readRDS("outputs/sensitivity_analysis/SA_UNI_indexed.rds")
health_exp <- readRDS("outputs/health_expenditure.rds")
UNI_indexed <- UNI_indexed %>% inner_join(health_exp, by = c("Country", "Year"))

UNI_indexed_Abaccess <- UNI_indexed %>% filter(Abaccess == 1)
UNI_indexed_Abdata <- UNI_indexed %>% filter(ab.data == 1)
UNI_indexed_both <- UNI_indexed %>% filter(Abaccess == 1, ab.data == 1)

LRI_health_exp_plot   <- readRDS("outputs/global_map_df/economic/SA_LRI_health_exp.rds")
Fever_health_exp_plot <- readRDS("outputs/global_map_df/economic/SA_Fever_health_exp.rds")

intercept_LRI_ideal <- group_df_LRI$med_AbLRI[1]
slope_LRI_ideal <- (group_df_LRI$med_AbLRI[2] - group_df_LRI$med_AbLRI[1]) / (group_df_LRI$seq_DDD[2] - group_df_LRI$seq_DDD[1])
xmax_LRI <- (1 - intercept_LRI_ideal)/slope_LRI_ideal + group_df_LRI$seq_DDD[1]

intercept_LRI_obs <- group_df_LRI$med_AbLRI[8]
slope_LRI_obs <- (group_df_LRI$med_AbLRI[9] - group_df_LRI$med_AbLRI[8]) / (group_df_LRI$seq_DDD[9] - group_df_LRI$seq_DDD[8])
med_AbLRI_at_xmax <- intercept_LRI_obs + slope_LRI_obs * (xmax_LRI - group_df_LRI$seq_DDD[8])

group_df_LRI_2 <- group_df_LRI[1:8,]
group_df_LRI_2[9,] <- NA
group_df_LRI_2$seq_DDD[9] <- xmax_LRI
group_df_LRI_2$ideal[9] <- 0
group_df_LRI_2$lack_med_AbLRI[9] <- 1-med_AbLRI_at_xmax

intercept_Fever_0_ideal <- group_df_Fever_0$med_AbFever[1]
slope_Fever_0_ideal <- (group_df_Fever_0$med_AbFever[2] - group_df_Fever_0$med_AbFever[1]) / (group_df_Fever_0$seq_DDD[2] - group_df_Fever_0$seq_DDD[1])
xmax_Fever_0 <- (1 - intercept_Fever_0_ideal)/slope_Fever_0_ideal + group_df_Fever_0$seq_DDD[1]

intercept_Fever_0_obs <- group_df_Fever_0$med_AbFever[9]
slope_Fever_0_obs <- (group_df_Fever_0$med_AbFever[10] - group_df_Fever_0$med_AbFever[9]) / (group_df_Fever_0$seq_DDD[10] - group_df_Fever_0$seq_DDD[9])
med_AbFever_at_xmax_0 <- intercept_Fever_0_obs + slope_Fever_0_obs * (xmax_Fever_0 - group_df_Fever_0$seq_DDD[9])

group_df_Fever_0_2 <- group_df_Fever_0
group_df_Fever_1_2 <- group_df_Fever_1

write.csv(group_df_LRI, "outputs/plot_df/SA_Figure3_LRI.csv")
write.csv(group_df_Fever_0, "outputs/plot_df/SA_Figure3_Fever0.csv")
write.csv(group_df_Fever_1, "outputs/plot_df/SA_Figure3_Fever1.csv")
write.csv(group_df_LRI_2, "outputs/plot_df/SA_Figure3_LRI_2.csv")
write.csv(group_df_Fever_0_2, "outputs/plot_df/SA_Figure3_Fever0_2.csv")
write.csv(group_df_Fever_1_2, "outputs/plot_df/SA_Figure3_Fever1_2.csv")




# ============================================================
# Export source data for sensitivity-analysis Figure 3
# ============================================================

dir.create(
  "outputs/plot_objects",
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# Panel A: LRTI by GDP
# ------------------------------------------------------------

fig3_A <- LRI_gdp_plot %>%
  transmute(
    Panel = "A",
    Syndrome = "LRTI",
    Malaria_status = "Not applicable",
    X_variable = "GDP per capita (PPP)",
    X_value = GDP,
    Median = med_AbLRI,
    Lower = lower_AbLRI,
    Upper = upper_AbLRI,
    Projected = NA_real_
  )


# ------------------------------------------------------------
# Panel B: AFI by GDP
# ------------------------------------------------------------

fig3_B <- Fever_gdp_plot %>%
  transmute(
    Panel = "B",
    Syndrome = "AFI",
    Malaria_status = "Not stratified",
    X_variable = "GDP per capita (PPP)",
    X_value = GDP,
    Median = med_AbFever,
    Lower = lower_AbFever,
    Upper = upper_AbFever,
    Projected = NA_real_
  )


# ------------------------------------------------------------
# Panel C: LRTI by health expenditure
# ------------------------------------------------------------

fig3_C <- LRI_health_exp_plot %>%
  transmute(
    Panel = "C",
    Syndrome = "LRTI",
    Malaria_status = "Not applicable",
    X_variable = "Health expenditure per capita (PPP)",
    X_value = health_exp,
    Median = med_AbLRI,
    Lower = lower_AbLRI,
    Upper = upper_AbLRI,
    Projected = NA_real_
  )


# ------------------------------------------------------------
# Panel D: AFI by health expenditure
# ------------------------------------------------------------

fig3_D <- Fever_health_exp_plot %>%
  transmute(
    Panel = "D",
    Syndrome = "AFI",
    Malaria_status = "Not stratified",
    X_variable = "Health expenditure per capita (PPP)",
    X_value = health_exp,
    Median = med_AbFever,
    Lower = lower_AbFever,
    Upper = upper_AbFever,
    Projected = NA_real_
  )


# ------------------------------------------------------------
# Panel E: LRTI by DDD per episode
# ------------------------------------------------------------

fig3_E <- group_df_LRI %>%
  transmute(
    Panel = "E",
    Syndrome = "LRTI",
    Malaria_status = "Not applicable",
    X_variable = "DDD per episode",
    X_value = seq_DDD,
    Median = med_AbLRI,
    Lower = lower_AbLRI,
    Upper = upper_AbLRI,
    Projected = NA_real_
  )


# ------------------------------------------------------------
# Panel F: AFI, malaria non-endemic
# ------------------------------------------------------------

fig3_F <- group_df_Fever_0 %>%
  transmute(
    Panel = "F",
    Syndrome = "AFI",
    Malaria_status = "Malaria non-endemic",
    X_variable = "DDD per episode",
    X_value = seq_DDD,
    Median = med_AbFever,
    Lower = lower_AbFever,
    Upper = upper_AbFever,
    Projected = NA_real_
  )


# ------------------------------------------------------------
# Panel G: AFI, malaria endemic
# ------------------------------------------------------------

fig3_G <- group_df_Fever_1 %>%
  transmute(
    Panel = "G",
    Syndrome = "AFI",
    Malaria_status = "Malaria endemic",
    X_variable = "DDD per episode",
    X_value = seq_DDD,
    Median = med_AbFever,
    Lower = lower_AbFever,
    Upper = upper_AbFever,
    Projected = NA_real_
  )


# ------------------------------------------------------------
# Panel H: LRTI observed vs projected
# ------------------------------------------------------------

fig3_H <- group_df_LRI_2 %>%
  transmute(
    Panel = "H",
    Syndrome = "LRTI",
    Malaria_status = "Not applicable",
    X_variable = "DDD per episode",
    X_value = seq_DDD,
    Median = 1 - lack_med_AbLRI,
    Lower = NA_real_,
    Upper = NA_real_,
    Projected = 1 - ideal
  )


# ------------------------------------------------------------
# Panel I: AFI, malaria non-endemic, observed vs projected
# ------------------------------------------------------------

fig3_I <- group_df_Fever_0_2 %>%
  transmute(
    Panel = "I",
    Syndrome = "AFI",
    Malaria_status = "Malaria non-endemic",
    X_variable = "DDD per episode",
    X_value = seq_DDD,
    Median = 1 - lack_med_AbFever,
    Lower = NA_real_,
    Upper = NA_real_,
    Projected = 1 - ideal
  )


# ------------------------------------------------------------
# Panel J: AFI, malaria endemic, observed vs projected
# ------------------------------------------------------------

fig3_J <- group_df_Fever_1_2 %>%
  transmute(
    Panel = "J",
    Syndrome = "AFI",
    Malaria_status = "Malaria endemic",
    X_variable = "DDD per episode",
    X_value = seq_DDD,
    Median = 1 - lack_med_AbFever,
    Lower = NA_real_,
    Upper = NA_real_,
    Projected = 1 - ideal
  )


# ------------------------------------------------------------
# Combine all panels
# ------------------------------------------------------------

SA_Figure_3_source_data <- bind_rows(
  fig3_A,
  fig3_B,
  fig3_C,
  fig3_D,
  fig3_E,
  fig3_F,
  fig3_G,
  fig3_H,
  fig3_I,
  fig3_J
) %>%
  arrange(Panel, X_value)


# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

write.csv(
  SA_Figure_3_source_data,
  "outputs/plot_objects/SA_Figure_6_source_data.csv",
  row.names = FALSE
)





# ------------------------1. P(Antibiotic treatment|infection) ~ Antibiotic consumption/incidence-----------------

pdf("figures/Figure_S6.pdf", width = 7.5, height = 7)


layout(matrix(c(
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
  0, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,10,10,10,10,10,10,10,10,10,10,10,11,11,11,11,11,11,11,11,11,11,11,
  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
  0,12,12,12,12,12,12,12,12,12,12,12,13,13,13,13,13,13,13,13,13,13,13,14,14,14,14,14,14,14,14,14,14,14,15,15,15,15,15,15,15,15,15,15,15), 
  nrow = 6, byrow = TRUE), 
  widths = c(1.2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), 
  heights = c(0.1,1,0.1,1,0.1,1))

#### 1. (title 1) ####
par(mar = c(0, 0, 0, 0))
plot.new()
text(x = -0.035, y = 0.5, labels = "Antibiotic access by GDP and health expenditure", 
     adj = 0, cex = 0.9, font = 2, xpd = TRUE)

#### 2. (title 2)  ####
par(mar = c(0, 0, 0, 0))
plot.new()
text(x = -0.035, y = 0.5, labels = "Antibiotic access and antibiotic consumption relative to incidence", 
     adj = 0, cex = 0.9, font = 2, xpd = TRUE)

#### 3. (title 3)  ####
par(mar = c(0, 0, 0, 0))
plot.new()
text(x = -0.035, y = 0.5, labels = "Observed and expected antibiotic access, and antibiotic consumption relative to incidence", 
     adj = 0, cex = 0.9, font = 2, xpd = TRUE)








#### 4 (figure A) ####

par(mar = c(3, 2.5, 1.5, 1))  
par(mgp = c(3, 0.3, 0))     
par(tck = -0.025)
par(lwd = 1)

xmin <- log10(1000)         ;     xmax <- log10(30000)       
ymin <- 0 ; ymax <- 1


plot(med_AbLRI ~ GDP, 
     data = LRI_gdp_plot,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.75, 0.75, 0.75, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")




xaxis <- log10(c(1000,10000,30000))
xticks <- log10(c(1000, 10000))
xlabels <- round(c(1000, 2000, 10000, 20000, 30000), digits = 0)
minor_xlabels <- as.vector(c(sapply(c(1000), function(v) v * 2:9),sapply(c(10000), function(v) v * 2:3)))

axis(1, at = xaxis, labels = FALSE, las = 2, pos = ymin, tck = 0)
axis(1, at = xticks, labels = FALSE, las = 2, pos = ymin)
axis(1, at = log10(minor_xlabels), labels = FALSE, las = 2, pos = ymin, tck = -0.015, lwd = 1)
text(x = log10(xlabels), y = ymin-(ymax-ymin)*0.065,  
     labels = xlabels/1000, srt = 0, adj = 0.5, xpd = TRUE, cex = 0.65)

axis(2,  at = seq(ymin, ymax, by = (ymax - ymin)/4), labels = FALSE, las = 2, pos = xmin)
text(y = seq(ymin, ymax, by = (ymax-ymin)/4), x = xmin-0.04*(xmax-xmin),  
     labels = format(round(seq(ymin, ymax, by = (ymax-ymin)/4), digits = 2), nsmall = 1), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines(med_AbLRI ~ log10(GDP), data = LRI_gdp_plot, col = "#0072B2", lwd = 1.5)

# Add shaded area for confidence interval
polygon(c(log10(LRI_gdp_plot$GDP), rev(log10(LRI_gdp_plot$GDP))),
        c(LRI_gdp_plot$lower_AbLRI, rev(LRI_gdp_plot$upper_AbLRI)),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)

mtext(expression("Antibiotic treatment probability"), side = 2, line = 2, cex = 0.50)  
mtext("GDP per capita (thousands, PPP)", side = 1, line = 1, cex = 0.50)
mtext("A. LRTI", at = xmin - 0.33*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.50, font = 1)




#### 5 (figure B) ####

xmin <- log10(1000)         ;     xmax <- log10(30000)       
ymin <- 0 ; ymax <- 1


plot(med_AbFever ~ GDP, 
     data = Fever_gdp_plot,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.75, 0.75, 0.75, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")



xaxis <- log10(c(1000,10000,30000))
xticks <- log10(c(1000, 10000))
xlabels <- round(c(1000, 2000, 10000, 20000, 30000), digits = 0)
minor_xlabels <- as.vector(c(sapply(c(1000), function(v) v * 2:9),sapply(c(10000), function(v) v * 2:3)))

axis(1, at = xaxis, labels = FALSE, las = 2, pos = ymin, tck = 0)
axis(1, at = xticks, labels = FALSE, las = 2, pos = ymin)
axis(1, at = log10(minor_xlabels), labels = FALSE, las = 2, pos = ymin, tck = -0.015, lwd = 1)
text(x = log10(xlabels), y = ymin-(ymax-ymin)*0.065,  
     labels = xlabels/1000, srt = 0, adj = 0.5, xpd = TRUE, cex = 0.65)


axis(2,  at = seq(ymin, ymax, by = (ymax - ymin)/4), labels = FALSE, las = 2, pos = xmin)
text(y = seq(ymin, ymax, by = (ymax-ymin)/4), x = xmin-0.04*(xmax-xmin),  
     labels = format(round(seq(ymin, ymax, by = (ymax-ymin)/4), digits = 2), nsmall = 1), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines(med_AbFever ~ log10(GDP), data = Fever_gdp_plot, col = "#0072B2", lwd = 1.5)

# Add shaded area for confidence interval
polygon(c(log10(Fever_gdp_plot$GDP), rev(log10(Fever_gdp_plot$GDP))),
        c(Fever_gdp_plot$lower_AbFever, rev(Fever_gdp_plot$upper_AbFever)),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)

mtext("GDP per capita (thousands, PPP)", side = 1, line = 1, cex = 0.50)
mtext("B. AFI", at = xmin - 0.33*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.50, font = 1)










#### 6 (figure C) ####

par(mar = c(3, 2.5, 1.5, 1))  
par(mgp = c(3, 0.3, 0))     
par(tck = -0.025)
par(lwd = 1)

xmin <- log10(50)         ;     xmax <- log10(2000)       
ymin <- 0 ; ymax <- 1


plot(med_AbLRI ~ health_exp, 
     data = LRI_health_exp_plot,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.75, 0.75, 0.75, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")



xaxis <- log10(c(50,100,1000, 2000))
xticks <- log10(c(100, 1000))
xlabels <- round(c(50, 100, 200, 1000, 2000), digits = 0)
minor_xlabels <- as.vector(c(sapply(c(10), function(v) v * 5:9),sapply(c(100), function(v) v * 2:9),sapply(c(1000), function(v) v * 2)))

axis(1, at = xaxis, labels = FALSE, las = 2, pos = ymin, tck = 0)
axis(1, at = xticks, labels = FALSE, las = 2, pos = ymin)
axis(1, at = log10(minor_xlabels), labels = FALSE, las = 2, pos = ymin, tck = -0.015, lwd = 1)
text(x = log10(xlabels), y = ymin-(ymax-ymin)*0.065,  
     labels = xlabels, srt = 0, adj = 0.5, xpd = TRUE, cex = 0.65)

axis(2,  at = seq(ymin, ymax, by = (ymax - ymin)/4), labels = FALSE, las = 2, pos = xmin)
text(y = seq(ymin, ymax, by = (ymax-ymin)/4), x = xmin-0.04*(xmax-xmin),  
     labels = format(round(seq(ymin, ymax, by = (ymax-ymin)/4), digits = 2), nsmall = 1), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines(med_AbLRI ~ log10(health_exp), data = LRI_health_exp_plot, col = "#0072B2", lwd = 1.5)

# Add shaded area for confidence interval
polygon(c(log10(LRI_health_exp_plot$health_exp), rev(log10(LRI_health_exp_plot$health_exp))),
        c(LRI_health_exp_plot$lower_AbLRI, rev(LRI_health_exp_plot$upper_AbLRI)),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)

mtext("Health expenditure per capita (PPP)", side = 1, line = 1, cex = 0.50)
mtext("C. LRTI", at = xmin - 0.33*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.50, font = 1)




#### 7 (figure D) ####

xmin <- log10(50)         ;     xmax <- log10(2000)       
ymin <- 0 ; ymax <- 1


plot(med_AbFever ~ health_exp, 
     data = Fever_health_exp_plot,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.75, 0.75, 0.75, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")





xaxis <- log10(c(50,100,1000, 2000))
xticks <- log10(c(100, 1000))
xlabels <- round(c(50, 100, 200, 1000, 2000), digits = 0)
minor_xlabels <- as.vector(c(sapply(c(10), function(v) v * 5:9),sapply(c(100), function(v) v * 2:9),sapply(c(1000), function(v) v * 2)))

axis(1, at = xaxis, labels = FALSE, las = 2, pos = ymin, tck = 0)
axis(1, at = xticks, labels = FALSE, las = 2, pos = ymin)
axis(1, at = log10(minor_xlabels), labels = FALSE, las = 2, pos = ymin, tck = -0.015, lwd = 1)
text(x = log10(xlabels), y = ymin-(ymax-ymin)*0.065,  
     labels = xlabels, srt = 0, adj = 0.5, xpd = TRUE, cex = 0.65)


axis(2,  at = seq(ymin, ymax, by = (ymax - ymin)/4), labels = FALSE, las = 2, pos = xmin)
text(y = seq(ymin, ymax, by = (ymax-ymin)/4), x = xmin-0.04*(xmax-xmin),  
     labels = format(round(seq(ymin, ymax, by = (ymax-ymin)/4), digits = 2), nsmall = 1), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines(med_AbFever ~ log10(health_exp), data = Fever_health_exp_plot, col = "#0072B2", lwd = 1.5)

# Add shaded area for confidence interval
polygon(c(log10(Fever_health_exp_plot$health_exp), rev(log10(Fever_health_exp_plot$health_exp))),
        c(Fever_health_exp_plot$lower_AbFever, rev(Fever_health_exp_plot$upper_AbFever)),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)

mtext("Health expenditure per capita (PPP)", side = 1, line = 1, cex = 0.50)
mtext("D. AFI", at = xmin - 0.33*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.50, font = 1)



#### 8. (Fig E)  ####
xmin <- log10(10)         ;     xmax <- log10(300)       
ymin <- 0 ; ymax <- 1

par(mar = c(3, 2.5, 1.5, 1))  
par(mgp = c(3, 0.3, 0))     
par(tck = -0.025)
par(lwd = 1)


plot(med_AbLRI ~ seq_DDD, 
     data = group_df_LRI,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.75, 0.75, 0.75, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")



xaxis <- log10(c(10,300))
xticks <- log10(c(10, 100))
xlabels <- round(c(10, 20, 100, 200, 300), digits = 0)
minor_xlabels <- as.vector(c(sapply(c(10), function(v) v * 2:9), sapply(c(100), function(v) v * 2:3)))

axis(1, at = xaxis, labels = FALSE, las = 2, pos = ymin, tck = 0)
axis(1, at = xticks, labels = FALSE, las = 2, pos = ymin)
axis(1, at = log10(minor_xlabels), labels = FALSE, las = 2, pos = ymin, tck = -0.015, lwd = 1)
text(x = log10(xlabels), y = ymin-(ymax-ymin)*0.065,  
     labels = xlabels, srt = 0, adj = 0.5, xpd = TRUE, cex = 0.65)

axis(2,  at = seq(ymin, ymax, by = (ymax - ymin)/4), labels = FALSE, las = 2, pos = xmin)
text(y = seq(ymin, ymax, by = (ymax-ymin)/4), x = xmin-0.04*(xmax-xmin),  
     labels = format(round(seq(ymin, ymax, by = (ymax-ymin)/4), digits = 2), nsmall = 1), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)



lines(med_AbLRI ~ log10(seq_DDD), data = group_df_LRI, col = "#0072B2", lwd = 1.5)

# Add shaded area for confidence interval
polygon(c(log10(group_df_LRI$seq_DDD), rev(log10(group_df_LRI$seq_DDD))),
        c(group_df_LRI$lower_AbLRI, rev(group_df_LRI$upper_AbLRI)),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)




# Probability of antibotics treatment by antibiotic consumption (DDD) per incidence
mtext("E. LRTI", at = xmin - 0.3*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.5, font = 1)
mtext("DDD per episode", side = 1, line = 1, cex = 0.5)
mtext(expression("Antibiotic treatment probability"), side = 2, line = 2, cex = 0.5)  



#### 9. (Fig F) ####
xmin <- log10(1)         ;     xmax <- log10(7)       
ymin <- 0 ; ymax <- 1



plot(med_AbFever ~ seq_DDD, 
     data = group_df_Fever_0,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.75, 0.75, 0.75, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")



xaxis <- log10(c(1,7))
xticks <- log10(c(1))
xlabels <- round(c(1, 2, 7), digits = 0)
minor_xlabels <- as.vector(sapply(c(1), function(v) v * 2:7))

axis(1, at = xaxis, labels = FALSE, las = 2, pos = ymin, tck = 0)
axis(1, at = xticks, labels = FALSE, las = 2, pos = ymin)
axis(1, at = log10(minor_xlabels), labels = FALSE, las = 2, pos = ymin, tck = -0.015, lwd = 1)
text(x = log10(xlabels), y = ymin-(ymax-ymin)*0.065,  
     labels = xlabels, srt = 0, adj = 0.5, xpd = TRUE, cex = 0.65)

axis(2,  at = seq(ymin, ymax, by = (ymax - ymin)/4), labels = FALSE, las = 2, pos = xmin)
text(y = seq(ymin, ymax, by = (ymax-ymin)/4), x = xmin-0.04*(xmax-xmin),  
     labels = format(round(seq(ymin, ymax, by = (ymax-ymin)/4), digits = 2), nsmall = 1), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)


lines(med_AbFever ~ log10(seq_DDD), data = group_df_Fever_0, col = "#0072B2", lwd = 1.5)

# Add shaded area for confidence interval
polygon(c(log10(group_df_Fever_0$seq_DDD), rev(log10(group_df_Fever_0$seq_DDD))),
        c(group_df_Fever_0$lower_AbFever, rev(group_df_Fever_0$upper_AbFever)),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)




mtext("F. AFI (malaria non-endemic countries)", at = xmin - 0.3*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.5, font = 1)
mtext("DDD per episode", side = 1, line = 1, cex = 0.5)





#### 10. (Fig G) ####
xmin <- log10(0.9)         ;     xmax <- log10(4)       
ymin <- 0 ; ymax <- 1



plot(med_AbFever ~ seq_DDD, 
     data = group_df_Fever_1,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.75, 0.75, 0.75, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")



xticks <- log10(c(1))
xlabels <- c(0.9, 1, 2, 4)
minor_xlabels <- as.vector(c(0.9,1,2,3,4))

axis(1, at = xticks, labels = FALSE, las = 2, pos = ymin)
axis(1, at = log10(minor_xlabels), labels = FALSE, las = 2, pos = ymin, tck = -0.015, lwd = 1)
text(x = log10(xlabels), y = ymin-(ymax-ymin)*0.065,  
     labels = xlabels, srt = 0, adj = 0.5, xpd = TRUE, cex = 0.65)

axis(2,  at = seq(ymin, ymax, by = (ymax - ymin)/4), labels = FALSE, las = 2, pos = xmin)
text(y = seq(ymin, ymax, by = (ymax-ymin)/4), x = xmin-0.04*(xmax-xmin),  
     labels = format(round(seq(ymin, ymax, by = (ymax-ymin)/4), digits = 2), nsmall = 1), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)


lines(med_AbFever ~ log10(seq_DDD), data = group_df_Fever_1, col = "#0072B2", lwd = 1.5)

# Add shaded area for confidence interval
polygon(c(log10(group_df_Fever_1$seq_DDD), rev(log10(group_df_Fever_1$seq_DDD))),
        c(group_df_Fever_1$lower_AbFever, rev(group_df_Fever_1$upper_AbFever)),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)




mtext("G. AFI (malaria-endemic countries)", at = xmin - 0.3*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.5, font = 1)
mtext("DDD per episode", side = 1, line = 1, cex = 0.5)





#### 11. (legend) ####
par(mar = c(0, 0, 0, 0))
plot.new()




#### 12. (Fig H)  ####
xmin <- log10(10)         ;     xmax <- log10(300)       
ymin <- 0 ; ymax <- 1

par(mar = c(3, 2.5, 1.5, 1))  
par(mgp = c(3, 0.3, 0))     
par(tck = -0.025)
par(lwd = 1)


plot(lack_med_AbLRI ~ seq_DDD, 
     data = group_df_LRI_2,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.75, 0.75, 0.75, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")




# Example: log10 scale

xaxis <- log10(c(10,200))
xticks <- log10(c(10, 100))
xlabels <- round(c(10, 20, 100, 200), digits = 0)
minor_xlabels <- as.vector(c(sapply(c(10), function(v) v * 2:9), sapply(c(100), function(v) v * 2)))

axis(1, at = xaxis, labels = FALSE, las = 2, pos = ymin, tck = 0)
axis(1, at = xticks, labels = FALSE, las = 2, pos = ymin)
axis(1, at = log10(minor_xlabels), labels = FALSE, las = 2, pos = ymin, tck = -0.015, lwd = 1)
text(x = log10(xlabels), y = ymin-(ymax-ymin)*0.065,  
     labels = xlabels, srt = 0, adj = 0.5, xpd = TRUE, cex = 0.65)

axis(2,  at = seq(ymin, ymax, by = (ymax - ymin)/4), labels = FALSE, las = 2, pos = xmin)
text(y = seq(ymin, ymax, by = (ymax-ymin)/4), x = xmin-0.04*(xmax-xmin),  
     labels = format(round(seq(ymin, ymax, by = (ymax-ymin)/4), digits = 2), nsmall = 1), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines((1-lack_med_AbLRI) ~ log10(seq_DDD), data = group_df_LRI_2, col = "#0072B2", lwd = 1.5)
lines((1-ideal) ~ log10(seq_DDD), data = group_df_LRI_2, col = "darkblue", lwd = 1.5)

polygon(c(log10(group_df_LRI_2$seq_DDD), rev(log10(group_df_LRI_2$seq_DDD))),
        c((1-group_df_LRI_2$ideal), rev((1-group_df_LRI_2$lack_med_AbLRI))),
        col = rgb(0/255, 0/255, 139/255, 0.25), border = NA)




# Probability of antibotics treatment by antibiotic consumption (DDD) per incidence
mtext("H. LRTI", at = xmin - 0.3*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.5, font = 1)
mtext("DDD per episode", side = 1, line = 1, cex = 0.5)
mtext(expression("Antibiotic treatment probability"), side = 2, line = 2, cex = 0.5)  






#### 13. (Fig I) ####
xmin <- log10(1)         ;     xmax <- log10(7)       
ymin <- 0 ; ymax <- 1


plot(med_AbFever ~ seq_DDD, 
     data = group_df_Fever_0_2,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.75, 0.75, 0.75, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1", 
     axes=F,ann=F,
     main = "")



xaxis <- log10(c(1,7))
xticks <- log10(c(1))
xlabels <- round(c(1, 2, 7), digits = 0)
minor_xlabels <- as.vector(sapply(c(1), function(v) v * 2:7))

axis(1, at = xaxis, labels = FALSE, las = 2, pos = ymin, tck = 0)
axis(1, at = xticks, labels = FALSE, las = 2, pos = ymin)
axis(1, at = log10(minor_xlabels), labels = FALSE, las = 2, pos = ymin, tck = -0.015, lwd = 1)
text(x = log10(xlabels), y = ymin-(ymax-ymin)*0.065,  
     labels = xlabels, srt = 0, adj = 0.5, xpd = TRUE, cex = 0.65)

axis(2,  at = seq(ymin, ymax, by = (ymax - ymin)/4), labels = FALSE, las = 2, pos = xmin)
text(y = seq(ymin, ymax, by = (ymax-ymin)/4), x = xmin-0.04*(xmax-xmin),  
     labels = format(round(seq(ymin, ymax, by = (ymax-ymin)/4), digits = 2), nsmall = 1), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines((1-lack_med_AbFever) ~ log10(seq_DDD), data = group_df_Fever_0_2, col = "#0072B2", lwd = 1.5)
lines((1-ideal) ~ log10(seq_DDD), data = group_df_Fever_0_2, col = "darkblue", lwd = 1.5)

# Add shaded area for confidence interval
polygon(c(log10(group_df_Fever_0_2$seq_DDD), rev(log10(group_df_Fever_0_2$seq_DDD))),
        c((1-group_df_Fever_0_2$ideal), rev((1-group_df_Fever_0_2$lack_med_AbFever))),
        col = rgb(0/255, 0/255, 139/255, 0.25), border = NA)



mtext("I. AFI (malaria non-endemic countries)", at = xmin - 0.3*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.5, font = 1)
mtext("DDD per episode", side = 1, line = 1, cex = 0.5)





#### 104 (Fig J) ####
xmin <- log10(0.9)         ;     xmax <- log10(4)       
ymin <- 0 ; ymax <- 1


plot(med_AbFever ~ seq_DDD, 
     data = group_df_Fever_1,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.75, 0.75, 0.75, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")



xticks <- log10(c(1))
xlabels <- c(0.9, 1, 2, 4)
minor_xlabels <- as.vector(c(0.9,1,2,3,4))

axis(1, at = xticks, labels = FALSE, las = 2, pos = ymin)
axis(1, at = log10(minor_xlabels), labels = FALSE, las = 2, pos = ymin, tck = -0.015, lwd = 1)
text(x = log10(xlabels), y = ymin-(ymax-ymin)*0.065,  
     labels = xlabels, srt = 0, adj = 0.5, xpd = TRUE, cex = 0.65)

axis(2,  at = seq(ymin, ymax, by = (ymax - ymin)/4), labels = FALSE, las = 2, pos = xmin)
text(y = seq(ymin, ymax, by = (ymax-ymin)/4), x = xmin-0.04*(xmax-xmin),  
     labels = format(round(seq(ymin, ymax, by = (ymax-ymin)/4), digits = 2), nsmall = 1), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines((1-lack_med_AbFever) ~ log10(seq_DDD), data = group_df_Fever_1_2, col = "#0072B2", lwd = 1.5)
lines((1-ideal) ~ log10(seq_DDD), data = group_df_Fever_1_2, col = "darkblue", lwd = 1.5)

# Add shaded area for confidence interval
polygon(c(log10(group_df_Fever_1_2$seq_DDD), rev(log10(group_df_Fever_1_2$seq_DDD))),
        c((1-group_df_Fever_1_2$ideal), rev((1-group_df_Fever_1_2$lack_med_AbFever))),
        col = rgb(0/255, 0/255, 139/255, 0.25), border = NA)



mtext("J. AFI (malaria-endemic countries)", at = xmin - 0.3*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.5, font = 1)
mtext("DDD per episode", side = 1, line = 1, cex = 0.5)








#### 15. (legend) ####
par(mar = c(0, 0, 0, 0))
plot.new()

legend_labels <- c(
  "Projected value from lowest consumption deciles", 
  "Observed value", 
  "Gap between observed and expected values"
)
legend_colors <- c( "darkblue","#0072B2", "darkblue")
n <- length(legend_labels)


# Coordinates for the custom legend
x0 <- 0.0
x1 <- 0.07
y_start <- 0.61
line_spacing <- 0.07


for (i in seq_len(n)) {
  y <- y_start - (i - 1) * line_spacing
  
  if (i == 3) {
    # Polygon only (no line)
    polygon(
      x = c(x0-0.01, x1+0.01, x1+0.01, x0-0.01),
      y = c(y - 0.02, y - 0.02, y + 0.02, y + 0.02),
      col = rgb(0/255, 0/255, 139/255, 0.25),
      border = NA
    )
  } else {
    segments(x0, y, x1, y, col = legend_colors[i], lwd = 1.5)
  }
  
  # Label
  text(x = x1 + 0.04, y = y, labels = legend_labels[i], adj = 0, cex = 0.56)
}






dev.off()



########## 11. Figure 4 prep ###########


# ------------------with 4 antibiotic classes: BSP, NSP, Cephalosporins, Aminoglycosides ---------

#### 1. prep

rm(list = ls())
gc()
set.seed(123)

requiredPackages = c("dplyr", "tidyr", "data.table", "MASS")
lib.path <- .libPaths()
for(p in requiredPackages){
  if(!require(p, character.only = TRUE)) install.packages(p, lib = lib.path)
  library(p, character.only = TRUE)
}

n_sim <- 10000

#### 2. AbLRI ~ DDD per LRTI

ab.data_LRI <- readRDS("outputs/sensitivity_analysis/SA_ab.data_LRI.UNI.rds")
matrix.LRI <- readRDS("outputs/sensitivity_analysis/SA_U5_DDD_LRI_matrix.rds")

matrix.meta <- ab.data_LRI %>% dplyr::select(Country, Year)
matrix.LRI <- matrix.meta %>% inner_join(matrix.LRI, by = c("Country", "Year"))
matrix.LRI <- matrix.LRI[, c(1:30002, 40003:50002)]

matrix.LRI.sum <- data.frame(matrix(NA_real_, nrow = nrow(matrix.LRI), ncol = (n_sim + 2)))
matrix.LRI.sum[, 1:2] <- as.matrix(matrix.LRI[, 1:2])

# sum antibiotic classes

for (i in 1:n_sim) {
  cols_to_sum <- seq(from = i + 2, to = (4 * n_sim + 2), by = n_sim)
  matrix.LRI.sum[, i + 2] <- as.numeric(rowSums(matrix.LRI[, cols_to_sum]))
}

dim(matrix.LRI.sum)
str(matrix.LRI.sum)

colnames(matrix.LRI.sum)[1:2] <- c("Country", "Year")
colnames(matrix.LRI.sum)[3:ncol(matrix.LRI.sum)] <- paste0("sim_", 1:n_sim)

ab.data_LRI$DDD_sum <- ab.data_LRI %>%
  dplyr::select(
    "U5_DDD_LRI_BSP",
    "U5_DDD_LRI_NSP",
    "U5_DDD_LRI_Cephalosporins",
    "U5_DDD_LRI_Aminoglycosides"
  ) %>%
  rowSums()

meta_LRI <- ab.data_LRI %>% dplyr::select(Country, Year)

meta_LRI$Year <- as.integer(meta_LRI$Year)
matrix.LRI.sum$Year <- as.integer(matrix.LRI.sum$Year)

matrix.LRI.sum <- meta_LRI %>%
  inner_join(matrix.LRI.sum, by = c("Country", "Year"))

### 1-1. DDD(LRI) ~ PC1 + PC2

Ab.colnames_LRI <- c(
  "U5_DDD_LRI_BSP",
  "U5_DDD_LRI_NSP",
  "U5_DDD_LRI_Cephalosporins",
  "U5_DDD_LRI_Aminoglycosides"
)

LRI.DDD_sim <- matrix.LRI %>%
  dplyr::select(-Country, -Year) %>%
  as.data.frame()

mod_LRI_list <- list()
coef_matrix_LRI_DDD <- list()

# Each model has 3 coefficients:
# (Intercept), PC1, PC2

for(i in 1:length(Ab.colnames_LRI)) {
  
  DDD_name <- Ab.colnames_LRI[i]
  
  ab.data_LRI$DDD_loop <- ab.data_LRI[[DDD_name]]
  
  mod <- glm(
    log(DDD_loop) ~ PC1 + PC2,
    data = ab.data_LRI
  )
  
  mod_LRI_list[[DDD_name]] <- mod
  
  coef_mat <- matrix(
    NA_real_,
    nrow = n_sim,
    ncol = length(coef(mod))
  )
  
  coef_names <- names(coef(mod))
  colnames(coef_mat) <- paste0(DDD_name, "_", coef_names)
  
  coef_matrix_LRI_DDD[[DDD_name]] <- coef_mat
}

coef_matrix_LRI_DDD <- do.call(cbind, coef_matrix_LRI_DDD)

# Simulate coefficients for each antibiotic class

for (j in 1:4) {
  
  for (i in 1:n_sim) {
    
    ab.data_LRI$LRI.DDD_sim <- LRI.DDD_sim[, n_sim * (j - 1) + i]
    
    mod <- glm(
      log(LRI.DDD_sim) ~ PC1 + PC2,
      data = ab.data_LRI
    )
    
    coef_mean <- coef(mod)
    coef_vcov <- vcov(mod)
    
    coef_matrix_LRI_DDD[
      i,
      3 * (j - 1) + (1:3)
    ] <- MASS::mvrnorm(
      1,
      mu = coef_mean,
      Sigma = coef_vcov
    )
    
    if (i %% 100 == 0) cat("Completed simulation:", j, "_", i, "\n")
    
  }
}

coef_df_LRI_DDD <- as.data.frame(coef_matrix_LRI_DDD)

head(coef_df_LRI_DDD)
str(coef_df_LRI_DDD)

# Check: 4 classes × 3 coefficients = 12 columns

stopifnot(ncol(coef_df_LRI_DDD) == 12)

saveRDS(
  coef_df_LRI_DDD,
  "outputs/sensitivity_analysis/SA_coef_df_LRI_DDD_by_4class.rds"
)



#### 2: time consuming (~50 min; jump to step 3 if modification is not required)

rm(list = ls())
gc()

requiredPackages = c("dplyr", "tidyr", "data.table", "MASS")
lib.path <- .libPaths()
for(p in requiredPackages){
  if(!require(p, character.only = TRUE)) install.packages(p, lib = lib.path)
  library(p, character.only = TRUE)
}

n_sim <- 10000

UNI_indexed <- readRDS("outputs/sensitivity_analysis/SA_UNI_indexed.rds")
coef_df_LRI <- readRDS("outputs/sensitivity_analysis/SA_coef_df_LRI.rds")
coef_df_LRI_DDD <- readRDS("outputs/sensitivity_analysis/SA_coef_df_LRI_DDD_by_4class.rds")

UNI_indexed_Abaccess <- UNI_indexed %>% filter(Abaccess == 1)
UNI_indexed_Abdata <- UNI_indexed %>% filter(ab.data == 1)
UNI_indexed_both <- UNI_indexed %>% filter(Abaccess == 1, ab.data == 1)

# Confirm treatment model uses intercept + PC1-PC2

stopifnot(
  identical(
    names(coef_df_LRI),
    c("(Intercept)", "PC1", "PC2")
  )
)

# Confirm class-specific DDD models use 3 coefficients per class

stopifnot(
  ncol(coef_df_LRI_DDD) == 12
)

# Common PC1-PC2 prediction matrix

X_pred <- model.matrix(
  ~ PC1 + PC2,
  data = UNI_indexed
)

X_pred_LRI <- X_pred[, names(coef_df_LRI), drop = FALSE]

# Antibiotic treatment linear predictor

LRI_matrix <- as.matrix(X_pred_LRI) %*%
  t(as.matrix(coef_df_LRI))

# DDD linear predictors for each antibiotic class

LRI_DDD_list <- list()

for(i in 1:4){
  
  coef_cols <- 3 * (i - 1) + (1:3)
  
  LRI_DDD_list[[i]] <- as.matrix(X_pred) %*%
    t(as.matrix(coef_df_LRI_DDD[, coef_cols, drop = FALSE]))
}

LRI_DDD_matrix <- do.call(cbind, LRI_DDD_list)

# Logistic model for antibiotic treatment

AbLRI_sim <- plogis(LRI_matrix)

# Log-linear model for DDD

LRI_DDD_sim <- exp(LRI_DDD_matrix)

max(AbLRI_sim)

n_row <- nrow(LRI_DDD_sim)

LRI_DDD_ratio_sim <- matrix(
  NA_real_,
  nrow = n_row,
  ncol = 4 * n_sim
)

# ratio between antibiotic classes (standard: BSP)

for (i in 1:n_sim) {
  
  sim_cols <- seq(i, 4 * n_sim, by = n_sim)
  
  sim_mat <- LRI_DDD_sim[, sim_cols, drop = FALSE]
  
  standard.ab <- sim_mat[, 1]
  
  LRI_DDD_ratio_sim[, sim_cols] <- sim_mat / standard.ab
}

saveRDS(
  LRI_DDD_ratio_sim,
  "outputs/plot_df/p_best/SA_ratio_LRI4.rds"
)

dim(LRI_DDD_ratio_sim)

# Estimate probability for BSP

p <- seq(0, 1, length.out = 1001)

p_best_LRI <- matrix(
  NA_real_,
  nrow = n_row,
  ncol = n_sim
)

for (i in 1:n_sim) {
  
  sim_cols <- seq(i, 4 * n_sim, by = n_sim)
  
  ratio_i <- LRI_DDD_ratio_sim[, sim_cols, drop = FALSE]
  
  error_i <- matrix(
    NA_real_,
    nrow = n_row,
    ncol = length(p)
  )
  
  for (k in seq_along(p)) {
    
    prod_term <- 1 - p[k] * ratio_i
    
    log_prod_term <- log(prod_term)
    
    pred_prob <- 1 - exp(rowSums(log_prod_term))
    
    error_i[, k] <- abs(
      AbLRI_sim[, i] - pred_prob
    )
    
  }
  
  best_k <- apply(error_i, 1, which.min)
  
  p_best_LRI[, i] <- p[best_k]
  
  cat("Completed sim:", i, "\n")
}

saveRDS(
  p_best_LRI,
  "outputs/plot_df/p_best/SA_p_LRI4.rds"
)



### 3. Matrix

rm(list = ls())
gc()

requiredPackages = c("dplyr", "tidyr", "data.table", "MASS")
lib.path <- .libPaths()
for(p in requiredPackages){
  if(!require(p, character.only = TRUE)) install.packages(p, lib = lib.path)
  library(p, character.only = TRUE)
}

n_sim <- 10000

UNI_indexed <- readRDS("outputs/sensitivity_analysis/SA_UNI_indexed.rds")
p_best_LRI <- readRDS("outputs/plot_df/p_best/SA_p_LRI4.rds")
LRI_DDD_ratio_sim <- readRDS("outputs/plot_df/p_best/SA_ratio_LRI4.rds")

LRI_Ab_probability.matrix <- list()

for (i in 1:4) {
  
  LRI_ratio_single_Ab <-
    LRI_DDD_ratio_sim[, (n_sim * (i - 1) + 1):(n_sim * i), drop = FALSE]
  
  LRI_Ab_probability.matrix[[i]] <-
    p_best_LRI * LRI_ratio_single_Ab
}

names(LRI_Ab_probability.matrix) <- c(
  "BSP",
  "NSP",
  "Cephalosporins",
  "Aminoglycosides"
)

n_row <- nrow(LRI_Ab_probability.matrix[[1]])

class_Ab_LRI <- matrix(
  NA_real_,
  nrow = n_row,
  ncol = 12
)

digit_number <- c(4, 4, 8, 4)

for(i in 1:4) {
  
  single_Ab_LRI <- LRI_Ab_probability.matrix[[i]]
  
  class_Ab_LRI[, 3 * i - 2] <-
    round(
      apply(single_Ab_LRI, 1, median),
      digits = digit_number[i]
    )
  
  class_Ab_LRI[, 3 * i - 1] <-
    round(
      apply(single_Ab_LRI, 1, quantile, probs = 0.025),
      digits = digit_number[i]
    )
  
  class_Ab_LRI[, 3 * i] <-
    round(
      apply(single_Ab_LRI, 1, quantile, probs = 0.975),
      digits = digit_number[i]
    )
  
  print(i)
}

class_Ab_LRI_df <- as_tibble(class_Ab_LRI)

Ab.names <- c(
  "BSP",
  "NSP",
  "Cephalosporins",
  "Aminoglycosides"
)

colnames(class_Ab_LRI_df) <- as.vector(
  rbind(
    paste0(Ab.names, "_Median"),
    paste0(Ab.names, "_Lower"),
    paste0(Ab.names, "_Upper")
  )
)

combine_quantiles <- function(point, lower, upper, digits = 4) {
  
  sprintf(
    paste0(
      "%.", digits, "f (%.", digits, "f, %.", digits, "f)"
    ),
    point,
    lower,
    upper
  )
}

combined_LRI <- tibble(.rows = n_row)

for (i in 1:4) {
  
  ab_group <- Ab.names[i]
  
  col_base <- 3 * (i - 1)
  
  point_LRI <- class_Ab_LRI_df[[col_base + 1]]
  lower_LRI <- class_Ab_LRI_df[[col_base + 2]]
  upper_LRI <- class_Ab_LRI_df[[col_base + 3]]
  
  digits <- digit_number[i]
  
  combined_LRI[[ab_group]] <-
    combine_quantiles(
      point_LRI,
      lower_LRI,
      upper_LRI,
      digits = digits
    )
}

combined_LRI <- as.data.frame(combined_LRI)

combined_LRI_result <- cbind(
  UNI_indexed[, 1:2],
  combined_LRI
)



##### 0 #####
rm(list = ls())
gc()

requiredPackages = c("dplyr", "tidyr", "tibble", "data.table", "stringr", "MASS")
lib.path <- .libPaths()
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p, lib = lib.path)
  library(p,character.only = TRUE)
}



##### 1 #####
rm(list = ls())
gc()

n_sim <- 10000
set.seed(123)

#---------------------------prepare data for figure------------------------
Ab_num <- 4

UNI_indexed <- readRDS("outputs/sensitivity_analysis/SA_UNI_indexed.rds")

UNI_indexed_Abaccess <- UNI_indexed %>% filter(Abaccess == 1)
UNI_indexed_Abdata <- UNI_indexed %>% filter(ab.data == 1)
UNI_indexed_both <- UNI_indexed %>% filter(Abaccess == 1, ab.data == 1)

# outputs generated from the current PC1-PC2 models
p_best_LRI <- readRDS("outputs/plot_df/p_best/SA_p_LRI4.rds")
LRI_DDD_ratio_sim <- readRDS("outputs/plot_df/p_best/SA_ratio_LRI4.rds")

stopifnot(
  nrow(p_best_LRI) == nrow(UNI_indexed),
  ncol(p_best_LRI) == n_sim,
  nrow(LRI_DDD_ratio_sim) == nrow(UNI_indexed),
  ncol(LRI_DDD_ratio_sim) == Ab_num * n_sim
)

LRI_Ab_probability.matrix <- list()

for (i in 1:Ab_num) {
  LRI_ratio_single_Ab <- LRI_DDD_ratio_sim[,(n_sim*(i-1)+1):(n_sim*i)]
  LRI_Ab_probability.matrix[[i]] <- p_best_LRI*LRI_ratio_single_Ab
}

Ab.names <- c("BSP","NSP","Cephalosporins","Aminoglycosides")
names(LRI_Ab_probability.matrix) <- Ab.names

for (i in 1:Ab_num) {
  ab_name <- Ab.names[i]
  AbLRI_long <- as.vector(LRI_Ab_probability.matrix[[i]])
  assign(paste0("AbLRI_long_", ab_name), AbLRI_long)
}

# Create matching indices
Year_index <- rep(UNI_indexed$Year, times = n_sim)
Country_index <- rep(UNI_indexed$Country, times = n_sim)
GDP_index <- rep(UNI_indexed$GDP, times = n_sim)
Income_index <- rep(UNI_indexed$Income, times = n_sim)
Malaria_index <- rep(UNI_indexed$Endemic, times = n_sim)
GBD_index <- rep(UNI_indexed$gbd_region, times = n_sim)


# Combine into dataframe
sim_df_LRI <- tibble(
  Country = Country_index,
  Year = Year_index,
  GDP = GDP_index,
  Income = Income_index,
  Malaria = Malaria_index,
  GBD_region = GBD_index,
  AbLRI_BSP = AbLRI_long_BSP,
  AbLRI_NSP = AbLRI_long_NSP,
  AbLRI_Cephalosporins = AbLRI_long_Cephalosporins,
  AbLRI_Aminoglycosides = AbLRI_long_Aminoglycosides,
  sim = rep(1:n_sim, each = nrow(UNI_indexed))
)


rm(p_best_LRI, LRI_DDD_ratio_sim)
rm_list <- c(paste0("AbLRI_long_", Ab.names))
rm(list = rm_list, envir = .GlobalEnv)

Country_Year_Antibiotics <- tibble(
  Country = sim_df_LRI$Country,
  Year = sim_df_LRI$Year,
  sim = sim_df_LRI$sim,
  GDP = sim_df_LRI$GDP,
  Income = sim_df_LRI$Income,
  Malaria = sim_df_LRI$Malaria,
  GBD_region = sim_df_LRI$GBD_region
)

U5_matrix <- readRDS("outputs/sensitivity_analysis/SA_matrix.result.ver.1.rds")

U5_matrix_meta <- U5_matrix[,1:2]
U5_matrix_LRI <- U5_matrix[,3:(n_sim+2)] + U5_matrix[,(2*n_sim+3):(3*n_sim+2)]

U5_matrix_LRI_index <- as.vector(as.matrix(U5_matrix_LRI))
U5_Year_index <- rep(U5_matrix_meta$Year, times = n_sim)
U5_Country_index <- rep(U5_matrix_meta$Country, times = n_sim)


sim_incidence <- tibble(
  Country = U5_Country_index,
  Year = U5_Year_index,
  Incidence_LRI = U5_matrix_LRI_index,
  sim = rep(1:n_sim, each = nrow(U5_matrix_meta))
)


Country_Year_Incidence <- tibble(
  Country = sim_incidence$Country,
  Year = sim_incidence$Year,
  sim = sim_incidence$sim
)

a <- unique(Country_Year_Antibiotics$Country)
b <- unique(Country_Year_Incidence$Country)
setdiff(a,b)
setdiff(b,a)

country_mapping <- c(
  "Bolivia" = "Bolivia (Plurinational State of)",
  "Congo, Dem. Rep." = "Democratic Republic of the Congo",
  "Congo, Rep." = "Congo",
  "Cote d'Ivoire" = "Côte d'Ivoire",
  "Egypt, Arab Rep." = "Egypt",
  "Gambia, The" = "Gambia",
  "Iran, Islamic Rep." = "Iran (Islamic Republic of)",
  "Kyrgyz Republic" = "Kyrgyzstan",
  "Lao PDR" = "Lao People's Democratic Republic",
  "Moldova" = "Republic of Moldova",
  "St. Vincent and the Grenadines" = "Saint Vincent and the Grenadines",
  "Tanzania" = "United Republic of Tanzania",
  "Turkiye" = "Türkiye",
  "Venezuela, RB" = "Venezuela (Bolivarian Republic of)",
  "West Bank and Gaza" = "Palestine",
  "Yemen, Rep." = "Yemen" 
)

Country_Year_Antibiotics <- Country_Year_Antibiotics %>% mutate(Country = ifelse(Country %in% names(country_mapping), country_mapping[Country], Country))
sim_df_LRI <- sim_df_LRI %>% mutate(Country = ifelse(Country %in% names(country_mapping), country_mapping[Country], Country))

Country_Year_meta <- inner_join(Country_Year_Incidence, Country_Year_Antibiotics, by = c("Country", "Year", "sim"))

rm(U5_matrix, U5_matrix_LRI, U5_matrix_meta, 
   UNI_indexed, UNI_indexed_Abaccess, UNI_indexed_Abdata, UNI_indexed_both,
   LRI_ratio_single_Ab, LRI_Ab_probability.matrix,
   U5_Country_index, U5_matrix_LRI_index, U5_Year_index, 
   Year_index, Malaria_index, Income_index, GBD_index, GDP_index, Country_index, AbLRI_long)
gc()



###### 2 ######

Antibiotics_LRI <- Country_Year_meta %>% inner_join(sim_df_LRI, by = c("Country","Year","sim","GDP","Income","Malaria","GBD_region"))
rm(sim_df_LRI)

colnames(sim_incidence)
colnames(Antibiotics_LRI)

rm(Country_Year_Incidence, Country_Year_Antibiotics)
gc()



###### 3 ########

treated_LRI_matrix <- tibble(matrix(NA_real_, nrow = nrow(sim_incidence), ncol = Ab_num))
for (i in 1:Ab_num) {
  treated_LRI_matrix[,i] <- Antibiotics_LRI[, 7 + i] * sim_incidence$Incidence_LRI
}

proportion_LRI_matrix <- tibble(matrix(NA_real_, nrow = nrow(sim_incidence), ncol = Ab_num))
for (i in 1:Ab_num) {
  proportion_LRI_matrix[,i] <- treated_LRI_matrix[,i] / sim_incidence$Incidence_LRI
}

colnames(treated_LRI_matrix) <- gsub("AbLRI_", "Treated_LRI_", colnames(Antibiotics_LRI[, 8:11]))
colnames(proportion_LRI_matrix) <- gsub("AbLRI_", "proportion_LRI_", colnames(Antibiotics_LRI[, 8:11]))

treated_LRI_matrix <- cbind(Country_Year_meta, treated_LRI_matrix)
proportion_LRI_matrix <- cbind(Country_Year_meta, proportion_LRI_matrix)

rm(Antibiotics_LRI)
gc()

saveRDS(Country_Year_meta,"outputs/sensitivity_analysis/SA_Country_Year_meta4.rds")
saveRDS(sim_incidence,"outputs/sensitivity_analysis/SA_sim_incidence4.rds")
saveRDS(treated_LRI_matrix,"outputs/sensitivity_analysis/SA_treated_LRI_matrix4.rds")
saveRDS(proportion_LRI_matrix, "outputs/sensitivity_analysis/SA_proportion_LRI_matrix4.rds")

rm(list = ls())
gc()

























##### 0 #####
rm(list = ls())
gc()

requiredPackages = c("dplyr", "tidyr", "tibble", "data.table", "stringr", "MASS")
lib.path <- .libPaths()
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p, lib = lib.path)
  library(p,character.only = TRUE)
}

##### 1 #####
n_sim <- 10000
set.seed(123)


Country_Year_meta <- readRDS("outputs/sensitivity_analysis/SA_Country_Year_meta4.rds")
sim_incidence <- readRDS("outputs/sensitivity_analysis/SA_sim_incidence4.rds")
treated_LRI_matrix <- readRDS("outputs/sensitivity_analysis/SA_treated_LRI_matrix4.rds")
proportion_LRI_matrix <- readRDS("outputs/sensitivity_analysis/SA_proportion_LRI_matrix4.rds")

digit_number <- c(3,4,3,3)

#----------------------------1.1 by country and year (treatment) ------------------------------


### LRI ###
treated_LRI_result <- treated_LRI_matrix %>%
  group_by(Country, Year, GDP, Income, Malaria, GBD_region) %>%
  summarise(
    Treated_LRI_BSP_med = median(Treated_LRI_BSP, na.rm = TRUE),
    Treated_LRI_BSP_lower = quantile(Treated_LRI_BSP, 0.025, na.rm = TRUE),
    Treated_LRI_BSP_upper = quantile(Treated_LRI_BSP, 0.975, na.rm = TRUE),
    
    Treated_LRI_NSP_med = median(Treated_LRI_NSP, na.rm = TRUE),
    Treated_LRI_NSP_lower = quantile(Treated_LRI_NSP, 0.025, na.rm = TRUE),
    Treated_LRI_NSP_upper = quantile(Treated_LRI_NSP, 0.975, na.rm = TRUE),
    
    Treated_LRI_Cephalosporins_med = median(Treated_LRI_Cephalosporins, na.rm = TRUE),
    Treated_LRI_Cephalosporins_lower = quantile(Treated_LRI_Cephalosporins, 0.025, na.rm = TRUE),
    Treated_LRI_Cephalosporins_upper = quantile(Treated_LRI_Cephalosporins, 0.975, na.rm = TRUE),
    
    Treated_LRI_Aminoglycosides_med = median(Treated_LRI_Aminoglycosides, na.rm = TRUE),
    Treated_LRI_Aminoglycosides_lower = quantile(Treated_LRI_Aminoglycosides, 0.025, na.rm = TRUE),
    Treated_LRI_Aminoglycosides_upper = quantile(Treated_LRI_Aminoglycosides, 0.975, na.rm = TRUE),
    
    .groups = "drop"
  )


abx_names <- c(
  "Treated_LRI_BSP", "Treated_LRI_NSP", "Treated_LRI_Cephalosporins", "Treated_LRI_Aminoglycosides")

for (i in 1:length(abx_names)) {
  med_col <- paste0(abx_names[i], "_med")
  lwr_col <- paste0(abx_names[i], "_lower")
  upr_col <- paste0(abx_names[i], "_upper")
  fmt_col <- paste0(abx_names[i])
  
  treated_LRI_result[[fmt_col]] <- sprintf(
    "%s (%s, %s)",
    formatC(treated_LRI_result[[med_col]], format = "f", big.mark = " ", digits = 0),
    formatC(treated_LRI_result[[lwr_col]], format = "f", big.mark = " ", digits = 0),
    formatC(treated_LRI_result[[upr_col]], format = "f", big.mark = " ", digits = 0)
  )
}

treated_LRI_result_summary <- treated_LRI_result[, c(1:6, 19:22)]
treated_LRI_result_df <- treated_LRI_result[, c(1:18)]

rm(treated_LRI_result)
gc()



proportion_LRI_result <- proportion_LRI_matrix %>%
  group_by(Country, Year, GDP, Income, Malaria, GBD_region) %>%
  summarise(
    proportion_LRI_BSP_med = median(proportion_LRI_BSP, na.rm = TRUE),
    proportion_LRI_BSP_lower = quantile(proportion_LRI_BSP, 0.025, na.rm = TRUE),
    proportion_LRI_BSP_upper = quantile(proportion_LRI_BSP, 0.975, na.rm = TRUE),
    
    proportion_LRI_NSP_med = median(proportion_LRI_NSP, na.rm = TRUE),
    proportion_LRI_NSP_lower = quantile(proportion_LRI_NSP, 0.025, na.rm = TRUE),
    proportion_LRI_NSP_upper = quantile(proportion_LRI_NSP, 0.975, na.rm = TRUE),
    
    proportion_LRI_Cephalosporins_med = median(proportion_LRI_Cephalosporins, na.rm = TRUE),
    proportion_LRI_Cephalosporins_lower = quantile(proportion_LRI_Cephalosporins, 0.025, na.rm = TRUE),
    proportion_LRI_Cephalosporins_upper = quantile(proportion_LRI_Cephalosporins, 0.975, na.rm = TRUE),
    
    proportion_LRI_Aminoglycosides_med = median(proportion_LRI_Aminoglycosides, na.rm = TRUE),
    proportion_LRI_Aminoglycosides_lower = quantile(proportion_LRI_Aminoglycosides, 0.025, na.rm = TRUE),
    proportion_LRI_Aminoglycosides_upper = quantile(proportion_LRI_Aminoglycosides, 0.975, na.rm = TRUE),
    
    .groups = "drop"
  )


abx_names <- c(
  "proportion_LRI_BSP", "proportion_LRI_NSP", "proportion_LRI_Cephalosporins", "proportion_LRI_Aminoglycosides"
)

for (i in 1:length(abx_names)) {
  med_col <- paste0(abx_names[i], "_med")
  lwr_col <- paste0(abx_names[i], "_lower")
  upr_col <- paste0(abx_names[i], "_upper")
  fmt_col <- paste0(abx_names[i])
  
  proportion_LRI_result[[fmt_col]] <- sprintf(
    "%s (%s, %s)",
    formatC(proportion_LRI_result[[med_col]], format = "f", big.mark = " ", digit_number[i]),
    formatC(proportion_LRI_result[[lwr_col]], format = "f", big.mark = " ", digit_number[i]),
    formatC(proportion_LRI_result[[upr_col]], format = "f", big.mark = " ", digit_number[i])
  )
}

proportion_LRI_result_summary <- proportion_LRI_result[, c(1:6, 19:22)]
proportion_LRI_result_df <- proportion_LRI_result[, c(1:18)]

rm(proportion_LRI_result)

gc()



LRI_result_summary <- cbind(treated_LRI_result_summary, proportion_LRI_result_summary[7:10])

years_to_save <- 2003:2019

#### 1.1 CSV file ####
for (yr in years_to_save) {
  write.csv(LRI_result_summary[LRI_result_summary$Year == yr, ], paste0("outputs/table/by_class_country_year/probability/SA_LRI_", yr, "4.csv"), row.names = FALSE)
}

rm(LRI_result_summary, treated_LRI_result, treated_LRI_result_df, treated_LRI_result_summary, 
   proportion_LRI_matrix, proportion_LRI_result, proportion_LRI_result_df, proportion_LRI_result_summary)

gc()






treated_LRI_matrix <- cbind(treated_LRI_matrix, sim_incidence[,3])




#----------------------------1.2 by year (treatment) ------------------------------

#---LRI---
treated_LRI_sim <- treated_LRI_matrix %>% 
  group_by(Year, sim) %>%
  summarise(
    Treated_LRI_BSP = sum(Treated_LRI_BSP, na.rm = TRUE),
    Treated_LRI_NSP = sum(Treated_LRI_NSP, na.rm = TRUE),
    Treated_LRI_Cephalosporins = sum(Treated_LRI_Cephalosporins, na.rm = TRUE),
    Treated_LRI_Aminoglycosides = sum(Treated_LRI_Aminoglycosides, na.rm = TRUE),
    Incidence = sum(Incidence_LRI, na.rm = TRUE),
    .groups = "drop"
  )

treated_LRI_result <- treated_LRI_sim %>%
  group_by(Year) %>%
  summarise(
    
    proportion_LRI_BSP_med = median(Treated_LRI_BSP/Incidence, na.rm = TRUE),
    proportion_LRI_BSP_lower = quantile(Treated_LRI_BSP/Incidence, 0.025, na.rm = TRUE),
    proportion_LRI_BSP_upper = quantile(Treated_LRI_BSP/Incidence, 0.975, na.rm = TRUE),
    
    proportion_LRI_NSP_med = median(Treated_LRI_NSP/Incidence, na.rm = TRUE),
    proportion_LRI_NSP_lower = quantile(Treated_LRI_NSP/Incidence, 0.025, na.rm = TRUE),
    proportion_LRI_NSP_upper = quantile(Treated_LRI_NSP/Incidence, 0.975, na.rm = TRUE),
    
    proportion_LRI_Cephalosporins_med = median(Treated_LRI_Cephalosporins/Incidence, na.rm = TRUE),
    proportion_LRI_Cephalosporins_lower = quantile(Treated_LRI_Cephalosporins/Incidence, 0.025, na.rm = TRUE),
    proportion_LRI_Cephalosporins_upper = quantile(Treated_LRI_Cephalosporins/Incidence, 0.975, na.rm = TRUE),
    
    proportion_LRI_Aminoglycosides_med = median(Treated_LRI_Aminoglycosides/Incidence, na.rm = TRUE),
    proportion_LRI_Aminoglycosides_lower = quantile(Treated_LRI_Aminoglycosides/Incidence, 0.025, na.rm = TRUE),
    proportion_LRI_Aminoglycosides_upper = quantile(Treated_LRI_Aminoglycosides/Incidence, 0.975, na.rm = TRUE),
    
    .groups = "drop"
  )


abx_names <- c(
  "proportion_LRI_BSP", "proportion_LRI_NSP", "proportion_LRI_Cephalosporins","proportion_LRI_Aminoglycosides"
)

for (i in 1:length(abx_names)) {
  med_col <- paste0(abx_names[i], "_med")
  lwr_col <- paste0(abx_names[i], "_lower")
  upr_col <- paste0(abx_names[i], "_upper")
  fmt_col <- paste0(abx_names[i])
  
  treated_LRI_result[[fmt_col]] <- sprintf(
    "%s (%s, %s)",
    formatC(treated_LRI_result[[med_col]], format = "f", big.mark = " ", digits = digit_number[i]),
    formatC(treated_LRI_result[[lwr_col]], format = "f", big.mark = " ", digits = digit_number[i]),
    formatC(treated_LRI_result[[upr_col]], format = "f", big.mark = " ", digits = digit_number[i])
  )
}

treated_LRI_summary <- treated_LRI_result[,c(1,14:17)]
treated_LRI_df <- treated_LRI_result[,c(1:13)]

#### 1.2 CSV file ####
write.csv(treated_LRI_summary, "outputs/table/by_class_year/probability/SA_LRI4.csv")
write.csv(treated_LRI_df, "outputs/by_class_df/by_year/probability/SA_LRI4.csv")


rm(treated_LRI_sim, treated_LRI_result, treated_LRI_summary, treated_LRI_df)
























#----------------------------1.3 by year and malaria (treatment) ------------------------------


#---LRI---
treated_LRI_sim <- treated_LRI_matrix %>% 
  group_by(Year, Malaria, sim) %>%
  summarise(
    Treated_LRI_BSP = sum(Treated_LRI_BSP, na.rm = TRUE),
    Treated_LRI_NSP = sum(Treated_LRI_NSP, na.rm = TRUE),
    Treated_LRI_Cephalosporins = sum(Treated_LRI_Cephalosporins, na.rm = TRUE),
    Treated_LRI_Aminoglycosides = sum(Treated_LRI_Aminoglycosides, na.rm = TRUE),
    Incidence = sum(Incidence_LRI, na.rm = TRUE),
    .groups = "drop"
  )

treated_LRI_result <- treated_LRI_sim %>%
  group_by(Year, Malaria) %>%
  summarise(
    
    proportion_LRI_BSP_med = median(Treated_LRI_BSP/Incidence, na.rm = TRUE),
    proportion_LRI_BSP_lower = quantile(Treated_LRI_BSP/Incidence, 0.025, na.rm = TRUE),
    proportion_LRI_BSP_upper = quantile(Treated_LRI_BSP/Incidence, 0.975, na.rm = TRUE),
    
    proportion_LRI_NSP_med = median(Treated_LRI_NSP/Incidence, na.rm = TRUE),
    proportion_LRI_NSP_lower = quantile(Treated_LRI_NSP/Incidence, 0.025, na.rm = TRUE),
    proportion_LRI_NSP_upper = quantile(Treated_LRI_NSP/Incidence, 0.975, na.rm = TRUE),
    
    proportion_LRI_Cephalosporins_med = median(Treated_LRI_Cephalosporins/Incidence, na.rm = TRUE),
    proportion_LRI_Cephalosporins_lower = quantile(Treated_LRI_Cephalosporins/Incidence, 0.025, na.rm = TRUE),
    proportion_LRI_Cephalosporins_upper = quantile(Treated_LRI_Cephalosporins/Incidence, 0.975, na.rm = TRUE),
    
    proportion_LRI_Aminoglycosides_med = median(Treated_LRI_Aminoglycosides/Incidence, na.rm = TRUE),
    proportion_LRI_Aminoglycosides_lower = quantile(Treated_LRI_Aminoglycosides/Incidence, 0.025, na.rm = TRUE),
    proportion_LRI_Aminoglycosides_upper = quantile(Treated_LRI_Aminoglycosides/Incidence, 0.975, na.rm = TRUE),
    
    .groups = "drop"
  )

abx_names <- c(
  "proportion_LRI_BSP", "proportion_LRI_NSP", "proportion_LRI_Cephalosporins", "proportion_LRI_Aminoglycosides")

for (i in 1:length(abx_names)) {
  med_col <- paste0(abx_names[i], "_med")
  lwr_col <- paste0(abx_names[i], "_lower")
  upr_col <- paste0(abx_names[i], "_upper")
  fmt_col <- paste0(abx_names[i])
  
  treated_LRI_result[[fmt_col]] <- sprintf(
    "%s (%s, %s)",
    formatC(treated_LRI_result[[med_col]], format = "f", big.mark = " ", digits = digit_number[i]),
    formatC(treated_LRI_result[[lwr_col]], format = "f", big.mark = " ", digits = digit_number[i]),
    formatC(treated_LRI_result[[upr_col]], format = "f", big.mark = " ", digits = digit_number[i])
  )
}

treated_LRI_summary <- treated_LRI_result[,c(1:2,15:18)]
treated_LRI_df <- treated_LRI_result[,c(1:14)]

#### 1.3 CSV file ####
write.csv(treated_LRI_summary, "outputs/table/by_class_malaria_year/probability/SA_LRI4.csv")

Malaria_name <- unique(treated_LRI_df$Malaria)

for (mal in Malaria_name) {
  write.csv(treated_LRI_summary[treated_LRI_summary$Malaria == mal, ], paste0("outputs/table/by_class_malaria_year/probability/SA_LRI_", mal, "4.csv"), row.names = FALSE)
}

write.csv(treated_LRI_df, "outputs/by_class_df/by_malaria_year/probability/SA_LRI_all4.csv")

for (mal in Malaria_name) {
  write.csv(treated_LRI_df[treated_LRI_df$Malaria == mal, ], paste0("outputs/by_class_df/by_malaria_year/probability/SA_LRI_", mal, "4.csv"), row.names = FALSE)
}

rm(treated_LRI_sim, treated_LRI_result, treated_LRI_summary, treated_LRI_df)




















#----------------------------1.4 by year and GBD region (treatment) ------------------------------


#---LRI---
treated_LRI_sim <- treated_LRI_matrix %>% 
  group_by(Year, GBD_region, sim) %>%
  summarise(
    Treated_LRI_BSP = sum(Treated_LRI_BSP, na.rm = TRUE),
    Treated_LRI_NSP = sum(Treated_LRI_NSP, na.rm = TRUE),
    Treated_LRI_Cephalosporins = sum(Treated_LRI_Cephalosporins, na.rm = TRUE),
    Treated_LRI_Aminoglycosides = sum(Treated_LRI_Aminoglycosides, na.rm = TRUE),
    Incidence = sum(Incidence_LRI, na.rm = TRUE),
    .groups = "drop"
  )

treated_LRI_result <- treated_LRI_sim %>%
  group_by(Year, GBD_region) %>%
  summarise(
    
    proportion_LRI_BSP_med = median(Treated_LRI_BSP/Incidence, na.rm = TRUE),
    proportion_LRI_BSP_lower = quantile(Treated_LRI_BSP/Incidence, 0.025, na.rm = TRUE),
    proportion_LRI_BSP_upper = quantile(Treated_LRI_BSP/Incidence, 0.975, na.rm = TRUE),
    
    proportion_LRI_NSP_med = median(Treated_LRI_NSP/Incidence, na.rm = TRUE),
    proportion_LRI_NSP_lower = quantile(Treated_LRI_NSP/Incidence, 0.025, na.rm = TRUE),
    proportion_LRI_NSP_upper = quantile(Treated_LRI_NSP/Incidence, 0.975, na.rm = TRUE),
    
    proportion_LRI_Cephalosporins_med = median(Treated_LRI_Cephalosporins/Incidence, na.rm = TRUE),
    proportion_LRI_Cephalosporins_lower = quantile(Treated_LRI_Cephalosporins/Incidence, 0.025, na.rm = TRUE),
    proportion_LRI_Cephalosporins_upper = quantile(Treated_LRI_Cephalosporins/Incidence, 0.975, na.rm = TRUE),
    
    proportion_LRI_Aminoglycosides_med = median(Treated_LRI_Aminoglycosides/Incidence, na.rm = TRUE),
    proportion_LRI_Aminoglycosides_lower = quantile(Treated_LRI_Aminoglycosides/Incidence, 0.025, na.rm = TRUE),
    proportion_LRI_Aminoglycosides_upper = quantile(Treated_LRI_Aminoglycosides/Incidence, 0.975, na.rm = TRUE),
    
    .groups = "drop"
  )


abx_names <- c(
  "proportion_LRI_BSP", "proportion_LRI_NSP", "proportion_LRI_Cephalosporins", "proportion_LRI_Aminoglycosides")

for (i in 1:length(abx_names)) {
  med_col <- paste0(abx_names[i], "_med")
  lwr_col <- paste0(abx_names[i], "_lower")
  upr_col <- paste0(abx_names[i], "_upper")
  fmt_col <- paste0(abx_names[i])
  
  treated_LRI_result[[fmt_col]] <- sprintf(
    "%s (%s, %s)",
    formatC(treated_LRI_result[[med_col]], format = "f", big.mark = " ", digits = digit_number[i]),
    formatC(treated_LRI_result[[lwr_col]], format = "f", big.mark = " ", digits = digit_number[i]),
    formatC(treated_LRI_result[[upr_col]], format = "f", big.mark = " ", digits = digit_number[i])
  )
}

treated_LRI_summary <- treated_LRI_result[,c(1:2,15:18)]
treated_LRI_df <- treated_LRI_result[,c(1:14)]



#### 1.4 CSV file ####
write.csv(treated_LRI_summary, "outputs/table/by_class_gbd_region_year/probability/SA_LRI4.csv")

write.csv(treated_LRI_df, "outputs/by_class_df/by_gbd_region_year/probability/SA_LRI_all_Region4.csv")

GBD_region_name <- unique(treated_LRI_df$GBD_region)
for (gbd in GBD_region_name) {
  write.csv(treated_LRI_summary[treated_LRI_summary$GBD_region == gbd, ], paste0("outputs/table/by_class_gbd_region_year/probability/SA_LRI_", gbd, "4.csv"), row.names = FALSE)
}

GBD_region_name <- unique(treated_LRI_df$GBD_region)
for (gbd in GBD_region_name) {
  write.csv(treated_LRI_df[treated_LRI_df$GBD_region == gbd, ], paste0("outputs/by_class_df/by_gbd_region_year/probability/SA_LRI_", gbd, "4.csv"), row.names = FALSE)
}





rm(treated_LRI_sim, treated_LRI_result, treated_LRI_summary, treated_LRI_df)














#----------------------------1.5 by year and Income (treatment-LRI) ------------------------------


first_sim <- treated_LRI_matrix[1:(nrow(treated_LRI_matrix)/10000),]
first_sim <- first_sim[,c(1,2,5)]
income_dist <- pivot_wider(first_sim, names_from = Year, values_from = Income)


first5_years <- as.character(2003:2007)
last5_years  <- as.character(2015:2019)

income_dist_summary <- income_dist %>%
  rowwise() %>%
  mutate(
    First5 = {
      first_vals <- c_across(all_of(first5_years))
      tab <- table(first_vals)
      if(length(tab) == 0) NA_character_ else names(sort(tab, decreasing = TRUE))[1]
    },
    Last5 = {
      last_vals <- c_across(all_of(last5_years))
      tab <- table(last_vals)
      if(length(tab) == 0) NA_character_ else names(sort(tab, decreasing = TRUE))[1]
    }
  ) %>%
  ungroup()

income_dist_summary[which_na(income_dist_summary$First5),]
income_dist_summary[income_dist_summary$Country == "Djibouti",]$First5 <- "LM" 
treated_LRI_matrix <- treated_LRI_matrix %>% inner_join(income_dist_summary[,c("Country","First5","Last5")], by = "Country") 
rm(income_dist, income_dist_summary, first_sim, first5_years, last5_years)




#---LRI---
treated_LRI_sim <- treated_LRI_matrix %>% 
  group_by(Year, First5, sim) %>%
  summarise(
    Treated_LRI_BSP = sum(Treated_LRI_BSP, na.rm = TRUE),
    Treated_LRI_NSP = sum(Treated_LRI_NSP, na.rm = TRUE),
    Treated_LRI_Cephalosporins = sum(Treated_LRI_Cephalosporins, na.rm = TRUE),
    Treated_LRI_Aminoglycosides = sum(Treated_LRI_Aminoglycosides, na.rm = TRUE),
    Incidence = sum(Incidence_LRI, na.rm = TRUE),
    .groups = "drop"
  )

treated_LRI_result <- treated_LRI_sim %>%
  group_by(Year, First5) %>%
  summarise(
    
    proportion_LRI_BSP_med = median(Treated_LRI_BSP/Incidence, na.rm = TRUE),
    proportion_LRI_BSP_lower = quantile(Treated_LRI_BSP/Incidence, 0.025, na.rm = TRUE),
    proportion_LRI_BSP_upper = quantile(Treated_LRI_BSP/Incidence, 0.975, na.rm = TRUE),
    
    proportion_LRI_NSP_med = median(Treated_LRI_NSP/Incidence, na.rm = TRUE),
    proportion_LRI_NSP_lower = quantile(Treated_LRI_NSP/Incidence, 0.025, na.rm = TRUE),
    proportion_LRI_NSP_upper = quantile(Treated_LRI_NSP/Incidence, 0.975, na.rm = TRUE),
    
    proportion_LRI_Cephalosporins_med = median(Treated_LRI_Cephalosporins/Incidence, na.rm = TRUE),
    proportion_LRI_Cephalosporins_lower = quantile(Treated_LRI_Cephalosporins/Incidence, 0.025, na.rm = TRUE),
    proportion_LRI_Cephalosporins_upper = quantile(Treated_LRI_Cephalosporins/Incidence, 0.975, na.rm = TRUE),
    
    proportion_LRI_Aminoglycosides_med = median(Treated_LRI_Aminoglycosides/Incidence, na.rm = TRUE),
    proportion_LRI_Aminoglycosides_lower = quantile(Treated_LRI_Aminoglycosides/Incidence, 0.025, na.rm = TRUE),
    proportion_LRI_Aminoglycosides_upper = quantile(Treated_LRI_Aminoglycosides/Incidence, 0.975, na.rm = TRUE),
    
    .groups = "drop"
  )



abx_names <- c(
  "proportion_LRI_BSP", "proportion_LRI_NSP", "proportion_LRI_Cephalosporins", "proportion_LRI_Aminoglycosides")

for (i in 1:length(abx_names)) {
  med_col <- paste0(abx_names[i], "_med")
  lwr_col <- paste0(abx_names[i], "_lower")
  upr_col <- paste0(abx_names[i], "_upper")
  fmt_col <- paste0(abx_names[i])
  
  treated_LRI_result[[fmt_col]] <- sprintf(
    "%s (%s, %s)",
    formatC(treated_LRI_result[[med_col]], format = "f", big.mark = " ", digits = digit_number[i]),
    formatC(treated_LRI_result[[lwr_col]], format = "f", big.mark = " ", digits = digit_number[i]),
    formatC(treated_LRI_result[[upr_col]], format = "f", big.mark = " ", digits = digit_number[i])
  )
}

treated_LRI_summary_first5 <- treated_LRI_result[,c(1:2,15:18)]
treated_LRI_df_first5 <- treated_LRI_result[,c(1:14)]



#### 1.5.1 CSV file ####
write.csv(treated_LRI_summary_first5, "outputs/table/by_class_income_year/probability/SA_LRI4_first5.csv")
write.csv(treated_LRI_df_first5, "outputs/by_class_df/by_income_year/probability/SA_LRI4_first5.csv")

Income <- unique(treated_LRI_df_first5$First5)
for (inc in Income) {
  write.csv(treated_LRI_summary_first5[treated_LRI_summary_first5$First5 == inc, ], paste0("outputs/table/by_class_income_year/probability/SA_LRI4_first5_", inc, ".csv"), row.names = FALSE)
}

Income <- unique(treated_LRI_df_first5$First5)
for (inc in Income) {
  write.csv(treated_LRI_df_first5[treated_LRI_df_first5$First5 == inc, ], paste0("outputs/by_class_df/by_income_year/probability/SA_LRI4_first5_", inc, ".csv"), row.names = FALSE)
}

rm(treated_LRI_result, treated_LRI_df_first5, treated_LRI_summary_first5, treated_LRI_sim)









#---LRI---
treated_LRI_sim <- treated_LRI_matrix %>% 
  group_by(Year, Last5, sim) %>%
  summarise(
    Treated_LRI_BSP = sum(Treated_LRI_BSP, na.rm = TRUE),
    Treated_LRI_NSP = sum(Treated_LRI_NSP, na.rm = TRUE),
    Treated_LRI_Cephalosporins = sum(Treated_LRI_Cephalosporins, na.rm = TRUE),
    Treated_LRI_Aminoglycosides = sum(Treated_LRI_Aminoglycosides, na.rm = TRUE),
    Incidence = sum(Incidence_LRI, na.rm = TRUE),
    .groups = "drop"
  )

treated_LRI_result <- treated_LRI_sim %>%
  group_by(Year, Last5) %>%
  summarise(
    
    proportion_LRI_BSP_med = median(Treated_LRI_BSP/Incidence, na.rm = TRUE),
    proportion_LRI_BSP_lower = quantile(Treated_LRI_BSP/Incidence, 0.025, na.rm = TRUE),
    proportion_LRI_BSP_upper = quantile(Treated_LRI_BSP/Incidence, 0.975, na.rm = TRUE),
    
    proportion_LRI_NSP_med = median(Treated_LRI_NSP/Incidence, na.rm = TRUE),
    proportion_LRI_NSP_lower = quantile(Treated_LRI_NSP/Incidence, 0.025, na.rm = TRUE),
    proportion_LRI_NSP_upper = quantile(Treated_LRI_NSP/Incidence, 0.975, na.rm = TRUE),
    
    proportion_LRI_Cephalosporins_med = median(Treated_LRI_Cephalosporins/Incidence, na.rm = TRUE),
    proportion_LRI_Cephalosporins_lower = quantile(Treated_LRI_Cephalosporins/Incidence, 0.025, na.rm = TRUE),
    proportion_LRI_Cephalosporins_upper = quantile(Treated_LRI_Cephalosporins/Incidence, 0.975, na.rm = TRUE),
    
    proportion_LRI_Aminoglycosides_med = median(Treated_LRI_Aminoglycosides/Incidence, na.rm = TRUE),
    proportion_LRI_Aminoglycosides_lower = quantile(Treated_LRI_Aminoglycosides/Incidence, 0.025, na.rm = TRUE),
    proportion_LRI_Aminoglycosides_upper = quantile(Treated_LRI_Aminoglycosides/Incidence, 0.975, na.rm = TRUE),
    
    .groups = "drop"
  )


abx_names <- c(
  "proportion_LRI_BSP", "proportion_LRI_NSP", "proportion_LRI_Cephalosporins", "proportion_LRI_Aminoglycosides")

for (i in 1:length(abx_names)) {
  med_col <- paste0(abx_names[i], "_med")
  lwr_col <- paste0(abx_names[i], "_lower")
  upr_col <- paste0(abx_names[i], "_upper")
  fmt_col <- paste0(abx_names[i])
  
  treated_LRI_result[[fmt_col]] <- sprintf(
    "%s (%s, %s)",
    formatC(treated_LRI_result[[med_col]], format = "f", big.mark = " ", digits = digit_number[i]),
    formatC(treated_LRI_result[[lwr_col]], format = "f", big.mark = " ", digits = digit_number[i]),
    formatC(treated_LRI_result[[upr_col]], format = "f", big.mark = " ", digits = digit_number[i])
  )
}

treated_LRI_summary_last5 <- treated_LRI_result[,c(1:2,15:18)]
treated_LRI_df_last5 <- treated_LRI_result[,c(1:14)]




#### 1.5.2 CSV file ####
write.csv(treated_LRI_summary_last5, "outputs/table/by_class_income_year/probability/SA_LRI4_last5.csv")
write.csv(treated_LRI_df_last5, "outputs/by_class_df/by_income_year/probability/SA_LRI4_last5.csv")

Income <- unique(treated_LRI_df_last5$Last5)
for (inc in Income) {
  write.csv(treated_LRI_summary_last5[treated_LRI_summary_last5$Last5 == inc, ], paste0("outputs/table/by_class_income_year/probability/SA_LRI4_last5_", inc, ".csv"), row.names = FALSE)
}

Income <- unique(treated_LRI_df_last5$Last5)
for (inc in Income) {
  write.csv(treated_LRI_df_last5[treated_LRI_df_last5$Last5 == inc, ], paste0("outputs/by_class_df/by_income_year/probability/SA_LRI4_last5_", inc, ".csv"), row.names = FALSE)
}

rm(treated_LRI_result, treated_LRI_df_last5, treated_LRI_summary_last5, treated_LRI_sim)

#
















##### 0 #####
rm(list = ls())
gc()

# Optimized strictly for dataframe manipulation
requiredPackages = c("dplyr", "tidyr", "tibble")
lib.path <- .libPaths()
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p, lib = lib.path)
  library(p,character.only = TRUE)
}


##### 1 #####
n_sim <- 10000
set.seed(123)

# confirm upstream antibiotic-treatment model uses PC1-PC2
coef_df_LRI <- readRDS("outputs/sensitivity_analysis/SA_coef_df_LRI.rds")

stopifnot(
  identical(
    names(coef_df_LRI),
    c("(Intercept)", "PC1", "PC2")
  )
)

rm(coef_df_LRI)

# The following objects were generated using the current PC1-PC2 models
Country_Year_meta <- readRDS("outputs/sensitivity_analysis/SA_Country_Year_meta4.rds")
sim_incidence <- readRDS("outputs/sensitivity_analysis/SA_sim_incidence4.rds")
treated_LRI_matrix <- readRDS("outputs/sensitivity_analysis/SA_treated_LRI_matrix4.rds")
proportion_LRI_matrix <- readRDS("outputs/sensitivity_analysis/SA_proportion_LRI_matrix4.rds")

treated_LRI_matrix$Treated_all <- rowSums(treated_LRI_matrix[,8:11])
proportion_LRI_matrix$Proportion_LRI_all <- rowSums(proportion_LRI_matrix[,8:11])

head(treated_LRI_matrix)
unique(treated_LRI_matrix$Year)
unique(treated_LRI_matrix$sim)
digit_number <- c(3,4,3,3)


#----------------------------1.2 by year (treatment) ------------------------------
#==============================
# Incremental antibiotic share (LRI)
# Baseline year: 2003
#==============================

baseline_year <- 2003

#--- 1. Aggregate by Year × sim ---
prop_LRI_sim <- proportion_LRI_matrix %>%
  group_by(Year, sim) %>%
  summarise(
    BSP   = mean(proportion_LRI_BSP, na.rm = TRUE),
    NSP   = mean(proportion_LRI_NSP, na.rm = TRUE),
    CEPH  = mean(proportion_LRI_Cephalosporins, na.rm = TRUE),
    AMINO = mean(proportion_LRI_Aminoglycosides, na.rm = TRUE),
    ALL   = mean(Proportion_LRI_all, na.rm = TRUE),
    .groups = "drop"
  )

#--- 2. Extract baseline (2003) per simulation ---
baseline_2003 <- prop_LRI_sim %>%
  filter(Year == baseline_year) %>%
  dplyr::select(
    sim,
    BSP_2003   = BSP,
    NSP_2003   = NSP,
    CEPH_2003  = CEPH,
    AMINO_2003 = AMINO,
    ALL_2003   = ALL
  )

#--- 3. Compute incremental shares per simulation ---
prop_LRI_incremental <- prop_LRI_sim %>%
  left_join(baseline_2003, by = "sim") %>%
  mutate(
    d_ALL = ALL - ALL_2003
  ) %>%
  mutate(
    inc_BSP   = (BSP   - BSP_2003)   / d_ALL,
    inc_NSP   = (NSP   - NSP_2003)   / d_ALL,
    inc_CEPH  = (CEPH  - CEPH_2003)  / d_ALL,
    inc_AMINO = (AMINO - AMINO_2003) / d_ALL
  ) %>%
  # protect against division by ~0
  mutate(
    across(starts_with("inc_"),
           ~ ifelse(abs(d_ALL) < 1e-6, NA, .))
  )

#--- 4. Summarise across simulations by year ---
treated_LRI_result <- prop_LRI_incremental %>%
  group_by(Year) %>%
  summarise(
    inc_LRI_BSP_med   = median(inc_BSP, na.rm = TRUE),
    inc_LRI_BSP_lower = quantile(inc_BSP, 0.025, na.rm = TRUE),
    inc_LRI_BSP_upper = quantile(inc_BSP, 0.975, na.rm = TRUE),
    
    inc_LRI_NSP_med   = median(inc_NSP, na.rm = TRUE),
    inc_LRI_NSP_lower = quantile(inc_NSP, 0.025, na.rm = TRUE),
    inc_LRI_NSP_upper = quantile(inc_NSP, 0.975, na.rm = TRUE),
    
    inc_LRI_Cephalosporins_med   = median(inc_CEPH, na.rm = TRUE),
    inc_LRI_Cephalosporins_lower = quantile(inc_CEPH, 0.025, na.rm = TRUE),
    inc_LRI_Cephalosporins_upper = quantile(inc_CEPH, 0.975, na.rm = TRUE),
    
    inc_LRI_Aminoglycosides_med   = median(inc_AMINO, na.rm = TRUE),
    inc_LRI_Aminoglycosides_lower = quantile(inc_AMINO, 0.025, na.rm = TRUE),
    inc_LRI_Aminoglycosides_upper = quantile(inc_AMINO, 0.975, na.rm = TRUE),
    
    .groups = "drop"
  )

#--- 5. Format for table output ---
abx_names <- c(
  "inc_LRI_BSP",
  "inc_LRI_NSP",
  "inc_LRI_Cephalosporins",
  "inc_LRI_Aminoglycosides"
)

for (i in seq_along(abx_names)) {
  med_col <- paste0(abx_names[i], "_med")
  lwr_col <- paste0(abx_names[i], "_lower")
  upr_col <- paste0(abx_names[i], "_upper")
  fmt_col <- abx_names[i]
  
  treated_LRI_result[[fmt_col]] <- sprintf(
    "%s (%s, %s)",
    formatC(treated_LRI_result[[med_col]], format = "f",
            digits = digit_number[i]),
    formatC(treated_LRI_result[[lwr_col]], format = "f",
            digits = digit_number[i]),
    formatC(treated_LRI_result[[upr_col]], format = "f",
            digits = digit_number[i])
  )
}

#--- 6. Final outputs ---
treated_LRI_summary <- treated_LRI_result %>%
  filter(Year != baseline_year) %>%
  dplyr::select(Year, all_of(abx_names))

treated_LRI_df <- treated_LRI_result %>%
  filter(Year != baseline_year)

#--- 7. Write CSVs ---
write.csv(
  treated_LRI_summary,
  "outputs/table/by_class_year/probability/SA_LRI_incremental_share_prop.csv",
  row.names = FALSE
)

write.csv(
  treated_LRI_df,
  "outputs/by_class_df/by_year/probability/SA_LRI_incremental_share_df_prop.csv",
  row.names = FALSE
)

#--- 8. Clean up ---
rm(
  prop_LRI_sim,
  prop_LRI_incremental,
  baseline_2003,
  treated_LRI_result,
  treated_LRI_summary,
  treated_LRI_df
)



#----------------------------1.4 by year and GBD region (treatment) ------------------------------

baseline_year <- 2003

prop_LRI_sim <- proportion_LRI_matrix %>%
  group_by(Year, GBD_region, sim) %>%
  summarise(
    BSP   = mean(proportion_LRI_BSP, na.rm = TRUE),
    NSP   = mean(proportion_LRI_NSP, na.rm = TRUE),
    CEPH  = mean(proportion_LRI_Cephalosporins, na.rm = TRUE),
    AMINO = mean(proportion_LRI_Aminoglycosides, na.rm = TRUE),
    ALL   = mean(Proportion_LRI_all, na.rm = TRUE),
    .groups = "drop"
  )

baseline_2003 <- prop_LRI_sim %>%
  filter(Year == baseline_year) %>%
  dplyr::select(
    GBD_region, sim,
    BSP_2003   = BSP,
    NSP_2003   = NSP,
    CEPH_2003  = CEPH,
    AMINO_2003 = AMINO,
    ALL_2003   = ALL
  )

prop_LRI_incremental <- prop_LRI_sim %>%
  left_join(baseline_2003, by = c("GBD_region", "sim")) %>%
  mutate(
    d_ALL = ALL - ALL_2003
  ) %>%
  mutate(
    inc_BSP   = (BSP   - BSP_2003)   / d_ALL,
    inc_NSP   = (NSP   - NSP_2003)   / d_ALL,
    inc_CEPH  = (CEPH  - CEPH_2003)  / d_ALL,
    inc_AMINO = (AMINO - AMINO_2003) / d_ALL
  ) %>%
  mutate(
    across(starts_with("inc_"),
           ~ ifelse(abs(d_ALL) < 1e-6, NA, .))
  )

treated_LRI_result <- prop_LRI_incremental %>%
  group_by(Year, GBD_region) %>%
  summarise(
    inc_LRI_BSP_med   = median(inc_BSP, na.rm = TRUE),
    inc_LRI_BSP_lower = quantile(inc_BSP, 0.025, na.rm = TRUE),
    inc_LRI_BSP_upper = quantile(inc_BSP, 0.975, na.rm = TRUE),
    
    inc_LRI_NSP_med   = median(inc_NSP, na.rm = TRUE),
    inc_LRI_NSP_lower = quantile(inc_NSP, 0.025, na.rm = TRUE),
    inc_LRI_NSP_upper = quantile(inc_NSP, 0.975, na.rm = TRUE),
    
    inc_LRI_Cephalosporins_med   = median(inc_CEPH, na.rm = TRUE),
    inc_LRI_Cephalosporins_lower = quantile(inc_CEPH, 0.025, na.rm = TRUE),
    inc_LRI_Cephalosporins_upper = quantile(inc_CEPH, 0.975, na.rm = TRUE),
    
    inc_LRI_Aminoglycosides_med   = median(inc_AMINO, na.rm = TRUE),
    inc_LRI_Aminoglycosides_lower = quantile(inc_AMINO, 0.025, na.rm = TRUE),
    inc_LRI_Aminoglycosides_upper = quantile(inc_AMINO, 0.975, na.rm = TRUE),
    
    .groups = "drop"
  )

abx_names <- c(
  "inc_LRI_BSP",
  "inc_LRI_NSP",
  "inc_LRI_Cephalosporins",
  "inc_LRI_Aminoglycosides"
)

for (i in seq_along(abx_names)) {
  med_col <- paste0(abx_names[i], "_med")
  lwr_col <- paste0(abx_names[i], "_lower")
  upr_col <- paste0(abx_names[i], "_upper")
  fmt_col <- abx_names[i]
  
  treated_LRI_result[[fmt_col]] <- sprintf(
    "%s (%s, %s)",
    formatC(treated_LRI_result[[med_col]], format = "f",
            digits = digit_number[i]),
    formatC(treated_LRI_result[[lwr_col]], format = "f",
            digits = digit_number[i]),
    formatC(treated_LRI_result[[upr_col]], format = "f",
            digits = digit_number[i])
  )
}

treated_LRI_summary <- treated_LRI_result %>%
  filter(Year != baseline_year) %>%
  dplyr::select(Year, GBD_region, all_of(abx_names))

treated_LRI_df <- treated_LRI_result %>%
  filter(Year != baseline_year)


#### 1.4 CSV file ####
write.csv(
  treated_LRI_summary,
  "outputs/table/by_class_gbd_region_year/probability/SA_LRI_incremental_prop.csv",
  row.names = FALSE
)

write.csv(
  treated_LRI_df,
  "outputs/by_class_df/by_gbd_region_year/probability/SA_LRI_incremental_df_prop.csv",
  row.names = FALSE
)

GBD_region_name <- unique(treated_LRI_df$GBD_region)
for (gbd in GBD_region_name) {
  write.csv(treated_LRI_summary[treated_LRI_summary$GBD_region == gbd, ], paste0("outputs/table/by_class_gbd_region_year/probability/SA_LRI_prop_", gbd, "4.csv"), row.names = FALSE)
}

GBD_region_name <- unique(treated_LRI_df$GBD_region)
for (gbd in GBD_region_name) {
  write.csv(treated_LRI_df[treated_LRI_df$GBD_region == gbd, ], paste0("outputs/by_class_df/by_gbd_region_year/probability/SA_LRI_prop_", gbd, "4.csv"), row.names = FALSE)
}

rm(
  prop_LRI_sim,
  prop_LRI_incremental,
  baseline_2003,
  treated_LRI_result,
  treated_LRI_summary,
  treated_LRI_df
)



#----------------------------1.5 by year and Income (treatment-LRI) ------------------------------

first_sim <- treated_LRI_matrix[1:(nrow(treated_LRI_matrix)/10000),]
first_sim <- first_sim[,c(1,2,5)]
income_dist <- pivot_wider(first_sim, names_from = Year, values_from = Income)

first5_years <- as.character(2003:2007)
last5_years  <- as.character(2015:2019)

income_dist_summary <- income_dist %>%
  rowwise() %>%
  mutate(
    First5 = {
      first_vals <- c_across(all_of(first5_years))
      tab <- table(first_vals)
      if(length(tab) == 0) NA_character_ else names(sort(tab, decreasing = TRUE))[1]
    },
    Last5 = {
      last_vals <- c_across(all_of(last5_years))
      tab <- table(last_vals)
      if(length(tab) == 0) NA_character_ else names(sort(tab, decreasing = TRUE))[1]
    }
  ) %>%
  ungroup()

income_dist_summary[which_na(income_dist_summary$First5),]
income_dist_summary[income_dist_summary$Country == "Djibouti",]$First5 <- "LM" 
treated_LRI_matrix <- treated_LRI_matrix %>% inner_join(income_dist_summary[,c("Country","First5","Last5")], by = "Country") 
rm(income_dist, income_dist_summary, first_sim, first5_years, last5_years)



#---LRI---
baseline_year <- 2003

prop_LRI_sim <- proportion_LRI_matrix %>%
  inner_join(
    treated_LRI_matrix %>% distinct(Country, Last5),
    by = "Country"
  ) %>%
  group_by(Year, Last5, sim) %>%
  summarise(
    BSP   = mean(proportion_LRI_BSP, na.rm = TRUE),
    NSP   = mean(proportion_LRI_NSP, na.rm = TRUE),
    CEPH  = mean(proportion_LRI_Cephalosporins, na.rm = TRUE),
    AMINO = mean(proportion_LRI_Aminoglycosides, na.rm = TRUE),
    ALL   = mean(Proportion_LRI_all, na.rm = TRUE),
    .groups = "drop"
  )

baseline_2003 <- prop_LRI_sim %>%
  filter(Year == baseline_year) %>%
  dplyr::select(
    Last5, sim,
    BSP_2003   = BSP,
    NSP_2003   = NSP,
    CEPH_2003  = CEPH,
    AMINO_2003 = AMINO,
    ALL_2003   = ALL
  )

prop_LRI_incremental <- prop_LRI_sim %>%
  left_join(baseline_2003, by = c("Last5", "sim")) %>%
  mutate(
    d_ALL = ALL - ALL_2003,
    inc_BSP   = (BSP   - BSP_2003)   / d_ALL,
    inc_NSP   = (NSP   - NSP_2003)   / d_ALL,
    inc_CEPH  = (CEPH  - CEPH_2003)  / d_ALL,
    inc_AMINO = (AMINO - AMINO_2003) / d_ALL
  ) %>%
  mutate(
    across(starts_with("inc_"),
           ~ ifelse(abs(d_ALL) < 1e-6, NA, .))
  )

treated_LRI_result <- prop_LRI_incremental %>%
  group_by(Year, Last5) %>%
  summarise(
    inc_LRI_BSP_med   = median(inc_BSP, na.rm = TRUE),
    inc_LRI_BSP_lower = quantile(inc_BSP, 0.025, na.rm = TRUE),
    inc_LRI_BSP_upper = quantile(inc_BSP, 0.975, na.rm = TRUE),
    
    inc_LRI_NSP_med   = median(inc_NSP, na.rm = TRUE),
    inc_LRI_NSP_lower = quantile(inc_NSP, 0.025, na.rm = TRUE),
    inc_LRI_NSP_upper = quantile(inc_NSP, 0.975, na.rm = TRUE),
    
    inc_LRI_Cephalosporins_med   = median(inc_CEPH, na.rm = TRUE),
    inc_LRI_Cephalosporins_lower = quantile(inc_CEPH, 0.025, na.rm = TRUE),
    inc_LRI_Cephalosporins_upper = quantile(inc_CEPH, 0.975, na.rm = TRUE),
    
    inc_LRI_Aminoglycosides_med   = median(inc_AMINO, na.rm = TRUE),
    inc_LRI_Aminoglycosides_lower = quantile(inc_AMINO, 0.025, na.rm = TRUE),
    inc_LRI_Aminoglycosides_upper = quantile(inc_AMINO, 0.975, na.rm = TRUE),
    
    .groups = "drop"
  )

abx_names <- c(
  "inc_LRI_BSP",
  "inc_LRI_NSP",
  "inc_LRI_Cephalosporins",
  "inc_LRI_Aminoglycosides"
)

for (i in seq_along(abx_names)) {
  med_col <- paste0(abx_names[i], "_med")
  lwr_col <- paste0(abx_names[i], "_lower")
  upr_col <- paste0(abx_names[i], "_upper")
  fmt_col <- abx_names[i]
  
  treated_LRI_result[[fmt_col]] <- sprintf(
    "%s (%s, %s)",
    formatC(treated_LRI_result[[med_col]], format = "f",
            digits = digit_number[i]),
    formatC(treated_LRI_result[[lwr_col]], format = "f",
            digits = digit_number[i]),
    formatC(treated_LRI_result[[upr_col]], format = "f",
            digits = digit_number[i])
  )
}

treated_LRI_summary_last5 <- treated_LRI_result %>%
  filter(Year != baseline_year) %>%
  dplyr::select(Year, Last5, all_of(abx_names))

treated_LRI_df_last5 <- treated_LRI_result %>%
  filter(Year != baseline_year)


#### 1.5.2 CSV file ####
write.csv(
  treated_LRI_summary_last5,
  "outputs/table/by_class_income_year/probability/SA_LRI_incremental_last5_prop.csv",
  row.names = FALSE
)

write.csv(
  treated_LRI_df_last5,
  "outputs/by_class_df/by_income_year/probability/SA_LRI_incremental_last5_df_prop.csv",
  row.names = FALSE
)

Income <- unique(treated_LRI_df_last5$Last5)
for (inc in Income) {
  write.csv(treated_LRI_summary_last5[treated_LRI_summary_last5$Last5 == inc, ], paste0("outputs/table/by_class_income_year/probability/SA_LRI4_last5_prop_", inc, ".csv"), row.names = FALSE)
}

Income <- unique(treated_LRI_df_last5$Last5)
for (inc in Income) {
  write.csv(treated_LRI_df_last5[treated_LRI_df_last5$Last5 == inc, ], paste0("outputs/by_class_df/by_income_year/probability/SA_LRI4_last5_prop_", inc, ".csv"), row.names = FALSE)
}

rm(
  prop_LRI_sim,
  prop_LRI_incremental,
  baseline_2003,
  treated_LRI_result,
  treated_LRI_df_last5,
  treated_LRI_summary_last5
)



#--- 2. Extract baseline (2003) per simulation ---
treated_LRI_2019 <- treated_LRI_matrix %>% filter(Year == 2019) 
treated_LRI_2019$proportion_NSP_Amino <- (treated_LRI_2019$Treated_LRI_NSP + treated_LRI_2019$Treated_LRI_Aminoglycosides) / treated_LRI_2019$Treated_all

quantile(treated_LRI_2019$proportion_NSP_Amino, 0.025, na.rm = TRUE)
quantile(treated_LRI_2019$proportion_NSP_Amino, 0.500, na.rm = TRUE)
quantile(treated_LRI_2019$proportion_NSP_Amino, 0.975, na.rm = TRUE)


#### 0. prep ####
rm(list = ls())
gc()

requiredPackages = c("dplyr")
lib.path <- .libPaths()
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p, lib = lib.path)
  library(p,character.only = TRUE)
}

options(scipen = 100, digits = 6) # 지수표시 안되게
n_sim <- 10000
set.seed(123)


# confirm upstream models use PC1-PC2

coef_df_LRI <- readRDS("outputs/sensitivity_analysis/SA_coef_df_LRI.rds")
coef_df_LRI_DDD_by_4class <- readRDS("outputs/sensitivity_analysis/SA_coef_df_LRI_DDD_by_4class.rds")

stopifnot(
  identical(
    names(coef_df_LRI),
    c("(Intercept)", "PC1", "PC2")
  ),
  ncol(coef_df_LRI_DDD_by_4class) == 12,
  !any(grepl("PC3", names(coef_df_LRI_DDD_by_4class)))
)

rm(
  coef_df_LRI,
  coef_df_LRI_DDD_by_4class
)


# plot data below were generated from the current PC1-PC2 models

gbd_region <- c("North Africa and Middle East","Central Europe, Eastern Europe, and Central Asia","Sub-Saharan Africa","South Asia","Latin America and Caribbean","Southeast Asia, East Asia, and Oceania")

for (i in 1:6) { # GBD_region
  path_LRI <- paste0("outputs/by_class_df/by_gbd_region_year/probability/SA_LRI_", gbd_region[i], "4.csv")
  obj_LRI   <- paste0("LRI_year_proportion_GBD_region_", i)
  assign(obj_LRI,   read.csv(path_LRI,   stringsAsFactors = FALSE))
}

for (i in 1:6) { # GBD_region
  path_LRI_prop <- paste0("outputs/by_class_df/by_gbd_region_year/probability/SA_LRI_prop_", gbd_region[i], "4.csv")
  obj_LRI_prop   <- paste0("LRI_year_proportion_GBD_region_prop_", i)
  assign(obj_LRI_prop,   read.csv(path_LRI_prop,   stringsAsFactors = FALSE))
}

Income <- c("L", "LM", "UM")

for (inc in Income) { # Income
  path_LRI_last   <- paste0("outputs/by_class_df/by_income_year/probability/SA_LRI4_last5_", inc, ".csv")
  obj_LRI_last    <- paste0("LRI_year_proportion_income_last_", inc)
  assign(obj_LRI_last,   read.csv(path_LRI_last,   stringsAsFactors = FALSE))
}

for (inc in Income) { # Income
  path_LRI_last_prop   <- paste0("outputs/by_class_df/by_income_year/probability/SA_LRI4_last5_prop_", inc, ".csv")
  obj_LRI_last_prop    <- paste0("LRI_year_proportion_income_last_prop_", inc)
  assign(obj_LRI_last_prop,   read.csv(path_LRI_last_prop,   stringsAsFactors = FALSE))
}

abx_names <- c(
  "BSP", "NSP", "Cephalosporins", "Aminoglycosides")

LRI_year_proportion_prop <- read.csv("outputs/by_class_df/by_year/probability/SA_LRI_incremental_share_df_prop.csv")

LRI_incr_2019       <- LRI_year_proportion_prop               [LRI_year_proportion_prop               $Year == 2019,]
LRI_incr_2019_L     <- LRI_year_proportion_income_last_prop_L [LRI_year_proportion_income_last_prop_L $Year == 2019,]
LRI_incr_2019_LM    <- LRI_year_proportion_income_last_prop_LM[LRI_year_proportion_income_last_prop_LM$Year == 2019,]
LRI_incr_2019_UM    <- LRI_year_proportion_income_last_prop_UM[LRI_year_proportion_income_last_prop_UM$Year == 2019,]
LRI_incr_2019_GBD_1 <- LRI_year_proportion_GBD_region_prop_1  [LRI_year_proportion_GBD_region_prop_1  $Year == 2019,]
LRI_incr_2019_GBD_2 <- LRI_year_proportion_GBD_region_prop_2  [LRI_year_proportion_GBD_region_prop_2  $Year == 2019,]
LRI_incr_2019_GBD_3 <- LRI_year_proportion_GBD_region_prop_3  [LRI_year_proportion_GBD_region_prop_3  $Year == 2019,]
LRI_incr_2019_GBD_4 <- LRI_year_proportion_GBD_region_prop_4  [LRI_year_proportion_GBD_region_prop_4  $Year == 2019,]
LRI_incr_2019_GBD_5 <- LRI_year_proportion_GBD_region_prop_5  [LRI_year_proportion_GBD_region_prop_5  $Year == 2019,]
LRI_incr_2019_GBD_6 <- LRI_year_proportion_GBD_region_prop_6  [LRI_year_proportion_GBD_region_prop_6  $Year == 2019,]

colnames(LRI_incr_2019_L)    [2] <- "Category"
colnames(LRI_incr_2019_LM   )[2] <- "Category"
colnames(LRI_incr_2019_UM   )[2] <- "Category"
colnames(LRI_incr_2019_GBD_1)[2] <- "Category"
colnames(LRI_incr_2019_GBD_2)[2] <- "Category"
colnames(LRI_incr_2019_GBD_3)[2] <- "Category"
colnames(LRI_incr_2019_GBD_4)[2] <- "Category"
colnames(LRI_incr_2019_GBD_5)[2] <- "Category"
colnames(LRI_incr_2019_GBD_6)[2] <- "Category"

LRI_incr_2019$Category <- "Overall"  
LRI_incr_2019 <- LRI_incr_2019 %>% relocate(Category, .after = Year)

incr_df <- rbind(LRI_incr_2019, LRI_incr_2019_L, LRI_incr_2019_LM, LRI_incr_2019_UM, LRI_incr_2019_GBD_6, LRI_incr_2019_GBD_2, LRI_incr_2019_GBD_5, LRI_incr_2019_GBD_4, LRI_incr_2019_GBD_1, LRI_incr_2019_GBD_3)

incr_df$Category <- factor(
  incr_df$Category,
  levels = c( "Overall", "L", "LM", "UM", "Sub-Saharan Africa", "South Asia", "Southeast Asia, East Asia, and Oceania", "North Africa and Middle East", "Latin America and Caribbean", "Central Europe, Eastern Europe, and Central Asia"))
xpos <- seq_along(incr_df$Category)

incr_df$col <- NA
incr_df$col <- c("#0072B2", "#CC79A7", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#E69F00", "#A7B97F", "#009E73", "#0072B2")



write.csv(incr_df, "outputs/plot_df/SA_increment_proportion.csv")





# ============================================================
# Export source data for sensitivity-analysis Figure 4
# ============================================================

dir.create(
  "outputs/plot_objects",
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# Panels A-D: antibiotic access by income
# ------------------------------------------------------------

income_names <- c(
  "L"  = "Low Income Countries",
  "LM" = "Lower Middle Income Countries",
  "UM" = "Upper Middle Income Countries"
)

make_income_panel <- function(abx, panel) {
  
  bind_rows(
    lapply(names(income_names), function(inc) {
      
      df <- get(
        paste0("LRI_year_proportion_income_last_", inc)
      )
      
      df %>%
        transmute(
          Panel = panel,
          Measure = "Antibiotic treatment probability",
          Antibiotic_class = abx,
          Stratification = "Income",
          Group = income_names[[inc]],
          Year = Year,
          Median = .data[[paste0(
            "proportion_LRI_", abx, "_med"
          )]],
          Lower = .data[[paste0(
            "proportion_LRI_", abx, "_lower"
          )]],
          Upper = .data[[paste0(
            "proportion_LRI_", abx, "_upper"
          )]]
        )
    })
  )
}


fig4_A <- make_income_panel("BSP", "A")
fig4_B <- make_income_panel("NSP", "B")
fig4_C <- make_income_panel("Cephalosporins", "C")
fig4_D <- make_income_panel("Aminoglycosides", "D")


# ------------------------------------------------------------
# Panels E-H: antibiotic access by GBD region
# ------------------------------------------------------------

GBD_region_names <- c(
  "1" = "North Africa and Middle East",
  "2" = "Central Europe, Eastern Europe, and Central Asia",
  "3" = "Sub-Saharan Africa",
  "4" = "South Asia",
  "5" = "Latin America and Caribbean",
  "6" = "Southeast Asia, East Asia, and Oceania"
)

make_GBD_panel <- function(abx, panel) {
  
  bind_rows(
    lapply(1:6, function(i) {
      
      df <- get(
        paste0("LRI_year_proportion_GBD_region_", i)
      )
      
      df %>%
        transmute(
          Panel = panel,
          Measure = "Antibiotic treatment probability",
          Antibiotic_class = abx,
          Stratification = "GBD region",
          Group = unname(
            GBD_region_names[as.character(i)]
          ),
          Year = Year,
          Median = .data[[paste0(
            "proportion_LRI_", abx, "_med"
          )]],
          Lower = .data[[paste0(
            "proportion_LRI_", abx, "_lower"
          )]],
          Upper = .data[[paste0(
            "proportion_LRI_", abx, "_upper"
          )]]
        )
    })
  )
}


fig4_E <- make_GBD_panel("BSP", "E")
fig4_F <- make_GBD_panel("NSP", "F")
fig4_G <- make_GBD_panel("Cephalosporins", "G")
fig4_H <- make_GBD_panel("Aminoglycosides", "H")


# ------------------------------------------------------------
# Panels I-L:
# class-specific contribution to increase in access, 2003-2019
# ------------------------------------------------------------

make_increment_panel <- function(abx, panel) {
  
  incr_df %>%
    transmute(
      Panel = panel,
      Measure = paste0(
        "Proportion of increase in antibiotic access ",
        "from 2003 to 2019 attributable to class"
      ),
      Antibiotic_class = abx,
      
      Stratification = case_when(
        as.character(Category) == "Overall" ~ "Overall",
        as.character(Category) %in% c(
          "L", "LM", "UM"
        ) ~ "Income",
        TRUE ~ "GBD region"
      ),
      
      Group = case_when(
        as.character(Category) == "L" ~
          "Low Income Countries",
        as.character(Category) == "LM" ~
          "Lower Middle Income Countries",
        as.character(Category) == "UM" ~
          "Upper Middle Income Countries",
        TRUE ~ as.character(Category)
      ),
      
      Year = Year,
      
      Median = .data[[paste0(
        "inc_LRI_", abx, "_med"
      )]],
      Lower = .data[[paste0(
        "inc_LRI_", abx, "_lower"
      )]],
      Upper = .data[[paste0(
        "inc_LRI_", abx, "_upper"
      )]]
    )
}


fig4_I <- make_increment_panel("BSP", "I")
fig4_J <- make_increment_panel("NSP", "J")
fig4_K <- make_increment_panel("Cephalosporins", "K")
fig4_L <- make_increment_panel("Aminoglycosides", "L")


# ------------------------------------------------------------
# Combine all panels
# ------------------------------------------------------------

SA_Figure_4_source_data <- bind_rows(
  fig4_A,
  fig4_B,
  fig4_C,
  fig4_D,
  fig4_E,
  fig4_F,
  fig4_G,
  fig4_H,
  fig4_I,
  fig4_J,
  fig4_K,
  fig4_L
) %>%
  arrange(
    Panel,
    Group,
    Year
  )


# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

write.csv(
  SA_Figure_4_source_data,
  "outputs/plot_objects/SA_Figure_7_source_data.csv",
  row.names = FALSE
)






pdf("figures/Figure_S7.pdf", width = 7.5, height = 7)

layout(matrix(c(1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
                0,  4,  4,  4,  4,  5,  5,  5,  5,  6,  6,  6,  6,  7,  7,  7,  7,  8,  8,  8,  8,
                0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
                2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,
                0,  9,  9,  9,  9, 10, 10, 10, 10, 11, 11, 11, 11, 12, 12, 12, 12, 13, 13, 13, 13,
                0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
                3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,
                0, 14, 14, 14, 14, 14, 15, 15, 15, 15, 15, 16, 16, 16, 16, 16, 17, 17, 17, 17, 17 
), nrow = 8, byrow = TRUE), 
widths = c(0.75 ,1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), 
heights = c(0.1,1,0.03,0.1,1,0.03,0.1,1.6))



#### 1. (title 1) ####

par(mar = c(0, 0, 0, 0))
plot.new()
text(x = -0.035, y = 0.5, labels = "Antibiotic access by class and income", 
     adj = 0, cex = 0.7, font = 2, xpd = TRUE)


#### 2. (title 2) ####

par(mar = c(0, 0, 0, 0))
plot.new()
text(x = -0.035, y = 0.5, labels = "Antibiotic access by class and GBD region", 
     adj = 0, cex = 0.7, font = 2, xpd = TRUE)

#### 3. (title 3) ####

par(mar = c(0, 0, 0, 0))
plot.new()
text(x = -0.035, y = 0.5, labels = "Class-specific contribution to the increase in antibiotic access, 2003\u20132019", 
     adj = 0, cex = 0.7, font = 2, xpd = TRUE)

#### 5. (Fig.A) #####
par(mar = c(2.5, 2.5, 2, 1))  
par(mgp = c(3, 0.3, 0))     
par(tck = -0.025)
par(lwd = 1)


xmin <- 2003; xmax <- 2019
ymin <- log10(0.001) ; ymax <- log10(1)


plot(proportion_LRI_BSP_med ~ Year, 
     data =  LRI_year_proportion_income_last_L,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")



yticks <- log10(c(0.001, 0.01, 0.1, 1))   
ylabels <- c(0.001, 0.01, 0.1, 1)         
minor_labels <- as.vector(sapply(c(0.001, 0.01, 0.1), function(v) v * 2:9))

axis(1, at = seq(xmin, xmax, by = 4), labels = FALSE, las = 2, pos = ymin)
text(x = seq(xmin, xmax, by = 4), y = ymin-(ymax-ymin)*0.08,  
     labels = format(round(seq(xmin, xmax, by = (xmax-xmin)/6*0+4), digits = 0), nsmall = 0), xpd = TRUE, cex = 0.65)

axis(2,  at = yticks, labels = FALSE, las = 2, pos = xmin)
axis(2, at = log10(minor_labels), labels = FALSE, las = 2, pos = xmin, tck = -0.015, lwd = 0.7)
text(y = yticks, x = xmin-0.04*(xmax-xmin),  
     labels = ylabels, srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines(log10(proportion_LRI_BSP_med) ~ Year, data = LRI_year_proportion_income_last_L,  col = "#CC79A7", lwd = 1.5)
lines(log10(proportion_LRI_BSP_med) ~ Year, data = LRI_year_proportion_income_last_LM, col = "#009E73", lwd = 1.5)
lines(log10(proportion_LRI_BSP_med) ~ Year, data = LRI_year_proportion_income_last_UM, col = "#0072B2", lwd = 1.5)



# Add shaded area for confidence interval
polygon(c(LRI_year_proportion_income_last_L$Year, rev(LRI_year_proportion_income_last_L$Year)),
        c(log10(LRI_year_proportion_income_last_L$proportion_LRI_BSP_lower), rev(log10(LRI_year_proportion_income_last_L$proportion_LRI_BSP_upper))),
        col = rgb(204/255, 121/255, 167/255, 0.25), border = NA)

polygon(c(LRI_year_proportion_income_last_LM$Year, rev(LRI_year_proportion_income_last_LM$Year)),
        c(log10(LRI_year_proportion_income_last_LM$proportion_LRI_BSP_lower), rev(log10(LRI_year_proportion_income_last_LM$proportion_LRI_BSP_upper))),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

polygon(c(LRI_year_proportion_income_last_UM$Year, rev(LRI_year_proportion_income_last_UM$Year)),
        c(log10(LRI_year_proportion_income_last_UM$proportion_LRI_BSP_lower), rev(log10(LRI_year_proportion_income_last_UM$proportion_LRI_BSP_upper))),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)
mtext(expression("Antibiotic treatment probability"), side = 2, line = 2.5, cex = 0.5)  
mtext("Year", side = 1, line = 1, cex = 0.5)
mtext("A. BSP", at = xmin - 0.54*(xmax-xmin), side = 3, line = 0.5, adj = -0.1, cex = 0.5, font = 1)











#### 6. (Fig.B) #####
par(mar = c(2.5, 2.5, 2, 1))  
par(mgp = c(3, 0.3, 0))     
par(tck = -0.025)
par(lwd = 1)


xmin <- 2003; xmax <- 2019
ymin <- log10(0.001) ; ymax <- log10(1)


plot(proportion_LRI_NSP_med ~ Year, 
     data =  LRI_year_proportion_income_last_L,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")



yticks <- log10(c(0.001, 0.01, 0.1, 1))   
ylabels <- c(0.001, 0.01, 0.1, 1)         
minor_labels <- as.vector(sapply(c(0.001, 0.01, 0.1), function(v) v * 2:9))

axis(1, at = seq(xmin, xmax, by = 4), labels = FALSE, las = 2, pos = ymin)
text(x = seq(xmin, xmax, by = 4), y = ymin-(ymax-ymin)*0.08,  
     labels = format(round(seq(xmin, xmax, by = (xmax-xmin)/6*0+4), digits = 0), nsmall = 0), xpd = TRUE, cex = 0.65)

axis(2,  at = yticks, labels = FALSE, las = 2, pos = xmin)
axis(2, at = log10(minor_labels), labels = FALSE, las = 2, pos = xmin, tck = -0.015, lwd = 0.7)
text(y = yticks, x = xmin-0.04*(xmax-xmin),  
     labels = ylabels, srt = 0, adj = 1, xpd = TRUE, cex = 0.65)  

lines(log10(proportion_LRI_NSP_med) ~ Year, data = LRI_year_proportion_income_last_L,  col = "#CC79A7", lwd = 1.5)
lines(log10(proportion_LRI_NSP_med) ~ Year, data = LRI_year_proportion_income_last_LM, col = "#009E73", lwd = 1.5)
lines(log10(proportion_LRI_NSP_med) ~ Year, data = LRI_year_proportion_income_last_UM, col = "#0072B2", lwd = 1.5)



# Add shaded area for confidence interval
polygon(c(LRI_year_proportion_income_last_L$Year, rev(LRI_year_proportion_income_last_L$Year)),
        c(log10(LRI_year_proportion_income_last_L$proportion_LRI_NSP_lower), rev(log10(LRI_year_proportion_income_last_L$proportion_LRI_NSP_upper))),
        col = rgb(204/255, 121/255, 167/255, 0.25), border = NA)

polygon(c(LRI_year_proportion_income_last_LM$Year, rev(LRI_year_proportion_income_last_LM$Year)),
        c(log10(LRI_year_proportion_income_last_LM$proportion_LRI_NSP_lower), rev(log10(LRI_year_proportion_income_last_LM$proportion_LRI_NSP_upper))),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

polygon(c(LRI_year_proportion_income_last_UM$Year, rev(LRI_year_proportion_income_last_UM$Year)),
        c(log10(LRI_year_proportion_income_last_UM$proportion_LRI_NSP_lower), rev(log10(LRI_year_proportion_income_last_UM$proportion_LRI_NSP_upper))),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)
mtext("Year", side = 1, line = 1, cex = 0.5)
mtext("B. NSP", at = xmin - 0.54*(xmax-xmin), side = 3, line = 0.5, adj = -0.1, cex = 0.5, font = 1)














#### 7. (Fig.C) #####
par(mar = c(2.5, 2.5, 2, 1))  
par(mgp = c(3, 0.3, 0))     
par(tck = -0.025)
par(lwd = 1)


xmin <- 2003; xmax <- 2019
ymin <- log10(0.001) ; ymax <- log10(1)


plot(proportion_LRI_Cephalosporins_med ~ Year, 
     data =  LRI_year_proportion_income_last_L,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")



yticks <- log10(c(0.001, 0.01, 0.1, 1))   
ylabels <- c(0.001, 0.01, 0.1, 1)         
minor_labels <- as.vector(sapply(c(0.001, 0.01, 0.1), function(v) v * 2:9))

axis(1, at = seq(xmin, xmax, by = 4), labels = FALSE, las = 2, pos = ymin)
text(x = seq(xmin, xmax, by = 4), y = ymin-(ymax-ymin)*0.08,  
     labels = format(round(seq(xmin, xmax, by = (xmax-xmin)/6*0+4), digits = 0), nsmall = 0), xpd = TRUE, cex = 0.65)

axis(2,  at = yticks, labels = FALSE, las = 2, pos = xmin)
axis(2, at = log10(minor_labels), labels = FALSE, las = 2, pos = xmin, tck = -0.015, lwd = 0.7)
text(y = yticks, x = xmin-0.04*(xmax-xmin),  
     labels = ylabels, srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines(log10(proportion_LRI_Cephalosporins_med) ~ Year, data = LRI_year_proportion_income_last_L,  col = "#CC79A7", lwd = 1.5)
lines(log10(proportion_LRI_Cephalosporins_med) ~ Year, data = LRI_year_proportion_income_last_LM, col = "#009E73", lwd = 1.5)
lines(log10(proportion_LRI_Cephalosporins_med) ~ Year, data = LRI_year_proportion_income_last_UM, col = "#0072B2", lwd = 1.5)



# Add shaded area for confidence interval
polygon(c(LRI_year_proportion_income_last_L$Year, rev(LRI_year_proportion_income_last_L$Year)),
        c(log10(LRI_year_proportion_income_last_L$proportion_LRI_Cephalosporins_lower), rev(log10(LRI_year_proportion_income_last_L$proportion_LRI_Cephalosporins_upper))),
        col = rgb(204/255, 121/255, 167/255, 0.25), border = NA)

polygon(c(LRI_year_proportion_income_last_LM$Year, rev(LRI_year_proportion_income_last_LM$Year)),
        c(log10(LRI_year_proportion_income_last_LM$proportion_LRI_Cephalosporins_lower), rev(log10(LRI_year_proportion_income_last_LM$proportion_LRI_Cephalosporins_upper))),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

polygon(c(LRI_year_proportion_income_last_UM$Year, rev(LRI_year_proportion_income_last_UM$Year)),
        c(log10(LRI_year_proportion_income_last_UM$proportion_LRI_Cephalosporins_lower), rev(log10(LRI_year_proportion_income_last_UM$proportion_LRI_Cephalosporins_upper))),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)
mtext("Year", side = 1, line = 1, cex = 0.5)
mtext("C. Cephalosporins", at = xmin - 0.54*(xmax-xmin), side = 3, line = 0.5, adj = -0.1, cex = 0.5, font = 1)











#### 8. (Fig.D) #####
par(mar = c(2.5, 2.5, 2, 1))  
par(mgp = c(3, 0.3, 0))     
par(tck = -0.025)
par(lwd = 1)


xmin <- 2003; xmax <- 2019
ymin <- log10(0.001) ; ymax <- log10(1)


plot(proportion_LRI_Aminoglycosides_med ~ Year, 
     data =  LRI_year_proportion_income_last_L,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")



yticks <- log10(c(0.001, 0.01, 0.1, 1))   
ylabels <- c(0.001, 0.01, 0.1, 1)         
minor_labels <- as.vector(sapply(c(0.001, 0.01, 0.1), function(v) v * 2:9))

axis(1, at = seq(xmin, xmax, by = 4), labels = FALSE, las = 2, pos = ymin)
text(x = seq(xmin, xmax, by = 4), y = ymin-(ymax-ymin)*0.08,  
     labels = format(round(seq(xmin, xmax, by = (xmax-xmin)/6*0+4), digits = 0), nsmall = 0), xpd = TRUE, cex = 0.65)

axis(2,  at = yticks, labels = FALSE, las = 2, pos = xmin)
axis(2, at = log10(minor_labels), labels = FALSE, las = 2, pos = xmin, tck = -0.015, lwd = 0.7)
text(y = yticks, x = xmin-0.04*(xmax-xmin),  
     labels = ylabels, srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines(log10(proportion_LRI_Aminoglycosides_med) ~ Year, data = LRI_year_proportion_income_last_L,  col = "#CC79A7", lwd = 1.5)
lines(log10(proportion_LRI_Aminoglycosides_med) ~ Year, data = LRI_year_proportion_income_last_LM, col = "#009E73", lwd = 1.5)
lines(log10(proportion_LRI_Aminoglycosides_med) ~ Year, data = LRI_year_proportion_income_last_UM, col = "#0072B2", lwd = 1.5)



# Add shaded area for confidence interval
polygon(c(LRI_year_proportion_income_last_L$Year, rev(LRI_year_proportion_income_last_L$Year)),
        c(log10(LRI_year_proportion_income_last_L$proportion_LRI_Aminoglycosides_lower), rev(log10(LRI_year_proportion_income_last_L$proportion_LRI_Aminoglycosides_upper))),
        col = rgb(204/255, 121/255, 167/255, 0.25), border = NA)

polygon(c(LRI_year_proportion_income_last_LM$Year, rev(LRI_year_proportion_income_last_LM$Year)),
        c(log10(LRI_year_proportion_income_last_LM$proportion_LRI_Aminoglycosides_lower), rev(log10(LRI_year_proportion_income_last_LM$proportion_LRI_Aminoglycosides_upper))),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

polygon(c(LRI_year_proportion_income_last_UM$Year, rev(LRI_year_proportion_income_last_UM$Year)),
        c(log10(LRI_year_proportion_income_last_UM$proportion_LRI_Aminoglycosides_lower), rev(log10(LRI_year_proportion_income_last_UM$proportion_LRI_Aminoglycosides_upper))),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)

mtext("D. Aminoglycosides", at = xmin - 0.54*(xmax-xmin), side = 3, line = 0.5, adj = -0.1, cex = 0.5, font = 1)







#### 9. (legend) ####
par(mar = c(0, 0, 0, 0))
plot.new()

legend_labels <- c("Upper Middle Income Countries", "Lower Middle Income Countries", "Low Income Countries")
legend_colors <- c("#0072B2", "#009E73", "#CC79A7")
n <- length(legend_labels)

# Coordinates for the custom legend
x0 <- 0.0  
x1 <- 0.07  
y_start <- 0.6
line_spacing <- 0.07

# Draw each line, label, and shaded polygon
for (i in seq_len(n)) {
  y <- y_start - (i - 1) * line_spacing
  

  polygon(
    x = c(x0-0.01, x1+0.01, x1+0.01, x0-0.01),
    y = c(y - 0.02, y - 0.02, y + 0.02, y + 0.02),
    col = adjustcolor(legend_colors[i], alpha.f = 0.25),
    border = NA
  )
  
 
  segments(x0, y, x1, y, col = legend_colors[i], lwd = 1.5)
  

  text(x = x1 + 0.04, y = y, labels = legend_labels[i], adj = 0, cex = 0.55)
}

















#### 10. (Fig.E) #####
par(mar = c(2.5, 2.5, 2, 1))  
par(mgp = c(3, 0.3, 0))     
par(tck = -0.025)
par(lwd = 1)


xmin <- 2003; xmax <- 2019
ymin <- log10(0.001) ; ymax <- log10(1)


plot(LRI_year_proportion_GBD_region_1$proportion_LRI_BSP_med ~ LRI_year_proportion_GBD_region_1$Year, 
     data =  LRI_year_proportion_GBD_region_1,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")



yticks <- log10(c(0.001, 0.01, 0.1, 1))   
ylabels <- c(0.001, 0.01, 0.1, 1)         
minor_labels <- as.vector(sapply(c(0.001, 0.01, 0.1), function(v) v * 2:9))

axis(1, at = seq(xmin, xmax, by = 4), labels = FALSE, las = 2, pos = ymin)
text(x = seq(xmin, xmax, by = 4), y = ymin-(ymax-ymin)*0.08,  
     labels = format(round(seq(xmin, xmax, by = (xmax-xmin)/6*0+4), digits = 0), nsmall = 0), xpd = TRUE, cex = 0.65)

axis(2,  at = yticks, labels = FALSE, las = 2, pos = xmin)
axis(2, at = log10(minor_labels), labels = FALSE, las = 2, pos = xmin, tck = -0.015, lwd = 0.7)
text(y = yticks, x = xmin-0.04*(xmax-xmin),  
     labels = ylabels, srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines(log10(proportion_LRI_BSP_med) ~ Year, data = LRI_year_proportion_GBD_region_6, col = "#D55E00", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_6$Year, rev(LRI_year_proportion_GBD_region_6$Year)),
        c(log10(LRI_year_proportion_GBD_region_6$proportion_LRI_BSP_lower), rev(log10(LRI_year_proportion_GBD_region_6$proportion_LRI_BSP_upper))),
        col = rgb(204/255, 121/255, 167/255, 0.25), border = NA)

lines(log10(proportion_LRI_BSP_med) ~ Year, data = LRI_year_proportion_GBD_region_2, col = "#CC79A7", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_2$Year, rev(LRI_year_proportion_GBD_region_2$Year)),
        c(log10(LRI_year_proportion_GBD_region_2$proportion_LRI_BSP_lower), rev(log10(LRI_year_proportion_GBD_region_2$proportion_LRI_BSP_upper))),
        col = rgb(213/255, 94/255, 0/255, 0.25), border = NA)

lines(log10(proportion_LRI_BSP_med) ~ Year, data = LRI_year_proportion_GBD_region_5, col = "#E69F00", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_5$Year, rev(LRI_year_proportion_GBD_region_5$Year)),
        c(log10(LRI_year_proportion_GBD_region_5$proportion_LRI_BSP_lower), rev(log10(LRI_year_proportion_GBD_region_5$proportion_LRI_BSP_upper))),
        col = rgb(230/255, 159/255, 0/255, 0.25), border = NA)

lines(log10(proportion_LRI_BSP_med) ~ Year, data = LRI_year_proportion_GBD_region_4, col = "#A7B97F", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_4$Year, rev(LRI_year_proportion_GBD_region_4$Year)),
        c(log10(LRI_year_proportion_GBD_region_4$proportion_LRI_BSP_lower), rev(log10(LRI_year_proportion_GBD_region_4$proportion_LRI_BSP_upper))),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

lines(log10(proportion_LRI_BSP_med) ~ Year, data = LRI_year_proportion_GBD_region_1, col = "#009E73", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_1$Year, rev(LRI_year_proportion_GBD_region_1$Year)),
        c(log10(LRI_year_proportion_GBD_region_1$proportion_LRI_BSP_lower), rev(log10(LRI_year_proportion_GBD_region_1$proportion_LRI_BSP_upper))),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

lines(log10(proportion_LRI_BSP_med) ~ Year, data = LRI_year_proportion_GBD_region_3, col = "#0072B2", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_3$Year, rev(LRI_year_proportion_GBD_region_3$Year)),
        c(log10(LRI_year_proportion_GBD_region_3$proportion_LRI_BSP_lower), rev(log10(LRI_year_proportion_GBD_region_3$proportion_LRI_BSP_upper))),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)

mtext(expression("Antibiotic treatment probability"), side = 2, line = 2.5, cex = 0.5)  
mtext("Year", side = 1, line = 1, cex = 0.5)
mtext("E. BSP", at = xmin - 0.54*(xmax-xmin), side = 3, line = 0.5, adj = -0.1, cex = 0.5, font = 1)




















#### 11. (Fig.F) #####
par(mar = c(2.5, 2.5, 2, 1))  
par(mgp = c(3, 0.3, 0))     
par(tck = -0.025)
par(lwd = 1)


xmin <- 2003; xmax <- 2019
ymin <- log10(0.001) ; ymax <- log10(1)


plot(LRI_year_proportion_GBD_region_1$proportion_LRI_NSP_med ~ LRI_year_proportion_GBD_region_1$Year, 
     data =  LRI_year_proportion_GBD_region_1,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")



yticks <- log10(c(0.001, 0.01, 0.1, 1))   
ylabels <- c(0.001, 0.01, 0.1, 1)         
minor_labels <- as.vector(sapply(c(0.001, 0.01, 0.1), function(v) v * 2:9))

axis(1, at = seq(xmin, xmax, by = 4), labels = FALSE, las = 2, pos = ymin)
text(x = seq(xmin, xmax, by = 4), y = ymin-(ymax-ymin)*0.08,  
     labels = format(round(seq(xmin, xmax, by = (xmax-xmin)/6*0+4), digits = 0), nsmall = 0), xpd = TRUE, cex = 0.65)

axis(2,  at = yticks, labels = FALSE, las = 2, pos = xmin)
axis(2, at = log10(minor_labels), labels = FALSE, las = 2, pos = xmin, tck = -0.015, lwd = 0.7)
text(y = yticks, x = xmin-0.04*(xmax-xmin),  
     labels = ylabels, srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines(log10(proportion_LRI_NSP_med) ~ Year, data = LRI_year_proportion_GBD_region_6, col = "#D55E00", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_6$Year, rev(LRI_year_proportion_GBD_region_6$Year)),
        c(log10(LRI_year_proportion_GBD_region_6$proportion_LRI_NSP_lower), rev(log10(LRI_year_proportion_GBD_region_6$proportion_LRI_NSP_upper))),
        col = rgb(204/255, 121/255, 167/255, 0.25), border = NA)

lines(log10(proportion_LRI_NSP_med) ~ Year, data = LRI_year_proportion_GBD_region_2, col = "#CC79A7", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_2$Year, rev(LRI_year_proportion_GBD_region_2$Year)),
        c(log10(LRI_year_proportion_GBD_region_2$proportion_LRI_NSP_lower), rev(log10(LRI_year_proportion_GBD_region_2$proportion_LRI_NSP_upper))),
        col = rgb(213/255, 94/255, 0/255, 0.25), border = NA)

lines(log10(proportion_LRI_NSP_med) ~ Year, data = LRI_year_proportion_GBD_region_5, col = "#E69F00", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_5$Year, rev(LRI_year_proportion_GBD_region_5$Year)),
        c(log10(LRI_year_proportion_GBD_region_5$proportion_LRI_NSP_lower), rev(log10(LRI_year_proportion_GBD_region_5$proportion_LRI_NSP_upper))),
        col = rgb(230/255, 159/255, 0/255, 0.25), border = NA)

lines(log10(proportion_LRI_NSP_med) ~ Year, data = LRI_year_proportion_GBD_region_4, col = "#A7B97F", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_4$Year, rev(LRI_year_proportion_GBD_region_4$Year)),
        c(log10(LRI_year_proportion_GBD_region_4$proportion_LRI_NSP_lower), rev(log10(LRI_year_proportion_GBD_region_4$proportion_LRI_NSP_upper))),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

lines(log10(proportion_LRI_NSP_med) ~ Year, data = LRI_year_proportion_GBD_region_1, col = "#009E73", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_1$Year, rev(LRI_year_proportion_GBD_region_1$Year)),
        c(log10(LRI_year_proportion_GBD_region_1$proportion_LRI_NSP_lower), rev(log10(LRI_year_proportion_GBD_region_1$proportion_LRI_NSP_upper))),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

lines(log10(proportion_LRI_NSP_med) ~ Year, data = LRI_year_proportion_GBD_region_3, col = "#0072B2", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_3$Year, rev(LRI_year_proportion_GBD_region_3$Year)),
        c(log10(LRI_year_proportion_GBD_region_3$proportion_LRI_NSP_lower), rev(log10(LRI_year_proportion_GBD_region_3$proportion_LRI_NSP_upper))),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)

mtext("Year", side = 1, line = 1, cex = 0.5)
mtext("F. NSP", at = xmin - 0.54*(xmax-xmin), side = 3, line = 0.5, adj = -0.1, cex = 0.5, font = 1)






















#### 12. (Fig.G) #####
par(mar = c(2.5, 2.5, 2, 1))  
par(mgp = c(3, 0.3, 0))     
par(tck = -0.025)
par(lwd = 1)


xmin <- 2003; xmax <- 2019
ymin <- log10(0.001) ; ymax <- log10(1)


plot(LRI_year_proportion_GBD_region_1$proportion_LRI_Cephalosporins_med ~ LRI_year_proportion_GBD_region_1$Year, 
     data =  LRI_year_proportion_GBD_region_1,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")



yticks <- log10(c(0.001, 0.01, 0.1, 1))   
ylabels <- c(0.001, 0.01, 0.1, 1)         
minor_labels <- as.vector(sapply(c(0.001, 0.01, 0.1), function(v) v * 2:9))

axis(1, at = seq(xmin, xmax, by = 4), labels = FALSE, las = 2, pos = ymin)
text(x = seq(xmin, xmax, by = 4), y = ymin-(ymax-ymin)*0.08,  
     labels = format(round(seq(xmin, xmax, by = (xmax-xmin)/6*0+4), digits = 0), nsmall = 0), xpd = TRUE, cex = 0.65)

axis(2,  at = yticks, labels = FALSE, las = 2, pos = xmin)
axis(2, at = log10(minor_labels), labels = FALSE, las = 2, pos = xmin, tck = -0.015, lwd = 0.7)
text(y = yticks, x = xmin-0.04*(xmax-xmin),  
     labels = ylabels, srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines(log10(proportion_LRI_Cephalosporins_med) ~ Year, data = LRI_year_proportion_GBD_region_6, col = "#D55E00", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_6$Year, rev(LRI_year_proportion_GBD_region_6$Year)),
        c(log10(LRI_year_proportion_GBD_region_6$proportion_LRI_Cephalosporins_lower), rev(log10(LRI_year_proportion_GBD_region_6$proportion_LRI_Cephalosporins_upper))),
        col = rgb(204/255, 121/255, 167/255, 0.25), border = NA)

lines(log10(proportion_LRI_Cephalosporins_med) ~ Year, data = LRI_year_proportion_GBD_region_2, col = "#CC79A7", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_2$Year, rev(LRI_year_proportion_GBD_region_2$Year)),
        c(log10(LRI_year_proportion_GBD_region_2$proportion_LRI_Cephalosporins_lower), rev(log10(LRI_year_proportion_GBD_region_2$proportion_LRI_Cephalosporins_upper))),
        col = rgb(213/255, 94/255, 0/255, 0.25), border = NA)

lines(log10(proportion_LRI_Cephalosporins_med) ~ Year, data = LRI_year_proportion_GBD_region_5, col = "#E69F00", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_5$Year, rev(LRI_year_proportion_GBD_region_5$Year)),
        c(log10(LRI_year_proportion_GBD_region_5$proportion_LRI_Cephalosporins_lower), rev(log10(LRI_year_proportion_GBD_region_5$proportion_LRI_Cephalosporins_upper))),
        col = rgb(230/255, 159/255, 0/255, 0.25), border = NA)

lines(log10(proportion_LRI_Cephalosporins_med) ~ Year, data = LRI_year_proportion_GBD_region_4, col = "#A7B97F", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_4$Year, rev(LRI_year_proportion_GBD_region_4$Year)),
        c(log10(LRI_year_proportion_GBD_region_4$proportion_LRI_Cephalosporins_lower), rev(log10(LRI_year_proportion_GBD_region_4$proportion_LRI_Cephalosporins_upper))),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

lines(log10(proportion_LRI_Cephalosporins_med) ~ Year, data = LRI_year_proportion_GBD_region_1, col = "#009E73", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_1$Year, rev(LRI_year_proportion_GBD_region_1$Year)),
        c(log10(LRI_year_proportion_GBD_region_1$proportion_LRI_Cephalosporins_lower), rev(log10(LRI_year_proportion_GBD_region_1$proportion_LRI_Cephalosporins_upper))),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

lines(log10(proportion_LRI_Cephalosporins_med) ~ Year, data = LRI_year_proportion_GBD_region_3, col = "#0072B2", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_3$Year, rev(LRI_year_proportion_GBD_region_3$Year)),
        c(log10(LRI_year_proportion_GBD_region_3$proportion_LRI_Cephalosporins_lower), rev(log10(LRI_year_proportion_GBD_region_3$proportion_LRI_Cephalosporins_upper))),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)

mtext("Year", side = 1, line = 1, cex = 0.5)
mtext("G. Cephalosporins", at = xmin - 0.54*(xmax-xmin), side = 3, line = 0.5, adj = -0.1, cex = 0.5, font = 1)
















#### 13. (Fig.H) #####
par(mar = c(2.5, 2.5, 2, 1))  
par(mgp = c(3, 0.3, 0))     
par(tck = -0.025)
par(lwd = 1)


xmin <- 2003; xmax <- 2019
ymin <- log10(0.001) ; ymax <- log10(1)


plot(LRI_year_proportion_GBD_region_1$proportion_LRI_Aminoglycosides_med ~ LRI_year_proportion_GBD_region_1$Year, 
     data =  LRI_year_proportion_GBD_region_1,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")



yticks <- log10(c(0.001, 0.01, 0.1, 1))   
ylabels <- c(0.001, 0.01, 0.1, 1)         
minor_labels <- as.vector(sapply(c(0.001, 0.01, 0.1), function(v) v * 2:9))

axis(1, at = seq(xmin, xmax, by = 4), labels = FALSE, las = 2, pos = ymin)
text(x = seq(xmin, xmax, by = 4), y = ymin-(ymax-ymin)*0.08,  
     labels = format(round(seq(xmin, xmax, by = (xmax-xmin)/6*0+4), digits = 0), nsmall = 0), xpd = TRUE, cex = 0.65)

axis(2,  at = yticks, labels = FALSE, las = 2, pos = xmin)
axis(2, at = log10(minor_labels), labels = FALSE, las = 2, pos = xmin, tck = -0.015, lwd = 0.7)
text(y = yticks, x = xmin-0.04*(xmax-xmin),  
     labels = ylabels, srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines(log10(proportion_LRI_Aminoglycosides_med) ~ Year, data = LRI_year_proportion_GBD_region_6, col = "#D55E00", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_6$Year, rev(LRI_year_proportion_GBD_region_6$Year)),
        c(log10(LRI_year_proportion_GBD_region_6$proportion_LRI_Aminoglycosides_lower), rev(log10(LRI_year_proportion_GBD_region_6$proportion_LRI_Aminoglycosides_upper))),
        col = rgb(204/255, 121/255, 167/255, 0.25), border = NA)

lines(log10(proportion_LRI_Aminoglycosides_med) ~ Year, data = LRI_year_proportion_GBD_region_2, col = "#CC79A7", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_2$Year, rev(LRI_year_proportion_GBD_region_2$Year)),
        c(log10(LRI_year_proportion_GBD_region_2$proportion_LRI_Aminoglycosides_lower), rev(log10(LRI_year_proportion_GBD_region_2$proportion_LRI_Aminoglycosides_upper))),
        col = rgb(213/255, 94/255, 0/255, 0.25), border = NA)

lines(log10(proportion_LRI_Aminoglycosides_med) ~ Year, data = LRI_year_proportion_GBD_region_5, col = "#E69F00", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_5$Year, rev(LRI_year_proportion_GBD_region_5$Year)),
        c(log10(LRI_year_proportion_GBD_region_5$proportion_LRI_Aminoglycosides_lower), rev(log10(LRI_year_proportion_GBD_region_5$proportion_LRI_Aminoglycosides_upper))),
        col = rgb(230/255, 159/255, 0/255, 0.25), border = NA)

lines(log10(proportion_LRI_Aminoglycosides_med) ~ Year, data = LRI_year_proportion_GBD_region_4, col = "#A7B97F", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_4$Year, rev(LRI_year_proportion_GBD_region_4$Year)),
        c(log10(LRI_year_proportion_GBD_region_4$proportion_LRI_Aminoglycosides_lower), rev(log10(LRI_year_proportion_GBD_region_4$proportion_LRI_Aminoglycosides_upper))),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

lines(log10(proportion_LRI_Aminoglycosides_med) ~ Year, data = LRI_year_proportion_GBD_region_1, col = "#009E73", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_1$Year, rev(LRI_year_proportion_GBD_region_1$Year)),
        c(log10(LRI_year_proportion_GBD_region_1$proportion_LRI_Aminoglycosides_lower), rev(log10(LRI_year_proportion_GBD_region_1$proportion_LRI_Aminoglycosides_upper))),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

lines(log10(proportion_LRI_Aminoglycosides_med) ~ Year, data = LRI_year_proportion_GBD_region_3, col = "#0072B2", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_3$Year, rev(LRI_year_proportion_GBD_region_3$Year)),
        c(log10(LRI_year_proportion_GBD_region_3$proportion_LRI_Aminoglycosides_lower), rev(log10(LRI_year_proportion_GBD_region_3$proportion_LRI_Aminoglycosides_upper))),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)

mtext("Year", side = 1, line = 1, cex = 0.5)
mtext("H. Aminoglycosides", at = xmin - 0.54*(xmax-xmin), side = 3, line = 0.5, adj = -0.1, cex = 0.5, font = 1)




#### 14. (legend) ####
par(mar = c(0, 0, 0, 0))
plot.new()



legend_labels <- c("Southeast Asia, East Asia, and Oceania",
                   "Central Europe, Eastern Europe, and Central Asia",
                   "Latin America and Caribbean",
                   "South Asia",
                   "North Africa and Middle East",
                   "Sub-Saharan Africa")
legend_colors <- c("#D55E00",
                   "#CC79A7",
                   "#E69F00",
                   "#A7B97F",
                   "#009E73",
                   "#0072B2")
n <- length(legend_labels)




# Coordinates for the custom legend
x0 <- 0.0  
x1 <- 0.07  
y_start <- 0.7
line_spacing <- 0.07

# Draw each line, label, and shaded polygon
for (i in seq_len(n)) {
  y <- y_start - (i - 1) * line_spacing
  

  polygon(
    x = c(x0-0.01, x1+0.01, x1+0.01, x0-0.01),
    y = c(y - 0.02, y - 0.02, y + 0.02, y + 0.02),
    col = adjustcolor(legend_colors[i], alpha.f = 0.25),
    border = NA
  )
  
 
  segments(x0, y, x1, y, col = legend_colors[i], lwd = 1.5)
  

  text(x = x1 + 0.04, y = y, labels = legend_labels[i], adj = 0, cex = 0.55)
}











#### 15. (Fig.I) #####
par(mar = c(9, 4, 2, 1))  
par(mgp = c(3, 0.3, 0))     
par(tck = -0.025)
par(lwd = 1)


xmin <- 0.5; xmax <- length(xpos) + 0.5
ymin <- log10(0.001) ; ymax <- log10(1)


plot(log10(inc_LRI_BSP_med) ~ xpos, 
     data =  incr_df,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")



axis(1, at = seq(xmin-0.1, xmax, by = 3), labels = FALSE, las = 2, pos = ymin, tck = 0)

axis(1, at = xpos, labels = FALSE, pos = ymin)
text(x = xpos, y = ymin - 0.06 * (ymax - ymin), labels = incr_df$Category, srt = 45, adj = 1, xpd = TRUE, cex = 0.6)

yticks <- log10(c(0.001, 0.01, 0.1, 1))   
ylabels <- c(0.001, 0.01, 0.1, 1)         
minor_labels <- as.vector(sapply(c(0.001, 0.01, 0.1), function(v) v * 2:9))

axis(2,  at = yticks, labels = FALSE, las = 2, pos = xmin)
axis(2, at = log10(minor_labels), labels = FALSE, las = 2, pos = xmin, tck = -0.015, lwd = 0.7)
text(y = yticks, x = xmin-0.04*(xmax-xmin),  
     labels = ylabels, srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

segments(x0 = xpos, x1 = xpos, y0 = log10(incr_df$inc_LRI_BSP_lower), y1 = log10(incr_df$inc_LRI_BSP_upper), col = incr_df$col, lwd = 1.2)
points(xpos, log10(incr_df$inc_LRI_BSP_med), pch = 16, col = incr_df$col, cex = 0.8)

mtext(expression("Proportion"), side = 2, line = 2.5, cex = 0.5)  
mtext("I. BSP", at = xmin - 0.45*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.5, font = 1)

div_x <- c(1.5, 4.5)
segments(x0 = div_x, y0 = ymin+0.01, 
         x1 = div_x, y1 = ymax, 
         col = "grey70", lwd = 1)







#### 16. (Fig.J) #####
par(mar = c(9, 4, 2, 1))  
par(mgp = c(3, 0.3, 0))     
par(tck = -0.025)
par(lwd = 1)


xmin <- 0.5; xmax <- length(xpos) + 0.5
ymin <- log10(0.001) ; ymax <- log10(1)


plot(log10(inc_LRI_NSP_med) ~ xpos, 
     data =  incr_df,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")



axis(1, at = seq(xmin-0.1, xmax, by = 3), labels = FALSE, las = 2, pos = ymin, tck = 0)

axis(1, at = xpos, labels = FALSE, pos = ymin)
text(x = xpos, y = ymin - 0.06 * (ymax - ymin), labels = incr_df$Category, srt = 45, adj = 1, xpd = TRUE, cex = 0.6)

yticks <- log10(c(0.001, 0.01, 0.1, 1))   
ylabels <- c(0.001, 0.01, 0.1, 1)         
minor_labels <- as.vector(sapply(c(0.001, 0.01, 0.1), function(v) v * 2:9))

axis(2,  at = yticks, labels = FALSE, las = 2, pos = xmin)
axis(2, at = log10(minor_labels), labels = FALSE, las = 2, pos = xmin, tck = -0.015, lwd = 0.7)
text(y = yticks, x = xmin-0.04*(xmax-xmin),  
     labels = ylabels, srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

incr_df$inc_LRI_NSP_lower[c(2,10)] <- 0.001

points(xpos, log10(incr_df$inc_LRI_NSP_med), pch = 16, col = incr_df$col, cex = 0.8)

arrow_indices <- c(2, 10)
segment_indices <- c(1, 3, 4, 5, 6, 7, 8, 9)

segments(x0 = xpos[segment_indices], x1 = xpos[segment_indices], 
         y0 = log10(incr_df$inc_LRI_NSP_lower[segment_indices]), 
         y1 = log10(incr_df$inc_LRI_NSP_upper[segment_indices]), 
         col = incr_df$col[segment_indices], lwd = 1.2)

arrows(x0 = xpos[arrow_indices], x1 = xpos[arrow_indices], 
       y1 = log10(incr_df$inc_LRI_NSP_upper[arrow_indices]), 
       y0 = ymin, 
       code = 1,
       length = 0.04, col = incr_df$col[arrow_indices], lwd = 1.2)

mtext("J. NSP", at = xmin - 0.45*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.5, font = 1)

div_x <- c(1.5, 4.5)
segments(x0 = div_x, y0 = ymin+0.01, 
         x1 = div_x, y1 = ymax, 
         col = "grey70", lwd = 1)













#### 17. (Fig.K) #####

par(mar = c(9, 4, 2, 1))  
par(mgp = c(3, 0.3, 0))     
par(tck = -0.025)
par(lwd = 1)


xmin <- 0.5; xmax <- length(xpos) + 0.5
ymin <- log10(0.001) ; ymax <- log10(1)


plot(log10(inc_LRI_Cephalosporins_med) ~ xpos, 
     data =  incr_df,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")



axis(1, at = seq(xmin-0.1, xmax, by = 3), labels = FALSE, las = 2, pos = ymin, tck = 0)

axis(1, at = xpos, labels = FALSE, pos = ymin)
text(x = xpos, y = ymin - 0.06 * (ymax - ymin), labels = incr_df$Category, srt = 45, adj = 1, xpd = TRUE, cex = 0.6)

yticks <- log10(c(0.001, 0.01, 0.1, 1))   
ylabels <- c(0.001, 0.01, 0.1, 1)         
minor_labels <- as.vector(sapply(c(0.001, 0.01, 0.1), function(v) v * 2:9))

axis(2,  at = yticks, labels = FALSE, las = 2, pos = xmin)
axis(2, at = log10(minor_labels), labels = FALSE, las = 2, pos = xmin, tck = -0.015, lwd = 0.7)
text(y = yticks, x = xmin-0.04*(xmax-xmin),  
     labels = ylabels, srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

segments(x0 = xpos, x1 = xpos, y0 = log10(incr_df$inc_LRI_Cephalosporins_lower), y1 = log10(incr_df$inc_LRI_Cephalosporins_upper), col = incr_df$col, lwd = 1.2)
points(xpos, log10(incr_df$inc_LRI_Cephalosporins_med), pch = 16, col = incr_df$col, cex = 0.8)

mtext("K. Cephalosporins", at = xmin - 0.45*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.5, font = 1)

div_x <- c(1.5, 4.5)
segments(x0 = div_x, y0 = ymin+0.01, 
         x1 = div_x, y1 = ymax, 
         col = "grey70", lwd = 1)











#### 18. (Fig.L) #####

par(mar = c(9, 4, 2, 1))  
par(mgp = c(3, 0.3, 0))     
par(tck = -0.025)
par(lwd = 1)


xmin <- 0.5; xmax <- length(xpos) + 0.5
ymin <- log10(0.001) ; ymax <- log10(1)


plot(log10(inc_LRI_Aminoglycosides_med) ~ xpos, 
     data =  incr_df,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")



axis(1, at = seq(xmin-0.1, xmax, by = 3), labels = FALSE, las = 2, pos = ymin, tck = 0)

axis(1, at = xpos, labels = FALSE, pos = ymin)
text(x = xpos, y = ymin - 0.06 * (ymax - ymin), labels = incr_df$Category, srt = 45, adj = 1, xpd = TRUE, cex = 0.6)

yticks <- log10(c(0.001, 0.01, 0.1, 1))   
ylabels <- c(0.001, 0.01, 0.1, 1)         
minor_labels <- as.vector(sapply(c(0.001, 0.01, 0.1), function(v) v * 2:9))

axis(2,  at = yticks, labels = FALSE, las = 2, pos = xmin)
axis(2, at = log10(minor_labels), labels = FALSE, las = 2, pos = xmin, tck = -0.015, lwd = 0.7)
text(y = yticks, x = xmin-0.04*(xmax-xmin),  
     labels = ylabels, srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

segments(x0 = xpos, x1 = xpos, y0 = log10(incr_df$inc_LRI_Aminoglycosides_lower), y1 = log10(incr_df$inc_LRI_Aminoglycosides_upper), col = incr_df$col, lwd = 1.2)
points(xpos, log10(incr_df$inc_LRI_Aminoglycosides_med), pch = 16, col = incr_df$col, cex = 0.8)

mtext("L. Aminoglycosides", at = xmin - 0.45*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.5, font = 1)

div_x <- c(1.5, 4.5)
segments(x0 = div_x, y0 = ymin+0.01, 
         x1 = div_x, y1 = ymax, 
         col = "grey70", lwd = 1)



dev.off()

