rm(list = ls())
gc()

requiredPackages = c("dplyr", "tidyr", "pbapply", "parallel")

lib.path <- .libPaths()
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p, lib = lib.path)
  library(p,character.only = TRUE)
}

options(scipen = 100, digits = 4)

set.seed(123)
n_sim <- 10000

gbd.data <- read.csv("data/All_inf2.csv")
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
  
  fitted_mean <- shape * scale   # Fitted values
  fitted_lower <- qgamma(0.025, shape = shape, scale = scale)
  fitted_upper <- qgamma(0.975, shape = shape, scale = scale)

  error_val <- abs(fitted_mean - val)/val  # errors
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

prop_fever_RI <- readRDS("outputs/prop.Fever.RI_logit.rds")
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
saveRDS(U5_matrix, "outputs/matrix.result.ver.1.rds")

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

saveRDS(U5_prop, "outputs/U5_prop.rds")
saveRDS(merged.data, "outputs/Incidence_LRI_Fever_Both.ver.1.rds")

# Incidence of LRI (w/ and w/o fever) and Fever
U5_matrix <- data.frame(
  Country = merged.data$Country,
  Year = merged.data$Year,
  U5_LRI = U5_LRI_matrix,
  U5_fever_only = U5_fever_only_matrix)

saveRDS(U5_matrix, "outputs/matrix.result.ver.2.rds")

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

saveRDS(merged.data, "outputs/Incidence_LRI_Fever_Both.ver.2.rds")








