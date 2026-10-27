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
ab.data_LRI <- readRDS("outputs/ab.data_LRI.UNI.rds")
matrix.LRI <- readRDS("outputs/U5_DDD_LRI_matrix.rds")

matrix.meta <- ab.data_LRI %>% dplyr::select(Country, Year)
matrix.LRI <- matrix.meta %>% inner_join(matrix.LRI, by = c("Country", "Year"))
matrix.LRI <- matrix.LRI[, c(1:30002, 40003:50002)]

matrix.LRI.sum <- data.frame(matrix(NA_real_, nrow = nrow(matrix.LRI), ncol = (n_sim + 2)))
matrix.LRI.sum[, 1:2] <- as.matrix(matrix.LRI[, 1:2])

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

dir.create(
  "outputs/antibiotic_class",
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  coef_df_LRI_DDD,
  "outputs/antibiotic_class/coef_df_LRI_DDD_by_4class.rds"
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

UNI_indexed <- readRDS("outputs/UNI_indexed.rds")
coef_df_LRI <- readRDS("outputs/coef_df_LRI.rds")
coef_df_LRI_DDD <- readRDS("outputs/antibiotic_class/coef_df_LRI_DDD_by_4class.rds")

UNI_indexed_Abaccess <- UNI_indexed %>% filter(Abaccess == 1)
UNI_indexed_Abdata <- UNI_indexed %>% filter(ab.data == 1)
UNI_indexed_both <- UNI_indexed %>% filter(Abaccess == 1, ab.data == 1)

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

dir.create(
  "outputs/plot_df/p_best",
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  LRI_DDD_ratio_sim,
  "outputs/plot_df/p_best/ratio_LRI4.rds"
)

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
  "outputs/plot_df/p_best/p_LRI4.rds"
)


