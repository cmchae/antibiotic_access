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
analysis_data <- readRDS("outputs/AbLRI.data.UNI.rds")
AbAccess <- analysis_data %>%
  dplyr::select(c("ChildSex","ChildAgeMonth","Year","GDP","Income","Country",
                  "AbLRI","AbFever.nonARI","AbLRI_SE","AbFever.nonARI_SE","PC1","PC2","Endemic","gbd_region","n_AbLRI","n_LRI","n_Ab_Fever.nonARI","n_Fever.nonARI"))

AbAccess <- na.omit(AbAccess)
AbAccess <- AbAccess %>% mutate(PersonTime = 1/26) 
AbAccess$Year <- as.integer(AbAccess$Year)
unique(AbAccess$Country)

ab.data_LRI <- readRDS("outputs/ab.data_LRI.UNI.rds")
ab.data_Fever <- readRDS("outputs/ab.data_Fever.UNI.rds")
ab.data_all <- readRDS("outputs/ab.data_all.UNI.rds")
matrix.LRI <- readRDS("outputs/U5_DDD_LRI_matrix.rds")
matrix.Fever <- readRDS("outputs/U5_DDD_Fever_only_matrix.rds")
matrix.all <- readRDS("outputs/U5_DDD_all_matrix.rds")

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

dir.create(
  "outputs/draft_df",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  AbAccess_Country_level,
  "outputs/draft_df/Suppl_1_DHS_model_data.csv",
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
pdf("outputs/(B) Suppl_Figure_logit_pois_Ab_LRTI_and_Fever.pdf", width = 6, height = 2.5)

par(mfrow = c(1, 2))
par(mar = c(3, 3.5, 1.5, 2.5))  # Bottom, Left, Top, Right
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

axis(1, at = seq(0.5, 0.9, by = 0.1), cex.axis = 0.8, mgp = c(3, 0.20, 0))
axis(2, at = seq(0.5, 0.9, by = 0.1), las = 1, cex.axis = 0.8)

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

axis(1, at = seq(0.4, 0.7, by = 0.1), cex.axis = 0.8, mgp = c(3, 0.20, 0))
axis(2, at = seq(0.4, 0.7, by = 0.1), las = 1, cex.axis = 0.8)

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
write.csv(tab_wide, "outputs/draft_df/Suppl_2_IQVIA_IHME_data.csv", row.names = FALSE)


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

saveRDS(UNI_indexed, "outputs/UNI_indexed.rds")
saveRDS(coef_df_AbLRI, "outputs/coef_df_LRI.rds")
saveRDS(coef_df_AbFever, "outputs/coef_df_Fever.rds")
saveRDS(coef_df_LRI_DDD, "outputs/coef_df_LRI_DDD.rds")
saveRDS(coef_df_Fever_DDD, "outputs/coef_df_Fever_DDD.rds")
saveRDS(coef_df_all_DDD, "outputs/coef_df_all_DDD.rds")