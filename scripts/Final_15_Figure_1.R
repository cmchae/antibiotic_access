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
coef_df_LRI <- readRDS("outputs/coef_df_LRI.rds")
coef_df_Fever <- readRDS("outputs/coef_df_Fever.rds")
coef_df_LRI_DDD <- readRDS("outputs/coef_df_LRI_DDD.rds")
coef_df_Fever_DDD <- readRDS("outputs/coef_df_Fever_DDD.rds")

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
  Probability_LRI[[as.character(yr)]] <- read.csv(paste0("outputs/global_map_df/1.1_Untreated_number_", yr, "_LRI.csv"))
  Probability_Fever_0[[as.character(yr)]] <- read.csv(paste0("outputs/global_map_df/1.1_Untreated_number_0_", yr, "_Fever.csv"))
  Probability_Fever_1[[as.character(yr)]] <- read.csv(paste0("outputs/global_map_df/1.1_Untreated_number_1_", yr, "_Fever.csv"))
  DDD_per_treated_LRI[[as.character(yr)]] <- read.csv(paste0("outputs/global_map_df/2.1_DDD_LRI_", yr, ".csv"))
  DDD_per_treated_Fever[[as.character(yr)]] <- read.csv(paste0("outputs/global_map_df/2.1_DDD_Fever_", yr, ".csv"))
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
      legend.position = c(0.03, 0.15),        # (x, y) coordinates (0 = left/bottom, 1 = right/top)
      legend.justification = c(0, 0),         # Anchor the legend box to its bottom-left corner
      legend.title = element_text(size = 6), 
      legend.text = element_text(size = 6),
      legend.key.height = unit(0.5, "cm"),
      legend.key.width = unit(0.3, "cm"),
      plot.title = element_text(size = 8, face = "bold"),
      plot.subtitle = element_text(size = 8)
    )
  
  plot_list_probability_LRI[[period_name]] <- p
}

tiff("outputs/table/by_year_malaria/probability_LRTI_2004_2006.tiff", width = 7.5, height = 5, units = "in", res = 300)
plot_list_probability_LRI[["2004_2006"]]
dev.off()
tiff("outputs/table/by_year_malaria/probability_LRTI_2009_2011.tiff", width = 7.5, height = 5, units = "in", res = 300)
plot_list_probability_LRI[["2009_2011"]]
dev.off()
tiff("outputs/table/by_year_malaria/probability_LRTI_2014_2016.tiff", width = 7.5, height = 5, units = "in", res = 300)
plot_list_probability_LRI[["2014_2016"]]
dev.off()
tiff("outputs/table/by_year_malaria/probability_LRTI_2019.tiff", width = 7.5, height = 5, units = "in", res = 300)
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
      legend.position = c(0.03, 0.15),        # (x, y) coordinates (0 = left/bottom, 1 = right/top)
      legend.justification = c(0, 0),         # Anchor the legend box to its bottom-left corner
      legend.title = element_text(size = 6), 
      legend.text = element_text(size = 6),
      legend.key.height = unit(0.5, "cm"),
      legend.key.width = unit(0.3, "cm"),
      plot.title = element_text(size = 8, face = "bold"),
      plot.subtitle = element_text(size = 8)
    )
  
  plot_list_probability_Fever[[period_name]] <- p
}

tiff("outputs/table/by_year_malaria/probability_Fever_2004_2006.tiff", width = 7.5, height = 5, units = "in", res = 300)
plot_list_probability_Fever[["2004_2006"]]
dev.off()
tiff("outputs/table/by_year_malaria/probability_Fever_2009_2011.tiff", width = 7.5, height = 5, units = "in", res = 300)
plot_list_probability_Fever[["2009_2011"]]
dev.off()
tiff("outputs/table/by_year_malaria/probability_Fever_2014_2016.tiff", width = 7.5, height = 5, units = "in", res = 300)
plot_list_probability_Fever[["2014_2016"]]
dev.off()
tiff("outputs/table/by_year_malaria/probability_Fever_2019.tiff", width = 7.5, height = 5, units = "in", res = 300)
plot_list_probability_Fever[["2019"]]
dev.off()

combined_plot <- (
  plot_list_probability_LRI[["2019"]] + # First plot
    ggtitle("A. Global antibiotic access for LRTI in 2019") + # Title for the first plot
    theme(plot.margin = unit(c(0, 0, 0, 0), "pt"),
          plot.title = element_text(hjust = +0.03),
          panel.grid = element_blank()) # Top, Right, Bottom, Left margins
  ) / plot_spacer() / (
  plot_list_probability_Fever[["2019"]] + # Second plot
    ggtitle("B. Global antibiotic access for AFI in 2019") + # Title for the second plot
    theme(plot.margin = unit(c(0, 0, 0, 0), "pt"),
          plot.title = element_text(hjust = +0.03),
          panel.grid = element_blank()) # Top, Right, Bottom, Left margins
  ) +
  plot_layout(heights = c(1, 0.04, 1))

dir.create(
  "figures",
  recursive = TRUE,
  showWarnings = FALSE
)

# 2. Save the combined plot to PDF
ggsave("figures/Figure_1.pdf", plot = combined_plot, width = 7.5, height = 7.5)
















# ============================================================
# Export country-level data used for plotting
# ============================================================

periods <- list(
  "2004_2006" = 2004:2006,
  "2009_2011" = 2009:2011,
  "2014_2016" = 2014:2016,
  "2019" = 2019
)

# ---- Antibiotic treatment probability ----

plot_data_probability_LRTI <- bind_rows(
  lapply(names(periods), function(period_name) {
    
    years <- periods[[period_name]]
    
    LRI_joint %>%
      filter(Year %in% years) %>%
      group_by(Country) %>%
      summarize(
        Value = sum(Treated_incidence) / sum(Total_incidence),
        .groups = "drop"
      ) %>%
      mutate(
        Syndrome = "LRTI",
        Measure = "Antibiotic treatment probability",
        Period = period_name
      )
  })
)

plot_data_probability_AFI <- bind_rows(
  lapply(names(periods), function(period_name) {
    
    years <- periods[[period_name]]
    
    Fever_joint %>%
      filter(Year %in% years) %>%
      group_by(Country) %>%
      summarize(
        Value = sum(Treated_incidence) / sum(Total_incidence),
        .groups = "drop"
      ) %>%
      mutate(
        Syndrome = "AFI",
        Measure = "Antibiotic treatment probability",
        Period = period_name
      )
  })
)

# ---- Combine and export ----

plot_data_all <- bind_rows(
  plot_data_probability_LRTI,
  plot_data_probability_AFI,
) %>%
  dplyr::select(
    Country,
    Syndrome,
    Measure,
    Period,
    Value
  ) %>%
  arrange(Measure, Syndrome, Period, Country)

write.csv(
  plot_data_all,
  "outputs/plot_objects/figure_1.csv",
  row.names = FALSE
)






