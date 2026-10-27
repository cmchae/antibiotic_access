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

coef_df_LRI <- readRDS("outputs/coef_df_LRI.rds")
coef_df_Fever <- readRDS("outputs/coef_df_Fever.rds")
coef_df_LRI_DDD <- readRDS("outputs/coef_df_LRI_DDD.rds")
coef_df_Fever_DDD <- readRDS("outputs/coef_df_Fever_DDD.rds")
coef_df_all_DDD <- readRDS("outputs/coef_df_all_DDD.rds")

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

LRI_year_proportion      <- read.csv("outputs/plot_df/1.3_plot_Ab_LRI_year_proportion.csv")
LRI_year_proportion_L    <- read.csv("outputs/plot_df/1.3_plot_Ab_LRI_year_proportion_L.csv")
LRI_year_proportion_LM   <- read.csv("outputs/plot_df/1.3_plot_Ab_LRI_year_proportion_LM.csv")
LRI_year_proportion_UM   <- read.csv("outputs/plot_df/1.3_plot_Ab_LRI_year_proportion_UM.csv")
LRI_year_proportion_malaria_0      <- read.csv("outputs/plot_df/1.4_plot_Ab_LRI_year_proportion_malaria_0.csv")
LRI_year_proportion_malaria_1      <- read.csv("outputs/plot_df/1.4_plot_Ab_LRI_year_proportion_malaria_1.csv")
Fever_year_proportion    <- read.csv("outputs/plot_df/1.3_plot_Ab_Fever_year_proportion.csv")
Fever_year_proportion_mal0_L  <- read.csv("outputs/plot_df/1.3_plot_Ab_Fever_year_proportion_mal0_L.csv")
Fever_year_proportion_mal0_LM <- read.csv("outputs/plot_df/1.3_plot_Ab_Fever_year_proportion_mal0_LM.csv")
Fever_year_proportion_mal0_UM <- read.csv("outputs/plot_df/1.3_plot_Ab_Fever_year_proportion_mal0_UM.csv")
Fever_year_proportion_mal1_L  <- read.csv("outputs/plot_df/1.3_plot_Ab_Fever_year_proportion_mal1_L.csv")
Fever_year_proportion_mal1_LM <- read.csv("outputs/plot_df/1.3_plot_Ab_Fever_year_proportion_mal1_LM.csv")
Fever_year_proportion_mal1_UM <- read.csv("outputs/plot_df/1.3_plot_Ab_Fever_year_proportion_mal1_UM.csv")
Fever_year_proportion_malaria_0      <- read.csv("outputs/plot_df/1.4_plot_Ab_Fever_year_proportion_malaria_0.csv")
Fever_year_proportion_malaria_1      <- read.csv("outputs/plot_df/1.4_plot_Ab_Fever_year_proportion_malaria_1.csv")

for (i in 1:7) { # GBD_region
  path_LRI <- paste0("outputs/plot_df/1.5_plot_Ab_LRI_year_proportion_GBD_region_", i, ".csv")
  obj_LRI <- paste0("LRI_year_proportion_GBD_region_", i)
  assign(obj_LRI, read.csv(path_LRI, stringsAsFactors = FALSE))
}

Malaria_level <- c(0,1)

for (mal in Malaria_level) {
  for (i in 1:7) {
    obj_Fever <- paste0("Fever_year_proportion_GBD_region_", mal,"_", i)
    path_Fever <- paste0("outputs/plot_df/1.5_plot_Ab_Fever_year_proportion_GBD_region_", mal,"_", i,".csv")
    assign(obj_Fever, read.csv(path_Fever, stringsAsFactors = FALSE))
  }
}

# read the plot data generated from the current PC1-PC2 models: DDD
LRI_year_DDD_per_treated      <- read.csv("outputs/plot_df/2.3_plot_Ab_LRI_year_DDD_per_treated.csv")
LRI_year_DDD_per_treated_L    <- read.csv("outputs/plot_df/2.3_plot_Ab_LRI_year_DDD_per_treated_L.csv")
LRI_year_DDD_per_treated_LM   <- read.csv("outputs/plot_df/2.3_plot_Ab_LRI_year_DDD_per_treated_LM.csv")
LRI_year_DDD_per_treated_UM   <- read.csv("outputs/plot_df/2.3_plot_Ab_LRI_year_DDD_per_treated_UM.csv")
Fever_year_DDD_per_treated    <- read.csv("outputs/plot_df/2.3_plot_Ab_Fever_year_DDD_per_treated.csv")
Fever_year_DDD_per_treated_L  <- read.csv("outputs/plot_df/2.3_plot_Ab_Fever_year_DDD_per_treated_L.csv")
Fever_year_DDD_per_treated_LM <- read.csv("outputs/plot_df/2.3_plot_Ab_Fever_year_DDD_per_treated_LM.csv")
Fever_year_DDD_per_treated_UM <- read.csv("outputs/plot_df/2.3_plot_Ab_Fever_year_DDD_per_treated_UM.csv")
all_year_DDD_per_treated      <- read.csv("outputs/plot_df/2.3_plot_Ab_all_year_DDD_per_treated.csv")
all_year_DDD_per_treated_L    <- read.csv("outputs/plot_df/2.3_plot_Ab_all_year_DDD_per_treated_L.csv")
all_year_DDD_per_treated_LM   <- read.csv("outputs/plot_df/2.3_plot_Ab_all_year_DDD_per_treated_LM.csv")
all_year_DDD_per_treated_UM   <- read.csv("outputs/plot_df/2.3_plot_Ab_all_year_DDD_per_treated_UM.csv")


# ============================================================
# Export source data for Figure 2
# ============================================================

# Helper function
format_fig2 <- function(df, panel, syndrome,
                        malaria_status = "Not applicable",
                        stratification = "Overall",
                        group = "Overall") {
  
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

fig2_A <- format_fig2(
  LRI_year_proportion,
  panel = "A",
  syndrome = "LRTI"
)

fig2_B <- format_fig2(
  Fever_year_proportion_malaria_0,
  panel = "B",
  syndrome = "AFI",
  malaria_status = "Malaria non-endemic"
)

fig2_C <- format_fig2(
  Fever_year_proportion_malaria_1,
  panel = "C",
  syndrome = "AFI",
  malaria_status = "Malaria endemic"
)


# ------------------------------------------------------------
# Panels D-F: antibiotic access by income
# ------------------------------------------------------------

fig2_D <- bind_rows(
  format_fig2(
    LRI_year_proportion_L,
    "D", "LRTI",
    stratification = "Income",
    group = "Low income"
  ),
  format_fig2(
    LRI_year_proportion_LM,
    "D", "LRTI",
    stratification = "Income",
    group = "Lower middle income"
  ),
  format_fig2(
    LRI_year_proportion_UM,
    "D", "LRTI",
    stratification = "Income",
    group = "Upper middle income"
  )
)

fig2_E <- bind_rows(
  format_fig2(
    Fever_year_proportion_mal0_L,
    "E", "AFI",
    malaria_status = "Malaria non-endemic",
    stratification = "Income",
    group = "Low income"
  ),
  format_fig2(
    Fever_year_proportion_mal0_LM,
    "E", "AFI",
    malaria_status = "Malaria non-endemic",
    stratification = "Income",
    group = "Lower middle income"
  ),
  format_fig2(
    Fever_year_proportion_mal0_UM,
    "E", "AFI",
    malaria_status = "Malaria non-endemic",
    stratification = "Income",
    group = "Upper middle income"
  )
)

fig2_F <- bind_rows(
  format_fig2(
    Fever_year_proportion_mal1_L,
    "F", "AFI",
    malaria_status = "Malaria endemic",
    stratification = "Income",
    group = "Low income"
  ),
  format_fig2(
    Fever_year_proportion_mal1_LM,
    "F", "AFI",
    malaria_status = "Malaria endemic",
    stratification = "Income",
    group = "Lower middle income"
  ),
  format_fig2(
    Fever_year_proportion_mal1_UM,
    "F", "AFI",
    malaria_status = "Malaria endemic",
    stratification = "Income",
    group = "Upper middle income"
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


# Panel G: LRTI
fig2_G <- bind_rows(
  lapply(c(1, 2, 3, 5, 6, 7), function(i) {
    
    df <- get(paste0("LRI_year_proportion_GBD_region_", i))
    
    format_fig2(
      df,
      panel = "G",
      syndrome = "LRTI",
      stratification = "GBD region",
      group = GBD_region_names[as.character(i)]
    )
  })
)


# Panel H: AFI, malaria non-endemic
fig2_H <- bind_rows(
  lapply(c(1, 2, 3, 5, 6, 7), function(i) {
    
    df <- get(
      paste0("Fever_year_proportion_GBD_region_0_", i)
    )
    
    format_fig2(
      df,
      panel = "H",
      syndrome = "AFI",
      malaria_status = "Malaria non-endemic",
      stratification = "GBD region",
      group = GBD_region_names[as.character(i)]
    )
  })
)


# Panel I: AFI, malaria endemic
fig2_I <- bind_rows(
  lapply(c(1, 2, 3, 5, 6, 7), function(i) {
    
    df <- get(
      paste0("Fever_year_proportion_GBD_region_1_", i)
    )
    
    format_fig2(
      df,
      panel = "I",
      syndrome = "AFI",
      malaria_status = "Malaria endemic",
      stratification = "GBD region",
      group = GBD_region_names[as.character(i)]
    )
  })
)


# ------------------------------------------------------------
# Combine and export
# ------------------------------------------------------------

Figure_2_source_data <- bind_rows(
  fig2_A,
  fig2_B,
  fig2_C,
  fig2_D,
  fig2_E,
  fig2_F,
  fig2_G,
  fig2_H,
  fig2_I
) %>%
  arrange(Panel, Group, Year)

write.csv(
  Figure_2_source_data,
  "outputs/plot_objects/Figure_2_source_data.csv",
  row.names = FALSE
)


pdf("figures/Figure_2.pdf", width = 7.5, height = 7)

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

# plot(type = "n", plot_data$GDP_Capita_PPP, plot_data$AbARI,
xmin <- 2003; xmax <- 2019
ymin <- 0; ymax <- 1

par(mar = c(2.5, 2.5, 2, 1))  # Bottom, Left, Top, Right
par(mgp = c(3, 0.3, 0))     # Axis label and tick distance
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

# title(main = "Probability of treatment for LRI ~ Antibiotics consumption per incidence", cex.main = 1.5)

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

# plot(type = "n", plot_data$GDP_Capita_PPP, plot_data$AbARI,
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

# title(main = "Probability of treatment for Fever ~ Antibiotics consumption per incidence", cex.main = 1.5)

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

# title(main = "Probability of treatment for Fever ~ Antibiotics consumption per incidence", cex.main = 1.5)
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

# plot(type = "n", plot_data$GDP_Capita_PPP, plot_data$AbARI,
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

# title(main = "Probability of treatment for LRI ~ Antibiotics consumption per incidence", cex.main = 1.5)

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

# plot(type = "n", plot_data$GDP_Capita_PPP, plot_data$AbARI,
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

# title(main = "Probability of treatment for Fever ~ Antibiotics consumption per incidence", cex.main = 1.5)

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

# plot(type = "n", plot_data$GDP_Capita_PPP, plot_data$AbARI,
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

# title(main = "Probability of treatment for Fever ~ Antibiotics consumption per incidence", cex.main = 1.5)

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
x0 <- 0.0  # start x for line
x1 <- 0.07  # end x for line
y_start <- 0.6
line_spacing <- 0.07

# Draw each line, label, and shaded polygon
for (i in seq_len(n)) {
  y <- y_start - (i - 1) * line_spacing
  
  # Polygon around line (±0.01 vertically for shading)
  polygon(
    x = c(x0-0.01, x1+0.01, x1+0.01, x0-0.01),
    y = c(y - 0.02, y - 0.02, y + 0.02, y + 0.02),
    col = adjustcolor(legend_colors[i], alpha.f = 0.25),
    border = NA
  )
  
  # Line
  segments(x0, y, x1, y, col = legend_colors[i], lwd = 1.5)
  
  # Text
  text(x = x1 + 0.04, y = y, labels = legend_labels[i], adj = 0, cex = 0.56)
}





#### 11. (Fig. G) ####

par(mar = c(2.5, 2.5, 2, 1))  # Bottom, Left, Top, Right
par(mgp = c(3, 0.3, 0))     # Axis label and tick distance
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
x0 <- 0.0  # start x for line
x1 <- 0.07  # end x for line
y_start <- 0.7
line_spacing <- 0.07

# Draw each line, label, and shaded polygon
for (i in seq_len(n)) {
  y <- y_start - (i - 1) * line_spacing

  # Polygon around line (±0.01 vertically for shading)
  polygon(
    x = c(x0-0.01, x1+0.01, x1+0.01, x0-0.01),
    y = c(y - 0.02, y - 0.02, y + 0.02, y + 0.02),
    col = adjustcolor(legend_colors[i], alpha.f = 0.25),
    border = NA
  )
  
  # Line
  segments(x0, y, x1, y, col = legend_colors[i], lwd = 1.5)
  
  # Text
  text(x = x1 + 0.04, y = y, labels = legend_labels[i], adj = 0, cex = 0.56)
}




dev.off()

#
