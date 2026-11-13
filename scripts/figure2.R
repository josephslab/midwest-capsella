# Generate Figure 2: Trait Data Points with Model LS Means and SE Error Bars

library(dplyr)
library(ggplot2)
library(glmmTMB)
library(emmeans)
library(gridExtra)

# Resolve paths whether run from project root or scripts/ directory
data_path <- if (file.exists("data/1stHarvRformatedSeeds_Nov27.csv")) {
  "data/1stHarvRformatedSeeds_Nov27.csv"
} else {
  "../data/1stHarvRformatedSeeds_Nov27.csv"
}

fig_dir <- if (dir.exists("figures")) "figures" else "../figures"

# 1. Load and prepare data ------------------------------------------------
Harv1Individuals <- read.csv(data_path)
Harv1Individuals$City <- as.factor(Harv1Individuals$City)

Harv1Individuals <- Harv1Individuals %>%
  mutate(
    Fecundity = (Seed_count / 5) * N_Fruits,
    across(c(HavHeight_Scaled, Rosbranches_Scaled, Caulibranches_Scaled), ~ as.vector(scale(.))),
    BushyScore = (Rosbranches_Scaled - Caulibranches_Scaled - HavHeight_Scaled) / 2
  )

Harv1GCMidWestInd <- Harv1Individuals %>%
  filter(GrowingConditions == "GC", City != "New York") %>%
  droplevels()

measured_traits <- c(
  "Rosette_D", "Height_Flower", "Height_Harvest", "N_Fruits",
  "Rosette_branches", "Cauline_branches", "DaysFridgetoGerm",
  "DaysPlantToBolt", "DaysBoltToFlower", "Seed_count", "Fecundity", "BushyScore"
)

remove_outliers_3sd <- function(x) {
  x <- as.numeric(x)
  z <- abs(x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  ifelse(!is.na(z) & z > 3.0, NA_real_, x)
}

Harv1GCMidWestInd <- Harv1GCMidWestInd %>%
  group_by(City) %>%
  mutate(across(all_of(measured_traits), remove_outliers_3sd)) %>%
  ungroup() %>%
  mutate(
    FlowerTime = DaysPlantToBolt + DaysBoltToFlower,
    across(c(City, Genotype, ParentGrowingConditions, F2_Rep), as.factor)
  )

# 2. Fit glmmTMB models for the 4 Figure 2 traits -------------------------
gamma_disp <- ~ City + ParentGrowingConditions

mod_Bolt <- glmmTMB(
  DaysPlantToBolt ~ ParentGrowingConditions + City + F2_Rep + (1 | Genotype),
  dispformula = gamma_disp,
  family      = Gamma(link = "log"),
  data        = Harv1GCMidWestInd,
  REML        = FALSE
)

mod_Rosette_D <- glmmTMB(
  Rosette_D ~ ParentGrowingConditions + City + F2_Rep + (1 | Genotype),
  family = gaussian(),
  data   = Harv1GCMidWestInd,
  REML   = FALSE
)

mod_FlowerTime <- glmmTMB(
  FlowerTime ~ ParentGrowingConditions + City + F2_Rep + (1 | Genotype),
  dispformula = gamma_disp,
  family      = Gamma(link = "log"),
  data        = Harv1GCMidWestInd,
  REML        = FALSE
)

mod_Fecundity <- glmmTMB(
  Fecundity ~ ParentGrowingConditions + City + F2_Rep + (1 | Genotype),
  family = gaussian(),
  data   = Harv1GCMidWestInd,
  REML   = FALSE
)

# 3. Plot helper: raw jittered data points + model LS means & SE ----------
city_labels <- c("Detroit", "Grand Rapids", "Kalamazoo", "Madison", "Milwaukee")
city_fills  <- c("blue", "#06d6a0", "purple", "gold", "#ff6f61")
city_shapes <- c(21, 24, 22, 23, 25)

plot_lsmeans_with_points <- function(data, y_var, model, y_label, plot_title) {
  # Obtain marginal means on the response scale so they match the raw trait values
  emm_df <- as.data.frame(emmeans(model, ~ City, type = "response"))
  if ("response" %in% names(emm_df)) {
    emm_df$emmean <- emm_df$response
  }

  ggplot() +
    geom_jitter(
      data   = data,
      aes(x = City, y = .data[[y_var]]),
      width  = 0.2,
      size   = 2,
      alpha  = 0.45,
      colour = "grey30"
    ) +
    geom_errorbar(
      data      = emm_df,
      aes(x = City, ymin = emmean - SE, ymax = emmean + SE),
      width     = 0.3,
      linewidth = 0.8,
      colour    = "black"
    ) +
    geom_point(
      data = emm_df,
      aes(x = City, y = emmean, shape = City, colour = City, fill = City),
      size = 7
    ) +
    scale_colour_manual(values = rep("black", 5), labels = city_labels) +
    scale_fill_manual(values = city_fills, labels = city_labels) +
    scale_shape_manual(values = city_shapes, labels = city_labels) +
    theme(
      panel.background = element_rect(fill = "white", colour = "white", linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position  = "none",
      axis.line        = element_line(linewidth = 0.5, colour = "black"),
      axis.text        = element_text(colour = "black", size = 20),
      axis.title       = element_text(size = 25),
      plot.title       = element_text(size = 30)
    ) +
    labs(x = "City", y = y_label, title = plot_title)
}

DaysPlantToBolt_plot <- plot_lsmeans_with_points(
  Harv1GCMidWestInd, "DaysPlantToBolt", mod_Bolt, "Bolting time (days)", "A"
)
Rosette_plot <- plot_lsmeans_with_points(
  Harv1GCMidWestInd, "Rosette_D", mod_Rosette_D, "Rosette diameter (cm)", "B"
)
Flowering_plot <- plot_lsmeans_with_points(
  Harv1GCMidWestInd, "FlowerTime", mod_FlowerTime, "Flowering time (days)", "C"
)
Fecundity_plot <- plot_lsmeans_with_points(
  Harv1GCMidWestInd, "Fecundity", mod_Fecundity, "Fecundity", "D"
)

combined_LSPLOTS <- arrangeGrob(
  DaysPlantToBolt_plot, Rosette_plot, Flowering_plot, Fecundity_plot,
  ncol = 2, padding = unit(2, "cm")
)

ggsave(
  filename = file.path(fig_dir, "figure2_lsmeans_with_points.png"),
  plot     = combined_LSPLOTS,
  width    = 20,
  height   = 16,
  dpi      = 600
)
