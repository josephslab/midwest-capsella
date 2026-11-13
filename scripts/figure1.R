# Generate Figure 1: Phenotypic Variance Boxplots Across Cities

library(dplyr)
library(ggplot2)
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
    City = as.factor(City)
  )

# 2. Plot helper and panel generation -------------------------------------
plot_variance_boxplot <- function(data, y_var, y_label, plot_title) {
  ggplot(data, aes(x = City, y = .data[[y_var]])) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 2) +
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

box_bolting   <- plot_variance_boxplot(Harv1GCMidWestInd, "DaysPlantToBolt", "Bolting time (days)", "A")
box_flowering <- plot_variance_boxplot(Harv1GCMidWestInd, "FlowerTime", "Flowering time (days)", "B")
box_fh        <- plot_variance_boxplot(Harv1GCMidWestInd, "Height_Flower", "Height at flowering (cm)", "C")
box_hh        <- plot_variance_boxplot(Harv1GCMidWestInd, "Height_Harvest", "Height at harvest (cm)", "D")

combined_Boxplots <- arrangeGrob(
  box_bolting, box_flowering, box_fh, box_hh,
  ncol = 2, padding = unit(2, "cm")
)

ggsave(
  filename = file.path(fig_dir, "figure1_variance_boxplots.png"),
  plot     = combined_Boxplots,
  width    = 20,
  height   = 16,
  dpi      = 600
)
