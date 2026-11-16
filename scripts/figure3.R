# Generate Multipanel Figure: Mean Pi (Left) and Genomic PCA (Right) with Shared Bottom Legend

library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)

# Resolve paths whether run from project root or scripts/ directory
data_dir <- if (dir.exists("data")) "data" else "../data"
fig_dir  <- if (dir.exists("figures")) "figures" else "../figures"

city_labels <- c("Detroit", "Grand Rapids", "Kalamazoo", "Madison", "Milwaukee")
city_fills  <- c("blue", "#06d6a0", "purple", "gold", "#ff6f61")
city_shapes <- c(21, 24, 22, 23, 25)

# =========================================================================
# 1. Left Panel: Mean Pi and 95% Bootstrap CIs (from chrissy_pi_4fold.R)
# =========================================================================
pi <- read.delim(file.path(data_dir, "Midwest_4-fold_10kb_pi.txt"))

set.seed(123)
n_boot <- 1000

bootstrap_city_pi <- function(df) {
  mysites <- unique(df$window_pos_1)
  df$window_pos_1 <- as.numeric(df$window_pos_1)

  boot_vals <- replicate(n_boot, {
    sampled_sites <- sample(mysites, replace = TRUE, size = length(mysites))
    boot_sample   <- df %>% filter(window_pos_1 %in% sampled_sites)
    if (nrow(boot_sample) > 0) {
      mean(boot_sample$avg_pi, na.rm = TRUE)
    } else {
      NA_real_
    }
  })
  boot_vals <- boot_vals[!is.na(boot_vals)]

  df %>%
    group_by(pop) %>%
    summarise(
      mean_pi   = mean(avg_pi, na.rm = TRUE),
      median_pi = median(avg_pi, na.rm = TRUE),
      lower_CI  = quantile(boot_vals, 0.025, na.rm = TRUE),
      upper_CI  = quantile(boot_vals, 0.975, na.rm = TRUE),
      .groups   = "drop"
    )
}

Detroit <- subset(pi, pop == "D")
GR      <- subset(pi, pop == "GR")
Kal     <- subset(pi, pop == "Kal")
Mil     <- subset(pi, pop == "Mil")
Mad     <- subset(pi, pop == "Mad")

Cities_Pi_CI <- bind_rows(
  bootstrap_city_pi(Detroit),
  bootstrap_city_pi(GR),
  bootstrap_city_pi(Kal),
  bootstrap_city_pi(Mil),
  bootstrap_city_pi(Mad)
)

Cities_Pi_CI$pop <- c("Detroit", "Grand Rapids", "Kalamazoo", "Milwaukee", "Madison")
Cities_Pi_CI_reordered <- Cities_Pi_CI %>%
  arrange(pop) %>%
  mutate(pop = factor(pop, levels = city_labels))

MeanCI_plot <- ggplot(Cities_Pi_CI_reordered, aes(x = pop, y = mean_pi, fill = pop)) +
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.2, linewidth = 1) +
  geom_point(aes(shape = pop, colour = pop, fill = pop), size = 7) +
  scale_colour_manual(name = "City", values = rep("black", 5), labels = city_labels) +
  scale_fill_manual(name = "City", values = city_fills, labels = city_labels) +
  scale_shape_manual(name = "City", values = city_shapes, labels = city_labels) +
  theme(
    panel.background = element_rect(fill = "white", colour = "white", linewidth = 0.5, linetype = "solid"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position  = "none",
    axis.line        = element_line(linewidth = 0.5, colour = "black"),
    axis.text.x      = element_text(colour = "black", size = 18, angle = 25, hjust = 1),
    axis.text.y      = element_text(colour = "black", size = 20),
    axis.title       = element_text(size = 25),
    plot.title       = element_text(size = 30)
  ) +
  labs(y = expression("Mean " * pi), x = "City", title = "A")

# =========================================================================
# 2. Right Panel: PCA Plot (from PCA-plot-Fig5.R)
# =========================================================================
midwestEv <- read.table(file.path(data_dir, "midwestCapsella-pca.eigenvec"))
names(midwestEv)[2] <- "id"

midwestPopInfo <- read.table(file.path(data_dir, "midwest_pixy_pops.txt"))
names(midwestPopInfo) <- c("id", "pop")
midwestPopInfo$pop <- factor(
  midwestPopInfo$pop,
  levels = c("D", "GR", "Kal", "Mad", "Mil"),
  labels = city_labels
)

midwestEv <- dplyr::left_join(midwestEv, midwestPopInfo, by = "id")

system("awk '{print $NR}' ../data/midwestCapsella-pca.rel > ../data/midwestCapsella-pca.rel.diag")
mydiag = read.table('../data/midwestCapsella-pca.rel.diag')
totalVar = sum(mydiag$V1)
evals    <- read.table(file.path(data_dir, "midwestCapsella-pca.eigenval"))
evals$perc <- evals$V1 / totalVar

scatterV3V4_with_legend <- ggplot(midwestEv, aes(x = V3, y = V4)) +
  geom_point(aes(shape = pop, colour = pop, fill = pop), size = 7) +
  scale_colour_manual(name = "City", values = rep("black", 5), labels = city_labels) +
  scale_fill_manual(name = "City", values = city_fills, labels = city_labels) +
  scale_shape_manual(name = "City", values = city_shapes, labels = city_labels) +
  theme(
    panel.background  = element_rect(fill = "white", colour = "white", linewidth = 0.5, linetype = "solid"),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    legend.position   = "bottom",
    legend.direction  = "horizontal",
    legend.background = element_rect(colour = NA, fill = "white"),
    legend.key        = element_rect(colour = "white", fill = NA),
    legend.key.width  = unit(1.2, "cm"),
    legend.spacing.x  = unit(0.6, "cm"),
    legend.text       = element_text(size = 20, margin = margin(r = 25, unit = "pt")),
    legend.title      = element_text(size = 22, margin = margin(r = 20, unit = "pt")),
    axis.line         = element_line(linewidth = 0.5, colour = "black"),
    axis.text         = element_text(colour = "black", size = 20),
    axis.title        = element_text(size = 25),
    plot.title        = element_text(size = 30)
  ) +
  guides(
    fill   = guide_legend(nrow = 1),
    shape  = guide_legend(nrow = 1),
    colour = guide_legend(nrow = 1)
  ) +
  labs(
    y     = paste0("PC 2 (", round(evals$perc[2], digits = 4) * 100, "%)"),
    x     = paste0("PC 1 (", round(evals$perc[1], digits = 4) * 100, "%)"),
    title = "B"
  )

# =========================================================================
# 3. Extract Shared Bottom Legend and Combine Panels (Aligned X-Axes)
# =========================================================================
extract_legend <- function(a_gplot) {
  tmp <- ggplot_gtable(ggplot_build(a_gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  tmp$grobs[[leg]]
}

shared_legend <- extract_legend(scatterV3V4_with_legend)
scatterV3V4   <- scatterV3V4_with_legend + theme(legend.position = "none")

# Convert both plots to gtables and equalize row heights so x-axes align exactly
g_left  <- ggplotGrob(MeanCI_plot + theme(plot.margin = margin(15, 25, 15, 15, "pt")))
g_right <- ggplotGrob(scatterV3V4 + theme(plot.margin = margin(15, 15, 15, 25, "pt")))

max_heights     <- grid::unit.pmax(g_left$heights, g_right$heights)
g_left$heights  <- max_heights
g_right$heights <- max_heights

top_panels <- arrangeGrob(
  g_left, g_right,
  ncol = 2
)

combined_pi_pca <- arrangeGrob(
  top_panels,
  shared_legend,
  nrow = 2,
  heights = c(10, 1.2)
)

ggsave(
  filename = file.path(fig_dir, "figure3_pi_and_pca.png"),
  plot     = combined_pi_pca,
  width    = 20,
  height   = 10,
  dpi      = 600,
  bg       = "white"
)

