###############################################################################
# FINAL PIPELINE — FIGURE 1 + FIGURE 2
# FIXES:
# - CDH correctly inverted (cooling_relief)
# - FIG 2 quadrant labels moved OUT of data space (no overlap)
# - extra left margin added
# - plots shown before saving
###############################################################################

library(tidyverse)
library(readxl)
library(ggrepel)
library(scales)

theme_set(theme_classic(base_size = 14))

# -------------------------------
# 0. LOAD + CLEAN
# -------------------------------
dat <- read_excel("Site_Summary_Master_OP.xlsx",
                  sheet = "Site_Year_OP") %>%
  rename_with(~ gsub(" ", "_", .)) %>%
  mutate(
    State = toupper(State),
    Mean_HDH_daily = as.numeric(Mean_HDH_daily),
    Mean_CDH_daily = as.numeric(Mean_CDH_daily),
    `%_abv_thres`  = as.numeric(`%_abv_thres`)
  ) %>%
  filter(!is.na(State))

# -------------------------------
# 1. SITE-LEVEL SUMMARY
# -------------------------------
site_level <- dat %>%
  group_by(State, Site_Name) %>%
  summarise(
    mean_HDH = mean(Mean_HDH_daily, na.rm = TRUE),
    mean_CDH = mean(Mean_CDH_daily, na.rm = TRUE),
    pct_hot  = mean(`%_abv_thres`, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    cooling_relief = -mean_CDH   # 🔥 FIXED INTERPRETATION
  )

# -------------------------------
# 2. GLOBAL k-MEANS SPLITS
# -------------------------------
kmeans_midpoint <- function(x){
  x <- x[is.finite(x)]
  if(length(x) < 2) return(NA_real_)
  km <- kmeans(x, centers = 2)
  mean(km$centers)
}

global_hdh_split  <- kmeans_midpoint(site_level$mean_HDH)
global_cool_split <- kmeans_midpoint(site_level$cooling_relief)

# -------------------------------
# 3. STATE RISK (75th percentile)
# -------------------------------
site_level <- site_level %>%
  group_by(State) %>%
  mutate(
    n_sites = n(),
    
    hdh_thresh = ifelse(n_sites >= 2,
                        quantile(mean_HDH, 0.75, na.rm = TRUE),
                        NA_real_),
    
    cdh_thresh = ifelse(n_sites >= 2,
                        quantile(mean_CDH, 0.75, na.rm = TRUE),
                        NA_real_),
    
    High_Stress = mean_HDH >= hdh_thresh,
    Low_Cooling = mean_CDH >= cdh_thresh,
    
    At_Risk = High_Stress & Low_Cooling
  ) %>%
  ungroup()

# =============================================================================
# FIGURE 1 — THERMAL RISK MAP
# =============================================================================
p1 <- ggplot(site_level,
             aes(x = mean_HDH,
                 y = cooling_relief,
                 color = State)) +
  
  geom_point(size = 3, alpha = 0.85) +
  
  geom_point(
    data = subset(site_level, At_Risk),
    color = "red",
    size = 4
  ) +
  
  geom_text_repel(
    aes(label = Site_Name),
    size = 3.5,
    box.padding = 0.5,
    point.padding = 0.3,
    max.overlaps = 100
  ) +
  
  labs(
    title = "Thermal Risk Across Genomic Sampling Sites",
    subtitle = "Red = high stress + reduced cooling (state-relative 75th percentile)",
    x = "Thermal Stress (Mean HDH)",
    y = "Cooling Relief (−CDH; higher = more cooling)",
    color = "State"
  )

# SHOW FIRST
print(p1)

ggsave("FIG1_thermal_risk.png", p1, width = 11, height = 8)

# =============================================================================
# FIGURE 2 — GLOBAL THERMAL REGIMES (FIXED LABEL OVERLAP)
# =============================================================================
p2 <- ggplot(site_level,
             aes(x = mean_HDH,
                 y = cooling_relief,
                 color = State)) +
  
  geom_point(size = 3, alpha = 0.85) +
  
  geom_vline(xintercept = global_hdh_split, linetype = "dashed") +
  geom_hline(yintercept = global_cool_split, linetype = "dashed") +
  
  geom_text_repel(
    aes(label = Site_Name),
    size = 3.5,
    box.padding = 0.5,
    point.padding = 0.3,
    max.overlaps = 100
  ) +
  
  # --------------------------------------------------------------------------
# 🔥 FIXED QUADRANT LABELS (NO OVERLAP USING Inf ANCHORS)
# --------------------------------------------------------------------------

annotate("text",
         x = -Inf, y = Inf,
         label = "Low Stress / High Cooling\nStable Systems",
         hjust = -0.1, vjust = 1.2,
         size = 4.5, fontface = "bold") +
  
  annotate("text",
           x = Inf, y = Inf,
           label = "High Stress / High Cooling\nBuffered Systems",
           hjust = 1.1, vjust = 1.2,
           size = 4.5, fontface = "bold") +
  
  annotate("text",
           x = -Inf, y = -Inf,
           label = "Low Stress / Low Cooling\nThermal Refugia",
           hjust = -0.1, vjust = -0.2,
           size = 4.5, fontface = "bold") +
  
  annotate("text",
           x = Inf, y = -Inf,
           label = "High Stress / Low Cooling\nChronic Heat Stress",
           hjust = 1.1, vjust = -0.2,
           size = 4.5, fontface = "bold") +
  
  labs(
    title = "Thermal Stress vs Cooling Relief (System Regimes)",
    subtitle = "Quadrants defined by global k-means splits",
    x = "Thermal Stress (Mean HDH)",
    y = "Cooling Relief (−CDH; higher = more cooling)",
    color = "State"
  ) +
  
  # -------------------------------
# SPACING FIX (LEFT SIDE)
# -------------------------------
scale_x_continuous(
  expand = expansion(mult = c(0.25, 0.05))
) +
  coord_cartesian(clip = "off") +
  theme(
    plot.margin = margin(20, 60, 20, 20)
  )

# SHOW FIRST
print(p2)

ggsave("FIG2_stress_vs_cooling.png", p2, width = 12, height = 9)

###############################################################################
# END
###############################################################################

eelgrass_df <- dat %>%
  group_by(State, Site_Name) %>%
  summarise(
    mean_events  = mean(Number_warming_events, na.rm = TRUE),
    mean_exposure = mean(`%_abv_thres`, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    eelgrass_loss = mean_exposure > 18.5
  )

p_eelgrass <- ggplot(eelgrass_df,
                     aes(x = mean_events,
                         y = mean_exposure,
                         color = State)) +
  
  geom_point(size = 3.5, alpha = 0.85) +
  
  # EELGRASS LOSS THRESHOLD
  geom_hline(yintercept = 18.5,
             linetype = "dashed",
             color = "darkgreen",
             linewidth = 1) +
  
  geom_text_repel(
    aes(label = Site_Name),
    size = 3.5,
    max.overlaps = 100
  ) +
  
  # highlight loss sites
  geom_point(
    data = subset(eelgrass_df, eelgrass_loss),
    color = "red",
    size = 4
  ) +
  
  labs(
    title = "Thermal Disturbance and Eelgrass Collapse Threshold",
    subtitle = "Red = sites exceeding 18.5% time above thermal threshold (eelgrass loss zone)",
    x = "Mean Number of Warming Events (Disturbance Frequency)",
    y = "Mean % Time Above Thermal Threshold (Exposure Duration)",
    color = "State"
  ) +
  
  theme_classic(base_size = 14)

# 👇 SHOW BEFORE SAVING
print(p_eelgrass)

ggsave("FIG_eelgrass_collapse_threshold.png",
       p_eelgrass,
       width = 11,
       height = 8)
###############################################################################
# 3-PANEL SPATIAL ECOLOGY FIGURE (CORRECTED INTERPRETATION)
#
# A = Spatial thermal exposure
# B = Latitudinal gradient
# C = Threshold exceedance (NOT eelgrass loss)
###############################################################################

library(tidyverse)
library(readxl)
library(ggrepel)
library(patchwork)

theme_set(theme_classic(base_size = 14))

# -------------------------------
# LOAD DATA
# -------------------------------
dat <- read_excel("Site_Summary_Master_OP.xlsx",
                  sheet = "Site_Year_OP") %>%
  rename_with(~ gsub(" ", "_", .)) %>%
  mutate(
    State = toupper(State),
    Mean_HDH_daily = as.numeric(Mean_HDH_daily),
    Mean_CDH_daily = as.numeric(Mean_CDH_daily),
    `%_abv_thres`  = as.numeric(`%_abv_thres`),
    Lat  = as.numeric(Lat),
    Long = as.numeric(Long)
  ) %>%
  filter(!is.na(State))

# -------------------------------
# SITE-LEVEL SUMMARY
# -------------------------------
site_level <- dat %>%
  group_by(State, Site_Name, Lat, Long) %>%
  summarise(
    mean_HDH = mean(Mean_HDH_daily, na.rm = TRUE),
    mean_CDH = mean(Mean_CDH_daily, na.rm = TRUE),
    pct_hot  = mean(`%_abv_thres`, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    cooling_relief = -mean_CDH,
    high_exposure = pct_hot > 18   # ecological threshold
  )

# =============================================================================
# PANEL A — SPATIAL THERMAL EXPOSURE (FIXED ZOOM)
# =============================================================================
pA <- ggplot(site_level,
             aes(x = Long, y = Lat)) +
  
  borders("state", colour = "grey80", fill = "grey95") +
  
  geom_point(aes(color = pct_hot,
                 size = mean_HDH),
             alpha = 0.85) +
  
  scale_color_viridis_c(option = "plasma") +
  
  geom_text_repel(aes(label = Site_Name),
                  size = 3,
                  max.overlaps = 100) +
  
  coord_cartesian(
    xlim = c(-77, -66),
    ylim = c(34, 45),
    clip = "off"
  ) +
  
  labs(
    title = "A. Spatial Distribution of Thermal Exposure",
    subtitle = "% time above thermal threshold across coastal sites",
    x = "Longitude",
    y = "Latitude",
    color = "% Time Above Threshold",
    size = "Thermal Stress (HDH)"
  )

# =============================================================================
# PANEL B — LATITUDINAL GRADIENT
# =============================================================================
pB <- ggplot(site_level,
             aes(x = Lat, y = mean_HDH, color = State)) +
  
  geom_point(size = 3) +
  
  geom_smooth(method = "loess", se = FALSE, color = "black") +
  
  labs(
    title = "B. Latitudinal Gradient in Thermal Stress",
    x = "Latitude",
    y = "Thermal Stress (Mean HDH)"
  )

# =============================================================================
# PANEL C — THRESHOLD EXCEEDANCE MAP (>18%)
# =============================================================================
pC <- ggplot(site_level,
             aes(x = Long, y = Lat)) +
  
  borders("state", colour = "grey80", fill = "grey95") +
  
  geom_point(aes(color = high_exposure,
                 size = pct_hot),
             alpha = 0.85) +
  
  scale_color_manual(values = c("FALSE" = "grey60",
                                "TRUE" = "red")) +
  
  geom_text_repel(aes(label = Site_Name),
                  size = 3,
                  max.overlaps = 100) +
  
  coord_cartesian(
    xlim = c(-77, -66),
    ylim = c(34, 45),
    clip = "off"
  ) +
  
  labs(
    title = "C. Sites Exceeding Thermal Exposure Threshold",
    subtitle = ">18% time above thermal threshold (ecological benchmark)",
    x = "Longitude",
    y = "Latitude",
    color = "High Exposure",
    size = "% Time Above Threshold"
  )

# -------------------------------
# COMBINE
# -------------------------------
final_fig <- pA + pB + pC +
  plot_layout(ncol = 1)

# -------------------------------
# SHOW BEFORE SAVE
# -------------------------------
print(final_fig)

ggsave("FIG_spatial_3panel_corrected.png",
       final_fig,
       width = 11,
       height = 14,
       dpi = 300)

###############################################################################
# END
###############################################################################
