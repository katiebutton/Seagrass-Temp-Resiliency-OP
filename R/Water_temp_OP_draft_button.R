###############################################################################
# FINAL CLEAN PIPELINE
# FIG 0 = NC (CAHA vs CALO)
# FIG 1 = Thermal risk across all sites
# FIG 2 = Spatial + latitudinal + threshold panels
# FIG 3 = Eelgrass disturbance threshold
###############################################################################

library(tidyverse)
library(readxl)
library(ggrepel)
library(scales)
library(patchwork)
library(here)

theme_set(theme_classic(base_size = 14))

# -------------------------------
# 0. LOAD + CLEAN (RUN ONCE)
# -------------------------------
dat <- read_excel(here::here("Data","Site_Summary_Master_OP.xlsx"),
                  sheet = "Site_Year_OP") %>%
  rename_with(~ gsub(" ", "_", .)) %>%
  mutate(
    State = toupper(State),
    Mean_HDH_daily = as.numeric(Mean_HDH_daily),
    Mean_CDH_daily = as.numeric(Mean_CDH_daily),
    `%_abv_thres`  = as.numeric(`%_abv_thres`),
    Number_warming_events = as.numeric(Number_warming_events),
    Lat = as.numeric(Lat),
    Long = as.numeric(Long)
  ) %>%
  filter(!is.na(State))

# -------------------------------
# 1. SITE-LEVEL SUMMARY
# -------------------------------
site_level <- dat %>%
  group_by(Park_Code, Site_Name, State, Lat, Long) %>%
  summarise(
    mean_HDH = mean(Mean_HDH_daily, na.rm = TRUE),
    mean_CDH = mean(Mean_CDH_daily, na.rm = TRUE),
    pct_hot  = mean(`%_abv_thres`, na.rm = TRUE),
    mean_events = mean(Number_warming_events, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    cooling_relief = -mean_CDH,
    high_exposure = pct_hot > 18
  )

# -------------------------------
# 1B. NC SUBSET (CAHA vs CALO)
# -------------------------------
nc_sites <- site_level %>%
  filter(Park_Code %in% c("CAHA", "CALO")) %>%
  group_by(Park_Code) %>%
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

# -------------------------------
# 2. EELGRASS ANALYSIS
# -------------------------------
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

# -------------------------------
# 3. GLOBAL THRESHOLDS
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
# 4. STATE-LEVEL RISK METRICS
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
# FIGURE 0 — NC ONLY (CAHA vs CALO)
# =============================================================================
p0 <- ggplot(nc_sites,
             aes(x = mean_HDH,
                 y = cooling_relief,
                 color = Park_Code)) +
  
  geom_point(size = 3.5, alpha = 0.9) +
  
  # 🔴 PARK-LEVEL RISK FLAG (same logic as state-level)
  geom_point(
    data = subset(nc_sites, At_Risk),
    shape = 21,
    fill = "red",
    color = "black",
    size = 4,
    stroke = 1
  ) +
  
  geom_text_repel(
    aes(label = Site_Name),
    size = 3,
    max.overlaps = 100
  ) +
  
  labs(
    title = "North Carolina Thermal Regimes",
    subtitle = "CAHA vs CALO — red = park-relative high stress + low cooling",
    x = "Thermal Stress (Mean HDH)",
    y = "Cooling Relief (−CDH)",
    color = "Park"
  ) +
  
  theme_classic(base_size = 14)

print(p0)
ggsave("FIG0_CAHA_CALO.png", p2, width = 12, height = 9)
# =============================================================================
# FIGURE 1 — THERMAL RISK (ALL SITES)
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
    max.overlaps = 100
  ) +
  
  labs(
    title = "Thermal Risk Across Sampling Sites",
    subtitle = "Red = high stress + low cooling (75th percentile threshold)",
    x = "Thermal Stress (Mean HDH)",
    y = "Cooling Relief (−CDH)",
    color = "State"
  )

print(p1)

ggsave("FIG1_thermal_risk.png", p1, width = 11, height = 8)
# =============================================================================
# FIGURE 2 — GLOBAL THERMAL REGIMES
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
  
  # -------------------------------
# QUADRANT LABELS (FIXED)
# -------------------------------

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
    subtitle = "Quadrants defined by global splits",
    x = "Thermal Stress (Mean HDH)",
    y = "Cooling Relief (−CDH; higher = more cooling)",
    color = "State"
  ) +
  
  scale_x_continuous(
    expand = expansion(mult = c(0.25, 0.05))
  ) +
  
  coord_cartesian(clip = "off") +
  
  theme(
    plot.margin = margin(20, 60, 20, 20)
  )

# SHOW
print(p2)

ggsave("FIG2_stress_vs_cooling.png", p2, width = 12, height = 9)

# =============================================================================
# FIGURE 3 — SPATIAL + LATITUDINAL + THRESHOLD
# =============================================================================

# -------------------------------
# CLEAN COORDINATES (IMPORTANT)
# -------------------------------
site_level <- site_level %>%
  filter(
    !is.na(Long),
    !is.na(Lat),
    Long > -85, Long < -60,
    Lat > 30, Lat < 50
  )

# =============================================================================
# PANEL A — SPATIAL THERMAL EXPOSURE (FIXED)
# =============================================================================
pA <- ggplot(site_level, aes(Long, Lat)) +
  
  borders("state",
          fill = "grey95",
          color = "grey80") +
  
  geom_point(
    aes(color = pct_hot,
        size = mean_HDH),
    alpha = 0.85
  ) +
  
  scale_color_viridis_c(option = "plasma") +
  
  geom_text_repel(
    aes(label = Site_Name),
    size = 2.8,
    max.overlaps = 50
  ) +
  
  coord_cartesian(
    xlim = c(-80, -66),
    ylim = c(34, 45),
    expand= FALSE,
    clip = "on"   # 🔥 IMPORTANT FIX
  ) +
  
  labs(
    title = "A. Spatial Distribution of Thermal Exposure",
    subtitle = "East Coast sites",
    x = "Longitude",
    y = "Latitude",
    color = "% Time Above Threshold",
    size = "Thermal Stress (HDH)"
  ) +
  
  theme_classic(base_size = 14)

# =============================================================================
# PANEL B — LATITUDINAL GRADIENT
# =============================================================================
pB <- ggplot(site_level, aes(Lat, mean_HDH, color = State)) +
  
  geom_point(size = 3) +
  
  geom_smooth(method = "loess", se = FALSE, color = "black") +
  
  labs(
    title = "B. Latitudinal Gradient",
    x = "Latitude",
    y = "Thermal Stress (Mean HDH)"
  ) +
  
  theme_classic(base_size = 14)

# =============================================================================
# PANEL C — THRESHOLD EXCEEDANCE (MATCHED ZOOM)
# =============================================================================
pC <- ggplot(site_level, aes(Long, Lat)) +
  
  borders("state",
          fill = "grey95",
          color = "grey80") +
  
  geom_point(
    aes(color = high_exposure,
        size = pct_hot),
    alpha = 0.85
  ) +
  
  scale_color_manual(values = c("FALSE" = "grey60",
                                "TRUE" = "red")) +
  
  geom_text_repel(
    aes(label = Site_Name),
    size = 2.8,
    max.overlaps = 50
  ) +
  
  coord_cartesian(
    xlim = c(-80, -66),
    ylim = c(34, 45),
    expand= FALSE,
    clip = "on"
  ) +
  
  labs(
    title = "C. Threshold Exceedance",
    x = "Longitude",
    y = "Latitude",
    color = "High Exposure",
    size = "% Time Above Threshold"
  ) +
  
  theme_classic(base_size = 14)

# =============================================================================
# COMBINE
# =============================================================================
final_fig <- pA / pB / pC +
  plot_layout(heights = c(1.4, 1, 1.4))

# SHOW
print(final_fig)

# SAVE
ggsave("FIG3_spatial_3panel.png",
       final_fig,
       width = 13,
       height = 14,
       dpi = 300)

# =============================================================================
# FIGURE 4 — EELGRASS DISTURBANCE THRESHOLD
# =============================================================================
p_eelgrass <- ggplot(eelgrass_df,
                     aes(x = mean_events,
                         y = mean_exposure,
                         color = State)) +
  
  geom_point(size = 3.5, alpha = 0.85) +
  
  geom_hline(yintercept = 18.5,
             linetype = "dashed",
             color = "darkgreen",
             linewidth = 1) +
  
  geom_point(
    data = subset(eelgrass_df, eelgrass_loss),
    color = "red",
    size = 4
  ) +
  
  geom_text_repel(
    aes(label = Site_Name),
    size = 3.5,
    max.overlaps = 100
  ) +
  
  labs(
    title = "Thermal Disturbance and Eelgrass Threshold",
    subtitle = "Red = exceedance of 18.5% exposure threshold",
    x = "Warming Events",
    y = "% Time Above Thermal Threshold",
    color = "State"
  )

print(p_eelgrass)

ggsave("FIG4_eelgrass.png", p_eelgrass, width = 11, height = 8)

###############################################################################
# END
###############################################################################
