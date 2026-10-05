# ==============================================================================
# STANDARD OXYGEN DISSOCIATION CURVE WITH O2 LIGATION AT THE INFLECTION POINT
# ==============================================================================
# Author: Dr. Holger H. Burchert (holger.burchert@unibas.ch)
# Affiliation: Department of Sport, Exercise and Health - University of Basel, Switzerland

# Main idea:
# This is the code (SHbO2CO2_Dash_2016_Adair.R) from Burchert (2025), producing 
# the standard oxygen dissocation curve with the hemoglobin ligation states at the 
# inflection point. Only some labels where changed for clarity. 

# References:
# Holger Burchert (DPhil). (2025). Oxygen dissociation curve inflection point 
# during incremental exercise: a trigger for the Bohr effect (Data & Analysis Code) 
# (Version v1.3) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.15614693


library(ggplot2)
library(tidyr)
library(pracma)
library(cowplot)

source("SHbO2CO2_Dash_2016_Adair.R", local = TRUE)


# CALCULATE O2 EQUILIBRIUM CURVE AND O2 LIGATION-STATE CURVES
# ------------------------------------------------------------------------------
states    <- c("Hb4", "Hb4O2_1", "Hb4O2_2", "Hb4O2_3", "Hb4O2_4")
PO2_sweep <- seq(0, 100, by = 0.1)

# Calculate O2 ligation states and HbO2 saturation for each PO2 value
rows <- lapply(PO2_sweep, function(p) {
  
  # Get model outputs for current PO2
  r <- SHbO2CO2(PO2 = p)
  
  # Return ligation-state percentages and HbO2 saturation
  data.frame(
    PO2     = p,
    Hb4     = 100 * r$Hb4,
    Hb4O2_1 = 100 * r$Hb4O2_1,
    Hb4O2_2 = 100 * r$Hb4O2_2,
    Hb4O2_3 = 100 * r$Hb4O2_3,
    Hb4O2_4 = 100 * r$Hb4O2_4,
    Mean    = 100 * r$SHbO2
  )
})

# Combine rows into one data frame
OEC_df <- do.call(rbind, rows)
rownames(OEC_df) <- NULL


# DETERMINE OEC INFLECTION POINT
# ------------------------------------------------------------------------------

# Inflection point = PO2 at which the slope of the saturation curve is maximal
dS          <- gradient(OEC_df$Mean, PO2_sweep)
i_inflect   <- which.max(dS)
ObsPO2      <- OEC_df$PO2[i_inflect]

# Reshape ligation-state data from wide to long format
OEC_long <- pivot_longer(
  OEC_df,
  cols      = all_of(states),
  names_to  = "State",
  values_to = "Percent"
)


# O2 LIGATION HISTOGRAM AT THE INFLECTION POINT
# ------------------------------------------------------------------------------

panel_max   <- 100
bin_width   <- 7.5
block_left  <- panel_max - length(states) * bin_width
bin_centers <- block_left + bin_width / 2 +
  (seq_along(states) - 1) * bin_width

names(bin_centers) <- states

# Extract ligation-state percentages at the inflection PO2
i_obs    <- which.min(abs(OEC_df$PO2 - ObsPO2))
hist_row <- OEC_df[i_obs, states, drop = FALSE]

ligation_bins <- data.frame(
  State   = states,
  Percent = as.numeric(hist_row[1, ]),
  xpos    = bin_centers[states],
  row.names = NULL
)

# Exact HbO2 saturation at the inflection point
Sat_exact <- SHbO2CO2(PO2 = ObsPO2)$SHbO2 * 100

annot <- data.frame(
  xpos  = min(ObsPO2 + 2, 99),
  ypos  = 95,
  label = sprintf("(%.1f, %.1f)", ObsPO2, Sat_exact)
)


# CREATE PLOT
# ------------------------------------------------------------------------------

standard_OEC <-
  ggplot() +
  
  # O2 ligation-state curves
  geom_line(
    data = OEC_long,
    aes(x = PO2, y = Percent, color = State),
    linewidth = 0.75
  ) +
  
  # HbO2 saturation curve
  geom_line(
    data = OEC_df,
    aes(x = PO2, y = Mean, color = "SHbO2"),
    linewidth = 0.75
  ) +
  
  # Vertical line at OEC inflection point
  geom_vline(
    xintercept = ObsPO2,
    linetype   = "dashed",
    color      = "black",
    linewidth  = 0.4
  ) +
  # Inflection point
  geom_point(
    aes(x = ObsPO2, y = Sat_exact),
    color = "black",
    size  = 2.5,
    inherit.aes = FALSE
  ) +
  
  # Ligation-state histogram at the inflection point
  geom_col(
    data = ligation_bins,
    aes(x = xpos, y = Percent, fill = State),
    width       = bin_width,
    inherit.aes = FALSE
  ) +
  
  # Histogram percentages
  geom_text(
    data = ligation_bins,
    aes(
      x = xpos,
      y = Percent + 5,
      label = sprintf("%.0f", Percent)
    ),
    size  = 2.5,
    vjust = 0
  ) +
  
  # Inflection-point annotation
  geom_text(
    data = annot,
    aes(x = xpos, y = ypos, label = label),
    hjust = 0,
    vjust = 0.5,
    size  = 3
  ) +
  
  # Curve colors and legend labels
  scale_color_manual(
    name = NULL,
    values = c(
      Hb4     = "#E64B35FF",
      Hb4O2_1 = "#4DBBD5FF",
      Hb4O2_2 = "#00A087FF",
      Hb4O2_3 = "#3C5488FF",
      Hb4O2_4 = "#F39B7FFF",
      SHbO2   = "#444444"
    ),
    labels = c(
      Hb4     = expression(0~O[2]~bound),
      Hb4O2_1 = expression(1~O[2]~bound),
      Hb4O2_2 = expression(2~O[2]~bound),
      Hb4O2_3 = expression(3~O[2]~bound),
      Hb4O2_4 = expression(4~O[2]~bound),
      SHbO2   = expression(HbO[2]~Saturation)
    ),
    breaks = c(
      "SHbO2",
      "Hb4",
      "Hb4O2_1",
      "Hb4O2_2",
      "Hb4O2_3",
      "Hb4O2_4"
    )
  ) +
  
  guides(
    color = guide_legend(
      nrow = 2,
      byrow = TRUE,
      override.aes = list(linewidth = 1.6),
      keywidth  = grid::unit(0.6, "lines"),
      keyheight = grid::unit(0.9, "lines")
    )
  ) +
  
  # Histogram colors
  scale_fill_manual(
    values = c(
      Hb4     = "#E64B35FF",
      Hb4O2_1 = "#4DBBD5FF",
      Hb4O2_2 = "#00A087FF",
      Hb4O2_3 = "#3C5488FF",
      Hb4O2_4 = "#F39B7FFF"
    ),
    guide = "none"
  ) +
  
  # Axes
  scale_x_continuous(
    breaks = seq(0, 100, 20),
    limits = c(0, 100),
    name   = expression(PO[2]~"(mmHg)")
  ) +
  
  scale_y_continuous(
    breaks = seq(0, 100, 20),
    name   = expression(HbO[2]~"Saturation & Hb Molecules (%)")
  )+
  
  # Theme
  theme_minimal() +
  theme(
    legend.position      = "top",
    legend.justification = "center",
    legend.title         = element_blank(),
    
    # Compact legend spacing so the longer labels fit
    legend.spacing.x     = grid::unit(0.02, "cm"),
    legend.key.width     = grid::unit(0.6, "lines"),
    
    panel.grid.major     = element_line(color = "gray80", linewidth = 0.2),
    panel.grid.minor     = element_line(color = "gray90", linewidth = 0.2)
  )

print(standard_OEC)


# PANEL B: HbO2 SATURATION CURVE WITH SWAPPED AXES
# ------------------------------------------------------------------------------
# Annotation for swapped axes
annot_swapped <- data.frame(
  xpos  = Sat_exact + 8,
  ypos  = ObsPO2 - 4,
  label = sprintf("(%.1f, %.1f)", Sat_exact, ObsPO2)
)
# Dummy data used only to create an invisible second legend row
legend_spacer <- data.frame(
  Mean = c(0, 0),
  PO2  = c(0, 0)
)

swapped_OEC <-
  ggplot() +
  
  # HbO2 saturation curve only
  geom_line(
    data = OEC_df,
    aes(x = Mean, y = PO2, color = "SHbO2"),
    linewidth = 0.75
  ) +
  
  # Invisible dummy legend entry
  # This forces Panel B's legend to occupy two rows,
  # matching the legend height of Panel A
  geom_line(
    data = legend_spacer,
    aes(x = Mean, y = PO2, color = "Spacer"),
    linewidth = 0.75,
    alpha = 0,
    show.legend = TRUE
  ) +
  
  # Inflection point
  geom_hline(
    yintercept = ObsPO2,
    linetype   = "dashed",
    color      = "black",
    linewidth  = 0.4
  ) +
  # Inflection point
  geom_point(
    aes(x = Sat_exact, y = ObsPO2),
    color = "black",
    size  = 2.5,
    inherit.aes = FALSE
  ) +
  
  # Inflection-point annotation
  geom_text(
    data = annot_swapped,
    aes(x = xpos, y = ypos, label = label),
    hjust = 0,
    vjust = 0.5,
    size  = 3
  ) +
  
  # HbO2 saturation curve color and invisible spacer
  scale_color_manual(
    name = NULL,
    values = c(
      SHbO2  = "#444444",
      Spacer = "#444444"
    ),
    labels = c(
      SHbO2  = expression(HbO[2]~Saturation),
      Spacer = ""
    ),
    breaks = c("SHbO2", "Spacer")
  ) +
  
  # Force two legend rows:
  # row 1 = HbO2 Saturation
  # row 2 = invisible spacer
  guides(
    color = guide_legend(
      nrow  = 2,
      byrow = TRUE,
      override.aes = list(linewidth = 1.6),
      keywidth  = grid::unit(0.6, "lines"),
      keyheight = grid::unit(0.9, "lines")
    )
  ) +
  
  # HbO2 saturation on x-axis
  scale_x_continuous(
    breaks = seq(0, 100, 20),
    limits = c(0, 100),
    name   = expression(HbO[2]~"Saturation (%)")
  ) +
  
  # PO2 on y-axis
  scale_y_continuous(
    breaks = seq(0, 100, 20),
    limits = c(0, 100),
    name   = expression(PO[2]~"(mmHg)")
  ) +
  
  theme_minimal() +
  
  theme(
    legend.position      = "top",
    legend.justification = "center",
    legend.title         = element_blank(),
    
    legend.spacing.x     = grid::unit(0.02, "cm"),
    legend.key.width     = grid::unit(0.6, "lines"),
    
    panel.grid.major = element_line(
      color = "gray80",
      linewidth = 0.2
    ),
    
    panel.grid.minor = element_line(
      color = "gray90",
      linewidth = 0.2
    )
  )


# COMBINE PANEL A AND PANEL B
# ------------------------------------------------------------------------------

two_panel_OEC <- cowplot::plot_grid(
  standard_OEC,
  swapped_OEC,
  ncol = 2,
  align = "h",
  axis = "tb",
  rel_widths = c(1, 1),
  labels = c("A", "B"),
  label_size = 12,
  label_fontface = "bold"
)

print(two_panel_OEC)


# EXPORT
# ------------------------------------------------------------------------------
# PDF
ggsave(
  filename = "OEC_ligation_and_swapped_axis.pdf",
  plot     = two_panel_OEC,
  device   = cairo_pdf,
  width    = 17.8,
  height   = 8.9,
  units    = "cm"
)

# High-resolution TIFF
ggsave(
  filename    = "OEC_ligation_and_swapped_axis.tiff",
  plot        = two_panel_OEC,
  device      = "tiff",
  width       = 17.8,
  height      = 8.9,
  units       = "cm",
  dpi         = 600,
  compression = "lzw"
)