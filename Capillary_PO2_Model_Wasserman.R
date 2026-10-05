# ==============================================================================
# Wasserman'-style's capillary PO2 as function of fractional capillary distance
# ==============================================================================
# Author: Dr. Holger H. Burchert (holger.burchert@unibas.ch)
# Affiliation: Department of Sport, Exercise and Health - University of Basel, Switzerland

# Main idea:
# 1. Start with arterial PO2
# 2. Convert PO2 to saturation using Severinghaus
# 3. Convert saturation to O2 content
# 4. Let O2 content fall linearly along the capillary
# 5. Convert local saturation back to PO2
# 6. Plot PO2 vs fractional capillary distance
# 7. Estimate inflection points using numerical derivatives

# Simplifications:
# - No Bohr effect
# - No critical PO2 cutoff
# - Dissolved O2 ignored for simplicity

# References:
# Wasserman K. Coupling of external to cellular respiration during exercise: 
# the wisdom of the body revisited. Am J Physiol-Endocrinol Metab 266: E519–E539, 
# 1994. doi: 10.1152/ajpendo.1994.266.4.E519.

# Severinghaus JW. Simple, accurate equations for human blood O2 dissociation 
# computations. J Appl Physiol 46: 599–602, 1979. doi: 10.1152/jappl.1979.46.3.599.

# Ellis RK. Determination of PO2 from saturation. J Appl Physiol 67: 902–902, 1989. 
# doi: 10.1152/jappl.1989.67.2.902.

# REQUIRED LIBRARIES
# ------------------------------------------------------------------------------
library(pracma)
library(ggplot2)                             

# SEVERINGHAUS / Ellis EQUATIONS
# ------------------------------------------------------------------------------
# Convert PO2 to hemoglobin saturation
Severinghaus_S <- function(PO2) {
  (((PO2^3 + 150 * PO2)^-1 * 23400) + 1)^-1
}

# Convert hemoglobin saturation to PO2
Severinghaus_PO2 <- function(S) {
  A <- 11700 * (S^-1 - 1)^-1
  B <- (50^3 + A^2)^0.5
  PO2 <- (B + A)^(1/3) - (B - A)^(1/3)
  return(PO2)
}

# MODEL ASSUMPTIONS 
# ------------------------------------------------------------------------------
Hb                 <- 15
PaO2               <- 90
Qm_VO2m_ratio      <- c(9, 7, 6, 5)
capillary_distance <- seq(0, 1, length.out = 300)

SaO2 <- Severinghaus_S(PaO2)
CaO2 <- 1.34 * Hb * SaO2

# COLLECT DATA INTO DATA FRAMES 
# ------------------------------------------------------------------------------
curve_data  <- data.frame()
inflections <- data.frame(
  Qm_VO2m_ratio = numeric(),
  capillary_distance_inflection = numeric(),
  PO2_inflection = numeric()
)

# CALCULATE AND COLLECT EACH CURVE 
# ------------------------------------------------------------------------------
for (ratio in Qm_VO2m_ratio) {
  
  # Define the blood-flow-to-O2-consumption ratio:
  #   Qm_VO2m_ratio = Qm / VO2m

  # From Fick:
  #   VO2m = Qm * (CaO2 - CvO2)

  # Therefore:
  #   CaO2 - CvO2 = VO2m / Qm
  #               = 1 / (Qm / VO2m)
  #               = 1 / Qm_VO2m_ratio

  # This gives L O2 / L blood.
  # Convert to mL O2 / dL blood:
  #   1 L O2 / L blood = 100 mL O2 / dL blood

  # Therefore:
  #   CaO2 - CvO2 = 100 / Qm_VO2m_ratio
  
  extraction <- 100 / ratio
  O2_content <- CaO2 - extraction * capillary_distance
  S          <- O2_content / (1.34 * Hb)
  PO2        <- Severinghaus_PO2(S)
  
  curve_data <- rbind(
    curve_data, 
    data.frame(
      capillary_distance = capillary_distance,
      PO2 = PO2,
      Qm_VO2m_ratio = factor(ratio)
    )
  )
  
  # First and second derivatives
  dPO2        <- gradient(PO2,  capillary_distance)
  d2PO2       <- gradient(dPO2, capillary_distance)
  sign_change <- which(diff(sign(d2PO2)) != 0)
  
  if (length(sign_change) > 0) {
    i <- sign_change[1]
    
    inflections <- rbind(
      inflections,
      data.frame(
        Qm_VO2m_ratio = ratio,
        capillary_distance_inflection = capillary_distance[i],
        PO2_inflection = PO2[i]
      )
    )
  }
}

# BUILD PLOT 
# ------------------------------------------------------------------------------
p <- ggplot(curve_data, aes(x = capillary_distance, y = PO2, group = Qm_VO2m_ratio)) +
  geom_line(linewidth = 0.8) +
  geom_point(
    data = inflections,
    aes(x = capillary_distance_inflection, y = PO2_inflection),
    inherit.aes = FALSE,
    size = 2
  ) +
  geom_hline(yintercept = 21.5, linetype = "dashed") +
  annotate(
    "text",
    size = 3,
    x = 0.25,
    y = 15,
    label = "'PO'[2]~'at HbO'[2]~'inflection'",
    parse = TRUE,
    vjust = -0.4
  ) +
  annotate("text", x = 0.98, y = 3, label = "5", hjust = -0.3) +
  geom_text(
    data = subset(curve_data, capillary_distance == max(capillary_distance)),
    aes(label = as.character(Qm_VO2m_ratio)),
    hjust = -0.3
  ) +
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.2), 
    expand = expansion(mult = c(0.02, 0.12))
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 20),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "Fractional capillary distance",
    y = expression("Capillary " * PO[2] * " (mmHg)")
  ) +
  theme_classic()

print(p)

ggsave("capillary_PO2.pdf", p, width = 90, height = 100, units = "mm")

ggsave(
  filename    = "capillary_PO2.tiff",
  plot        = p,
  device      = "tiff",
  width       = 90,
  height      = 100,
  units       = "mm",
  dpi         = 600,
  compression = "lzw"
)

print(inflections)


# Session information for reproducibility
session_info <- c(capture.output(sessionInfo()),"", capture.output(rstudioapi::versionInfo()))
writeLines(session_info, "session_info.txt") 
