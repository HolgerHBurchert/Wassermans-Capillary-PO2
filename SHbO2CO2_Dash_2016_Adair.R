#---------------------------------------------------------------------------------------
# Simulation of oxyhemoglobin (HbO2) and carbaminohemoglobin (HbCO2) dissociation
# curves and computation of total O2 and CO2 contents in whole blood, revised from
# the original model of Dash and Bassingthwaighte, ABME 38(4):1683-1701, 2010. The
# revision makes the model further simplified, as it bypasses the computations of
# the indices n1, n2, n3 and n4, which are complex expressions. Rather the revision
# necessiates the computations of K4p in terms of P50. Also the calculations of P50
# in terms of pH is enhanced based on a 3rd degree polynomial interpolation. 
#---------------------------------------------------------------------------------------
# Developed by: Ranjan Dash, PhD (Last modified: 2/16/2017; structured input and output)
# Department of Physiology and Biotechnology and Bioengineering Center
# Medical College of Wisconsin, Milwaukee, WI-53226
#---------------------------------------------------------------------------------------

# Input physiological variables for calculations of SHbO2CO2 and blood O2CO2 contents
SHbO2CO2 <- function(PO2 = 100, PCO2 = 40, pHrbc = 7.24, DPGrbc = 0.00465, Temp = 37, Hbrbc = 0.00528, Hct = 0.45) {
  
  # This model can provide unrealistic P50 values if pHrbc > 10 (highly unphysiological).
  if (pHrbc > 10) {
    #disp('nonrealistic pHrbc > 10; so set it as 10'); disp(pHrbc);
    pHrbc <- 10
  }
  
  # Parameters those are fixed in the model (i.e. water fractions, RBCs hemoglobin
  # concentration, equilibrium constants, and Hill coefficient)
  Wpl <- 0.94                           # fractional water space of plasma; unitless
  Wrbc <- 0.65                          # fractional water space of RBCs; unitless
  Wbl <- (1 - Hct) * Wpl + Hct * Wrbc   # fractional water space of blood; unitless
  K1 <- 10^(-6.12)                      # CO2 hydration reaction equilibrium constant K
  K2 <- 21.5e-6                         # CO2 + HbNH2 equilibrium constant; unitless
  K2dp <- 1e-6                          # HbNHCOOH dissociation constant; M
  K2p <- K2 / K2dp                      # kf2p/kb2p; 1/M
  K3 <- 11.3e-6                         # CO2 + O2HbNH2 equilibrium constant; unitless
  K3dp <- 1e-6                          # O2HbNHCOOH dissociation constant; M
  K3p <- K3 / K3dp                      # kf3p/kb3p; 1/M
  K5dp <- 2.4e-8                        # HbNH3+ dissociation constant; M
  K6dp <- 1.2e-8                        # O2HbNH3+ dissociation constant; M
  Rrbc <- 0.69                          # Gibbs-Donnan ratio across the RBC membrane
  mol2ml <- 22267.4                     # Conversion factor from mol of gas to ml of gas at STP
  
  # Variables those are fixed in the model with values at standard physiological conditions
  # (i.e. PO20, PCO20, pHpl0, pHrbc0, DPGrbc0, Temp0, P500,aO20, aCO20)
  PO20 <- 100                           # standard O2 partial pressure in blood; mmHg
  PCO20 <- 40                           # standard CO2 partial pressure in blood; mmHg
  pHrbc0 <- 7.24                        # standard pH in RBCs; unitless
  pHpl0 <- pHrbc0 - log10(Rrbc)         # standard pH in plsama; unitless
  DPGrbc0 <- 4.65e-3                    # standard 2,3-DPG concentration in RBCs; M
  Temp0 <- 37                           # standard temperature in blood; degC
  P500 <- 26.8                          # standard PO2 for 50% SHbO2; mmHg
  aO20 <- 1.46e-6                       # solubility of O2 in water at 37 C; M/mmHg
  aCO20 <- 32.66e-6                     # solubility of CO2 in water at 37 C; M/mmHg
  # aO20 = 1.1129e-6;                   # solubility of O2 in water at 37 C; M/mmHg
  # aCO20 = 2.6715e-5;                  # solubility of CO2 in water at 37 C; M/mmHg
  
  # Calculation of intermediate variables in the computations of SHbO2 and SHbCO2
  pHpl <- pHrbc - log10(Rrbc)           # Rrbc = Hpl/Hrbc = 10^-(pHpl-pHrbc)
  delpHrbc <- pHrbc - pHrbc0; delpHpl <- pHpl - pHpl0
  delPCO2 <- PCO2 - PCO20; delDPGrbc <- DPGrbc - DPGrbc0; delTemp <- Temp - Temp0
  aO2 <- aO20 * (1 - 1e-2 * delTemp + 4.234e-4 * delTemp^2) # Corrected solubility of O2
  aCO2 <- aCO20 * (1 - 1.86e-2 * delTemp + 6.515e-4 * delTemp^2) # Corrected solubility of CO2
  O2 <- aO2 * PO2
  CO2 <- aCO2 * PCO2
  Hrbc <- 10^(-pHrbc)
  Hpl <- 10^(-pHpl)
  
  P501 <- P500 + 1.2 * (-21.279 * delpHrbc + 8.872 * delpHrbc^2 - 1.47 * delpHrbc^3) # all standard conditions, except pH
  P502 <- P500 + 1.7 * (4.28e-2 * delPCO2 + 3.64e-5 * delPCO2^2) # all standard conditions, except CO2
  P503 <- P500 + 1.0 * (795.633533 * delDPGrbc - 19660.8947 * delDPGrbc^2) # all standard conditions, except DPG
  P504 <- P500 + 0.98 * (1.4945 * delTemp + 4.335e-2 * delTemp^2 + 7e-4 * delTemp^3) # all standard conditions, except T
  P50 <- P500 * (P501 / P500) * (P502 / P500) * (P503 / P500) * (P504 / P500)
  C50 <- aO2 * P50
  
  # NEW IN THIS VERSION OF THE CODE: PO2 dependent variable Hill coefficient
  # nH = 2.7; Hill coefficient; unitless (redefined as a function PO2)
  alpha <- 2.8    #  Roughton et al data
  beta  <- 1.2    #  Roughton et al data
  gamma <- 29.2   #  Roughton et al data
  nH <- alpha - beta * 10^(-PO2 / gamma)    #PO2 dependent variable nH
  
  # Compute the apparent equilibrium constant of Hb with O2 and CO2 (KHbO2 and KHbCO2); 
  # O2 and CO2 saturations of Hb (SHbO2 and SHbCO2); and O2 and CO2 contents in blood. 
  BPH1 <- 1 + K2dp / Hrbc # Binding polynomial involving K2dp and Hrbc 
  BPH2 <- 1 + K3dp / Hrbc # Binding polynomial involving K3dp and Hrbc 
  BPH3 <- 1 + Hrbc / K5dp # Binding polynomial involving K5dp and Hrbc 
  BPH4 <- 1 + Hrbc / K6dp # Binding polynomial involving K6dp and Hrbc 
  K4p <- (O2^(nH - 1) * (K2p * BPH1 * CO2 + BPH3)) / (C50^nH * (K3p * BPH2 * CO2 + BPH4))
  KHbO2 <- K4p * (K3p * BPH2 * CO2 + BPH4) / (K2p * BPH1 * CO2 + BPH3)
  KHbCO2 <- (K2p * BPH1 + K3p * K4p * BPH2 * O2) / (BPH3 + K4p * BPH4 * O2)
  SHbO2 <- KHbO2 * O2 / (1 + KHbO2 * O2)
  SHbCO2 <- KHbCO2 * CO2 / (1 + KHbCO2 * CO2)
  
  O2free <- Wbl * O2 #  M (mol O2 per L blood)
  O2bound <- 4 * Hct * Hbrbc * SHbO2 # M (mol O2 per L blood)
  O2tot <- O2free + O2bound # M (mol O2 per L blood)
  O2cont <- mol2ml * O2tot / 10 # mL O2/100 mL blood
  CO2free <- Wbl * CO2 # M (mol CO2 per L blood)
  CO2bicarb <- ((1 - Hct) * Wpl + Hct * Wrbc * Rrbc) * (K1 * CO2 / Hpl) # M (mol CO2 per L blood)
  CO2bound <- 4 * Hct * Hbrbc * SHbCO2 # M (mol CO2 per L blood)
  CO2totfb <- CO2free + CO2bound # M (mol CO2 per L blood)
  CO2contfb <- mol2ml * CO2totfb / 10 # mL CO2/100 mL blood
  CO2tot <- CO2free + CO2bicarb + CO2bound # M (mol CO2 per L blood)
  CO2cont <- mol2ml * CO2tot / 10 # mL CO2/100 mL blood
  
  # ALTERNATIVE CALCULATIONS FOR SHBO2 AND SHBCO2 BASED ON DIFFERENT HB-BOUND SPECIES
  HbNH2 <- Hbrbc / ((K2p * CO2 * BPH1 + BPH3) + K4p * O2 * (K3p * CO2 * BPH2 + BPH4))
  HbNH3p <- HbNH2 * Hrbc / K5dp
  O2HbNH2 <- K4p * O2 * HbNH2
  O2HbNH3p <- O2HbNH2 * Hrbc / K6dp
  HbNHCOOH <- K2p * CO2 * HbNH2
  HbNHCOOm <- K2dp * HbNHCOOH / Hrbc
  O2HbNHCOOH <- K3p * CO2 * O2HbNH2
  O2HbNHCOOm <- K3dp * O2HbNHCOOH / Hrbc
  SHbO2kin <- (O2HbNH2 + O2HbNH3p + O2HbNHCOOH + O2HbNHCOOm) / Hbrbc
  SHbCO2kin <- (HbNHCOOH + HbNHCOOm + O2HbNHCOOH + O2HbNHCOOm) / Hbrbc
  

  # ADAIR ADD-ON: PO2-shift Adair ligation states Hb4..Hb4(O2)4
  # ---------------------------------------------------------------------------
  # Fitted PO2-shift parameters from the Dash_Adair_Bridge.R file 
  k_pH  <-  0.41958876  
  k_CO2 <-  0.11986640 
  k_T   <- -0.02413729
  k_DPG <-  0.11146272
  
  # Adair anchor coefficients (Roughton 1972 set)
  a1_std_ad <- 2.18e-2
  a2_std_ad <- 9.12e-4
  a3_std_ad <- 3.75e-6
  a4_std_ad <- 2.47e-6
  
  # Anchor conditions (must match training anchor)
  pHrbc_std_ad  <- 7.24
  PCO2_std_ad   <- 40
  Temp_std_ad   <- 37
  DPGrbc_std_ad <- 0.00465
  
  # Deltas from anchor (MATCHING THE FIT)
  dpH_ad      <- pHrbc - pHrbc_std_ad
  dlogCO2_ad  <- log10(PCO2_std_ad / PCO2)         # std/PCO2
  dT_ad       <- Temp - Temp_std_ad
  dlogDPG_ad  <- log10(DPGrbc_std_ad / DPGrbc)     # std/DPGrbc
  
  # PO2 scaling (MATCHING THE FIT)
  PO2scale_adair <- 10^(k_pH*dpH_ad + k_CO2*dlogCO2_ad + k_T*dT_ad + k_DPG*dlogDPG_ad)
  PO2vir_adair   <- PO2 * PO2scale_adair
  
  # Adair polynomial at virtual PO2
  term1_ad <- a1_std_ad * PO2vir_adair
  term2_ad <- a2_std_ad * PO2vir_adair^2
  term3_ad <- a3_std_ad * PO2vir_adair^3
  term4_ad <- a4_std_ad * PO2vir_adair^4
  
  Z_ad <- 1 + term1_ad + term2_ad + term3_ad + term4_ad
  
  Hb4     <- 1 / Z_ad
  Hb4O2_1 <- term1_ad / Z_ad
  Hb4O2_2 <- term2_ad / Z_ad
  Hb4O2_3 <- term3_ad / Z_ad
  Hb4O2_4 <- term4_ad / Z_ad
  
  # Saturation implied by Adair ligation states (diagnostic)
  SHbO2_adair <- (Hb4O2_1 + 2*Hb4O2_2 + 3*Hb4O2_3 + 4*Hb4O2_4) / 4
  
  # Output physiological variables that are computed (e.g. SHbO2CO2 and blood O2CO2 contents)
  return(list(aO2 = aO2, 
              aCO2 = aCO2, 
              nH = nH, 
              P50 = P50, 
              K4p = K4p, 
              KHbO2 = KHbO2, 
              KHbCO2 = KHbCO2, 
              SHbO2 = SHbO2, 
              SHbCO2 = SHbCO2, 
              O2tot = O2tot, 
              O2cont = O2cont, 
              CO2tot = CO2tot, 
              CO2cont = CO2cont, 
              CO2totfb = CO2totfb,
              CO2contfb = CO2contfb, 
              HbNH2 = HbNH2, 
              HbNH3p = HbNH3p, 
              O2HbNH2 = O2HbNH2, 
              O2HbNH3p = O2HbNH3p, 
              HbNHCOOH = HbNHCOOH, 
              HbNHCOOm = HbNHCOOm, 
              O2HbNHCOOH = O2HbNHCOOH, 
              O2HbNHCOOm = O2HbNHCOOm, 
              SHbO2kin = SHbO2kin, 
              SHbCO2kin = SHbCO2kin, 
              
              # Added Adair/PO2-shift outputs (suffix _adair)
              PO2scale_adair = PO2scale_adair,
              PO2vir_adair = PO2vir_adair,
              SHbO2_adair = SHbO2_adair,
              Hb4 = Hb4,
              Hb4O2_1 = Hb4O2_1,
              Hb4O2_2 = Hb4O2_2,
              Hb4O2_3 = Hb4O2_3,
              Hb4O2_4 = Hb4O2_4))
}

#---------------------------------------------------------------------------------------
# BIOCHEMICAL REACTIONS FOR DERIVATION OF THE SHBO2 AND SHBCO2 EQUATIONS
# The equations for O2 and CO2 saturations of hemoglobin (SHbO2 and SHbCO2) are  
# derived by considering the various kinetic reactions involving the binding of
# O2 and CO2 with hemoglobin in RBCs:
#
#            kf1p       K1dp
# 1. CO2+H2O <--> H2CO3 <--> HCO3- + H+;  K1=(kf1p/kb1p)*K1dp
#            kb1p		K1 = 7.43e-7 M; K1dp = 5.5e-4 M
#
#              kf2p          K2dp
# 2. CO2+HbNH2 <--> HbNHCOOH <--> HbNHCOO- + H+;  K2=(kf2p/kb2p)*K2dp
#              kb2p		K2 = 21.5e-6; K2dp = 1.0e-6 M
#
#                kf3p            K3dp
# 3. CO2+O2HbNH2 <--> O2HbNHCOOH <--> O2HbNHCOO- + H+; K3=(kf3p/kb3p)*K3dp
#                kb3p		K3 = 11.3e-6; K3dp = 1.0e-6 M
#
#              kf4p          
# 4. O2+HbNH2 <--> O2HbNH2;  K4p=func([O2];[H+];[CO2];[DPG];T)
#              kb4p
#
#           K5dp
# 5. HbNH3+ <--> HbNH2 + H+; K5dp = 2.4e-8 M
#
#             K6dp
# 6. O2HbNH3+ <--> O2HbNH2 + H+; K6dp = 1.2e-8 M
#---------------------------------------------------------------------------------------
