################################################################################
# Written by Loïc Labache, Ph.D.                                               #
# Holmes Lab, Department of Psychiatry - Rutgers University                    #
# May 12, 2025                                                                 #
################################################################################

################################################################################
# Script to test language lateralization effects;
# As long as the data file have the same columns names.
#
#       Lateralization_binary (factor TYP / ATYP)
#       Manual_Preference (factor G / D)
#       MainFUN (factor MG / MD)
#       Age_years
#       Educational_years
#       Total_Intracranial_Volume_eTIV_mm3
#       Gender (factor F / H)
################################################################################

options(contrasts = c("contr.sum", "contr.poly"))

# Libraries --------------------------------------------------------------------
library(psych)
library(car)
library(broom)
library(emmeans)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggdist)
library(prismatic)
library(effectsize)
library(purrr)

# Function ---------------------------------------------------------------------
run_language_effect_single_network = function(data,
                                              outcome_var,
                                              network_label = outcome_var) {
  data = data %>%
    mutate(Lateralization_binary = factor(Lateralization_binary,
                                          levels = c("TYP", "ATYP")),
           Manual_Preference = factor(Manual_Preference, levels = c("G", "D")),
           MainFUN = factor(MainFUN, levels = c("MG", "MD")),
           Gender = factor(Gender, levels = c("F", "H")))
  form = as.formula(paste0(outcome_var,
                           " ~ Lateralization_binary * Manual_Preference + ",
                           "MainFUN + Age_years + Educational_years + ",
                           "Total_Intracranial_Volume_eTIV_mm3 + Gender"))
  mod = lm(form, data = data)
  
  gl = glance(mod)
  p_model = gl$p.value
  R2_adj = gl$adj.r.squared
  beta_lat = coef(mod)["Lateralization_binary1"] %||% NA_real_
  t_lat = summary(mod)$coefficients["Lateralization_binary1", "t value"] %||% NA_real_
  
  anova_tbl = Anova(mod, type = "III")
  
  p_anova_lat = anova_tbl["Lateralization_binary", "Pr(>F)"]
  
  eta_tbl = eta_squared(anova_tbl, partial = TRUE)
  eta2_lat = eta_tbl[eta_tbl$Parameter == "Lateralization_binary", "Eta2_partial"]
  
  emm_lat = emmeans(mod, ~ Lateralization_binary)
  pair_lat = pairs(emm_lat, adjust = "tukey")
  
  emm_df = as.data.frame(emm_lat)
  diff_df = as.data.frame(summary(pair_lat, infer = TRUE))
  
  lsm_typ = emm_df$emmean[emm_df$Lateralization_binary == "TYP"]
  lsm_atyp = emm_df$emmean[emm_df$Lateralization_binary == "ATYP"]
  
  ci_typ = emm_df[emm_df$Lateralization_binary == "TYP",
                  c("lower.CL", "upper.CL")]
  ci_atyp = emm_df[emm_df$Lateralization_binary == "ATYP",
                   c("lower.CL", "upper.CL")]
  
  diff_est = diff_df$estimate[1]
  diff_ci = diff_df[1, c("lower.CL", "upper.CL")]
  diff_p = diff_df$p.value[1]
  
  summary_row = tibble(Network = network_label,
                       p_model = p_model,
                       R2_adj = R2_adj,
                       beta_lat = beta_lat,
                       t_value_lat = t_lat,
                       p_main_effect_lat = p_anova_lat,
                       partial_eta2_lat = eta2_lat,
                       lsm_typ = lsm_typ,
                       lsm_typ_lowerCI = ci_typ$lower.CL,
                       lsm_typ_upperCI = ci_typ$upper.CL,
                       lsm_atyp = lsm_atyp,
                       lsm_atyp_lowerCI = ci_atyp$lower.CL,
                       lsm_atyp_upperCI = ci_atyp$upper.CL,
                       diff_typ_minus_atyp = diff_est,
                       diff_lowerCI = diff_ci$lower.CL,
                       diff_upperCI = diff_ci$upper.CL,
                       diff_p_value = diff_p)
  
  pal_col = c("#ff8c00", "#9f34f0")
  
  plot_data = data %>%
    select(Lateralization_binary, !!sym(outcome_var)) %>%
    rename(Outcome = !!sym(outcome_var)) %>%
    filter(!is.na(Outcome))
  
  plot_data$Lateralization_binary = factor(plot_data$Lateralization_binary,
                                           levels = c("TYP", "ATYP"))
  
  lsm_plot = emm_df %>%
    transmute(Lateralization_binary,
              emmean,
              low_CI = lower.CL,
              upp_CI = upper.CL )
  
  fig_net = ggplot() +
    geom_hline(yintercept = 0,
               linetype = "dashed",
               color = "grey50",
               alpha = 0.25) +
    stat_slab(data = plot_data,
              aes(x = as.numeric(Lateralization_binary),
                  y = Outcome,
                  colour = Lateralization_binary,
                  fill = Lateralization_binary),
              side = "left",
              justification = 1,
              width = 0.4,
              expand = FALSE,
              trim = FALSE,
              normalize = "panels",
              alpha = 0.3) +
    geom_point(data = plot_data,
               aes(x = as.numeric(Lateralization_binary) + 0.35,
                   y = Outcome,
                   colour = Lateralization_binary,
                   fill = after_scale(clr_alpha(color, 0.4))),
               position = position_jitter(width = 0.10, height = 0),
               size = 0.5,
               stroke = 1,
               alpha = 0.3,
               shape = 21) +
    geom_errorbar(data = lsm_plot,
                  aes(x = as.numeric(Lateralization_binary) + 0.117,
                      ymin = low_CI,
                      ymax = upp_CI,
                      colour = Lateralization_binary),
                  width = 0.15,
                  linewidth = 0.5) +
    geom_point(data = lsm_plot,
               aes(x = as.numeric(Lateralization_binary) + 0.117,
                   y = emmean,
                   colour = Lateralization_binary,
                   fill = Lateralization_binary),
               shape = 21,
               size = 1) +
    scale_color_manual(values = pal_col) +
    scale_fill_manual(values  = pal_col) +
    scale_x_continuous(breaks = 1:2,
                       labels = levels(plot_data$Lateralization_binary)) +
    labs(x = "Language Lateralization",
         y = paste0("Asymmetry (", outcome_var, ")"),
         title = paste0("Asymmetry by Language Phenotype – ", network_label)) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none",
          strip.background = element_blank())
  
  out = list(network = network_label,
             model = mod,
             anova = anova_tbl,
             emm = emm_lat,
             pairwise = pair_lat,
             summary_row = summary_row,
             plot = fig_net)
  return(out)
}

#------------------------------------------------------------------------------
# Example: Task - Parieto-Frontal network (ALANs)
#------------------------------------------------------------------------------
my_project = "my_PATH"
dataNet = read.csv(file.path(my_project, "Data", "BIL&GIN_287participants_ALANsLUCA_TaskRest.txt"))
res_PF = run_language_effect_single_network(data = dataNet,
                                            outcome_var = "ALANs_LBJ_weighted_BOLD_ParietoFrontal_Asym",
                                            network_label = "ALANs_ParietoFrontal")
res_PF$summary_row
res_PF$plot
