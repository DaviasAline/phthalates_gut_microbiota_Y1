## Aline Davias
## 20/02/2024
## Analyses de la taxonomie 


# Packages loading ----
library(tidyverse)
library(phyloseq)
library(expss)
library(gtsummary)
library(Maaslin2)
library(labelled)
library(questionr)
library(sjlabelled)
library(openxlsx)
source("~/5. R projects/phthalates_gut_microbiota_Y1/3_programs/4_functions_AD_gumme.R")
rm(comp_effectifs, heatmap_cor_pairwise, model_covar, model_multi, model_summary, model_univ_multi, 
   table_cor, table_cor_sg, test_sensi_sg)
library(corrplot)
library(see)
library(psych)
library(compositions)
library(writexl)
library(broom)
library(qgcomp)

# Data reading ----
## Gut microbiota 
# taxa_table <- read_csv("0_source_data/taxa_table_ASVbased_Y1_AD_20220504_8.csv")
# bdd_microbiota <- 
#   read_labelled_csv(
#     "0_source_data/gut_microbiota_ASVbased_Y1_labelled_AD_20220504_7239.csv")%>%
#   select(
#     ident, 
#     starts_with("ch_feces_rel_")) %>%
#   filter(!is.na(ch_feces_rel_g1_Y1)) 
# 
# input_data <- 
#   read_labelled_csv(
#     "0_source_data/gut_microbiota_ASVbased_Y1_labelled_AD_20220504_7239.csv")%>%
#   select(
#     ident, 
#     starts_with("ch_feces_rel_g")) %>%
#   filter(!is.na(ch_feces_rel_g1_Y1)) %>%
#   column_to_rownames("ident")
# 
# var_label(input_data) <-  str_replace(
#   var_label(input_data),                                    # set correct variable names                   
#   "One year child feces relative abundance of ", "")
# colnames(input_data) <- var_label(input_data)
# 
# load("2_final_data/bdd_alpha.RData")
load("4_output/results_alpha_phyla.RData")
load("4_output/results_genera.RData")
rm(cor_mixed_table_expo, cor_mixed_table_outcome, df, forest_plot, 
   heatmap_genera, mahatan_plot, test_expo, test_outcome, test_outcome_2, 
   covar_vec_i, not_signi_vec, signi_vec, test, var, variables, pollutants_vec, 
   results_alpha_corrected_expo, results_alpha_corrected_outcome, 
   results_M0_corrected_expo, results_M0_corrected_outcome)

bdd_review <- left_join(bdd_alpha[,c("ident", alpha_vec, covariates, pollutants)],
                        bdd_taxa[,c("ident", 
                                    "ch_feces_rel_p1_Y1", "ch_feces_rel_p2_Y1",
                                    "ch_feces_rel_p3_Y1", "ch_feces_rel_p4_Y1")], 
                        by = "ident")
bdd_review <- bdd_review %>% mutate(ident = as.character(ident))
bdd_review <- left_join(bdd_review,
                        input_data_log, 
                        by = "ident")
metadata <- metadata %>% mutate(ident = as.character(ident))
bdd_review <- left_join(bdd_review, 
                        metadata[, c("ident", "ch_hospit_2M", "ch_hospit_2_12M")], 
                        by = "ident")

bdd_review <- bdd_review %>% 
  filter(!is.na(ch_feces_RUN_Y1)) %>% 
  mutate(
    ch_hospit_Y1 = ifelse(ch_hospit_2M | ch_hospit_2_12M > 0, 1, 0))

bdd_review <- bdd_review %>% filter(!is.na(ch_feces_RUN_Y1))

# Vectors ----
# 
# genera_linear_complet <- bdd_microbiota %>%                                     # création du vecteur genera long                         
#   select(contains("ch_feces_rel_g")) %>%
#   filter(!is.na(ch_feces_rel_g1_Y1)) %>%
#   select_if(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.3) %>%
#   colnames()
# 
# genera_linear <- input_data %>%                                                 # création du vecteur genera court   
#   select_if(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.3) %>%
#   colnames()
# 
# input_data_log <- input_data %>%                                                # imputation et log transformation
#   select(all_of(genera_linear)) %>%
#   mutate_all(., ~ ifelse(. == 0, 1/5000, .)) %>%        # remplacement des valeurs 0 par 1/5000
#   mutate_all(~ log(.)) %>%                              # transformation logarithmique
#   rename_with(~gsub("genus ", "", .), everything()) %>% # changement des noms de colonnes pour qu'ils n'aient pas d'espace
#   rownames_to_column(var = "ident")
# 
# genera_linear <- str_replace_all(genera_linear, "genus ", "")                   # création du vecteur genera court   


covariates_sensi_14 <- tibble(value = covariates) %>%                 # création d'un vecteur sans les covariables suspectées d'overadjustement
  filter(!value %in% c("po_w_kg_3cat", "po_he_3cat_i", "po_gd")) %>%
  pull(value)

covariates_sensi_15 <- paste(covariates, "ch_hospit_Y1")

phyla_vec <- c("ch_feces_rel_p1_Y1", "ch_feces_rel_p2_Y1", "ch_feces_rel_p3_Y1", "ch_feces_rel_p4_Y1")

# New functions ----
model <- function(data,                 # pas de changement 
                  outcome, 
                  exposure_vec, 
                  digit_beta_IC) {
  data %>%                  
    select(
      {{outcome}},
      all_of({{exposure_vec}}),
      all_of(covariates)) %>% 
    
    tbl_uvregression(
      method = lm ,
      y = {{outcome}},
      formula = "{y} ~ 
         {x} +
         ch_feces_RUN_Y1  +
         ch_feces_age_w_Y1_i +
         po_delmod +
         ch_food_intro_Y1_3cat_i +
         ch_antibio_Y1_2cat_i +
         mo_par_2cat +
         mo_pets_i +
         ch_sex +
         mo_tob_gr_anyt_yn_n2_i +
         Mo_ETS_anyT_yn1_opt_i +
         ch_ETS_12m_opt36m +
         mo_interpreg_3cat +
         mo_dipl_3cat_i +
         po_w_kg_3cat +
         po_he_3cat_i +
         ch_w_Y1_3cat_i +
         ch_he_Y1_3cat_i +
         po_gd +
         mo_age +
         mo_bmi_bepr_3cat_i +
         bf_duration_till48w_4cat_i",
      hide_n = TRUE,
      pvalue_fun = ~ style_pvalue(.x, digits = 2),
      estimate_fun = ~ style_sigfig(.x, digits = {{digit_beta_IC}})
    ) %>%
    add_global_p(keep = TRUE, singular.ok = TRUE) %>%
    bold_labels() }


model_sensi_14 <- function(data, 
                        outcome, 
                        exposure_vec, 
                        digit_beta_IC) {
  data %>%                  
    select(
      {{outcome}},
      all_of({{exposure_vec}}),
      all_of(covariates_sensi_14)) %>% 
    
    tbl_uvregression(
      method = lm ,
      y = {{outcome}},
      formula = "{y} ~ 
         {x} +
         ch_feces_RUN_Y1  +
         ch_feces_age_w_Y1_i +
         po_delmod +
         ch_food_intro_Y1_3cat_i +
         ch_antibio_Y1_2cat_i +
         mo_par_2cat +
         mo_pets_i +
         ch_sex +
         mo_tob_gr_anyt_yn_n2_i +
         Mo_ETS_anyT_yn1_opt_i +
         ch_ETS_12m_opt36m +
         mo_interpreg_3cat +
         mo_dipl_3cat_i +
         ch_w_Y1_3cat_i +
         ch_he_Y1_3cat_i +
         mo_age +
         mo_bmi_bepr_3cat_i +
         bf_duration_till48w_4cat_i",
      hide_n = TRUE,
      pvalue_fun = ~ style_pvalue(.x, digits = 2),
      estimate_fun = ~ style_sigfig(.x, digits = {{digit_beta_IC}})
    ) %>%
    add_global_p(keep = TRUE, singular.ok = TRUE) %>%
    bold_labels() }



model_sensi_15 <- function(data, 
                           outcome, 
                           exposure_vec, 
                           digit_beta_IC) {
  data %>%                  
    select(
      {{outcome}},
      all_of({{exposure_vec}}),
      all_of(covariates_sensi_15)) %>% 
    
    tbl_uvregression(
      method = lm ,
      y = {{outcome}},
      formula = "{y} ~ 
         {x} +
          ch_feces_RUN_Y1  +
         ch_feces_age_w_Y1_i +
         po_delmod +
         ch_food_intro_Y1_3cat_i +
         ch_antibio_Y1_2cat_i +
         mo_par_2cat +
         mo_pets_i +
         ch_sex +
         mo_tob_gr_anyt_yn_n2_i +
         Mo_ETS_anyT_yn1_opt_i +
         ch_ETS_12m_opt36m +
         mo_interpreg_3cat +
         mo_dipl_3cat_i +
         po_w_kg_3cat +
         po_he_3cat_i +
         ch_w_Y1_3cat_i +
         ch_he_Y1_3cat_i +
         po_gd +
         mo_age +
         mo_bmi_bepr_3cat_i +
         bf_duration_till48w_4cat_i +
         ch_hospit_Y1",
      hide_n = TRUE,
      pvalue_fun = ~ style_pvalue(.x, digits = 2),
      estimate_fun = ~ style_sigfig(.x, digits = {{digit_beta_IC}})
    ) %>%
    add_global_p(keep = TRUE, singular.ok = TRUE) %>%
    bold_labels() }

extract_qgcomp_info <- function(qgcomp_result) {                                # Création fonction pour extraire les informations d'intérêt 
  data.frame(
    Estimate = qgcomp_result$psi,
    Lower_CI = qgcomp_result$ci[1],
    Upper_CI = qgcomp_result$ci[2],
    p_value = qgcomp_result$pval[2]
  )
}



lm_func_sensi_15 <- function(outcome, exposure, data){   
  model <- lm({{outcome}} ~
                exposure +
                ch_feces_RUN_Y1 +
                ch_feces_age_w_Y1_i +
                po_delmod +
                ch_food_intro_Y1_3cat_i +
                ch_antibio_Y1_2cat_i +
                mo_par_2cat + 
                mo_pets_i +
                ch_sex +
                mo_tob_gr_anyt_yn_n2_i +
                Mo_ETS_anyT_yn1_opt_i +
                ch_ETS_12m_opt36m +
                mo_interpreg_3cat +
                mo_dipl_3cat_i +
                po_w_kg_3cat +
                po_he_3cat_i +
                ch_w_Y1_3cat_i +
                ch_he_Y1_3cat_i +
                po_gd +
                mo_age +
                mo_bmi_bepr_3cat_i +
                bf_duration_till48w_4cat_i +
                ch_hospit_Y1,
              data = data)
  return(model)
}

## Metadata
# load("2_final_data/metadata.RData")
# input_metadata <- metadata %>%
#   select(ident,                      
#          all_of(covariates), 
#          all_of(pollutants)) %>%
#   filter(!is.na(ch_feces_RUN_Y1)) %>%
#   column_to_rownames("ident")
# rm(metadata)
# 
# input_metadata <- input_metadata %>%
#   mutate(across(c("ch_feces_RUN_Y1",                   
#                    "po_delmod",
#                    "ch_food_intro_Y1_3cat_i",
#                    "ch_antibio_Y1_2cat_i",
#                    "mo_par_2cat",
#                    "mo_pets_i",
#                    "ch_sex",
#                    "mo_tob_gr_anyt_yn_n2_i",
#                    "Mo_ETS_anyT_yn1_opt_i",
#                    "ch_ETS_12m_opt36m",
#                    "mo_interpreg_3cat",
#                    "mo_dipl_3cat_i",
#                    "po_w_kg_3cat",
#                    "po_he_3cat_i",
#                    "ch_w_Y1_3cat_i",
#                    "ch_he_Y1_3cat_i",
#                    "mo_bmi_bepr_3cat_i",
#                    "bf_duration_till48w_4cat_i"), as.factor))



# Table A14: overadjustment? ----
## Table A14-a (alpha diversity) ----
table_A14_a <- tbl_merge(                         
  tbls = list(
    effectif_column(data = bdd_alpha, 
                    outcome = ch_feces_SpecRich_5000_ASV_Y1, 
                    exposure_vec = pollutants),
    model(data = bdd_alpha, 
          outcome = ch_feces_SpecRich_5000_ASV_Y1,
          exposure_vec = pollutants, 
          digit_beta_IC = 1),
    model_sensi_14(data = bdd_alpha, 
                outcome = ch_feces_SpecRich_5000_ASV_Y1,
                exposure_vec = pollutants, 
                digit_beta_IC = 1),
    model(data = bdd_alpha,
          outcome = ch_feces_Shannon_5000_ASV_Y1,
          exposure_vec = pollutants, 
          digit_beta_IC = 2), 
    model_sensi_14(data = bdd_alpha,
                outcome = ch_feces_Shannon_5000_ASV_Y1,
                exposure_vec = pollutants, 
                digit_beta_IC = 2)), 
  tab_spanner = c("", 
                  "**Specific richness, main analysis**", 
                  "**Specific richness, sensitivity analysis**", 
                  "**Shannon diversity, main analysis**", 
                  "**Shannon diversity, sensitivity analysis**"))

## Table A14-b (phyla) ----
table_A14_b <- tbl_merge(
  tbls = list(
    effectif_column(
      data = bdd_taxa, 
      outcome = ch_feces_rel_p1_Y1, 
      exposure_vec = pollutants),
    model(
      data = bdd_taxa, 
      outcome = ch_feces_rel_p1_Y1,
      exposure_vec = pollutants, 
      digit_beta_IC = 1),
    model_sensi_14(
      data = bdd_taxa, 
      outcome = ch_feces_rel_p1_Y1,
      exposure_vec = pollutants, 
      digit_beta_IC = 1),
    model(
      data = bdd_taxa,
      outcome = ch_feces_rel_p2_Y1,
      exposure_vec = pollutants, 
      digit_beta_IC = 1),
    model_sensi_14(
      data = bdd_taxa,
      outcome = ch_feces_rel_p2_Y1,
      exposure_vec = pollutants, 
      digit_beta_IC = 1),
    model(
      data = bdd_taxa, 
      outcome = ch_feces_rel_p3_Y1,
      exposure_vec = pollutants, 
      digit_beta_IC = 1), 
    model_sensi_14(
      data = bdd_taxa, 
      outcome = ch_feces_rel_p3_Y1,
      exposure_vec = pollutants, 
      digit_beta_IC = 1), 
    model(
      data = bdd_taxa, 
      outcome = ch_feces_rel_p4_Y1,
      exposure_vec = pollutants, 
      digit_beta_IC = 1), 
    model_sensi_14(
      data = bdd_taxa, 
      outcome = ch_feces_rel_p4_Y1,
      exposure_vec = pollutants, 
      digit_beta_IC = 1)), 
  tab_spanner = c("", 
                  "**Phylum Firmicutes, main analysis**",
                  "**Phylum Firmicutes, sensitivity analysis**",
                  "**Phylum Actinobacteria, main analysis**",
                  "**Phylum Actinobacteria, sensitivity analysis**",
                  "**Phylum Bacteroidetes, main analysis**",
                  "**Phylum Bacteroidetes, sensitivity analysis**",
                  "**Phylum Proteobacteria, main analysis**",
                  "**Phylum Proteobacteria, sensitivity analysis**"))


## Table A14-c (genera) ----
# input_metadata <- input_metadata %>% rownames_to_column(var = "ident")          # Rassemblement des données en 1 seul dataframe 
# data_df_log <- left_join(input_data_log, input_metadata, by = "ident")
# corres <- 
#   taxa_table %>% 
#   select(Phyla_corres = ch_feces_phylum_ASVbased_Y1, 
#          Class_corres = ch_feces_class_ASVbased_Y1, 
#          Order_corres = ch_feces_order_ASVbased_Y1, 
#          Family_corres = ch_feces_family_ASVbased_Y1, 
#          Outcome = ch_feces_genus_ASVbased_Y1) %>%
#   filter(Outcome %in% genera_linear) %>%
#   distinct(Outcome, .keep_all = TRUE)
# 
# tbls_by_outcome_log <- vector("list", length(genera_linear))                    # création liste pour stocker les tableaux par outcome
# names(tbls_by_outcome_log) <- genera_linear
# 
# for (outcome in genera_linear) {
#   tbls_for_outcome_log <- vector("list", length(pollutants))
#   names(tbls_for_outcome_log) <- pollutants
#   
#   for (exposure in pollutants) {                                                # formula setting
#     terms <- c(exposure, covariates)
#     formula <- reformulate(terms, response = outcome)
#     model <- lm(formula, data = data_df_log)
#     
#     tbl <-                                                                      # running linear regression
#       tbl_regression(
#         model, 
#         include = exposure,
#         estimate_fun = scales::label_number(accuracy = .01, decimal.mark = "."),
#         pvalue_fun = custom_pvalue_fun,
#         exponentiate = FALSE) %>%
#       bold_p() %>%
#       bold_labels() %>%
#       add_global_p(include = exposure, singular.ok = TRUE, keep = TRUE)
#     
#     tbls_for_outcome_log[[exposure]] <- tbl
#   }
#   tbls_by_outcome_log[[outcome]] <- tbls_for_outcome_log
# }


tbls_by_outcome_log_sensi <- vector("list", length(genera_linear))              # création liste pour stocker les tableaux par outcome
names(tbls_by_outcome_log_sensi) <- genera_linear

for (outcome in genera_linear) {
  tbls_for_outcome_log_sensi <- vector("list", length(pollutants))
  names(tbls_for_outcome_log_sensi) <- pollutants
  
  for (exposure in pollutants) {                                                # formula setting
    terms <- c(exposure, covariates_sensi)
    formula <- reformulate(terms, response = outcome)
    model <- lm(formula, data = data_df_log)
    
    tbl <-                                                                      # running linear regression
      tbl_regression(
        model, 
        include = exposure,
        estimate_fun = scales::label_number(accuracy = .01, decimal.mark = "."),
        pvalue_fun = custom_pvalue_fun,
        exponentiate = FALSE) %>%
      bold_p() %>%
      bold_labels() %>%
      add_global_p(include = exposure, singular.ok = TRUE, keep = TRUE)
    
    tbls_for_outcome_log_sensi[[exposure]] <- tbl
  }
  tbls_by_outcome_log_sensi[[outcome]] <- tbls_for_outcome_log_sensi
}
rm(terms, formula, model, tbls_for_outcome_log_sensi, outcome, exposure, tbl)

prep_table_A14_c <- vector("list", length = 24)                                   # Création liste vide pour stocker les résultats finaux

for (i in 1:24) {                                                               # Boucle sur chaque expo (1 à 24)
  
  tbls_to_merge <- vector("list", length = 92)                                  # Création liste temporaire pour stocker les tbl_regressions à fusionner pour l'index i
  spanner_labels <- vector("character", length = 92)
  
  for (j in 1:46) {                                                             # Boucle sur chaque outcome (1 à 46)
    
    tbls_to_merge[[2*j-1]] <- tbls_by_outcome_log[[j]][[i]]                     # Ajouter le tbl_regression de tbls_by_outcome_log et tbls_by_outcome_log_sensi à la liste temporaire
    tbls_to_merge[[2*j]]   <- tbls_by_outcome_log_sensi[[j]][[i]]
    
    spanner_labels[2*j-1] <- paste(genera_linear[j], "main analysis")           # Attribuer les entetes correspondants en utilisant le vecteur genera_linear
    spanner_labels[2*j]   <- paste(genera_linear[j], "sensitivity analysis")
  }
  
  prep_table_A14_c[[i]] <- tbl_merge(tbls_to_merge, tab_spanner = spanner_labels) # Faire un tbl_merge pour le i-ème groupe de tbl_regressions et le stocker dans la liste finale
}
table_A14_c <- tbl_stack(prep_table_A14_c)                                          # Assemblage final des résultats


# ### Tableau pour générer des figures
# table_log <- tibble()                                                           # Initialisation d'un tibble vide pour stocker les résultats finaux
# for (i in seq_along(tbls_by_outcome_log_sensi)) {                               # Nom de l'outcome pour cette itération
#   outcome_name <- names(tbls_by_outcome_log_sensi)[i]
#   
#   for (j in seq_along(tbls_by_outcome_log_sensi[[i]])) {                        # Itération sur chaque tbl_regression dans la liste courante
#     exposure_name <- names(tbls_by_outcome_log_sensi[[i]])[j]                   # Nom de la variable d'exposition pour cette itération
#     tbl_data <- tbls_by_outcome_log_sensi[[i]][[j]] %>%                         # Extraction des données du tableau tbl_regression
#       as_tibble() %>%
#       mutate(Outcome = outcome_name, Exposure = exposure_name)
#     
#     table_log <- bind_rows(table_log, tbl_data)                                 # Ajout des données extraites au tibble final
#   }
# }
# # rm(terms, formula, model, 
# #    tbl, tbl_data, tbls_for_this_outcome_log, 
# #    stacked_tbl_log, tbls_for_outcome_log,
# #    stacked_tbls_by_outcome_log,
# #    exposure, exposure_name, i, j, outcome, outcome_name)
# 
# table_log_sensi <-                                                                    # Ajout variable de la correspondance en phyla
#   left_join(table_log_sensi, corres, by = "Outcome") %>% 
#   select(Phyla_corres, Class_corres, Order_corres, Family_corres, everything())
# 
# table_log_sensi <- table_log_sensi %>%
#   select(Phyla_corres, 
#          Outcome, 
#          Pollutants = Exposure, 
#          Beta = "**Beta**", 
#          "95% CI" = "**95% CI**",
#          "p-value" = "**p-value**", 
#          "Characteristic" = "**Characteristic**") %>%
#   mutate(
#     Phyla_corres = as.factor(Phyla_corres), 
#     Phyla_corres = fct_relevel(Phyla_corres,
#                                "Firmicutes", "Actinobacteria", 
#                                "Bacteroidetes", "Proteobacteria", 
#                                "Verrucomicrobia", "Candidatus_Saccharibacteria"),
#     Time_window = case_when(grepl("t2", Pollutants) ~ "Mother, pregnancy trim. 2", 
#                             grepl("t3", Pollutants) ~ "Mother, pregnancy trim. 3",
#                             grepl("Y1", Pollutants) ~ "Child, 12 months", 
#                             .default = "Mother, pregnancy trim. 2"),
#     Pollutants = str_replace_all(Pollutants,
#                                  c(
#                                    "mo_" = "",
#                                    "ch_" = "",
#                                    "DEHP" = "ΣDEHP",
#                                    "DiNP" = "ΣDiNP",
#                                    "DINCH" = "ΣDINCH",
#                                    "_ms_i_cor_t2_ln" = "", 
#                                    "_ms_i_cor_t3_ln" = "", 
#                                    "_ms_i_cor_Y1_ln" = "", 
#                                    "_i_cor_t2_ln" = "", 
#                                    "_i_cor_t3_ln" = "", 
#                                    "_i_cor_Y1_ln" = ""
#                                    
#                                  )),
#     Pollutants_Time_window = case_when(Time_window == "Mother, pregnancy trim. 2" ~ paste(Pollutants, "trim.2", sep = " "), 
#                                        Time_window == "Mother, pregnancy trim. 3" ~ paste(Pollutants, "trim.3", sep = " "), 
#                                        Time_window == "Child, 12 months" ~ paste(Pollutants, "12 months", sep = " ")), 
#     Pollutants_Time_window = 
#       fct_relevel(Pollutants_Time_window, 
#                   "ΣDINCH 12 months", "ΣDINCH trim.3", "ΣDINCH trim.2", "ohMPHP 12 months",
#                   "ohMPHP trim.3", "ohMPHP trim.2", "MEP 12 months", "MEP trim.3",
#                   "MEP trim.2", "MBzP 12 months", "MBzP trim.3", "MBzP trim.2",
#                   "MiBP 12 months", "MiBP trim.3", "MiBP trim.2", "ΣDiNP 12 months",
#                   "ΣDiNP trim.3", "ΣDiNP trim.2", "MnBP 12 months", "MnBP trim.3",
#                   "MnBP trim.2", "ΣDEHP 12 months", "ΣDEHP trim.3", "ΣDEHP trim.2"),
#     Pollutants_Time_window_rec = str_replace_all(Pollutants_Time_window, 
#                                                  c("trim.2" = "t2", 
#                                                    "trim.3" = "t3", 
#                                                    "12 months" = "Y1")), 
#     `p-value` = gsub("__", "", `p-value`),
#     `p-value` = as.numeric(`p-value`),
#     `q-value` = `p-value`/(29*31), 
#     p_value_shape = ifelse(`p-value`<0.05, "p-value<0.05", "p-value≥0.05"),
#     q_value_shape = ifelse(`q-value`<0.05, "q-value<0.05", "q-value≥0.05"), 
#     sens_beta = ifelse(Beta < 0, "Beta<0", "Beta≥0"), 
#     sens_beta = fct_relevel(sens_beta, "Beta≥0", "Beta<0"))  %>% 
#   separate(col = "95% CI", into = c("lower_CI", "upper_CI"), sep = ",", remove = FALSE) %>%
#   mutate(
#     lower_CI = as.numeric(lower_CI),
#     upper_CI = as.numeric(upper_CI)
#   ) %>%
#   select(
#     Phyla_corres,
#     Outcome, 
#     Pollutants, 
#     Time_window, 
#     Pollutants_Time_window, Pollutants_Time_window_rec, 
#     Beta, sens_beta, 
#     "95% CI", lower_CI, upper_CI, 
#     "p-value", p_value_shape, 
#     "q-value", q_value_shape)
# 
# table_log_sensi$Outcome_rec <- table_log_sensi$Outcome %>%
#   fct_recode(
#     "Clostridium IV" = "Clostridium_IV",
#     "Clostridium sensu stricto" = "Clostridium_sensu_stricto",
#     "Clostridium XlVa" = "Clostridium_XlVa",
#     "Clostridium XVIII" = "Clostridium_XVIII",
#     "Erysipelotrichaceae incertae sedis" = "Erysipelotrichaceae_incertae_sedis",
#     "Escherichia and Shigella" = "Escherichia_Shigella",
#     "Lachnospiracea incertae sedis" = "Lachnospiracea_incertae_sedis",
#     "Ruminococcus 2" = "Ruminococcus2",
#     "Saccharibacteria genera incertae sedis" = "Saccharibacteria_genera_incertae_sedis"
#   )

# Table A15: need adjustment on hospitalization? ----
prep_table_A15 <- vector("list", length(outcomes_vec))              # création liste pour stocker les tableaux par outcome
names(prep_table_A15) <- outcomes_vec

for (outcome in outcomes_vec) {
  tbls_for_outcome_log_sensi <- vector("list", length(pollutants))
  names(tbls_for_outcome_log_sensi) <- pollutants
  
  for (exposure in pollutants) {                                                # formula setting
    terms <- c(exposure, covariates_sensi_15)
    formula <- reformulate(terms, response = outcome)
    model <- lm(formula, data = bdd_review)
    
    tbl <-                                                                      # running linear regression
      tbl_regression(
        model, 
        include = exposure,
        estimate_fun = scales::label_number(accuracy = .01, decimal.mark = "."),
        pvalue_fun = custom_pvalue_fun,
        exponentiate = FALSE) %>%
      bold_p() %>%
      bold_labels() %>%
      add_global_p(include = exposure, singular.ok = TRUE, keep = TRUE)
    
    tbls_for_outcome_log_sensi[[exposure]] <- tbl
  }
  prep_table_A15[[outcome]] <- tbls_for_outcome_log_sensi
}
rm(terms, formula, model, tbls_for_outcome_log_sensi, outcome, exposure, tbl)

prep_table_A14_c <- vector("list", length = 24)                                   # Création liste vide pour stocker les résultats finaux

for (i in 1:24) {                                                               # Boucle sur chaque expo (1 à 24)
  
  tbls_to_merge <- vector("list", length = 92)                                  # Création liste temporaire pour stocker les tbl_regressions à fusionner pour l'index i
  spanner_labels <- vector("character", length = 92)
  
  for (j in 1:46) {                                                             # Boucle sur chaque outcome (1 à 46)
    
    tbls_to_merge[[2*j-1]] <- tbls_by_outcome_log[[j]][[i]]                     # Ajouter le tbl_regression de tbls_by_outcome_log et prep_table_A15 à la liste temporaire
    tbls_to_merge[[2*j]]   <- prep_table_A15[[j]][[i]]
    
    spanner_labels[2*j-1] <- paste(outcomes_vec[j], "main analysis")           # Attribuer les entetes correspondants en utilisant le vecteur outcomes_vec
    spanner_labels[2*j]   <- paste(outcomes_vec[j], "sensitivity analysis")
  }
  
  prep_table_A14_c[[i]] <- tbl_merge(tbls_to_merge, tab_spanner = spanner_labels) # Faire un tbl_merge pour le i-ème groupe de tbl_regressions et le stocker dans la liste finale
}
table_A14_c <- tbl_stack(prep_table_A14_c)   

# Table A16: qgcomp? ----
## Table A16-a (t2) ----
outcomes_vec <- c(alpha_vec, phyla_vec, genera_linear)
prep_table_A16_t2 <- vector("list",                                             # Création d'une liste pour stocker les résultats
                            length = 
                              length(outcomes_vec))             
names(prep_table_A16_t2) <- outcomes_vec
pollutants_t2 <- tibble(value = pollutants) %>%
  filter(str_detect(value, "t2")) %>%
  pull(value)

for (i in 1:length(outcomes_vec)) {                                            # Boucle sur chaque outcome pour effectuer l'analyse qgcomp
  outcome_name <- outcomes_vec[i]
  
  formula <- as.formula(paste(outcome_name, "~",                                # Formule pour le modèle avec covariables
                              paste(c(pollutants_t2, covariates), collapse = "+")))
  
  qgcomp_fit <- qgcomp.noboot(formula,                                          # réalisation de l'analyse qgcomp
                              data = bdd_review, 
                              expnms = pollutants_t2, 
                              q = 4)                                            # choix de 4 quantiles
  
  prep_table_A16_t2[[i]] <- qgcomp_fit                                          # Stockage le résultat dans la liste
}
rm(outcome_name, formula, qgcomp_fit, pollutants_t2, i)

## Table A16-b (t3) ----
prep_table_A16_t3 <- vector("list", length = length(outcomes_vec))             # Création d'une liste pour stocker les résultats
names(prep_table_A16_t3) <- outcomes_vec
pollutants_t3 <- tibble(value = pollutants) %>%
  filter(str_detect(value, "t3")) %>%
  pull(value)

for (i in 1:length(outcomes_vec)) {                                            # Boucle sur chaque outcome pour effectuer l'analyse qgcomp
  outcome_name <- outcomes_vec[i]
  
  formula <- as.formula(paste(outcome_name, "~",                                # Formule pour le modèle avec covariables
                              paste(c(pollutants_t3, covariates), collapse = "+")))
  
  qgcomp_fit <- qgcomp.noboot(formula,                                          # réalisation de l'analyse qgcomp
                              data = bdd_review, 
                              expnms = pollutants_t3, 
                              q = 4)                                            # choix de 4 quantiles
  
  prep_table_A16_t3[[i]] <- qgcomp_fit                                          # Stockage le résultat dans la liste
}
rm(outcome_name, formula, qgcomp_fit, pollutants_t3, i)

## Table A16-c (Y1) ----
prep_table_A16_Y1 <- vector("list", length = length(outcomes_vec))             # Création d'une liste pour stocker les résultats
names(prep_table_A16_Y1) <- outcomes_vec
pollutants_Y1 <- tibble(value = pollutants) %>%
  filter(str_detect(value, "Y1")) %>%
  pull(value)

for (i in 1:length(outcomes_vec)) {                                            # Boucle sur chaque outcome pour effectuer l'analyse qgcomp
  outcome_name <- outcomes_vec[i]
  
  formula <- as.formula(paste(outcome_name, "~",                                # Formule pour le modèle avec covariables
                              paste(c(pollutants_Y1, covariates), collapse = "+")))
  
  qgcomp_fit <- qgcomp.noboot(formula,                                          # réalisation de l'analyse qgcomp
                              data = bdd_review, 
                              expnms = pollutants_Y1, 
                              q = 4)                                            # choix de 4 quantiles
  
  prep_table_A16_Y1[[i]] <- qgcomp_fit                                          # Stockage le résultat dans la liste
}

prep_table_A16 <- list(prep_table_A16_t2, prep_table_A16_t3, prep_table_A16_Y1) # regroupement des résultats
names(prep_table_A16) <- c("t2", "t3", "Y1")
rm(outcome_name, formula, qgcomp_fit, pollutants_Y1, 
   prep_table_A16_t2, prep_table_A16_t3, prep_table_A16_Y1)



combined_results <- list()                                                      # Initialiser une liste pour stocker les résultats combinés

for (i in 1:length(outcomes_vec)) {                                            # Boucle pour extraire et combiner les résultats de t2, t3 et Y1
  outcome_name <- outcomes_vec[i]
  
  result_t2 <- extract_qgcomp_info(prep_table_A16$t2[[i]])                      # Extraire les informations pour chaque temps (t2, t3, Y1)
  result_t3 <- extract_qgcomp_info(prep_table_A16$t3[[i]])
  result_Y1 <- extract_qgcomp_info(prep_table_A16$Y1[[i]])
  
  result_t2$Time <- "t2"                                                        # Ajouter une colonne pour identifier le temps
  result_t3$Time <- "t3"
  result_Y1$Time <- "Y1"
  
  combined_result <- rbind(result_t2, result_t3, result_Y1)                     # Combiner les résultats pour cet outcome
  combined_result$Outcome <- outcome_name
  
  combined_results[[i]] <- combined_result                                      # Ajouter le résultat combiné à la liste des résultats finaux
}

table_A16 <- do.call(rbind, combined_results)                                   # Combiner tous les résultats en un seul data frame

table_A16 <- table_A16[, c("Outcome", "Time",                                   # Réorganiser les colonnes pour avoir Outcome en première colonne
                           "Estimate", "Lower_CI", 
                           "Upper_CI", "p_value")]

rm(combined_result, combined_results, result_t2, result_t3, result_Y1)

table_A16 <- table_A16 %>%
  mutate(
    Lower_CI = round(Lower_CI, 2) %>% as.character() %>% paste0(","),
    Upper_CI = round(Upper_CI, 2) %>% as.character(),
    "95% CI" = paste0(Lower_CI, Upper_CI)) %>%
  select(-Lower_CI, -Upper_CI) %>%
  select(Outcome, Time, Estimate, `95% CI`, p_value)

write.xlsx(table_A16, "4_output/Table_qgcomp.xlsx")



# stacked_tbls_by_outcome_log <- vector("list", length(genera_linear))            # Création liste pour stocker les tableaux empilés par outcome
# names(stacked_tbls_by_outcome_log) <- genera_linear
# 
# for (outcome in names(tbls_by_outcome_log)) {                                   # Récupérer les tableaux de régression pour cet outcome
#   tbls_for_this_outcome_log <- tbls_by_outcome_log[[outcome]]
#   stacked_tbl_log <- do.call(tbl_stack, list(tbls = tbls_for_this_outcome_log)) # Empiler les tableaux en un seul tableau
#   stacked_tbls_by_outcome_log[[outcome]] <- stacked_tbl_log                     # Ajouter le tableau empilé à la liste des tableaux empilés
# }
# 
# results_tbl_log <- tbl_merge(tbls = stacked_tbls_by_outcome_log,                # Fusionner les tableaux empilés en un seul tableau
#                              tab_spanner = genera_linear)




# Figures ----
### Mahatan plot Fig.4 ----
# mahatan_plot <- table_log  %>%
#   mutate(
#     Outcome_rec = 
#       fct_relevel(Outcome_rec, 
#                   "Alistipes", "Anaerostipes", "Bacteroides", "Anaerotruncus",
#                   "Bifidobacterium", "Blautia", "Escherichia and Shigella", "Butyricicoccus", 
#                   "Cellulosibacter", "Clostridium IV", "Clostridium sensu stricto",
#                   "Clostridium XlVa", "Clostridium XVIII", "Collinsella", "Coprococcus",
#                   "Dialister", "Dorea", "Eggerthella", "Eisenbergiella", "Enterobacter",
#                   "Enterococcus", "Erysipelotrichaceae incertae sedis", 
#                   "Faecalibacterium", "Flavonifractor", "Fusicatenibacter", "Gemmiger",
#                   "Granulicatella", "Haemophilus", "Hungatella", "Intestinibacter",
#                   "Klebsiella", "Lachnospiracea incertae sedis", "Lactococcus",
#                   "Oscillibacter", "Parabacteroides", "Peptoniphilus", "Romboutsia",
#                   "Roseburia", "Ruminococcus", "Ruminococcus 2", "Saccharibacteria genera incertae sedis",
#                   "Terrisporobacter","Streptococcus", "Subdoligranulum", "Veillonella", "Akkermansia")) %>%
#   ggplot(aes(x = -log10(`p-value`), y = Outcome_rec)) +
#   geom_point(aes(shape = sens_beta), size = 2) +
#   geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red") +
#   geom_vline(xintercept = -log10(0.05/(14*33)), linetype = "dashed", color = "blue") +
#   theme_lucid() +
#   labs(x = "-log10(P-value)", 
#        y = "Genera", 
#        shape = "") +
#   geom_text(aes(label = ifelse(`p-value` < 0.008, as.character(Pollutants_Time_window_rec), "")), hjust = -0.05, vjust = -0.3, angle = 35, size = 3.5) +
#   scale_shape_manual(values = c("Beta<0" = 15, "Beta≥0" = 17)) +# 15: carré plein, 17: triangle plein
#   theme(
#     legend.position = "right",
#     legend.box = "vertical", 
#     legend.justification = "right", 
#     axis.text.y = element_text(face = "italic"))

table_log$categorie <- factor(table_log$Outcome_rec, levels = c(levels(table_log$Outcome_rec), "  ", "   ", "    "))

df <- tibble(
  Phyla_corres = rep(NA, 3),
  Outcome = rep(NA, 3),
  Pollutants = rep(NA, 3),
  Time_window = rep(NA, 3),
  Pollutants_Time_window = rep(NA, 3),
  Pollutants_Time_window_rec = rep(NA, 3),
  Beta = rep(NA, 3),
  sens_beta = rep(NA, 3),
  `95% CI` = rep(NA, 3),
  lower_CI = rep(NA, 3),
  upper_CI = rep(NA, 3),
  `p-value` = rep(NA, 3),
  p_value_shape = rep(NA, 3),
  `q-value` = rep(NA, 3),
  q_value_shape = rep(NA, 3),
  Outcome_rec = rep(NA, 3),
  categorie = c("  ", "   ", "    ") # Assignation des valeurs avec des espaces
)


table_log <- rbind(table_log, df)

mahatan_plot <- table_log  %>%
  mutate(
    categorie = 
      fct_relevel(categorie, 
                  "Saccharibacteria genera incertae sedis", "Peptoniphilus",
                  "Granulicatella", "Anaerotruncus", "Lactococcus", "Terrisporobacter",
                  "Oscillibacter", "Haemophilus", "Coprococcus", "Erysipelotrichaceae incertae sedis",
                  "Butyricicoccus", "Dialister", "Subdoligranulum", "Intestinibacter",
                  "Klebsiella", "Eisenbergiella", "Hungatella", "Dorea", "Eggerthella",
                  "Romboutsia", "Clostridium IV", "Ruminococcus 2", "Flavonifractor",
                  "Alistipes", "Collinsella", "Parabacteroides", "Fusicatenibacter",
                  "Veillonella", "Roseburia", "Enterobacter", "Cellulosibacter",
                  "Enterococcus", "Clostridium sensu stricto", "Clostridium XVIII",
                  "Ruminococcus", "Gemmiger", "Anaerostipes", "Lachnospiracea incertae sedis",
                  "Clostridium XlVa", "Streptococcus", "Faecalibacterium", "Akkermansia",
                  "Escherichia and Shigella", "Blautia", "Bacteroides", "Bifidobacterium", "  ", "   ")) %>%
  ggplot(aes(x = -log10(`p-value`), y = categorie)) +
  geom_point(aes(shape = sens_beta), size = 2) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = -log10(0.05/(14*33)), linetype = "dashed", color = "blue") +
  theme_lucid() +
  labs(x = "-log10(P-value)", 
       y = "Genera", 
       shape = "") +
  geom_text(aes(label = ifelse(`p-value` < 0.008, as.character(Pollutants_Time_window_rec), "")), hjust = -0.07, vjust = -0.2, angle = 45, size = 3.5) +
  scale_shape_manual(values = c("Beta<0" = 15, "Beta≥0" = 17)) +# 15: carré plein, 17: triangle plein
  theme(
    legend.position = "bottom",
    legend.box = "vertical", 
    legend.justification = "center", 
    axis.text.y = element_text(face = "italic", size = 12)) +
  xlim(0, 5.8) 


mahatan_plot
ggsave("4_output/manhattan_plot_sensi.tiff", 
       mahatan_plot, 
       device = "tiff",
       units = "cm",
       dpi = 300,
       height = 27, 
       width = 27)

### Forestplot final ----
forest_plot <- table_log %>% 
  filter(`p-value`<0.00011) %>% 
  mutate(Beta = as.numeric(Beta)) %>%
  ggplot(aes(x = Outcome_rec, 
             y = Beta, 
             min = lower_CI, 
             ymax = upper_CI, 
             color = Pollutants_Time_window_rec)) +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_pointrange(
    position = position_dodge(width = 0.5), 
    size = 0.4) +
  labs(x = "Genera", y = "") +
  theme_lucid() +
  coord_flip()  +
  # facet_wrap(vars(Phyla_corres), scales = "free_y", ncol = 1) +
  # scale_shape_manual(values = c(19, 8),
  #                    name = "p-value") +
  guides(color = guide_legend(title = "Exposure and time window"))+
  theme(
    # axis.title = element_text(size = 7),
    #     axis.text = element_text(size = 6),
    #     legend.text = element_text(size = 7),
    #     legend.title = element_text(size = 7), 
    legend.position = "right",
    legend.box = "vertical", 
    legend.justification = "right", 
    axis.text.y = element_text(face = "italic")
    # legend.spacing.y = unit(0, "cm"), 
    # legend.spacing.x = unit(0, "cm"), 
    # legend.box.margin = margin(0,0,0,0, "cm"), 
    # legend.margin = margin(0,0,0,0, "cm")
  ) 

forest_plot
ggsave("4_output/forest_plot_genera.tiff", 
       forest_plot, 
       device = "tiff",
       units = "cm",
       dpi = 300,
       height = 15, 
       width = 25)

# Significant associations ----
## p<0.00011 ----
results_signi <-  table_log %>%
  filter(`p-value`<0.00011) 

write_xlsx(results_signi,
           path = "4_output/results_genera.xlsx")   # penser à copier coller le tbl_regression complet en plus

save.image("4_output/results_genera.RData")

## p<0.0035 ----
table_log %>% 
  select(-Pollutants, 
         -Time_window, 
         -Pollutants_Time_window, 
         -sens_beta, 
         -p_value_shape, 
         -q_value_shape, 
         -`q-value`, 
         -Outcome_rec, 
         -categorie ) %>%
  filter(`p-value`<0.0035) %>% 
  arrange(Phyla_corres, Class_corres, Order_corres, Family_corres, Outcome) %>%
  View()
