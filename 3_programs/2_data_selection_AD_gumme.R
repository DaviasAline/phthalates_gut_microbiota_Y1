# 2_data_selection
# A. Davias
# 19/10/2021


# Chargement des packages ----
library(tidyverse)
library(haven)
library(reshape2)
library(openxlsx)
library(GGally)
library(gtsummary)
library(summarytools)
library(patchwork)
library(AICcmodavg)
library(ggpubr)
library(grid)
library(questionr)
library(Hmisc)
library(rmarkdown)
library(knitr)
library(labelled)
library(distill)
library(rmdformats)
library(parameters)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(lazyeval)
library(car)
library(expss)
theme_gtsummary_language("en", decimal.mark = ".", big.mark = " ")
theme_gtsummary_compact(set_theme = TRUE)


# Chargement des données ----
load("1_intermediate_data/1_data_cleaning_AD_gumme.RData")

# Création des vecteurs ----
## vecteur phthalates ----
phthalates <- bdd %>% select(
  mo_DEHP_ms_i_cor_t2_ln, mo_DEHP_ms_i_cor_t3_ln, ch_DEHP_ms_i_cor_Y1_ln,       # Métabolites DEHP
  mo_MnBP_i_cor_t2_ln, mo_MnBP_i_cor_t3_ln, ch_MnBP_i_cor_Y1_ln,                # Métabolite DBP
  mo_DiNP_ms_i_cor_t2_ln, mo_DiNP_ms_i_cor_t3_ln, ch_DiNP_ms_i_cor_Y1_ln,       # Métabolites DiNP
  mo_MiBP_i_cor_t2_ln, mo_MiBP_i_cor_t3_ln, ch_MiBP_i_cor_Y1_ln,                # Métabolite DiBP
  mo_MBzP_i_cor_t2_ln, mo_MBzP_i_cor_t3_ln, ch_MBzP_i_cor_Y1_ln,                # Metabolite MBzP
  mo_MEP_i_cor_t2_ln, mo_MEP_i_cor_t3_ln, ch_MEP_i_cor_Y1_ln,                   # Metabolite DEP
  mo_ohMPHP_i_cor_t2_ln, mo_ohMPHP_i_cor_t3_ln, ch_ohMPHP_i_cor_Y1_ln,          # Metabolite DEHP
  mo_DINCH_ms_i_cor_t2_ln, mo_DINCH_ms_i_cor_t3_ln, ch_DINCH_ms_i_cor_Y1_ln     # Metabolites DINCH
) %>% colnames()

phthalates_pre <- bdd %>% select(all_of(phthalates)) %>% select(contains("t2"), contains("t3")) %>% colnames()
phthalates_post <- bdd %>% select(all_of(phthalates)) %>% select(contains("Y1")) %>% colnames()

rm(phthalates_total, 
   phthalates_total_ln, 
   phthalates_total_sg, 
   phthalates_total_sg_ln)

## vecteur covariables ----
covariates_pre <- bdd %>% select(                       
  "ch_feces_RUN_Y1",                   
  "ch_feces_age_w_Y1_i",
  "po_delmod",
  "ch_food_intro_Y1_3cat_i",
  "ch_antibio_Y1_2cat_i",
  "mo_par_2cat",
  "mo_pets_i",
  "ch_sex",
  "mo_tob_gr_anyt_yn_n2_i",
  "Mo_ETS_anyT_yn1_opt_i",
  "ch_ETS_12m_opt36m",
  "mo_interpreg_3cat",
  "mo_dipl_3cat_i",
  # "po_w_kg_3cat",
  # "po_he_3cat_i",
  # "ch_w_Y1_3cat_i",
  # "ch_he_Y1_3cat_i",
  # "po_gd",
  "mo_age",
  "mo_bmi_bepr_3cat_i",
  "bf_duration_till48w_4cat_i") %>% colnames()

covariates_post <- bdd %>% select(                       
  "ch_feces_RUN_Y1",                   
  "ch_feces_age_w_Y1_i",
  "po_delmod",
  "ch_food_intro_Y1_3cat_i",
  "ch_antibio_Y1_2cat_i",
  "mo_par_2cat",
  "mo_pets_i",
  "ch_sex",
  "mo_tob_gr_anyt_yn_n2_i",
  "Mo_ETS_anyT_yn1_opt_i",
  "ch_ETS_12m_opt36m",
  "mo_interpreg_3cat",
  "mo_dipl_3cat_i",
  "po_w_kg_3cat",
  "po_he_3cat_i",
  # "ch_w_Y1_3cat_i",
  # "ch_he_Y1_3cat_i",
  "po_gd",
  "mo_age",
  "mo_bmi_bepr_3cat_i",
  "bf_duration_till48w_4cat_i") %>% colnames()

## vecteur outcomes ----
alpha_vec <- c("ch_feces_SpecRich_5000_ASV_Y1", 
               "ch_feces_Shannon_5000_ASV_Y1") 
phyla_vec <- c("ch_feces_rel_p1_Y1", "ch_feces_rel_p2_Y1", 
               "ch_feces_rel_p3_Y1", "ch_feces_rel_p4_Y1")
genera_vec <- genera
genera_names <- gsub("Escherichia_Shigella", "Escherichia and Shigella", genera_vec)
genera_names <- gsub("_", " ", genera_names)
genera_names <- gsub("Ruminococcus2", "Ruminococcus 2", genera_names)
rm(genera)
outcomes <- c(alpha_vec, phyla_vec, genera_vec)
outcomes_names <- c("Specific richness", "Shannon diversity", 
                    "Firmicutes", "Actinobacteria", 
                    "Bacteroidetes", "Proteobacteria", 
                    genera_names)

# Labels de variables ----
bdd = modify(bdd,{
  var_lab(ident) = "SEPAGES dentity"
  var_lab(ch_feces_ID_Y1) = "Child age at the one-year stool collection"
  
  var_lab(ch_feces_age_w_Y1) = "Child age at stool collection (weeks)"    # variables dispo continues et catégorielles
  var_lab(ch_feces_age_w_Y1_4cat) = "Child age at stool collection"
  
  var_lab(ch_feces_age_w_Y1_i) = "Child age at stool collection (weeks)"   
  var_lab(ch_feces_age_w_Y1_4cat_i) = "Child age at stool collection"
  
  var_lab(po_gd) = "Gestational duration, completed weeks"    
  var_lab(po_gd_4cat) = "Gestational term, completed weeks" 
  
  var_lab(mo_age) = "Maternal age at conception (years)"                               
  var_lab(mo_age_4cat) = "Maternal age at conception"
  
  var_lab(mo_bmi_bepr) = "Maternal BMI before pregnancy (kg/m2)" 
  var_lab(mo_bmi_bepr_3cat) = "Maternal BMI before pregnancy" 
  var_lab(mo_bmi_bepr_i) = "Maternal BMI before pregnancy (kg/m2)" 
  var_lab(mo_bmi_bepr_3cat_i) = "Maternal BMI before pregnancy" 
  
  var_lab(mo_par) = "Maternal parity" 
  var_lab(mo_par_2cat) = "Maternal parity" 
  
  var_lab(po_w) = "Birth weight (g)" 
  var_lab(po_w_kg) = "Birth weight (kg)"
  var_lab(po_w_kg_3cat) = "Birth weight" 
  
  var_lab(po_he) = "Birth length (cm)" 
  var_lab(po_he_3cat) = "Birth length" 
  var_lab(po_he_i) = "Birth length (cm)" 
  var_lab(po_he_3cat_i) = "Birth length" 
  
  var_lab(ch_w_Y1) = "Weight at one year (kg)" 
  var_lab(ch_he_Y1) = "Length at one year (cm)"
  var_lab(ch_w_Y1_i) = "Weight at one year (kg)" 
  var_lab(ch_he_Y1_i) = "Length at one year (cm)"
  
  var_lab(ch_w_Y1_3cat) = "Weight at one year"
  var_lab(ch_he_Y1_3cat) = "Length at one year"
  var_lab(ch_w_Y1_3cat_i) = "Weight at one year"
  var_lab(ch_he_Y1_3cat_i) = "Length at one year"
  
  var_lab(bf_duration_till48w) = "Breastfeeding duration (weeks)"
  var_lab(bf_duration_till48w_4cat) = "Breastfeeding duration"
  var_lab(bf_duration_till48w_i) = "Breastfeeding duration (weeks)"
  var_lab(bf_duration_till48w_4cat_i) = "Breastfeeding duration"
  
  var_lab(ch_antibio_Y1) = "Number of antibiotics use between 0-12 months old"
  var_lab(ch_antibio_Y1_3cat) = "Number of antibiotics use between 0-12 months old"
  var_lab(ch_antibio_Y1_2cat) = "Antibiotics use between 0-12 months old"
  var_lab(ch_antibio_Y1_i) = "Number of antibiotics use between 0-12 months old"
  var_lab(ch_antibio_Y1_3cat_i) = "Number of antibiotics use between 0-12 months old"
  var_lab(ch_antibio_Y1_2cat_i) = "Antibiotics use between 0-12 months old"
  var_lab(ch_hospit_Y1) = "Hospitalization durin between 0-12 months old"
  
  var_lab(po_delmod) = "Delivery mode"              # variables forcément catégorielles 
  var_lab(ch_sex) = "Child sex"
  var_lab(mo_dipl) = "Maternal education" 
  var_lab(mo_dipl_3cat) = "Maternal education"
  var_lab(mo_dipl_3cat_i) = "Maternal education"
  var_lab(mo_pets) = "Presence of pets"  
  var_lab(mo_pets_i) = "Presence of pets"  
  var_lab(ch_food_intro_Y1) = "Period of introduction of solid food 0-1Y"
  var_lab(ch_food_intro_Y1_3cat) = "Period of introduction of solid food 0-1Y"
  var_lab(ch_food_intro_Y1_i) = "Period of introduction of solid food 0-1Y"
  var_lab(ch_food_intro_Y1_3cat_i) = "Period of introduction of solid food 0-1Y"
  var_lab(mo_tob_gr_anyt_yn_n2) = "Maternal active smoking during pregnancy"
  var_lab(Mo_ETS_anyT_yn1_opt) = "Maternal passive smoking during pregnancy"
  var_lab(ch_ETS_12m_opt36m) = "Child passive smoking during pregnancy"
  var_lab(mo_tob_gr_anyt_yn_n2_i) = "Maternal active smoking during pregnancy"
  var_lab(Mo_ETS_anyT_yn1_opt_i) = "Maternal passive smoking during pregnancy"
  var_lab(mo_interpreg) = "Interpregnancy interval (years)"
  var_lab(mo_interpreg_5cat) = "Interpregnancy interval"
  var_lab(mo_interpreg_3cat) = "Interpregnancy interval"
  
  var_lab(mo_ethnicity) = "Maternal ethnicity"                                    # varibales que l'on utilisera pas 
  var_lab(ch_sibling) = "Number of other children living at home (sibling or not)" 
  var_lab(fa_dipl) = "Paternal education"  
  var_lab(mo_cats) = "Presence of cats" 
  var_lab(mo_dogs) = "Presence of dogs"
  var_lab(mo_birds) = "Presence of birds" 
  var_lab(mo_rodents) = "Presence of rodents" 
  var_lab(mo_other_pets) = "Presence of other pets"
})

var_label(bdd[, phthalates]) <- colnames(bdd[, phthalates]) %>%
  str_replace_all(
    c("mo_" = "Maternal exposure to ",
      "ch_" = "Child exposure to ",
      "_ms_i_cor_" = " ", 
      "_i_cor_" = " ", 
      "t2" = "at trim.2",
      "t3" = "at trim.3",
      "Y1" = "at 12 months", 
      "_ln" = ""))

# Vérification codage ----
bdd %>% filter(statut == "inclu") %>% select(all_of(covariates_post))  %>% tbl_summary() 
bdd %>% filter(statut == "inclu") %>% select(all_of(phthalates))  %>% tbl_summary() 

save.image("1_intermediate_data/2_data_selection_AD_gumme.RData")
