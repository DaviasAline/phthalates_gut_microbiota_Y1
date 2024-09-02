# 1_data_cleaning
# Aline Davias 
# 2024_08_27

# 0_data_reading ----
## Chargement des packages ----
library(haven)
library(tidyverse)
library(questionr)
library(lubridate)
library(mice)
library(labelled)
library(readr)
library(expss)
library(gtsummary)
library(sjlabelled)

## Chargement des données ----
load("1_intermediate_data/0_source_data_reading_AD_gumme.RData")

# 1_data_cleaning covariates ----
## Data cleaning ----
bdd <- bdd %>%
  mutate(
    statut = ifelse(!is.na(ch_feces_OrderSeq_ASVbased_Y1), "inclu", "exclu"), 
    statut = as.factor(statut)) %>%
  rename(
    mo_ethnicity = mt1saa1_q02,
    ch_sibling = mt1saa1_q04,
    fa_dipl = ft2sac1_q03,
    mo_pets = mt3eaf1_q01,       # Variables animaux (questionnaire MT3EAF1_V1)
    mo_dogs = mt3eaf1_q01p1,
    mo_cats = mt3eaf1_q01p5,
    mo_birds = mt3eaf1_q01p9,
    mo_rodents = mt3eaf1_q01p13,
    mo_other_pets =  mt3eaf1_q01p17)


covar_vec_cat <- c("ch_sex", "po_delmod","mo_par", "mo_dipl", "mo_ethnicity","ch_sibling","fa_dipl","mo_pets",
                   "mo_dogs","mo_cats","mo_birds","mo_rodents","mo_other_pets","mo_tob_gr_anyt_yn_n2","mo_tob_gr_anyt_yn1_n2",
                   "ch_ETS_12m_opt36m","Mo_ETS_anyT_yn1_opt")
bdd[, covar_vec_cat] <- lapply(bdd[, covar_vec_cat], as.character)
bdd[, covar_vec_cat] <- lapply(bdd[, covar_vec_cat], as.factor)

bdd <- bdd %>%
  mutate(
    mo_tob_gr_anyt_yn_n2 = fct_recode(mo_tob_gr_anyt_yn_n2,            # tabagisme actif de la mère pdt la grossesse
                                      "No" = "0",
                                      "Yes" = "1"), 
    mo_tob_gr_anyt_yn1_n2 = fct_recode(mo_tob_gr_anyt_yn1_n2,
                                       "No" = "0",
                                       "Yes" = "1"), 
    ch_ETS_12m_opt36m = fct_recode(ch_ETS_12m_opt36m,                  # tabagisme passif de l'enfant 
                                   "No" = "0",
                                   "Yes" = "1"), 
    Mo_ETS_anyT_yn1_opt = fct_recode(Mo_ETS_anyT_yn1_opt,              # tabagisme passif de la mère pdt la grossesse
                                     "No" = "0",
                                     "Yes" = "1"), 
    ch_sex = fct_recode(ch_sex,                                        # Labellisation de ch_sex
                        "Male" = "1", 
                        "Female" = "2"), 
    po_delmod = fct_recode(po_delmod,                                  # Labellisation de po_delmod
                           "Vaginal delivery" = "1", 
                           "C-section" = "2"), 
    mo_par_2cat = fct_recode(mo_par,                                        # Labellisation de mo_par
                             "None" = "0", 
                             "1 child or more" = "1", 
                             "1 child or more" = "2"), 
    mo_ethnicity = fct_recode(mo_ethnicity,                            # Labellisation de mo_ethnicity
                              "Afrique" = "1",
                              "Amériques" = "2",
                              "Asie du Sud-Est" = "3",
                              "Europe" = "4",
                              "Méditéranée Orientale" = "5",
                              "Pacifique Occidental" = "6",            
                              "Autre" = "7", 
                              "Ne sait pas/ne souhaite pas répondre" = "99"), 
    mo_dipl = fct_recode(mo_dipl,                                      # Labellisation des variables éducation mère mo_dipl
                         "BEP/CAP/Highschool" = "2", 
                         "1-2 years after graduation" = "3", 
                         "3-4 years after graduation" = "4", 
                         "≥5 years after graduation" = "5"), 
    fa_dipl = fct_recode(fa_dipl,                                      # Labellisation des variables éducation père fa_dipl              
                         "BEP/CAP/Highschool" = "1", 
                         "BEP/CAP/Highschool" = "2", 
                         "BEP/CAP/Highschool" = "3", 
                         "BEP/CAP/Highschool" = "4", 
                         "BEP/CAP/Highschool" = "5", 
                         "1-2years after graduation" = "6",
                         "3-4years after graduation" = "7",
                         "≥5years after graduation" = "8"), 
    mo_pets = fct_recode(mo_pets,                                      # Labellisation des variables animaux 
                         "No" = "0", 
                         "One or more" = "1"))

# warning message mo_ethnicity : Unknown levels in `f`: 6  : normal car absence d'observation pour cette modalité 
# warning message fa_dipl : Unknown levels in `f`: 1, 2  : normal car absence d'observation pour ces modalités 

## Création de la variable hospit 0 - 1 an ----
bdd <- bdd %>% 
  mutate(
    ch_hospit_Y1 = ifelse(ch_hospit_2M | ch_hospit_2_12M > 0, 1, 0))

## Création variable intervalle intergrossesse ----
bdd <- bdd %>% 
  rename(
    mo_number_previous_pregnancy = mt1haa1_q22,                   # Nombre de grossesse antérieures 
    mo_first_pregnancy = mt1haa1_q22p1,                           # Infos intervalle avec 1ère grossesse
    mo_month_first_pregnancy = mt1haa1_q22p3, 
    mo_year_first_pregnancy = mt1haa1_q22p4,        
    mo_second_pregnancy = mt1haa1_q22p19,                        # Infos intervalle avec 2nd grossesse
    mo_month_second_pregnancy = mt1haa1_q22p21, 
    mo_year_second_pregnancy = mt1haa1_q22p22, 
    mo_third_pregnancy = mt1haa1_q22p37,                         # Infos intervalle avec 3eme grossesse
    mo_month_third_pregnancy = mt1haa1_q22p39, 
    mo_year_third_pregnancy = mt1haa1_q22p40,
    mo_fourth_pregnancy = mt1haa1_q22p55,                        # Infos intervalle avec 4eme grossesse
    mo_month_fourth_pregnancy = mt1haa1_q22p57, 
    mo_year_fourth_pregnancy = mt1haa1_q22p58,
    mo_fifth_pregnancy = mt1haa1_q22p73,                         # Infos intervalle avec 5eme grossesse
    mo_month_fifth_pregnancy = mt1haa1_q22p75, 
    mo_year_fifth_pregnancy = mt1haa1_q22p76,
    mo_sixth_pregnancy = mt1haa1_q22p91,                         # Infos intervalle avec 6eme grossesse
    mo_month_sixth_pregnancy = mt1haa1_q22p93, 
    mo_year_sixth_pregnancy = mt1haa1_q22p94) %>%
  mutate(
    mo_number_previous_pregnancy = as.double(mo_number_previous_pregnancy), 
    mo_first_pregnancy = as.double(mo_first_pregnancy), 
    mo_month_first_pregnancy = as.double(mo_month_first_pregnancy), 
    mo_year_first_pregnancy = as.double(mo_year_first_pregnancy), 
    mo_second_pregnancy = as.double(mo_second_pregnancy), 
    mo_month_second_pregnancy = as.double(mo_month_second_pregnancy), 
    mo_year_second_pregnancy = as.double(mo_year_second_pregnancy),
    mo_third_pregnancy = as.double(mo_third_pregnancy), 
    mo_month_third_pregnancy = as.double(mo_month_third_pregnancy), 
    mo_year_third_pregnancy = as.double(mo_year_third_pregnancy), 
    mo_fourth_pregnancy = as.double(mo_fourth_pregnancy), 
    mo_month_fourth_pregnancy = as.double(mo_month_fourth_pregnancy), 
    mo_year_fourth_pregnancy = as.double(mo_year_fourth_pregnancy), 
    mo_fifth_pregnancy = as.double(mo_fifth_pregnancy), 
    mo_month_fifth_pregnancy = as.double(mo_month_fifth_pregnancy), 
    mo_year_fifth_pregnancy = as.double(mo_year_fifth_pregnancy), 
    mo_sixth_pregnancy = as.double(mo_sixth_pregnancy), 
    mo_month_sixth_pregnancy = as.double(mo_month_sixth_pregnancy), 
    mo_year_sixth_pregnancy = as.double(mo_year_sixth_pregnancy)
  )


## Pour chaque variable mois codé en NA ou en 99 alors que l'on a une valeur pour l'année, on impute par la valeur 1
bdd <- bdd %>% 
  mutate(
    mo_month_first_pregnancy = ifelse(mo_month_first_pregnancy < 13, mo_month_first_pregnancy, 1),
    mo_month_second_pregnancy = ifelse(mo_month_second_pregnancy < 13, mo_month_second_pregnancy, 1),
    mo_month_third_pregnancy = ifelse(mo_month_third_pregnancy < 13, mo_month_third_pregnancy, 1),
    mo_month_fourth_pregnancy = ifelse(mo_month_fourth_pregnancy < 13, mo_month_fourth_pregnancy, 1),
    mo_month_fifth_pregnancy = ifelse(mo_month_fifth_pregnancy < 13, mo_month_fifth_pregnancy, 1),
    mo_month_sixth_pregnancy = ifelse(mo_month_sixth_pregnancy < 13, mo_month_sixth_pregnancy, 1)) %>% 
  mutate(                                                       
    mo_month_first_pregnancy = ifelse((!is.na(mo_year_first_pregnancy)) & (is.na(mo_month_first_pregnancy)),
                                       1, 
                                       mo_month_first_pregnancy), 
    mo_month_second_pregnancy = ifelse((!is.na(mo_year_second_pregnancy)) & (is.na(mo_month_second_pregnancy)),
                                        1, 
                                        mo_month_second_pregnancy), 
    mo_month_third_pregnancy = ifelse((!is.na(mo_year_third_pregnancy)) & (is.na(mo_month_third_pregnancy)),
                                       1, 
                                       mo_month_third_pregnancy), 
    mo_month_fourth_pregnancy = ifelse((!is.na(mo_year_fourth_pregnancy)) & (is.na(mo_month_fourth_pregnancy)),
                                        1, 
                                        mo_month_fourth_pregnancy), 
    mo_month_fifth_pregnancy = ifelse((!is.na(mo_year_fifth_pregnancy)) & (is.na(mo_month_fifth_pregnancy)),
                                       1, 
                                       mo_month_fifth_pregnancy), 
    mo_month_sixth_pregnancy = ifelse((!is.na(mo_year_sixth_pregnancy)) & (is.na(mo_month_sixth_pregnancy)),
                                       1, mo_month_sixth_pregnancy)) %>%
  
  mutate(                               # Rassembler le mois/l'année pour chaque grossesse
    mo_date_first_pregnancy = str_c(mo_year_first_pregnancy, mo_month_first_pregnancy, sep="-"), 
    mo_date_second_pregnancy = str_c(mo_year_second_pregnancy, mo_month_second_pregnancy, sep="-"),
    mo_date_third_pregnancy = str_c(mo_year_third_pregnancy, mo_month_third_pregnancy, sep="-"), 
    mo_date_fourth_pregnancy = str_c(mo_year_fourth_pregnancy, mo_month_fourth_pregnancy, sep="-"), 
    mo_date_fifth_pregnancy = str_c(mo_year_fifth_pregnancy, mo_month_fifth_pregnancy, sep="-"), 
    mo_date_sixth_pregnancy = str_c(mo_year_sixth_pregnancy, mo_month_sixth_pregnancy, sep="-")) %>% 
  
  mutate(                                                      # Conversion en variable date 
    mo_date_first_pregnancy = ym(mo_date_first_pregnancy), 
    mo_date_second_pregnancy = ym(mo_date_second_pregnancy), 
    mo_date_third_pregnancy = ym(mo_date_third_pregnancy), 
    mo_date_fourth_pregnancy = ym(mo_date_fourth_pregnancy), 
    mo_date_fifth_pregnancy = ym(mo_date_fifth_pregnancy), 
    mo_date_sixth_pregnancy = ym(mo_date_sixth_pregnancy))


## Création de la variable date de la dernière grossesse ayant conduit à un enfant vivant 
bdd <- bdd %>%
  mutate(
    mo_date_last_pregnancy = NA
    # , 
    # mo_date_last_pregnancy = as.numeric(mo_date_last_pregnancy)
    ) %>%
  mutate(mo_date_last_pregnancy = if_else(mo_sixth_pregnancy == 1, 
                                          mo_date_sixth_pregnancy, 
                                          mo_date_last_pregnancy))%>%
  mutate(mo_date_last_pregnancy = if_else(mo_fifth_pregnancy == 1 & (is.na(mo_date_last_pregnancy)), 
                                          mo_date_fifth_pregnancy, 
                                          mo_date_last_pregnancy)) %>%
  mutate(mo_date_last_pregnancy = if_else(mo_fourth_pregnancy == 1 & (is.na(mo_date_last_pregnancy)), 
                                          mo_date_fourth_pregnancy, 
                                          mo_date_last_pregnancy)) %>%
  mutate(mo_date_last_pregnancy = if_else(mo_third_pregnancy == 1 & (is.na(mo_date_last_pregnancy)), 
                                          mo_date_third_pregnancy, 
                                          mo_date_last_pregnancy)) %>%
  mutate(mo_date_last_pregnancy = if_else(mo_second_pregnancy == 1 & (is.na(mo_date_last_pregnancy)), 
                                          mo_date_second_pregnancy, 
                                          mo_date_last_pregnancy)) %>%
  mutate(mo_date_last_pregnancy = if_else(mo_first_pregnancy == 1 & (is.na(mo_date_last_pregnancy)), 
                                          mo_date_first_pregnancy, 
                                          mo_date_last_pregnancy))
bdd <- bdd %>%
  mutate(
    mo_number_previous_pregnancy = as.factor(mo_number_previous_pregnancy),
    mo_interpreg = difftime(po_datedel, 
                            mo_date_last_pregnancy, 
                            units = "weeks"),
    mo_interpreg = as.numeric(mo_interpreg), 
    mo_interpreg = mo_interpreg / 52.1429,             # Transformation en année
    mo_interpreg_5cat = cut(mo_interpreg,              # création variable catégorielle 
                            include.lowest = TRUE,
                            right = FALSE,
                            dig.lab = 4,
                            breaks = c(1, 2, 3, 4, 11.9945106894433)), 
    mo_interpreg_5cat = fct_recode(mo_interpreg_5cat,
                                   "Under 2 years" = "[1,2)",
                                   "Between 2 and 3 years" = "[2,3)",
                                   "Between 3 and 4 years" = "[3,4)",
                                   "Over 4 years" = "[4,11.99]"), 
    mo_interpreg_5cat = as.character(mo_interpreg_5cat), 
    mo_par_2cat = as.character(mo_par_2cat), 
    mo_interpreg_5cat = ifelse(mo_number_previous_pregnancy == 0 & (is.na(mo_interpreg_5cat)),   # Prise en compte des NA qui sont en fait des femmes primipares
                                "Primiparous", 
                                mo_interpreg_5cat), 
    mo_interpreg_5cat = ifelse( mo_par_2cat == "None" & (is.na(mo_interpreg_5cat)), 
                                 "Primiparous", 
                                 mo_interpreg_5cat), 
    mo_interpreg_5cat = as.factor(mo_interpreg_5cat), 
    mo_par_2cat = as.factor( mo_par_2cat)) %>%
  mutate(
    mo_interpreg_5cat = fct_relevel(mo_interpreg_5cat,
                                    "Under 2 years", 
                                    "Between 2 and 3 years", 
                                    "Between 3 and 4 years", 
                                    "Over 4 years", 
                                    "Primiparous"),
    mo_interpreg_3cat = fct_recode(mo_interpreg_5cat,
                                   "2 years and more" = "Between 2 and 3 years",
                                   "2 years and more" = "Between 3 and 4 years",
                                   "2 years and more" = "Over 4 years"))


## Choix du codage des covariables ----
bdd <- bdd %>%                                    # création de variables catégorielles à partir de variables numériques 
  mutate(
    ch_feces_age_w_Y1 = as.numeric(ch_feces_age_w_Y1),
    ch_feces_age_w_Y1_4cat = cut(
      ch_feces_age_w_Y1,
      include.lowest = TRUE,
      right = FALSE,
      dig.lab = 4,
      breaks = c(42, 51, 53, 55, 71)),
    po_w_kg = po_w / 1000,
    po_w_kg_3cat = cut(po_w_kg,
                       include.lowest = TRUE,
                       right = FALSE,
                       dig.lab = 4,
                       breaks = c(0.9, 3, 3.5, 4.7)), 
    po_he_3cat = cut(po_he,
                     include.lowest = TRUE,
                     right = FALSE,
                     dig.lab = 4,
                     breaks = c(30, 50, 52, 60)),
    ch_w_Y1_3cat = cut(ch_w_Y1,
                       include.lowest = TRUE,
                       right = FALSE,
                       dig.lab = 4,
                       breaks = c(6.0, 8.5, 10, 14)),
    ch_he_Y1_3cat = cut(ch_he_Y1,
                        include.lowest = TRUE,
                        right = FALSE,
                        dig.lab = 4,
                        breaks = c(67.75, 75, 78, 83.75)),
    mo_bmi_bepr_3cat = cut(mo_bmi_bepr,
                           include.lowest = TRUE, 
                           right = FALSE,
                           dig.lab = 4,
                           breaks = c(16, 19, 24, 42)),
    po_gd_4cat = cut(po_gd,                                                 
                     include.lowest = TRUE,
                     right = FALSE,
                     dig.lab = 4,
                     breaks = c(28, 38, 40, 41, 42)), 
    mo_age_4cat = cut(mo_age,                                                
                      include.lowest = TRUE,
                      right = FALSE,
                      dig.lab = 4,
                      breaks = c(20, 27, 33, 36, 46)), 
    bf_duration_till48w_4cat = cut(bf_duration_till48w,                         
                                   include.lowest = TRUE,
                                   right = FALSE,
                                   dig.lab = 4,
                                   breaks = c(0, 1, 24, 47, 48))) %>%
  
  
  mutate(
    # création d'étiquettes de catégories pour les nouvelles variables
    ch_feces_age_w_Y1_4cat = fct_recode(
      ch_feces_age_w_Y1_4cat,
      "<51 weeks" = "[42,51)",
      "51-52 weeks" = "[51,53)",
      "53-54 weeks" = "[53,55)",
      "≥55 weeks" = "[55,71]"),
    ch_antibio_Y1 = as.character(ch_antibio_Y1),
    ch_antibio_Y1_3cat = fct_recode(
      ch_antibio_Y1,
      "2 and more" = "2",
      "2 and more" = "3",
      "2 and more" = "4",
      "2 and more" = "5"),
    ch_antibio_Y1_2cat = fct_recode(
      ch_antibio_Y1,
      "No" = "0",
      "Yes" = "1",
      "Yes" = "2",
      "Yes" = "3",
      "Yes" = "4",
      "Yes" = "5"),
    ch_food_intro_Y1_3cat = fct_recode(
      ch_food_intro_Y1,
      "Between 0 and 6 months old" = "Between 3 and 6 months old",
      "Between 0 and 6 months old" = "Between 0 and 3 months old"),
    po_w_kg_3cat = fct_recode(
      po_w_kg_3cat,
      "<3 kg" = "[0.9,3)",
      "3-3.4 kg" = "[3,3.5)",
      "≥3.5 kg" = "[3.5,4.7]"),
    po_he_3cat = fct_recode(
      po_he_3cat,
      "<50 cm" = "[30,50)",
      "50-51 cm" = "[50,52)",
      "≥52 cm" = "[52,60]"),
    ch_w_Y1_3cat = fct_recode(
      ch_w_Y1_3cat,
      "<8.5 kg" = "[6,8.5)",
      "8.5-9.9 kg" = "[8.5,10)",
      "≥10 kg" = "[10,14]"),
    ch_he_Y1_3cat = fct_recode(
      ch_he_Y1_3cat,
      "<75 cm" = "[67.75,75)",
      "75-77.9 cm" = "[75,78)",
      "≥78 cm" = "[78,83.75]"),
    mo_bmi_bepr_3cat = fct_recode(
      mo_bmi_bepr_3cat,
      "<19 kg/m2" = "[16,19)",
      "19-23.9 kg/m2" = "[19,24)",
      "≥24 kg/m2" = "[24,42]"),
    po_gd_4cat = fct_recode(
      po_gd_4cat,
      "<38 weeks of amenorrhea" = "[28,38)",
      "38-39 weeks of amenorrhea" = "[38,40)",
      "40 weeks of amenorrhea" = "[40,41)",
      ">40 weeks of amenorrhea" = "[41,42]"),
    mo_age_4cat = fct_recode(
      mo_age_4cat,
      "20-26 years" = "[20,27)",
      "27-32 years" = "[27,33)",
      "33-35 years" = "[33,36)",
      ">35 years" = "[36,46]"),
    bf_duration_till48w_4cat = fct_recode(
      bf_duration_till48w_4cat,
      "Not breastfed" = "[0,1)",
      "<24 weeks" = "[1,24)",
      "24-47 weeks" = "[24,47)",
      "Still breastfeed at 48 weeks" = "[47,48]"),
    mo_interpreg_3cat = fct_relevel(
      mo_interpreg_3cat,
      "Primiparous",
      "2 years and more",
      "Under 2 years"),
    mo_dipl = as.character(mo_dipl),
    mo_dipl_3cat = fct_recode(
      mo_dipl,
      "2 years or less after graduation" = "BEP/CAP/Highschool",
      "2 years or less after graduation" = "1-2 years after graduation"))


## Imputation des valeurs manquantes covariables ----
bdd %>%
  filter(statut == "inclu") %>%
  select(                 
    ch_feces_age_w_Y1,     # 7 NA /356
    po_gd,                 # 0 NA /356
    mo_age,                # 0 NA /356
    mo_bmi_bepr,           # 1 NA /356
    mo_par,                # 0 NA /356
    po_w_kg,               # 0 NA /356
    po_he,                 # 1 NA /356
    ch_w_Y1,               # 5 NA /356
    ch_he_Y1,              # 13 NA /356
    bf_duration_till48w,   # 9 NA /356
    ch_antibio_Y1,         # 1 NA /356
    po_delmod,             # 0 NA /356
    ch_sex,                # 0 NA /356
    mo_dipl_3cat,          # 2 NA /356
    mo_pets,               # 45 NA /356
    mo_interpreg_3cat,     # 0 NA /356
    mo_tob_gr_anyt_yn_n2,  # 26 NA /356
    Mo_ETS_anyT_yn1_opt,   # 17 NA /356
    ch_ETS_12m_opt36m,     # 0 NA /356
    ch_food_intro_Y1,      # 37 NA /356
    ch_hospit_Y1,          # 9 NA / 356
    ch_hospit_2M,          # 3 NA / 356
    ch_hospit_2_12M        # 8 NA / 356
  ) %>%
  tbl_summary()


### Imputations (package mice)
bdd_imput <- bdd %>%
  mutate(
    ch_antibio_Y1 = as.numeric(ch_antibio_Y1), 
    ch_hospit_Y1 = as.factor(ch_hospit_Y1),
    ch_hospit_2M = as.factor(ch_hospit_2M),
    ch_hospit_2_12M = as.factor(ch_hospit_2_12M),
    ch_food_intro_Y1 = as.factor(ch_food_intro_Y1)) %>%
  select(ident, 
         ch_feces_age_w_Y1,                                                     # variables à imputer en continu
         mo_bmi_bepr,
         po_he, 
         ch_w_Y1, 
         ch_he_Y1,
         bf_duration_till48w,
         ch_antibio_Y1,
         
         mo_dipl_3cat,                                                          # variables à imputer en catégoriel (predictive mean matching)
         mo_pets,
         mo_tob_gr_anyt_yn_n2,
         Mo_ETS_anyT_yn1_opt, 
         ch_food_intro_Y1, 
         ch_hospit_Y1, 
         ch_hospit_2M, 
         ch_hospit_2_12M)

str(bdd_imput)                                                                  # vérifier le codage des variables 

imput <-  mice(bdd_imput, m=1, maxit=50, seed=500)                              # imputation

imput$method                                                                    # voir quelles méthodes ont été utilisés pour chaque variable à imputer

imput$imp$ch_feces_age_w_Y1                                                     # voir les imputations 
imput$imp$mo_bmi_bepr
imput$imp$po_he 
imput$imp$ch_w_Y1
imput$imp$ch_he_Y1
imput$imp$bf_duration_till48w
imput$imp$ch_antibio_Y1
imput$imp$mo_dipl_3cat  
imput$imp$mo_pets
imput$imp$mo_tob_gr_anyt_yn_n2
imput$imp$Mo_ETS_anyT_yn1_opt 
imput$imp$ch_food_intro_Y1

bdd_imput <- complete(imput)                                                    # extrair la bdd imputée
bdd_imput <- as.data.frame(bdd_imput) 
colnames(bdd_imput) <-                                                          # les nouvelles variables imputées se terminent en _i
  c("ident", paste(colnames(bdd_imput[, 2:16]), "i", sep="_"))  

bdd <- left_join(bdd,                                                           # merger les variables imputées au reste de la base de données
                 bdd_imput,
                 by = "ident")  

# on suprime les valeurs imputées créées chez des enfants qui n'ont pas eu de prélévement de selles à 1 an
bdd <- bdd %>%
  mutate(ch_feces_age_w_Y1_i = 
           ifelse(statut == "inclu", ch_feces_age_w_Y1_i, ch_feces_age_w_Y1))   
# parmis les "exclu" de l'analyse, il y a quand même 4 valeurs pour la variable age au prélévelement de la selle 
# --> correspond à des échantillons avec une faible quantité d'adn donc exclu des analyses

rm(imput, bdd_imput)

# on re créé les variables catégorielles imputées
bdd <- bdd %>%
  mutate(
    ch_feces_age_w_Y1_4cat_i = cut(
      ch_feces_age_w_Y1_i,
      include.lowest = TRUE,
      right = FALSE,
      dig.lab = 4,
      breaks = c(42, 51, 53, 55, 71),
    ),
    mo_bmi_bepr_3cat_i = cut(
      mo_bmi_bepr_i,
      include.lowest = TRUE,
      right = FALSE,
      dig.lab = 4,
      breaks = c(16, 19, 24, 42)
    ),
    po_he_3cat_i = cut(
      po_he_i,
      include.lowest = TRUE,
      right = FALSE,
      dig.lab = 4,
      breaks = c(30, 50, 52, 60)
    ),
    ch_w_Y1_3cat_i = cut(
      ch_w_Y1_i,
      include.lowest = TRUE,
      right = FALSE,
      dig.lab = 4,
      breaks = c(6.0, 8.5, 10, 14)
    ),
    ch_he_Y1_3cat_i = cut(
      ch_he_Y1_i,
      include.lowest = TRUE,
      right = FALSE,
      dig.lab = 4,
      breaks = c(67.75, 75, 78, 83.75)
    ),
    bf_duration_till48w_4cat_i = cut(
      bf_duration_till48w_i,
      include.lowest = TRUE,
      right = FALSE,
      dig.lab = 4,
      breaks = c(0, 1, 24, 47, 48)
    )
  ) %>%
  
  mutate(                                     # création d'étiquettes de catégories pour les nouvelles variables
    ch_feces_age_w_Y1_4cat_i = fct_recode(
      ch_feces_age_w_Y1_4cat_i,
      "<51 weeks" = "[42,51)",
      "51-52 weeks" = "[51,53)",
      "53-54 weeks" = "[53,55)",
      "≥55 weeks" = "[55,71]"
    ),
    mo_bmi_bepr_3cat_i = fct_recode(
      mo_bmi_bepr_3cat_i,
      "<19 kg/m2" = "[16,19)",
      "19-23.9 kg/m2" = "[19,24)",
      "≥24 kg/m2" = "[24,42]"
    ),
    po_he_3cat_i = fct_recode(
      po_he_3cat_i,
      "<50 cm" = "[30,50)",
      "50-51 cm" = "[50,52)",
      "≥52 cm" = "[52,60]"
    ),
    ch_w_Y1_3cat_i = fct_recode(
      ch_w_Y1_3cat_i,
      "<8.5 kg" = "[6,8.5)",
      "8.5-9.9 kg" = "[8.5,10)",
      "≥10 kg" = "[10,14]"
    ),
    ch_he_Y1_3cat_i = fct_recode(
      ch_he_Y1_3cat_i,
      "<75 cm" = "[67.75,75)",
      "75-77.9 cm" = "[75,78)",
      "≥78 cm" = "[78,83.75]"
    ),
    bf_duration_till48w_4cat_i = fct_recode(
      bf_duration_till48w_4cat_i,
      "Not breastfed" = "[0,1)",
      "<24 weeks" = "[1,24)",
      "24-47 weeks" = "[24,47)",
      "Still breastfeed at 48 weeks" = "[47,48]"
    ),
    ch_antibio_Y1_3cat_i = fct_recode(
      as.factor(as.character(ch_antibio_Y1_i)),
      "2 and more" = "2",
      "2 and more" = "3",
      "2 and more" = "4",
      "2 and more" = "5"
    ),
    ch_antibio_Y1_2cat_i = fct_recode(
      as.factor(as.character(ch_antibio_Y1_i)),
      "No" = "0",
      "Yes" = "1",
      "Yes" = "2",
      "Yes" = "3",
      "Yes" = "4",
      "Yes" = "5"
    ),
    ch_food_intro_Y1_3cat_i = fct_recode(
      ch_food_intro_Y1_i,
      "Between 0 and 6 months old" = "Between 3 and 6 months old",
      "Between 0 and 6 months old" = "Between 0 and 3 months old"
    )
  )

## Ordre des catégories des variables catégorielles ----
bdd[, covar_vec_cat] <- lapply(bdd[, covar_vec_cat], as.factor) 

bdd <- bdd %>% 
  mutate(
    ch_food_intro_Y1_3cat_i = fct_relevel(
      ch_food_intro_Y1_3cat_i,
      "Between 0 and 6 months old",
      "Between 6 and 12 months old",
      "Not introduced at 12 months old"), 
    ch_food_intro_Y1_3cat = fct_relevel(
      ch_food_intro_Y1_3cat,
      "Between 0 and 6 months old",
      "Between 6 and 12 months old",
      "Not introduced at 12 months old"), 
    
    mo_par_2cat = fct_relevel(
      mo_par_2cat,
      "None", "1 child or more"),
    
    mo_interpreg_3cat = fct_relevel(
      mo_interpreg_3cat,
      "Under 2 years", "2 years and more", "Primiparous"),
    
    mo_dipl_3cat_i = fct_relevel(
      mo_dipl_3cat_i,
      "2 years or less after graduation", 
      "3-4 years after graduation", 
      "≥5 years after graduation"),
    mo_dipl_3cat = fct_relevel(
      mo_dipl_3cat,
      "2 years or less after graduation", 
      "3-4 years after graduation", 
      "≥5 years after graduation"),
    
    po_w_kg_3cat = fct_relevel(
      po_w_kg_3cat,
      "<3 kg", "3-3.4 kg", "≥3.5 kg"),
    
    po_he_3cat_i = fct_relevel(
      po_he_3cat_i,
      "<50 cm", "50-51 cm", "≥52 cm"),
    po_he_3cat = fct_relevel(
      po_he_3cat,
      "<50 cm", "50-51 cm", "≥52 cm"),
    
    ch_w_Y1_3cat_i = fct_relevel(
      ch_w_Y1_3cat_i,
      "<8.5 kg", "8.5-9.9 kg", "≥10 kg"),
    ch_w_Y1_3cat = fct_relevel(
      ch_w_Y1_3cat,
      "<8.5 kg", "8.5-9.9 kg", "≥10 kg"),
    
    ch_he_Y1_3cat_i = fct_relevel(
      ch_he_Y1_3cat_i,
      "<75 cm", "75-77.9 cm", "≥78 cm"),
    ch_he_Y1_3cat = fct_relevel(
      ch_he_Y1_3cat,
      "<75 cm", "75-77.9 cm", "≥78 cm"),
    
    ch_w_Y1_3cat_i = fct_relevel(
      ch_w_Y1_3cat_i,
      "<8.5 kg", "8.5-9.9 kg", "≥10 kg"),
    
    mo_bmi_bepr_3cat_i = fct_relevel(
      mo_bmi_bepr_3cat_i,
      "<19 kg/m2", "19-23.9 kg/m2", "≥24 kg/m2"),
    mo_bmi_bepr_3cat = fct_relevel(
      mo_bmi_bepr_3cat,
      "<19 kg/m2", "19-23.9 kg/m2", "≥24 kg/m2"),
    
    bf_duration_till48w_4cat_i = fct_relevel(
      bf_duration_till48w_4cat_i,
      "Not breastfed", 
      "<24 weeks", 
      "24-47 weeks", 
      "Still breastfeed at 48 weeks"), 
    bf_duration_till48w_4cat = fct_relevel(
      bf_duration_till48w_4cat,
      "Not breastfed", 
      "<24 weeks", 
      "24-47 weeks", 
      "Still breastfeed at 48 weeks"))

rm(covar_vec_cat)

# 2_data_cleaning expo ----
phthalates_total <- bdd %>% select(
  
  mo_MEOHP_i_cor_t2, mo_MEOHP_i_cor_t3, ch_MEOHP_i_cor_Y1, 
  mo_MECPP_i_cor_t2, mo_MECPP_i_cor_t3, ch_MECPP_i_cor_Y1, 
  mo_MEHHP_i_cor_t2, mo_MEHHP_i_cor_t3, ch_MEHHP_i_cor_Y1,    
  mo_MEHP_i_cor_t2, mo_MEHP_i_cor_t3, ch_MEHP_i_cor_Y1, 
  mo_MMCHP_i_cor_t2, mo_MMCHP_i_cor_t3, ch_MMCHP_i_cor_Y1, 
  mo_DEHP_ms_i_cor_t2, mo_DEHP_ms_i_cor_t3, ch_DEHP_ms_i_cor_Y1,                # Métabolites DEHP
  
  mo_MnBP_i_cor_t2, mo_MnBP_i_cor_t3, ch_MnBP_i_cor_Y1,                         # Métabolite DBP
  
  mo_ohMiNP_i_cor_t2, mo_ohMiNP_i_cor_t3, ch_ohMiNP_i_cor_Y1, 
  mo_oxoMiNP_i_cor_t2, mo_oxoMiNP_i_cor_t3, ch_oxoMiNP_i_cor_Y1, 
  mo_cxMiNP_i_cor_t2, mo_cxMiNP_i_cor_t3, ch_cxMiNP_i_cor_Y1, 
  mo_DiNP_ms_i_cor_t2, mo_DiNP_ms_i_cor_t3, ch_DiNP_ms_i_cor_Y1,                # Métabolites DiNP
  
  mo_MiBP_i_cor_t2, mo_MiBP_i_cor_t3, ch_MiBP_i_cor_Y1,                         # Métabolite DiBP
  mo_MBzP_i_cor_t2, mo_MBzP_i_cor_t3, ch_MBzP_i_cor_Y1,                         # Metabolite MBzP
  mo_MEP_i_cor_t2, mo_MEP_i_cor_t3, ch_MEP_i_cor_Y1,                            # Metabolite DEP
  mo_ohMPHP_i_cor_t2, mo_ohMPHP_i_cor_t3, ch_ohMPHP_i_cor_Y1,                   # Metabolite DPHP
  
  mo_ohMINCH_i_cor_t2, mo_ohMINCH_i_cor_t3, ch_ohMINCH_i_cor_Y1, 
  mo_oxoMINCH_i_cor_t2, mo_oxoMINCH_i_cor_t3, ch_oxoMINCH_i_cor_Y1,
  mo_DINCH_ms_i_cor_t2, mo_DINCH_ms_i_cor_t3, ch_DINCH_ms_i_cor_Y1) %>%         # Metabolites DINCH
  
  colnames()

phthalates_total_ln <- paste(phthalates_total, "ln", sep = "_")
phthalates_total_sg <- gsub("_i_cor_", "_i_cor_sg_", phthalates_total)
phthalates_total_sg_ln <- gsub("_i_cor_", "_i_cor_sg_", phthalates_total_ln)

## Création phthalates ln ----
bdd[, phthalates_total] <-                                                      # on fait en sorte d'etre sur que toutes les variables sont codées correctement
  lapply(bdd[, phthalates_total], as.numeric)
data_phthalates <- bdd %>%                                                      # on créé une bdd spéciale phthalates
  select(ident, all_of(phthalates_total))
data_phthalates[, phthalates_total] <-                                          # on log transforme
  log(data_phthalates[, phthalates_total]) 
colnames(data_phthalates) <-                                                    # on ajoute "_ln" à chaque fin de nom de variable
  c("ident", paste(colnames(data_phthalates[, 2:55]), "ln", sep="_")) 
bdd <- left_join(bdd,                                                           # on merge avec la bdd principale
                 data_phthalates, 
                 by ="ident")
rm(data_phthalates)

## Création phthalates ln ajustées sur gravité spécifique ----
bdd[, phthalates_total_sg] <-                                                   # on fait en sorte d'etre sur que toutes les variables sont codées correctement
  lapply(bdd[, phthalates_total_sg], as.numeric)
data_phthalates_sg <- bdd %>%                                                   # on créé une bdd spéciale phthalates
  select(ident, all_of(phthalates_total_sg))
data_phthalates_sg[, phthalates_total_sg] <-                                    # on log transforme
  log(data_phthalates_sg[, phthalates_total_sg]) 
colnames(data_phthalates_sg) <-                                                 # on ajoute "_ln" à chaque fin de nom de variable
  c("ident", paste(colnames(data_phthalates_sg[, 2:55]), "ln", sep="_")) 
bdd <- left_join(bdd,                                                           # on merge avec la bdd principale
                 data_phthalates_sg, 
                 by ="ident")
rm(data_phthalates_sg)

# bdd_sg <- read_sas(
#   "0_source_data/old/metadata_220901/base_aline_220901.sas7bdat",                  # base de données SEPAGES
#   catalog_file = "0_source_data/formats.sas7bcat") %>%
#   select(ident, 
#          mo_pool_sg_T1, mo_pool_sg_T3, ch_pool_sg_M2, ch_pool_sg_Y1) 

# # Nettoyage gravité spécifique xx_pool_sg_xx 
# # A actualiser dans les codes de nettoyage de données au moment du changement de projet
# sg_vec <- bdd_sg %>% select(mo_pool_sg_T1, mo_pool_sg_T3, ch_pool_sg_M2, ch_pool_sg_Y1) %>% colnames()
# sg_vec_2 <- sg_vec %>% paste("2", sep = "_")
# sg_vec_ter <- sg_vec_2 %>% paste("ter", sep = "_")
# bdd_sg[, sg_vec] <- lapply(bdd_sg[, sg_vec], as.numeric)
# 
# bdd_sg <- bdd_sg %>% 
#   mutate(ch_pool_sg_M2_2 = ifelse(ch_pool_sg_M2 >900, NA, ch_pool_sg_M2)) %>%   # Exclusion de ident outlier 
#   
#   mutate(mo_pool_sg_T1_2 = mo_pool_sg_T1*1000,                                  # Changement de l'unité de la gravité spécifique 
#          mo_pool_sg_T3_2 = mo_pool_sg_T3*1000, 
#          ch_pool_sg_M2_2 = ch_pool_sg_M2_2*1000, 
#          ch_pool_sg_Y1_2 = ch_pool_sg_Y1*1000)
# 
# bdd_sg_ter <- bdd_sg %>%                          # Création variables sg sous formes de tertiles 
#   select(ident, all_of(sg_vec_2))%>% 
#   select(ident, everything()) 
# bdd_sg_ter[, sg_vec_2] <- lapply(
#   bdd_sg_ter[, sg_vec_2],
#   quant.cut,
#   nbclass = 3,
#   include.lowest = TRUE,
#   right = FALSE,
#   dig.lab = 2)
# 
# colnames(bdd_sg_ter) <- c("ident", paste(colnames(bdd_sg_ter[, 2:5]), "ter", sep="_"))
# bdd_sg_ter <- bdd_sg_ter %>%
#   mutate(
#     mo_pool_sg_T1_2_ter = fct_recode(mo_pool_sg_T1_2_ter,
#                                      "1st tertile" = "[1006,1016)",
#                                      "2nd tertile" = "[1016,1020)",
#                                      "3rd tertile" = "[1020,1035]"), 
#     mo_pool_sg_T3_2_ter = fct_recode(mo_pool_sg_T3_2_ter,
#                                      "1st tertile" = "[1e+03,1.01e+03)",
#                                      "2nd tertile" = "[1.01e+03,1.02e+03)",
#                                      "3rd tertile" = "[1.02e+03,1.04e+03]"), 
#     ch_pool_sg_M2_2_ter = fct_recode(ch_pool_sg_M2_2_ter,
#                                      "1st tertile" = "[1000,1004)",
#                                      "2nd tertile" = "[1004,1005)",
#                                      "3rd tertile" = "[1005,1011]"), 
#     ch_pool_sg_Y1_2_ter = fct_recode(ch_pool_sg_Y1_2_ter,
#                                      "1st tertile" = "[1003,1010)",
#                                      "2nd tertile" = "[1010,1013)",
#                                      "3rd tertile" = "[1013,1029]"))
# bdd_sg <- left_join(bdd_sg, 
#                     bdd_sg_ter, 
#                     by = "ident")
# metadata <- left_join(metadata, 
#                       bdd_sg, 
#                       by = "ident")


# 3_data cleaning outcomes ----
data_genera <- bdd %>%
  select(ident, 
         contains("ch_feces_rel_g"))
genera_total <- data_genera %>%
  select(contains("ch_feces_rel_g")) %>%
  filter(!is.na(ch_feces_rel_g1_Y1)) %>%
  select_if(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.3) %>%
  colnames()

data_genera <- data_genera %>% 
  column_to_rownames("ident")

var_label(data_genera) <-  str_replace(
  var_label(data_genera),                                    # set correct variable names                   
  "One year child feces relative abundance of ", "")
colnames(data_genera) <- var_label(data_genera)

genera <- data_genera %>%
  filter(!is.na(`genus Bifidobacterium`)) %>%
  select_if(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.3) %>%
  colnames()

descrip_genera_linear <- tbl_merge(
  tbls = 
    list(
      tbl_1 = 
        data_genera %>% 
        filter(!is.na(`genus Bifidobacterium`)) %>%
        select_if(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.3) %>%
        tbl_summary(
          type = list(everything() ~ "continuous"), 
          statistic = list(everything() ~ "{median} ({p25}, {p75})"), 
          digits = list(all_continuous() ~ c(2, 1, 1))), 
      tbl_2 = 
        data_genera %>% 
        filter(!is.na(`genus Bifidobacterium`)) %>%
        select_if(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.3) %>%
        mutate_all(~ ifelse(.>0, "Yes", "No")) %>%
        set_label(genera) %>%
        tbl_summary(
          type = list(everything() ~ "categorical"))), 
  tab_spanner = c("**Continuous**", "**Categorical (Y/N)**"))

data_genera_log <- data_genera %>% 
  select(all_of(genera)) %>%
  mutate_all(., ~ ifelse(. == 0, 1/5000, .)) %>%        # remplacement des valeurs 0 par 1/5000
  mutate_all(~ log(.)) %>%                              # transformation logarithmique
  rename_with(~gsub("genus ", "", .), everything()) %>% # changement des noms de colonnes pour qu'ils n'aient pas d'espace
  rownames_to_column(var = "ident")

genera <- str_replace_all(genera, "genus ", "")

data_genera_log$ident <- as.numeric(data_genera_log$ident)
bdd <- left_join(bdd, data_genera_log, by = "ident")

rm(data_genera, data_genera_log, descrip_genera_linear)

save.image("1_intermediate_data/1_data_cleaning_AD_gumme.RData")
