# A. Davias
# 08/02/2023

## Chargement des données et des fonctions existantes ----
load("2_final_data/metadata.RData")
load("2_final_data/bdd_alpha.RData")
load("2_final_data/bdd_taxa.RData")
source("3_programs/4_functions_AD_gumme.R", encoding = 'UTF-8')      # fonctions
source("3_programs/4_vectors_AD_gumme.R", echo=TRUE)
library(psych)
library(writexl)


## Création de fonctions ----

# Fonction pour obtenir les résultats bruts
lm_func <- function(outcome, exposure, data){   
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
                bf_duration_till48w_4cat_i,
              data = data)
  return(model)
}

# Fonction pour le nettoyage de résultats bruts
data_prep <- function(results_list, outcome_name) {
  results_list <- map(results_list, broom::tidy, conf.int = TRUE) %>%
    bind_rows(.id = "model_id") %>%
    filter(!term %in%
             c("(Intercept)",
               "ch_feces_RUN_Y1R3",
               "ch_feces_age_w_Y1_i",
               "po_delmodVaginal delivery",
               "ch_food_intro_Y1_3cat_iBetween 6 and 12 months old",
               "ch_food_intro_Y1_3cat_iNot introduced at 12 months old",
               "ch_antibio_Y1_2cat_iYes",
               "mo_par_2cat1 child or more",
               "mo_pets_iOne or more",
               "ch_sexMale",
               "mo_tob_gr_anyt_yn_n2_iYes",
               "Mo_ETS_anyT_yn1_opt_iYes",
               "ch_ETS_12m_opt36mYes",
               "mo_interpreg_3cat2 years and more",
               "mo_interpreg_3catPrimiparous",
               "mo_dipl_3cat_i3-4years after graduation",
               "mo_dipl_3cat_i>=5years after graduation",
               "po_w_kg_3cat3-3.4 Kg",
               "po_w_kg_3cat>= 3.5 Kg",
               "po_he_3cat_i50-51 cm",
               "po_he_3cat_i>= 52 cm",
               "ch_w_Y1_3cat_i8.5-9.9 Kg",
               "ch_w_Y1_3cat_i>=10 Kg",
               "ch_he_Y1_3cat_i75-77.9 cm",
               "ch_he_Y1_3cat_i>=78 cm",
               "po_gd",
               "mo_age",
               "mo_bmi_bepr_3cat_i19-23.9 Kg/m2",
               "mo_bmi_bepr_3cat_i>=24 Kg/m2",
               "bf_duration_till48w_4cat_i<24 weeks",
               "bf_duration_till48w_4cat_i24-47 weeks",
               "bf_duration_till48w_4cat_iStill breastfeed at 48 weeks"))%>%
    mutate(
      model_type = as.factor("adjusted"),
      outcome_name = outcome_name,
      term = str_replace_all(term, "exposure", "")) %>% 
    rename(outcome = outcome_name) %>%
    rename(exposure = model_id) %>%
    select(model_type,
           outcome,
           exposure,
           term,
           everything())
  
  return(results_list)
}


# Fonction pour tableaux de résultats article
model <- function(data, 
                  outcome, 
                  exposure_vec, 
                  digit_beta_IC) {
  data %>%                  
    select(
      {{outcome}},
      all_of({{exposure_vec}}),
      all_of(covar_vec_i)) %>% 
    
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

effectif_column <- function(data, outcome, exposure_vec) {
  data %>%                   
    filter(!is.na({{outcome}}))%>%
    select(all_of({{exposure_vec}})) %>%
    tbl_summary(missing = "no", 
                statistic = all_continuous() ~ "{N_nonmiss}") %>%
    bold_labels()
}

# Fonction pour figures forestplots à partir du tableau bruts des résultats 
forestplot <- function(results_list, outcome_name) {
  results <- results_list %>%
    ggplot(aes(x = exposure, 
               y = estimate, 
               min = conf.low, 
               ymax = conf.high, 
               #color = interaction(exposure_window, term_2), 
               color = term_rec, 
               shape = p_value_shape)) +
    geom_hline(yintercept = 0, linetype="dashed") +
    geom_pointrange(position = position_dodge(width = 0.7), size = 0.3) +
    labs(x = "Exposures", y = outcome_name) +
    theme_bw() +
    coord_flip()  +
    scale_shape_manual(values = c(19, 21),
                       name = "p-value") +
    guides(color = guide_legend(title = ""))+
    theme(axis.title = element_text(size = 7),
          axis.text = element_text(size = 6),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7), 
          legend.position = "bottom",
          legend.box = "vertical", 
          legend.justification = "right", 
          legend.spacing.y = unit(0, "cm"), 
          legend.spacing.x = unit(0, "cm"), 
          legend.box.margin = margin(0,0,0,0, "cm"), 
          legend.margin = margin(0,0,0,0, "cm"))
  
  return(results)
}


# Fonction pour voir la table de correlation avant de faire la correction pour compararaison multiple selon correlation entre les tests
cor_mixed <- function(data,method="spearman"){
  ## Purpose: 2 -  use the function mixedCor from psych to calculate a 
  #                correlation matrix on mixed type variables (continuous and categorical).
  #                Note: mixedCor consider categorical variable as ordered factors.
  ## Inputs: - data: dataframe (not mice) with only the variables to consider 
  ##                 (mixed type allowed)
  ## Output: - correlation matrix
  
  ## STEP 1 : TYPE OF VARIABLES
  continuous_var = which(sapply(data, class) == "numeric")
  names(continuous_var)=NULL
  # Categorical var
  categorical_var = which(sapply(data, class) == "factor")
  #  - 2 levels only
  binary_var = categorical_var[sapply(data[,categorical_var],nlevels)==2] 
  binary_var = binary_var[!is.na(binary_var)]
  names(binary_var)=NULL
  #  - More than 2 levels (but less than 8)
  poly_var = categorical_var[(sapply(data[,categorical_var],nlevels)>2 & sapply(data[,categorical_var],nlevels)<8)] %>% na.exclude()
  names(poly_var)=NULL
  
  ## STEP 2 : CORRELATION MATRIX USING MIXEDCOR FUNCTION (FROM PSYCH)
  # data converted in numeric (necessary)
  data[,] = lapply(data[,],as.numeric)
  # Correlation matrix
  cor = data  %>% 
    mixedCor(c=continuous_var,p=poly_var,d=binary_var,use="pairwise.complete.obs",method=method)%>% pluck('rho')
  return(cor)
}

# Fonction pour voir l'alpha corrigé pour la correction pour comparaison multiple 
alpha_corrected <- function(data, alpha=0.05) {
  # Purpose: Alpha-risk correction (FWER), inspired from https://www-ncbi-nlm-nih-gov.gate2.inist.fr/pmc/articles/PMC3325408/
  # Inputs
  # - data = dataset with all exposures to consider
  # - alpha = risk to correct
  # Output
  # - alpha corrected
  
  M <- ncol(data)
  
  ## FIRST STEP: CORRELATION MATRIX (see function "cor_mixed" in the path "general_functions.R")
  cor = cor_mixed(data) #handle mixed typed variables
  
  ## SECOND STEP : EIGENVALUES OF CORRELATION MATRIX
  lambdas <- base::eigen(cor)$values
  
  ## THIRD STEP: CORRECT ALPHA
  M0 <- M - sum( (lambdas > 1) * (lambdas - 1))
  M0 <- floor(M0) # always round down estimated value
  alpha_corrected <- alpha/M0
  
  return(alpha_corrected)
  
}

# Fonction pour voir le nombre de tests rééls après prise en compte de la corrélation entre les expositions 
# Peut aussi s'appliquer à la corrélation entre les outcomes si besoin 
# puis faire nombre d'expo x nombre d'outcome
# multiplier ce nombre à la p-value
# = cette correction pour comparaison multiple adapté de Li et al 2012 est une correction de Bonferroni qui prend en compte la corrélation entre les tests effectués
M0_corrected <- function(data, alpha=0.05) {
  # Purpose: Alpha-risk correction (FWER), inspired from https://www-ncbi-nlm-nih-gov.gate2.inist.fr/pmc/articles/PMC3325408/
  # Inputs
  # - data = dataset with all exposures to consider
  # - alpha = risk to correct
  # Output
  # - alpha corrected
  
  M <- ncol(data)
  
  ## FIRST STEP: CORRELATION MATRIX (see function "cor_mixed" in the path "general_functions.R")
  cor = cor_mixed(data, method = "pearson") #handle mixed typed variables
  
  ## SECOND STEP : EIGENVALUES OF CORRELATION MATRIX
  lambdas <- base::eigen(cor)$values
  
  ## THIRD STEP: CORRECT ALPHA
  M0 <- M - sum( (lambdas > 1) * (lambdas - 1))
  M0 <- floor(M0) # always round down estimated value
  alpha_corrected <- alpha/M0
  
  return(M0)
  
}



## Création des vecteurs ----
conf_vec <- bdd_alpha %>% 
  select(all_of(phthalates_vec))%>% 
  select(contains(c("DEHP", "MnBP"))) %>% 
  colnames()
explo_vec <- bdd_alpha %>% 
  select(all_of(phthalates_vec)) %>% 
  select(!contains(c("MEOHP", "MECPP", "MEHHP", "MEHP", "MMCHP",    # on ne met pas en analyse principale les métabolite quand on peut étudier la masse molaire
                     "DEHP",                                        # analyse confirmatoire
                     "MnBP",                                        # analyse confirmatoire
                     "ohMiNP", "oxoMiNP", "cxMiNP",                 # on ne met pas en analyse principale les métabolite quand on peut étudier la masse molaire
                     "ohMINCH", "oxoMINCH"))) %>%                   # on ne met pas en analyse principale les métabolite quand on peut étudier la masse molaire
  colnames()
pollutants_vec <- bdd_alpha %>% 
  select(all_of(phthalates_vec)) %>% 
  select(!contains(c("MEOHP", "MECPP", "MEHHP", "MEHP", "MMCHP", "ohMiNP", "oxoMiNP", "cxMiNP", "ohMINCH", "oxoMINCH"))) %>% 
  select(!contains("M2")) %>%
  colnames()
alpha_vec <- bdd_alpha %>% 
  select(all_of(alpha_vec)) %>% 
  select(contains("5000")) %>% 
  colnames() 



## Results - Unipolluants analysis ----
### Tableaux bruts ----
#### Colonnes principales ----
results_list_1 <- lapply(bdd_alpha[, pollutants_vec], lm_func, outcome = bdd_alpha$ch_feces_SpecRich_5000_ASV_Y1, data = bdd_alpha)
results_list_2 <- lapply(bdd_alpha[, pollutants_vec], lm_func, outcome = bdd_alpha$ch_feces_Shannon_5000_ASV_Y1, data = bdd_alpha)
results_list_3 <- lapply(bdd_alpha[, pollutants_vec], lm_func, outcome = bdd_alpha$ch_feces_Faith_5000_ASV_Y1, data = bdd_alpha)

results_list_4 <- lapply(bdd_taxa[, pollutants_vec], lm_func, outcome = bdd_taxa$ch_feces_rel_p1_Y1, data = bdd_taxa)
results_list_5 <- lapply(bdd_taxa[, pollutants_vec], lm_func, outcome = bdd_taxa$ch_feces_rel_p2_Y1, data = bdd_taxa)
results_list_6 <- lapply(bdd_taxa[, pollutants_vec], lm_func, outcome = bdd_taxa$ch_feces_rel_p3_Y1, data = bdd_taxa)
results_list_7 <- lapply(bdd_taxa[, pollutants_vec], lm_func, outcome = bdd_taxa$ch_feces_rel_p4_Y1, data = bdd_taxa)

results_list_8 <- lapply(bdd_taxa[, pollutants_vec], lm_func, outcome = bdd_taxa$ch_feces_rel_g1_Y1, data = bdd_taxa)
results_list_9 <- lapply(bdd_taxa[, pollutants_vec], lm_func, outcome = bdd_taxa$ch_feces_rel_g2_Y1, data = bdd_taxa)
results_list_10 <- lapply(bdd_taxa[, pollutants_vec], lm_func, outcome = bdd_taxa$ch_feces_rel_g3_Y1, data = bdd_taxa)
results_list_11 <- lapply(bdd_taxa[, pollutants_vec], lm_func, outcome = bdd_taxa$ch_feces_rel_g4_Y1, data = bdd_taxa)

results_multi <- list(
  
  data_prep(results_list = results_list_1, outcome_name = "Specific richness"),
  data_prep(results_list = results_list_2, outcome_name = "Shannon diversity"),
  data_prep(results_list = results_list_3, outcome_name = "Faith phylogenetic diversity"),
  
  data_prep(results_list = results_list_4, outcome_name = "Firmicutes"),
  data_prep(results_list = results_list_5, outcome_name = "Actinobacteria"),
  data_prep(results_list = results_list_6, outcome_name = "Bacteroidetes"),
  data_prep(results_list = results_list_7, outcome_name = "Proteobacteria"),
  
  data_prep(results_list = results_list_8, outcome_name = "Bifidobacterium"),
  data_prep(results_list = results_list_9, outcome_name = "Bacteroides"),
  data_prep(results_list = results_list_10, outcome_name = "Blautia"),
  data_prep(results_list = results_list_11, outcome_name = "Escherichia and Shigella")) %>%
  
  bind_rows() %>%
  mutate(
    exposure = as.factor(exposure)) %>%
  mutate(
    # exposure_type = case_when(
    #   exposure == "ch_ohMPHP_cat_M2_2" ~ "categorical_2",
    #   exposure == "ch_ohMINCH_cat_M2_2" ~ "categorical_2",
    #   exposure == "ch_oxoMINCH_cat_M2_2" ~ "categorical_2",
    #   TRUE ~ "continuous"), 
    exposure = str_replace_all(exposure,
                               c("mo_" = "",
                                 "ch_" = "",
                                 "_total_i_cor_" = " ",
                                 "_i_cor_" = " ", 
                                 "_ln" = "",
                                 "_conj" = "",
                                 "ln" = "",
                                 "_cor" = "", 
                                 "_t2" = " t2", 
                                 "_t3" = " t3", 
                                 #"_M2" = " M2", 
                                 "_Y1" = " Y1", 
                                 "_ms" = ""
                                 #, 
                                 #"_cat M2_2" = " M2"
                                 )), 

    p_value_shape = case_when(p.value < 0.1 ~ "p.value <0.1",
                              p.value > 0.1 ~ "p.value >0.1"), 
    p_value_shape = fct_relevel(p_value_shape,
                                "p.value >0.1", 
                                "p.value <0.1"),
    exposure_window = case_when(grepl("t2", exposure) ~ "Trim.2", 
                                grepl("t3", exposure) ~ "Trim.3", 
                                #grepl("M2", exposure) ~ "2 months", 
                                grepl("Y1", exposure) ~ "12 months", 
                                TRUE ~ "Trim.2"), 
    exposure_window = fct_relevel(exposure_window, 
                                  "Trim.2", 
                                  "Trim.3", 
                                  #"2 months", 
                                  "12 months"),
    analysis = case_when(grepl("DEHP", exposure) ~ "confirmatory", 
                         grepl("MnBP", exposure) ~ "confirmatory", 
                         TRUE ~ "exploratory")
    # , 
    # term_rec = case_when(exposure_type == "continuous" ~ "Continuous", 
    #                      exposure_type == "categorical_2" & term == ">LOD" ~ ">LOD, compared to <LOD"),
    # term_rec = fct_relevel(term_rec,
    #                        "Continuous", ">LOD, compared to <LOD")
    ) %>%
  select("model_type", 
         "analysis",
         "outcome", 
         "exposure",
         #"exposure_type",
         "exposure_window",
         "term",
         #"term_rec",
         "estimate",
         "std.error",
         "statistic",
         "conf.low",
         "conf.high",
         "p.value", 
         "p_value_shape")

#### Tableau pour Marion ---- 
bdd_marion <- results_multi %>%
  select(outcome, 
         Pollutant = exposure, 
         Variable_type = exposure_type, 
         #Variable_categories = term_rec, 
         estimate, 
         p.value) %>%
  separate(Pollutant, c("Pollutant", "Timing")) %>%
  mutate(
    Timing = fct_recode(Timing, "T2" = "t2", "T3" = "t3"), 
    Variable_type = fct_recode(Variable_type, "cat_2" = "categorical_2"),
    ch_feces_specrich_5000_Y1_b = ifelse(outcome == "Specific richness", estimate, NA), 
    ch_feces_specrich_5000_Y1_p = ifelse(outcome == "Specific richness", p.value, NA), 
    ch_feces_shannon_5000_Y1_b = ifelse(outcome == "Shannon diversity", estimate, NA), 
    ch_feces_shannon_5000_Y1_p = ifelse(outcome == "Shannon diversity", p.value, NA), 
    ch_feces_faith_5000_Y1_b = ifelse(outcome == "Faith phylogenetic diversity", estimate, NA), 
    ch_feces_faith_5000_Y1_p = ifelse(outcome == "Faith phylogenetic diversity", p.value, NA), 
    
    ch_feces_firmicutes_5000_Y1_b = ifelse(outcome == "Firmicutes", estimate, NA), 
    ch_feces_firmicutes_5000_Y1_p = ifelse(outcome == "Firmicutes", p.value, NA), 
    ch_feces_actinobacteria_5000_Y1_b = ifelse(outcome == "Actinobacteria", estimate, NA), 
    ch_feces_actinobacteria_5000_Y1_p = ifelse(outcome == "Actinobacteria", p.value, NA),
    ch_feces_bacteroidetes_5000_Y1_b = ifelse(outcome == "Bacteroidetes", estimate, NA), 
    ch_feces_bacteroidetes_5000_Y1_p = ifelse(outcome == "Bacteroidetes", p.value, NA),
    ch_feces_proteobacteria_5000_Y1_b = ifelse(outcome == "Proteobacteria", estimate, NA), 
    ch_feces_proteobacteria_5000_Y1_p = ifelse(outcome == "Proteobacteria", p.value, NA),
    
    ch_feces_blautia_5000_Y1_b = ifelse(outcome == "Blautia", estimate, NA), 
    ch_feces_blautia_5000_Y1_p = ifelse(outcome == "Blautia", p.value, NA),
    ch_feces_bifidobacterium_5000_Y1_b = ifelse(outcome == "Bifidobacterium", estimate, NA), 
    ch_feces_bifidobacterium_5000_Y1_p = ifelse(outcome == "Bifidobacterium", p.value, NA),
    ch_feces_bacteroides_5000_Y1_b = ifelse(outcome == "Bacteroides", estimate, NA), 
    ch_feces_bacteroides_5000_Y1_p = ifelse(outcome == "Bacteroides", p.value, NA),
    ch_feces_escherichia_shigella_5000_Y1_b = ifelse(outcome == "Escherichia and Shigella", estimate, NA), 
    ch_feces_escherichia_shigella_5000_Y1_p = ifelse(outcome == "Escherichia and Shigella", p.value, NA)) %>%
  select(-outcome, -estimate, -p.value)

write_xlsx(
  bdd_marion, 
  "C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/7. Présentations écrites/9. Figure_overview_SEPAGES/SEPAGES_overview_phthalates.xlsx")


#### Correction pour comparaison multiple (tous les tests pris en compte) ----
# test_expo <- bdd_alpha %>% 
#   filter(!is.na(ch_feces_Shannon_5000_ASV_Y1)) %>%
#   select(all_of(pollutants_vec)) %>% na.omit()
# 
# test_outcome <- bdd_alpha %>% 
#   filter(!is.na(ch_feces_Shannon_5000_ASV_Y1)) %>%
#   select(ident, 
#          ch_feces_Shannon_5000_ASV_Y1, 
#          ch_feces_Faith_5000_ASV_Y1, 
#          ch_feces_SpecRich_5000_ASV_Y1) %>% 
#   na.omit() %>%
#   select(-ident)
# 
# # test_outcome_2 <- bdd_taxa %>% 
# #   filter(!is.na(ch_feces_rel_p1_Y1)) %>%
# #   select(ident, 
# #          ch_feces_rel_p1_Y1:ch_feces_rel_p4_Y1,
# #          ch_feces_rel_g1_Y1:ch_feces_rel_g4_Y1) %>% na.omit()
# # 
# # test_outcome <- left_join(test_outcome, test_outcome_2, by = "ident") %>% select(-ident)
# 
# var_lab(test_outcome$ch_feces_Shannon_5000_ASV_Y1) <- NULL 
# var_lab(test_outcome$ch_feces_Faith_5000_ASV_Y1) <- NULL 
# var_lab(test_outcome$ch_feces_SpecRich_5000_ASV_Y1) <- NULL 
# 
# var_lab(test_outcome$ch_feces_rel_p1_Y1) <- NULL 
# var_lab(test_outcome$ch_feces_rel_p2_Y1) <- NULL 
# var_lab(test_outcome$ch_feces_rel_p3_Y1) <- NULL 
# var_lab(test_outcome$ch_feces_rel_p4_Y1) <- NULL 
# var_lab(test_outcome$ch_feces_rel_g1_Y1) <- NULL 
# var_lab(test_outcome$ch_feces_rel_g2_Y1) <- NULL 
# var_lab(test_outcome$ch_feces_rel_g3_Y1) <- NULL 
# var_lab(test_outcome$ch_feces_rel_g4_Y1) <- NULL 
# 
# cor_mixed_table_expo <- cor_mixed(data = test_expo, method = "pearson")
# results_alpha_corrected_expo <- alpha_corrected(data = test_expo, alpha = 0.05)
# results_M0_corrected_expo <- M0_corrected(data = test_expo, alpha = 0.05)
# results_M0_corrected_expo
# ## on passe de 24 expositions à 14 expositions après prise en compte de leur corrélation
# 
# cor_mixed_table_outcome <- cor_mixed(data = test_outcome, method = "pearson")
# results_alpha_corrected_outcome <- alpha_corrected(data = test_outcome, alpha = 0.05)
# results_M0_corrected_outcome <- M0_corrected(data = test_outcome, alpha = 0.05)
# results_M0_corrected_outcome
# ## on passe de 3 outcomes à 2 outcomes après prise en compte de leur corrélation
# 
# # 
# # results_multi <- results_multi %>%
# #   mutate(
# #     FWER.p.value = p.value * 17 * 6,
# #     FWER.p.value = ifelse(FWER.p.value > 1, 1, FWER.p.value)) 
# # 
# # results_multi <- results_multi %>%
# #   mutate(
# #     FWER.p.value_2 = round(FWER.p.value, 2))


#### Correction pour comparaison multiple (seulement les tests alpha diversité pris en compte) ----
test_expo <- bdd_alpha %>% 
  filter(!is.na(ch_feces_Shannon_5000_ASV_Y1)) %>%
  select(all_of(pollutants_vec)) %>% na.omit()

test_outcome <- bdd_alpha %>% 
  filter(!is.na(ch_feces_Shannon_5000_ASV_Y1)) %>%
  select(ch_feces_Shannon_5000_ASV_Y1, 
         ch_feces_Faith_5000_ASV_Y1, 
         ch_feces_SpecRich_5000_ASV_Y1) %>% 
  na.omit() 

var_lab(test_outcome$ch_feces_Shannon_5000_ASV_Y1) <- NULL 
var_lab(test_outcome$ch_feces_Faith_5000_ASV_Y1) <- NULL 
var_lab(test_outcome$ch_feces_SpecRich_5000_ASV_Y1) <- NULL 

cor_mixed_table_expo <- cor_mixed(data = test_expo, method = "pearson")
results_alpha_corrected_expo <- alpha_corrected(data = test_expo, alpha = 0.05)
results_M0_corrected_expo <- M0_corrected(data = test_expo, alpha = 0.05)
results_M0_corrected_expo
## on passe de 24 expositions à 14 expositions après prise en compte de leur corrélation

cor_mixed_table_outcome <- cor_mixed(data = test_outcome, method = "pearson")
results_alpha_corrected_outcome <- alpha_corrected(data = test_outcome, alpha = 0.05)
results_M0_corrected_outcome <- M0_corrected(data = test_outcome, alpha = 0.05)
results_M0_corrected_outcome
## on passe de 3 outcomes à 2 outcomes après prise en compte de leur corrélation

results_multi <- results_multi %>%
  mutate(
    FWER.p.value_alpha = ifelse(outcome %in% c("Specific richness", 
                                         "Shannon diversity", 
                                         "Faith phylogenetic diversity"), 
                          p.value * 14 * 2, NA),
    FWER.p.value_alpha = ifelse(FWER.p.value_alpha > 1, ">0.99", FWER.p.value_alpha)) 

results_multi <- results_multi %>%
  mutate(
    FWER.p.value_shape_alpha = case_when(p.value< 0.0012 & 
                                     outcome %in% c("Specific richness", 
                                                    "Shannon diversity", 
                                                    "Faith phylogenetic diversity")~ "p.value <0.0012",
                                   p.value > 0.0012 & 
                                     outcome %in% c("Specific richness", 
                                                    "Shannon diversity", 
                                                    "Faith phylogenetic diversity")~ "p.value >0.0012"), 
    FWER.p.value_shape_alpha = fct_relevel(FWER.p.value_shape_alpha,
                                "p.value >0.0012", 
                                "p.value <0.0012"))


#### Correction pour comparaison multiple (taxonomie pour les polluants associés à l'alpha div) ----
test_expo <- bdd_alpha %>% 
  filter(!is.na(ch_feces_Shannon_5000_ASV_Y1)) %>%
  select(ch_DEHP_ms_i_cor_Y1_ln, ch_MEP_i_cor_Y1_ln, ch_ohMPHP_i_cor_Y1_ln) %>% na.omit()

test_outcome <- bdd_taxa %>% 
  filter(!is.na(ch_feces_rel_p1_Y1)) %>%
  select(ch_feces_rel_p1_Y1, ch_feces_rel_p2_Y1, ch_feces_rel_p3_Y1, ch_feces_rel_p4_Y1,
         ch_feces_rel_g1_Y1, ch_feces_rel_g2_Y1, ch_feces_rel_g3_Y1, ch_feces_rel_g4_Y1) %>% 
  na.omit() 

var_lab(test_outcome$ch_feces_rel_p1_Y1) <- NULL 
var_lab(test_outcome$ch_feces_rel_p2_Y1) <- NULL 
var_lab(test_outcome$ch_feces_rel_p3_Y1) <- NULL 
var_lab(test_outcome$ch_feces_rel_p4_Y1) <- NULL 
var_lab(test_outcome$ch_feces_rel_g1_Y1) <- NULL 
var_lab(test_outcome$ch_feces_rel_g2_Y1) <- NULL 
var_lab(test_outcome$ch_feces_rel_g3_Y1) <- NULL 
var_lab(test_outcome$ch_feces_rel_g4_Y1) <- NULL 

cor_mixed_table_expo <- cor_mixed(data = test_expo, method = "pearson")
results_alpha_corrected_expo <- alpha_corrected(data = test_expo, alpha = 0.05)
results_M0_corrected_expo <- M0_corrected(data = test_expo, alpha = 0.05)
results_M0_corrected_expo
## on passe de 3 expositions à 2 expositions après prise en compte de leur corrélation

cor_mixed_table_outcome <- cor_mixed(data = test_outcome, method = "pearson")
results_alpha_corrected_outcome <- alpha_corrected(data = test_outcome, alpha = 0.05)
results_M0_corrected_outcome <- M0_corrected(data = test_outcome, alpha = 0.05)
results_M0_corrected_outcome
## on passe de 8 outcomes à 4 outcomes après prise en compte de leur corrélation

results_multi <- results_multi %>%
  mutate(
    FWER.p.value_taxa = ifelse(exposure %in% c("DEHP Y1", 
                                            "MEP Y1", 
                                            "ohMPHP Y1") & 
                            !outcome %in% c("Specific richness", 
                                           "Shannon diversity", 
                                           "Faith phylogenetic diversity"), 
                          p.value * 4 * 2,NA), 
    FWER.p.value_taxa = ifelse(FWER.p.value_taxa > 1, ">0.99", FWER.p.value_taxa)) 

results_multi <- results_multi %>%
  mutate(
    FWER.p.value_shape_taxa = 
      case_when(p.value< 0.006 & exposure %in% c("DEHP Y1", "MEP Y1", "ohMPHP Y1") & 
                  !outcome %in% c("Specific richness", "Shannon diversity", "Faith phylogenetic diversity") ~ 
                  "p<0.006 (pollutants associated with α-div and corrected for multiple testing)",
                p.value> 0.006 & exposure %in% c("DEHP Y1", "MEP Y1", "ohMPHP Y1") & 
                  !outcome %in% c("Specific richness", "Shannon diversity", "Faith phylogenetic diversity") ~ 
                  "p>0.006 (pollutants associated with α-div and corrected for multiple testing)",
                
                p.value< 0.05 & !exposure %in% c("DEHP Y1", "MEP Y1", "ohMPHP Y1") & 
                  !outcome %in% c("Specific richness", "Shannon diversity", "Faith phylogenetic diversity") ~ 
                  "p<0.05 (pollutants not associated with α-div and not corrected for multiple testing)",
                p.value> 0.05 & !exposure %in% c("DEHP Y1", "MEP Y1", "ohMPHP Y1") & 
                  !outcome %in% c("Specific richness", "Shannon diversity", "Faith phylogenetic diversity") ~ 
                  "p>0.05 (pollutants not associated with α-div and not corrected for multiple testing)"))


### Tableaux article - alpha diversité ----
# analyses confirmatoires sans correction pour les polluants avec hypothèses à priori
results_multi_alpha_conf <- tbl_merge(                         
  tbls = list(
    effectif_column(data = bdd_alpha, 
                    outcome = ch_feces_SpecRich_5000_ASV_Y1, 
                    exposure_vec = conf_vec),
    model(data = bdd_alpha, 
          outcome = ch_feces_SpecRich_5000_ASV_Y1,
          exposure_vec = conf_vec, 
          digit_beta_IC = 1),
    model(data = bdd_alpha,
          outcome = ch_feces_Shannon_5000_ASV_Y1,
          exposure_vec = conf_vec, 
          digit_beta_IC = 2),
    model(data = bdd_alpha, 
          outcome = ch_feces_Faith_5000_ASV_Y1,
          exposure_vec = conf_vec, 
          digit_beta_IC = 1)), 
  tab_spanner = c("", 
                  "**Specific richness**", 
                  "**Shannon diversity**", 
                  "**Faith's phylogenetic diversity**"))

# analyses exploratoires avec correction pour les polluants sans hypothèses à priori
results_multi_alpha_explo <- tbl_merge(
  tbls = list(
    effectif_column(
      data = bdd_alpha, 
      outcome = ch_feces_SpecRich_5000_ASV_Y1, 
      exposure_vec = explo_vec),
    model(
      data = bdd_alpha, 
      outcome = ch_feces_SpecRich_5000_ASV_Y1,
      exposure_vec = explo_vec, 
      digit_beta_IC = 1),
    model(
      data = bdd_alpha,
      outcome = ch_feces_Shannon_5000_ASV_Y1,
      exposure_vec = explo_vec, 
      digit_beta_IC = 2),
    model(
      data = bdd_alpha, 
      outcome = ch_feces_Faith_5000_ASV_Y1,
      exposure_vec = explo_vec, 
      digit_beta_IC = 1)), 
  tab_spanner = c("", 
                  "**Specific richness**", 
                  "**Shannon diversity**", 
                  "**Faith's phylogenetic diversity**"))


### Tableaux article - taxonomie ----
# Analyses "explicatives" non corrigées pour comparaison multiple pour la taxonomie 
# pour les associations significatives pour l'alpha diversité
# pour les autres polluants non significatifs --> on met tout en analyse sup. 
signi_vec <- c("ch_DEHP_ms_i_cor_Y1_ln", "ch_MnBP_i_cor_Y1_ln", "ch_MEP_i_cor_Y1_ln", "ch_ohMPHP_i_cor_Y1_ln")
not_signi_vec <- bdd_alpha %>% 
  select(all_of(conf_vec), all_of(explo_vec)) %>% 
  select(-c("ch_DEHP_ms_i_cor_Y1_ln", "ch_MnBP_i_cor_Y1_ln", "ch_MEP_i_cor_Y1_ln", "ch_ohMPHP_i_cor_Y1_ln")) %>%
  colnames()


results_multi_phyla_signi <- tbl_merge(
  tbls = list(
    effectif_column(
      data = bdd_taxa, 
      outcome = ch_feces_rel_p1_Y1, 
      exposure_vec = signi_vec),
    model(
      data = bdd_taxa, 
      outcome = ch_feces_rel_p1_Y1,
      exposure_vec = signi_vec, 
      digit_beta_IC = 1),
    model(
      data = bdd_taxa,
      outcome = ch_feces_rel_p2_Y1,
      exposure_vec = signi_vec, 
      digit_beta_IC = 1),
    model(
      data = bdd_taxa, 
      outcome = ch_feces_rel_p3_Y1,
      exposure_vec = signi_vec, 
      digit_beta_IC = 1), 
    model(
      data = bdd_taxa, 
      outcome = ch_feces_rel_p4_Y1,
      exposure_vec = signi_vec, 
      digit_beta_IC = 1)), 
  tab_spanner = c("", 
                  "**Phylum Firmicutes**",
                  "**Phylum Actinobacteria**",
                  "**Phylum Bacteroidetes**",
                  "**Phylum Proteobacteria**"))

results_multi_genera_signi <- tbl_merge(
  tbls = list(
    effectif_column(
      data = bdd_taxa, 
      outcome = ch_feces_rel_g1_Y1, 
      exposure_vec = signi_vec),
    model(
      data = bdd_taxa, 
      outcome = ch_feces_rel_g1_Y1,
      exposure_vec = signi_vec, 
      digit_beta_IC = 1),
    model(
      data = bdd_taxa,
      outcome = ch_feces_rel_g2_Y1,
      exposure_vec = signi_vec, 
      digit_beta_IC = 1),
    model(
      data = bdd_taxa, 
      outcome = ch_feces_rel_g3_Y1,
      exposure_vec = signi_vec, 
      digit_beta_IC = 1), 
    model(
      data = bdd_taxa, 
      outcome = ch_feces_rel_g4_Y1,
      exposure_vec = signi_vec, 
      digit_beta_IC = 1)), 
  tab_spanner = c("", 
                  "**Genus Bifidobacterium**",
                  "**Genus Bacteroides**",
                  "**Genus Blautia**",
                  "**Genera Eschericha / Shigella**"))


results_multi_phyla_not_signi <- tbl_merge(
  tbls = list(
    effectif_column(
      data = bdd_taxa, 
      outcome = ch_feces_rel_p1_Y1, 
      exposure_vec = not_signi_vec),
    model(
      data = bdd_taxa, 
      outcome = ch_feces_rel_p1_Y1,
      exposure_vec = not_signi_vec, 
      digit_beta_IC = 1),
    model(
      data = bdd_taxa,
      outcome = ch_feces_rel_p2_Y1,
      exposure_vec = not_signi_vec, 
      digit_beta_IC = 1),
    model(
      data = bdd_taxa, 
      outcome = ch_feces_rel_p3_Y1,
      exposure_vec = not_signi_vec, 
      digit_beta_IC = 1), 
    model(
      data = bdd_taxa, 
      outcome = ch_feces_rel_p4_Y1,
      exposure_vec = not_signi_vec, 
      digit_beta_IC = 1)), 
  tab_spanner = c("", 
                  "**Phylum Firmicutes**",
                  "**Phylum Actinobacteria**",
                  "**Phylum Bacteroidetes**",
                  "**Phylum Proteobacteria**"))

results_multi_genera_not_signi <- tbl_merge(
  tbls = list(
    effectif_column(
      data = bdd_taxa, 
      outcome = ch_feces_rel_g1_Y1, 
      exposure_vec = not_signi_vec),
    model(
      data = bdd_taxa, 
      outcome = ch_feces_rel_g1_Y1,
      exposure_vec = not_signi_vec, 
      digit_beta_IC = 1),
    model(
      data = bdd_taxa,
      outcome = ch_feces_rel_g2_Y1,
      exposure_vec = not_signi_vec, 
      digit_beta_IC = 1),
    model(
      data = bdd_taxa, 
      outcome = ch_feces_rel_g3_Y1,
      exposure_vec = not_signi_vec, 
      digit_beta_IC = 1), 
    model(
      data = bdd_taxa, 
      outcome = ch_feces_rel_g4_Y1,
      exposure_vec = not_signi_vec, 
      digit_beta_IC = 1)), 
  tab_spanner = c("", 
                  "**Genus Bifidobacterium**",
                  "**Genus Bacteroides**",
                  "**Genus Blautia**",
                  "**Genera Eschericha / Shigella**"))

### Tableaux article - investigation des genres firmicutes ----
# Analyses "explicatives" non corrigées pour comparaison multiple pour la taxonomie 
# pour les associations significatives pour l'alpha diversité
# pour les autres polluants non significatifs --> on met tout en analyse sup. 
taxa <- read.csv("~/5. R projects/pollutants_gut_microbiota_Y1/0_source_data/taxa_table_ASVbased_Y1_AD_20220504_8.csv")

taxa %>% filter(ch_feces_phylum_ASVbased_Y1 == "Firmicutes") %>% select(ch_feces_phylum_ASVbased_Y1, ch_feces_genus_ASVbased_Y1) %>% View()
bdd_taxa %>% select(ident, contains("ch_feces_rel_g")) %>% View()
descrip_num(data = bdd_taxa, vars = c("ch_feces_rel_g1_Y1", "ch_feces_rel_g2_Y1", "ch_feces_rel_g3_Y1", "ch_feces_rel_g4_Y1"))
densityplot(data = bdd_taxa, vars = c("ch_feces_rel_g1_Y1", "ch_feces_rel_g2_Y1", "ch_feces_rel_g3_Y1", "ch_feces_rel_g4_Y1"))
densityplot(data = bdd_taxa, vars = c("ch_feces_rel_p1_Y1", "ch_feces_rel_p2_Y1", "ch_feces_rel_p3_Y1", "ch_feces_rel_p4_Y1"))

descrip_num(data = bdd_taxa, 
            vars = c("ch_feces_rel_g6_Y1", "ch_feces_rel_g7_Y1", "ch_feces_rel_g8_Y1", "ch_feces_rel_g9_Y1", 
                     "ch_feces_rel_g10_Y1", "ch_feces_rel_g12_Y1", "ch_feces_rel_g13_Y1", "ch_feces_rel_g14_Y1", "ch_feces_rel_g15_Y1"))

bdd_taxa <- bdd_taxa %>%
  mutate(
    ch_feces_rel_g6_Y1_ln = ch_feces_rel_g6_Y1^2,
    ch_feces_rel_g7_Y1_ln = ch_feces_rel_g7_Y1^2, 
    ch_feces_rel_g8_Y1_ln = ch_feces_rel_g8_Y1^2, 
    ch_feces_rel_g9_Y1_ln = ch_feces_rel_g9_Y1^2,
    ch_feces_rel_g10_Y1_ln = ch_feces_rel_g10_Y1^2, 
    ch_feces_rel_g12_Y1_ln = ch_feces_rel_g12_Y1^2, 
    ch_feces_rel_g13_Y1_ln = ch_feces_rel_g13_Y1^2, 
    ch_feces_rel_g14_Y1_ln = ch_feces_rel_g14_Y1^2, 
    ch_feces_rel_g15_Y1_ln = ch_feces_rel_g15_Y1^2)

densityplot(data = bdd_taxa, 
            vars = c("ch_feces_rel_g6_Y1", "ch_feces_rel_g7_Y1", "ch_feces_rel_g8_Y1", "ch_feces_rel_g9_Y1", 
                     "ch_feces_rel_g10_Y1", "ch_feces_rel_g12_Y1", "ch_feces_rel_g13_Y1", "ch_feces_rel_g14_Y1", "ch_feces_rel_g15_Y1", 
                     "ch_feces_rel_g6_Y1_ln", "ch_feces_rel_g7_Y1_ln", "ch_feces_rel_g8_Y1_ln", "ch_feces_rel_g9_Y1_ln", 
                     "ch_feces_rel_g10_Y1_ln", "ch_feces_rel_g12_Y1_ln", "ch_feces_rel_g13_Y1_ln", "ch_feces_rel_g14_Y1_ln", "ch_feces_rel_g15_Y1_ln"))

firmi_vec <- c("ch_feces_rel_g6_Y1_ln", "ch_feces_rel_g7_Y1_ln", "ch_feces_rel_g8_Y1_ln", "ch_feces_rel_g9_Y1_ln", 
               "ch_feces_rel_g10_Y1_ln", "ch_feces_rel_g12_Y1_ln", "ch_feces_rel_g13_Y1_ln", "ch_feces_rel_g14_Y1_ln", "ch_feces_rel_g15_Y1_ln")


results_list_12 <- lapply(bdd_taxa[, pollutants_vec], lm_func, outcome = bdd_taxa$ch_feces_rel_g6_Y1_ln, data = bdd_taxa)
results_list_13 <- lapply(bdd_taxa[, pollutants_vec], lm_func, outcome = bdd_taxa$ch_feces_rel_g7_Y1_ln, data = bdd_taxa)
results_list_14 <- lapply(bdd_taxa[, pollutants_vec], lm_func, outcome = bdd_taxa$ch_feces_rel_g8_Y1_ln, data = bdd_taxa)
results_list_15 <- lapply(bdd_taxa[, pollutants_vec], lm_func, outcome = bdd_taxa$ch_feces_rel_g9_Y1_ln, data = bdd_taxa)
results_list_16 <- lapply(bdd_taxa[, pollutants_vec], lm_func, outcome = bdd_taxa$ch_feces_rel_g10_Y1_ln, data = bdd_taxa)
results_list_17 <- lapply(bdd_taxa[, pollutants_vec], lm_func, outcome = bdd_taxa$ch_feces_rel_g12_Y1_ln, data = bdd_taxa)
results_list_18 <- lapply(bdd_taxa[, pollutants_vec], lm_func, outcome = bdd_taxa$ch_feces_rel_g13_Y1_ln, data = bdd_taxa)
results_list_19 <- lapply(bdd_taxa[, pollutants_vec], lm_func, outcome = bdd_taxa$ch_feces_rel_g14_Y1_ln, data = bdd_taxa)
results_list_20 <- lapply(bdd_taxa[, pollutants_vec], lm_func, outcome = bdd_taxa$ch_feces_rel_g15_Y1_ln, data = bdd_taxa)

lm(ch_feces_rel_g6_Y1_ln ~
     mo_DiNP_ms_i_cor_t2_ln +
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
     bf_duration_till48w_4cat_i,
   data = bdd_taxa)




### Tableaux article - mixture effects ----
### Figures article - forestplots ----
#### Alpha diversity ----
results_multi <- results_multi %>%
  mutate(
    exposure =  fct_relevel(
      exposure, 
      "DINCH Y1", "DINCH t3", "DINCH t2", 
      "ohMPHP Y1", "ohMPHP t3", "ohMPHP t2", 
      "MEP Y1", "MEP t2", "MEP t3", 
      "MBzP Y1", "MBzP t3", "MBzP t2", 
      "MiBP Y1", "MiBP t3", "MiBP t2", 
      "DiNP Y1", "DiNP t3", "DiNP t2", 
      "MnBP Y1", "MnBP t3", "MnBP t2", 
      "DEHP Y1", "DEHP t3", "DEHP t2"), 
    exposure = fct_recode(
      exposure,
      "ΣDINCH Y1" = "DINCH Y1",
      "ΣDINCH t3" = "DINCH t3",
      "ΣDINCH t2" = "DINCH t2",
      "ΣDiNP Y1" = "DiNP Y1",
      "ΣDiNP t3" = "DiNP t3",
      "ΣDiNP t2" = "DiNP t2",
      "ΣDEHP Y1" = "DEHP Y1",
      "ΣDEHP t3" = "DEHP t3",
      "ΣDEHP t2" = "DEHP t2"))

# forestplot <- function(results_list, outcome_name) {
#   results <- results_list %>%
#     ggplot(aes(x = exposure, 
#                y = estimate, 
#                min = conf.low, 
#                ymax = conf.high, 
#                #color = interaction(exposure_window, term_2), 
#                #color = term_rec, 
#                shape = FWER.p.value_shape)) +
#     geom_hline(yintercept = 0, linetype="dashed") +
#     geom_pointrange(position = position_dodge(width = 0.7), size = 0.4) +
#     labs(x = "Exposures", y = outcome_name) +
#     theme_bw() +
#     coord_flip()  +
#     scale_shape_manual(values = c(19, 21),
#                        name = "p-value corrected") +
#     guides(color = "none")+
#     theme(axis.title = element_text(size = 7),
#           axis.text = element_text(size = 7),
#           legend.text = element_text(size = 7),
#           legend.title = element_text(size = 7), 
#           legend.position = "bottom",
#           legend.box = "vertical", 
#           legend.justification = "center", 
#           legend.spacing.y = unit(0, "cm"), 
#           legend.spacing.x = unit(0, "cm"), 
#           legend.box.margin = margin(0,0,0,0, "cm"), 
#           legend.margin = margin(0,0,0,0, "cm"))
#   
#   return(results)
# }



forestplot <- function(results_list, outcome_name) {
  results <- results_list %>%
    ggplot(aes(x = exposure,
               y = estimate,
               min = conf.low,
               ymax = conf.high,
               #color = interaction(exposure_window, term_2),
               #color = term_rec,
               color = FWER.p.value_shape)) +
    geom_hline(yintercept = 0, linetype="dashed") +
    geom_pointrange(position = position_dodge(width = 0.7), size = 0.4,
                    aes(color = ifelse(FWER.p.value_shape_alpha == "p.value >0.0012", "p.value >0.0012", "p.value <0.0012"))) +
    labs(x = "Exposures", y = outcome_name) +
    theme_bw() +
    coord_flip()  +
    scale_color_manual(values = c("black", "red"),
                       name = "") +
    #guides(color = "none")+
    theme(axis.title = element_text(size = 9),
          axis.text = element_text(size = 9),
          legend.text = element_text(size = 9),
          legend.title = element_text(size = 9),
          legend.position = "bottom",
          legend.box = "vertical",
          legend.justification = "center",
          legend.spacing.y = unit(0, "cm"),
          legend.spacing.x = unit(0, "cm"),
          legend.box.margin = margin(0,0,0,0, "cm"),
          legend.margin = margin(0,0,0,0, "cm"))

  return(results)
}

forestplot_shannon <- function(results_list, outcome_name) {
  results <- results_list %>%
    ggplot(aes(x = exposure,
               y = estimate,
               min = conf.low,
               ymax = conf.high,
               #color = interaction(exposure_window, term_2),
               #color = term_rec,
               color = FWER.p.value_shape)) +
    geom_hline(yintercept = 0, linetype="dashed") +
    geom_pointrange(position = position_dodge(width = 0.7), size = 0.4,
                    aes(color = ifelse(FWER.p.value_shape_alpha == "p.value >0.0012", "p.value >0.0012", "p.value <0.0012"))) +
    labs(x = "Exposures", y = outcome_name) +
    theme_bw() +
    coord_flip()  +
    scale_color_manual(values = c( "red", "black"),
                       name = "") +
    #guides(color = "none")+
    theme(axis.title = element_text(size = 9),
          axis.text = element_text(size = 9),
          legend.text = element_text(size = 9),
          legend.title = element_text(size = 9),
          legend.position = "bottom",
          legend.box = "vertical",
          legend.justification = "center",
          legend.spacing.y = unit(0, "cm"),
          legend.spacing.x = unit(0, "cm"),
          legend.box.margin = margin(0,0,0,0, "cm"),
          legend.margin = margin(0,0,0,0, "cm"))
  
  return(results)
}

leg <- results_multi %>%
  filter(outcome == "Shannon diversity") %>%
  #filter(analysis == "confirmatory") %>%
  filter(model_type == "adjusted") %>%
  forestplot_shannon(outcome_name = "Shannon diversity") +
  theme(axis.text.y = element_blank(), 
        axis.title.y = element_blank())
leg <- get_legend(leg) %>% as_ggplot()

forestplot_alpha_confirmatory_1 <-
  results_multi %>%
  filter(outcome == "Specific richness") %>%
  #filter(analysis == "confirmatory") %>%
  filter(model_type == "adjusted") %>%
  forestplot(outcome_name = "Specific richness") +
  theme(legend.position = "none")

forestplot_alpha_confirmatory_2 <-
  results_multi %>%
  filter(outcome == "Shannon diversity") %>%
  #filter(analysis == "confirmatory") %>%
  filter(model_type == "adjusted") %>%
  forestplot_shannon(outcome_name = "Shannon diversity") +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank())+
  theme(legend.position = "none")

forestplot_alpha_confirmatory_3 <- 
  results_multi %>%
  filter(outcome == "Faith phylogenetic diversity") %>%
  #filter(analysis == "confirmatory") %>%
  filter(model_type == "adjusted") %>%
  forestplot(outcome_name = "Faith's phylogenetic diversity") +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank())+
  theme(legend.position = "none")


forestplot_alpha_confirmatory <- 
  (forestplot_alpha_confirmatory_1 + forestplot_alpha_confirmatory_2 + forestplot_alpha_confirmatory_3) / leg + 
  plot_layout(heights = c(14, 1))
forestplot_alpha_confirmatory

ggsave("4_output/Review Coauteurs/forestplot_alpha_phthalates.tiff", 
       plot = forestplot_alpha_confirmatory, 
       device = "tiff",
       units = "mm",
       width = 180, 
       height = 160,
       dpi = 300,
       limitsize = FALSE)

#### Taxonomy -----
forestplot <- function(results_list, outcome_name) {
  results <- results_list %>%
    ggplot(aes(x = exposure, 
               y = estimate, 
               min = conf.low, 
               ymax = conf.high, 
               #color = interaction(exposure_window, term_2), 
               color = FWER.p.value_shape_taxa
               )) +
    geom_hline(yintercept = 0, linetype="dashed") +
    geom_pointrange(position = position_dodge(width = 0.4), size = 0.4) +
    labs(x = "Exposures", y = outcome_name) +
    theme_bw() +
    coord_flip()  +
    scale_color_manual(values = c(  "orange","red", "grey", "black"),
                       name = "") +
    guides(color = guide_legend(title = ""))+
    theme(axis.title = element_text(size = 7),
          axis.text = element_text(size = 6),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 6), 
          legend.position = "bottom",
          legend.box = "vertical", 
          legend.justification = "center", 
          legend.spacing.y = unit(0, "cm"), 
          legend.spacing.x = unit(0, "cm"), 
          legend.box.margin = margin(0,0,0,0, "cm"), 
          legend.margin = margin(0,0,0,0, "cm"))
  
  return(results)
}

##### Confirmatory ----

leg_taxonomy <-  
  results_multi %>%
  filter(outcome == "Firmicutes") %>%
  filter(exposure %in% c("MEP Y1", "DEHP Y1", "ohMPHP Y1")) %>%
  filter(model_type == "adjusted") %>%
  mutate(
    FWER.p.value_shape_taxa = fct_relevel(FWER.p.value_shape_taxa, 
                                          "p>0.006 (pollutants associated with α-div and corrected for multiple testing)",
                                          "p<0.006 (pollutants associated with α-div and corrected for multiple testing)",
                                          "p<0.05 (pollutants not associated with α-div and not corrected for multiple testing)",
                                          "p>0.05 (pollutants not associated with α-div and not corrected for multiple testing)"
    )
  ) %>% 
  forestplot(outcome_name = "Firmicutes") +
  guides(color = guide_legend(direction = "vertical"))
leg_taxonomy <- get_legend(leg_taxonomy) %>% as_ggplot()

forestplot_alpha_taxonomy_1 <-
  results_multi %>%
  filter(outcome == "Firmicutes") %>%
  filter(exposure %in% c("MEP Y1", "DEHP Y1", "ohMPHP Y1")) %>%
  filter(model_type == "adjusted") %>%
  mutate(
    FWER.p.value_shape_taxa = fct_relevel(FWER.p.value_shape_taxa, 
                                          "p>0.006 (pollutants associated with α-div and corrected for multiple testing)",
                                          "p<0.006 (pollutants associated with α-div and corrected for multiple testing)",
                                          "p<0.05 (pollutants not associated with α-div and not corrected for multiple testing)",
                                          "p>0.05 (pollutants not associated with α-div and not corrected for multiple testing)"
    )
  ) %>% 
  forestplot(outcome_name = "Firmicutes")+
  theme(legend.position = "none")



forestplot_alpha_taxonomy_2 <-
  results_multi %>%
  filter(outcome == "Actinobacteria") %>%
  filter(exposure %in% c("MEP Y1", "DEHP Y1", "ohMPHP Y1")) %>%
  filter(model_type == "adjusted") %>%
  mutate(
    FWER.p.value_shape_taxa = fct_relevel(FWER.p.value_shape_taxa, 
        "p>0.006 (pollutants associated with α-div and corrected for multiple testing)",
        "p<0.006 (pollutants associated with α-div and corrected for multiple testing)",
        "p<0.05 (pollutants not associated with α-div and not corrected for multiple testing)",
        "p>0.05 (pollutants not associated with α-div and not corrected for multiple testing)"
      )
  ) %>% 
  forestplot(outcome_name = "Actinobacteria") +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank())+
  theme(legend.position = "none")

forestplot_alpha_taxonomy_3 <- 
  results_multi %>%
  filter(outcome == "Bacteroidetes") %>%
  filter(exposure %in% c("MEP Y1", "DEHP Y1", "ohMPHP Y1")) %>%
  filter(model_type == "adjusted") %>%
  mutate(
    FWER.p.value_shape_taxa = fct_relevel(FWER.p.value_shape_taxa, 
                                          "p>0.006 (pollutants associated with α-div and corrected for multiple testing)",
                                          "p<0.006 (pollutants associated with α-div and corrected for multiple testing)",
                                          "p<0.05 (pollutants not associated with α-div and not corrected for multiple testing)",
                                          "p>0.05 (pollutants not associated with α-div and not corrected for multiple testing)"
    )
  ) %>% 
  forestplot(outcome_name = "Bacteroidetes") +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank())+
  theme(legend.position = "none")

forestplot_alpha_taxonomy_4 <- 
  results_multi %>%
  filter(outcome == "Proteobacteria") %>%
  filter(exposure %in% c("MEP Y1", "DEHP Y1", "ohMPHP Y1")) %>%
  filter(model_type == "adjusted") %>%
  mutate(
    FWER.p.value_shape_taxa = fct_relevel(FWER.p.value_shape_taxa, 
                                          "p>0.006 (pollutants associated with α-div and corrected for multiple testing)",
                                          "p<0.006 (pollutants associated with α-div and corrected for multiple testing)",
                                          "p<0.05 (pollutants not associated with α-div and not corrected for multiple testing)",
                                          "p>0.05 (pollutants not associated with α-div and not corrected for multiple testing)"
    )
  ) %>% 
  forestplot(outcome_name = "Proteobacteria") +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank())+
  theme(legend.position = "none")


forestplot_taxonomy <- 
  (forestplot_alpha_taxonomy_1 + forestplot_alpha_taxonomy_2 + forestplot_alpha_taxonomy_3 + forestplot_alpha_taxonomy_4) + 
  plot_layout(ncol = 4)

forestplot_taxonomy <- forestplot_taxonomy + leg_taxonomy + plot_layout(nrow = 2, ncol = 4, heights = c(4, 2))

ggsave("4_output/Review Coauteurs/forestplot_taxonomy.tiff", 
       plot = forestplot_taxonomy, 
       device = "tiff",
       units = "mm",
       width = 180, 
       height = 60,
       dpi = 300,
       limitsize = FALSE)

##### Explo ----


### Figures article - Heatmap correlation entre pollutants ----
bdd_phthalates <- bdd_alpha %>% 
  select(all_of(phthalates_vec)) %>%
  select(!contains(c("MEOHP", "MECPP", "MEHHP", "MEHP", "MMCHP",     
                     "ohMiNP", "oxoMiNP", "cxMiNP",                 
                     "ohMINCH", "oxoMINCH"))) %>%  
  select(!contains("cat"))

colnames(bdd_phthalates) <- colnames(bdd_phthalates) %>%
  str_replace_all(c("mo_" = "",
                    "ch_" = "",
                    "_total_i_cor_" = " ",
                    "_i_cor_" = " ", 
                    "_ln" = "",
                    "_conj" = "",
                    "ln" = "",
                    "_cor" = "", 
                    "_t2" = " t2", 
                    "_t3" = " t3", 
                    "_M2" = " M2", 
                    "_Y1" = " Y1", 
                    "_ms" = "", 
                    "_cat M2_2" = " M2", 
                    "DINCH Y1" = "ΣDINCH Y1",
                    "DINCH t3" = "ΣDINCH t3",
                    "DINCH t2" = "ΣDINCH t2",
                    "DiNP Y1" = "ΣDiNP Y1",
                    "DiNP M2" = "ΣDiNP M2",
                    "DiNP t3" = "ΣDiNP t3",
                    "DiNP t2" = "ΣDiNP t2",
                    "DEHP Y1" = "ΣDEHP Y1",
                    "DEHP M2" = "ΣDEHP M2",
                    "DEHP t3" = "ΣDEHP t3",
                    "DEHP t2" = "ΣDEHP t2")) 

heatmap_phthalates <- heatmap_cor(bdd_phthalates, decimal = 1)
ggsave(filename = "4_output/phthalates/heatmap_cor_phthalates.tiff",
       plot = heatmap_phthalates, 
       units = "mm",
       width = 250, 
       height = 250,
       dpi = 300,
       limitsize = FALSE)

### Figures article - Heatmap correlation pollutants and covariates  ----
vars_1 <- metadata %>% 
  select(all_of(phthalates_vec_num)) %>%
  select(!contains(c("MEOHP", "MECPP", "MEHHP", "MEHP", "MMCHP",     
                     "ohMiNP", "oxoMiNP", "cxMiNP",                 
                     "ohMINCH", "oxoMINCH"))) %>%  
  colnames()
cormat <- metadata %>% 
  select(all_of(vars_1),
         all_of(covar_vec_num))
colnames(cormat) <- colnames(cormat) %>%
  str_replace_all(
    c("ch_feces_age_w_Y1" = "Child age (weeks)",
      "ch_antibio_Y1" = "Antibiotics use 0-12 months",
      "mo_par" = "Maternal parity",
      "po_w_kg" = "Birth weight (kg)",
      "po_he"= "Birth length (cm)", 
      "ch_w_Y1"="Weight at one year (Kg)", 
      "ch_he_Y1"="Length at one year (cm)", 
      "po_gd"= "Gestational age (weeks)", 
      "mo_age"="Maternal age before pregnancy", 
      "mo_bmi_bepr"="Maternal BMI before pregnancy", 
      "bf_duration_till48w"="Breastfeeding duration (weeks)"))


heatmap <- heatmap_cor_pairwise(data = cormat, 
                     vars_1 = vars_1, 
                     vars_2 = c("Child age (weeks)",
                                "Antibiotics use 0-12 months",
                                "Maternal parity",
                                "Birth weight (kg)",
                                "Birth length (cm)", 
                                "Weight at one year (Kg)", 
                                "Length at one year (cm)", 
                                "Gestational age (weeks)", 
                                "Maternal age before pregnancy", 
                                "Maternal BMI before pregnancy", 
                                "Breastfeeding duration (weeks)"), 
                     decimal = 1)

ggsave(filename = "4_output/heatmap_cor_phthalates_covar.tiff",
       plot = heatmap, 
       units = "mm",
       width = 100, 
       height = 150,
       dpi = 300,
       limitsize = FALSE)


## Analyses de sensibilité ----
### Rarefaction ----
phthalates_signi_vec <- 
  bdd_alpha %>% 
  select(all_of(phthalates_vec)) %>%
  select(ch_DEHP_ms_i_cor_Y1_ln, 
         ch_MnBP_i_cor_Y1_ln, 
         mo_DiNP_ms_i_cor_t2_ln,
         ch_DiNP_ms_i_cor_Y1_ln, 
         ch_MiBP_i_cor_Y1_ln, 
         ch_MBzP_i_cor_Y1_ln, 
         ch_MEP_i_cor_Y1_ln, 
         ch_ohMPHP_cat_M2_2, 
         ch_ohMPHP_i_cor_Y1_ln) %>%
  colnames()

bdd_col_1 <- bdd_alpha %>% filter(!is.na(ch_feces_Shannon_5000_ASV_Y1)) 
bdd_col_2_3 <- bdd_alpha %>% filter(!is.na(ch_feces_Shannon_10000_ASV_Y1))

sensi_obs_rare_phthalates <-    
  tbl_merge(
    tbls = list(model(exposure_vec = phthalates_signi_vec, outcome = ch_feces_SpecRich_5000_ASV_Y1, data = bdd_col_1, digit_beta_IC = 1), 
                model(exposure_vec = phthalates_signi_vec, outcome = ch_feces_SpecRich_5000_ASV_Y1, data = bdd_col_2_3, digit_beta_IC = 1), 
                model(exposure_vec = phthalates_signi_vec, outcome = ch_feces_SpecRich_10000_ASV_Y1, data = bdd_col_2_3, digit_beta_IC = 1)),       
    tab_spanner = c("**Threshold 5,000 (n=350)**", "**Threshold 5,000 (n=339)**", "**Threshold 10,000 (n=339)**"))

sensi_sha_rare_phthalates <-    
  tbl_merge(
    tbls = list(model(exposure_vec = phthalates_signi_vec, outcome = ch_feces_Shannon_5000_ASV_Y1, data = bdd_col_1, digit_beta_IC = 2), 
                model(exposure_vec = phthalates_signi_vec, outcome = ch_feces_Shannon_5000_ASV_Y1, data = bdd_col_2_3, digit_beta_IC = 2), 
                model(exposure_vec = phthalates_signi_vec, outcome = ch_feces_Shannon_10000_ASV_Y1, data = bdd_col_2_3, digit_beta_IC = 2)),       
    tab_spanner = c("**Threshold 5,000 (n=350)**", "**Threshold 5,000 (n=339)**", "**Threshold 10,000 (n=339)**"))

sensi_fai_rare_phthalates <-    
  tbl_merge(
    tbls = list(model(exposure_vec = phthalates_signi_vec, outcome = ch_feces_Faith_5000_ASV_Y1, data = bdd_col_1, digit_beta_IC = 1), 
                model(exposure_vec = phthalates_signi_vec, outcome = ch_feces_Faith_5000_ASV_Y1, data = bdd_col_2_3, digit_beta_IC = 1), 
                model(exposure_vec = phthalates_signi_vec, outcome = ch_feces_Faith_10000_ASV_Y1, data = bdd_col_2_3, digit_beta_IC = 1)), 
    tab_spanner = c("**Threshold 5,000 (n=350)**", "**Threshold 5,000 (n=339)**", "**Threshold 10,000 (n=339)**"))


### Gravité spécifique ----
phthalates_signi_vec <- 
  bdd_alpha %>% 
  select(all_of(phthalates_vec)) %>%
  select(ch_DEHP_ms_i_cor_Y1_ln, 
         ch_MnBP_i_cor_Y1_ln, 
         mo_DiNP_ms_i_cor_t2_ln,
         ch_DiNP_ms_i_cor_Y1_ln, 
         ch_MiBP_i_cor_Y1_ln, 
         ch_MBzP_i_cor_Y1_ln, 
         ch_MEP_i_cor_Y1_ln, 
         ch_ohMPHP_i_cor_Y1_ln) %>%
  colnames()


phthalates_signi_sg_vec <- 
  bdd_alpha %>% 
  select(all_of(phthalates_vec_num_sg_ln)) %>%
  select(ch_DEHP_ms_i_cor_sg_Y1_ln,
         ch_MnBP_i_cor_sg_Y1_ln, 
         mo_DiNP_ms_i_cor_sg_t2_ln,
         ch_DiNP_ms_i_cor_sg_Y1_ln, 
         ch_MiBP_i_cor_sg_Y1_ln, 
         ch_MBzP_i_cor_sg_Y1_ln, 
         ch_MEP_i_cor_sg_Y1_ln, 
         ch_ohMPHP_i_cor_sg_Y1_ln) %>%
  colnames()

sensi_sg_alpha_phthalates <-
  tbl_merge(
    tbls = list(
      model(data = bdd_alpha, outcome = ch_feces_SpecRich_5000_ASV_Y1, exposure_vec = phthalates_signi_vec, digit_beta_IC = 1), 
      model(data = bdd_alpha, outcome = ch_feces_SpecRich_5000_ASV_Y1, exposure_vec = phthalates_signi_sg_vec, digit_beta_IC = 1), 
      
      model(data = bdd_alpha, outcome = ch_feces_Shannon_5000_ASV_Y1, exposure_vec = phthalates_signi_vec, digit_beta_IC = 2), 
      model(data = bdd_alpha, outcome = ch_feces_Shannon_5000_ASV_Y1, exposure_vec = phthalates_signi_sg_vec, digit_beta_IC = 2), 
      
      model(data = bdd_alpha, outcome = ch_feces_Faith_5000_ASV_Y1, exposure_vec = phthalates_signi_vec, digit_beta_IC = 1), 
      model(data = bdd_alpha, outcome = ch_feces_Faith_5000_ASV_Y1, exposure_vec = phthalates_signi_sg_vec, digit_beta_IC = 1)), 
    tab_spanner = c(
      "**Specific richness, principal analysis, fully adjusted**", 
      "**Specific richness, sensitivity analysis, fully adjusted + sg**", 
      "**Shannon diversity, principal analysis, fully adjusted**", 
      "**Shannon diversity, sensitivity analysis, fully adjusted + sg**", 
      "**Faith diversity, principal analysis, fully adjusted**", 
      "**Faith diversity, sensitivity analysis, fully adjusted + sg**"))

sensi_sg_phyla_phthalates <- tbl_merge(
  tbls = list(
    model(data = bdd_taxa, outcome = ch_feces_rel_p1_Y1, exposure_vec = phthalates_signi_vec, digit_beta_IC = 1), 
    model(data = bdd_taxa, outcome = ch_feces_rel_p1_Y1, exposure_vec = phthalates_signi_sg_vec, digit_beta_IC = 1), 
    
    model(data = bdd_taxa, outcome = ch_feces_rel_p2_Y1, exposure_vec = phthalates_signi_vec, digit_beta_IC = 1), 
    model(data = bdd_taxa, outcome = ch_feces_rel_p2_Y1, exposure_vec = phthalates_signi_sg_vec, digit_beta_IC = 1), 
    
    model(data = bdd_taxa, outcome = ch_feces_rel_p3_Y1, exposure_vec = phthalates_signi_vec, digit_beta_IC = 1), 
    model(data = bdd_taxa, outcome = ch_feces_rel_p3_Y1, exposure_vec = phthalates_signi_sg_vec, digit_beta_IC = 1), 
    
    model(data = bdd_taxa, outcome = ch_feces_rel_p4_Y1, exposure_vec = phthalates_signi_vec, digit_beta_IC = 1), 
    model(data = bdd_taxa, outcome = ch_feces_rel_p4_Y1, exposure_vec = phthalates_signi_sg_vec, digit_beta_IC = 1)), 
  tab_spanner = c(
    "**Phylum Firmicutes, principal analysis, fully adjusted**",
    "**Phylum Firmicutes, sensitivity analysis, fully adjusted + sg**", 
    "**Phylum Actinobacteria, principal analysis, fully adjusted**",
    "**Phylum Actinobacteria, sensitivity analysis, fully adjusted + sg**", 
    "**Phylum Bacteroidetes, principal analysis, fully adjusted**",
    "**Phylum Bacteroidetes, sensitivity analysis, fully adjusted + sg**", 
    "**Phylum Proteobacteria, principal analysis, fully adjusted**",
    "**Phylum Proteobacteria, sensitivity analysis, fully adjusted + sg**"))


sensi_sg_genera_phthalates <- tbl_merge(
  tbls = list(
    model(data = bdd_taxa, outcome = ch_feces_rel_g1_Y1, exposure_vec = phthalates_signi_vec, digit_beta_IC = 1), 
    model(data = bdd_taxa, outcome = ch_feces_rel_g1_Y1, exposure_vec = phthalates_signi_sg_vec, digit_beta_IC = 1), 
    
    model(data = bdd_taxa, outcome = ch_feces_rel_g2_Y1, exposure_vec = phthalates_signi_vec, digit_beta_IC = 1), 
    model(data = bdd_taxa, outcome = ch_feces_rel_g2_Y1, exposure_vec = phthalates_signi_sg_vec, digit_beta_IC = 1), 
    
    model(data = bdd_taxa, outcome = ch_feces_rel_g3_Y1, exposure_vec = phthalates_signi_vec, digit_beta_IC = 1), 
    model(data = bdd_taxa, outcome = ch_feces_rel_g3_Y1, exposure_vec = phthalates_signi_sg_vec, digit_beta_IC = 1), 
    
    model(data = bdd_taxa, outcome = ch_feces_rel_g4_Y1, exposure_vec = phthalates_signi_vec, digit_beta_IC = 1), 
    model(data = bdd_taxa, outcome = ch_feces_rel_g4_Y1, exposure_vec = phthalates_signi_sg_vec, digit_beta_IC = 1)), 
  tab_spanner = c(
    "**Genus Bifidobacterium, principal analysis, fully adjusted**",
    "**Genus Bifidobacterium, sensitivity analysis, fully adjusted + sg**", 
    "**Genus Bacteroides, principal analysis, fully adjusted**",
    "**Genus Bacteroides, sensitivity analysis, fully adjusted + sg**", 
    "**Genus Blautia, principal analysis, fully adjusted**",
    "**Genus Blautia, sensitivity analysis, fully adjusted + sg**", 
    "**Genera Eschericha / Shigella, principal analysis, fully adjusted**",
    "**Genera Eschericha / Shigella, sensitivity analysis, fully adjusted + sg**"
  ))


## Additional file 2 modif pour somme molaire ----
test <- metadata %>% 
  select(statut, 
         mo_DEHP_ms_i_cor_t2, mo_DEHP_ms_i_cor_t3, ch_DEHP_ms_i_cor_M2, ch_DEHP_ms_i_cor_Y1, 
         mo_DiNP_ms_i_cor_t2, mo_DiNP_ms_i_cor_t3, ch_DiNP_ms_i_cor_M2, ch_DiNP_ms_i_cor_Y1, 
         mo_DINCH_ms_i_cor_t2, mo_DINCH_ms_i_cor_t3, ch_DINCH_ms_i_cor_Y1) %>% 
  tbl_summary(by = "statut") %>%
  add_p() %>%
  bold_labels()

test <- metadata %>% filter(statut == "inclu")
test <- descrip_num(data = test, vars = c("mo_DEHP_ms_i_cor_t2", "mo_DEHP_ms_i_cor_t3", "ch_DEHP_ms_i_cor_M2", "ch_DEHP_ms_i_cor_Y1", 
                                  "mo_DiNP_ms_i_cor_t2", "mo_DiNP_ms_i_cor_t3", "ch_DiNP_ms_i_cor_M2", "ch_DiNP_ms_i_cor_Y1", 
                                  "mo_DINCH_ms_i_cor_t2", "mo_DINCH_ms_i_cor_t3", "ch_DINCH_ms_i_cor_Y1"))
writexl::write_xlsx(test, "C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/7. Présentations écrites/5. Article_phthalates_microbiote_Y1/tableauS1.xlsx")
boxplot(data = test, vars = c("mo_DEHP_ms_i_cor_t2", "mo_DEHP_ms_i_cor_t3", "ch_DEHP_ms_i_cor_M2", "ch_DEHP_ms_i_cor_Y1", 
                                  "mo_DiNP_ms_i_cor_t2", "mo_DiNP_ms_i_cor_t3", "ch_DiNP_ms_i_cor_M2", "ch_DiNP_ms_i_cor_Y1", 
                                  "mo_DINCH_ms_i_cor_t2", "mo_DINCH_ms_i_cor_t3", "ch_DINCH_ms_i_cor_Y1"))

library(haven)
base_aline_211115 <- read_sas("0_source_data/base_aline_211115.sas7bdat", 
                              NULL)
test <- base_aline_211115 %>% filter(statut == "inclu")
descrip_num(data = base_aline_211115, vars = c("mo_DEHP_ms_i_cor_t2", "mo_DEHP_ms_i_cor_t3", "ch_DEHP_ms_i_cor_M2", "ch_DEHP_ms_i_cor_Y1", 
                                  "mo_DiNP_ms_i_cor_t2", "mo_DiNP_ms_i_cor_t3", "ch_DiNP_ms_i_cor_M2", "ch_DiNP_ms_i_cor_Y1", 
                                  "mo_DINCH_ms_i_cor_t2", "mo_DINCH_ms_i_cor_t3", "ch_DINCH_ms_i_cor_Y1"))
boxplot(data = test, vars = c("mo_DEHP_ms_i_cor_t2", "mo_DEHP_ms_i_cor_t3", "ch_DEHP_ms_i_cor_M2", "ch_DEHP_ms_i_cor_Y1", 
                              "mo_DiNP_ms_i_cor_t2", "mo_DiNP_ms_i_cor_t3", "ch_DiNP_ms_i_cor_M2", "ch_DiNP_ms_i_cor_Y1", 
                              "mo_DINCH_ms_i_cor_t2", "mo_DINCH_ms_i_cor_t3", "ch_DINCH_ms_i_cor_Y1"))


a_convertir <- c("mo_DEHP_ms_i_cor_t2", "mo_DEHP_ms_i_cor_t3", "ch_DEHP_ms_i_cor_M2", "ch_DEHP_ms_i_cor_Y1", 
                 "mo_DiNP_ms_i_cor_t2", "mo_DiNP_ms_i_cor_t3", "ch_DiNP_ms_i_cor_M2", "ch_DiNP_ms_i_cor_Y1", 
                 "mo_DINCH_ms_i_cor_t2", "mo_DINCH_ms_i_cor_t3", "ch_DINCH_ms_i_cor_Y1")
base_aline_211115[a_convertir] <- lapply(base_aline_211115[a_convertir], as.numeric)


base_aline_211115 %>% select("mo_DEHP_ms_i_cor_t2", "mo_DEHP_ms_i_cor_t3", "ch_DEHP_ms_i_cor_M2", "ch_DEHP_ms_i_cor_Y1", 
                             "mo_DiNP_ms_i_cor_t2", "mo_DiNP_ms_i_cor_t3", "ch_DiNP_ms_i_cor_M2", "ch_DiNP_ms_i_cor_Y1", 
                             "mo_DINCH_ms_i_cor_t2", "mo_DINCH_ms_i_cor_t3", "ch_DINCH_ms_i_cor_Y1") %>% 
  tbl_summary()


## Additional file 2 ----
bdd_test <- metadata %>% filter(statut == "inclu")
additional_file_2.1 <- descrip_num(data = bdd_test, vars = phthalates_vec_num)
writexl::write_xlsx(additional_file_2.1, 
                    "C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/5. R projects/pollutants_gut_microbiota_Y1/4_output/phthalates/AdditionalFile2.xlsx")


phthalates_descrip_vec <- metadata %>% select(all_of(phthalates_vec)) %>% colnames() %>% str_replace("_ln", "")
additional_file_2.2 <- metadata %>%
  select(all_of(phthalates_descrip_vec),
         statut) %>%
  mutate(statut = fct_relevel(statut, "inclu", "exclu")) %>%
  tbl_summary(by = statut,
              missing = "no", 
              digits = list(all_continuous() ~ 1, 
                            all_categorical()~ 0)) %>%
  add_p()
