# A. Davias
# 08/02/2023

# Chargement des données et des fonctions existantes ----
load("2_final_data/metadata.RData")
load("2_final_data/bdd_alpha.RData")
load("2_final_data/bdd_taxa.RData")
source("3_programs/4_functions_AD_gumme.R", encoding = 'UTF-8')      # fonctions
rm(comp_effectifs, heatmap_cor_pairwise, model_covar, model_multi, model_summary, model_univ_multi, 
   table_cor, table_cor_sg, test_sensi_sg)
source("3_programs/4_vectors_AD_gumme.R", echo=TRUE)
rm(pollutant_vec_t2, pollutant_vec_t3, pollutant_vec_M2, pollutant_vec_Y1, 
   covar_vec, covar_vec_cat, covar_vec_cat_i,
   phthalates_vec_cat, phthalates_vec_ln, phthalates_vec_ter, 
   taxa_vec, 
   list = ls()[grep("sg", ls())])
rm(list = ls()[grep("phenols", ls())])
rm(list = ls()[grep("pfas", ls())])
rm(list = ls()[grep("num", ls())])
library(psych)
library(writexl)


# Création de fonctions ----

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



# Création des vecteurs ----
pollutants_vec <- bdd_alpha %>% 
  select(all_of(phthalates_vec)) %>% 
  select(!contains(c("MEOHP", "MECPP", "MEHHP", "MEHP", "MMCHP", "ohMiNP", "oxoMiNP", "cxMiNP", "ohMINCH", "oxoMINCH"))) %>% 
  select(!contains("M2")) %>%
  colnames()
alpha_vec <- bdd_alpha %>% 
  select("ch_feces_SpecRich_5000_ASV_Y1", "ch_feces_Shannon_5000_ASV_Y1") %>% 
  colnames() 
rm(phthalates_vec)



# Results - Unipolluants analysis ----
## Tableaux bruts ----
results_list_1 <- lapply(bdd_alpha[, pollutants_vec], lm_func, outcome = bdd_alpha$ch_feces_SpecRich_5000_ASV_Y1, data = bdd_alpha)
results_list_2 <- lapply(bdd_alpha[, pollutants_vec], lm_func, outcome = bdd_alpha$ch_feces_Shannon_5000_ASV_Y1, data = bdd_alpha)

results_list_3 <- lapply(bdd_taxa[, pollutants_vec], lm_func, outcome = bdd_taxa$ch_feces_rel_p1_Y1, data = bdd_taxa)
results_list_4 <- lapply(bdd_taxa[, pollutants_vec], lm_func, outcome = bdd_taxa$ch_feces_rel_p2_Y1, data = bdd_taxa)
results_list_5 <- lapply(bdd_taxa[, pollutants_vec], lm_func, outcome = bdd_taxa$ch_feces_rel_p3_Y1, data = bdd_taxa)
results_list_6 <- lapply(bdd_taxa[, pollutants_vec], lm_func, outcome = bdd_taxa$ch_feces_rel_p4_Y1, data = bdd_taxa)

results_multi <- list(
  
  data_prep(results_list = results_list_1, outcome_name = "Specific richness"),
  data_prep(results_list = results_list_2, outcome_name = "Shannon diversity"),
  
  data_prep(results_list = results_list_3, outcome_name = "Firmicutes"),
  data_prep(results_list = results_list_4, outcome_name = "Actinobacteria"),
  data_prep(results_list = results_list_5, outcome_name = "Bacteroidetes"),
  data_prep(results_list = results_list_6, outcome_name = "Proteobacteria")) %>%
  
  bind_rows() %>%
  mutate(
    exposure = as.factor(exposure)) %>%
  mutate(
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
                                 "_Y1" = " Y1", 
                                 "_ms" = ""
                                 )), 

    p_value_shape = case_when(p.value < 0.1 ~ "p.value <0.1",
                              p.value > 0.1 ~ "p.value >0.1"), 
    p_value_shape = fct_relevel(p_value_shape,
                                "p.value >0.1", 
                                "p.value <0.1"),
    exposure_window = case_when(grepl("t2", exposure) ~ "Trim.2", 
                                grepl("t3", exposure) ~ "Trim.3", 
                                grepl("Y1", exposure) ~ "12 months", 
                                TRUE ~ "Trim.2"), 
    exposure_window = fct_relevel(exposure_window, 
                                  "Trim.2", 
                                  "Trim.3", 
                                  "12 months")) %>%
  select("model_type", 
         "outcome", 
         "exposure",
         "exposure_window",
         "term",
         "estimate",
         "std.error",
         "statistic",
         "conf.low",
         "conf.high",
         "p.value", 
         "p_value_shape")

rm(results_list_1, results_list_2, results_list_3, results_list_4, results_list_5, results_list_6)


## Correction pour comparaison multiple (seulement les tests alpha diversité pris en compte) ----
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
## on passe de 2 outcomes à 2 outcomes après prise en compte de leur corrélation

results_multi <- results_multi %>%
  mutate(
    FWER.p.value_alpha = ifelse(outcome %in% c("Specific richness", "Shannon diversity"), 
                          p.value * 14 * 3, NA),
    FWER.p.value_alpha = ifelse(FWER.p.value_alpha > 1, ">0.99", FWER.p.value_alpha),
    FWER.p.value_shape_alpha = case_when(p.value< 0.0012 & 
                                     outcome %in% c("Specific richness", 
                                                    "Shannon diversity")~ "p.value <0.0012",
                                   p.value > 0.0012 & 
                                     outcome %in% c("Specific richness", 
                                                    "Shannon diversity")~ "p.value >0.0012"), 
    FWER.p.value_shape_alpha = fct_relevel(FWER.p.value_shape_alpha,
                                "p.value >0.0012", 
                                "p.value <0.0012"))
rm(test_expo, test_outcome, 
   cor_mixed_table_expo, cor_mixed_table_outcome, 
   results_alpha_corrected_expo, results_alpha_corrected_outcome, 
   results_M0_corrected_expo, results_M0_corrected_outcome)

## Tableaux article -----
### Table 1 : alpha diversité en régressions linéaires ----
table_1 <- tbl_merge(                         
  tbls = list(
    effectif_column(data = bdd_alpha, 
                    outcome = ch_feces_SpecRich_5000_ASV_Y1, 
                    exposure_vec = pollutants_vec),
    model(data = bdd_alpha, 
          outcome = ch_feces_SpecRich_5000_ASV_Y1,
          exposure_vec = pollutants_vec, 
          digit_beta_IC = 1),
    model(data = bdd_alpha,
          outcome = ch_feces_Shannon_5000_ASV_Y1,
          exposure_vec = pollutants_vec, 
          digit_beta_IC = 2)), 
  tab_spanner = c("", 
                  "**Specific richness**", 
                  "**Shannon diversity**"))


### Table 2 : taxonomie en régressions linéaires ----
# Analyses "explicatives" non corrigées pour comparaison multiple pour la taxonomie 
# pour les associations significatives pour l'alpha diversité
# pour les autres polluants non significatifs --> on met tout en analyse sup. 
signi_vec <- c("ch_DEHP_ms_i_cor_Y1_ln","ch_MEP_i_cor_Y1_ln", "ch_ohMPHP_i_cor_Y1_ln")
not_signi_vec <- bdd_alpha %>% 
  select(all_of(pollutants_vec)) %>% 
  select(-c("ch_DEHP_ms_i_cor_Y1_ln", "ch_MEP_i_cor_Y1_ln", "ch_ohMPHP_i_cor_Y1_ln")) %>%
  colnames()


table_2 <- tbl_merge(
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

### Table A.9 : taxonomie en régressions linéaires ----
table_A9 <- tbl_merge(
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

## Figures article ----
### Figure 1 : Alpha diversity ----
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
  filter(model_type == "adjusted") %>%
  forestplot_shannon(outcome_name = "Shannon diversity") +
  theme(axis.text.y = element_blank(), 
        axis.title.y = element_blank())
leg <- get_legend(leg) %>% as_ggplot()

forestplot_alpha_1 <-
  results_multi %>%
  filter(outcome == "Specific richness") %>%
  filter(model_type == "adjusted") %>%
  forestplot(outcome_name = "Specific richness") +
  theme(legend.position = "none")

forestplot_alpha_2 <-
  results_multi %>%
  filter(outcome == "Shannon diversity") %>%
  filter(model_type == "adjusted") %>%
  forestplot_shannon(outcome_name = "Shannon diversity") +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank())+
  theme(legend.position = "none")


fig_1 <- 
  (forestplot_alpha_1 + forestplot_alpha_2) / leg + 
  plot_layout(heights = c(14, 1))

rm(forestplot_alpha_1, forestplot_alpha_2, leg)

ggsave("4_output/Figure 1 (forestplot_alpha_phthalates).tiff", 
       plot = fig_1, 
       device = "tiff",
       units = "mm",
       width = 180, 
       height = 160,
       dpi = 300,
       limitsize = FALSE)


## Analyses de sensibilité ----
### Rarefaction ----
phthalates_signi_vec <- 
  bdd_alpha %>% 
  select(all_of(pollutants_vec)) %>%
  select(ch_DEHP_ms_i_cor_Y1_ln, 
         ch_MnBP_i_cor_Y1_ln, 
         mo_DiNP_ms_i_cor_t2_ln,
         ch_DiNP_ms_i_cor_Y1_ln, 
         ch_MiBP_i_cor_Y1_ln, 
         ch_MBzP_i_cor_Y1_ln, 
         ch_MEP_i_cor_Y1_ln, 
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
      model(data = bdd_alpha, outcome = ch_feces_Shannon_5000_ASV_Y1, exposure_vec = phthalates_signi_sg_vec, digit_beta_IC = 2)), 
    tab_spanner = c(
      "**Specific richness, principal analysis, fully adjusted**", 
      "**Specific richness, sensitivity analysis, fully adjusted + sg**", 
      "**Shannon diversity, principal analysis, fully adjusted**", 
      "**Shannon diversity, sensitivity analysis, fully adjusted + sg**"))

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
