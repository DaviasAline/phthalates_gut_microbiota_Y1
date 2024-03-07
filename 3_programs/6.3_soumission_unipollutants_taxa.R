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

# Data reading ----
## Gut microbiota ----
taxa_table <- read_csv("0_source_data/taxa_table_ASVbased_Y1_AD_20220504_8.csv")
bdd_microbiota <- 
  read_labelled_csv(
    "0_source_data/gut_microbiota_ASVbased_Y1_labelled_AD_20220504_7239.csv")%>%
  select(
    ident, 
    starts_with("ch_feces_rel_")) %>%
  filter(!is.na(ch_feces_rel_g1_Y1)) 

input_data <- 
  read_labelled_csv(
    "0_source_data/gut_microbiota_ASVbased_Y1_labelled_AD_20220504_7239.csv")%>%
  select(
    ident, 
    starts_with("ch_feces_rel_g")) %>%
  filter(!is.na(ch_feces_rel_g1_Y1)) %>%
  column_to_rownames("ident")

var_label(input_data) <-  str_replace(
  var_label(input_data),                                    # set correct variable names                   
  "One year child feces relative abundance of ", "")
colnames(input_data) <- var_label(input_data)



## Vectors ----
covariates <- c("ch_feces_RUN_Y1",                   
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
                "ch_w_Y1_3cat_i",
                "ch_he_Y1_3cat_i",
                "po_gd",
                "mo_age",
                "mo_bmi_bepr_3cat_i",
                "bf_duration_till48w_4cat_i")
pollutants <- c("mo_DEHP_ms_i_cor_t2_ln", "mo_DEHP_ms_i_cor_t3_ln", "ch_DEHP_ms_i_cor_Y1_ln",   # Métabolites DEHP
                "mo_MnBP_i_cor_t2_ln", "mo_MnBP_i_cor_t3_ln", "ch_MnBP_i_cor_Y1_ln",            # Métabolite DBP
                "mo_DiNP_ms_i_cor_t2_ln", "mo_DiNP_ms_i_cor_t3_ln", "ch_DiNP_ms_i_cor_Y1_ln",   # Métabolites DiNP
                "mo_MiBP_i_cor_t2_ln", "mo_MiBP_i_cor_t3_ln", "ch_MiBP_i_cor_Y1_ln",            # Métabolite DiBP
                "mo_MBzP_i_cor_t2_ln", "mo_MBzP_i_cor_t3_ln", "ch_MBzP_i_cor_Y1_ln",            # Metabolite MBzP
                "mo_MEP_i_cor_t2_ln", "mo_MEP_i_cor_t3_ln", "ch_MEP_i_cor_Y1_ln",               # Metabolite DEP
                "mo_ohMPHP_i_cor_t2_ln", "mo_ohMPHP_i_cor_t3_ln", "ch_ohMPHP_i_cor_Y1_ln",      # Metabolite DEHP
                "mo_DINCH_ms_i_cor_t2_ln", "mo_DINCH_ms_i_cor_t3_ln", "ch_DINCH_ms_i_cor_Y1_ln")# Metabolites DINCH

genera_linear_complet <- bdd_microbiota %>%
  select(contains("ch_feces_rel_g")) %>%
  filter(!is.na(ch_feces_rel_g1_Y1)) %>%
  select_if(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.3) %>%
  colnames()

## Metadata ----
load("2_final_data/metadata.RData")
input_metadata <- metadata %>%
  select(ident,                      
         all_of(covariates), 
         all_of(pollutants)) %>%
  filter(!is.na(ch_feces_RUN_Y1)) %>%
  column_to_rownames("ident")
rm(metadata)

input_metadata <- input_metadata %>%
  mutate_at(vars(c("ch_feces_RUN_Y1",                   
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
                   "ch_w_Y1_3cat_i",
                   "ch_he_Y1_3cat_i",
                   "mo_bmi_bepr_3cat_i",
                   "bf_duration_till48w_4cat_i")), as.factor)


# Functions coding ----
# Fonction pour voir la table de correlation avant de faire la correction pour comparaison multiple selon correlation entre les tests
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

# Fonction pour le choix du nombre de décimales
custom_pvalue_fun <- function(x) {
  sapply(x, function(p) {
    if (is.na(p)) {
      return(NA) # Retourner NA si p est NA
    } else if (p < 0.001) {
      # Pour p < 0.001, utiliser la notation scientifique pour afficher toutes les décimales
      return(format(p, scientific = TRUE))
    } else if (p >= 0.001 & p < 0.01) {
      # Pour 0.001 <= p < 0.01, afficher avec 3 décimales
      return(sprintf("%.3f", p))
    } else {
      # Pour p >= 0.01, afficher avec 2 décimales
      return(sprintf("%.2f", p))
    }
  })
}


# Genera selection ----
## Description des genres à étudier en régression linéaire (n=46) 
input_data %>%
  select_if(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.3) %>%
  View()

genera_linear <- input_data %>%
  select_if(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.3) %>%
  colnames()

descrip_genera_linear <- tbl_merge(
  tbls = 
    list(
      tbl_1 = 
        input_data %>% 
        select_if(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.3) %>%
        tbl_summary(
          type = list(everything() ~ "continuous"), 
          statistic = list(everything() ~ "{median} ({p25}, {p75})"), 
          digits = list(all_continuous() ~ c(2, 1, 1))), 
      tbl_2 = 
        input_data %>% 
        select_if(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.3) %>%
        mutate_all(~ ifelse(.>0, "Yes", "No")) %>%
        set_label(genera_linear) %>%
        tbl_summary(
          type = list(everything() ~ "categorical"))), 
  tab_spanner = c("**Continuous**", "**Categorical (Y/N)**"))


## Description des genres à exclure (n=146) 
# input_data %>%
#   select_if(~ sum(. != 0, na.rm = TRUE) / length(.) < 0.3) %>%
#   View()
# 
# genera_excluded <- input_data %>%
#   select_if(~ sum(. != 0, na.rm = TRUE) / length(.) < 0.3) %>%
#   colnames()
# 
# descrip_genera_excluded <- tbl_merge(
#   tbls = 
#     list(
#       tbl_1 = 
#         input_data %>% 
#         select_if(~ sum(. != 0, na.rm = TRUE) / length(.) < 0.3) %>%
#         tbl_summary(
#           type = list(everything() ~ "continuous"), 
#           statistic = list(everything() ~ "{median} ({p25}, {p75})"), 
#           digits = list(all_continuous() ~ c(2, 1, 1))), 
#       tbl_2 = 
#         input_data %>% 
#         select_if(~ sum(. != 0, na.rm = TRUE) / length(.) < 0.3) %>%
#         mutate_all(~ ifelse(.>0, "Yes", "No")) %>%
#         set_label(genera_excluded) %>%
#         tbl_summary(
#           type = list(everything() ~ "categorical"))), 
#   tab_spanner = c("**Continuous**", "**Categorical (Y/N)**"))


# Data cleaning ----
## Selection des genres à analyser
# input_data_raw <- input_data %>% 
#   select(all_of(genera_linear)) %>%
#   mutate_all(., ~ ifelse(. == 0, 1/5000, .)) %>%        # remplacement des valeurs 0 par 1/5000
#   rename_with(~gsub("genus ", "", .), everything()) %>% # changement des noms de colonnes pour qu'ils n'aient pas d'espace
#   rownames_to_column(var = "ident")

input_data_log <- input_data %>% 
  select(all_of(genera_linear)) %>%
  mutate_all(., ~ ifelse(. == 0, 1/5000, .)) %>%        # remplacement des valeurs 0 par 1/5000
  mutate_all(~ log(.)) %>%                              # transformation logarithmique
  rename_with(~gsub("genus ", "", .), everything()) %>% # changement des noms de colonnes pour qu'ils n'aient pas d'espace
  rownames_to_column(var = "ident")

# input_data_clr <- input_data %>% 
#   select(all_of(genera_linear)) %>%
#   mutate_all(., ~ ifelse(. == 0, 1/5000, .)) %>%     # remplacement des valeurs 0 par 1/5000
#   mutate_all(~ log(.))
# rowMeans <- rowMeans(input_data_clr)
# input_data_clr <- sweep(input_data_clr, 1, rowMeans, FUN = "-")   # on soustrait chaque valeur par la log mean de la lign
# input_data_clr <- input_data_clr %>%
#   rename_with(~gsub("genus ", "", .), everything()) %>% # changement des noms de colonnes pour qu'ils n'aient pas d'espace
#   rownames_to_column(var = "ident")
# rm(rowMeans)

genera_linear <- str_replace_all(genera_linear, "genus ", "")

# Data description ----
## Covariates ----
descrip_covar <- 
  input_metadata %>% 
  select(all_of(covariates)) %>% 
  tbl_summary(type = list(all_categorical() ~ "categorical"))

## Included genera (n=46) ----
### Tables ----
descrip_genera_linear

# comp_transfo <- tbl_merge(
#   tbls = list(
#     input_data_raw %>% select(-ident) %>% tbl_summary(),
#     input_data_log %>% select(-ident) %>% tbl_summary(),
#     input_data_clr %>% select(-ident) %>% tbl_summary()), 
#   tab_spanner = c("**Raw**", "**Log**", "**Clr**"))
# comp_transfo

### Density plots ----
# densityplot(data = input_data_raw, vars = genera_linear[1:23])
# densityplot(data = input_data_raw, vars = genera_linear[24:46])
densityplot(data = input_data_log, vars = genera_linear[1:23])
densityplot(data = input_data_log, vars = genera_linear[24:46])
# densityplot(data = input_data_clr, vars = genera_linear[1:23])
# densityplot(data = input_data_clr, vars = genera_linear[24:46])


### Heatmap of correlation ----
cormat_genera <- input_data %>% 
  select_if(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.3) %>%
  rename_with(~gsub("genus ", "", .), everything())
cormat_genera <- round(cor(cormat_genera, 
                           use = "pairwise.complete.obs", 
                           method = "pearson"), 1)
heatmap_genera <-                                       # heatmap
  reshape2::melt(cormat_genera, na.rm = TRUE) %>% # passer en df long rapidement 
  ggplot(aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1, 1),
    space = "Lab",
    name = "Pearson\nCorrelation"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(
    angle = 45,
    vjust = 1,
    #size = 12,
    hjust = 1
  )) +
  coord_fixed() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.4, 0.7),
    legend.direction = "horizontal"
  ) +
  guides(fill = guide_colorbar(
    barwidth = 7,
    barheight = 1,
    title.position = "top",
    title.hjust = 0.5
  ))
heatmap_genera
rm(cormat_genera)

# Calcul du nombre de tests ----
## Pour les analydes de taxonomie : 50 outcomes (4 phyla et 46 genres) et 24 expositions 
test_outcome <- input_data_log %>%
  select(ident, all_of(genera_linear)) %>%
  na.omit()
test_outcome_2 <- bdd_microbiota %>% 
  select(ident, ch_feces_rel_p1_Y1, ch_feces_rel_p2_Y1, ch_feces_rel_p3_Y1, ch_feces_rel_p4_Y1) %>%
  mutate(ident = as.character(ident)) %>%
  na.omit()
test_outcome <- left_join(test_outcome, test_outcome_2, by = "ident") %>%
  select(-ident)
variables <- names(test_outcome)
for(var in variables) {
  var_lab(test_outcome[[var]]) <- NULL}
cor_mixed_table_outcome <- cor_mixed(data = test_outcome, method = "pearson")
results_alpha_corrected_outcome <- alpha_corrected(data = test_outcome, alpha = 0.05)
results_M0_corrected_outcome <- M0_corrected(data = test_outcome, alpha = 0.05)
results_M0_corrected_outcome
## on passe de 50 outcomes à 33 outcomes après prise en compte de leur corrélation

test_expo <- input_metadata %>%
  select(all_of(pollutants)) %>%
  na.omit()
variables <- names(test_expo)
for(var in variables) {
  var_lab(test_expo[[var]]) <- NULL}
cor_mixed_table_expo <- cor_mixed(data = test_expo, method = "pearson")
results_alpha_corrected_expo <- alpha_corrected(data = test_expo, alpha = 0.05)
results_M0_corrected_expo <- M0_corrected(data = test_expo, alpha = 0.05)
results_M0_corrected_expo
## on passe de 24 expo à 14 expo après prise en compte de leur corrélation


# Running linear regressions ----
## Rassemblement des données en 1 seul dataframe 
input_metadata <- input_metadata %>% rownames_to_column(var = "ident")
data_df_log <- left_join(input_data_log, input_metadata, by = "ident")
corres <- 
  taxa_table %>% 
  select(Phyla_corres = ch_feces_phylum_ASVbased_Y1, 
         Outcome = ch_feces_genus_ASVbased_Y1) %>%
  filter(Outcome %in% genera_linear) %>%
  distinct(Outcome, .keep_all = TRUE)

## Version log ----
tbls_by_outcome_log <- vector("list", length(genera_linear))                    # création liste pour stocker les tableaux par outcome
names(tbls_by_outcome_log) <- genera_linear

for (outcome in genera_linear) {
  tbls_for_outcome_log <- vector("list", length(pollutants))
  names(tbls_for_outcome_log) <- pollutants
  
  for (exposure in pollutants) {                                                # formula setting
    terms <- c(exposure, covariates)
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
    
    tbls_for_outcome_log[[exposure]] <- tbl
  }
  tbls_by_outcome_log[[outcome]] <- tbls_for_outcome_log
}

### Tableau pour l'article ----
stacked_tbls_by_outcome_log <- vector("list", length(genera_linear))            # Création liste pour stocker les tableaux empilés par outcome
names(stacked_tbls_by_outcome_log) <- genera_linear

for (outcome in names(tbls_by_outcome_log)) {                                   # Récupérer les tableaux de régression pour cet outcome
  tbls_for_this_outcome_log <- tbls_by_outcome_log[[outcome]]
  stacked_tbl_log <- do.call(tbl_stack, list(tbls = tbls_for_this_outcome_log)) # Empiler les tableaux en un seul tableau
  stacked_tbls_by_outcome_log[[outcome]] <- stacked_tbl_log                     # Ajouter le tableau empilé à la liste des tableaux empilés
}

results_tbl_log <- tbl_merge(tbls = stacked_tbls_by_outcome_log,                # Fusionner les tableaux empilés en un seul tableau
                             tab_spanner = genera_linear)


### Tableau pour générer des figures ----
table_log <- tibble()                                                           # Initialisation d'un tibble vide pour stocker les résultats finaux
for (i in seq_along(tbls_by_outcome_log)) {                                     # Nom de l'outcome pour cette itération
  outcome_name <- names(tbls_by_outcome_log)[i]
  
  for (j in seq_along(tbls_by_outcome_log[[i]])) {                              # Itération sur chaque tbl_regression dans la liste courante
    exposure_name <- names(tbls_by_outcome_log[[i]])[j]                         # Nom de la variable d'exposition pour cette itération
    tbl_data <- tbls_by_outcome_log[[i]][[j]] %>%                               # Extraction des données du tableau tbl_regression
      as_tibble() %>%
      mutate(Outcome = outcome_name, Exposure = exposure_name)
    
    table_log <- bind_rows(table_log, tbl_data)                                 # Ajout des données extraites au tibble final
  }
}
rm(terms, formula, model, 
   tbl, tbl_data, tbls_for_this_outcome_log, 
   stacked_tbl_log, tbls_for_outcome_log,
   stacked_tbls_by_outcome_log,
   exposure, exposure_name, i, j, outcome, outcome_name)

table_log <-                                                                    # Ajout variable de la correspondance en phyla
  left_join(table_log, corres, by = "Outcome") %>% 
  select(Phyla_corres, everything())

table_log <- table_log %>%
  select(Phyla_corres, 
         Outcome, 
         Pollutants = Exposure, 
         Beta = "**Beta**", 
         "95% CI" = "**95% CI**",
         "p-value" = "**p-value**", 
         "Characteristic" = "**Characteristic**") %>%
  mutate(
    Phyla_corres = as.factor(Phyla_corres), 
    Phyla_corres = fct_relevel(Phyla_corres,
                               "Firmicutes", "Actinobacteria", 
                               "Bacteroidetes", "Proteobacteria", 
                               "Verrucomicrobia", "Candidatus_Saccharibacteria"),
    Time_window = case_when(grepl("t2", Pollutants) ~ "Mother, pregnancy trim. 2", 
                            grepl("t3", Pollutants) ~ "Mother, pregnancy trim. 3",
                            grepl("Y1", Pollutants) ~ "Child, 12 months", 
                            .default = "Mother, pregnancy trim. 2"),
    Pollutants = str_replace_all(Pollutants,
                                 c(
                                   "mo_" = "",
                                   "ch_" = "",
                                   "DEHP" = "ΣDEHP",
                                   "DiNP" = "ΣDiNP",
                                   "DINCH" = "ΣDINCH",
                                   "_ms_i_cor_t2_ln" = "", 
                                   "_ms_i_cor_t3_ln" = "", 
                                   "_ms_i_cor_Y1_ln" = "", 
                                   "_i_cor_t2_ln" = "", 
                                   "_i_cor_t3_ln" = "", 
                                   "_i_cor_Y1_ln" = ""
                                   
                                 )),
    Pollutants_Time_window = case_when(Time_window == "Mother, pregnancy trim. 2" ~ paste(Pollutants, "trim.2", sep = " "), 
                                       Time_window == "Mother, pregnancy trim. 3" ~ paste(Pollutants, "trim.3", sep = " "), 
                                       Time_window == "Child, 12 months" ~ paste(Pollutants, "12 months", sep = " ")), 
    Pollutants_Time_window = 
      fct_relevel(Pollutants_Time_window, 
                  "ΣDINCH 12 months", "ΣDINCH trim.3", "ΣDINCH trim.2", "ohMPHP 12 months",
                  "ohMPHP trim.3", "ohMPHP trim.2", "MEP 12 months", "MEP trim.3",
                  "MEP trim.2", "MBzP 12 months", "MBzP trim.3", "MBzP trim.2",
                  "MiBP 12 months", "MiBP trim.3", "MiBP trim.2", "ΣDiNP 12 months",
                  "ΣDiNP trim.3", "ΣDiNP trim.2", "MnBP 12 months", "MnBP trim.3",
                  "MnBP trim.2", "ΣDEHP 12 months", "ΣDEHP trim.3", "ΣDEHP trim.2"),
    Pollutants_Time_window_rec = str_replace_all(Pollutants_Time_window, 
                                                 c("trim.2" = "t2", 
                                                   "trim.3" = "t3", 
                                                   "12 months" = "Y1")), 
    `p-value` = gsub("__", "", `p-value`),
    `p-value` = as.numeric(`p-value`),
    `q-value` = `p-value`/(29*31), 
    p_value_shape = ifelse(`p-value`<0.05, "p-value<0.05", "p-value≥0.05"),
    q_value_shape = ifelse(`q-value`<0.05, "q-value<0.05", "q-value≥0.05"), 
    sens_beta = ifelse(Beta < 0, "Beta<0", "Beta≥0"), 
    sens_beta = fct_relevel(sens_beta, "Beta≥0", "Beta<0"))  %>% 
  separate(col = "95% CI", into = c("lower_CI", "upper_CI"), sep = ",", remove = FALSE) %>%
  mutate(
    lower_CI = as.numeric(lower_CI),
    upper_CI = as.numeric(upper_CI)
  ) %>%
  select(
    Phyla_corres,
    Outcome, 
    Pollutants, 
    Time_window, 
    Pollutants_Time_window, Pollutants_Time_window_rec, 
    Beta, sens_beta, 
    "95% CI", lower_CI, upper_CI, 
    "p-value", p_value_shape, 
    "q-value", q_value_shape)

table_log$Outcome_rec <- table_log$Outcome %>%
  fct_recode(
    "Clostridium IV" = "Clostridium_IV",
    "Clostridium sensu stricto" = "Clostridium_sensu_stricto",
    "Clostridium XlVa" = "Clostridium_XlVa",
    "Clostridium XVIII" = "Clostridium_XVIII",
    "Erysipelotrichaceae incertae sedis" = "Erysipelotrichaceae_incertae_sedis",
    "Escherichia and Shigella" = "Escherichia_Shigella",
    "Lachnospiracea incertae sedis" = "Lachnospiracea_incertae_sedis",
    "Ruminococcus 2" = "Ruminococcus2",
    "Saccharibacteria genera incertae sedis" = "Saccharibacteria_genera_incertae_sedis"
  )



# Figures ----
### Mahatan plot Fig.4 ----
mahatan_plot <- table_log  %>%
  mutate(
    Outcome_rec = 
      fct_relevel(Outcome_rec, 
                  "Alistipes", "Anaerostipes", "Bacteroides", "Anaerotruncus",
                  "Bifidobacterium", "Blautia", "Escherichia and Shigella", "Butyricicoccus", 
                  "Cellulosibacter", "Clostridium IV", "Clostridium sensu stricto",
                  "Clostridium XlVa", "Clostridium XVIII", "Collinsella", "Coprococcus",
                  "Dialister", "Dorea", "Eggerthella", "Eisenbergiella", "Enterobacter",
                  "Enterococcus", "Erysipelotrichaceae incertae sedis", 
                  "Faecalibacterium", "Flavonifractor", "Fusicatenibacter", "Gemmiger",
                  "Granulicatella", "Haemophilus", "Hungatella", "Intestinibacter",
                  "Klebsiella", "Lachnospiracea incertae sedis", "Lactococcus",
                  "Oscillibacter", "Parabacteroides", "Peptoniphilus", "Romboutsia",
                  "Roseburia", "Ruminococcus", "Ruminococcus 2", "Saccharibacteria genera incertae sedis",
                  "Terrisporobacter","Streptococcus", "Subdoligranulum", "Veillonella", "Akkermansia")) %>%
  ggplot(aes(x = -log10(`p-value`), y = Outcome_rec)) +
  geom_point(aes(shape = sens_beta), size = 2) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = -log10(0.05/(14*33)), linetype = "dashed", color = "blue") +
  theme_lucid() +
  labs(x = "-log10(P-value)", 
       y = "Genera", 
       shape = "") +
  geom_text(aes(label = ifelse(`p-value` < 0.008, as.character(Pollutants_Time_window_rec), "")), hjust = -0.05, vjust = -0.3, angle = 35, size = 3.5) +
  scale_shape_manual(values = c("Beta<0" = 15, "Beta≥0" = 17)) +# 15: carré plein, 17: triangle plein
  theme(
    legend.position = "right",
    legend.box = "vertical", 
    legend.justification = "right", 
    axis.text.y = element_text(face = "italic"))

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
  geom_text(aes(label = ifelse(`p-value` < 0.008, as.character(Pollutants_Time_window_rec), "")), hjust = -0.07, vjust = -0.2, angle = 22.5, size = 3.5) +
  scale_shape_manual(values = c("Beta<0" = 15, "Beta≥0" = 17)) +# 15: carré plein, 17: triangle plein
  theme(
    legend.position = "right",
    legend.box = "vertical", 
    legend.justification = "right", 
    axis.text.y = element_text(face = "italic"), 
    plot.margin = margin(t = 30))

mahatan_plot
ggsave("4_output/manhattan_plot.tiff", 
       mahatan_plot, 
       device = "tiff",
       units = "cm",
       dpi = 300,
       height = 20, 
       width = 44)

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
