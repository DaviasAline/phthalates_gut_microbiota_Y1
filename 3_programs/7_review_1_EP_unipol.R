# Revisions EP - article phthalates - analyses univariées
# Aline Davias 
# 27.08.2024


# Chargement des packages ----
library(tidyverse)
library(haven)
library(reshape2)
library(GGally)
library(gtsummary)
library(summarytools)
library(patchwork)
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
library(rstatix)
library(grDevices)
library(lazyeval)
library(mice)
library(car)
library(expss)
library(phyloseq)
library(vegan)
library(SRS)
library(see)
library(corrplot)
library(bkmr)
library(fields)
library(future)
library(future.apply)
library(writexl)
library(ggrepel)

# Chargement des données ----
load("1_intermediate_data/2_data_selection_AD_gumme.RData")
taxa_table <- read_csv("0_source_data/taxa_table_ASVbased_Y1_AD_20220504_8.csv")
corres <- 
  taxa_table %>% 
  select(Phyla_corres = ch_feces_phylum_ASVbased_Y1, 
         Class_corres = ch_feces_class_ASVbased_Y1, 
         Order_corres = ch_feces_order_ASVbased_Y1, 
         Family_corres = ch_feces_family_ASVbased_Y1, 
         Outcome = ch_feces_genus_ASVbased_Y1) %>%
  filter(Outcome %in% genera_vec) %>%
  distinct(Outcome, .keep_all = TRUE)
rm(taxa_table)

asv_raw_not_rarefied <-              # pour les analyses betadiv
  read_labelled_csv(
    "0_source_data/gut_microbiota_ASVbased_Y1_labelled_AD_20220504_7239.csv") 

# Création des fonctions ----
source("3_programs/3_functions_AD_gumme.R", echo=TRUE)
rm(model_covar, model_multi, model_summary, model_univ_multi)

# Fonction pour générer un tableau de régression
generate_tbl <- function(outcome, var, covariates, data) {
  formula <- as.formula(paste(outcome, "~", var, "+", paste(covariates, collapse = "+")))
  model <- lm(formula, data = data)
  tbl <- 
    tbl_regression(
      model, 
      pvalue_fun = custom_pvalue_fun,
      estimate_fun = scales::label_number(accuracy = .01, decimal.mark = ".")) %>%
    bold_p() %>%
    bold_labels() %>%
    add_global_p(singular.ok = TRUE, keep = TRUE)
  return(tbl)
}

# Fonction pour exécuter la régression linéaire (analyses principales)
main_model <- function(outcome, data) {
  pre_tbls <- map(phthalates_pre, ~ generate_tbl(outcome, .x, covariates_pre, data))
  post_tbls <- map(phthalates_post, ~ generate_tbl(outcome, .x, covariates_post, data))

  names(pre_tbls) <- phthalates_pre          # Combiner les résultats pré et post
  names(post_tbls) <- phthalates_post
  c(pre_tbls, post_tbls)
}

# Fonction pour exécuter la régression linéaire (analyses de sensibilité gravité spécifique)
sensi_12_a <- function(outcome, exposure, covariates) {
  formula_str <- paste(outcome, "~", exposure, "+", paste(covariates, collapse = "+"))
  model <- lm(as.formula(formula_str), data = bdd)
  summary_model <- summary(model)
  
  # Extraction des résultats pour l'exposition principale seulement
  beta <- summary_model$coefficients[exposure, "Estimate"]
  conf_int <- confint(model, level = 0.95)[exposure, ]
  p_value <- summary_model$coefficients[exposure, "Pr(>|t|)"]
  
  return(list(beta_sg_std = beta, conf.low_sg_std = conf_int[1], conf.high_sg_std = conf_int[2], p.value_sg_std = p_value))
}

sensi_12_b <- function(outcome, exposure, covariates) {
  formula_str <- paste(outcome, "~", exposure, "+", paste(covariates, collapse = "+"))
  model <- lm(as.formula(formula_str), data = bdd)
  summary_model <- summary(model)
  
  beta <- summary_model$coefficients[exposure, "Estimate"]                      # Extraction des résultats pour l'exposition principale seulement
  conf_int <- confint(model, level = 0.95)[exposure, ]
  p_value <- summary_model$coefficients[exposure, "Pr(>|t|)"]
  
  return(list(beta_sg_adj = beta, 
              conf.low_sg_adj = conf_int[1], 
              conf.high_sg_adj = conf_int[2], 
              p.value_sg_adj = p_value))
}

# Fonction pour exécuter la régression linéaire (analyses de sensibilité hospitalisation)
sensi_13 <- function(outcome, exposure, covariates) {
  formula_str <- paste(outcome, "~", exposure, "+", paste(covariates, collapse = "+"))
  model <- lm(as.formula(formula_str), data = bdd)
  summary_model <- summary(model)
  
  # Extraction des résultats pour l'exposition principale seulement
  beta <- summary_model$coefficients[exposure, "Estimate"]
  conf_int <- confint(model, level = 0.95)[exposure, ]
  p_value <- summary_model$coefficients[exposure, "Pr(>|t|)"]
  
  return(list(beta_hospit = beta, conf.low_hospit = conf_int[1], conf.high_hospit = conf_int[2], p.value_hospit = p_value))
}

# Fonction pour exécuter la régression linéaire (analyses de sensibilité age gestationel, poids et taille à la naissance et à un an)
sensi_14 <- function(outcome, exposure, covariates) {
  formula_str <- paste(outcome, "~", exposure, "+", paste(covariates, collapse = "+"))
  model <- lm(as.formula(formula_str), data = bdd)
  summary_model <- summary(model)
  
  # Extraction des résultats pour l'exposition principale seulement
  beta <- summary_model$coefficients[exposure, "Estimate"]
  conf_int <- confint(model, level = 0.95)[exposure, ]
  p_value <- summary_model$coefficients[exposure, "Pr(>|t|)"]
  
  return(list(beta_w_he = beta, conf.low_w_he = conf_int[1], conf.high_w_he = conf_int[2], p.value_w_he = p_value))
}

# Fonction pour extraire les résultats liés à la variable explicative d'intérêt
extract_exposure_results <- function(tbl, exposure_var) {
  tbl %>%
    as_tibble() %>%
    filter(term == exposure_var) %>%
    select(beta = estimate, conf.low, conf.high, p.value)
}

# Fonction pour traiter une outcome et extraire les résultats
process_outcome_results <- function(outcome, exposure_list) {
  map_df(names(exposure_list), function(exposure_var) {
    tbl <- exposure_list[[exposure_var]]$table_body
    filtered_results <- extract_exposure_results(tbl, exposure_var)
    filtered_results %>%
      mutate(Outcome = outcome, Exposure = exposure_var)
  })
}

# Fonctions pour mise en forme tableaux articles 
filter_exposure_tbl <- function(tbl, exposure_var) {
  tbl %>%
    modify_table_body(~ .x %>% dplyr::filter(variable == exposure_var))
}

create_filtered_tbl_merge <- function(exposure_var, outcome_list, spanner_names) {
  tbls <- map(outcome_list, ~ filter_exposure_tbl(.x[[exposure_var]], exposure_var))  # Filtrer chaque tbl_regression
  tbl_merge(tbls, tab_spanner = spanner_names)  # Fusionner les tbl_regressions filtrés
}

create_tbl_for_range <- function(range, spanner_names) {
  map(phthalates, ~ create_filtered_tbl_merge(.x, results_main_list[range], spanner_names)) %>%
    tbl_stack()
}

# Fonction pour adapter les décimales des p-values
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

format_p_value <- function(p) {
  if (p >= 0.01) {
    # Si p est supérieur ou égal à 0.01, on garde 2 décimales
    formatted_p <- sprintf("%.2f", p)
  } else {
    # Si p est inférieur à 0.01, on augmente le nombre de décimales jusqu'à trouver un chiffre significatif
    formatted_p <- sprintf("%.1g", p)
  }
  return(formatted_p)
}

# Analyses principales ----
## Analyses univarivées (alpha diversité et taxonomie) ----
results_main_list <- map(outcomes, ~ main_model(.x, bdd))                       # Création de la liste contenant tous les tbl_reg
names(results_main_list) <- outcomes                                            # Ajout des noms des outcomes à chaque sous liste

results_main <- map_df(names(results_main_list), function(outcome) {            # Création tableau brut des resultats (pour figures)
  process_outcome_results(outcome, results_main_list[[outcome]])
})

results_main <-                                                                 # Ajout variable de la correspondance en phyla
  left_join(results_main, corres, by = "Outcome") %>% 
  select(Phyla_corres, Class_corres, Order_corres, Family_corres, everything())

results_main <- results_main %>%                                                # Réorganisation tableau brut des resultats (pour figures)                   
  arrange(Outcome, Exposure) %>%
  rename(Exposure_window = Exposure) %>%
  mutate(
    Exposure_window = factor(Exposure_window, levels = phthalates), 
    Outcome = factor(Outcome, levels = outcomes), 
    Exposure_window = as.factor(Exposure_window),
    Exposure_window = str_replace_all(Exposure_window,
                               c("mo_" = "",
                                 "ch_" = "",
                                 "_i_cor_" = " ", 
                                 "_ln" = "",
                                 "_conj" = "",
                                 "ln" = "",
                                 "_cor" = "", 
                                 "_t2" = " t2", 
                                 "_t3" = " t3", 
                                 "_Y1" = " Y1", 
                                 "_ms" = "")), 
    Outcome_name = str_replace_all(Outcome, 
                                   c("ch_feces_SpecRich_5000_ASV_Y1" = "Specific richness", 
                                     "ch_feces_Shannon_5000_ASV_Y1" = "Shannon diversity", 
                                     "ch_feces_rel_p1_Y1" = "Firmicutes", 
                                     "ch_feces_rel_p2_Y1" = "Actinobacteria", 
                                     "ch_feces_rel_p3_Y1" = "Bacteroidetes", 
                                     "ch_feces_rel_p4_Y1" = "Proteobacteria")),
    Window = case_when(grepl("t2", Exposure_window) ~ "Trim.2", 
                                grepl("t3", Exposure_window) ~ "Trim.3", 
                                grepl("Y1", Exposure_window) ~ "12 months", 
                                TRUE ~ "Trim.2"), 
    Window = fct_relevel(Window, 
                                  "Trim.2", 
                                  "Trim.3", 
                                  "12 months"), 
    FWER.p.value_alpha = ifelse(Outcome_name %in% c("Specific richness", "Shannon diversity"), 
                                p.value * 14 * 3, NA),
    FWER.p.value_alpha = ifelse(FWER.p.value_alpha > 1, ">0.99", FWER.p.value_alpha),
    FWER.p.value_shape_alpha = case_when(p.value< 0.0012 & 
                                           Outcome_name %in% c("Specific richness", 
                                                          "Shannon diversity")~ "p.value <0.0012",
                                         p.value > 0.0012 & 
                                           Outcome_name %in% c("Specific richness", 
                                                          "Shannon diversity")~ "p.value >0.0012"), 
    FWER.p.value_shape_alpha = fct_relevel(FWER.p.value_shape_alpha,
                                           "p.value >0.0012", 
                                           "p.value <0.0012"), 
    Outcome_name = fct_recode(
      Outcome_name, 
        "Clostridium IV" = "Clostridium_IV",
        "Clostridium sensu stricto" = "Clostridium_sensu_stricto",
        "Clostridium XlVa" = "Clostridium_XlVa",
        "Clostridium XVIII" = "Clostridium_XVIII",
        "Erysipelotrichaceae incertae sedis" = "Erysipelotrichaceae_incertae_sedis",
        "Escherichia and Shigella" = "Escherichia_Shigella",
        "Lachnospiracea incertae sedis" = "Lachnospiracea_incertae_sedis",
        "Ruminococcus 2" = "Ruminococcus2",
        "Saccharibacteria genera incertae sedis" = "Saccharibacteria_genera_incertae_sedis"),
    Phyla_corres = as.factor(Phyla_corres), 
    Phyla_corres = fct_relevel(Phyla_corres,
                               "Firmicutes", "Actinobacteria", 
                               "Bacteroidetes", "Proteobacteria", 
                               "Verrucomicrobia", "Candidatus_Saccharibacteria"),
    Exposure_window = str_replace_all(Exposure_window,
                               c("DEHP" = "ΣDEHP",
                                 "DiNP" = "ΣDiNP",
                                 "DINCH" = "ΣDINCH")),
    Exposure = str_replace_all(Exposure_window,
                                   c(" t2" = "",
                                     " t3" = "",
                                     " Y1" = "")),
    Exposure_rec = gsub("T", "t", Window),
    Exposure_window_rec = paste(Exposure, Exposure_rec, sep = " "),
    Exposure_window_rec = 
      fct_relevel(Exposure_window_rec, 
                  "ΣDINCH 12 months", "ΣDINCH trim.3", "ΣDINCH trim.2", "ohMPHP 12 months",
                  "ohMPHP trim.3", "ohMPHP trim.2", "MEP 12 months", "MEP trim.3",
                  "MEP trim.2", "MBzP 12 months", "MBzP trim.3", "MBzP trim.2",
                  "MiBP 12 months", "MiBP trim.3", "MiBP trim.2", "ΣDiNP 12 months",
                  "ΣDiNP trim.3", "ΣDiNP trim.2", "MnBP 12 months", "MnBP trim.3",
                  "MnBP trim.2", "ΣDEHP 12 months", "ΣDEHP trim.3", "ΣDEHP trim.2"),
    sens_beta = ifelse(beta < 0, "Beta<0", "Beta≥0"), 
    sens_beta = fct_relevel(sens_beta, "Beta≥0", "Beta<0"))  %>% 
  select(
    Phyla_corres, Class_corres, Order_corres, Family_corres, 
    Outcome, Outcome_name, 
    Exposure, Exposure_rec,
    Window,
    Exposure_window, Exposure_window_rec,
    beta, conf.low, conf.high, p.value, 
    sens_beta,
    FWER.p.value_alpha, FWER.p.value_shape_alpha)


## Analyses univariées (beta diversité) ----
### Rarefaction 
#### Data preparation 
asv_raw_not_rarefied <- 
  asv_raw_not_rarefied %>%
  select(ident, 
         starts_with("ch_feces_raw_asv")) %>%
  na.omit()
row.names(asv_raw_not_rarefied) <- NULL 
asv_raw_not_rarefied <- asv_raw_not_rarefied %>%
  column_to_rownames("ident") %>%
  t() %>%
  as.data.frame()%>% 
  rownames_to_column("ch_feces_ASV_ID_Y1") %>%
  mutate(
    ch_feces_ASV_ID_Y1 = str_replace_all(ch_feces_ASV_ID_Y1, c("ch_feces_raw_" = "", "_Y1"=""))) %>%
  column_to_rownames("ch_feces_ASV_ID_Y1") %>%
  otu_table(taxa_are_rows = TRUE)     

#### from the raw ASV table, choose a threshold for the sequencing depth 
#### define a subset of samples to keep = samples with a sequencing depth > the chosen threshold  
#### samples with a sequencing depth < the chosen threshold become missing data 
keep_5000 <- 
  names(which(sample_sums(asv_raw_not_rarefied)>= 5000)) %>%
  prune_samples(asv_raw_not_rarefied) %>% 
  as.data.frame()                # Loss of 6 samples 

#### reduce the sequencing depth of the samples to the chosen threshold
#### the sequences kept within each sample are randomly selected 
ASV_rarefied_5000_Y1 <- keep_5000 %>%
  SRS(5000, set_seed = TRUE, seed = 1)
rownames(ASV_rarefied_5000_Y1)<- rownames(keep_5000)

#### put the dataframe with rarefied ASVs in columns and samples in rows
ASV_rarefied_5000_Y1 <- 
  ASV_rarefied_5000_Y1 %>%
  t() %>%
  as.data.frame()
rm(keep_5000)

#### check if the rarefaction worked properly 
#### ok the rowsums are equal to 5000, the threshold we chose
rowSums(ASV_rarefied_5000_Y1)   

### Metric calculation (Bray Curtis) 
set.seed(1996)
metric_bray_curtis <- vegan::vegdist(ASV_rarefied_5000_Y1, method = "bray")
### warning message: Plus d’une classe "dist" est trouvée en cache : Utilisation de la première, depuis l’espace de noms 'BiocGenerics'. Aussi défini par ‘spam’
### ok, phyloseq package is supposed to use BiocGeneric 

metric_bray_curtis <- 
  metric_bray_curtis %>% 
  as.matrix %>% 
  as.data.frame() %>%
  rownames_to_column(var = "ident")

### Data preparation
bdd_betadiv <- bdd %>%
  select(ident, 
         all_of(covariates_post), 
         all_of(phthalates)) %>%
  mutate(ident = as.character(ident))

for (var in phthalates) {
  new_var_name <- str_replace_all(var, "_ln", "_ter")                          # Création d'expo tertiles
  bdd_betadiv[[new_var_name]] <- 
    cut(bdd_betadiv[[var]],
        breaks = quantile(bdd_betadiv[[var]], 
                          probs = seq(0, 1, by = 1/3), 
                          na.rm = TRUE), 
        include.lowest = TRUE, 
        labels = c("1st tertile", "2nd tertile", "3rd tertile"))
}
rm(var, new_var_name)

phthalates_ter <- str_replace_all(phthalates, "_ln", "_ter")
phthalates_pre_ter <- str_replace_all(phthalates_pre, "_ln", "_ter")
phthalates_post_ter <- str_replace_all(phthalates_post, "_ln", "_ter")

#### merge betadiversity data and metadata
bdd_betadiv <- inner_join(metric_bray_curtis, bdd_betadiv, by = "ident")
bdd_betadiv[phthalates_ter] <- lapply(bdd_betadiv[phthalates_ter], as.factor)
covariates_pre_beta <- bdd_betadiv %>% select(all_of(covariates_pre)) %>% select(-mo_interpreg_3cat) %>% colnames()
covariates_post_beta <- bdd_betadiv %>% select(all_of(covariates_post)) %>% select(-mo_interpreg_3cat) %>% colnames()

### Trim 2. ----
#### filter t2 (because of the NA on the pollutants)
bdd_betadiv %>% select(ident, contains("t2")) %>% filter_all(any_vars(is.na(.)))
bdd_betadiv_t2 <- bdd_betadiv %>% 
  select(ident, contains("t2"), everything()) %>%
  filter(ident != 15804) %>%
  select(-"15804")
all_dist_t2 <- bdd_betadiv_t2 %>%
  select(all_of(.[["ident"]])) %>%
  as.dist()

explanatory_vars_t2 <- 
  bdd_betadiv_t2 %>% 
  select(all_of(phthalates_ter)) %>%
  select(contains("t2")) %>%
  colnames()

results_betadiv_multivar_bray_curtis_t2 <-
  lapply(explanatory_vars_t2, function(x) {
    formula <- reformulate(c(x, covariates_pre_beta), response = "all_dist_t2")
    adonis2(formula, data = bdd_betadiv_t2, permutations = 999)
  })
results_betadiv_multivar_bray_curtis_t2 <- 
  do.call(rbind, results_betadiv_multivar_bray_curtis_t2) %>%
  rownames_to_column(var = "Explanatory variables") %>% 
  mutate(
    Pollutants = c(rep("mo_DEHP_ms_i_cor_t2_ter", times = 18), 
                   rep("mo_MnBP_i_cor_t2_ter", times = 18), 
                   rep("mo_DiNP_ms_i_cor_t2_ter", times = 18),
                   rep("mo_MiBP_i_cor_t2_ter", times = 18),
                   rep("mo_MBzP_i_cor_t2_ter", times = 18),
                   rep("mo_MEP_i_cor_t2_ter", times = 18),
                   rep("mo_ohMPHP_i_cor_t2_ter", times = 18),
                   rep("mo_DINCH_ms_i_cor_t2_ter", times = 18))) %>%
  select(Pollutants, everything())


results_betadiv_multivar_bray_curtis_t2 %>%                      # Visualisation des résultats significatifs multivarié
  filter(`Pr(>F)` < 0.05) %>%
  filter(`Explanatory variables` %in% c(explanatory_vars_t2)) %>%
  View()


### Trim 3. ----
#### filter t3 (because of the NA on the pollutants)
bdd_betadiv %>% select(ident, contains("t3")) %>% filter_all(any_vars(is.na(.)))
bdd_betadiv_t3 <- bdd_betadiv %>% 
  select(ident, contains("t3"), everything()) %>%
  filter(!ident %in% c(17827, 15929, 26891, 25668, 23330, 28199)) %>%
  select(-"15929", -"17827", -"23330", -"25668", -"26891", -"28199")
all_dist_t3 <- bdd_betadiv_t3 %>%
  select(all_of(.[["ident"]])) %>%
  as.dist()

explanatory_vars_t3 <- 
  bdd_betadiv_t3 %>% 
  select(all_of(phthalates_ter)) %>%
  select(contains("t3")) %>%
  colnames()


results_betadiv_multivar_bray_curtis_t3 <-
  lapply(explanatory_vars_t3, function(x) {
    formula <- reformulate(c(x, covariates_pre_beta), response = "all_dist_t3")
    adonis2(formula, data = bdd_betadiv_t3, permutations = 999)
  })
results_betadiv_multivar_bray_curtis_t3 <-
  do.call(rbind, results_betadiv_multivar_bray_curtis_t3) %>%
  rownames_to_column(var = "Explanatory variables") %>%
  mutate(
    Pollutants = c(rep("mo_DEHP_ms_i_cor_t3_ter", times = 18), 
                   rep("mo_MnBP_i_cor_t3_ter", times = 18), 
                   rep("mo_DiNP_ms_i_cor_t3_ter", times = 18),
                   rep("mo_MiBP_i_cor_t3_ter", times = 18),
                   rep("mo_MBzP_i_cor_t3_ter", times = 18),
                   rep("mo_MEP_i_cor_t3_ter", times = 18),
                   rep("mo_ohMPHP_i_cor_t3_ter", times = 18),
                   rep("mo_DINCH_ms_i_cor_t3_ter", times = 18))) %>%
  select(Pollutants, everything())

# results_betadiv_univar_bray_curtis_t3 %>%                         # Visualisation des résultats significatifs univarié
#   filter(`Pr(>F)` < 0.05) %>%
#   filter(`Explanatory variables` %in% c(explanatory_vars_t3)) %>%
#   View()

results_betadiv_multivar_bray_curtis_t3 %>%                     # Visualisation des résultats significatifs multivarié
  filter(`Pr(>F)` < 0.05) %>%
  filter(`Explanatory variables` %in% c(explanatory_vars_t3)) %>%
  View()


### Y1 ----
#### filter Y1 (because of the NA on the pollutants)
bdd_betadiv %>% select(ident, contains("Y1")) %>% filter_all(any_vars(is.na(.)))
bdd_betadiv_Y1 <- bdd_betadiv %>% 
  select(ident, contains("Y1"), everything()) %>%
  filter(!ident %in% c(23994, 25166, 26766, 14668, 26923)) %>%
  select(-"23994", -"25166", -"26766", -"14668", -"26923")
all_dist_Y1 <- bdd_betadiv_Y1 %>%
  select(all_of(.[["ident"]])) %>%
  as.dist()

explanatory_vars_Y1 <- 
  bdd_betadiv_Y1 %>% 
  select(all_of(phthalates_ter))  %>%
  select(contains("Y1")) %>%
  colnames()


results_betadiv_multivar_bray_curtis_Y1 <-
  lapply(explanatory_vars_Y1, function(x) {
    formula <- reformulate(c(x, covariates_post_beta), response = "all_dist_Y1")
    adonis2(formula, data = bdd_betadiv_Y1, permutations = 999)
  })
results_betadiv_multivar_bray_curtis_Y1 <-
  do.call(rbind, results_betadiv_multivar_bray_curtis_Y1) %>%
  rownames_to_column(var = "Explanatory variables") %>%
  mutate(
    Pollutants = c(rep("ch_DEHP_ms_i_cor_Y1_ter", times = 21), 
                   rep("ch_MnBP_i_cor_Y1_ter", times = 21), 
                   rep("ch_DiNP_ms_i_cor_Y1_ter", times = 21),
                   rep("ch_MiBP_i_cor_Y1_ter", times = 21),
                   rep("ch_MBzP_i_cor_Y1_ter", times = 21),
                   rep("ch_MEP_i_cor_Y1_ter", times = 21),
                   rep("ch_ohMPHP_i_cor_Y1_ter", times = 21),
                   rep("ch_DINCH_ms_i_cor_Y1_ter", times = 21))) %>%
  select(Pollutants, everything())


results_betadiv_multivar_bray_curtis_Y1 %>%         # Visualisation des résultats significatifs multivarié
  filter(`Pr(>F)` < 0.05) %>%
  filter(`Explanatory variables` %in% c(explanatory_vars_Y1)) %>%
  View()

results_betadiv <- list(results_betadiv_multivar_bray_curtis_t2, 
                        results_betadiv_multivar_bray_curtis_t3, 
                        results_betadiv_multivar_bray_curtis_Y1)

rm(ASV_rarefied_5000_Y1, metric_bray_curtis, 
   phthalates_pre_ter, phthalates_post_ter, 
   covariates_pre_beta, covariates_post_beta, 
   bdd_betadiv,
   bdd_betadiv_t2, bdd_betadiv_t3, bdd_betadiv_Y1, 
   all_dist_t2,all_dist_t3, all_dist_Y1, 
   explanatory_vars_t2, explanatory_vars_t3, explanatory_vars_Y1, 
   results_betadiv_multivar_bray_curtis_t2, 
   results_betadiv_multivar_bray_curtis_t3, 
   results_betadiv_multivar_bray_curtis_Y1)

# Article ----
## Table 1 : alpha div - univar ----
Table_1 <- results_main %>% 
  mutate(
    Exposure = factor(Exposure, levels = c("ΣDEHP", "MnBP", "ΣDiNP", "MiBP", "MBzP", "MEP", "ohMPHP", "ΣDINCH")), 
    Window = factor(Window, levels = c("Trim.2", "Trim.3", "12 months")), 
    p.value = case_when(p.value < 0.001 ~ "<0.001",
                        p.value < 0.01 ~ format(round(p.value, 3), nsmall = 3),
                        p.value > 0.01 ~ format(round(p.value, 2), nsmall = 2)), 
    beta = case_when(Outcome_name == "Specific richness" ~ format(round(beta, 1), nsmall = 1), 
                     Outcome_name == "Shannon diversity" ~format(round(beta, 2), nsmall  = 2)),
    conf.low = case_when(Outcome_name == "Specific richness" ~ format(round(conf.low, 1), nsmall = 1), 
                         Outcome_name == "Shannon diversity" ~format(round(conf.low, 2), nsmall  = 2)),
    conf.high = case_when(Outcome_name == "Specific richness" ~ format(round(conf.high, 1), nsmall = 1), 
                          Outcome_name == "Shannon diversity" ~format(round(conf.high, 2), nsmall  = 2)),
    CI = paste(conf.low, conf.high, sep = ", ")) %>%
  arrange(Outcome, Exposure, Window) %>%
  select(Outcome_name, Exposure_window_rec, beta, CI, p.value) %>% 
  filter(Outcome_name %in% c("Specific richness", "Shannon diversity")) %>% 
  pivot_wider(names_from = Outcome_name, values_from = c("beta", "CI", "p.value")) %>%
  select(Exposure_window_rec, contains("richness"), contains("Shannon")) 

write_xlsx(Table_1, "4_output/review/Table_1.xlsx")

## Table 2 : phyla - univar ----
Table_2 <- results_main %>% 
  mutate(
    Exposure = factor(Exposure, levels = c("ΣDEHP", "MnBP", "ΣDiNP", "MiBP", "MBzP", "MEP", "ohMPHP", "ΣDINCH")), 
    Window = factor(Window, levels = c("Trim.2", "Trim.3", "12 months")), 
    beta = format(round(beta, 1), nsmall = 1), 
    conf.low = format(round(conf.low, 1), nsmall = 1), 
    conf.high = format(round(conf.high, 1), nsmall = 1),
    p.value = case_when(p.value < 0.001 ~ "<0.001",
                        p.value < 0.01 ~ format(round(p.value, 3), nsmall = 3),
                        p.value > 0.01 ~ format(round(p.value, 2), nsmall = 2)), 
    CI = paste(conf.low, conf.high, sep = ", ")) %>%
  arrange(Outcome, Exposure, Window) %>%
  select(Outcome_name, Exposure_window_rec, beta, CI, p.value) %>% 
  filter(Outcome_name %in% c("Firmicutes", "Actinobacteria", "Bacteroidetes", "Proteobacteria")) %>% 
  pivot_wider(names_from = Outcome_name, values_from = c("beta", "CI", "p.value")) %>%
  select(Exposure_window_rec, 
         contains("Firmicutes"), 
         contains("Actinobacteria"), 
         contains("Bacteroidetes"), 
         contains("Proteobacteria")) %>%
  filter(Exposure_window_rec %in% c("ΣDEHP 12 months", 
                         "MEP 12 months", 
                         "ohMPHP 12 months")) 
write_xlsx(Table_2, "4_output/review/Table_2.xlsx")

## Figure 1 : alpha div - univar ----
forestplot_rich <- function(results_main, Outcome_name) {
  results <- results_main %>%
    ggplot(aes(x = Exposure_window,
               y = beta,
               min = conf.low,
               ymax = conf.high,
               #color = interaction(exposure_window, term_2),
               #color = term_rec,
               color = FWER.p.value_shape)) +
    geom_hline(yintercept = 0, linetype="dashed") +
    geom_pointrange(position = position_dodge(width = 0.7), size = 0.4,
                    aes(color = ifelse(FWER.p.value_shape_alpha == "p.value >0.0012", "p.value >0.0012", "p.value <0.0012"))) +
    labs(x = "Exposures", y = Outcome_name) +
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

forestplot_shannon <- function(results_main, Outcome_name) {
  results <- results_main %>%
    ggplot(aes(x = Exposure_window,
               y = beta,
               min = conf.low,
               ymax = conf.high,
               #color = interaction(exposure_window, term_2),
               #color = term_rec,
               color = FWER.p.value_shape)) +
    geom_hline(yintercept = 0, linetype="dashed") +
    geom_pointrange(position = position_dodge(width = 0.7), size = 0.4,
                    aes(color = ifelse(FWER.p.value_shape_alpha == "p.value >0.0012", "p.value >0.0012", "p.value <0.0012"))) +
    labs(x = "Exposures", y = Outcome_name) +
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


leg <- results_main %>%
  mutate(
    Exposure_window =  fct_relevel(
      Exposure_window, 
      "ΣDINCH Y1", "ΣDINCH t3", "ΣDINCH t2", 
      "ohMPHP Y1", "ohMPHP t3", "ohMPHP t2", 
      "MEP Y1", "MEP t2", "MEP t3", 
      "MBzP Y1", "MBzP t3", "MBzP t2", 
      "MiBP Y1", "MiBP t3", "MiBP t2", 
      "ΣDiNP Y1", "ΣDiNP t3", "ΣDiNP t2", 
      "MnBP Y1", "MnBP t3", "MnBP t2", 
      "ΣDEHP Y1", "ΣDEHP t3", "ΣDEHP t2")) %>%
  filter(Outcome_name == "Shannon diversity") %>%
  forestplot_shannon(Outcome_name = "Shannon diversity") +
  theme(axis.text.y = element_blank(), 
        axis.title.y = element_blank())
leg <- get_legend(leg) %>% as_ggplot()

forestplot_alpha_1 <-
  results_main %>%
  mutate(
    Exposure_window =  fct_relevel(
      Exposure_window, 
      "ΣDINCH Y1", "ΣDINCH t3", "ΣDINCH t2", 
      "ohMPHP Y1", "ohMPHP t3", "ohMPHP t2", 
      "MEP Y1", "MEP t2", "MEP t3", 
      "MBzP Y1", "MBzP t3", "MBzP t2", 
      "MiBP Y1", "MiBP t3", "MiBP t2", 
      "ΣDiNP Y1", "ΣDiNP t3", "ΣDiNP t2", 
      "MnBP Y1", "MnBP t3", "MnBP t2", 
      "ΣDEHP Y1", "ΣDEHP t3", "ΣDEHP t2")) %>%
  filter(Outcome_name == "Specific richness") %>%
  forestplot_rich(Outcome_name = "Specific richness") +
  theme(legend.position = "none")

forestplot_alpha_2 <-
  results_main %>%
  mutate(
    Exposure_window =  fct_relevel(
      Exposure_window, 
      "ΣDINCH Y1", "ΣDINCH t3", "ΣDINCH t2", 
      "ohMPHP Y1", "ohMPHP t3", "ohMPHP t2", 
      "MEP Y1", "MEP t2", "MEP t3", 
      "MBzP Y1", "MBzP t3", "MBzP t2", 
      "MiBP Y1", "MiBP t3", "MiBP t2", 
      "ΣDiNP Y1", "ΣDiNP t3", "ΣDiNP t2", 
      "MnBP Y1", "MnBP t3", "MnBP t2", 
      "ΣDEHP Y1", "ΣDEHP t3", "ΣDEHP t2")) %>%
  filter(Outcome_name == "Shannon diversity") %>%
  forestplot_shannon(Outcome_name = "Shannon diversity") +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank())+
  theme(legend.position = "none")


Figure_1 <- 
  (forestplot_alpha_1 + forestplot_alpha_2) / leg + 
  plot_layout(heights = c(14, 1))

rm(forestplot_alpha_1, forestplot_alpha_2, leg, 
   forestplot_rich, forestplot_shannon)

ggsave("4_output/review/Figure 1 (forestplot_alpha_phthalates).tiff", 
       plot = Figure_1, 
       device = "tiff",
       units = "mm",
       width = 180, 
       height = 160,
       dpi = 300,
       limitsize = FALSE)

## Figure 2 : alpha div - multivar ----

## Figure 3 : phyla - multivar ----

## Figure 4 : genera - univar ----
Figure_4 <- results_main  %>%
  filter(Outcome %in% genera_vec) %>% 
  mutate(
    Outcome_name = 
      fct_relevel(Outcome_name, 
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
                  "Escherichia and Shigella", "Blautia", "Bacteroides", "Bifidobacterium")) %>%
  ggplot(aes(x = -log10(`p.value`), y = Outcome_name)) +
  geom_point(aes(shape = sens_beta), size = 2) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = -log10(0.05/(14*33)), linetype = "dashed", color = "blue") +
  theme_lucid() +
  labs(x = "-log10(P-value)", 
       y = "Genera", 
       shape = "") +
  geom_text_repel(aes(label = ifelse(`p.value` < 0.008, as.character(Exposure_window), "")), 
                  hjust = -0.07,
                  vjust = -0.2,
                  angle = 42,
                  size = 3.5,
                  #segment.color = NA,   # supprimes les traits qui relient les points aux textes 
                  direction = "y",      # argument pour faire les ajustements de position seulement dans le sens vertical
                  #angle = 20, 
                  seed = 1996,            # seed pour obtenir toujours le même ajustement 
                  #box.padding = 0.5,
                  # hjust = -0.05,        # place le texte un tout petit peu au dessus des points
                  # vjust = 0.5,
                  nudge_x = 0.005, 
                  min.segment.length = 3,  # ne met que les traits de plus de 3
                  max.overlaps = 20
                  ) +
  scale_shape_manual(values = c("Beta<0" = 15, "Beta≥0" = 17)) +# 15: carré plein, 17: triangle plein
  theme(
    legend.position = "bottom",
    legend.box = "vertical", 
    legend.justification = "center", 
    axis.text.y = element_text(face = "italic", size = 12)) +
  xlim(0, 5.8) 

Figure_4
ggsave("4_output/review/Figure_4 (manhatanplot genera).tiff", 
       plot = Figure_4, 
       device = "tiff",
       units = "mm",
       width = 340, 
       height = 200,
       dpi = 300,
       limitsize = FALSE)

# Supplementary analysis ----
## Table S1 : methodo expo ----
# ok

## Table S2 : descript covar mère (cf papier phenols) ----
Table_S2 <- bdd %>% 
  select("mo_age",
         "mo_par_2cat",
         "mo_dipl_3cat",
         "mo_bmi_bepr_3cat", 
         "po_delmod",   
         "mo_interpreg_3cat",
         "mo_tob_gr_anyt_yn_n2",
         "Mo_ETS_anyT_yn1_opt",
         "statut") %>%
  mutate(statut = fct_relevel(statut, "inclu", "exclu")) %>%
  tbl_summary(
    by = statut, 
    type = list(mo_tob_gr_anyt_yn_n2 ~ "categorical", 
                Mo_ETS_anyT_yn1_opt ~ "categorical"), 
    missing = "no") %>%
  add_p() %>%
  bold_labels()

## Table S3 : descript covar enfant (cf papier phenols) ----
Table_S3 <- bdd %>%
  select("ch_sex",
         "po_gd",
         "po_w_kg_3cat",
         "po_he_3cat",
         "ch_food_intro_Y1_3cat",
         "ch_antibio_Y1_2cat",
         "ch_ETS_12m_opt36m",
         "mo_pets",
         "ch_w_Y1_3cat", 
         "ch_he_Y1_3cat", 
         "bf_duration_till48w_4cat_i",
         "statut") %>%
  mutate(statut = fct_relevel(statut, "inclu", "exclu")) %>%
  tbl_summary(
    by = statut, 
    type = list(ch_ETS_12m_opt36m ~ "categorical"), 
    missing = "no") %>%
  add_p() %>%
  bold_labels()

## Table S4 : descrip expo ----
phthalates_S4 <- phthalates %>% str_replace_all("_ln", "")

Table_S4 <- bdd %>%
  select(all_of(phthalates_S4), 
         statut) %>%
  mutate(statut = fct_relevel(statut, "inclu", "exclu")) %>%
  tbl_summary(by = "statut", 
              missing = "no") %>%
  bold_labels() %>%
  add_n()

rm(phthalates_S4)

## Table S5 : descrip microbiote ----
Table_S5 <- bdd %>%
  select(ident, 
         all_of(outcomes))

Table_S5 <- descrip_num(data = bdd, vars = outcomes)
Table_S5 <- Table_S5 %>%
  rename("Min." = Min,
         "1st quartile" = Q1, 
         "3rd quartile" = Q3, 
         "Max." = Max) %>%
  rename(Outcome = `Variable names`)

Table_S5 <- left_join(Table_S5, corres, by = "Outcome")
Table_S5 <- Table_S5 %>%
  mutate(
    `Variable labels` = fct_recode(`Variable labels`,
                                   "Shannon diversity" = "Shannon diversity in 1Y child gut microbiota ASV based (seq. depth = 5 000",
                                   "Specific richness" = "Specific richness diversity in 1Y child gut microbiota ASV based (seq. depth = 5 000",
                                   "Escherichia and Shigella" = "Escherichia_Shigella"),
    `Variable labels` = gsub("Escherichia_Shigella", "Escherichia and Shigella", `Variable labels`),
    `Variable labels` = gsub("_", " ", `Variable labels`),
    Phyla_corres = fct_recode(Phyla_corres, "Candidatus Saccharibacteria" = "Candidatus_Saccharibacteria"),
    Phyla_corres = fct_relevel(Phyla_corres, 
                               "Firmicutes", "Actinobacteria", "Bacteroidetes", "Proteobacteria", "Verrucomicrobia", "Candidatus Saccharibacteria"))%>%
  arrange(Phyla_corres, desc(Median))
write.xlsx(Table_S5, file = "4_output/review/Table_S5 (descrip outcomes).xlsx")

## Table S6 : alpha div - multivar (PIP) ----
# cf code 7.1_revision_1_EP_multipol.R

## Table S7 : beta div - univar ----
Table_S7 <- do.call(rbind, results_betadiv, quote = FALSE)
Table_S7 <-
  Table_S7 %>%
  mutate(
    `Explanatory variables` = if_else(`Explanatory variables` %in% phthalates_ter,
                                      `Explanatory variables`,
                                      str_remove(`Explanatory variables`, "\\d+$"))) %>%
  filter(`Explanatory variables` %in% phthalates_ter) %>%
  select(-Df, -R2) %>%
  mutate(
    Pollutants = factor(Pollutants, 
                        levels = phthalates_ter)) %>%
  arrange(Pollutants)
rm(phthalates_ter)

write_xlsx(Table_S7, "4_output/review/Table_S7.xlsx")

## Table S8 : phyla - univar ----
Table_S8 <- results_main %>% 
  mutate(
    Exposure = factor(Exposure, levels = c("ΣDEHP", "MnBP", "ΣDiNP", "MiBP", "MBzP", "MEP", "ohMPHP", "ΣDINCH")), 
    Window = factor(Window, levels = c("Trim.2", "Trim.3", "12 months")), 
    beta = format(round(beta, 1), nsmall = 1), 
    conf.low = format(round(conf.low, 1), nsmall = 1), 
    conf.high = format(round(conf.high, 1), nsmall = 1),
    p.value = case_when(p.value < 0.001 ~ "<0.001",
                        p.value < 0.01 ~ format(round(p.value, 3), nsmall = 3),
                        p.value > 0.01 ~ format(round(p.value, 2), nsmall = 2)), 
    CI = paste(conf.low, conf.high, sep = ", ")) %>%
  arrange(Outcome, Exposure, Window) %>%
  select(Outcome_name, Exposure_window_rec, beta, CI, p.value) %>% 
  filter(Outcome_name %in% c("Firmicutes", "Actinobacteria", "Bacteroidetes", "Proteobacteria")) %>% 
  pivot_wider(names_from = Outcome_name, values_from = c("beta", "CI", "p.value")) %>%
  select(Exposure_window_rec, 
         contains("Firmicutes"), 
         contains("Actinobacteria"), 
         contains("Bacteroidetes"), 
         contains("Proteobacteria")) %>%
  filter(!Exposure_window_rec %in% c("ΣDEHP 12 months", 
                                    "MEP 12 months", 
                                    "ohMPHP 12 months")) 

write_xlsx(Table_S8, "4_output/review/Table_S8.xlsx")

## Table S9 : phyla - multivar (PIP) ----
# cf code 7.1_revision_1_EP_multipol.R

## Table S10 : genera - univar ----
Table_S10 <- results_main %>% 
  mutate(
    Exposure = factor(Exposure, levels = c("ΣDEHP", "MnBP", "ΣDiNP", "MiBP", "MBzP", "MEP", "ohMPHP", "ΣDINCH")), 
    Window = factor(Window, levels = c("Trim.2", "Trim.3", "12 months")), 
    beta = format(round(beta, 1), nsmall = 1), 
    conf.low = format(round(conf.low, 1), nsmall = 1), 
    conf.high = format(round(conf.high, 1), nsmall = 1),
    p.value = case_when(p.value < 0.001 ~ "<0.001",
                        p.value < 0.01 ~ format(round(p.value, 3), nsmall = 3),
                        p.value > 0.01 ~ format(round(p.value, 2), nsmall = 2)), 
    CI = paste(conf.low, conf.high, sep = ", ")) %>%
  arrange(Outcome, Exposure, Window) %>%
  select(Outcome_name, Exposure_window_rec, beta, CI, p.value) %>% 
  filter(Outcome_name %in% genera_names) %>% 
  pivot_wider(names_from = Outcome_name, values_from = c("beta", "CI", "p.value")) %>% 
  select(
    Exposure_window_rec, 
    unlist(lapply(genera_names, function(name) {
      paste0(c("beta_", "CI_", "p.value_"), name)
    }))
  ) 
write_xlsx(Table_S10, "4_output/review/Table_S10.xlsx")

## Table S11 : sensitivity - rarefaction threshold ----
Table_S11_a <- results_main %>%                                                  # Colonnes analyses principales
  filter(Outcome %in% alpha_vec) %>%
  filter(p.value <0.05) %>% 
  select(Outcome_name, Exposure_window_rec, beta, conf.low, conf.high, p.value) %>%
  mutate(
    Outcome_name = factor(Outcome_name, levels = c("Specific richness", "Shannon diversity"))) %>%
  arrange(Outcome_name, desc(Exposure_window_rec))

bdd_Table_S11 <- bdd %>% filter(!is.na(ch_feces_SpecRich_10000_ASV_Y1))          
Table_S11_b <- map(alpha_vec, ~ main_model(.x, bdd_Table_S11))               # Colonnes analyses de sensi 5000 à l'effectif 339   
names(Table_S11_b) <- alpha_vec                                                  
Table_S11_b <- map_df(names(Table_S11_b), function(outcome) {            
  process_outcome_results(outcome, Table_S11_b[[outcome]])
})
Table_S11_b <- Table_S11_b %>% 
  filter(p.value <0.05 | 
         Outcome == "ch_feces_Shannon_5000_ASV_Y1" & Exposure == "ch_MBzP_i_cor_Y1_ln") %>%
  rename(beta_5000 = beta, 
         conf.low_5000 = conf.low, 
         conf.high_5000 = conf.high, 
         p.value_5000 = p.value,
         Exposure_window_rec = Exposure, 
         Outcome_name = Outcome) %>%
  mutate(
    Exposure_window_rec = str_replace_all(Exposure_window_rec,
                                      c("mo_" = "",
                                        "ch_" = "",
                                        "_i_cor_" = " ", 
                                        "_ln" = "",
                                        "ln" = "",
                                        "_cor" = "", 
                                        "t2" = "trim.2", 
                                        "t3" = "trim.3", 
                                        "Y1" = "12 months", 
                                        "_ms" = "", 
                                        "DEHP" = "ΣDEHP",
                                        "DiNP" = "ΣDiNP",
                                        "DINCH" = "ΣDINCH")), 
    Outcome_name = str_replace_all(Outcome_name, 
                                   c("ch_feces_SpecRich_5000_ASV_Y1" = "Specific richness", 
                                     "ch_feces_Shannon_5000_ASV_Y1" = "Shannon diversity"))) %>%
  select(Outcome_name, Exposure_window_rec, 
         beta_5000, conf.low_5000, conf.high_5000, p.value_5000)

Table_S11_c <- map(c("ch_feces_SpecRich_10000_ASV_Y1",                           # Colonnes analyses de sensi 10000 à l'effectif 339
                    "ch_feces_Shannon_10000_ASV_Y1" ), 
                  ~ main_model(.x, bdd_Table_S11))               
names(Table_S11_c) <- c("ch_feces_SpecRich_10000_ASV_Y1", "ch_feces_Shannon_10000_ASV_Y1" )                                                  
Table_S11_c <- map_df(names(Table_S11_c), function(outcome) {            
  process_outcome_results(outcome, Table_S11_c[[outcome]])
})
Table_S11_c <- Table_S11_c %>% 
  filter(p.value <0.05 | 
         Outcome == "ch_feces_Shannon_10000_ASV_Y1" & Exposure == "ch_MBzP_i_cor_Y1_ln")%>%
  rename(beta_10000 = beta, 
         conf.low_10000 = conf.low, 
         conf.high_10000 = conf.high, 
         p.value_10000 = p.value,
         Exposure_window_rec = Exposure, 
         Outcome_name = Outcome) %>%
  mutate(
    Exposure_window_rec = str_replace_all(Exposure_window_rec,
                                          c("mo_" = "",
                                            "ch_" = "",
                                            "_i_cor_" = " ", 
                                            "_ln" = "",
                                            "ln" = "",
                                            "_cor" = "", 
                                            "t2" = "trim.2", 
                                            "t3" = "trim.3", 
                                            "Y1" = "12 months", 
                                            "_ms" = "", 
                                            "DEHP" = "ΣDEHP",
                                            "DiNP" = "ΣDiNP",
                                            "DINCH" = "ΣDINCH")), 
    Outcome_name = str_replace_all(Outcome_name, 
                                   c("ch_feces_SpecRich_10000_ASV_Y1" = "Specific richness", 
                                     "ch_feces_Shannon_10000_ASV_Y1" = "Shannon diversity"))) %>%
  select(Outcome_name, Exposure_window_rec, 
         beta_10000, conf.low_10000, conf.high_10000, p.value_10000)

Table_S11 <- Table_S11_a %>%
  left_join(Table_S11_b, by = c("Outcome_name", "Exposure_window_rec")) %>%
  left_join(Table_S11_c, by = c("Outcome_name", "Exposure_window_rec")) 
rm(bdd_Table_S11, Table_S11_a, Table_S11_b, Table_S11_c)

Table_S11 <- Table_S11 %>%
  mutate(
    beta = format(round(beta, 2), nsmall = 2),
    beta_5000 = format(round(beta_5000, 2), nsmall = 2), 
    beta_10000 = format(round(beta_10000, 2), nsmall = 2),
    conf.low = format(round(conf.low, 1), nsmall = 1), 
    conf.high = format(round(conf.high, 1), nsmall = 1), 
    conf.low_5000 = format(round(conf.low_5000, 1), nsmall = 1),  
    conf.high_5000 = format(round(conf.high_5000, 1), nsmall = 1), 
    conf.low_10000 = format(round(conf.low_10000, 1), nsmall = 1), 
    conf.high_10000 = format(round(conf.high_10000, 1), nsmall = 1), 
    CI = paste(conf.low, conf.high, sep = ", "),
    CI_5000 = paste(conf.low_5000, conf.high_5000, sep = ", "),
    CI_10000 = paste(conf.low_10000, conf.high_10000, sep = ", "), 
    p.value = case_when(p.value < 0.001 ~ "<0.001",
                        p.value < 0.01 ~ format(round(p.value, 3), nsmall = 3),
                        p.value > 0.01 ~ format(round(p.value, 2), nsmall = 2)), 
    p.value_5000 = case_when(p.value_5000 < 0.001 ~ "<0.001",
                             p.value_5000 < 0.01 ~ format(round(p.value_5000, 3), nsmall = 3),
                             p.value_5000 > 0.01 ~ format(round(p.value_5000, 2), nsmall = 2)), 
    p.value_10000 = case_when(p.value_10000 < 0.001 ~ "<0.001",
                              p.value_10000 < 0.01 ~ format(round(p.value_10000, 3), nsmall = 3),
                              p.value_10000 > 0.01 ~ format(round(p.value_10000, 2), nsmall = 2))) %>%
  select(Outcome_name, 
         Exposure_window_rec, 
         beta, CI, p.value,
         beta_5000, CI_5000, p.value_5000, 
         beta_10000, CI_10000, p.value_10000)

write_xlsx(Table_S11, "4_output/review/Table_S11.xlsx")
  

## Table S12 : sensitivity - specific gravity ----
### utilisation des variables d'espo standardisés sur la gravité spécifique 
phthalates_sensi_sg_std <- phthalates %>% str_replace_all("_i_cor", "_i_cor_sg")
phthalates_sensi_sg_std_pre <- bdd %>% select(all_of(phthalates_sensi_sg_std)) %>% select(contains("t2"), contains("t3")) %>% colnames()
phthalates_sensi_sg_std_post <- bdd %>% select(all_of(phthalates_sensi_sg_std)) %>% select(contains("Y1")) %>% colnames()

Table_S12_a <- list()

for (i in 1:length(outcomes)) {                                                 # Boucle sur chaque outcome
  outcome <- outcomes[i]
  outcome_results <- list()
  
  for (j in 1:length(phthalates_sensi_sg_std_pre)) {                                # Boucle sur les variables explicatives "phthalates_sensi_sg_std_pre"
    exposure <- phthalates_sensi_sg_std_pre[j]
    regression_result <- sensi_12_a(outcome, exposure, covariates_pre)
    regression_result$Outcome_name <- outcome
    regression_result$Exposure_window_rec <- exposure
    outcome_results[[length(outcome_results) + 1]] <- regression_result
  }
  
  for (k in 1:length(phthalates_sensi_sg_std_post)) {                               # Boucle sur les variables explicatives "phthalates_sensi_sg_std_post"
    exposure <- phthalates_sensi_sg_std_post[k]
    regression_result <- sensi_12_a(outcome, exposure, covariates_post)
    regression_result$Outcome_name <- outcome
    regression_result$Exposure_window_rec <- exposure
    outcome_results[[length(outcome_results) + 1]] <- regression_result
  }
  
  Table_S12_a[[i]] <- outcome_results                                             # Stocker les résultats pour cet outcome
}

Table_S12_a <- do.call(rbind, lapply(Table_S12_a, function(outcome_list) {          # Convertir la liste de résultats en un seul tableau
  do.call(rbind, lapply(outcome_list, function(x) {
    as.data.frame(x)
  }))
}))

rm(i, j, k, outcome_results, regression_result, outcome, exposure)

Table_S12_a <- Table_S12_a %>%                                                      # Réorganisation tableau brut des resultats (pour figures)                   
  mutate(
    Exposure_window_rec = factor(Exposure_window_rec, levels = phthalates_sensi_sg_std), 
    Outcome_name = factor(Outcome_name, levels = outcomes), 
    Exposure_window_rec = str_replace_all(Exposure_window_rec,
                                          c("mo_" = "",
                                            "ch_" = "",
                                            "_i_cor_sg_" = " ", 
                                            "_ln" = "",
                                            "ln" = "",
                                            "_cor" = "", 
                                            "t2" = "trim.2", 
                                            "t3" = "trim.3", 
                                            "Y1" = "12 months", 
                                            "_ms" = "", 
                                            "DEHP" = "ΣDEHP",
                                            "DiNP" = "ΣDiNP",
                                            "DINCH" = "ΣDINCH")), 
    Outcome_name = str_replace_all(Outcome_name, 
                                   c("ch_feces_SpecRich_5000_ASV_Y1" = "Specific richness", 
                                     "ch_feces_Shannon_5000_ASV_Y1" = "Shannon diversity", 
                                     "ch_feces_rel_p1_Y1" = "Firmicutes", 
                                     "ch_feces_rel_p2_Y1" = "Actinobacteria", 
                                     "ch_feces_rel_p3_Y1" = "Bacteroidetes", 
                                     "ch_feces_rel_p4_Y1" = "Proteobacteria", 
                                     "Escherichia_Shigella" = "Escherichia and Shigella",
                                     "_" = " "))) %>%
  select(Outcome_name, Exposure_window_rec, 
         beta_sg_std, conf.low_sg_std, conf.high_sg_std, p.value_sg_std) %>%
  arrange(Outcome_name, Exposure_window_rec)

Table_S12_a <- left_join(
  results_main[, c("Outcome_name", "Exposure_window_rec", "beta",
                   "conf.low", "conf.high", "p.value")], 
  Table_S12_a, 
  by = c("Outcome_name", "Exposure_window_rec"))

rm(phthalates_sensi_sg_std, 
   phthalates_sensi_sg_std_pre, 
   phthalates_sensi_sg_std_post)

### utilisation de la standardisation sur les variables de la gravité spécifique 
covariates_t2_sensi_12_b <- c(covariates_pre, "mo_pool_sg_T1")
covariates_t3_sensi_12_b <- c(covariates_pre, "mo_pool_sg_T3")
covariates_post_sensi_12_b <- c(covariates_post, "ch_pool_sg_Y1")

phthalates_t2 <- bdd %>% select(all_of(phthalates)) %>% select(contains("t2")) %>% colnames()
phthalates_t3 <- bdd %>% select(all_of(phthalates)) %>% select(contains("t3")) %>% colnames()

Table_S12_b <- list()

for (i in 1:length(outcomes)) {                                                 # Boucle sur chaque outcome
  outcome <- outcomes[i]
  outcome_results <- list()
  
  for (j in 1:length(phthalates_t2)) {                                          # Boucle sur les variables explicatives "phthalates_t2"
    exposure <- phthalates_t2[j]
    regression_result <- sensi_12_b(outcome, exposure, covariates_t2_sensi_12_b)
    regression_result$Outcome_name <- outcome
    regression_result$Exposure_window_rec <- exposure
    outcome_results[[length(outcome_results) + 1]] <- regression_result
  }
  
  for (k in 1:length(phthalates_t3)) {                                          # Boucle sur les variables explicatives "phthalates_t3"
    exposure <- phthalates_t3[k]
    regression_result <- sensi_12_b(outcome, exposure, covariates_t3_sensi_12_b)
    regression_result$Outcome_name <- outcome
    regression_result$Exposure_window_rec <- exposure
    outcome_results[[length(outcome_results) + 1]] <- regression_result
  }
  
  for (m in 1:length(phthalates_post)) {                                        # Boucle sur les variables explicatives "phthalates_post"
    exposure <- phthalates_post[m]
    regression_result <- sensi_12_b(outcome, exposure, covariates_post_sensi_12_b)
    regression_result$Outcome_name <- outcome
    regression_result$Exposure_window_rec <- exposure
    outcome_results[[length(outcome_results) + 1]] <- regression_result
  }
  
  Table_S12_b[[i]] <- outcome_results                                           # Stocker les résultats pour cet outcome
}

Table_S12_b <- do.call(rbind, lapply(Table_S12_b, function(outcome_list) {      # Convertir la liste de résultats en un seul tableau
  do.call(rbind, lapply(outcome_list, function(x) {
    as.data.frame(x)
  }))
}))

rm(i, j, k, m, outcome_results, regression_result, outcome, exposure)

Table_S12_b <- Table_S12_b %>%                                                  # Réorganisation tableau brut des resultats (pour figures)                   
  mutate(
    Exposure_window_rec = factor(Exposure_window_rec, levels = phthalates), 
    Outcome_name = factor(Outcome_name, levels = outcomes), 
    Exposure_window_rec = str_replace_all(Exposure_window_rec,
                                          c("mo_" = "",
                                            "ch_" = "",
                                            "_i_cor_" = " ", 
                                            "_ln" = "",
                                            "ln" = "",
                                            "_cor" = "", 
                                            "t2" = "trim.2", 
                                            "t3" = "trim.3", 
                                            "Y1" = "12 months", 
                                            "_ms" = "", 
                                            "DEHP" = "ΣDEHP",
                                            "DiNP" = "ΣDiNP",
                                            "DINCH" = "ΣDINCH")), 
    Outcome_name = str_replace_all(Outcome_name, 
                                   c("ch_feces_SpecRich_5000_ASV_Y1" = "Specific richness", 
                                     "ch_feces_Shannon_5000_ASV_Y1" = "Shannon diversity", 
                                     "ch_feces_rel_p1_Y1" = "Firmicutes", 
                                     "ch_feces_rel_p2_Y1" = "Actinobacteria", 
                                     "ch_feces_rel_p3_Y1" = "Bacteroidetes", 
                                     "ch_feces_rel_p4_Y1" = "Proteobacteria", 
                                     "Escherichia_Shigella" = "Escherichia and Shigella",
                                     "_" = " "))) %>%
  select(Outcome_name, Exposure_window_rec, 
         beta_sg_adj, conf.low_sg_adj, conf.high_sg_adj, p.value_sg_adj) %>%
  arrange(Outcome_name, Exposure_window_rec)

Table_S12_b <- left_join(
  results_main[, c("Outcome_name", "Exposure_window_rec", "beta",
                   "conf.low", "conf.high", "p.value")], 
  Table_S12_b, 
  by = c("Outcome_name", "Exposure_window_rec"))

rm(covariates_t2_sensi_12_b, 
   covariates_t3_sensi_12_b, 
   covariates_post_sensi_12_b, 
   phthalates_t2, 
   phthalates_t3)

### Assemblage de Table_S12_a et Table_S12_b
Table_S12 <- right_join(Table_S12_a, Table_S12_b, 
                        by = c("Outcome_name", 
                               "Exposure_window_rec",
                               "beta", "conf.low", "conf.high", "p.value" ))

Table_S12 <- Table_S12 %>%
  filter(p.value <0.05 | p.value_sg_std <0.05 | p.value_sg_adj <0.05) %>%
  mutate(
    Outcome_name = factor(Outcome_name, levels = c("Specific richness", "Shannon diversity", 
                                                   "Firmicutes", "Actinobacteria", 
                                                   "Bacteroidetes", "Proteobacteria",
                                                   genera_names)), 
    
    conf.low = format(round(conf.low, 1), digits = 1),
    conf.high = format(round(conf.high, 1), digits = 1),
    conf.low_sg_std = format(round(conf.low_sg_std, 1), digits = 1),
    conf.high_sg_std = format(round(conf.high_sg_std, 1), digits = 1),
    conf.low_sg_adj = format(round(conf.low_sg_adj, 1), digits = 1),
    conf.high_sg_adj = format(round(conf.high_sg_adj, 1), digits = 1),
    
    "95% CI" = paste(conf.low, conf.high, sep = ", "), 
    "95% CI sg std" = paste(conf.low_sg_std, conf.high_sg_std, sep = ", "), 
    "95% CI sg adj" = paste(conf.low_sg_adj, conf.high_sg_adj, sep = ", "), 
    
    beta = ifelse(Outcome_name == "Shannon diversity", 
                  format(round(beta, 2), digits = 2),
                  format(round(beta, 1), digits = 1)), 
    beta_sg_std = ifelse(Outcome_name == "Shannon diversity", 
                         format(round(beta_sg_std, 2), digits = 2),
                         format(round(beta_sg_std, 1), digits = 1)), 
    beta_sg_adj = ifelse(Outcome_name == "Shannon diversity", 
                         format(round(beta_sg_adj, 2), digits = 2),
                         format(round(beta_sg_adj, 1), digits = 1)), 
    
    p.value = case_when(p.value < 0.001 ~ "<0.001",
                        p.value < 0.01 ~ format(round(p.value, 3), nsmall = 3),
                        p.value > 0.01 ~ format(round(p.value, 2), nsmall = 2)), 
    p.value_sg_std = case_when(p.value_sg_std < 0.001 ~ "<0.001",
                               p.value_sg_std < 0.01 ~ format(round(p.value_sg_std, 3), nsmall = 3),
                               p.value_sg_std > 0.01 ~ format(round(p.value_sg_std, 2), nsmall = 2)),
    p.value_sg_adj = case_when(p.value_sg_adj < 0.001 ~ "<0.001",
                               p.value_sg_adj < 0.01 ~ format(round(p.value_sg_adj, 3), nsmall = 3),
                               p.value_sg_adj > 0.01 ~ format(round(p.value_sg_adj, 2), nsmall = 2))) %>%
  arrange(Outcome_name) %>%
  select("Outcome_name", "Exposure_window_rec", 
         "beta", "95% CI", "p.value",
         "beta_sg_std", "95% CI sg std","p.value_sg_std",
         "beta_sg_adj", "95% CI sg adj","p.value_sg_adj") %>%
  arrange(Outcome_name)

write_xlsx(Table_S12, "4_output/review/Table_S12.xlsx")

rm(Table_S12_a, Table_S12_b)

## Table S13 : sensitivity - hospitalization ----
covariates_pre_sensi_13 <- c(covariates_pre, "ch_hospit_Y1")
covariates_post_sensi_13 <- c(covariates_post, "ch_hospit_Y1")

Table_S13 <- list()

for (i in 1:length(outcomes)) {                                                 # Boucle sur chaque outcome
  outcome <- outcomes[i]
  outcome_results <- list()
  
  for (j in 1:length(phthalates_pre)) {                                         # Boucle sur les variables explicatives "phthalates_pre"
    exposure <- phthalates_pre[j]
    regression_result <- sensi_13(outcome, exposure, covariates_pre_sensi_13)
    regression_result$Outcome_name <- outcome
    regression_result$Exposure_window_rec <- exposure
    outcome_results[[length(outcome_results) + 1]] <- regression_result
  }
  
  for (k in 1:length(phthalates_post)) {                                        # Boucle sur les variables explicatives "phthalates_post"
    exposure <- phthalates_post[k]
    regression_result <- sensi_13(outcome, exposure, covariates_post_sensi_13)
    regression_result$Outcome_name <- outcome
    regression_result$Exposure_window_rec <- exposure
    outcome_results[[length(outcome_results) + 1]] <- regression_result
  }
  
  Table_S13[[i]] <- outcome_results                                             # Stocker les résultats pour cet outcome
}

Table_S13 <- do.call(rbind, lapply(Table_S13, function(outcome_list) {          # Convertir la liste de résultats en un seul tableau
  do.call(rbind, lapply(outcome_list, function(x) {
    as.data.frame(x)
  }))
}))

rm(i, j, k, outcome_results, regression_result, outcome, exposure)

Table_S13 <- Table_S13 %>%                                        # Réorganisation tableau brut des resultats (pour figures)                   
  mutate(
    Exposure_window_rec = factor(Exposure_window_rec, levels = phthalates), 
    Outcome_name = factor(Outcome_name, levels = outcomes), 
    Exposure_window_rec = str_replace_all(Exposure_window_rec,
                                          c("mo_" = "",
                                            "ch_" = "",
                                            "_i_cor_" = " ", 
                                            "_ln" = "",
                                            "ln" = "",
                                            "_cor" = "", 
                                            "t2" = "trim.2", 
                                            "t3" = "trim.3", 
                                            "Y1" = "12 months", 
                                            "_ms" = "", 
                                            "DEHP" = "ΣDEHP",
                                            "DiNP" = "ΣDiNP",
                                            "DINCH" = "ΣDINCH")), 
    Outcome_name = str_replace_all(Outcome_name, 
                                   c("ch_feces_SpecRich_5000_ASV_Y1" = "Specific richness", 
                                     "ch_feces_Shannon_5000_ASV_Y1" = "Shannon diversity", 
                                     "ch_feces_rel_p1_Y1" = "Firmicutes", 
                                     "ch_feces_rel_p2_Y1" = "Actinobacteria", 
                                     "ch_feces_rel_p3_Y1" = "Bacteroidetes", 
                                     "ch_feces_rel_p4_Y1" = "Proteobacteria", 
                                     "Escherichia_Shigella" = "Escherichia and Shigella",
                                     "_" = " "))) %>%
  select(Outcome_name, Exposure_window_rec, 
         beta_hospit, conf.low_hospit, conf.high_hospit, p.value_hospit) %>%
  arrange(Outcome_name, Exposure_window_rec)

Table_S13 <- left_join(
  results_main[, c("Outcome_name", "Exposure_window_rec", "beta",
                   "conf.low", "conf.high", "p.value")], 
  Table_S13, 
  by = c("Outcome_name", "Exposure_window_rec"))

Table_S13 <- Table_S13 %>%
  filter(p.value <0.05 | p.value_hospit <0.05) %>%
  mutate(
    Outcome_name = factor(Outcome_name, levels = c("Specific richness", "Shannon diversity", 
                                                   "Firmicutes", "Actinobacteria", 
                                                   "Bacteroidetes", "Proteobacteria",
                                                   genera_names)), 
    conf.low = format(round(conf.low, 1), digits = 1),
    conf.high = format(round(conf.high, 1), digits = 1),
    conf.low_hospit = format(round(conf.low_hospit, 1), digits = 1),
    conf.high_hospit = format(round(conf.high_hospit, 1), digits = 1),
    "95% CI" = paste(conf.low, conf.high, sep = ", "), 
    "95% CI hospit" = paste(conf.low_hospit, conf.high_hospit, sep = ", "), 
    beta = ifelse(Outcome_name == "Shannon diversity", 
                  format(round(beta, 2), digits = 2),
                  format(round(beta, 1), digits = 1)), 
    beta_hospit = ifelse(Outcome_name == "Shannon diversity", 
                     format(round(beta_hospit, 2), digits = 2),
                     format(round(beta_hospit, 1), digits = 1)), 
    p.value = case_when(p.value < 0.001 ~ "<0.001",
                        p.value < 0.01 ~ format(round(p.value, 3), nsmall = 3),
                        p.value > 0.01 ~ format(round(p.value, 2), nsmall = 2)), 
    p.value_hospit = case_when(p.value_hospit < 0.001 ~ "<0.001",
                        p.value_hospit < 0.01 ~ format(round(p.value_hospit, 3), nsmall = 3),
                        p.value_hospit > 0.01 ~ format(round(p.value_hospit, 2), nsmall = 2))) %>%
  arrange(Outcome_name) %>%
  select("Outcome_name", "Exposure_window_rec", 
         "beta", "95% CI", "p.value",
         "beta_hospit", "95% CI hospit","p.value_hospit") %>%
  arrange(Outcome_name)

rm(covariates_pre_sensi_13, 
   covariates_post_sensi_13)

write_xlsx(Table_S13, "4_output/review/Table_S13.xlsx")


## Table S14 : sensitivity - gestational age, weight and length at birth and one year ----
covariates_pre_sensi_14 <- c(covariates_pre, "po_w_kg_3cat", "po_he_3cat_i", "ch_w_Y1_3cat_i", "ch_he_Y1_3cat_i", "po_gd")
covariates_post_sensi_14 <- c(covariates_post, "ch_w_Y1_3cat_i", "ch_he_Y1_3cat_i")

Table_S14 <- list()

for (i in 1:length(outcomes)) {                                                 # Boucle sur chaque outcome
  outcome <- outcomes[i]
  outcome_results <- list()
  
  for (j in 1:length(phthalates_pre)) {                                         # Boucle sur les variables explicatives "phthalates_pre"
    exposure <- phthalates_pre[j]
    regression_result <- sensi_14(outcome, exposure, covariates_pre_sensi_14)
    regression_result$Outcome_name <- outcome
    regression_result$Exposure_window_rec <- exposure
    outcome_results[[length(outcome_results) + 1]] <- regression_result
  }
  
  for (k in 1:length(phthalates_post)) {                                        # Boucle sur les variables explicatives "phthalates_post"
    exposure <- phthalates_post[k]
    regression_result <- sensi_14(outcome, exposure, covariates_post_sensi_14)
    regression_result$Outcome_name <- outcome
    regression_result$Exposure_window_rec <- exposure
    outcome_results[[length(outcome_results) + 1]] <- regression_result
  }
  
  Table_S14[[i]] <- outcome_results                                             # Stocker les résultats pour cet outcome
}

Table_S14 <- do.call(rbind, lapply(Table_S14, function(outcome_list) {          # Convertir la liste de résultats en un seul tableau
  do.call(rbind, lapply(outcome_list, function(x) {
    as.data.frame(x)
  }))
}))

rm(i, j, k, outcome_results, regression_result, outcome, exposure)


Table_S14 <- Table_S14 %>%                                                      # Réorganisation tableau brut des resultats (pour figures)                   
  mutate(
    Exposure_window_rec = factor(Exposure_window_rec, levels = phthalates), 
    Outcome_name = factor(Outcome_name, levels = outcomes), 
    Exposure_window_rec = str_replace_all(Exposure_window_rec,
                                          c("mo_" = "",
                                            "ch_" = "",
                                            "_i_cor_" = " ", 
                                            "_ln" = "",
                                            "ln" = "",
                                            "_cor" = "", 
                                            "t2" = "trim.2", 
                                            "t3" = "trim.3", 
                                            "Y1" = "12 months", 
                                            "_ms" = "", 
                                            "DEHP" = "ΣDEHP",
                                            "DiNP" = "ΣDiNP",
                                            "DINCH" = "ΣDINCH")), 
    Outcome_name = str_replace_all(Outcome_name, 
                                   c("ch_feces_SpecRich_5000_ASV_Y1" = "Specific richness", 
                                     "ch_feces_Shannon_5000_ASV_Y1" = "Shannon diversity", 
                                     "ch_feces_rel_p1_Y1" = "Firmicutes", 
                                     "ch_feces_rel_p2_Y1" = "Actinobacteria", 
                                     "ch_feces_rel_p3_Y1" = "Bacteroidetes", 
                                     "ch_feces_rel_p4_Y1" = "Proteobacteria", 
                                     "Escherichia_Shigella" = "Escherichia and Shigella",
                                     "_" = " "))) %>%
  select(Outcome_name, Exposure_window_rec, 
         beta_w_he, conf.low_w_he, conf.high_w_he, p.value_w_he) %>%
  arrange(Outcome_name, Exposure_window_rec)

Table_S14 <- left_join(
  results_main[, c("Outcome_name", "Exposure_window_rec", "beta",
                   "conf.low", "conf.high", "p.value")], 
  Table_S14, 
  by = c("Outcome_name", "Exposure_window_rec"))

Table_S14 <- Table_S14 %>%
  filter(p.value <0.05 | p.value_w_he <0.05) %>%
  mutate(
    Outcome_name = factor(Outcome_name, levels = c("Specific richness", "Shannon diversity", 
                                                   "Firmicutes", "Actinobacteria", 
                                                   "Bacteroidetes", "Proteobacteria",
                                                   genera_names)), 
    conf.low = format(round(conf.low, 1), digits = 1),
    conf.high = format(round(conf.high, 1), digits = 1),
    conf.low_w_he = format(round(conf.low_w_he, 1), digits = 1),
    conf.high_w_he = format(round(conf.high_w_he, 1), digits = 1),
    "95% CI" = paste(conf.low, conf.high, sep = ", "), 
    "95% CI w_he" = paste(conf.low_w_he, conf.high_w_he, sep = ", "), 
    beta = ifelse(Outcome_name == "Shannon diversity", 
                  format(round(beta, 2), digits = 2),
                  format(round(beta, 1), digits = 1)), 
    beta_w_he = ifelse(Outcome_name == "Shannon diversity", 
                       format(round(beta_w_he, 2), digits = 2),
                       format(round(beta_w_he, 1), digits = 1)), 
    p.value = case_when(p.value < 0.001 ~ "<0.001",
                        p.value < 0.01 ~ format(round(p.value, 3), nsmall = 3),
                        p.value > 0.01 ~ format(round(p.value, 2), nsmall = 2)), 
    p.value_w_he = case_when(p.value_w_he < 0.001 ~ "<0.001",
                        p.value_w_he < 0.01 ~ format(round(p.value_w_he, 3), nsmall = 3),
                        p.value_w_he > 0.01 ~ format(round(p.value_w_he, 2), nsmall = 2))) %>%
  arrange(Outcome_name) %>%
  select("Outcome_name", "Exposure_window_rec", 
         "beta", "95% CI", "p.value",
         "beta_w_he", "95% CI w_he","p.value_w_he") %>%
  arrange(Outcome_name)


rm(covariates_pre_sensi_14, 
   covariates_post_sensi_14)

write_xlsx(Table_S14, "4_output/review/Table_S14.xlsx")


## Figure S1 (dag) ----
# ok

## Figure S2 (heatmap cor expo) ----
bdd_figure_S2 <- bdd %>% 
  select(all_of(phthalates), 
         ch_feces_rel_p1_Y1) %>%
  filter(!is.na(ch_feces_rel_p1_Y1)) %>%
  select(-ch_feces_rel_p1_Y1)

colnames(bdd_figure_S2) <- colnames(bdd_figure_S2) %>%
  str_replace_all(c("mo_" = "",
                    "ch_" = "",
                    "_i_cor_" = " ", 
                    "_ln" = "",
                    "_t2" = " t2", 
                    "_t3" = " t3", 
                    "_Y1" = " Y1", 
                    "_ms" = "", 
                    "DINCH Y1" = "ΣDINCH Y1",
                    "DINCH t3" = "ΣDINCH t3",
                    "DINCH t2" = "ΣDINCH t2",
                    "DiNP Y1" = "ΣDiNP Y1",
                    "DiNP t3" = "ΣDiNP t3",
                    "DiNP t2" = "ΣDiNP t2",
                    "DEHP Y1" = "ΣDEHP Y1",
                    "DEHP t3" = "ΣDEHP t3",
                    "DEHP t2" = "ΣDEHP t2")) 

bdd_figure_S2 <- bdd_figure_S2 %>% na.omit()
for (nom in names(bdd_figure_S2)) {                   # enleve les étiquettes des variables 
  var_lab(bdd_figure_S2[[nom]]) <- NULL
}

bdd_figure_S2 <- bdd_figure_S2 %>%
  select("ΣDEHP t2","ΣDEHP t3", "ΣDEHP Y1",
         "MnBP t2", "MnBP t3", "MnBP Y1",
         "ΣDiNP t2", "ΣDiNP t3",  "ΣDiNP Y1",
         "MiBP t2", "MiBP t3", "MiBP Y1",
         "MBzP t2", "MBzP t3", "MBzP Y1",
         "MEP t2", "MEP t3", "MEP Y1",
         "ohMPHP t2", "ohMPHP t3", "ohMPHP Y1",
         "ΣDINCH t2", "ΣDINCH t3", "ΣDINCH Y1")

Figure_S2 <- cor(bdd_figure_S2, 
                          use = "pairwise.complete.obs", 
                          method = "pearson")

tiff(filename = "4_output/review/Figure_S2 (heatmap cor expo).tiff", 
     units = "mm", width = 250, height = 250, res = 300)
corrplot(Figure_S2, 
         method = 'color', 
         type = "lower", 
         tl.col = 'black', 
         tl.srt = 45, 
         addCoef.col = "black",
         number.cex = 0.8,
         tl.cex = 0.8,
         number.digits = 1,
         col = rev(COL2(diverging = "RdYlBu")))
dev.off()
rm(bdd_figure_S2)

## Figure S3 (heatmap cor expo and covariates) ----
Figure_S3 <- bdd %>%
  filter(!is.na(ch_feces_rel_p1_Y1)) %>%
  mutate(mo_par = as.numeric(as.character(mo_par))) %>%
  select(ident, 
         all_of(phthalates),
         "ch_feces_age_w_Y1_i",
         "ch_antibio_Y1_i",
         "mo_par",
         "po_w_kg", "po_he_i", 
         "ch_w_Y1_i", "po_he_i", 
         "po_gd", "mo_age",
         "mo_bmi_bepr_i", 
         "bf_duration_till48w_i") %>%
  na.omit()

Figure_S3 <- round(cor(Figure_S3, 
                    use = "pairwise.complete.obs", 
                    method = "spearman"), 1)

Figure_S3 <- Figure_S3 %>% 
  as.data.frame() %>% 
  select("ch_feces_age_w_Y1_i",
         "ch_antibio_Y1_i",
         "mo_par",
         "po_w_kg", "po_he_i", 
         "ch_w_Y1_i", "po_he_i", 
         "po_gd", "mo_age",
         "mo_bmi_bepr_i", 
         "bf_duration_till48w_i") %>% 
  t() %>%
  as.data.frame() %>%
  select(all_of(phthalates)) %>%
  as.matrix()

rownames(Figure_S3) <- rownames(Figure_S3) %>%
  str_replace_all(
    c("ch_feces_age_w_Y1_i" = "Child age (weeks)",
      "ch_antibio_Y1_i" = "Antibiotics use 0-12 months",
      "mo_par" = "Maternal parity",
      "po_w_kg" = "Birth weight (kg)",
      "po_he_i"= "Birth length (cm)", 
      "ch_w_Y1_i"="Weight at one year (kg)", 
      "ch_he_Y1_i"="Length at one year (cm)", 
      "po_gd"= "Gestational age (weeks)", 
      "mo_age"="Maternal age before pregnancy", 
      "mo_bmi_bepr_i"="Maternal BMI before pregnancy", 
      "bf_duration_till48w_i"="Breastfeeding duration (weeks)"))

colnames(Figure_S3) <- colnames(Figure_S3) %>%
  str_replace_all(
    c("mo_" = "",
      "ch_" = "",
      "_i_cor_" = " ", 
      "_ms" = "", 
      "DINCH" = "ΣDINCH",
      "DiNP" = "ΣDiNP",
      "DEHP" = "ΣDEHP", 
      "_ln" =""))
Figure_S3 <- Figure_S3 %>%
  as.data.frame() %>%
  select("ΣDEHP t2", "ΣDEHP t3", "ΣDEHP Y1",
         "MnBP t2", "MnBP t3", "MnBP Y1",
         "ΣDiNP t2", "ΣDiNP t3", "ΣDiNP Y1",
         "MiBP t2", "MiBP t3", "MiBP Y1",
         "MBzP t2", "MBzP t3", "MBzP Y1",
         "MEP t2", "MEP t3", "MEP Y1", 
         "ohMPHP t2", "ohMPHP t3", "ohMPHP Y1",
         "ΣDINCH t2", "ΣDINCH t3", "ΣDINCH Y1") %>%
  as.matrix()

tiff(filename = "4_output/review/Figure_S3 (heatmap cor expo covar).tiff", 
     units = "mm", width = 300, height = 150, res = 300)
corrplot(Figure_S3, 
         method = 'color', 
         tl.col = 'black', 
         tl.srt = 45, 
         addCoef.col = "black",
         number.cex = 0.8,
         tl.cex = 0.8,
         number.digits = 1,
         col = rev(COL2(diverging = "RdYlBu")))
dev.off()

## Figure S4 (heatmap cor outcomes) ----
Figure_S4 <- bdd %>%
  select(all_of(alpha_vec), 
         all_of(phyla_vec)) %>%
  rename("Specific richness" = "ch_feces_SpecRich_5000_ASV_Y1",
         "Shannon diversity" = "ch_feces_Shannon_5000_ASV_Y1", 
         "Firmicutes" = "ch_feces_rel_p1_Y1", 
         "Actinobacteria" = "ch_feces_rel_p2_Y1", 
         "Bacteroidetes" = "ch_feces_rel_p3_Y1", 
         "Proteobacteria" = "ch_feces_rel_p4_Y1")

Figure_S4 <- cor(Figure_S4, 
                 use = "pairwise.complete.obs", 
                 method = "pearson")

tiff(filename = "4_output/review/Figure_S4 (heatmap outcomes).tiff", 
     units = "mm", width = 150, height = 150, res = 300)
corrplot(Figure_S4, 
         method = 'color', 
         type = "lower", 
         tl.col = 'black', 
         tl.srt = 45, 
         addCoef.col = "black",
         number.cex = 0.8,
         tl.cex = 0.8,
         number.digits = 1,
         col = rev(COL2(diverging = "RdYlBu")))
dev.off()

## Figure S5 (alpha - multivar) ----
# cf code 7.1_revision_1_EP_multipol.R

## Figure S6 (phyla - multivar) ----
# cf code 7.1_revision_1_EP_multipol.R


# Exporter les résultats ----
rm(asv_raw_not_rarefied, 
   barplot, 
   boxplot, 
   comp_effectifs, 
   create_filtered_tbl_merge, 
   create_tbl_for_range, 
   custom_pvalue_fun, 
   densityplot, 
   descrip_num, 
   extract_exposure_results, 
   filter_exposure_tbl, 
   format_p_value, 
   generate_tbl, 
   heatmap_cor, 
   heatmap_cor_pairwise, 
   histogram, 
   main_model, 
   outliers, 
   process_outcome, 
   process_outcome_results, 
   scatterplot, 
   sensi_12_a, 
   sensi_12_b, 
   sensi_13, 
   sensi_14, 
   table_cor, 
   table_cor_sg,
   test_sensi_sg, 
   verif_distrib)                        
save.image("4_output/review/results_review_unipol.RData")

# Récupérer les résultats ----
load("4_output/review/results_review_unipol.RData")

# Autre points non enregistré dans les résultats ----
## réponse au reviewer EP lié au choix de la méthode de comparaison multiple ----
results_q_value_BH_div <- results_main %>%
  select(Outcome, Exposure_window, p.value) %>%
  filter(Outcome %in% alpha_vec) %>%
  mutate(
    q.value_div = p.adjust(p.value, method = "BH"))

results_q_value_BH_taxa <- results_main %>%
  select(Outcome, Exposure_window, p.value) %>%
  filter(!Outcome %in% alpha_vec) %>%
  mutate(
    q.value_taxa = p.adjust(p.value, method = "BH"))

results_main <- left_join(results_main, results_q_value_BH_div, by = c("Outcome", "Exposure_window", "p.value"))
results_main <- left_join(results_main, results_q_value_BH_taxa, by = c("Outcome", "Exposure_window", "p.value"))

justif_review_1 <- results_main %>% 
  select(Outcome, Exposure_window, beta, conf.low, conf.high, p.value, q.value_div) %>%
  filter(Outcome %in% alpha_vec) %>%
  filter(q.value_div <0.05) %>% 
  mutate(
    beta = format(round(beta, 2), nsmall = 2),
    conf.low = format(round(conf.low, 2), nsmall = 2),
    conf.high = format(round(conf.high, 2), nsmall = 2),
    `95%CI` = paste(conf.low, conf.high, sep = ", "),
    p.value = case_when(p.value < 0.001 ~ "<0.001",
                        p.value < 0.01 ~ format(round(p.value, 3), nsmall = 3),
                        p.value > 0.01 ~ format(round(p.value, 2), nsmall = 2)),
    q.value_div = case_when(q.value_div < 0.001 ~ "<0.001",
                            q.value_div < 0.01 ~ format(round(q.value_div, 3), nsmall = 3),
                            q.value_div > 0.01 ~ format(round(q.value_div, 2), nsmall = 2))) %>%
  select(Outcome, Exposure = Exposure_window, beta, `95%CI`, p.value, q.value = q.value_div)

justif_review_1_bis <- results_main %>% 
  select(Outcome, Exposure_window, beta, conf.low, conf.high, p.value, q.value_taxa) %>%
  filter(!Outcome %in% alpha_vec) %>%
  filter(q.value_taxa <0.05) %>% 
  mutate(
    beta = format(round(beta, 1), nsmall = 1),
    conf.low = format(round(conf.low, 1), nsmall = 1),
    conf.high = format(round(conf.high, 1), nsmall = 1),
    `95%CI` = paste(conf.low, conf.high, sep = ", "),
    p.value = case_when(p.value < 0.001 ~ "<0.001",
                        p.value < 0.01 ~ format(round(p.value, 3), nsmall = 3),
                        p.value > 0.01 ~ format(round(p.value, 2), nsmall = 2)),
    q.value_taxa = case_when(q.value_taxa < 0.001 ~ "<0.001",
                             q.value_taxa < 0.01 ~ format(round(q.value_taxa, 3), nsmall = 3),
                             q.value_taxa > 0.01 ~ format(round(q.value_taxa, 2), nsmall = 2))) %>%
  select(Outcome, Exposure = Exposure_window, beta, `95%CI`, p.value, q.value = q.value_taxa)


justif_review_1 <- rbind(justif_review_1, justif_review_1_bis)
rm(justif_review_1_bis)
write_xlsx(justif_review_1, "4_output/review/justif_review_1.xlsx")

## figure pour le résumé graphique ----
test <- results_main %>%
  mutate(
    FWER.p.value_taxa = ifelse(!Outcome_name %in% c("Specific richness", "Shannon diversity"), 
                                p.value * 14 * 33, NA),
    FWER.p.value_taxa = ifelse(FWER.p.value_taxa > 1, ">0.99", FWER.p.value_taxa),
    FWER.p.value_shape_taxa = case_when(p.value< 0.00011 & 
                                           !Outcome_name %in% c("Specific richness", 
                                                               "Shannon diversity")~ "p.value <0.00011",
                                         p.value > 0.00011 & 
                                           !Outcome_name %in% c("Specific richness", 
                                                               "Shannon diversity")~ "p.value >0.00011"), 
    FWER.p.value_shape_taxa = fct_relevel(FWER.p.value_shape_taxa,
                                           "p.value >0.00011", 
                                           "p.value <0.00011"), 
    
    FWER.p.value_shape = ifelse(is.na(FWER.p.value_alpha), as.character(FWER.p.value_shape_taxa), as.character(FWER.p.value_shape_alpha)), 
    FWER.p.value_shape = as.factor(FWER.p.value_shape)) %>%
  select(Outcome_name, Exposure_window, beta, conf.low, conf.high, p.value, FWER.p.value_shape) %>%
  filter(FWER.p.value_shape %in% c("p.value <0.00011", "p.value <0.0012"))


plot_abstract <- test %>%
  ggplot(aes(x = Exposure_window,
             y = beta,
             min = conf.low,
             ymax = conf.high,
             #color = interaction(exposure_window, term_2),
             #color = term_rec,
             color = FWER.p.value_shape)) +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_pointrange(position = position_dodge(width = 0.7), size = 0.4,
                  aes(color = FWER.p.value_shape)) +
  # labs(x = "Exposures", y = Outcome_name) +
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


## Associations entre la gravité spécifique et les paramètres de microbiote (non ajusté) ----
sg_vec <- c("mo_pool_sg_T1_rec", "mo_pool_sg_T3_rec", "ch_pool_sg_Y1_rec")

results_sg_outcome <- data.frame(
  Outcome = character(),
  Explicative = character(),
  Beta = numeric(),
  CI_Lower = numeric(),
  CI_Upper = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)

bdd <- bdd %>% 
  mutate(
    mo_pool_sg_T1_rec = mo_pool_sg_T1 * 10, 
    mo_pool_sg_T3_rec = mo_pool_sg_T3 * 10, 
    ch_pool_sg_Y1_rec = ch_pool_sg_Y1 * 10)

# Boucle pour faire les régressions linéaires
for (outcome in outcomes) {
  for (explicative in sg_vec) {
    formule <- as.formula(paste(outcome, "~", explicative))
    
    modele <- lm(formule, data = bdd)
    
    beta <- coef(modele)[2]  # Coefficient de la variable explicative
    p_value <- summary(modele)$coefficients[2, 4]  # p-value de la variable explicative
    ci <- confint(modele, level = 0.95)[2, ]  # CI pour le coefficient
    ci_lower <- ci[1]
    ci_upper <- ci[2]
    
    results_sg_outcome <- rbind(results_sg_outcome, data.frame(
      Outcome = outcome,
      Explicative = explicative,
      Beta = beta,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      P_value = p_value
    ))
  }
}

results_sg_outcome <- results_sg_outcome %>%
  mutate(
    Beta = ifelse(Outcome == "ch_feces_Shannon_5000_ASV_Y1", 
                  format(round(Beta, 2), digits = 2),
                  format(round(Beta, 1), digits = 1)), 
    CI_Lower = ifelse(Outcome == "ch_feces_Shannon_5000_ASV_Y1", 
                      format(round(CI_Lower, 2), digits = 2),
                      format(round(CI_Lower, 1), digits = 1)), 
    CI_Upper = ifelse(Outcome == "ch_feces_Shannon_5000_ASV_Y1", 
                      format(round(CI_Upper, 2), digits = 2),
                      format(round(CI_Upper, 1), digits = 1)), 
    CI = paste(CI_Lower, CI_Upper, sep = ","), 
    P_value = case_when(P_value < 0.0001 ~ format(round(P_value, 5), nsmall = 5),
                        P_value < 0.001 ~ format(round(P_value, 4), nsmall = 4),
                        P_value < 0.01 ~ format(round(P_value, 3), nsmall = 3),
                        P_value > 0.01 ~ format(round(P_value, 2), nsmall = 2)), 
    Outcome = fct_recode(
      Outcome, 
      "Specific richness" = "ch_feces_SpecRich_5000_ASV_Y1", 
      "Shannon diversity" = "ch_feces_Shannon_5000_ASV_Y1",
      "Firmicutes" = "ch_feces_rel_p1_Y1",
      "Actinobacteria" = "ch_feces_rel_p2_Y1", 
      "Bacteroidetes" = "ch_feces_rel_p3_Y1", 
      "Proteobacteria" = "ch_feces_rel_p4_Y1",
      "Clostridium IV" = "Clostridium_IV",
      "Clostridium sensu stricto" = "Clostridium_sensu_stricto",
      "Clostridium XlVa" = "Clostridium_XlVa",
      "Clostridium XVIII" = "Clostridium_XVIII",
      "Erysipelotrichaceae incertae sedis" = "Erysipelotrichaceae_incertae_sedis",
      "Escherichia and Shigella" = "Escherichia_Shigella",
      "Lachnospiracea incertae sedis" = "Lachnospiracea_incertae_sedis",
      "Saccharibacteria genera incertae sedis" = "Saccharibacteria_genera_incertae_sedis"),
    Outcome = factor(Outcome, levels = outcomes_names)) %>%
  select(Outcome, Explicative, Beta, CI, P_value)

rm(outcome, formule, modele, sg_vec, explicative, beta, p_value, ci, ci_lower, ci_upper)



## Associations entre la gravité spécifique et les paramètres de microbiote (ajusté) ----
sg_vec <- c("mo_pool_sg_T1_rec", "mo_pool_sg_T3_rec", "ch_pool_sg_Y1_rec")
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
               "po_gd",
               "mo_age",
               "mo_bmi_bepr_3cat_i",
               "bf_duration_till48w_4cat_i")

results_sg_outcome_adj <- data.frame(
  Outcome = character(),
  Explicative = character(),
  Beta = numeric(),
  CI_Lower = numeric(),
  CI_Upper = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)

bdd <- bdd %>% 
  mutate(
    mo_pool_sg_T1_rec = mo_pool_sg_T1 * 10, 
    mo_pool_sg_T3_rec = mo_pool_sg_T3 * 10, 
    ch_pool_sg_Y1_rec = ch_pool_sg_Y1 * 10)

# Boucle pour faire les régressions linéaires
for (outcome in outcomes) {
  for (explicative in sg_vec) {
    formule <- as.formula(paste(outcome, "~", explicative, "+", paste(covariates, collapse = " + ")))
    
    modele <- lm(formule, data = bdd)
    
    beta <- coef(modele)[2]  # Coefficient de la variable explicative
    p_value <- summary(modele)$coefficients[2, 4]  # p-value de la variable explicative
    ci <- confint(modele, level = 0.95)[2, ]  # CI pour le coefficient
    ci_lower <- ci[1]
    ci_upper <- ci[2]
    
    results_sg_outcome_adj <- rbind(results_sg_outcome_adj, data.frame(
      Outcome = outcome,
      Explicative = explicative,
      Beta = beta,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      P_value = p_value
    ))
  }
}

results_sg_outcome_adj <- results_sg_outcome_adj %>%
  mutate(
    Beta = ifelse(Outcome == "ch_feces_Shannon_5000_ASV_Y1", 
                  format(round(Beta, 2), digits = 2),
                  format(round(Beta, 1), digits = 1)), 
    CI_Lower = ifelse(Outcome == "ch_feces_Shannon_5000_ASV_Y1", 
                      format(round(CI_Lower, 2), digits = 2),
                      format(round(CI_Lower, 1), digits = 1)), 
    CI_Upper = ifelse(Outcome == "ch_feces_Shannon_5000_ASV_Y1", 
                      format(round(CI_Upper, 2), digits = 2),
                      format(round(CI_Upper, 1), digits = 1)), 
    CI = paste(CI_Lower, CI_Upper, sep = ","), 
    P_value = case_when(P_value < 0.0001 ~ format(round(P_value, 5), nsmall = 5),
                        P_value < 0.001 ~ format(round(P_value, 4), nsmall = 4),
                        P_value < 0.01 ~ format(round(P_value, 3), nsmall = 3),
                        P_value > 0.01 ~ format(round(P_value, 2), nsmall = 2)), 
    Outcome = fct_recode(
      Outcome, 
      "Specific richness" = "ch_feces_SpecRich_5000_ASV_Y1", 
      "Shannon diversity" = "ch_feces_Shannon_5000_ASV_Y1",
      "Firmicutes" = "ch_feces_rel_p1_Y1",
      "Actinobacteria" = "ch_feces_rel_p2_Y1", 
      "Bacteroidetes" = "ch_feces_rel_p3_Y1", 
      "Proteobacteria" = "ch_feces_rel_p4_Y1",
      "Clostridium IV" = "Clostridium_IV",
      "Clostridium sensu stricto" = "Clostridium_sensu_stricto",
      "Clostridium XlVa" = "Clostridium_XlVa",
      "Clostridium XVIII" = "Clostridium_XVIII",
      "Erysipelotrichaceae incertae sedis" = "Erysipelotrichaceae_incertae_sedis",
      "Escherichia and Shigella" = "Escherichia_Shigella",
      "Lachnospiracea incertae sedis" = "Lachnospiracea_incertae_sedis",
      "Saccharibacteria genera incertae sedis" = "Saccharibacteria_genera_incertae_sedis"),
    Outcome = factor(Outcome, levels = outcomes_names)) %>%
  select(Outcome, Explicative, Beta, CI, P_value)

rm(outcome, formule, modele, sg_vec, explicative, beta, p_value, ci, ci_lower, ci_upper, covariates)


results_sg_outcome_adj <- results_sg_outcome_adj %>% 
  rename(Beta_adj = Beta, 
         CI_adj = CI, 
         P_value_adj = P_value)
results_sg_outcome <- left_join(results_sg_outcome, 
                                results_sg_outcome_adj, 
                                by = c("Outcome", "Explicative"))
rm(results_sg_outcome_adj)

results_sg_outcome %>% filter(P_value <0.05) %>% View()
results_sg_outcome %>% filter(Explicative == "mo_pool_sg_T1_rec") %>% filter(P_value <0.05 |  P_value_adj <0.05) %>% View()
results_sg_outcome %>% filter(Explicative == "mo_pool_sg_T3_rec") %>% filter(P_value <0.05 |  P_value_adj <0.05) %>% View()
results_sg_outcome %>% filter(Explicative == "ch_pool_sg_Y1_rec") %>% filter(P_value <0.05 |  P_value_adj <0.05) %>% View()


# bdd <- bdd %>%
#   mutate(rowSums_genera = rowSums(bdd[, genera_total]))
# 
# cor.test(bdd$ch_feces_SpecRich_5000_ASV_Y1, bdd$rowSums_genera)
# cor.test(bdd$ch_feces_Shannon_5000_ASV_Y1, bdd$rowSums_genera)


## Associations GS --> pththalates ----
phthalates_t2 <- bdd %>% select(all_of(phthalates)) %>% select(contains("t2")) %>% colnames()
phthalates_t3 <- bdd %>% select(all_of(phthalates)) %>% select(contains("t3")) %>% colnames()

bdd <- bdd %>% 
  mutate(
    mo_pool_sg_T1_rec = mo_pool_sg_T1 * 10, 
    mo_pool_sg_T3_rec = mo_pool_sg_T3 * 10, 
    ch_pool_sg_Y1_rec = ch_pool_sg_Y1 * 10)

list_exp_outcomes <- list(
  mo_pool_sg_T1_rec = phthalates_t2,
  mo_pool_sg_T3_rec = phthalates_t3,
  ch_pool_sg_Y1_rec = phthalates_post)


run_lm <- function(explicative, outcome) {
  model <- lm(as.formula(paste(outcome, "~", explicative)), data = bdd)
  
  tidy_results <- tidy(model) %>% filter(term == explicative) 
  confint_results <- confint(model, level = 0.95)[explicative, ]

  tidy_results %>%
    mutate(conf.low = confint_results[1],  
           conf.high = confint_results[2],
           Outcome = outcome,             
           Explicative = explicative)  
}

results_sg_exposure <- bind_rows(
  lapply(names(list_exp_outcomes), function(explicative) {
    outcomes <- list_exp_outcomes[[explicative]]
    bind_rows(lapply(outcomes, function(outcome) {
      run_lm(explicative, outcome)
    }))
  })
)

results_sg_exposure <- results_sg_exposure %>%
  select(Outcome, Explicative, estimate, conf.low, conf.high, p.value) %>%
  rename(Beta = estimate,
         P_value = p.value) %>%
  mutate(
      Beta = format(round(Beta, 2), digits = 2), 
      conf.low = format(round(conf.low, 2), digits = 2),
      conf.high = format(round(conf.high, 2), digits = 2),
      CI = paste(conf.low, conf.high, sep = ","), 
      P_value = as.numeric(P_value),
      P_value = format(round(P_value, 2), digits = 2)) %>%
  select(Outcome, Explicative, Beta, CI, P_value)

rm(phthalates_t2, phthalates_t3, list_exp_outcomes, run_lm)


write_xlsx(
  list("sensi_sg" = Table_S12, 
       "results_sg_outcome" = results_sg_outcome, 
       "results_sg_exposure" = results_sg_exposure), 
  "review_1_specific_gravity.xlsx")
