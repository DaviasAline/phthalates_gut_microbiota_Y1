# Analyses BKMR exposition aux polluants et microbiote 1 an 
# A. Davias
# 23/04/2023

# Chargement des données ----
load("2_final_data/metadata.RData")
load("2_final_data/bdd_alpha.RData")
load("2_final_data/bdd_taxa.RData")

# Chargement des packages ----
library(bkmr)
library(fields)
library(future)
library(future.apply)
library(writexl)
library(ggtext)
source("3_programs/4_functions_AD_gumme.R", encoding = 'UTF-8')    

# Chargement des vecteurs ----
source("3_programs/4_vectors_AD_gumme.R", echo=FALSE)
rm(phenols_vec, phenols_vec_2, phenols_vec_cat, phenols_vec_ln, phenols_vec_num, phenols_vec_num_2, phenols_vec_num_M2, 
   phenols_vec_num_sg, phenols_vec_num_sg_ln, phenols_vec_num_sg_M2, phenols_vec_num_sg_t2, phenols_vec_num_sg_t3, 
   phenols_vec_num_sg_Y1, phenols_vec_num_t2, phenols_vec_num_t3, phenols_vec_num_Y1, phenols_vec_ter, 
   pfas_vec, pfas_vec_2, pfas_vec_3, pfas_vec_cat, pfas_vec_exclu17673, pfas_vec_ln, pfas_vec_num, pfas_vec_ter, 
   alpha_vec, 
   phthalates_vec_num_M2, phthalates_vec_num_sg, phthalates_vec_num_sg_ln, phthalates_vec_num_sg_M2, phthalates_vec_num_sg_t2, 
   phthalates_vec_num_sg_t3, phthalates_vec_num_sg_Y1, 
   pollutant_vec_M2, pollutant_vec_t2, pollutant_vec_t3, pollutant_vec_Y1, 
   comp_effectifs, heatmap_cor_pairwise, model_covar, model_multi, model_summary, model_univ_multi, 
   table_cor, table_cor_sg, test_sensi_sg, 
   bdd_alpha)

phthalates_vec_t2 <- bdd_taxa %>% 
  select(all_of(phthalates_vec)) %>%
  select(contains("t2")) %>%
  select(!contains(c("MEOHP", "MECPP", "MEHHP", "MEHP", "MMCHP",     
                     "ohMiNP", "oxoMiNP", "cxMiNP",                 
                     "ohMINCH", "oxoMINCH"))) %>%  
  select(!contains("cat")) %>%
  colnames()

phthalates_vec_t3 <- bdd_taxa %>% 
  select(all_of(phthalates_vec)) %>%
  select(contains("t3"))%>%
  select(!contains(c("MEOHP", "MECPP", "MEHHP", "MEHP", "MMCHP",     
                     "ohMiNP", "oxoMiNP", "cxMiNP",                 
                     "ohMINCH", "oxoMINCH"))) %>%  
  select(!contains("cat")) %>%
  colnames()

phthalates_vec_Y1 <- bdd_taxa %>% 
  select(all_of(phthalates_vec)) %>%
  select(contains("Y1"))%>%
  select(!contains(c("MEOHP", "MECPP", "MEHHP", "MEHP", "MMCHP",     
                     "ohMiNP", "oxoMiNP", "cxMiNP",                 
                     "ohMINCH", "oxoMINCH"))) %>%  
  select(!contains("cat")) %>%
  colnames()

# Chargement des fonctions ----
bkmr_t2 <- function(outcome) {
  set.seed(111)
  results_t2 <- kmbayes(    # modèle sans les variables catégorielles "numérisées" / toutes fenetres confondues
    y = outcome, 
    Z = mixture_taxa_t2, 
    X = covariates_taxa_t2, 
    iter = 50000,          # mettre 50 000
    verbose = FALSE,       # if TRUE, la sortie intermédiaire résumant la progression de l'ajustement du modèle est imprimée
    varsel = TRUE)         # if TRUE, we can fit the model with variable selection and estimate the posterior inclusion probability (PIP) for each of the exposures zim
  return(results_t2)
}

bkmr_t3 <- function(outcome) {
  set.seed(111)
  results_t3 <- kmbayes(    # modèle avec les variables catégorielles "numérisées" / toutes fenetres confondues
    y = outcome, 
    Z = mixture_taxa_t3, 
    X = covariates_taxa_t3, 
    iter = 50000,          # mettre 50 000
    verbose = FALSE,       # if TRUE, la sortie intermédiaire résumant la progression de l'ajustement du modèle est imprimée
    varsel = TRUE)         # if TRUE, we can fit the model with variable selection and estimate the posterior inclusion probability (PIP) for each of the exposures zim
  return(results_t3)
}

bkmr_Y1 <- function(outcome) {
  set.seed(111)
  results_Y1 <- kmbayes(    # modèle avec les variables catégorielles "numérisées" / fenetres séparées
    y = outcome, 
    Z = mixture_taxa_Y1, 
    X = covariates_taxa_Y1, 
    iter = 50000,          # mettre 50 000
    verbose = FALSE,       # if TRUE, la sortie intermédiaire résumant la progression de l'ajustement du modèle est imprimée
    varsel = TRUE)         # if TRUE, we can fit the model with variable selection and estimate the posterior inclusion probability (PIP) for each of the exposures zim
  return(results_Y1)
}


TracePlot_group <- function(model_p1, model_p2, model_p3, model_p4, titre){
  par(mfrow=c(4,3))
  TracePlot(fit = model_p1, par = "beta", sel = TRUE) 
  TracePlot(fit = model_p1, par = "sigsq.eps", sel = TRUE) 
  TracePlot(fit = model_p1, par = "r", comp = 1, sel = TRUE) 
  TracePlot(fit = model_p2, par = "beta", sel = TRUE) 
  TracePlot(fit = model_p2, par = "sigsq.eps", sel = TRUE) 
  TracePlot(fit = model_p2, par = "r", comp = 1, sel = TRUE) 
  TracePlot(fit = model_p3, par = "beta", sel = TRUE) 
  TracePlot(fit = model_p3, par = "sigsq.eps", sel = TRUE) 
  TracePlot(fit = model_p3, par = "r", comp = 1, sel = TRUE) 
  TracePlot(fit = model_p4, par = "beta", sel = TRUE) 
  TracePlot(fit = model_p4, par = "sigsq.eps", sel = TRUE) 
  TracePlot(fit = model_p4, par = "r", comp = 1, sel = TRUE) 
  mtext(titre, outer=TRUE, line = -1.5,font=2, cex=1, padj = 0)
}

pip_results <- function(bkmr_p1, bkmr_p2, bkmr_p3, bkmr_p4,
                        bkmr_g1, bkmr_g2, bkmr_g3, bkmr_g4) {
  
  bkmr_pip_p1 <- ExtractPIPs(bkmr_p1) %>% as.data.frame() %>% rename(PIP_p1 = PIP) 
  bkmr_pip_p2 <- ExtractPIPs(bkmr_p2) %>% as.data.frame() %>% rename(PIP_p2 = PIP) 
  bkmr_pip_p3 <- ExtractPIPs(bkmr_p3) %>% as.data.frame() %>% rename(PIP_p3 = PIP) 
  bkmr_pip_p4 <- ExtractPIPs(bkmr_p4) %>% as.data.frame() %>% rename(PIP_p4 = PIP) 

  bkmr_pip_g1 <- ExtractPIPs(bkmr_g1) %>% as.data.frame() %>% rename(PIP_g1 = PIP) 
  bkmr_pip_g2 <- ExtractPIPs(bkmr_g2) %>% as.data.frame() %>% rename(PIP_g2 = PIP) 
  bkmr_pip_g3 <- ExtractPIPs(bkmr_g3) %>% as.data.frame() %>% rename(PIP_g3 = PIP) 
  bkmr_pip_g4 <- ExtractPIPs(bkmr_g4) %>% as.data.frame() %>% rename(PIP_g4 = PIP) 
  
  results <- list(bkmr_pip_p1, bkmr_pip_p2, bkmr_pip_p3, bkmr_pip_p4,
                  bkmr_pip_g1, bkmr_pip_g2, bkmr_pip_g3, bkmr_pip_g4)
  results <- reduce(results, left_join, by = "variable")
  return(results)
  
}

risks_overall <- function(fit, y, Z, covariates) {
  results <- OverallRiskSummaries(
    fit = fit, 
    y = y,
    Z = Z,
    X = covariates,
    qs = seq(0.25, 0.95, by = 0.10),
    q.fixed = 0.25,
    method = "exact")  # approx utilisé en run 5, exact utilisé en run 6
  return(results)
}     

risks_singvar <- function(fit, y, Z, covariates) {
  results <- SingVarRiskSummaries(
    fit = fit, 
    y = y,
    Z = Z,   
    X = covariates,
    qs.diff = c(0.25, 0.75),
    q.fixed = c(0.25, 0.50, 0.75),
    method = "exact")
  return(results)
}

plot_risks.overall <- function(risks.overall, title, y_title, x_title){
  ggplot(risks.overall,
         aes(
           quantile,
           est,
           ymin = est - 1.96 * sd,
           ymax = est + 1.96 * sd
         )) +
    geom_pointrange() +
    ylab(y_title) +           
    xlab(x_title) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
    theme_bw() +
    theme(plot.title = element_markdown(size = 12, hjust = 0.5))+ 
    ggtitle(title) 
    
}

plot_risks.singvar <- function(risks.singvar, window, taxa, 
                               title, y_title, x_title, 
                               legend.position, axis.y) {
  
  risks.singvar %>%
    mutate(
      pollutant = sub("...$", "", pollutant)) %>%
    filter({{window}} == window) %>%
    filter({{taxa}} == taxa) %>%
    ggplot(
      aes(
        pollutant,
        est,
        ymin = est - 1.96 * sd,
        ymax = est + 1.96 * sd,
        col = q.fixed)
    ) +
    geom_pointrange(position = position_dodge(width = 0.75), 
                    size = 0.15) +
    geom_hline(yintercept = 0, linetype="dashed") +
    labs(y = y_title, 
         x =x_title) + 
    coord_flip()+
    theme_bw() +
    theme(legend.position = legend.position, 
          plot.title = element_markdown(size = 12, hjust = 0.5), 
          axis.text.y = axis.y)+ 
    ggtitle(title)}

# Data cleaning ----
## Outcome variables ----
bdd_outcomes_taxa <- 
  bdd_taxa %>% 
  select(
    ident, 
    all_of(taxa_vec)) %>%
  na.omit()

## Covariates ----
bdd_covariates_i <-          # on force les variables catégorielles en variables continues 
  metadata %>%
  select(
    ident, 
    all_of(covar_vec), 
    all_of(covar_vec_i), 
    all_of(covar_vec_num_i)) %>%
  select(
    ident,  
    ch_feces_RUN_Y1,              # covariable à 2 catégories : ok 0/1
    ch_feces_age_w_Y1_i,  
    po_delmod,                    # covariable à 2 catégories : ok 0/1
    ch_food_intro_Y1_3cat_i,      # covariable à transformer 
    ch_antibio_Y1_i,         
    mo_par,  
    mo_pets_i,                    # covariable à 2 catégories : ok 0/1
    ch_sex,                       # covariable à 2 catégories : ok 0/1
    mo_tob_gr_anyt_yn_n2_i,       # covariable à 2 catégories : ok 0/1
    Mo_ETS_anyT_yn1_opt_i,        # covariable à 2 catégories : ok 0/1
    ch_ETS_12m_opt36m,            # covariable à 2 catégories : ok 0/1
    mo_interpreg_3cat,            # covariable à transformer 
    mo_dipl_3cat_i,               # covariable à transformer
    po_w_kg,
    po_he_i,
    ch_w_Y1_i,
    ch_he_Y1_i,
    po_gd,
    mo_age, 
    mo_bmi_bepr_i,
    bf_duration_till48w_i) %>%
  mutate(
    ch_feces_RUN_Y1 = as.numeric(fct_recode(ch_feces_RUN_Y1,
                                            "0" = "R2",
                                            "1" = "R3")), 
    po_delmod = as.numeric(fct_recode(po_delmod,
                                      "0" = "C-section",
                                      "1" = "Vaginal delivery")),
    mo_pets_i = as.numeric(fct_recode(mo_pets_i,
                                      "0" = "No",
                                      "1" = "One or more")),
    ch_sex = as.numeric(fct_recode(ch_sex,
                                   "0" = "Female",
                                   "1" = "Male")),
    mo_tob_gr_anyt_yn_n2_i = as.numeric(fct_recode(mo_tob_gr_anyt_yn_n2_i,
                                                   "0" = "No",
                                                   "1" = "Yes")),
    Mo_ETS_anyT_yn1_opt_i = as.numeric(fct_recode(Mo_ETS_anyT_yn1_opt_i,
                                                  "0" = "No",
                                                  "1" = "Yes")), 
    ch_ETS_12m_opt36m = as.numeric(fct_recode(ch_ETS_12m_opt36m,
                                              "0" = "No",
                                              "1" = "Yes")), 
    
    mo_interpreg_3cat_2years_and_more = as.numeric(fct_recode(mo_interpreg_3cat,
                                                              "0" = "Under 2 years",
                                                              "1" = "2 years and more",
                                                              "0" = "Primiparous")),
    mo_interpreg_3cat_primiparous = as.numeric(fct_recode(mo_interpreg_3cat,
                                                          "0" = "Under 2 years",
                                                          "0" = "2 years and more",
                                                          "1" = "Primiparous")), 
    
    ch_food_intro_Y1_3cat_i_6_12m = as.numeric(fct_recode(ch_food_intro_Y1_3cat_i,
                                                          "0" = "Between 0 and 6 months old",
                                                          "1" = "Between 6 and 12 months old",
                                                          "0" = "Not introduced at 12 months old")),
    ch_food_intro_Y1_3cat_i_not_intro = as.numeric(fct_recode(ch_food_intro_Y1_3cat_i,
                                                              "0" = "Between 0 and 6 months old",
                                                              "0" = "Between 6 and 12 months old",
                                                              "1" = "Not introduced at 12 months old")),
    mo_dipl_3cat_i_3_4y = as.numeric(fct_recode(mo_dipl_3cat_i,
                                                "0" = "2years or less after graduation",
                                                "1" = "3-4years after graduation",
                                                "0" = ">=5years after graduation")), 
    mo_dipl_3cat_i_5y = as.numeric(fct_recode(mo_dipl_3cat_i,
                                              "0" = "2years or less after graduation",
                                              "0" = "3-4years after graduation",
                                              "1" = ">=5years after graduation"))) %>%
  select(-c("mo_interpreg_3cat", "mo_dipl_3cat_i", "ch_food_intro_Y1_3cat_i"))
keep_covar <- bdd_covariates_i[, 2:25] %>% colnames()


## Exposure variables ----
bdd_expo_t2 <- metadata %>%          
  select(
    ident,
    all_of(phthalates_vec_t2))%>%           # variables d'exposition (déjà logtransformé)
  na.omit()                             # conserver que les lignes avec toutes les observations completes sans données manquantes 
keep_expo_t2 <- colnames(bdd_expo_t2[, 2:9])
bdd_expo_t2[, keep_expo_t2] <- lapply(bdd_expo_t2[, keep_expo_t2], scale)  # standardiser les expositions (diviser par sd)
colnames(bdd_expo_t2) <- c("ident", keep_expo_t2)

bdd_expo_t3 <- metadata %>%          
  select(
    ident,
    all_of(phthalates_vec_t3))%>%           # variables d'exposition (déjà logtransformé)
  na.omit()                             # conserver que les lignes avec toutes les observations completes sans données manquantes 
keep_expo_t3 <- colnames(bdd_expo_t3[, 2:9])
bdd_expo_t3[, keep_expo_t3] <- lapply(bdd_expo_t3[, keep_expo_t3], scale)  # standardiser les expositions (diviser par sd)
colnames(bdd_expo_t3) <- c("ident", keep_expo_t3)

bdd_expo_Y1 <- metadata %>%          
  select(
    ident,
    all_of(phthalates_vec_Y1))%>%           # variables d'exposition (déjà logtransformé)
  na.omit()                             # conserver que les lignes avec toutes les observations completes sans données manquantes 
keep_expo_Y1 <- colnames(bdd_expo_Y1[, 2:9])
bdd_expo_Y1[, keep_expo_Y1] <- lapply(bdd_expo_Y1[, keep_expo_Y1], scale)  # standardiser les expositions (diviser par sd)
colnames(bdd_expo_Y1) <- c("ident", keep_expo_Y1)

# Préparation des matrices ----
bdd_bkmr_taxa_t2 <- 
  list(bdd_outcomes_taxa, bdd_covariates_i, bdd_expo_t2) %>%
  reduce(left_join, by = "ident") %>%
  na.omit()

mixture_taxa_t2 <- bdd_bkmr_taxa_t2 %>% select(all_of(keep_expo_t2)) %>% as.matrix()
covariates_taxa_t2 <- bdd_bkmr_taxa_t2 %>% select(all_of(keep_covar)) %>% as.matrix()
outcome_p1_t2 <- bdd_bkmr_taxa_t2 %>% select(ch_feces_rel_p1_Y1) %>% as.matrix()
outcome_p2_t2 <- bdd_bkmr_taxa_t2 %>% select(ch_feces_rel_p2_Y1) %>% as.matrix()
outcome_p3_t2 <- bdd_bkmr_taxa_t2 %>% select(ch_feces_rel_p3_Y1) %>% as.matrix()
outcome_p4_t2 <- bdd_bkmr_taxa_t2 %>% select(ch_feces_rel_p4_Y1) %>% as.matrix()
outcome_g1_t2 <- bdd_bkmr_taxa_t2 %>% select(ch_feces_rel_g1_Y1) %>% as.matrix()
outcome_g2_t2 <- bdd_bkmr_taxa_t2 %>% select(ch_feces_rel_g2_Y1) %>% as.matrix()
outcome_g3_t2 <- bdd_bkmr_taxa_t2 %>% select(ch_feces_rel_g3_Y1) %>% as.matrix()
outcome_g4_t2 <- bdd_bkmr_taxa_t2 %>% select(ch_feces_rel_g4_Y1) %>% as.matrix()

bdd_bkmr_taxa_t3 <- 
  list(bdd_outcomes_taxa, bdd_covariates_i, bdd_expo_t3) %>%
  reduce(left_join, by = "ident") %>%
  na.omit()

mixture_taxa_t3 <- bdd_bkmr_taxa_t3 %>% select(all_of(keep_expo_t3)) %>% as.matrix()
covariates_taxa_t3 <- bdd_bkmr_taxa_t3 %>% select(all_of(keep_covar)) %>% as.matrix()
outcome_p1_t3 <- bdd_bkmr_taxa_t3 %>% select(ch_feces_rel_p1_Y1) %>% as.matrix()
outcome_p2_t3 <- bdd_bkmr_taxa_t3 %>% select(ch_feces_rel_p2_Y1) %>% as.matrix()
outcome_p3_t3 <- bdd_bkmr_taxa_t3 %>% select(ch_feces_rel_p3_Y1) %>% as.matrix()
outcome_p4_t3 <- bdd_bkmr_taxa_t3 %>% select(ch_feces_rel_p4_Y1) %>% as.matrix()
outcome_g1_t3 <- bdd_bkmr_taxa_t3 %>% select(ch_feces_rel_g1_Y1) %>% as.matrix()
outcome_g2_t3 <- bdd_bkmr_taxa_t3 %>% select(ch_feces_rel_g2_Y1) %>% as.matrix()
outcome_g3_t3 <- bdd_bkmr_taxa_t3 %>% select(ch_feces_rel_g3_Y1) %>% as.matrix()
outcome_g4_t3 <- bdd_bkmr_taxa_t3 %>% select(ch_feces_rel_g4_Y1) %>% as.matrix()

bdd_bkmr_taxa_Y1 <- 
  list(bdd_outcomes_taxa, bdd_covariates_i, bdd_expo_Y1) %>%
  reduce(left_join, by = "ident") %>%
  na.omit()

mixture_taxa_Y1 <- bdd_bkmr_taxa_Y1 %>% select(all_of(keep_expo_Y1)) %>% as.matrix()
covariates_taxa_Y1 <- bdd_bkmr_taxa_Y1 %>% select(all_of(keep_covar)) %>% as.matrix()
outcome_p1_Y1 <- bdd_bkmr_taxa_Y1 %>% select(ch_feces_rel_p1_Y1) %>% as.matrix()
outcome_p2_Y1 <- bdd_bkmr_taxa_Y1 %>% select(ch_feces_rel_p2_Y1) %>% as.matrix()
outcome_p3_Y1 <- bdd_bkmr_taxa_Y1 %>% select(ch_feces_rel_p3_Y1) %>% as.matrix()
outcome_p4_Y1 <- bdd_bkmr_taxa_Y1 %>% select(ch_feces_rel_p4_Y1) %>% as.matrix()
outcome_g1_Y1 <- bdd_bkmr_taxa_Y1 %>% select(ch_feces_rel_g1_Y1) %>% as.matrix()
outcome_g2_Y1 <- bdd_bkmr_taxa_Y1 %>% select(ch_feces_rel_g2_Y1) %>% as.matrix()
outcome_g3_Y1 <- bdd_bkmr_taxa_Y1 %>% select(ch_feces_rel_g3_Y1) %>% as.matrix()
outcome_g4_Y1 <- bdd_bkmr_taxa_Y1 %>% select(ch_feces_rel_g4_Y1) %>% as.matrix()


# Running BKMR ----
## Run 5 ----  
# Spécifier le planificateur future et le nombre de noyaux à utiliser
plan(multisession)
ncores <- availableCores() # ou spécifier le nombre de noyaux à utiliser

# Utiliser future_lapply pour exécuter les appels bkmr en parallèle pour les 3 outcomes
results_bkmr_phthalates_taxa_t2_run5 <- 
  future_lapply(list(outcome_p1_t2, outcome_p2_t2, outcome_p3_t2, outcome_p4_t2, 
                     outcome_g1_t2, outcome_g2_t2, outcome_g3_t2, outcome_g4_t2), 
                bkmr_t2,
                future.seed = TRUE)
results_bkmr_phthalates_taxa_t3_run5 <- 
  future_lapply(list(outcome_p1_t3, outcome_p2_t3, outcome_p3_t3, outcome_p4_t3, 
                     outcome_g1_t3, outcome_g2_t3, outcome_g3_t3, outcome_g4_t3), 
                bkmr_t3, 
                future.seed = TRUE)
results_bkmr_phthalates_taxa_Y1_run5 <- 
  future_lapply(list(outcome_p1_Y1, outcome_p2_Y1, outcome_p3_Y1, outcome_p4_Y1, 
                     outcome_g1_Y1, outcome_g2_Y1, outcome_g3_Y1, outcome_g4_Y1), 
                bkmr_Y1, 
                future.seed = TRUE)

results_bkmr_phthalates_taxa_run5 <- 
  list(results_bkmr_phthalates_taxa_t2_run5, 
       results_bkmr_phthalates_taxa_t3_run5, 
       results_bkmr_phthalates_taxa_Y1_run5)
rm(results_bkmr_phthalates_taxa_t2_run5, results_bkmr_phthalates_taxa_t3_run5, results_bkmr_phthalates_taxa_Y1_run5)
names(results_bkmr_phthalates_taxa_run5) <- c("T2", "T3", "Y1")
names(results_bkmr_phthalates_taxa_run5$T2) <- c("p1", "p2", "p3", "p4", "g1", "g2", "g3", "g4")
names(results_bkmr_phthalates_taxa_run5$T3) <- c("p1", "p2", "p3", "p4", "g1", "g2", "g3", "g4")
names(results_bkmr_phthalates_taxa_run5$Y1) <- c("p1", "p2", "p3", "p4", "g1", "g2", "g3", "g4")

## Run 6 ----
# Spécifier le planificateur future et le nombre de noyaux à utiliser
plan(multisession)
ncores <- availableCores() # ou spécifier le nombre de noyaux à utiliser

# Utiliser future_lapply pour exécuter les appels bkmr en parallèle pour les 3 outcomes
results_bkmr_phthalates_taxa_t2_run6 <- 
  future_lapply(list(outcome_p1_t2, outcome_p2_t2, outcome_p3_t2, outcome_p4_t2, 
                     outcome_g1_t2, outcome_g2_t2, outcome_g3_t2, outcome_g4_t2), 
                bkmr_t2, 
                future.seed = TRUE)
results_bkmr_phthalates_taxa_t3_run6 <- 
  future_lapply(list(outcome_p1_t3, outcome_p2_t3, outcome_p3_t3, outcome_p4_t3, 
                     outcome_g1_t3, outcome_g2_t3, outcome_g3_t3, outcome_g4_t3), 
                bkmr_t3, 
                future.seed = TRUE)
results_bkmr_phthalates_taxa_Y1_run6 <- 
  future_lapply(list(outcome_p1_Y1, outcome_p2_Y1, outcome_p3_Y1, outcome_p4_Y1, 
                     outcome_g1_Y1, outcome_g2_Y1, outcome_g3_Y1, outcome_g4_Y1), 
                bkmr_Y1, 
                future.seed = TRUE)

results_bkmr_phthalates_taxa_run6 <- 
  list(results_bkmr_phthalates_taxa_t2_run6, 
       results_bkmr_phthalates_taxa_t3_run6, 
       results_bkmr_phthalates_taxa_Y1_run6)
rm(results_bkmr_phthalates_taxa_t2_run6, results_bkmr_phthalates_taxa_t3_run6, results_bkmr_phthalates_taxa_Y1_run6)
names(results_bkmr_phthalates_taxa_run6) <- c("T2", "T3", "Y1")
names(results_bkmr_phthalates_taxa_run6$T2) <- c("p1", "p2", "p3", "p4", "g1", "g2", "g3", "g4")
names(results_bkmr_phthalates_taxa_run6$T3) <- c("p1", "p2", "p3", "p4", "g1", "g2", "g3", "g4")
names(results_bkmr_phthalates_taxa_run6$Y1) <- c("p1", "p2", "p3", "p4", "g1", "g2", "g3", "g4")

# Model convergence ----
## Run 5 ----
TracePlot_group(results_bkmr_phthalates_taxa_run5$T2$p1, 
                results_bkmr_phthalates_taxa_run5$T2$p2, 
                results_bkmr_phthalates_taxa_run5$T2$p3,
                results_bkmr_phthalates_taxa_run5$T2$p4,
                titre = "Model convergence Phthalates BKMR t2, Firmicutes, Actinobacteria, Bacteroidetes, Proteobacteria (de haut en bas)")
TracePlot_group(results_bkmr_phthalates_taxa_run5$T3$p1, 
                results_bkmr_phthalates_taxa_run5$T3$p2, 
                results_bkmr_phthalates_taxa_run5$T3$p3,
                results_bkmr_phthalates_taxa_run5$T3$p4,
                titre = "Model convergence Phthalates BKMR t3, Firmicutes, Actinobacteria, Bacteroidetes, Proteobacteria (de haut en bas)")
TracePlot_group(results_bkmr_phthalates_taxa_run5$Y1$p1, 
                results_bkmr_phthalates_taxa_run5$Y1$p2, 
                results_bkmr_phthalates_taxa_run5$Y1$p3,
                results_bkmr_phthalates_taxa_run5$Y1$p4,
                titre = "Model convergence Phthalates BKMR Y1, Firmicutes, Actinobacteria, Bacteroidetes, Proteobacteria (de haut en bas)")
TracePlot_group(results_bkmr_phthalates_taxa_run5$T2$g1, 
                results_bkmr_phthalates_taxa_run5$T2$g2, 
                results_bkmr_phthalates_taxa_run5$T2$g3,
                results_bkmr_phthalates_taxa_run5$T2$g4,
                titre = "Model convergence Phthalates BKMR t2, Bifidobacterium, Bacteroides, Blautia, Escherichia and Shigella (de haut en bas)")
TracePlot_group(results_bkmr_phthalates_taxa_run5$T3$g1, 
                results_bkmr_phthalates_taxa_run5$T3$g2, 
                results_bkmr_phthalates_taxa_run5$T3$g3,
                results_bkmr_phthalates_taxa_run5$T3$g4,
                titre = "Model convergence Phthalates BKMR t3, Bifidobacterium, Bacteroides, Blautia, Escherichia and Shigella (de haut en bas)")
TracePlot_group(results_bkmr_phthalates_taxa_run5$Y1$g1, 
                results_bkmr_phthalates_taxa_run5$Y1$g2, 
                results_bkmr_phthalates_taxa_run5$Y1$g3,
                results_bkmr_phthalates_taxa_run5$Y1$g4,
                titre = "Model convergence Phthalates BKMR Y1, Bifidobacterium, Bacteroides, Blautia, Escherichia and Shigella (de haut en bas)")
par(mfrow=c(1,1))

## Run 6 ----
TracePlot_group(results_bkmr_phthalates_taxa_run6$T2$p1, 
                results_bkmr_phthalates_taxa_run6$T2$p2, 
                results_bkmr_phthalates_taxa_run6$T2$p3,
                results_bkmr_phthalates_taxa_run6$T2$p4,
                titre = "Model convergence Phthalates BKMR t2, Firmicutes, Actinobacteria, Bacteroidetes, Proteobacteria (de haut en bas)")
TracePlot_group(results_bkmr_phthalates_taxa_run6$T3$p1, 
                results_bkmr_phthalates_taxa_run6$T3$p2, 
                results_bkmr_phthalates_taxa_run6$T3$p3,
                results_bkmr_phthalates_taxa_run6$T3$p4,
                titre = "Model convergence Phthalates BKMR t3, Firmicutes, Actinobacteria, Bacteroidetes, Proteobacteria (de haut en bas)")
TracePlot_group(results_bkmr_phthalates_taxa_run6$Y1$p1, 
                results_bkmr_phthalates_taxa_run6$Y1$p2, 
                results_bkmr_phthalates_taxa_run6$Y1$p3,
                results_bkmr_phthalates_taxa_run6$Y1$p4,
                titre = "Model convergence Phthalates BKMR Y1, Firmicutes, Actinobacteria, Bacteroidetes, Proteobacteria (de haut en bas)")
TracePlot_group(results_bkmr_phthalates_taxa_run6$T2$g1, 
                results_bkmr_phthalates_taxa_run6$T2$g2, 
                results_bkmr_phthalates_taxa_run6$T2$g3,
                results_bkmr_phthalates_taxa_run6$T2$g4,
                titre = "Model convergence Phthalates BKMR t2, Bifidobacterium, Bacteroides, Blautia, Escherichia and Shigella (de haut en bas)")
TracePlot_group(results_bkmr_phthalates_taxa_run6$T3$g1, 
                results_bkmr_phthalates_taxa_run6$T3$g2, 
                results_bkmr_phthalates_taxa_run6$T3$g3,
                results_bkmr_phthalates_taxa_run6$T3$g4,
                titre = "Model convergence Phthalates BKMR t3, Bifidobacterium, Bacteroides, Blautia, Escherichia and Shigella (de haut en bas)")
TracePlot_group(results_bkmr_phthalates_taxa_run6$Y1$g1, 
                results_bkmr_phthalates_taxa_run6$Y1$g2, 
                results_bkmr_phthalates_taxa_run6$Y1$g3,
                results_bkmr_phthalates_taxa_run6$Y1$g4,
                titre = "Model convergence Phthalates BKMR Y1, Bifidobacterium, Bacteroides, Blautia, Escherichia and Shigella (de haut en bas)")
par(mfrow=c(1,1))


# PIP ----
## Run 5 ----
pip_phthalates_taxa_run5 <- list(
  T2 = pip_results(
    results_bkmr_phthalates_taxa_run5$T2$p1,
    results_bkmr_phthalates_taxa_run5$T2$p2,
    results_bkmr_phthalates_taxa_run5$T2$p3,
    results_bkmr_phthalates_taxa_run5$T2$p4,
    results_bkmr_phthalates_taxa_run5$T2$g1,
    results_bkmr_phthalates_taxa_run5$T2$g2,
    results_bkmr_phthalates_taxa_run5$T2$g3,
    results_bkmr_phthalates_taxa_run5$T2$g4), 
  T3 = pip_results(
    results_bkmr_phthalates_taxa_run5$T3$p1,
    results_bkmr_phthalates_taxa_run5$T3$p2,
    results_bkmr_phthalates_taxa_run5$T3$p3,
    results_bkmr_phthalates_taxa_run5$T3$p4,
    results_bkmr_phthalates_taxa_run5$T3$g1,
    results_bkmr_phthalates_taxa_run5$T3$g2,
    results_bkmr_phthalates_taxa_run5$T3$g3,
    results_bkmr_phthalates_taxa_run5$T3$g4),
  Y1 = pip_results(
    results_bkmr_phthalates_taxa_run5$Y1$p1,
    results_bkmr_phthalates_taxa_run5$Y1$p2,
    results_bkmr_phthalates_taxa_run5$Y1$p3,
    results_bkmr_phthalates_taxa_run5$Y1$p4,
    results_bkmr_phthalates_taxa_run5$Y1$g1,
    results_bkmr_phthalates_taxa_run5$Y1$g2,
    results_bkmr_phthalates_taxa_run5$Y1$g3,
    results_bkmr_phthalates_taxa_run5$Y1$g4)) %>%
  bind_rows()

## Run 6 ----
pip_phthalates_taxa_run6 <- list(
  T2 = pip_results(
    results_bkmr_phthalates_taxa_run6$T2$p1,
    results_bkmr_phthalates_taxa_run6$T2$p2,
    results_bkmr_phthalates_taxa_run6$T2$p3,
    results_bkmr_phthalates_taxa_run6$T2$p4,
    results_bkmr_phthalates_taxa_run6$T2$g1,
    results_bkmr_phthalates_taxa_run6$T2$g2,
    results_bkmr_phthalates_taxa_run6$T2$g3,
    results_bkmr_phthalates_taxa_run6$T2$g4), 
  T3 = pip_results(
    results_bkmr_phthalates_taxa_run6$T3$p1,
    results_bkmr_phthalates_taxa_run6$T3$p2,
    results_bkmr_phthalates_taxa_run6$T3$p3,
    results_bkmr_phthalates_taxa_run6$T3$p4,
    results_bkmr_phthalates_taxa_run6$T3$g1,
    results_bkmr_phthalates_taxa_run6$T3$g2,
    results_bkmr_phthalates_taxa_run6$T3$g3,
    results_bkmr_phthalates_taxa_run6$T3$g4),
  Y1 = pip_results(
    results_bkmr_phthalates_taxa_run6$Y1$p1,
    results_bkmr_phthalates_taxa_run6$Y1$p2,
    results_bkmr_phthalates_taxa_run6$Y1$p3,
    results_bkmr_phthalates_taxa_run6$Y1$p4,
    results_bkmr_phthalates_taxa_run6$Y1$g1,
    results_bkmr_phthalates_taxa_run6$Y1$g2,
    results_bkmr_phthalates_taxa_run6$Y1$g3,
    results_bkmr_phthalates_taxa_run6$Y1$g4)) %>% 
  bind_rows()

# Overall ----
## Run 5 ----
results_bkmr_overall_phthalates_taxa_run5 <- 
  list(
    T2 = list(
      p1 = risks_overall(results_bkmr_phthalates_taxa_run5$T2$p1, outcome_p1_t2, mixture_taxa_t2, covariates_taxa_t2),
      p2 = risks_overall(results_bkmr_phthalates_taxa_run5$T2$p2, outcome_p2_t2, mixture_taxa_t2, covariates_taxa_t2),
      p3 = risks_overall(results_bkmr_phthalates_taxa_run5$T2$p3, outcome_p3_t2, mixture_taxa_t2, covariates_taxa_t2),
      p4 = risks_overall(results_bkmr_phthalates_taxa_run5$T2$p4, outcome_p4_t2, mixture_taxa_t2, covariates_taxa_t2),
      g1 = risks_overall(results_bkmr_phthalates_taxa_run5$T2$g1, outcome_g1_t2, mixture_taxa_t2, covariates_taxa_t2),
      g2 = risks_overall(results_bkmr_phthalates_taxa_run5$T2$g2, outcome_g2_t2, mixture_taxa_t2, covariates_taxa_t2),
      g3 = risks_overall(results_bkmr_phthalates_taxa_run5$T2$g3, outcome_g3_t2, mixture_taxa_t2, covariates_taxa_t2),
      g4 = risks_overall(results_bkmr_phthalates_taxa_run5$T2$g4, outcome_g4_t2, mixture_taxa_t2, covariates_taxa_t2)), 
    T3 = list(
      p1 = risks_overall(results_bkmr_phthalates_taxa_run5$T3$p1, outcome_p1_t3, mixture_taxa_t3, covariates_taxa_t3),
      p2 = risks_overall(results_bkmr_phthalates_taxa_run5$T3$p2, outcome_p2_t3, mixture_taxa_t3, covariates_taxa_t3),
      p3 = risks_overall(results_bkmr_phthalates_taxa_run5$T3$p3, outcome_p3_t3, mixture_taxa_t3, covariates_taxa_t3),
      p4 = risks_overall(results_bkmr_phthalates_taxa_run5$T3$p4, outcome_p4_t3, mixture_taxa_t3, covariates_taxa_t3),
      g1 = risks_overall(results_bkmr_phthalates_taxa_run5$T3$g1, outcome_g1_t3, mixture_taxa_t3, covariates_taxa_t3),
      g2 = risks_overall(results_bkmr_phthalates_taxa_run5$T3$g2, outcome_g2_t3, mixture_taxa_t3, covariates_taxa_t3),
      g3 = risks_overall(results_bkmr_phthalates_taxa_run5$T3$g3, outcome_g3_t3, mixture_taxa_t3, covariates_taxa_t3),
      g4 = risks_overall(results_bkmr_phthalates_taxa_run5$T3$g4, outcome_g4_t3, mixture_taxa_t3, covariates_taxa_t3)),
    Y1 = list(
      p1 = risks_overall(results_bkmr_phthalates_taxa_run5$Y1$p1, outcome_p1_Y1, mixture_taxa_Y1, covariates_taxa_Y1),
      p2 = risks_overall(results_bkmr_phthalates_taxa_run5$Y1$p2, outcome_p2_Y1, mixture_taxa_Y1, covariates_taxa_Y1),
      p3 = risks_overall(results_bkmr_phthalates_taxa_run5$Y1$p3, outcome_p3_Y1, mixture_taxa_Y1, covariates_taxa_Y1),
      p4 = risks_overall(results_bkmr_phthalates_taxa_run5$Y1$p4, outcome_p4_Y1, mixture_taxa_Y1, covariates_taxa_Y1),
      g1 = risks_overall(results_bkmr_phthalates_taxa_run5$Y1$g1, outcome_g1_Y1, mixture_taxa_Y1, covariates_taxa_Y1),
      g2 = risks_overall(results_bkmr_phthalates_taxa_run5$Y1$g2, outcome_g2_Y1, mixture_taxa_Y1, covariates_taxa_Y1),
      g3 = risks_overall(results_bkmr_phthalates_taxa_run5$Y1$g3, outcome_g3_Y1, mixture_taxa_Y1, covariates_taxa_Y1),
      g4 = risks_overall(results_bkmr_phthalates_taxa_run5$Y1$g4, outcome_g4_Y1, mixture_taxa_Y1, covariates_taxa_Y1)))

overall_phthalates_taxa_run5 <- bind_rows(
  list_t2 = do.call(rbind, results_bkmr_overall_phthalates_taxa_run5$T2) %>% 
    as.data.frame() %>% 
    rownames_to_column("genus") %>%
    mutate(genus = gsub("\\..$", "", genus), 
           window = "T2"), 
  list_t3 = do.call(rbind, results_bkmr_overall_phthalates_taxa_run5$T3) %>% 
    as.data.frame() %>% 
    rownames_to_column("genus") %>%
    mutate(genus = gsub("\\..$", "", genus), 
           window = "T3"),
  list_Y1 = do.call(rbind, results_bkmr_overall_phthalates_taxa_run5$Y1) %>% 
    as.data.frame() %>% 
    rownames_to_column("genus") %>%
    mutate(genus = gsub("\\..$", "", genus), 
           window = "Y1")) %>%
  select(window, everything())


## Run 6 ----
results_bkmr_overall_phthalates_taxa_run6 <- 
  list(
    T2 = list(
      p1 = risks_overall(results_bkmr_phthalates_taxa_run6$T2$p1, outcome_p1_t2, mixture_taxa_t2, covariates_taxa_t2),
      p2 = risks_overall(results_bkmr_phthalates_taxa_run6$T2$p2, outcome_p2_t2, mixture_taxa_t2, covariates_taxa_t2),
      p3 = risks_overall(results_bkmr_phthalates_taxa_run6$T2$p3, outcome_p3_t2, mixture_taxa_t2, covariates_taxa_t2),
      p4 = risks_overall(results_bkmr_phthalates_taxa_run6$T2$p4, outcome_p4_t2, mixture_taxa_t2, covariates_taxa_t2),
      g1 = risks_overall(results_bkmr_phthalates_taxa_run6$T2$g1, outcome_g1_t2, mixture_taxa_t2, covariates_taxa_t2),
      g2 = risks_overall(results_bkmr_phthalates_taxa_run6$T2$g2, outcome_g2_t2, mixture_taxa_t2, covariates_taxa_t2),
      g3 = risks_overall(results_bkmr_phthalates_taxa_run6$T2$g3, outcome_g3_t2, mixture_taxa_t2, covariates_taxa_t2),
      g4 = risks_overall(results_bkmr_phthalates_taxa_run6$T2$g4, outcome_g4_t2, mixture_taxa_t2, covariates_taxa_t2)), 
    T3 = list(
      p1 = risks_overall(results_bkmr_phthalates_taxa_run6$T3$p1, outcome_p1_t3, mixture_taxa_t3, covariates_taxa_t3),
      p2 = risks_overall(results_bkmr_phthalates_taxa_run6$T3$p2, outcome_p2_t3, mixture_taxa_t3, covariates_taxa_t3),
      p3 = risks_overall(results_bkmr_phthalates_taxa_run6$T3$p3, outcome_p3_t3, mixture_taxa_t3, covariates_taxa_t3),
      p4 = risks_overall(results_bkmr_phthalates_taxa_run6$T3$p4, outcome_p4_t3, mixture_taxa_t3, covariates_taxa_t3),
      g1 = risks_overall(results_bkmr_phthalates_taxa_run6$T3$g1, outcome_g1_t3, mixture_taxa_t3, covariates_taxa_t3),
      g2 = risks_overall(results_bkmr_phthalates_taxa_run6$T3$g2, outcome_g2_t3, mixture_taxa_t3, covariates_taxa_t3),
      g3 = risks_overall(results_bkmr_phthalates_taxa_run6$T3$g3, outcome_g3_t3, mixture_taxa_t3, covariates_taxa_t3),
      g4 = risks_overall(results_bkmr_phthalates_taxa_run6$T3$g4, outcome_g4_t3, mixture_taxa_t3, covariates_taxa_t3)),
    Y1 = list(
      p1 = risks_overall(results_bkmr_phthalates_taxa_run6$Y1$p1, outcome_p1_Y1, mixture_taxa_Y1, covariates_taxa_Y1),
      p2 = risks_overall(results_bkmr_phthalates_taxa_run6$Y1$p2, outcome_p2_Y1, mixture_taxa_Y1, covariates_taxa_Y1),
      p3 = risks_overall(results_bkmr_phthalates_taxa_run6$Y1$p3, outcome_p3_Y1, mixture_taxa_Y1, covariates_taxa_Y1),
      p4 = risks_overall(results_bkmr_phthalates_taxa_run6$Y1$p4, outcome_p4_Y1, mixture_taxa_Y1, covariates_taxa_Y1),
      g1 = risks_overall(results_bkmr_phthalates_taxa_run6$Y1$g1, outcome_g1_Y1, mixture_taxa_Y1, covariates_taxa_Y1),
      g2 = risks_overall(results_bkmr_phthalates_taxa_run6$Y1$g2, outcome_g2_Y1, mixture_taxa_Y1, covariates_taxa_Y1),
      g3 = risks_overall(results_bkmr_phthalates_taxa_run6$Y1$g3, outcome_g3_Y1, mixture_taxa_Y1, covariates_taxa_Y1),
      g4 = risks_overall(results_bkmr_phthalates_taxa_run6$Y1$g4, outcome_g4_Y1, mixture_taxa_Y1, covariates_taxa_Y1)))

overall_phthalates_taxa_run6 <- bind_rows(
  list_t2 = do.call(rbind, results_bkmr_overall_phthalates_taxa_run6$T2) %>% 
    as.data.frame() %>% 
    rownames_to_column("genus") %>%
    mutate(genus = gsub("\\..$", "", genus), 
           window = "T2"), 
  list_t3 = do.call(rbind, results_bkmr_overall_phthalates_taxa_run6$T3) %>% 
    as.data.frame() %>% 
    rownames_to_column("genus") %>%
    mutate(genus = gsub("\\..$", "", genus), 
           window = "T3"),
  list_Y1 = do.call(rbind, results_bkmr_overall_phthalates_taxa_run6$Y1) %>% 
    as.data.frame() %>% 
    rownames_to_column("genus") %>%
    mutate(genus = gsub("\\..$", "", genus), 
           window = "Y1")) %>%
  select(window, everything())

# Singvar ----
## Run 5 ----
results_bkmr_singvar_phthalates_taxa_run5 <- 
  list(
    T2 = list(
      p1 = risks_singvar(results_bkmr_phthalates_taxa_run5$T2$p1, outcome_p1_t2, mixture_taxa_t2, covariates_taxa_t2),
      p2 = risks_singvar(results_bkmr_phthalates_taxa_run5$T2$p2, outcome_p2_t2, mixture_taxa_t2, covariates_taxa_t2),
      p3 = risks_singvar(results_bkmr_phthalates_taxa_run5$T2$p3, outcome_p3_t2, mixture_taxa_t2, covariates_taxa_t2),
      p4 = risks_singvar(results_bkmr_phthalates_taxa_run5$T2$p4, outcome_p4_t2, mixture_taxa_t2, covariates_taxa_t2),
      g1 = risks_singvar(results_bkmr_phthalates_taxa_run5$T2$g1, outcome_g1_t2, mixture_taxa_t2, covariates_taxa_t2),
      g2 = risks_singvar(results_bkmr_phthalates_taxa_run5$T2$g2, outcome_g2_t2, mixture_taxa_t2, covariates_taxa_t2),
      g3 = risks_singvar(results_bkmr_phthalates_taxa_run5$T2$g3, outcome_g3_t2, mixture_taxa_t2, covariates_taxa_t2),
      g4 = risks_singvar(results_bkmr_phthalates_taxa_run5$T2$g4, outcome_g4_t2, mixture_taxa_t2, covariates_taxa_t2)), 
    T3 = list(
      p1 = risks_singvar(results_bkmr_phthalates_taxa_run5$T3$p1, outcome_p1_t3, mixture_taxa_t3, covariates_taxa_t3),
      p2 = risks_singvar(results_bkmr_phthalates_taxa_run5$T3$p2, outcome_p2_t3, mixture_taxa_t3, covariates_taxa_t3),
      p3 = risks_singvar(results_bkmr_phthalates_taxa_run5$T3$p3, outcome_p3_t3, mixture_taxa_t3, covariates_taxa_t3),
      p4 = risks_singvar(results_bkmr_phthalates_taxa_run5$T3$p4, outcome_p4_t3, mixture_taxa_t3, covariates_taxa_t3),
      g1 = risks_singvar(results_bkmr_phthalates_taxa_run5$T3$g1, outcome_g1_t3, mixture_taxa_t3, covariates_taxa_t3),
      g2 = risks_singvar(results_bkmr_phthalates_taxa_run5$T3$g2, outcome_g2_t3, mixture_taxa_t3, covariates_taxa_t3),
      g3 = risks_singvar(results_bkmr_phthalates_taxa_run5$T3$g3, outcome_g3_t3, mixture_taxa_t3, covariates_taxa_t3),
      g4 = risks_singvar(results_bkmr_phthalates_taxa_run5$T3$g4, outcome_g4_t3, mixture_taxa_t3, covariates_taxa_t3)),
    Y1 = list(
      p1 = risks_singvar(results_bkmr_phthalates_taxa_run5$Y1$p1, outcome_p1_Y1, mixture_taxa_Y1, covariates_taxa_Y1),
      p2 = risks_singvar(results_bkmr_phthalates_taxa_run5$Y1$p2, outcome_p2_Y1, mixture_taxa_Y1, covariates_taxa_Y1),
      p3 = risks_singvar(results_bkmr_phthalates_taxa_run5$Y1$p3, outcome_p3_Y1, mixture_taxa_Y1, covariates_taxa_Y1),
      p4 = risks_singvar(results_bkmr_phthalates_taxa_run5$Y1$p4, outcome_p4_Y1, mixture_taxa_Y1, covariates_taxa_Y1),
      g1 = risks_singvar(results_bkmr_phthalates_taxa_run5$Y1$g1, outcome_g1_Y1, mixture_taxa_Y1, covariates_taxa_Y1),
      g2 = risks_singvar(results_bkmr_phthalates_taxa_run5$Y1$g2, outcome_g2_Y1, mixture_taxa_Y1, covariates_taxa_Y1),
      g3 = risks_singvar(results_bkmr_phthalates_taxa_run5$Y1$g3, outcome_g3_Y1, mixture_taxa_Y1, covariates_taxa_Y1),
      g4 = risks_singvar(results_bkmr_phthalates_taxa_run5$Y1$g4, outcome_g4_Y1, mixture_taxa_Y1, covariates_taxa_Y1)))

valeurs <- c("p1", "p2", "p3", "p4", "g1", "g2", "g3", "g4")
for (i in 1:8) {
  results_bkmr_singvar_phthalates_taxa_run5$T2[[i]] <- results_bkmr_singvar_phthalates_taxa_run5$T2[[i]] %>% mutate(taxa = valeurs[i])}
for (i in 1:8) {
  results_bkmr_singvar_phthalates_taxa_run5$T3[[i]] <- results_bkmr_singvar_phthalates_taxa_run5$T3[[i]] %>% mutate(taxa = valeurs[i])}
for (i in 1:8) {
  results_bkmr_singvar_phthalates_taxa_run5$Y1[[i]] <- results_bkmr_singvar_phthalates_taxa_run5$Y1[[i]] %>% mutate(taxa = valeurs[i])}
rm(valeurs)

singvar_phthalates_taxa_run5 <- bind_rows(
  list_t2 = do.call(rbind, results_bkmr_singvar_phthalates_taxa_run5$T2) %>% 
    as.data.frame() %>% 
    mutate(window = "T2"), 
  list_t3 = do.call(rbind, results_bkmr_singvar_phthalates_taxa_run5$T3) %>% 
    as.data.frame() %>% 
    mutate(window = "T3"),
  list_Y1 = do.call(rbind, results_bkmr_singvar_phthalates_taxa_run5$Y1) %>% 
    as.data.frame() %>% 
    mutate(window = "Y1"))  %>%
  select(window, taxa, variable, everything())%>%
  rename(pollutant = variable) %>%
  arrange(window, taxa, pollutant)%>%
  mutate(
    pollutant = str_replace_all(pollutant, 
                                c("_i_cor_" = " ", 
                                  "mo_" = "", 
                                  "_ln" = "", 
                                  "ch_" = "", 
                                  "_ms" = "", 
                                  "DEHP" = "ΣDEHP", 
                                  "DiNP" = "ΣDiNP", 
                                  "DINCH" = "ΣDINCH")), 
    pollutant = fct_relevel(pollutant,
                            "ΣDINCH Y1", "ΣDINCH t3", "ΣDINCH t2", 
                            "ohMPHP Y1", "ohMPHP t3", "ohMPHP t2", 
                            "MEP Y1", "MEP t3", "MEP t2", 
                            "MBzP Y1", "MBzP t3", "MBzP t2", 
                            "MiBP Y1", "MiBP t3", "MiBP t2",
                            "MnBP Y1", "MnBP t3", "MnBP t2", 
                            "ΣDiNP Y1", "ΣDiNP t3", "ΣDiNP t2", 
                            "ΣDEHP Y1", "ΣDEHP t3", "ΣDEHP t2"))

## Run 6 ----
results_bkmr_singvar_phthalates_taxa_run6 <- 
  list(
    T2 = list(
      p1 = risks_singvar(results_bkmr_phthalates_taxa_run6$T2$p1, outcome_p1_t2, mixture_taxa_t2, covariates_taxa_t2),
      p2 = risks_singvar(results_bkmr_phthalates_taxa_run6$T2$p2, outcome_p2_t2, mixture_taxa_t2, covariates_taxa_t2),
      p3 = risks_singvar(results_bkmr_phthalates_taxa_run6$T2$p3, outcome_p3_t2, mixture_taxa_t2, covariates_taxa_t2),
      p4 = risks_singvar(results_bkmr_phthalates_taxa_run6$T2$p4, outcome_p4_t2, mixture_taxa_t2, covariates_taxa_t2),
      g1 = risks_singvar(results_bkmr_phthalates_taxa_run6$T2$g1, outcome_g1_t2, mixture_taxa_t2, covariates_taxa_t2),
      g2 = risks_singvar(results_bkmr_phthalates_taxa_run6$T2$g2, outcome_g2_t2, mixture_taxa_t2, covariates_taxa_t2),
      g3 = risks_singvar(results_bkmr_phthalates_taxa_run6$T2$g3, outcome_g3_t2, mixture_taxa_t2, covariates_taxa_t2),
      g4 = risks_singvar(results_bkmr_phthalates_taxa_run6$T2$g4, outcome_g4_t2, mixture_taxa_t2, covariates_taxa_t2)), 
    T3 = list(
      p1 = risks_singvar(results_bkmr_phthalates_taxa_run6$T3$p1, outcome_p1_t3, mixture_taxa_t3, covariates_taxa_t3),
      p2 = risks_singvar(results_bkmr_phthalates_taxa_run6$T3$p2, outcome_p2_t3, mixture_taxa_t3, covariates_taxa_t3),
      p3 = risks_singvar(results_bkmr_phthalates_taxa_run6$T3$p3, outcome_p3_t3, mixture_taxa_t3, covariates_taxa_t3),
      p4 = risks_singvar(results_bkmr_phthalates_taxa_run6$T3$p4, outcome_p4_t3, mixture_taxa_t3, covariates_taxa_t3),
      g1 = risks_singvar(results_bkmr_phthalates_taxa_run6$T3$g1, outcome_g1_t3, mixture_taxa_t3, covariates_taxa_t3),
      g2 = risks_singvar(results_bkmr_phthalates_taxa_run6$T3$g2, outcome_g2_t3, mixture_taxa_t3, covariates_taxa_t3),
      g3 = risks_singvar(results_bkmr_phthalates_taxa_run6$T3$g3, outcome_g3_t3, mixture_taxa_t3, covariates_taxa_t3),
      g4 = risks_singvar(results_bkmr_phthalates_taxa_run6$T3$g4, outcome_g4_t3, mixture_taxa_t3, covariates_taxa_t3)),
    Y1 = list(
      p1 = risks_singvar(results_bkmr_phthalates_taxa_run6$Y1$p1, outcome_p1_Y1, mixture_taxa_Y1, covariates_taxa_Y1),
      p2 = risks_singvar(results_bkmr_phthalates_taxa_run6$Y1$p2, outcome_p2_Y1, mixture_taxa_Y1, covariates_taxa_Y1),
      p3 = risks_singvar(results_bkmr_phthalates_taxa_run6$Y1$p3, outcome_p3_Y1, mixture_taxa_Y1, covariates_taxa_Y1),
      p4 = risks_singvar(results_bkmr_phthalates_taxa_run6$Y1$p4, outcome_p4_Y1, mixture_taxa_Y1, covariates_taxa_Y1),
      g1 = risks_singvar(results_bkmr_phthalates_taxa_run6$Y1$g1, outcome_g1_Y1, mixture_taxa_Y1, covariates_taxa_Y1),
      g2 = risks_singvar(results_bkmr_phthalates_taxa_run6$Y1$g2, outcome_g2_Y1, mixture_taxa_Y1, covariates_taxa_Y1),
      g3 = risks_singvar(results_bkmr_phthalates_taxa_run6$Y1$g3, outcome_g3_Y1, mixture_taxa_Y1, covariates_taxa_Y1),
      g4 = risks_singvar(results_bkmr_phthalates_taxa_run6$Y1$g4, outcome_g4_Y1, mixture_taxa_Y1, covariates_taxa_Y1)))

valeurs <- c("p1", "p2", "p3", "p4", "g1", "g2", "g3", "g4")
for (i in 1:8) {
  results_bkmr_singvar_phthalates_taxa_run6$T2[[i]] <- results_bkmr_singvar_phthalates_taxa_run6$T2[[i]] %>% mutate(taxa = valeurs[i])}
for (i in 1:8) {
  results_bkmr_singvar_phthalates_taxa_run6$T3[[i]] <- results_bkmr_singvar_phthalates_taxa_run6$T3[[i]] %>% mutate(taxa = valeurs[i])}
for (i in 1:8) {
  results_bkmr_singvar_phthalates_taxa_run6$Y1[[i]] <- results_bkmr_singvar_phthalates_taxa_run6$Y1[[i]] %>% mutate(taxa = valeurs[i])}
rm(valeurs)

singvar_phthalates_taxa_run6 <- bind_rows(
  list_t2 = do.call(rbind, results_bkmr_singvar_phthalates_taxa_run6$T2) %>% 
    as.data.frame() %>% 
    mutate(window = "T2"), 
  list_t3 = do.call(rbind, results_bkmr_singvar_phthalates_taxa_run6$T3) %>% 
    as.data.frame() %>% 
    mutate(window = "T3"),
  list_Y1 = do.call(rbind, results_bkmr_singvar_phthalates_taxa_run6$Y1) %>% 
    as.data.frame() %>% 
    mutate(window = "Y1"))  %>%
  select(window, taxa, variable, everything())%>%
  rename(pollutant = variable) %>%
  arrange(window, taxa, pollutant)%>%
  mutate(
    pollutant = str_replace_all(pollutant, 
                                c("_i_cor_" = " ", 
                                  "mo_" = "", 
                                  "_ln" = "", 
                                  "ch_" = "", 
                                  "_ms" = "", 
                                  "DEHP" = "ΣDEHP", 
                                  "DiNP" = "ΣDiNP", 
                                  "DINCH" = "ΣDINCH")), 
    pollutant = fct_relevel(pollutant,
                            "ΣDINCH Y1", "ΣDINCH t3", "ΣDINCH t2", 
                            "ohMPHP Y1", "ohMPHP t3", "ohMPHP t2", 
                            "MEP Y1", "MEP t3", "MEP t2", 
                            "MBzP Y1", "MBzP t3", "MBzP t2", 
                            "MiBP Y1", "MiBP t3", "MiBP t2",
                            "MnBP Y1", "MnBP t3", "MnBP t2", 
                            "ΣDiNP Y1", "ΣDiNP t3", "ΣDiNP t2", 
                            "ΣDEHP Y1", "ΣDEHP t3", "ΣDEHP t2"))

# Assemblage et export ----
write_xlsx(
  list(                           
    pip_phthalates_taxa_run5 = pip_phthalates_taxa_run5,
    overall_phthalates_taxa_run5 = overall_phthalates_taxa_run5, 
    singvar_phthalates_taxa_run5 = singvar_phthalates_taxa_run5), 
  path = "4_output/bkmr/Run5 (ms taxa)/results_bkmr_phthalates_taxa_run5ms.xlsx")

write_xlsx(
  list(                           
    pip_phthalates_taxa_run6 = pip_phthalates_taxa_run6,
    overall_phthalates_taxa_run6 = overall_phthalates_taxa_run6, 
    singvar_phthalates_taxa_run6 = singvar_phthalates_taxa_run6), 
  path = "4_output/bkmr/Run5 (ms taxa)/results_bkmr_phthalates_taxa_run6ms.xlsx")



# Plot overall ----
## Run 5 ----
dev.off()
plot_risks.overall_phthalates_taxa_p_run5 <- 
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run5$T2$p1, title = "Phylum Firmicutes", y_title = bquote("2"^{nd}~trim.~exposure), x_title = "") +
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run5$T2$p2, title = "Phylum Actinobacteria", y_title = "", x_title = "") + 
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run5$T2$p3, title = "Phylum Bacteroidetes", y_title = "", x_title = "") + 
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run5$T2$p4, title = "Phylum Proteobacteria", y_title = "", x_title = "") + 
  
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run5$T3$p1, title = "", y_title = bquote("3"^{rd}~trim.~exposure), x_title = "") +
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run5$T3$p2, title = "", y_title = "", x_title = "") + 
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run5$T3$p3, title = "", y_title = "", x_title = "") + 
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run5$T3$p4, title = "", y_title = "", x_title = "") + 
  
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run5$Y1$p1, title = "", y_title = "12-month exposure", x_title = "quantile") +
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run5$Y1$p2, title = "", y_title = "", x_title = "quantile") + 
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run5$Y1$p3, title = "", y_title = "", x_title = "quantile") + 
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run5$Y1$p4, title = "", y_title = "", x_title = "quantile") +   
  
  plot_layout(ncol = 4, nrow = 3) 

plot_risks.overall_phthalates_taxa_g_run5 <- 
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run5$T2$g3, title = "Genus *Blautia*", y_title = bquote("2"^{nd}~trim.~exposure), x_title = "") +  # on met g3 en 1er pour garder l'ordre des phyla
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run5$T2$g1, title = "Genus *Bifidobacterium*", y_title = "", x_title = "") + 
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run5$T2$g2, title = "Genus *Bacteroides*", y_title = "", x_title = "") + 
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run5$T2$g4, title = "Genera *Escherichia*<br>and *Shigella*", y_title = "", x_title = "") + 
  
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run5$T3$g3, title = "", y_title = bquote("3"^{rd}~trim.~exposure), x_title = "") +  # on met g3 en 1er pour garder l'ordre des phyla
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run5$T3$g1, title = "", y_title = "", x_title = "") + 
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run5$T3$g2, title = "", y_title = "", x_title = "") + 
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run5$T3$g4, title = "", y_title = "", x_title = "") + 
  
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run5$Y1$g3, title = "", y_title = "12-month exposure", x_title = "quantile") +  # on met g3 en 1er pour garder l'ordre des phyla
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run5$Y1$g1, title = "", y_title = "", x_title = "quantile") + 
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run5$Y1$g2, title = "", y_title = "", x_title = "quantile") + 
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run5$Y1$g4, title = "", y_title = "", x_title = "quantile") +  
  
  plot_layout(ncol = 4, nrow = 3) 

ggsave("4_output/bkmr/Run5 (ms taxa)/plot_risks.overall_phthalates_taxa_p_run5ms.tiff", 
       plot_risks.overall_phthalates_taxa_p_run5, 
       device = "tiff",
       units = "mm",
       dpi = 300, 
       width = 300,
       height = 220)

ggsave("4_output/bkmr/Run5 (ms taxa)/plot_risks.overall_phthalates_taxa_g_run5ms.tiff", 
       plot_risks.overall_phthalates_taxa_g_run5, 
       device = "tiff",
       units = "mm",
       dpi = 300, 
       width = 300,
       height = 220)


## Run 6 ----
plot_risks.overall_phthalates_taxa_p_run6 <- 
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run6$T2$p1, title = "Phylum Firmicutes", y_title = bquote("2"^{nd}~trim.~exposure), x_title = "") +
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run6$T2$p2, title = "Phylum Actinobacteria", y_title = "", x_title = "") + 
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run6$T2$p3, title = "Phylum Bacteroidetes", y_title = "", x_title = "") + 
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run6$T2$p4, title = "Phylum Proteobacteria", y_title = "", x_title = "") + 
  
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run6$T3$p1, title = "", y_title = bquote("3"^{rd}~trim.~exposure), x_title = "") +
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run6$T3$p2, title = "", y_title = "", x_title = "") + 
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run6$T3$p3, title = "", y_title = "", x_title = "") + 
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run6$T3$p4, title = "", y_title = "", x_title = "") + 
  
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run6$Y1$p1, title = "", y_title = "12-month exposure", x_title = "quantile") +
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run6$Y1$p2, title = "", y_title = "", x_title = "quantile") + 
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run6$Y1$p3, title = "", y_title = "", x_title = "quantile") + 
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run6$Y1$p4, title = "", y_title = "", x_title = "quantile") +   
  
  plot_layout(ncol = 4, nrow = 3) 

plot_risks.overall_phthalates_taxa_g_run6 <- 
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run6$T2$g3, title = "Genus *Blautia*", y_title = bquote("2"^{nd}~trim.~exposure), x_title = "") +  # on met g3 en 1er pour garder l'ordre des phyla
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run6$T2$g1, title = "Genus *Bifidobacterium*", y_title = "", x_title = "") + 
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run6$T2$g2, title = "Genus *Bacteroides*", y_title = "", x_title = "") + 
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run6$T2$g4, title = "Genera *Escherichia*<br>and *Shigella*", y_title = "", x_title = "") + 
  
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run6$T3$g3, title = "", y_title = bquote("3"^{rd}~trim.~exposure), x_title = "") +  # on met g3 en 1er pour garder l'ordre des phyla
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run6$T3$g1, title = "", y_title = "", x_title = "") + 
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run6$T3$g2, title = "", y_title = "", x_title = "") + 
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run6$T3$g4, title = "", y_title = "", x_title = "") + 
  
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run6$Y1$g3, title = "", y_title = "12-month exposure", x_title = "quantile") +  # on met g3 en 1er pour garder l'ordre des phyla
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run6$Y1$g1, title = "", y_title = "", x_title = "quantile") + 
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run6$Y1$g2, title = "", y_title = "", x_title = "quantile") + 
  plot_risks.overall(results_bkmr_overall_phthalates_taxa_run6$Y1$g4, title = "", y_title = "", x_title = "quantile") +  
  
  plot_layout(ncol = 4, nrow = 3) 

ggsave("4_output/bkmr/run6 (ms taxa)/plot_risks.overall_phthalates_taxa_p_run6ms.tiff", 
       plot_risks.overall_phthalates_taxa_p_run6, 
       device = "tiff",
       units = "mm",
       dpi = 300, 
       width = 300,
       height = 220)

ggsave("4_output/bkmr/run6 (ms taxa)/plot_risks.overall_phthalates_taxa_g_run6ms.tiff", 
       plot_risks.overall_phthalates_taxa_g_run6, 
       device = "tiff",
       units = "mm",
       dpi = 300, 
       width = 300,
       height = 220)

# Plot Singvar ----
## Run 5 ----
plot_risks.singvar_phthalates_taxa_p_run5 <- 
  plot_risks.singvar(singvar_phthalates_taxa_run5, window = "T2", taxa = "p1", 
                     title = "Phylum Firmicutes", x_title = bquote("2"^{nd}~trim.~exposure), y_title = "", 
                     legend.position = "none", axis.y = element_text(size = 12)) +
  plot_risks.singvar(singvar_phthalates_taxa_run5, window = "T2", taxa = "p2",  
                     title = "Phylum Actinobacteria", x_title = "", y_title = "", 
                     legend.position = "none", axis.y = element_blank()) + 
  plot_risks.singvar(singvar_phthalates_taxa_run5, window = "T2", taxa = "p3",  
                     title = "Phylum Bacteroidetes", x_title = "", y_title = "", 
                     legend.position = "none", axis.y = element_blank()) + 
  plot_risks.singvar(singvar_phthalates_taxa_run5, window = "T2", taxa = "p4",  
                     title = "Phylum Proteobacteria", x_title = "", y_title = "", 
                     legend.position = "right", axis.y = element_blank()) + 
  
  plot_risks.singvar(singvar_phthalates_taxa_run5, window = "T3", taxa = "p1", 
                     title = "", x_title = bquote("3"^{rd}~trim.~exposure), y_title = "", 
                     legend.position = "none", axis.y = element_text(size = 12)) +
  plot_risks.singvar(singvar_phthalates_taxa_run5, window = "T3", taxa = "p2",
                     title = "", x_title = "", y_title = "", 
                     legend.position = "none", axis.y = element_blank()) + 
  plot_risks.singvar(singvar_phthalates_taxa_run5, window = "T3", taxa = "p3", 
                     title = "", x_title = "", y_title = "", 
                     legend.position = "none", axis.y = element_blank()) + 
  plot_risks.singvar(singvar_phthalates_taxa_run5, window = "T3", taxa = "p4", 
                     title = "", x_title = "", y_title = "", 
                     legend.position = "right", axis.y = element_blank()) + 
  
  plot_risks.singvar(singvar_phthalates_taxa_run5, window = "Y1", taxa = "p1", 
                     title = "", x_title = "12-month exposure", y_title = "Relative abundance (%)", 
                     legend.position = "none", axis.y = element_text(size = 12)) +
  plot_risks.singvar(singvar_phthalates_taxa_run5, window = "Y1", taxa = "p2", 
                     title = "", x_title = "", y_title = "Relative abundance (%)", 
                     legend.position = "none", axis.y = element_blank()) + 
  plot_risks.singvar(singvar_phthalates_taxa_run5, window = "Y1", taxa = "p3", 
                     title = "", x_title = "", y_title = "Relative abundance (%)", 
                     legend.position = "none", axis.y = element_blank()) + 
  plot_risks.singvar(singvar_phthalates_taxa_run5, window = "Y1", taxa = "p4", 
                     title = "", x_title = "", y_title = "Relative abundance (%)", 
                     legend.position = "right", axis.y = element_blank()) +   
  
  plot_layout(ncol = 4, nrow = 3) 

plot_risks.singvar_phthalates_taxa_g_run5 <- 
  plot_risks.singvar(singvar_phthalates_taxa_run5, window = "T2", taxa = "g3",  
                     title = "Genus *Blautia*", x_title = bquote("2"^{nd}~trim.~exposure),  y_title = "", 
                     legend.position = "none", axis.y = element_text(size = 12)) + 
  plot_risks.singvar(singvar_phthalates_taxa_run5, window = "T2", taxa = "g1", 
                     title = "Genus *Bifidobacterium*", x_title = "", y_title = "", 
                     legend.position = "none", axis.y = element_blank()) +
  plot_risks.singvar(singvar_phthalates_taxa_run5, window = "T2", taxa = "g2",  
                     title = "Genus *Bacteroides*", x_title = "", y_title = "", 
                     legend.position = "none", axis.y = element_blank()) + 
  plot_risks.singvar(singvar_phthalates_taxa_run5, window = "T2", taxa = "g4",  
                     title = "Genera *Escherichia*<br>and *Shigella*", x_title = "", y_title = "", 
                     legend.position = "right", axis.y = element_blank()) + 
  
  plot_risks.singvar(singvar_phthalates_taxa_run5, window = "T3", taxa = "g3", 
                     title = "", x_title = bquote("3"^{rd}~trim.~exposure), y_title = "", 
                     legend.position = "none", axis.y = element_text(size = 12)) + 
  plot_risks.singvar(singvar_phthalates_taxa_run5, window = "T3", taxa = "g1", 
                     title = "", x_title = "", y_title = "", 
                     legend.position = "none", axis.y = element_blank()) +
  plot_risks.singvar(singvar_phthalates_taxa_run5, window = "T3", taxa = "g2",
                     title = "", x_title = "", y_title = "", 
                     legend.position = "none", axis.y = element_blank()) + 
  plot_risks.singvar(singvar_phthalates_taxa_run5, window = "T3", taxa = "g4", 
                     title = "", x_title = "", y_title = "", 
                     legend.position = "right", axis.y = element_blank()) + 
  
  plot_risks.singvar(singvar_phthalates_taxa_run5, window = "Y1", taxa = "g3", 
                     title = "", x_title = "12-month exposure", y_title = "Relative abundance (%)", 
                     legend.position = "none", axis.y = element_text(size = 12)) + 
  plot_risks.singvar(singvar_phthalates_taxa_run5, window = "Y1", taxa = "g1", 
                     title = "", x_title = "", y_title = "Relative abundance (%)", 
                     legend.position = "none", axis.y = element_blank()) +
  plot_risks.singvar(singvar_phthalates_taxa_run5, window = "Y1", taxa = "g2", 
                     title = "", x_title = "", y_title = "Relative abundance (%)", 
                     legend.position = "none", axis.y = element_blank()) + 
  plot_risks.singvar(singvar_phthalates_taxa_run5, window = "Y1", taxa = "g4", 
                     title = "", x_title = "", y_title = "Relative abundance (%)", 
                     legend.position = "right", axis.y = element_blank()) +   
  
  plot_layout(ncol = 4, nrow = 3) 


ggsave("4_output/bkmr/Run5 (ms taxa)/plot_risks.singvar_phthalates_taxa_p_run5ms.tiff", 
       plot_risks.singvar_phthalates_taxa_p_run5, 
       device = "tiff",
       units = "mm",
       dpi = 300, 
       width = 300,
       height = 220)

ggsave("4_output/bkmr/Run5 (ms taxa)/plot_risks.singvar_phthalates_taxa_g_run5ms.tiff", 
       plot_risks.singvar_phthalates_taxa_g_run5, 
       device = "tiff",
       units = "mm",
       dpi = 300, 
       width = 300,
       height = 220)



## Run 6 ----
plot_risks.singvar_phthalates_taxa_p_run6 <- 
  plot_risks.singvar(singvar_phthalates_taxa_run6, window = "T2", taxa = "p1", 
                     title = "Phylum Firmicutes", x_title = bquote("2"^{nd}~trim.~exposure), y_title = "", 
                     legend.position = "none", axis.y = element_text(size = 12)) +
  plot_risks.singvar(singvar_phthalates_taxa_run6, window = "T2", taxa = "p2",  
                     title = "Phylum Actinobacteria", x_title = "", y_title = "", 
                     legend.position = "none", axis.y = element_blank()) + 
  plot_risks.singvar(singvar_phthalates_taxa_run6, window = "T2", taxa = "p3",  
                     title = "Phylum Bacteroidetes", x_title = "", y_title = "", 
                     legend.position = "none", axis.y = element_blank()) + 
  plot_risks.singvar(singvar_phthalates_taxa_run6, window = "T2", taxa = "p4",  
                     title = "Phylum Proteobacteria", x_title = "", y_title = "", 
                     legend.position = "right", axis.y = element_blank()) + 
  
  plot_risks.singvar(singvar_phthalates_taxa_run6, window = "T3", taxa = "p1", 
                     title = "", x_title = bquote("3"^{rd}~trim.~exposure), y_title = "", 
                     legend.position = "none", axis.y = element_text(size = 12)) +
  plot_risks.singvar(singvar_phthalates_taxa_run6, window = "T3", taxa = "p2",
                     title = "", x_title = "", y_title = "", 
                     legend.position = "none", axis.y = element_blank()) + 
  plot_risks.singvar(singvar_phthalates_taxa_run6, window = "T3", taxa = "p3", 
                     title = "", x_title = "", y_title = "", 
                     legend.position = "none", axis.y = element_blank()) + 
  plot_risks.singvar(singvar_phthalates_taxa_run6, window = "T3", taxa = "p4", 
                     title = "", x_title = "", y_title = "", 
                     legend.position = "right", axis.y = element_blank()) + 
  
  plot_risks.singvar(singvar_phthalates_taxa_run6, window = "Y1", taxa = "p1", 
                     title = "", x_title = "12-month exposure", y_title = "Relative abundance (%)", 
                     legend.position = "none", axis.y = element_text(size = 12)) +
  plot_risks.singvar(singvar_phthalates_taxa_run6, window = "Y1", taxa = "p2", 
                     title = "", x_title = "", y_title = "Relative abundance (%)", 
                     legend.position = "none", axis.y = element_blank()) + 
  plot_risks.singvar(singvar_phthalates_taxa_run6, window = "Y1", taxa = "p3", 
                     title = "", x_title = "", y_title = "Relative abundance (%)", 
                     legend.position = "none", axis.y = element_blank()) + 
  plot_risks.singvar(singvar_phthalates_taxa_run6, window = "Y1", taxa = "p4", 
                     title = "", x_title = "", y_title = "Relative abundance (%)", 
                     legend.position = "right", axis.y = element_blank()) +   
  
  plot_layout(ncol = 4, nrow = 3) 

plot_risks.singvar_phthalates_taxa_g_run6 <- 
  plot_risks.singvar(singvar_phthalates_taxa_run6, window = "T2", taxa = "g3",  
                     title = "Genus *Blautia*", x_title = bquote("2"^{nd}~trim.~exposure),  y_title = "", 
                     legend.position = "none", axis.y = element_text(size = 12)) + 
  plot_risks.singvar(singvar_phthalates_taxa_run6, window = "T2", taxa = "g1", 
                     title = "Genus *Bifidobacterium*", x_title = "", y_title = "", 
                     legend.position = "none", axis.y = element_blank()) +
  plot_risks.singvar(singvar_phthalates_taxa_run6, window = "T2", taxa = "g2",  
                     title = "Genus *Bacteroides*", x_title = "", y_title = "", 
                     legend.position = "none", axis.y = element_blank()) + 
  plot_risks.singvar(singvar_phthalates_taxa_run6, window = "T2", taxa = "g4",  
                     title = "Genera *Escherichia*<br>and *Shigella*", x_title = "", y_title = "", 
                     legend.position = "right", axis.y = element_blank()) + 
  
  plot_risks.singvar(singvar_phthalates_taxa_run6, window = "T3", taxa = "g3", 
                     title = "", x_title = bquote("3"^{rd}~trim.~exposure), y_title = "", 
                     legend.position = "none", axis.y = element_text(size = 12)) + 
  plot_risks.singvar(singvar_phthalates_taxa_run6, window = "T3", taxa = "g1", 
                     title = "", x_title = "", y_title = "", 
                     legend.position = "none", axis.y = element_blank()) +
  plot_risks.singvar(singvar_phthalates_taxa_run6, window = "T3", taxa = "g2",
                     title = "", x_title = "", y_title = "", 
                     legend.position = "none", axis.y = element_blank()) + 
  plot_risks.singvar(singvar_phthalates_taxa_run6, window = "T3", taxa = "g4", 
                     title = "", x_title = "", y_title = "", 
                     legend.position = "right", axis.y = element_blank()) + 
  
  plot_risks.singvar(singvar_phthalates_taxa_run6, window = "Y1", taxa = "g3", 
                     title = "", x_title = "12-month exposure", y_title = "Relative abundance (%)", 
                     legend.position = "none", axis.y = element_text(size = 12)) + 
  plot_risks.singvar(singvar_phthalates_taxa_run6, window = "Y1", taxa = "g1", 
                     title = "", x_title = "", y_title = "Relative abundance (%)", 
                     legend.position = "none", axis.y = element_blank()) +
  plot_risks.singvar(singvar_phthalates_taxa_run6, window = "Y1", taxa = "g2", 
                     title = "", x_title = "", y_title = "Relative abundance (%)", 
                     legend.position = "none", axis.y = element_blank()) + 
  plot_risks.singvar(singvar_phthalates_taxa_run6, window = "Y1", taxa = "g4", 
                     title = "", x_title = "", y_title = "Relative abundance (%)", 
                     legend.position = "right", axis.y = element_blank()) +   
  
  plot_layout(ncol = 4, nrow = 3) 


ggsave("4_output/bkmr/run6 (ms taxa)/plot_risks.singvar_phthalates_taxa_p_run6ms.tiff", 
       plot_risks.singvar_phthalates_taxa_p_run6, 
       device = "tiff",
       units = "mm",
       dpi = 300, 
       width = 300,
       height = 220)

ggsave("4_output/bkmr/run6 (ms taxa)/plot_risks.singvar_phthalates_taxa_g_run6ms.tiff", 
       plot_risks.singvar_phthalates_taxa_g_run6, 
       device = "tiff",
       units = "mm",
       dpi = 300, 
       width = 300,
       height = 220)
