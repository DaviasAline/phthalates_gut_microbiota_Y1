# Revisions EP - article phthalates - analyses multivariées
# Aline Davias 
# 27.08.2024

# Chargement des packages ----
library(tidyverse)
library(bkmr)
library(fields)
library(future)
library(future.apply)
library(writexl)

# Chargement des données ----
load("1_intermediate_data/2_data_selection_AD_gumme.RData")
asv_raw_not_rarefied <-              # pour les analyses betadiv
  read_labelled_csv(
    "0_source_data/gut_microbiota_ASVbased_Y1_labelled_AD_20220504_7239.csv") 

# Chargement des fonctions ----
bkmr_t2_alpha <- function(outcome) {
  set.seed(111)
  results_t2 <- kmbayes(    # modèle sans les variables catégorielles "numérisées" / toutes fenetres confondues
    y = outcome, 
    Z = mixture_alpha_t2, 
    X = covariates_alpha_t2, 
    iter = 50000,          # mettre 50 000
    verbose = FALSE,       # if TRUE, la sortie intermédiaire résumant la progression de l'ajustement du modèle est imprimée
    varsel = TRUE)         # if TRUE, we can fit the model with variable selection and estimate the posterior inclusion probability (PIP) for each of the exposures zim
  return(results_t2)
}

bkmr_t3_alpha <- function(outcome) {
  set.seed(111)
  results_t3 <- kmbayes(    # modèle avec les variables catégorielles "numérisées" / toutes fenetres confondues
    y = outcome, 
    Z = mixture_alpha_t3, 
    X = covariates_alpha_t3, 
    iter = 50000,          # mettre 50 000
    verbose = FALSE,       # if TRUE, la sortie intermédiaire résumant la progression de l'ajustement du modèle est imprimée
    varsel = TRUE)         # if TRUE, we can fit the model with variable selection and estimate the posterior inclusion probability (PIP) for each of the exposures zim
  return(results_t3)
}


bkmr_Y1_alpha <- function(outcome) {
  set.seed(111)
  results_Y1 <- kmbayes(    # modèle avec les variables catégorielles "numérisées" / fenetres séparées
    y = outcome, 
    Z = mixture_alpha_Y1, 
    X = covariates_alpha_Y1, 
    iter = 50000,          # mettre 50 000
    verbose = FALSE,       # if TRUE, la sortie intermédiaire résumant la progression de l'ajustement du modèle est imprimée
    varsel = TRUE)         # if TRUE, we can fit the model with variable selection and estimate the posterior inclusion probability (PIP) for each of the exposures zim
  return(results_Y1)
}


bkmr_t2_phyla <- function(outcome) {
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

bkmr_t3_phyla <- function(outcome) {
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

bkmr_Y1_phyla <- function(outcome) {
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


TracePlot_group_alpha <- function(model_specrich, model_shannon, titre){
  par(mfrow=c(2,3))
  TracePlot(fit = model_specrich, par = "beta", sel = TRUE) 
  TracePlot(fit = model_specrich, par = "sigsq.eps", sel = TRUE) 
  TracePlot(fit = model_specrich, par = "r", comp = 1, sel = TRUE) 
  TracePlot(fit = model_shannon, par = "beta", sel = TRUE) 
  TracePlot(fit = model_shannon, par = "sigsq.eps", sel = TRUE) 
  TracePlot(fit = model_shannon, par = "r", comp = 1, sel = TRUE) 
  mtext(titre, outer=TRUE, line = -1.5,font=2, cex=1, padj = 0)
}

TracePlot_group_phyla <- function(model_p1, model_p2, model_p3, model_p4, titre){
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

pip_results_alpha <- function(bkmr_specrich, bkmr_shannon) {
  bkmr_pip_spechrich <- 
    ExtractPIPs(bkmr_specrich) %>% 
    as.data.frame() %>% 
    rename(PIP_specrich = PIP) 
  
  bkmr_pip_shannon <- 
    ExtractPIPs(bkmr_shannon) %>% 
    as.data.frame() %>% 
    rename(PIP_shannon = PIP) 
  
  results <- left_join(bkmr_pip_spechrich, bkmr_pip_shannon, by = "variable")
  return(results)
  
}

pip_results_phyla <- function(bkmr_p1, bkmr_p2, bkmr_p3, bkmr_p4) {
  
  bkmr_pip_p1 <- ExtractPIPs(bkmr_p1) %>% as.data.frame() %>% rename(PIP_p1 = PIP) 
  bkmr_pip_p2 <- ExtractPIPs(bkmr_p2) %>% as.data.frame() %>% rename(PIP_p2 = PIP) 
  bkmr_pip_p3 <- ExtractPIPs(bkmr_p3) %>% as.data.frame() %>% rename(PIP_p3 = PIP) 
  bkmr_pip_p4 <- ExtractPIPs(bkmr_p4) %>% as.data.frame() %>% rename(PIP_p4 = PIP) 
  
  results <- list(bkmr_pip_p1, bkmr_pip_p2, bkmr_pip_p3, bkmr_pip_p4)
  results <- reduce(results, left_join, by = "variable")
  return(results)
  
}

risks_overall <- function(fit, y, Z, covariates) {
  results <- OverallRiskSummaries(
    fit = fit, 
    y = y,
    Z = Z,
    X = covariates,
    qs = seq(0.10, 0.90, by = 0.10),
    q.fixed = 0.10,
    method = "exact")
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

plot_risks.overall_sperich <- function(risks.overall, est, sd, title){
  ggplot(risks.overall,
         aes(
           quantile,
           {{est}},
           ymin = {{est}} - 1.96 * {{sd}},
           ymax = {{est}} + 1.96 * {{sd}}
         )) +
    geom_pointrange() +
    labs(y = "Specific richness") + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
    theme_bw() +
    theme(plot.title = element_text(size = 12))+ 
    ggtitle(title)
}

plot_risks.overall_shannon <-function(risks.overall, est, sd){   
  ggplot(risks.overall,
         aes(
           quantile,
           {{est}},
           ymin = {{est}} - 1.96 * {{sd}},
           ymax = {{est}} + 1.96 * {{sd}}
         )) +
    geom_pointrange() + 
    labs(y = "Shannon diversity") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
    theme_bw() }


plot_risks.singvar_sperich <- function(risks.singvar, plot_title, window) {
  
  plot_risks.singvar_sperich <- risks.singvar %>%
    filter(grepl(window, variable)) %>%
    ggplot(
      aes(
        variable,
        est_rich,
        ymin = est_rich - 1.96 * sd_rich,
        ymax = est_rich + 1.96 * sd_rich,
        col = q.fixed)
    ) +
    geom_pointrange(position = position_dodge(width = 0.75), 
                    size = 0.15) +
    geom_hline(yintercept = 0, linetype="dashed") +
    labs(y = "Specific richness", 
         x ="Exposure") + 
    coord_flip()+
    theme_bw() +
    theme(legend.position = "none", 
          plot.title = element_text(size = 12))+ 
    ggtitle(plot_title)}

plot_risks.singvar_shannon <- function(risks.singvar, window) {
  plot_risks.singvar_shannon <- risks.singvar %>%
    filter(grepl(window, variable)) %>%
    ggplot(
      aes(
        variable,
        est_sha,
        ymin = est_sha - 1.96 * sd_sha,
        ymax = est_sha + 1.96 * sd_sha,
        col = q.fixed)
    ) +
    geom_pointrange(position = position_dodge(width = 0.75), 
                    size = 0.15) +
    geom_hline(yintercept = 0, linetype="dashed") + 
    labs(y = "Shannon diversity", 
         x ="Exposure")+
    coord_flip() +
    theme_bw() +
    theme(axis.text.y = element_blank(), axis.title.y = element_blank())+
    theme(legend.position = "none")}

# Nettoyage des données ----
## variables outcomes ----
bdd_outcomes_alpha_bkmr <- 
  bdd %>% 
  select(
    ident, 
    all_of(alpha_vec)) 

bdd_outcomes_phyla_bkmr <- 
  bdd %>% 
  select(
    ident, 
    all_of(phyla_vec)) 

## covariates ----
bdd_covariates_post_bkmr <-       # on force les variables catégorielles en variables continues 
  bdd %>%
  select(
    ident,  
    ch_feces_RUN_Y1,              # covariable à 2 catégories : ok 0/1
    ch_feces_age_w_Y1_i,  
    po_delmod,                    # covariable à 2 catégories : ok 0/1
    ch_food_intro_Y1_3cat_i,      # covariable à transformer 
    ch_antibio_Y1_2cat_i,         
    mo_par_2cat,  
    mo_pets_i,                    # covariable à 2 catégories : ok 0/1
    ch_sex,                       # covariable à 2 catégories : ok 0/1
    mo_tob_gr_anyt_yn_n2_i,       # covariable à 2 catégories : ok 0/1
    Mo_ETS_anyT_yn1_opt_i,        # covariable à 2 catégories : ok 0/1
    ch_ETS_12m_opt36m,            # covariable à 2 catégories : ok 0/1
    mo_interpreg_3cat,            # covariable à transformer 
    mo_dipl_3cat_i,               # covariable à transformer
    po_w_kg,
    po_he_i,
    # ch_w_Y1_i,
    # ch_he_Y1_i,
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
    ch_antibio_Y1_2cat_i = as.numeric(fct_recode(ch_antibio_Y1_2cat_i,
                                              "0" = "No",
                                              "1" = "Yes")), 
    mo_par_2cat = as.numeric(fct_recode(mo_par_2cat,
                                                 "0" = "None",
                                                 "1" = "1 child or more")), 
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
                                                "0" = "2 years or less after graduation",
                                                "1" = "3-4 years after graduation",
                                                "0" = "≥5 years after graduation")), 
    mo_dipl_3cat_i_5y = as.numeric(fct_recode(mo_dipl_3cat_i,
                                              "0" = "2 years or less after graduation",
                                              "0" = "3-4 years after graduation",
                                              "1" = "≥5 years after graduation"))) %>%
  select(-c("mo_interpreg_3cat", "mo_dipl_3cat_i", "ch_food_intro_Y1_3cat_i"))

bdd_covariates_pre_bkmr <- bdd_covariates_post_bkmr %>% select(-po_w_kg, -po_gd, -po_he_i)

covariates_pre_bkmr <- bdd_covariates_pre_bkmr %>% select(-ident) %>% colnames()
covariates_post_bkmr <- bdd_covariates_post_bkmr %>% select(-ident) %>% colnames()


## variables d'exposition ----
phthalates_t2 <- bdd %>% 
  select(all_of(phthalates)) %>%
  select(contains("t2")) %>%
  colnames()

phthalates_t3 <- bdd %>% 
  select(all_of(phthalates)) %>%
  select(contains("t3"))%>%
  colnames()

phthalates_Y1 <- bdd %>% 
  select(all_of(phthalates)) %>%
  select(contains("Y1"))%>%
  colnames()

bdd_expo_t2_bkmr <- bdd %>%          
  select(
    ident,
    all_of(phthalates_t2))%>%           # variables d'exposition (déjà logtransformé)
  na.omit()                             # conserver que les lignes avec toutes les observations completes sans données manquantes 
bdd_expo_t2_bkmr[, phthalates_t2] <- lapply(bdd_expo_t2_bkmr[, phthalates_t2], scale)  # standardiser les expositions (diviser par sd)
colnames(bdd_expo_t2_bkmr) <- c("ident", phthalates_t2)

bdd_expo_t3_bkmr <- bdd %>%          
  select(
    ident,
    all_of(phthalates_t3))%>%           # variables d'exposition (déjà logtransformé)
  na.omit()                             # conserver que les lignes avec toutes les observations completes sans données manquantes 
bdd_expo_t3_bkmr[, phthalates_t3] <- lapply(bdd_expo_t3_bkmr[, phthalates_t3], scale)  # standardiser les expositions (diviser par sd)
colnames(bdd_expo_t3_bkmr) <- c("ident", phthalates_t3)

bdd_expo_Y1_bkmr <- bdd %>%          
  select(
    ident,
    all_of(phthalates_Y1))%>%           # variables d'exposition (déjà logtransformé)
  na.omit()                             # conserver que les lignes avec toutes les observations completes sans données manquantes 
bdd_expo_Y1_bkmr[, phthalates_Y1] <- lapply(bdd_expo_Y1_bkmr[, phthalates_Y1], scale)  # standardiser les expositions (diviser par sd)
colnames(bdd_expo_Y1_bkmr) <- c("ident", phthalates_Y1)

# Préparation des matrices ----
bdd_bkmr_alpha_t2 <- 
  list(bdd_outcomes_alpha_bkmr, bdd_covariates_pre_bkmr, bdd_expo_t2_bkmr) %>%
  reduce(left_join, by = "ident") %>%
  na.omit()

mixture_alpha_t2 <- bdd_bkmr_alpha_t2 %>% select(all_of(phthalates_t2)) %>% as.matrix()
covariates_alpha_t2 <- bdd_bkmr_alpha_t2 %>% select(all_of(covariates_pre_bkmr)) %>% as.matrix()
outcome_specrich_t2 <- bdd_bkmr_alpha_t2 %>% select(ch_feces_SpecRich_5000_ASV_Y1) %>% as.matrix()
outcome_shannon_t2 <- bdd_bkmr_alpha_t2 %>% select(ch_feces_Shannon_5000_ASV_Y1) %>% as.matrix()

bdd_bkmr_alpha_t3 <- 
  list(bdd_outcomes_alpha_bkmr, bdd_covariates_pre_bkmr, bdd_expo_t3_bkmr) %>%
  reduce(left_join, by = "ident") %>%
  na.omit()

mixture_alpha_t3 <- bdd_bkmr_alpha_t3 %>% select(all_of(phthalates_t3)) %>% as.matrix()
covariates_alpha_t3 <- bdd_bkmr_alpha_t3 %>% select(all_of(covariates_pre_bkmr)) %>% as.matrix()
outcome_specrich_t3 <- bdd_bkmr_alpha_t3 %>% select(ch_feces_SpecRich_5000_ASV_Y1) %>% as.matrix()
outcome_shannon_t3 <- bdd_bkmr_alpha_t3 %>% select(ch_feces_Shannon_5000_ASV_Y1) %>% as.matrix()

bdd_bkmr_alpha_Y1 <- 
  list(bdd_outcomes_alpha_bkmr, bdd_covariates_post_bkmr, bdd_expo_Y1_bkmr) %>%
  reduce(left_join, by = "ident") %>%
  na.omit()

mixture_alpha_Y1 <- bdd_bkmr_alpha_Y1 %>% select(all_of(phthalates_Y1)) %>% as.matrix()
covariates_alpha_Y1 <- bdd_bkmr_alpha_Y1 %>% select(all_of(covariates_post_bkmr)) %>% as.matrix()
outcome_specrich_Y1 <- bdd_bkmr_alpha_Y1 %>% select(ch_feces_SpecRich_5000_ASV_Y1) %>% as.matrix()
outcome_shannon_Y1 <- bdd_bkmr_alpha_Y1 %>% select(ch_feces_Shannon_5000_ASV_Y1) %>% as.matrix()

bdd_bkmr_taxa_t2 <- 
  list(bdd_outcomes_phyla_bkmr, bdd_covariates_pre_bkmr, bdd_expo_t2_bkmr) %>%
  reduce(left_join, by = "ident") %>%
  na.omit()

mixture_taxa_t2 <- bdd_bkmr_taxa_t2 %>% select(all_of(phthalates_t2)) %>% as.matrix()
covariates_taxa_t2 <- bdd_bkmr_taxa_t2 %>% select(all_of(covariates_pre_bkmr)) %>% as.matrix()
outcome_p1_t2 <- bdd_bkmr_taxa_t2 %>% select(ch_feces_rel_p1_Y1) %>% as.matrix()
outcome_p2_t2 <- bdd_bkmr_taxa_t2 %>% select(ch_feces_rel_p2_Y1) %>% as.matrix()
outcome_p3_t2 <- bdd_bkmr_taxa_t2 %>% select(ch_feces_rel_p3_Y1) %>% as.matrix()
outcome_p4_t2 <- bdd_bkmr_taxa_t2 %>% select(ch_feces_rel_p4_Y1) %>% as.matrix()

bdd_bkmr_taxa_t3 <- 
  list(bdd_outcomes_phyla_bkmr, bdd_covariates_pre_bkmr, bdd_expo_t3_bkmr) %>%
  reduce(left_join, by = "ident") %>%
  na.omit()

mixture_taxa_t3 <- bdd_bkmr_taxa_t3 %>% select(all_of(phthalates_t3)) %>% as.matrix()
covariates_taxa_t3 <- bdd_bkmr_taxa_t3 %>% select(all_of(covariates_pre_bkmr)) %>% as.matrix()
outcome_p1_t3 <- bdd_bkmr_taxa_t3 %>% select(ch_feces_rel_p1_Y1) %>% as.matrix()
outcome_p2_t3 <- bdd_bkmr_taxa_t3 %>% select(ch_feces_rel_p2_Y1) %>% as.matrix()
outcome_p3_t3 <- bdd_bkmr_taxa_t3 %>% select(ch_feces_rel_p3_Y1) %>% as.matrix()
outcome_p4_t3 <- bdd_bkmr_taxa_t3 %>% select(ch_feces_rel_p4_Y1) %>% as.matrix()

bdd_bkmr_taxa_Y1 <- 
  list(bdd_outcomes_phyla_bkmr, bdd_covariates_post_bkmr, bdd_expo_Y1_bkmr) %>%
  reduce(left_join, by = "ident") %>%
  na.omit()

mixture_taxa_Y1 <- bdd_bkmr_taxa_Y1 %>% select(all_of(phthalates_Y1)) %>% as.matrix()
covariates_taxa_Y1 <- bdd_bkmr_taxa_Y1 %>% select(all_of(covariates_post_bkmr)) %>% as.matrix()
outcome_p1_Y1 <- bdd_bkmr_taxa_Y1 %>% select(ch_feces_rel_p1_Y1) %>% as.matrix()
outcome_p2_Y1 <- bdd_bkmr_taxa_Y1 %>% select(ch_feces_rel_p2_Y1) %>% as.matrix()
outcome_p3_Y1 <- bdd_bkmr_taxa_Y1 %>% select(ch_feces_rel_p3_Y1) %>% as.matrix()
outcome_p4_Y1 <- bdd_bkmr_taxa_Y1 %>% select(ch_feces_rel_p4_Y1) %>% as.matrix()

# Fit BKMR ----
# Spécifier le planificateur future et le nombre de noyaux à utiliser
plan(multisession)
ncores <- availableCores() # ou spécifier le nombre de noyaux à utiliser

# Utiliser future_lapply pour exécuter les appels bkmr en parallèle pour les 3 outcomes
results_bkmr_phthalates_alpha_t2 <- future_lapply(list(outcome_specrich_t2, outcome_shannon_t2),
                                                  bkmr_t2_alpha,
                                                  future.seed = TRUE)
results_bkmr_phthalates_alpha_t3 <- future_lapply(list(outcome_specrich_t3, outcome_shannon_t3),
                                                  bkmr_t3_alpha,
                                                  future.seed = TRUE)
results_bkmr_phthalates_alpha_Y1 <- future_lapply(list(outcome_specrich_Y1, outcome_shannon_Y1),
                                                  bkmr_Y1_alpha,
                                                  future.seed = TRUE)

results_bkmr_phthalates_taxa_t2 <- future_lapply(list(outcome_p1_t2, outcome_p2_t2, outcome_p3_t2, outcome_p4_t2), 
                                                 bkmr_t2_phyla,
                                                 future.seed = TRUE)
results_bkmr_phthalates_taxa_t3 <- future_lapply(list(outcome_p1_t3, outcome_p2_t3, outcome_p3_t3, outcome_p4_t3), 
                                                 bkmr_t3_phyla, 
                                                 future.seed = TRUE)
results_bkmr_phthalates_taxa_Y1 <- future_lapply(list(outcome_p1_Y1, outcome_p2_Y1, outcome_p3_Y1, outcome_p4_Y1), 
                                                 bkmr_Y1_phyla, 
                                                 future.seed = TRUE)

results_bkmr <- list(
  alpha_t2 = results_bkmr_phthalates_alpha_t2,
  alpha_t3 = results_bkmr_phthalates_alpha_t3,
  alpha_Y1 = results_bkmr_phthalates_alpha_Y1,
  taxa_t2 = results_bkmr_phthalates_taxa_t2,
  taxa_t3 = results_bkmr_phthalates_taxa_t3,
  taxa_Y1 = results_bkmr_phthalates_taxa_Y1
)

names(results_bkmr$alpha_t2) <- c("Specific richness", "Shannon diversity")
names(results_bkmr$alpha_t3) <- c("Specific richness", "Shannon diversity") 
names(results_bkmr$alpha_Y1) <- c("Specific richness", "Shannon diversity")
names(results_bkmr$taxa_t2) <- c("Firmicutes", "Actinobacteria", "Bacteroidetes", "Proteobacteria") 
names(results_bkmr$taxa_t3) <- c("Firmicutes", "Actinobacteria", "Bacteroidetes", "Proteobacteria") 
names(results_bkmr$taxa_Y1) <- c("Firmicutes", "Actinobacteria", "Bacteroidetes", "Proteobacteria") 


# extraire les résultats
# results_bkmr_phthalates_alpha_t2 <- list[[1]]
# results_bkmr_phthalates_alpha_t3 <- list[[2]]
# results_bkmr_phthalates_alpha_Y1 <- list[[4]]
# 
# results_bkmr_phthalates_rich_t2 <- results_bkmr_phthalates_alpha_t2[[1]]
# results_bkmr_phthalates_sha_t2 <- results_bkmr_phthalates_alpha_t2[[2]]
# 
# results_bkmr_phthalates_rich_t3 <- results_bkmr_phthalates_alpha_t3[[1]]
# results_bkmr_phthalates_sha_t3 <- results_bkmr_phthalates_alpha_t3[[2]]
# 
# results_bkmr_phthalates_rich_Y1 <- results_bkmr_phthalates_alpha_Y1[[1]]
# results_bkmr_phthalates_sha_Y1 <- results_bkmr_phthalates_alpha_Y1[[2]]


# model convergence ----
TracePlot_group_alpha(results_bkmr$alpha_t2$`Specific richness`, 
                      results_bkmr$alpha_t2$`Shannon diversity`, 
                titre = "Model convergence Phthalates BKMR t2, Specific richness, Shannon (de haut en bas)")
TracePlot_group_alpha(results_bkmr$alpha_t3$`Specific richness`, 
                      results_bkmr$alpha_t3$`Shannon diversity`, 
                titre = "Model convergence Phthalates BKMR t3, Specific richness, Shannon (de haut en bas)")
TracePlot_group_alpha(results_bkmr$alpha_Y1$`Specific richness`, 
                      results_bkmr$alpha_Y1$`Shannon diversity`, 
                      titre = "Model convergence Phthalates BKMR t3, Specific richness, Shannon (de haut en bas)")
par(mfrow=c(1,1))

TracePlot_group_phyla(results_bkmr$taxa_t2$Firmicutes, 
                      results_bkmr$taxa_t2$Actinobacteria,
                      results_bkmr$taxa_t2$Bacteroidetes,
                      results_bkmr$taxa_t2$Proteobacteria,
                      titre = "Model convergence Phthalates BKMR t2, Firmicutes, Actinobacteria, Bacteroidetes, Proteobacteria (de haut en bas)")
TracePlot_group_phyla(results_bkmr$taxa_t3$Firmicutes, 
                      results_bkmr$taxa_t3$Actinobacteria,
                      results_bkmr$taxa_t3$Bacteroidetes,
                      results_bkmr$taxa_t3$Proteobacteria,
                      titre = "Model convergence Phthalates BKMR t3, Firmicutes, Actinobacteria, Bacteroidetes, Proteobacteria (de haut en bas)")
TracePlot_group_phyla(results_bkmr$taxa_Y1$Firmicutes, 
                      results_bkmr$taxa_Y1$Actinobacteria,
                      results_bkmr$taxa_Y1$Bacteroidetes,
                      results_bkmr$taxa_Y1$Proteobacteria,
                      titre = "Model convergence Phthalates BKMR Y1, Firmicutes, Actinobacteria, Bacteroidetes, Proteobacteria (de haut en bas)")
par(mfrow=c(1,1))


# PIP ----
Table_S6 <- list(
  T2 = pip_results_alpha(
    results_bkmr$alpha_t2$`Specific richness`, 
    results_bkmr$alpha_t2$`Shannon diversity`), 
  T3 = pip_results_alpha(
    results_bkmr$alpha_t3$`Specific richness`, 
    results_bkmr$alpha_t3$`Shannon diversity`),
  Y1 = pip_results_alpha(
    results_bkmr$alpha_Y1$`Specific richness`, 
    results_bkmr$alpha_Y1$`Shannon diversity`)) %>% 
  bind_rows()

Table_S9 <- list(
  T2 = pip_results_phyla(
    results_bkmr$taxa_t2$Firmicutes,
    results_bkmr$taxa_t2$Actinobacteria,
    results_bkmr$taxa_t2$Bacteroidetes,
    results_bkmr$taxa_t2$Proteobacteria), 
  T3 = pip_results_phyla(
    results_bkmr$taxa_t3$Firmicutes,
    results_bkmr$taxa_t3$Actinobacteria,
    results_bkmr$taxa_t3$Bacteroidetes,
    results_bkmr$taxa_t3$Proteobacteria),
  Y1 = pip_results_phyla(
    results_bkmr$taxa_Y1$Firmicutes,
    results_bkmr$taxa_Y1$Actinobacteria,
    results_bkmr$taxa_Y1$Bacteroidetes,
    results_bkmr$taxa_Y1$Proteobacteria)) %>%
  bind_rows()

# Overall ----
# results_bkmr_overall_phthalates_rich_t2 <- 
# results_bkmr_overall_phthalates_sha_t2 <- 
# 
# results_bkmr_overall_phthalates_rich_t3 <- 
# 
# results_bkmr_overall_phthalates_rich_Y1 <- 
# 
# list_overall <- list(results_bkmr_overall_phthalates_rich_t2 = results_bkmr_overall_phthalates_rich_t2, 
#                      results_bkmr_overall_phthalates_sha_t2 = results_bkmr_overall_phthalates_sha_t2, 
#                      results_bkmr_overall_phthalates_rich_t3 = results_bkmr_overall_phthalates_rich_t3, 
#                      results_bkmr_overall_phthalates_sha_t3 = results_bkmr_overall_phthalates_sha_t3, 
#                      results_bkmr_overall_phthalates_rich_Y1 = results_bkmr_overall_phthalates_rich_Y1, 
#                      results_bkmr_overall_phthalates_sha_Y1 =results_bkmr_overall_phthalates_sha_Y1)
# save(list_overall, 
#      file = "4_output/bkmr/Run3 (ms)/resuts_bkmr_overall_phthalates_alpha_run3ms.RData")
# load("4_output/bkmr/Run3 (ms)/resuts_bkmr_overall_phthalates_alpha_run3ms.RData")
# list2env(list_overall, envir = .GlobalEnv)

results_bkmr_overall <- 
  list(
    T2 = list(
      rich = risks_overall(results_bkmr$alpha_t2$`Specific richness`, outcome_specrich_t2, mixture_alpha_t2, covariates_alpha_t2),
      shan = risks_overall(results_bkmr$alpha_t2$`Shannon diversity`, outcome_shannon_t2, mixture_alpha_t2, covariates_alpha_t2),
      p1 = risks_overall(results_bkmr$taxa_t2$Firmicutes, outcome_p1_t2, mixture_taxa_t2, covariates_taxa_t2),
      p2 = risks_overall(results_bkmr$taxa_t2$Actinobacteria, outcome_p2_t2, mixture_taxa_t2, covariates_taxa_t2),
      p3 = risks_overall(results_bkmr$taxa_t2$Bacteroidetes, outcome_p3_t2, mixture_taxa_t2, covariates_taxa_t2),
      p4 = risks_overall(results_bkmr$taxa_t2$Proteobacteria, outcome_p4_t2, mixture_taxa_t2, covariates_taxa_t2)), 
    T3 = list(
      rich = risks_overall(results_bkmr$alpha_t3$`Specific richness`, outcome_specrich_t3, mixture_alpha_t3, covariates_alpha_t3), 
      shan = risks_overall(results_bkmr$alpha_t3$`Shannon diversity`, outcome_shannon_t3, mixture_alpha_t3, covariates_alpha_t3),
      p1 = risks_overall(results_bkmr$taxa_t3$Firmicutes, outcome_p1_t3, mixture_taxa_t3, covariates_taxa_t3),
      p2 = risks_overall(results_bkmr$taxa_t3$Actinobacteria, outcome_p2_t3, mixture_taxa_t3, covariates_taxa_t3),
      p3 = risks_overall(results_bkmr$taxa_t3$Bacteroidetes, outcome_p3_t3, mixture_taxa_t3, covariates_taxa_t3),
      p4 = risks_overall(results_bkmr$taxa_t3$Proteobacteria, outcome_p4_t3, mixture_taxa_t3, covariates_taxa_t3)),
    Y1 = list(
      rich = risks_overall(results_bkmr$alpha_Y1$`Specific richness`, outcome_specrich_Y1, mixture_alpha_Y1, covariates_alpha_Y1), 
      shan = risks_overall(results_bkmr$alpha_Y1$`Shannon diversity`, outcome_shannon_Y1, mixture_alpha_Y1, covariates_alpha_Y1),
      p1 = risks_overall(results_bkmr$taxa_Y1$Firmicutes, outcome_p1_Y1, mixture_taxa_Y1, covariates_taxa_Y1),
      p2 = risks_overall(results_bkmr$taxa_Y1$Actinobacteria, outcome_p2_Y1, mixture_taxa_Y1, covariates_taxa_Y1),
      p3 = risks_overall(results_bkmr$taxa_Y1$Bacteroidetes, outcome_p3_Y1, mixture_taxa_Y1, covariates_taxa_Y1),
      p4 = risks_overall(results_bkmr$taxa_Y1$Proteobacteria, outcome_p4_Y1, mixture_taxa_Y1, covariates_taxa_Y1)))

overall_phthalates_phyla <- bind_rows(
  list_t2 = do.call(rbind, results_bkmr_overall$T2) %>% 
    as.data.frame() %>% 
    rownames_to_column("outcome") %>%
    mutate(outcome = gsub("\\..$", "", outcome), 
           window = "T2"), 
  list_t3 = do.call(rbind, results_bkmr_overall$T3) %>% 
    as.data.frame() %>% 
    rownames_to_column("outcome") %>%
    mutate(outcome = gsub("\\..$", "", outcome), 
           window = "T3"),
  list_Y1 = do.call(rbind, results_bkmr_overall$Y1) %>% 
    as.data.frame() %>% 
    rownames_to_column("outcome") %>%
    mutate(outcome = gsub("\\..$", "", outcome), 
           window = "Y1")) %>%
  select(window, everything())

# Singvar ----
results_bkmr_singvar_phthalates_rich_t2 <- 

results_bkmr_singvar_phthalates_rich_t3 <- 

results_bkmr_singvar_phthalates_rich_Y1 <- 

list_singvar <- list(results_bkmr_singvar_phthalates_rich_t2, results_bkmr_singvar_phthalates_sha_t2, results_bkmr_singvar_phthalates_fai_t2,
                     results_bkmr_singvar_phthalates_rich_t3, results_bkmr_singvar_phthalates_sha_t3, results_bkmr_singvar_phthalates_fai_t3)
save(list_singvar,
     file = "4_output/phthalates/bkmr_phthalates/Run3 (ms)/resuts_bkmr_singvar_phthalates_alpha_run3ms.RData")

load("4_output/bkmr/Run3 (ms)/resuts_bkmr_singvar_phthalates_alpha_run3ms.RData")
names(list_singvar) <- c("results_bkmr_singvar_phthalates_rich_t2", 
                         "results_bkmr_singvar_phthalates_sha_t2", 
                         "results_bkmr_singvar_phthalates_fai_t2", 
                         "results_bkmr_singvar_phthalates_rich_t3", 
                         "results_bkmr_singvar_phthalates_sha_t3", 
                         "results_bkmr_singvar_phthalates_fai_t3",  
                         "results_bkmr_singvar_phthalates_rich_M2", 
                         "results_bkmr_singvar_phthalates_sha_M2", 
                         "results_bkmr_singvar_phthalates_fai_M2",
                         "results_bkmr_singvar_phthalates_rich_Y1", 
                         "results_bkmr_singvar_phthalates_sha_Y1", 
                         "results_bkmr_singvar_phthalates_fai_Y1")
list2env(list_singvar, envir = .GlobalEnv)


results_bkmr_singvar <- 
  list(
    T2 = list(
      rich = risks_singvar(results_bkmr$alpha_t2$`Specific richness`, outcome_specrich_t2, mixture_alpha_t2, covariates_alpha_t2),
      shan = risks_singvar(results_bkmr$alpha_t2$`Shannon diversity`, outcome_shannon_t2, mixture_alpha_t2, covariates_alpha_t2),
      p1 = risks_singvar(results_bkmr$taxa_t2$Firmicutes, outcome_p1_t2, mixture_taxa_t2, covariates_taxa_t2),
      p2 = risks_singvar(results_bkmr$taxa_t2$Actinobacteria, outcome_p2_t2, mixture_taxa_t2, covariates_taxa_t2),
      p3 = risks_singvar(results_bkmr$taxa_t2$Bacteroidetes, outcome_p3_t2, mixture_taxa_t2, covariates_taxa_t2),
      p4 = risks_singvar(results_bkmr$taxa_t2$Proteobacteria, outcome_p4_t2, mixture_taxa_t2, covariates_taxa_t2)), 
    T3 = list(
      rich = risks_singvar(results_bkmr$alpha_t3$`Specific richness`, outcome_specrich_t3, mixture_alpha_t3, covariates_alpha_t3),
      shan = risks_singvar(results_bkmr$alpha_t3$`Shannon diversity`, outcome_shannon_t3, mixture_alpha_t3, covariates_alpha_t3),
      p1 = risks_singvar(results_bkmr$taxa_t3$Firmicutes, outcome_p1_t3, mixture_taxa_t3, covariates_taxa_t3),
      p2 = risks_singvar(results_bkmr$taxa_t3$Actinobacteria, outcome_p2_t3, mixture_taxa_t3, covariates_taxa_t3),
      p3 = risks_singvar(results_bkmr$taxa_t3$Bacteroidetes, outcome_p3_t3, mixture_taxa_t3, covariates_taxa_t3),
      p4 = risks_singvar(results_bkmr$taxa_t3$Proteobacteria, outcome_p4_t3, mixture_taxa_t3, covariates_taxa_t3)),
    Y1 = list(
      rich = risks_singvar(results_bkmr$alpha_Y1$`Specific richness`, outcome_specrich_Y1, mixture_alpha_Y1, covariates_alpha_Y1), 
      shan = risks_singvar(results_bkmr$alpha_Y1$`Shannon diversity`, outcome_shannon_Y1, mixture_alpha_Y1, covariates_alpha_Y1),
      p1 = risks_singvar(results_bkmr$taxa_Y1$Firmicutes, outcome_p1_Y1, mixture_taxa_Y1, covariates_taxa_Y1),
      p2 = risks_singvar(results_bkmr$taxa_Y1$Actinobacteria, outcome_p2_Y1, mixture_taxa_Y1, covariates_taxa_Y1),
      p3 = risks_singvar(results_bkmr$taxa_Y1$Bacteroidetes, outcome_p3_Y1, mixture_taxa_Y1, covariates_taxa_Y1),
      p4 = risks_singvar(results_bkmr$taxa_Y1$Proteobacteria, outcome_p4_Y1, mixture_taxa_Y1, covariates_taxa_Y1)))

valeurs <- c("rich", "shan", "p1", "p2", "p3", "p4")
for (i in 1:6) {
  results_bkmr_singvar$T2[[i]] <- results_bkmr_singvar$T2[[i]] %>% mutate(taxa = valeurs[i])}
for (i in 1:6) {
  results_bkmr_singvar$T3[[i]] <- results_bkmr_singvar$T3[[i]] %>% mutate(taxa = valeurs[i])}
for (i in 1:6) {
  results_bkmr_singvar$Y1[[i]] <- results_bkmr_singvar$Y1[[i]] %>% mutate(taxa = valeurs[i])}
rm(valeurs)

singvar_phthalates_phyla <- bind_rows(
  list_t2 = do.call(rbind, results_bkmr_singvar$T2) %>% 
    as.data.frame() %>% 
    mutate(window = "T2"), 
  list_t3 = do.call(rbind, results_bkmr_singvar$T3) %>% 
    as.data.frame() %>% 
    mutate(window = "T3"),
  list_Y1 = do.call(rbind, results_bkmr_singvar$Y1) %>% 
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

# Sauvegarde ----
## Overall results ----
### Run 3 ----
results_bkmr_overall_phthalates_rich_t2 <- results_bkmr_overall_phthalates_rich_t2 %>% rename(est_rich_t2 = est, sd_rich_t2 = sd)
results_bkmr_overall_phthalates_sha_t2 <- results_bkmr_overall_phthalates_sha_t2 %>% rename(est_sha_t2 = est, sd_sha_t2 = sd)

results_bkmr_overall_phthalates_rich_t3 <- results_bkmr_overall_phthalates_rich_t3 %>% rename(est_rich_t3 = est, sd_rich_t3 = sd)
results_bkmr_overall_phthalates_sha_t3 <- results_bkmr_overall_phthalates_sha_t3 %>% rename(est_sha_t3 = est, sd_sha_t3 = sd)

results_bkmr_overall_phthalates_rich_Y1 <- results_bkmr_overall_phthalates_rich_Y1 %>% rename(est_rich_Y1 = est, sd_rich_Y1 = sd)
results_bkmr_overall_phthalates_sha_Y1 <- results_bkmr_overall_phthalates_sha_Y1 %>% rename(est_sha_Y1 = est, sd_sha_Y1 = sd)


results_bkmr_overall_phthalates_alpha <- 
  results_bkmr_overall_phthalates_rich_t2 %>%
  left_join(results_bkmr_overall_phthalates_sha_t2, by = "quantile") %>%
  
  left_join(results_bkmr_overall_phthalates_rich_t3, by =  "quantile") %>%
  left_join(results_bkmr_overall_phthalates_sha_t3, by =  "quantile") %>%
  
  left_join(results_bkmr_overall_phthalates_rich_Y1, by =  "quantile") %>%
  left_join(results_bkmr_overall_phthalates_sha_Y1, by =  "quantile") 


## Singvar results ----
results_bkmr_singvar_phthalates_rich_t2 <- results_bkmr_singvar_phthalates_rich_t2 %>% rename(est_rich = est, sd_rich = sd)
results_bkmr_singvar_phthalates_sha_t2 <- results_bkmr_singvar_phthalates_sha_t2 %>% rename(est_sha = est, sd_sha = sd)

results_bkmr_singvar_phthalates_rich_t3 <- results_bkmr_singvar_phthalates_rich_t3 %>% rename(est_rich = est, sd_rich = sd)
results_bkmr_singvar_phthalates_sha_t3 <- results_bkmr_singvar_phthalates_sha_t3 %>% rename(est_sha = est, sd_sha = sd)

results_bkmr_singvar_phthalates_rich_Y1 <- results_bkmr_singvar_phthalates_rich_Y1 %>% rename(est_rich = est, sd_rich = sd)
results_bkmr_singvar_phthalates_sha_Y1 <- results_bkmr_singvar_phthalates_sha_Y1 %>% rename(est_sha = est, sd_sha = sd)

results_bkmr_singvar_phthalates_alpha_t2 <- 
  results_bkmr_singvar_phthalates_rich_t2 %>%
  left_join(results_bkmr_singvar_phthalates_sha_t2, by = c("q.fixed", "variable")) 

results_bkmr_singvar_phthalates_alpha_t3 <- 
  results_bkmr_singvar_phthalates_rich_t3 %>%
  left_join(results_bkmr_singvar_phthalates_sha_t3, by = c("q.fixed", "variable"))

results_bkmr_singvar_phthalates_alpha_Y1 <- 
  results_bkmr_singvar_phthalates_rich_Y1 %>%
  left_join(results_bkmr_singvar_phthalates_sha_Y1, by = c("q.fixed", "variable"))

results_bkmr_singvar_phthalates_alpha <- 
  results_bkmr_singvar_phthalates_alpha_t2 %>%
  bind_rows(results_bkmr_singvar_phthalates_alpha_t3) %>%
  bind_rows(results_bkmr_singvar_phthalates_alpha_Y1) %>%
  mutate(
    variable = str_replace_all(variable, 
                               c("_i_cor_" = " ", 
                                 "mo_" = "", 
                                 "_ln" = "", 
                                 "ch_" = "", 
                                 "_ms" = "", 
                                 "DEHP" = "ΣDEHP", 
                                 "DiNP" = "ΣDiNP", 
                                 "DINCH" = "ΣDINCH")), 
    variable = fct_relevel(variable,
                           "ΣDINCH Y1", "ΣDINCH t3", "ΣDINCH t2", "ohMPHP Y1", "ohMPHP t3",
                           "ohMPHP t2", "MEP Y1", "MEP t3", "MEP t2", "MBzP Y1",
                           "MBzP t3", "MBzP t2", "MiBP Y1", "MiBP t3",
                           "MiBP t2", "MnBP Y1", "MnBP t3", "MnBP t2", "ΣDiNP Y1",
                           "ΣDiNP t3", "ΣDiNP t2", "ΣDEHP Y1", "ΣDEHP t3",
                           "ΣDEHP t2"))


## Assemblage et export ----
# results_bkmr_phthalates <- list(                           
#   pip_phthalates_alpha_t2 = pip_phthalates_alpha_t2,
#   pip_phthalates_alpha_t3 = pip_phthalates_alpha_t3,
#   pip_phthalates_alpha_Y1 = pip_phthalates_alpha_Y1, 
#   overall_phthalates_alpha = results_bkmr_overall_phthalates_alpha, 
#   singvar_phthalates_alpha = results_bkmr_singvar_phthalates_alpha)
# write_xlsx(results_bkmr_phthalates, path = "4_output/phthalates/bkmr_phthalates/Run3 (ms)/results_bkmr_phthalates_run3ms.xlsx")
# 
# results_bkmr_phthalates_run4 <- list(                           
#   pip_phthalates_alpha_t2 = pip_phthalates_alpha_t2_run4,
#   pip_phthalates_alpha_t3 = pip_phthalates_alpha_t3_run4,
#   pip_phthalates_alpha_Y1 = pip_phthalates_alpha_Y1_run4, 
#   overall_phthalates_alpha = results_bkmr_overall_phthalates_alpha_run4, 
#   singvar_phthalates_alpha = results_bkmr_singvar_phthalates_alpha_run4)
# write_xlsx(results_bkmr_phthalates_run4, path = "4_output/phthalates/bkmr_phthalates/Run4 (ms)/results_bkmr_phthalates_run4ms.xlsx")
# 


# Plots ----
## overall ----
dev.off()
plot_risks.overall_phthalates_alpha_run3 <- 
  plot_risks.overall_sperich(results_bkmr_overall_phthalates_alpha, est = est_rich_t2, sd = sd_rich_t2, title = bquote("2"^{nd}~trim.~exposure)) +
  plot_risks.overall_shannon(results_bkmr_overall_phthalates_alpha, est = est_sha_t2, sd = sd_sha_t2) + 
  
  plot_risks.overall_sperich(results_bkmr_overall_phthalates_alpha, est = est_rich_t3, sd = sd_rich_t3, title = bquote("3"^{rd}~trim.~exposure)) + 
  plot_risks.overall_shannon(results_bkmr_overall_phthalates_alpha, est = est_sha_t3, sd = sd_sha_t3) + 
  
  plot_risks.overall_sperich(results_bkmr_overall_phthalates_alpha, est = est_rich_Y1, sd = sd_rich_Y1, title = "12-month exposure") + 
  plot_risks.overall_shannon(results_bkmr_overall_phthalates_alpha, est = est_sha_Y1, sd = sd_sha_Y1) + 
  
  plot_layout(ncol = 2, nrow = 3) 


ggsave("4_output/bkmr/Run3 (ms)/plot_risks.overall_phthalates_alpha_run3ms.tiff", 
       plot_risks.overall_phthalates_alpha_run3, 
       device = "tiff",
       units = "mm",
       dpi = 300, 
       width = 180,
       height = 187)



## Singvar ----
plot_risks.singvar_phthalates <- 
  plot_risks.singvar_sperich(risks.singvar = results_bkmr_singvar_phthalates_alpha, window = "t2", plot_title = bquote("2"^{nd}~trim.~exposure)) + 
  plot_risks.singvar_shannon(risks.singvar = results_bkmr_singvar_phthalates_alpha, window = "t2") + 
  
  plot_risks.singvar_sperich(risks.singvar = results_bkmr_singvar_phthalates_alpha, window = "t3", plot_title = bquote("3"^{rd}~trim.~exposure)) + 
  plot_risks.singvar_shannon(risks.singvar = results_bkmr_singvar_phthalates_alpha, window = "t3") + 
  
  plot_risks.singvar_sperich(risks.singvar = results_bkmr_singvar_phthalates_alpha, window = "Y1", plot_title = "12-month exposure") + 
  plot_risks.singvar_shannon(risks.singvar = results_bkmr_singvar_phthalates_alpha, window = "Y1") + 
  plot_layout(ncol = 2, nrow = 3) 


ggsave("4_output/bkmr/Run3 (ms)/plot_risks.singvar_phthalates_run3ms.tiff", 
       plot_risks.singvar_phthalates, 
       device = "tiff",
       units = "mm",
       dpi = 300, 
       width = 180,
       height = 250)

rm(bkmr_t2_alpha, bkmr_t3_alpha, bkmr_Y1_alpha, TracePlot_group_alpha, pip_results_alpha, risks_overall, risks_singvar,
   plot_risks.overall_shannon, plot_risks.overall_sperich, plot_risks.singvar_shannon, plot_risks.singvar_sperich, 
   bdd_outcomes_alpha_bkmr, covariates_alpha_t2, covariates_alpha_t3, covariates_alpha_Y1, 
   mixture_alpha_t2, mixture_alpha_t3, mixture_alpha_Y1, 
   outcome_shannon_t2, outcome_shannon_t3, outcome_shannon_Y1, outcome_specrich_t2, outcome_specrich_t3, outcome_specrich_Y1, 
   ncores, nom)
