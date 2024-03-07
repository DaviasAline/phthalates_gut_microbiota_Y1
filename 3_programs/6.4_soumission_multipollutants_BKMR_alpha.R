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
source("3_programs/4_functions_AD_gumme.R", encoding = 'UTF-8')    

# Chargement des vecteurs ----
source("3_programs/4_vectors_AD_gumme.R", echo=TRUE)

phthalates_vec_t2 <- bdd_alpha %>% 
  select(all_of(phthalates_vec)) %>%
  select(contains("t2")) %>%
  select(!contains(c("MEOHP", "MECPP", "MEHHP", "MEHP", "MMCHP",     
                     "ohMiNP", "oxoMiNP", "cxMiNP",                 
                     "ohMINCH", "oxoMINCH"))) %>%  
  select(!contains("cat")) %>%
  colnames()

phthalates_vec_t3 <- bdd_alpha %>% 
  select(all_of(phthalates_vec)) %>%
  select(contains("t3"))%>%
  select(!contains(c("MEOHP", "MECPP", "MEHHP", "MEHP", "MMCHP",     
                     "ohMiNP", "oxoMiNP", "cxMiNP",                 
                     "ohMINCH", "oxoMINCH"))) %>%  
  select(!contains("cat")) %>%
  colnames()

# phthalates_vec_M2 <- bdd_alpha %>% 
#   select(all_of(phthalates_vec)) %>%
#   select(contains("M2"))%>%
#   select(!contains(c("MEOHP", "MECPP", "MEHHP", "MEHP", "MMCHP",     
#                      "ohMiNP", "oxoMiNP", "cxMiNP"))) %>%  
#   select(!contains("cat")) %>%
#   colnames()

phthalates_vec_Y1 <- bdd_alpha %>% 
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
    Z = mixture_alpha_t2, 
    X = covariates_alpha_t2, 
    iter = 50000,          # mettre 50 000
    verbose = FALSE,       # if TRUE, la sortie intermédiaire résumant la progression de l'ajustement du modèle est imprimée
    varsel = TRUE)         # if TRUE, we can fit the model with variable selection and estimate the posterior inclusion probability (PIP) for each of the exposures zim
  return(results_t2)
}

bkmr_t3 <- function(outcome) {
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

# bkmr_M2 <- function(outcome) {
#   set.seed(111)
#   results_M2 <- kmbayes(    # modèle sans les variables catégorielles "numérisées" / fenetres séparées
#     y = outcome, 
#     Z = mixture_alpha_M2, 
#     X = covariates_alpha_M2, 
#     iter = 50000,          # mettre 50 000
#     verbose = FALSE,       # if TRUE, la sortie intermédiaire résumant la progression de l'ajustement du modèle est imprimée
#     varsel = TRUE)         # if TRUE, we can fit the model with variable selection and estimate the posterior inclusion probability (PIP) for each of the exposures zim
#   return(results_M2)
# }

bkmr_Y1 <- function(outcome) {
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


TracePlot_group <- function(model_specrich, model_shannon, model_faith, titre){
  par(mfrow=c(3,3))
  TracePlot(fit = model_specrich, par = "beta", sel = TRUE) 
  TracePlot(fit = model_specrich, par = "sigsq.eps", sel = TRUE) 
  TracePlot(fit = model_specrich, par = "r", comp = 1, sel = TRUE) 
  TracePlot(fit = model_shannon, par = "beta", sel = TRUE) 
  TracePlot(fit = model_shannon, par = "sigsq.eps", sel = TRUE) 
  TracePlot(fit = model_shannon, par = "r", comp = 1, sel = TRUE) 
  TracePlot(fit = model_faith, par = "beta", sel = TRUE) 
  TracePlot(fit = model_faith, par = "sigsq.eps", sel = TRUE) 
  TracePlot(fit = model_faith, par = "r", comp = 1, sel = TRUE) 
  mtext(titre, outer=TRUE, line = -1.5,font=2, cex=1, padj = 0)
}

pip_results <- function(bkmr_specrich, bkmr_shannon, bkmr_faith) {
  bkmr_pip_spechrich <- 
    ExtractPIPs(bkmr_specrich) %>% 
    as.data.frame() %>% 
    rename(PIP_specrich = PIP) 
  
  bkmr_pip_shannon <- 
    ExtractPIPs(bkmr_shannon) %>% 
    as.data.frame() %>% 
    rename(PIP_shannon = PIP) 
  
  bkmr_pip_faith <- 
    ExtractPIPs(bkmr_faith) %>% 
    as.data.frame() %>% 
    rename(PIP_faith = PIP) 
  
  results <- left_join(bkmr_pip_spechrich, bkmr_pip_shannon, bkmr_pip_faith, by = "variable")
  results <- left_join(results, bkmr_pip_faith, by = "variable")
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
    method = "approx")
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
    method = "approx")
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

plot_risks.overall_faith <- function(risks.overall, est, sd){
  ggplot(risks.overall,
         aes(
           quantile,
           {{est}},
           ymin = {{est}} - 1.96 * {{sd}},
           ymax = {{est}} + 1.96 * {{sd}}
         )) +
    geom_pointrange() +
    labs(y = "Faith phylogenetic diversity")+
    geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
    theme_bw()} 



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

plot_risks.singvar_faith <- function(risks.singvar, window) {
  plot_risks.singvar_faith <- risks.singvar %>%
    filter(grepl(window, variable)) %>%
    ggplot(
      aes(
        variable,
        est_fai,
        ymin = est_fai - 1.96 * sd_fai,
        ymax = est_fai + 1.96 * sd_fai,
        col = q.fixed)
    ) +
    geom_pointrange(position = position_dodge(width = 0.75), 
                    size = 0.15) +
    geom_hline(yintercept = 0, linetype="dashed") +
    labs(y = "Faith phylogenetic diversity", 
         x ="Exposure") + 
    coord_flip()+
    theme_bw() +
    theme(axis.text.y = element_blank(), axis.title.y = element_blank())}



# Nettoyage des données ----
## variables outcomes ----
bdd_outcomes_alpha <- 
  bdd_alpha %>% 
  select(
    ident, 
    ch_feces_SpecRich_5000_ASV_Y1,   
    ch_feces_Shannon_5000_ASV_Y1,
    ch_feces_Faith_5000_ASV_Y1) %>%
  na.omit()

bdd_outcomes_taxa <- bdd_taxa %>%
  select(ident, 
         all_of(taxa_vec)) %>%
  na.omit

## variables d'ajustement ----
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


## variables d'exposition ----
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
# 
# bdd_expo_M2 <- metadata %>%          
#   select(
#     ident,
#     all_of(phthalates_vec_M2))%>%           # variables d'exposition (déjà logtransformé)
#   select(
#     ident, 
#     where(is.numeric))%>%               # conserver que les expositions codées en numérique
#   na.omit()                             # conserver que les lignes avec toutes les observations completes sans données manquantes 
# keep_expo_M2 <- colnames(bdd_expo_M2[, 2:7])
# bdd_expo_M2[, keep_expo_M2] <- lapply(bdd_expo_M2[, keep_expo_M2], scale)  # standardiser les expositions (diviser par sd)
# colnames(bdd_expo_M2) <- c("ident", keep_expo_M2)

bdd_expo_Y1 <- metadata %>%          
  select(
    ident,
    all_of(phthalates_vec_Y1))%>%           # variables d'exposition (déjà logtransformé)
  na.omit()                             # conserver que les lignes avec toutes les observations completes sans données manquantes 
keep_expo_Y1 <- colnames(bdd_expo_Y1[, 2:9])
bdd_expo_Y1[, keep_expo_Y1] <- lapply(bdd_expo_Y1[, keep_expo_Y1], scale)  # standardiser les expositions (diviser par sd)
colnames(bdd_expo_Y1) <- c("ident", keep_expo_Y1)

# Préparation des matrices ----
bdd_bkmr_alpha_t2 <- 
  list(bdd_outcomes_alpha, bdd_covariates_i, bdd_expo_t2) %>%
  reduce(left_join, by = "ident") %>%
  na.omit()

mixture_alpha_t2 <- bdd_bkmr_alpha_t2 %>% select(all_of(keep_expo_t2)) %>% as.matrix()
covariates_alpha_t2 <- bdd_bkmr_alpha_t2 %>% select(all_of(keep_covar)) %>% as.matrix()
outcome_specrich_t2 <- bdd_bkmr_alpha_t2 %>% select(ch_feces_SpecRich_5000_ASV_Y1) %>% as.matrix()
outcome_shannon_t2 <- bdd_bkmr_alpha_t2 %>% select(ch_feces_Shannon_5000_ASV_Y1) %>% as.matrix()
outcome_faith_t2 <- bdd_bkmr_alpha_t2 %>% select(ch_feces_Faith_5000_ASV_Y1) %>% as.matrix()

bdd_bkmr_alpha_t3 <- 
  list(bdd_outcomes_alpha, bdd_covariates_i, bdd_expo_t3) %>%
  reduce(left_join, by = "ident") %>%
  na.omit()

mixture_alpha_t3 <- bdd_bkmr_alpha_t3 %>% select(all_of(keep_expo_t3)) %>% as.matrix()
covariates_alpha_t3 <- bdd_bkmr_alpha_t3 %>% select(all_of(keep_covar)) %>% as.matrix()
outcome_specrich_t3 <- bdd_bkmr_alpha_t3 %>% select(ch_feces_SpecRich_5000_ASV_Y1) %>% as.matrix()
outcome_shannon_t3 <- bdd_bkmr_alpha_t3 %>% select(ch_feces_Shannon_5000_ASV_Y1) %>% as.matrix()
outcome_faith_t3 <- bdd_bkmr_alpha_t3 %>% select(ch_feces_Faith_5000_ASV_Y1) %>% as.matrix()

# bdd_bkmr_alpha_M2 <- 
#   list(bdd_outcomes_alpha, bdd_covariates_i, bdd_expo_M2) %>%
#   reduce(left_join, by = "ident") %>%
#   na.omit()
# 
# mixture_alpha_M2 <- bdd_bkmr_alpha_M2 %>% select(all_of(keep_expo_M2)) %>% as.matrix()
# covariates_alpha_M2 <- bdd_bkmr_alpha_M2 %>% select(all_of(keep_covar)) %>% as.matrix()
# outcome_specrich_M2 <- bdd_bkmr_alpha_M2 %>% select(ch_feces_SpecRich_5000_ASV_Y1) %>% as.matrix()
# outcome_shannon_M2 <- bdd_bkmr_alpha_M2 %>% select(ch_feces_Shannon_5000_ASV_Y1) %>% as.matrix()
# outcome_faith_M2 <- bdd_bkmr_alpha_M2 %>% select(ch_feces_Faith_5000_ASV_Y1) %>% as.matrix()

bdd_bkmr_alpha_Y1 <- 
  list(bdd_outcomes_alpha, bdd_covariates_i, bdd_expo_Y1) %>%
  reduce(left_join, by = "ident") %>%
  na.omit()

mixture_alpha_Y1 <- bdd_bkmr_alpha_Y1 %>% select(all_of(keep_expo_Y1)) %>% as.matrix()
covariates_alpha_Y1 <- bdd_bkmr_alpha_Y1 %>% select(all_of(keep_covar)) %>% as.matrix()
outcome_specrich_Y1 <- bdd_bkmr_alpha_Y1 %>% select(ch_feces_SpecRich_5000_ASV_Y1) %>% as.matrix()
outcome_shannon_Y1 <- bdd_bkmr_alpha_Y1 %>% select(ch_feces_Shannon_5000_ASV_Y1) %>% as.matrix()
outcome_faith_Y1 <- bdd_bkmr_alpha_Y1 %>% select(ch_feces_Faith_5000_ASV_Y1) %>% as.matrix()


# Fit BKMR ----
## Run 3 ----  
# Spécifier le planificateur future et le nombre de noyaux à utiliser
# plan(multisession)
# ncores <- availableCores() # ou spécifier le nombre de noyaux à utiliser
# 
# # Utiliser future_lapply pour exécuter les appels bkmr en parallèle pour les 3 outcomes
# results_bkmr_phthalates_alpha_t2 <- future_lapply(list(outcome_specrich_t2, outcome_shannon_t2, outcome_faith_t2), bkmr_t2, future.seed = TRUE)
# results_bkmr_phthalates_alpha_t3 <- future_lapply(list(outcome_specrich_t3, outcome_shannon_t3, outcome_faith_t3), bkmr_t3, future.seed = TRUE)
# results_bkmr_phthalates_alpha_M2 <- future_lapply(list(outcome_specrich_M2, outcome_shannon_M2, outcome_faith_M2), bkmr_M2, future.seed = TRUE)
# results_bkmr_phthalates_alpha_Y1 <- future_lapply(list(outcome_specrich_Y1, outcome_shannon_Y1, outcome_faith_Y1), bkmr_Y1, future.seed = TRUE)
# 
# list <- list(results_bkmr_phthalates_alpha_t2, results_bkmr_phthalates_alpha_t3, results_bkmr_phthalates_alpha_M2, results_bkmr_phthalates_alpha_Y1)
# save(list, 
#      file = "4_output/phthalates/bkmr_phthalates/Run3 (ms)/resuts_bkmr_phthalates_alpha_run3ms.RData")
load("4_output/bkmr/Run3 (ms)/resuts_bkmr_phthalates_alpha_run3ms.RData")
names(list) <- c("results_bkmr_phthalates_alpha_t2", 
                 "results_bkmr_phthalates_alpha_t3", 
                 "results_bkmr_phthalates_alpha_M2", 
                 "results_bkmr_phthalates_alpha_Y1")
list2env(list, envir = .GlobalEnv)
## Run 4 ----
# Spécifier le planificateur future et le nombre de noyaux à utiliser
# plan(multisession)
# ncores <- availableCores() # ou spécifier le nombre de noyaux à utiliser
# 
# # Utiliser future_lapply pour exécuter les appels bkmr en parallèle pour les 3 outcomes
# results_bkmr_phthalates_alpha_t2_run4 <- future_lapply(list(outcome_specrich_t2, outcome_shannon_t2, outcome_faith_t2), bkmr_t2, future.seed = TRUE)
# results_bkmr_phthalates_alpha_t3_run4 <- future_lapply(list(outcome_specrich_t3, outcome_shannon_t3, outcome_faith_t3), bkmr_t3, future.seed = TRUE)
# results_bkmr_phthalates_alpha_M2_run4 <- future_lapply(list(outcome_specrich_M2, outcome_shannon_M2, outcome_faith_M2), bkmr_M2, future.seed = TRUE)
# results_bkmr_phthalates_alpha_Y1_run4 <- future_lapply(list(outcome_specrich_Y1, outcome_shannon_Y1, outcome_faith_Y1), bkmr_Y1, future.seed = TRUE)
# 
# list_run4 <- list(results_bkmr_phthalates_alpha_t2_run4, 
#                   results_bkmr_phthalates_alpha_t3_run4, 
#                   results_bkmr_phthalates_alpha_M2_run4, 
#                   results_bkmr_phthalates_alpha_Y1_run4)
# save(list_run4, 
#      file = "4_output/phthalates/bkmr_phthalates/Run4 (ms)/resuts_bkmr_phthalates_alpha_run4.RData")

# extraire les résultats ----
results_bkmr_phthalates_rich_t2 <- results_bkmr_phthalates_alpha_t2[[1]]
results_bkmr_phthalates_sha_t2 <- results_bkmr_phthalates_alpha_t2[[2]]
results_bkmr_phthalates_fai_t2 <- results_bkmr_phthalates_alpha_t2[[3]]

results_bkmr_phthalates_rich_t3 <- results_bkmr_phthalates_alpha_t3[[1]]
results_bkmr_phthalates_sha_t3 <- results_bkmr_phthalates_alpha_t3[[2]]
results_bkmr_phthalates_fai_t3 <- results_bkmr_phthalates_alpha_t3[[3]]

# results_bkmr_phthalates_rich_M2 <- results_bkmr_phthalates_alpha_M2[[1]]
# results_bkmr_phthalates_sha_M2 <- results_bkmr_phthalates_alpha_M2[[2]]
# results_bkmr_phthalates_fai_M2 <- results_bkmr_phthalates_alpha_M2[[3]]

results_bkmr_phthalates_rich_Y1 <- results_bkmr_phthalates_alpha_Y1[[1]]
results_bkmr_phthalates_sha_Y1 <- results_bkmr_phthalates_alpha_Y1[[2]]
results_bkmr_phthalates_fai_Y1 <- results_bkmr_phthalates_alpha_Y1[[3]]

# load("C:/Users/Aline/OneDrive - etu.univ-grenoble-alpes.fr/Documents/5. R projects/pollutants_gut_microbiota_Y1/4_output/phthalates/bkmr_phthalates/Run4 (ms)/resuts_bkmr_phthalates_alpha_run4.RData")
# results_bkmr_phthalates_rich_t2_run4 <- list_run4[[1]][[1]]
# results_bkmr_phthalates_sha_t2_run4 <- list_run4[[1]][[2]]
# results_bkmr_phthalates_fai_t2_run4 <- list_run4[[1]][[3]]
# 
# results_bkmr_phthalates_rich_t3_run4 <- list_run4[[2]][[1]]
# results_bkmr_phthalates_sha_t3_run4 <- list_run4[[2]][[2]]
# results_bkmr_phthalates_fai_t3_run4 <- list_run4[[2]][[3]]
# 
# results_bkmr_phthalates_rich_M2_run4 <- list_run4[[3]][[1]]
# results_bkmr_phthalates_sha_M2_run4 <- list_run4[[3]][[2]]
# results_bkmr_phthalates_fai_M2_run4 <- list_run4[[3]][[3]]
# 
# results_bkmr_phthalates_rich_Y1_run4 <- list_run4[[4]][[1]]
# results_bkmr_phthalates_sha_Y1_run4 <- list_run4[[4]][[2]]
# results_bkmr_phthalates_fai_Y1_run4 <- list_run4[[4]][[3]]
# 
# results_bkmr_phthalates_rich_t2_run4 <- results_bkmr_phthalates_alpha_t2_run4[[1]]
# results_bkmr_phthalates_sha_t2_run4 <- results_bkmr_phthalates_alpha_t2_run4[[2]]
# results_bkmr_phthalates_fai_t2_run4 <- results_bkmr_phthalates_alpha_t2_run4[[3]]
# 
# results_bkmr_phthalates_rich_t3_run4 <- results_bkmr_phthalates_alpha_t3_run4[[1]]
# results_bkmr_phthalates_sha_t3_run4 <- results_bkmr_phthalates_alpha_t3_run4[[2]]
# results_bkmr_phthalates_fai_t3_run4 <- results_bkmr_phthalates_alpha_t3_run4[[3]]
# 
# results_bkmr_phthalates_rich_M2_run4 <- results_bkmr_phthalates_alpha_M2_run4[[1]]
# results_bkmr_phthalates_sha_M2_run4 <- results_bkmr_phthalates_alpha_M2_run4[[2]]
# results_bkmr_phthalates_fai_M2_run4 <- results_bkmr_phthalates_alpha_M2_run4[[3]]
# 
# results_bkmr_phthalates_rich_Y1_run4 <- results_bkmr_phthalates_alpha_Y1_run4[[1]]
# results_bkmr_phthalates_sha_Y1_run4 <- results_bkmr_phthalates_alpha_Y1_run4[[2]]
# results_bkmr_phthalates_fai_Y1_run4 <- results_bkmr_phthalates_alpha_Y1_run4[[3]]


# model convergence ----
TracePlot_group(results_bkmr_phthalates_rich_t2, results_bkmr_phthalates_sha_t2, results_bkmr_phthalates_fai_t2, 
                titre = "Model convergence Phthalates BKMR t2, Specific richness, Shannon, Faith (de haut en bas)")
TracePlot_group(results_bkmr_phthalates_rich_t3, results_bkmr_phthalates_sha_t3, results_bkmr_phthalates_fai_t3, 
                titre = "Model convergence Phthalates BKMR t3, Specific richness, Shannon, Faith (de haut en bas)")
TracePlot_group(results_bkmr_phthalates_rich_M2, results_bkmr_phthalates_sha_M2, results_bkmr_phthalates_fai_M2, 
                titre = "Model convergence Phthalates BKMR M2, Specific richness, Shannon, Faith (de haut en bas)")
TracePlot_group(results_bkmr_phthalates_rich_Y1, results_bkmr_phthalates_sha_Y1, results_bkmr_phthalates_fai_Y1, 
                titre = "Model convergence Phthalates BKMR Y1, Specific richness, Shannon, Faith (de haut en bas)")
par(mfrow=c(1,1))


# TracePlot_group(results_bkmr_phthalates_rich_t2_run4, results_bkmr_phthalates_sha_t2_run4, results_bkmr_phthalates_fai_t2_run4, 
#                 titre = "Model convergence Phthalates BKMR t2, Specific richness, Shannon, Faith (de haut en bas)")
# TracePlot_group(results_bkmr_phthalates_rich_t3_run4, results_bkmr_phthalates_sha_t3_run4, results_bkmr_phthalates_fai_t3_run4, 
#                 titre = "Model convergence Phthalates BKMR t3, Specific richness, Shannon, Faith (de haut en bas)")
# TracePlot_group(results_bkmr_phthalates_rich_M2_run4, results_bkmr_phthalates_sha_M2_run4, results_bkmr_phthalates_fai_M2_run4, 
#                 titre = "Model convergence Phthalates BKMR M2, Specific richness, Shannon, Faith (de haut en bas)")
# TracePlot_group(results_bkmr_phthalates_rich_Y1_run4, results_bkmr_phthalates_sha_Y1_run4, results_bkmr_phthalates_fai_Y1_run4, 
#                 titre = "Model convergence Phthalates BKMR Y1, Specific richness, Shannon, Faith (de haut en bas)")
# par(mfrow=c(1,1))


# PIP ----
pip_phthalates_alpha_t2 <- pip_results(results_bkmr_phthalates_rich_t2, results_bkmr_phthalates_sha_t2, results_bkmr_phthalates_fai_t2)
pip_phthalates_alpha_t3 <- pip_results(results_bkmr_phthalates_rich_t3, results_bkmr_phthalates_sha_t3, results_bkmr_phthalates_fai_t3)
pip_phthalates_alpha_M2 <- pip_results(results_bkmr_phthalates_rich_M2, results_bkmr_phthalates_sha_M2, results_bkmr_phthalates_fai_M2)
pip_phthalates_alpha_Y1 <- pip_results(results_bkmr_phthalates_rich_Y1, results_bkmr_phthalates_sha_Y1, results_bkmr_phthalates_fai_Y1)
# 
# pip_phthalates_alpha_t2_run4 <- pip_results(results_bkmr_phthalates_rich_t2_run4, results_bkmr_phthalates_sha_t2_run4, results_bkmr_phthalates_fai_t2_run4)
# pip_phthalates_alpha_t3_run4 <- pip_results(results_bkmr_phthalates_rich_t3_run4, results_bkmr_phthalates_sha_t3_run4, results_bkmr_phthalates_fai_t3_run4)
# pip_phthalates_alpha_M2_run4 <- pip_results(results_bkmr_phthalates_rich_M2_run4, results_bkmr_phthalates_sha_M2_run4, results_bkmr_phthalates_fai_M2_run4)
# pip_phthalates_alpha_Y1_run4 <- pip_results(results_bkmr_phthalates_rich_Y1_run4, results_bkmr_phthalates_sha_Y1_run4, results_bkmr_phthalates_fai_Y1_run4)


# Overall ----
# results_bkmr_overall_phthalates_rich_t2 <- risks_overall(results_bkmr_phthalates_rich_t2, outcome_specrich_t2, mixture_alpha_t2, covariates_alpha_t2)
# results_bkmr_overall_phthalates_sha_t2 <- risks_overall(results_bkmr_phthalates_sha_t2, outcome_shannon_t2, mixture_alpha_t2, covariates_alpha_t2)
# results_bkmr_overall_phthalates_fai_t2 <- risks_overall(results_bkmr_phthalates_fai_t2, outcome_faith_t2, mixture_alpha_t2, covariates_alpha_t2)
# 
# results_bkmr_overall_phthalates_rich_t3 <- risks_overall(results_bkmr_phthalates_rich_t3, outcome_specrich_t3, mixture_alpha_t3, covariates_alpha_t3)
# results_bkmr_overall_phthalates_sha_t3 <- risks_overall(results_bkmr_phthalates_sha_t3, outcome_shannon_t3, mixture_alpha_t3, covariates_alpha_t3)
# results_bkmr_overall_phthalates_fai_t3 <- risks_overall(results_bkmr_phthalates_fai_t3, outcome_faith_t3, mixture_alpha_t3, covariates_alpha_t3)
# 
# results_bkmr_overall_phthalates_rich_M2 <- risks_overall(results_bkmr_phthalates_rich_M2, outcome_specrich_M2, mixture_alpha_M2, covariates_alpha_M2)
# results_bkmr_overall_phthalates_sha_M2 <- risks_overall(results_bkmr_phthalates_sha_M2, outcome_shannon_M2, mixture_alpha_M2, covariates_alpha_M2)
# results_bkmr_overall_phthalates_fai_M2 <- risks_overall(results_bkmr_phthalates_fai_M2, outcome_faith_M2, mixture_alpha_M2, covariates_alpha_M2)
# 
# results_bkmr_overall_phthalates_rich_Y1 <- risks_overall(results_bkmr_phthalates_rich_Y1, outcome_specrich_Y1, mixture_alpha_Y1, covariates_alpha_Y1)
# results_bkmr_overall_phthalates_sha_Y1 <- risks_overall(results_bkmr_phthalates_sha_Y1, outcome_shannon_Y1, mixture_alpha_Y1, covariates_alpha_Y1)
# results_bkmr_overall_phthalates_fai_Y1 <- risks_overall(results_bkmr_phthalates_fai_Y1, outcome_faith_Y1, mixture_alpha_Y1, covariates_alpha_Y1)
# 
# list_overall <- list(results_bkmr_overall_phthalates_rich_t2, results_bkmr_overall_phthalates_sha_t2, results_bkmr_overall_phthalates_fai_t2, 
#                      results_bkmr_overall_phthalates_rich_t3, results_bkmr_overall_phthalates_sha_t3, results_bkmr_overall_phthalates_fai_t3, 
#                      results_bkmr_overall_phthalates_rich_M2, results_bkmr_overall_phthalates_sha_M2, results_bkmr_overall_phthalates_fai_M2, 
#                      results_bkmr_overall_phthalates_rich_Y1, results_bkmr_overall_phthalates_sha_Y1, results_bkmr_overall_phthalates_fai_Y1)
# save(list_overall, 
#      file = "4_output/phthalates/bkmr_phthalates/Run3 (ms)/resuts_bkmr_overall_phthalates_alpha_run3ms.RData")
load("4_output/bkmr/Run3 (ms)/resuts_bkmr_overall_phthalates_alpha_run3ms.RData")
names(list_overall) <- c("results_bkmr_overall_phthalates_rich_t2", 
                         "results_bkmr_overall_phthalates_sha_t2", 
                         "results_bkmr_overall_phthalates_fai_t2", 
                         "results_bkmr_overall_phthalates_rich_t3", 
                         "results_bkmr_overall_phthalates_sha_t3", 
                         "results_bkmr_overall_phthalates_fai_t3", 
                         "results_bkmr_overall_phthalates_rich_M2", 
                         "results_bkmr_overall_phthalates_sha_M2", 
                         "results_bkmr_overall_phthalates_fai_M2", 
                         "results_bkmr_overall_phthalates_rich_Y1", 
                         "results_bkmr_overall_phthalates_sha_Y1", 
                         "results_bkmr_overall_phthalates_fai_Y1")
list2env(list_overall, envir = .GlobalEnv)

# results_bkmr_overall_phthalates_rich_t2_run4 <- risks_overall(results_bkmr_phthalates_rich_t2_run4, outcome_specrich_t2, mixture_alpha_t2, covariates_alpha_t2)
# results_bkmr_overall_phthalates_sha_t2_run4 <- risks_overall(results_bkmr_phthalates_sha_t2_run4, outcome_shannon_t2, mixture_alpha_t2, covariates_alpha_t2)
# results_bkmr_overall_phthalates_fai_t2_run4 <- risks_overall(results_bkmr_phthalates_fai_t2_run4, outcome_faith_t2, mixture_alpha_t2, covariates_alpha_t2)
# 
# results_bkmr_overall_phthalates_rich_t3_run4 <- risks_overall(results_bkmr_phthalates_rich_t3_run4, outcome_specrich_t3, mixture_alpha_t3, covariates_alpha_t3)
# results_bkmr_overall_phthalates_sha_t3_run4 <- risks_overall(results_bkmr_phthalates_sha_t3_run4, outcome_shannon_t3, mixture_alpha_t3, covariates_alpha_t3)
# results_bkmr_overall_phthalates_fai_t3_run4 <- risks_overall(results_bkmr_phthalates_fai_t3_run4, outcome_faith_t3, mixture_alpha_t3, covariates_alpha_t3)
# 
# results_bkmr_overall_phthalates_rich_M2_run4 <- risks_overall(results_bkmr_phthalates_rich_M2_run4, outcome_specrich_M2, mixture_alpha_M2, covariates_alpha_M2)
# results_bkmr_overall_phthalates_sha_M2_run4 <- risks_overall(results_bkmr_phthalates_sha_M2_run4, outcome_shannon_M2, mixture_alpha_M2, covariates_alpha_M2)
# results_bkmr_overall_phthalates_fai_M2_run4 <- risks_overall(results_bkmr_phthalates_fai_M2_run4, outcome_faith_M2, mixture_alpha_M2, covariates_alpha_M2)
# 
# results_bkmr_overall_phthalates_rich_Y1_run4 <- risks_overall(results_bkmr_phthalates_rich_Y1_run4, outcome_specrich_Y1, mixture_alpha_Y1, covariates_alpha_Y1)
# results_bkmr_overall_phthalates_sha_Y1_run4 <- risks_overall(results_bkmr_phthalates_sha_Y1_run4, outcome_shannon_Y1, mixture_alpha_Y1, covariates_alpha_Y1)
# results_bkmr_overall_phthalates_fai_Y1_run4 <- risks_overall(results_bkmr_phthalates_fai_Y1_run4, outcome_faith_Y1, mixture_alpha_Y1, covariates_alpha_Y1)
# 
# list_overall_run4 <- list(
#   results_bkmr_overall_phthalates_rich_t2_run4, results_bkmr_overall_phthalates_sha_t2_run4, results_bkmr_overall_phthalates_fai_t2_run4, 
#   results_bkmr_overall_phthalates_rich_t3_run4, results_bkmr_overall_phthalates_sha_t3_run4, results_bkmr_overall_phthalates_fai_t3_run4,
#   results_bkmr_overall_phthalates_rich_M2_run4, results_bkmr_overall_phthalates_sha_M2_run4, results_bkmr_overall_phthalates_fai_M2_run4,
#   results_bkmr_overall_phthalates_rich_Y1_run4, results_bkmr_overall_phthalates_sha_Y1_run4, results_bkmr_overall_phthalates_fai_Y1_run4)
# save(list_overall_run4, 
#      file = "4_output/phthalates/bkmr_phthalates/Run4 (ms)/resuts_bkmr_overall_phthalates_alpha_run4.RData")


# Singvar ----
# results_bkmr_singvar_phthalates_rich_t2 <- risks_singvar(results_bkmr_phthalates_rich_t2, outcome_specrich_t2, mixture_alpha_t2, covariates_alpha_t2)
# results_bkmr_singvar_phthalates_sha_t2 <- risks_singvar(results_bkmr_phthalates_sha_t2, outcome_shannon_t2, mixture_alpha_t2, covariates_alpha_t2)
# results_bkmr_singvar_phthalates_fai_t2 <- risks_singvar(results_bkmr_phthalates_fai_t2, outcome_faith_t2, mixture_alpha_t2, covariates_alpha_t2)
# 
# results_bkmr_singvar_phthalates_rich_t3 <- risks_singvar(results_bkmr_phthalates_rich_t3, outcome_specrich_t3, mixture_alpha_t3, covariates_alpha_t3)
# results_bkmr_singvar_phthalates_sha_t3 <- risks_singvar(results_bkmr_phthalates_sha_t3, outcome_shannon_t3, mixture_alpha_t3, covariates_alpha_t3)
# results_bkmr_singvar_phthalates_fai_t3 <- risks_singvar(results_bkmr_phthalates_fai_t3, outcome_faith_t3, mixture_alpha_t3, covariates_alpha_t3)
# 
# results_bkmr_singvar_phthalates_rich_M2 <- risks_singvar(results_bkmr_phthalates_rich_M2, outcome_specrich_M2, mixture_alpha_M2, covariates_alpha_M2)
# results_bkmr_singvar_phthalates_sha_M2 <- risks_singvar(results_bkmr_phthalates_sha_M2, outcome_shannon_M2, mixture_alpha_M2, covariates_alpha_M2)
# results_bkmr_singvar_phthalates_fai_M2 <- risks_singvar(results_bkmr_phthalates_fai_M2, outcome_faith_M2, mixture_alpha_M2, covariates_alpha_M2)
# 
# results_bkmr_singvar_phthalates_rich_Y1 <- risks_singvar(results_bkmr_phthalates_rich_Y1, outcome_specrich_Y1, mixture_alpha_Y1, covariates_alpha_Y1)
# results_bkmr_singvar_phthalates_sha_Y1 <- risks_singvar(results_bkmr_phthalates_sha_Y1, outcome_shannon_Y1, mixture_alpha_Y1, covariates_alpha_Y1)
# results_bkmr_singvar_phthalates_fai_Y1 <- risks_singvar(results_bkmr_phthalates_fai_Y1, outcome_faith_Y1, mixture_alpha_Y1, covariates_alpha_Y1)
# 
# list_singvar <- list(results_bkmr_singvar_phthalates_rich_t2, results_bkmr_singvar_phthalates_sha_t2, results_bkmr_singvar_phthalates_fai_t2, 
#                      results_bkmr_singvar_phthalates_rich_t3, results_bkmr_singvar_phthalates_sha_t3, results_bkmr_singvar_phthalates_fai_t3, 
#                      results_bkmr_singvar_phthalates_rich_M2, results_bkmr_singvar_phthalates_sha_M2, results_bkmr_singvar_phthalates_fai_M2, 
#                      results_bkmr_singvar_phthalates_rich_Y1, results_bkmr_singvar_phthalates_sha_Y1, results_bkmr_singvar_phthalates_fai_Y1)
# save(list_singvar, 
#      file = "4_output/phthalates/bkmr_phthalates/Run3 (ms)/resuts_bkmr_singvar_phthalates_alpha_run3ms.RData")

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
# 
# results_bkmr_singvar_phthalates_rich_t2_run4 <- risks_singvar(results_bkmr_phthalates_rich_t2_run4, outcome_specrich_t2, mixture_alpha_t2, covariates_alpha_t2)
# results_bkmr_singvar_phthalates_sha_t2_run4 <- risks_singvar(results_bkmr_phthalates_sha_t2_run4, outcome_shannon_t2, mixture_alpha_t2, covariates_alpha_t2)
# results_bkmr_singvar_phthalates_fai_t2_run4 <- risks_singvar(results_bkmr_phthalates_fai_t2_run4, outcome_faith_t2, mixture_alpha_t2, covariates_alpha_t2)
# 
# results_bkmr_singvar_phthalates_rich_t3_run4 <- risks_singvar(results_bkmr_phthalates_rich_t3_run4, outcome_specrich_t3, mixture_alpha_t3, covariates_alpha_t3)
# results_bkmr_singvar_phthalates_sha_t3_run4 <- risks_singvar(results_bkmr_phthalates_sha_t3_run4, outcome_shannon_t3, mixture_alpha_t3, covariates_alpha_t3)
# results_bkmr_singvar_phthalates_fai_t3_run4 <- risks_singvar(results_bkmr_phthalates_fai_t3_run4, outcome_faith_t3, mixture_alpha_t3, covariates_alpha_t3)
# 
# results_bkmr_singvar_phthalates_rich_M2_run4 <- risks_singvar(results_bkmr_phthalates_rich_M2_run4, outcome_specrich_M2, mixture_alpha_M2, covariates_alpha_M2)
# results_bkmr_singvar_phthalates_sha_M2_run4 <- risks_singvar(results_bkmr_phthalates_sha_M2_run4, outcome_shannon_M2, mixture_alpha_M2, covariates_alpha_M2)
# results_bkmr_singvar_phthalates_fai_M2_run4 <- risks_singvar(results_bkmr_phthalates_fai_M2_run4, outcome_faith_M2, mixture_alpha_M2, covariates_alpha_M2)
# 
# results_bkmr_singvar_phthalates_rich_Y1_run4 <- risks_singvar(results_bkmr_phthalates_rich_Y1_run4, outcome_specrich_Y1, mixture_alpha_Y1, covariates_alpha_Y1)
# results_bkmr_singvar_phthalates_sha_Y1_run4 <- risks_singvar(results_bkmr_phthalates_sha_Y1_run4, outcome_shannon_Y1, mixture_alpha_Y1, covariates_alpha_Y1)
# results_bkmr_singvar_phthalates_fai_Y1_run4 <- risks_singvar(results_bkmr_phthalates_fai_Y1_run4, outcome_faith_Y1, mixture_alpha_Y1, covariates_alpha_Y1)
# 
# list_singvar_run4 <- list(
#   results_bkmr_singvar_phthalates_rich_t2_run4, results_bkmr_singvar_phthalates_sha_t2_run4, results_bkmr_singvar_phthalates_fai_t2_run4, 
#   results_bkmr_singvar_phthalates_rich_t3_run4, results_bkmr_singvar_phthalates_sha_t3_run4, results_bkmr_singvar_phthalates_fai_t3_run4,
#   results_bkmr_singvar_phthalates_rich_M2_run4, results_bkmr_singvar_phthalates_sha_M2_run4, results_bkmr_singvar_phthalates_fai_M2_run4,
#   results_bkmr_singvar_phthalates_rich_Y1_run4, results_bkmr_singvar_phthalates_sha_Y1_run4, results_bkmr_singvar_phthalates_fai_Y1_run4)
# save(list_singvar_run4, 
#      file = "4_output/phthalates/bkmr_phthalates/Run4 (ms)/resuts_bkmr_singvar_phthalates_alpha_run4ms.RData")



# Sauvegarde ----
## Overall results ----
### Run 3 ----
results_bkmr_overall_phthalates_rich_t2 <- results_bkmr_overall_phthalates_rich_t2 %>% rename(est_rich_t2 = est, sd_rich_t2 = sd)
results_bkmr_overall_phthalates_sha_t2 <- results_bkmr_overall_phthalates_sha_t2 %>% rename(est_sha_t2 = est, sd_sha_t2 = sd)
results_bkmr_overall_phthalates_fai_t2 <- results_bkmr_overall_phthalates_fai_t2 %>% rename(est_fai_t2 = est, sd_fai_t2 = sd)

results_bkmr_overall_phthalates_rich_t3 <- results_bkmr_overall_phthalates_rich_t3 %>% rename(est_rich_t3 = est, sd_rich_t3 = sd)
results_bkmr_overall_phthalates_sha_t3 <- results_bkmr_overall_phthalates_sha_t3 %>% rename(est_sha_t3 = est, sd_sha_t3 = sd)
results_bkmr_overall_phthalates_fai_t3 <- results_bkmr_overall_phthalates_fai_t3 %>% rename(est_fai_t3 = est, sd_fai_t3 = sd)

results_bkmr_overall_phthalates_rich_M2 <- results_bkmr_overall_phthalates_rich_M2 %>% rename(est_rich_M2 = est, sd_rich_M2 = sd)
results_bkmr_overall_phthalates_sha_M2 <- results_bkmr_overall_phthalates_sha_M2 %>% rename(est_sha_M2 = est, sd_sha_M2 = sd)
results_bkmr_overall_phthalates_fai_M2 <- results_bkmr_overall_phthalates_fai_M2 %>% rename(est_fai_M2 = est, sd_fai_M2 = sd)

results_bkmr_overall_phthalates_rich_Y1 <- results_bkmr_overall_phthalates_rich_Y1 %>% rename(est_rich_Y1 = est, sd_rich_Y1 = sd)
results_bkmr_overall_phthalates_sha_Y1 <- results_bkmr_overall_phthalates_sha_Y1 %>% rename(est_sha_Y1 = est, sd_sha_Y1 = sd)
results_bkmr_overall_phthalates_fai_Y1 <- results_bkmr_overall_phthalates_fai_Y1 %>% rename(est_fai_Y1 = est, sd_fai_Y1 = sd)


results_bkmr_overall_phthalates_alpha <- 
  results_bkmr_overall_phthalates_rich_t2 %>%
  left_join(results_bkmr_overall_phthalates_sha_t2, by = "quantile") %>%
  left_join(results_bkmr_overall_phthalates_fai_t2, by = "quantile") %>%
  
  left_join(results_bkmr_overall_phthalates_rich_t3, by =  "quantile") %>%
  left_join(results_bkmr_overall_phthalates_sha_t3, by =  "quantile") %>%
  left_join(results_bkmr_overall_phthalates_fai_t3, by = "quantile") %>%
  
  left_join(results_bkmr_overall_phthalates_rich_M2, by =  "quantile") %>%
  left_join(results_bkmr_overall_phthalates_sha_M2, by =  "quantile") %>%
  left_join(results_bkmr_overall_phthalates_fai_M2, by = "quantile") %>%
  
  left_join(results_bkmr_overall_phthalates_rich_Y1, by =  "quantile") %>%
  left_join(results_bkmr_overall_phthalates_sha_Y1, by =  "quantile") %>%
  left_join(results_bkmr_overall_phthalates_fai_Y1, by = "quantile")


### Run 4 ----
# results_bkmr_overall_phthalates_rich_t2_run4 <- results_bkmr_overall_phthalates_rich_t2_run4 %>% rename(est_rich_t2 = est, sd_rich_t2 = sd)
# results_bkmr_overall_phthalates_sha_t2_run4 <- results_bkmr_overall_phthalates_sha_t2_run4 %>% rename(est_sha_t2 = est, sd_sha_t2 = sd)
# results_bkmr_overall_phthalates_fai_t2_run4 <- results_bkmr_overall_phthalates_fai_t2_run4 %>% rename(est_fai_t2 = est, sd_fai_t2 = sd)
# 
# results_bkmr_overall_phthalates_rich_t3_run4 <- results_bkmr_overall_phthalates_rich_t3_run4 %>% rename(est_rich_t3 = est, sd_rich_t3 = sd)
# results_bkmr_overall_phthalates_sha_t3_run4 <- results_bkmr_overall_phthalates_sha_t3_run4 %>% rename(est_sha_t3 = est, sd_sha_t3 = sd)
# results_bkmr_overall_phthalates_fai_t3_run4 <- results_bkmr_overall_phthalates_fai_t3_run4 %>% rename(est_fai_t3 = est, sd_fai_t3 = sd)
# 
# results_bkmr_overall_phthalates_rich_M2_run4 <- results_bkmr_overall_phthalates_rich_M2_run4 %>% rename(est_rich_M2 = est, sd_rich_M2 = sd)
# results_bkmr_overall_phthalates_sha_M2_run4 <- results_bkmr_overall_phthalates_sha_M2_run4 %>% rename(est_sha_M2 = est, sd_sha_M2 = sd)
# results_bkmr_overall_phthalates_fai_M2_run4 <- results_bkmr_overall_phthalates_fai_M2_run4 %>% rename(est_fai_M2 = est, sd_fai_M2 = sd)
# 
# results_bkmr_overall_phthalates_rich_Y1_run4 <- results_bkmr_overall_phthalates_rich_Y1_run4 %>% rename(est_rich_Y1 = est, sd_rich_Y1 = sd)
# results_bkmr_overall_phthalates_sha_Y1_run4 <- results_bkmr_overall_phthalates_sha_Y1_run4 %>% rename(est_sha_Y1 = est, sd_sha_Y1 = sd)
# results_bkmr_overall_phthalates_fai_Y1_run4 <- results_bkmr_overall_phthalates_fai_Y1_run4 %>% rename(est_fai_Y1 = est, sd_fai_Y1 = sd)
# 
# 
# results_bkmr_overall_phthalates_alpha_run4 <- 
#   results_bkmr_overall_phthalates_rich_t2_run4 %>%
#   left_join(results_bkmr_overall_phthalates_sha_t2_run4, by = "quantile") %>%
#   left_join(results_bkmr_overall_phthalates_fai_t2_run4, by = "quantile") %>%
#   
#   left_join(results_bkmr_overall_phthalates_rich_t3_run4, by =  "quantile") %>%
#   left_join(results_bkmr_overall_phthalates_sha_t3_run4, by =  "quantile") %>%
#   left_join(results_bkmr_overall_phthalates_fai_t3_run4, by = "quantile") %>%
#   
#   left_join(results_bkmr_overall_phthalates_rich_M2_run4, by =  "quantile") %>%
#   left_join(results_bkmr_overall_phthalates_sha_M2_run4, by =  "quantile") %>%
#   left_join(results_bkmr_overall_phthalates_fai_M2_run4, by = "quantile") %>%
#   
#   left_join(results_bkmr_overall_phthalates_rich_Y1_run4, by =  "quantile") %>%
#   left_join(results_bkmr_overall_phthalates_sha_Y1_run4, by =  "quantile") %>%
#   left_join(results_bkmr_overall_phthalates_fai_Y1_run4, by = "quantile")


## Singvar results ----
### Run 3 ----
results_bkmr_singvar_phthalates_rich_t2 <- results_bkmr_singvar_phthalates_rich_t2 %>% rename(est_rich = est, sd_rich = sd)
results_bkmr_singvar_phthalates_sha_t2 <- results_bkmr_singvar_phthalates_sha_t2 %>% rename(est_sha = est, sd_sha = sd)
results_bkmr_singvar_phthalates_fai_t2 <- results_bkmr_singvar_phthalates_fai_t2 %>% rename(est_fai = est, sd_fai = sd)

results_bkmr_singvar_phthalates_rich_t3 <- results_bkmr_singvar_phthalates_rich_t3 %>% rename(est_rich = est, sd_rich = sd)
results_bkmr_singvar_phthalates_sha_t3 <- results_bkmr_singvar_phthalates_sha_t3 %>% rename(est_sha = est, sd_sha = sd)
results_bkmr_singvar_phthalates_fai_t3 <- results_bkmr_singvar_phthalates_fai_t3 %>% rename(est_fai = est, sd_fai = sd)

results_bkmr_singvar_phthalates_rich_M2 <- results_bkmr_singvar_phthalates_rich_M2 %>% rename(est_rich = est, sd_rich = sd)
results_bkmr_singvar_phthalates_sha_M2 <- results_bkmr_singvar_phthalates_sha_M2 %>% rename(est_sha = est, sd_sha = sd)
results_bkmr_singvar_phthalates_fai_M2 <- results_bkmr_singvar_phthalates_fai_M2 %>% rename(est_fai = est, sd_fai = sd)

results_bkmr_singvar_phthalates_rich_Y1 <- results_bkmr_singvar_phthalates_rich_Y1 %>% rename(est_rich = est, sd_rich = sd)
results_bkmr_singvar_phthalates_sha_Y1 <- results_bkmr_singvar_phthalates_sha_Y1 %>% rename(est_sha = est, sd_sha = sd)
results_bkmr_singvar_phthalates_fai_Y1 <- results_bkmr_singvar_phthalates_fai_Y1 %>% rename(est_fai = est, sd_fai = sd)

results_bkmr_singvar_phthalates_alpha_t2 <- 
  results_bkmr_singvar_phthalates_rich_t2 %>%
  left_join(results_bkmr_singvar_phthalates_sha_t2, by = c("q.fixed", "variable")) %>%
  left_join(results_bkmr_singvar_phthalates_fai_t2, by =c("q.fixed", "variable"))

results_bkmr_singvar_phthalates_alpha_t3 <- 
  results_bkmr_singvar_phthalates_rich_t3 %>%
  left_join(results_bkmr_singvar_phthalates_sha_t3, by = c("q.fixed", "variable")) %>%
  left_join(results_bkmr_singvar_phthalates_fai_t3, by =c("q.fixed", "variable"))

results_bkmr_singvar_phthalates_alpha_M2 <- 
  results_bkmr_singvar_phthalates_rich_M2 %>%
  left_join(results_bkmr_singvar_phthalates_sha_M2, by = c("q.fixed", "variable")) %>%
  left_join(results_bkmr_singvar_phthalates_fai_M2, by =c("q.fixed", "variable"))

results_bkmr_singvar_phthalates_alpha_Y1 <- 
  results_bkmr_singvar_phthalates_rich_Y1 %>%
  left_join(results_bkmr_singvar_phthalates_sha_Y1, by = c("q.fixed", "variable")) %>%
  left_join(results_bkmr_singvar_phthalates_fai_Y1, by =c("q.fixed", "variable"))

results_bkmr_singvar_phthalates_alpha <- 
  results_bkmr_singvar_phthalates_alpha_t2 %>%
  bind_rows(results_bkmr_singvar_phthalates_alpha_t3) %>%
  bind_rows(results_bkmr_singvar_phthalates_alpha_M2) %>%
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
                           "ohMPHP t2", "MEP Y1", "MEP M2", "MEP t3", "MEP t2", "MBzP Y1",
                           "MBzP M2", "MBzP t3", "MBzP t2", "MiBP Y1", "MiBP M2", "MiBP t3",
                           "MiBP t2", "MnBP Y1", "MnBP M2", "MnBP t3", "MnBP t2", "ΣDiNP Y1",
                           "ΣDiNP M2", "ΣDiNP t3", "ΣDiNP t2", "ΣDEHP Y1", "ΣDEHP M2", "ΣDEHP t3",
                           "ΣDEHP t2"))


### Run 4 ----
# results_bkmr_singvar_phthalates_rich_t2_run4 <- results_bkmr_singvar_phthalates_rich_t2_run4 %>% rename(est_rich = est, sd_rich = sd)
# results_bkmr_singvar_phthalates_sha_t2_run4 <- results_bkmr_singvar_phthalates_sha_t2_run4 %>% rename(est_sha = est, sd_sha = sd)
# results_bkmr_singvar_phthalates_fai_t2_run4 <- results_bkmr_singvar_phthalates_fai_t2_run4 %>% rename(est_fai = est, sd_fai = sd)
# 
# results_bkmr_singvar_phthalates_rich_t3_run4 <- results_bkmr_singvar_phthalates_rich_t3_run4 %>% rename(est_rich = est, sd_rich = sd)
# results_bkmr_singvar_phthalates_sha_t3_run4 <- results_bkmr_singvar_phthalates_sha_t3_run4 %>% rename(est_sha = est, sd_sha = sd)
# results_bkmr_singvar_phthalates_fai_t3_run4 <- results_bkmr_singvar_phthalates_fai_t3_run4 %>% rename(est_fai = est, sd_fai = sd)
# 
# results_bkmr_singvar_phthalates_rich_M2_run4 <- results_bkmr_singvar_phthalates_rich_M2_run4 %>% rename(est_rich = est, sd_rich = sd)
# results_bkmr_singvar_phthalates_sha_M2_run4 <- results_bkmr_singvar_phthalates_sha_M2_run4 %>% rename(est_sha = est, sd_sha = sd)
# results_bkmr_singvar_phthalates_fai_M2_run4 <- results_bkmr_singvar_phthalates_fai_M2_run4 %>% rename(est_fai = est, sd_fai = sd)
# 
# results_bkmr_singvar_phthalates_rich_Y1_run4 <- results_bkmr_singvar_phthalates_rich_Y1_run4 %>% rename(est_rich = est, sd_rich = sd)
# results_bkmr_singvar_phthalates_sha_Y1_run4 <- results_bkmr_singvar_phthalates_sha_Y1_run4 %>% rename(est_sha = est, sd_sha = sd)
# results_bkmr_singvar_phthalates_fai_Y1_run4 <- results_bkmr_singvar_phthalates_fai_Y1_run4 %>% rename(est_fai = est, sd_fai = sd)
# 
# results_bkmr_singvar_phthalates_alpha_t2_run4 <- 
#   results_bkmr_singvar_phthalates_rich_t2_run4 %>%
#   left_join(results_bkmr_singvar_phthalates_sha_t2_run4, by = c("q.fixed", "variable")) %>%
#   left_join(results_bkmr_singvar_phthalates_fai_t2_run4, by =c("q.fixed", "variable"))
# 
# results_bkmr_singvar_phthalates_alpha_t3_run4 <- 
#   results_bkmr_singvar_phthalates_rich_t3_run4 %>%
#   left_join(results_bkmr_singvar_phthalates_sha_t3_run4, by = c("q.fixed", "variable")) %>%
#   left_join(results_bkmr_singvar_phthalates_fai_t3_run4, by =c("q.fixed", "variable"))
# 
# results_bkmr_singvar_phthalates_alpha_M2_run4 <- 
#   results_bkmr_singvar_phthalates_rich_M2_run4 %>%
#   left_join(results_bkmr_singvar_phthalates_sha_M2_run4, by = c("q.fixed", "variable")) %>%
#   left_join(results_bkmr_singvar_phthalates_fai_M2_run4, by =c("q.fixed", "variable"))
# 
# results_bkmr_singvar_phthalates_alpha_Y1_run4 <- 
#   results_bkmr_singvar_phthalates_rich_Y1_run4 %>%
#   left_join(results_bkmr_singvar_phthalates_sha_Y1_run4, by = c("q.fixed", "variable")) %>%
#   left_join(results_bkmr_singvar_phthalates_fai_Y1_run4, by =c("q.fixed", "variable"))
# 
# results_bkmr_singvar_phthalates_alpha_run4 <- 
#   results_bkmr_singvar_phthalates_alpha_t2_run4 %>%
#   bind_rows(results_bkmr_singvar_phthalates_alpha_t3_run4) %>%
#   bind_rows(results_bkmr_singvar_phthalates_alpha_M2_run4) %>%
#   bind_rows(results_bkmr_singvar_phthalates_alpha_Y1_run4) %>%
#   mutate(
#     variable = str_replace_all(variable, 
#                                c("_i_cor_" = " ", 
#                                  "mo_" = "", 
#                                  "_ln" = "", 
#                                  "ch_" = "", 
#                                  "_ms" = "", 
#                                  "DEHP" = "ΣDEHP", 
#                                  "DiNP" = "ΣDiNP", 
#                                  "DINCH" = "ΣDINCH")), 
#     variable = fct_relevel(variable,
#                            "ΣDINCH Y1", "ΣDINCH t3", "ΣDINCH t2", "ohMPHP Y1", "ohMPHP t3",
#                            "ohMPHP t2", "MEP Y1", "MEP M2", "MEP t3", "MEP t2", "MBzP Y1",
#                            "MBzP M2", "MBzP t3", "MBzP t2", "MiBP Y1", "MiBP M2", "MiBP t3",
#                            "MiBP t2", "MnBP Y1", "MnBP M2", "MnBP t3", "MnBP t2", "ΣDiNP Y1",
#                            "ΣDiNP M2", "ΣDiNP t3", "ΣDiNP t2", "ΣDEHP Y1", "ΣDEHP M2", "ΣDEHP t3",
#                            "ΣDEHP t2"))

## Assemblage et export ----
# results_bkmr_phthalates <- list(                           
#   pip_phthalates_alpha_t2 = pip_phthalates_alpha_t2,
#   pip_phthalates_alpha_t3 = pip_phthalates_alpha_t3,
#   pip_phthalates_alpha_M2 = pip_phthalates_alpha_M2,
#   pip_phthalates_alpha_Y1 = pip_phthalates_alpha_Y1, 
#   overall_phthalates_alpha = results_bkmr_overall_phthalates_alpha, 
#   singvar_phthalates_alpha = results_bkmr_singvar_phthalates_alpha)
# write_xlsx(results_bkmr_phthalates, path = "4_output/phthalates/bkmr_phthalates/Run3 (ms)/results_bkmr_phthalates_run3ms.xlsx")
# 
# results_bkmr_phthalates_run4 <- list(                           
#   pip_phthalates_alpha_t2 = pip_phthalates_alpha_t2_run4,
#   pip_phthalates_alpha_t3 = pip_phthalates_alpha_t3_run4,
#   pip_phthalates_alpha_M2 = pip_phthalates_alpha_M2_run4,
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
  plot_risks.overall_faith(results_bkmr_overall_phthalates_alpha, est = est_fai_t2, sd = sd_fai_t2) + 
  
  plot_risks.overall_sperich(results_bkmr_overall_phthalates_alpha, est = est_rich_t3, sd = sd_rich_t3, title = bquote("3"^{rd}~trim.~exposure)) + 
  plot_risks.overall_shannon(results_bkmr_overall_phthalates_alpha, est = est_sha_t3, sd = sd_sha_t3) + 
  plot_risks.overall_faith(results_bkmr_overall_phthalates_alpha, est = est_fai_t3, sd = sd_fai_t3) + 
  
  # plot_risks.overall_sperich(results_bkmr_overall_phthalates_alpha, est = est_rich_M2, sd = sd_rich_M2, title = "2-month exposure") + 
  # plot_risks.overall_shannon(results_bkmr_overall_phthalates_alpha, est = est_sha_M2, sd = sd_sha_M2) + 
  # plot_risks.overall_faith(results_bkmr_overall_phthalates_alpha, est = est_fai_M2, sd = sd_fai_M2) +  
  
  plot_risks.overall_sperich(results_bkmr_overall_phthalates_alpha, est = est_rich_Y1, sd = sd_rich_Y1, title = "12-month exposure") + 
  plot_risks.overall_shannon(results_bkmr_overall_phthalates_alpha, est = est_sha_Y1, sd = sd_sha_Y1) + 
  plot_risks.overall_faith(results_bkmr_overall_phthalates_alpha, est = est_fai_Y1, sd = sd_fai_Y1) +  
  
  plot_layout(ncol = 3, nrow = 3) 


ggsave("4_output/bkmr/Run3 (ms)/plot_risks.overall_phthalates_alpha_run3ms.tiff", 
       plot_risks.overall_phthalates_alpha_run3, 
       device = "tiff",
       units = "mm",
       dpi = 300, 
       width = 180,
       height = 187)
# 
# plot_risks.overall_phthalates_alpha_run4 <- 
#   plot_risks.overall_sperich(results_bkmr_overall_phthalates_alpha_run4, est = est_rich_t2, sd = sd_rich_t2, title = "2nd trim. exposure") +
#   plot_risks.overall_shannon(results_bkmr_overall_phthalates_alpha_run4, est = est_sha_t2, sd = sd_sha_t2) + 
#   plot_risks.overall_faith(results_bkmr_overall_phthalates_alpha_run4, est = est_fai_t2, sd = sd_fai_t2) + 
#   
#   plot_risks.overall_sperich(results_bkmr_overall_phthalates_alpha_run4, est = est_rich_t3, sd = sd_rich_t3, title = "3rd trim. exposure") + 
#   plot_risks.overall_shannon(results_bkmr_overall_phthalates_alpha_run4, est = est_sha_t3, sd = sd_sha_t3) + 
#   plot_risks.overall_faith(results_bkmr_overall_phthalates_alpha_run4, est = est_fai_t3, sd = sd_fai_t3) + 
#   
#   plot_risks.overall_sperich(results_bkmr_overall_phthalates_alpha_run4, est = est_rich_M2, sd = sd_rich_M2, title = "2 months exposure") + 
#   plot_risks.overall_shannon(results_bkmr_overall_phthalates_alpha_run4, est = est_sha_M2, sd = sd_sha_M2) + 
#   plot_risks.overall_faith(results_bkmr_overall_phthalates_alpha_run4, est = est_fai_M2, sd = sd_fai_M2) +  
#   
#   plot_risks.overall_sperich(results_bkmr_overall_phthalates_alpha_run4, est = est_rich_Y1, sd = sd_rich_Y1, title = "12 months exposure") + 
#   plot_risks.overall_shannon(results_bkmr_overall_phthalates_alpha_run4, est = est_sha_Y1, sd = sd_sha_Y1) + 
#   plot_risks.overall_faith(results_bkmr_overall_phthalates_alpha_run4, est = est_fai_Y1, sd = sd_fai_Y1) +  
#   
#   plot_layout(ncol = 3, nrow = 4) 
# 
# 
# ggsave("4_output/phthalates/bkmr_phthalates/Run4 (ms)/plot_risks.overall_phthalates_alpha_run4.tiff", 
#        plot_risks.overall_phthalates_alpha_run4, 
#        device = "tiff",
#        units = "mm",
#        dpi = 300, 
#        width = 180,
#        height = 250)


## Singvar ----
# plot_risks.singvar_phthalates_t2 <- 
#   plot_risks.singvar_sperich(risks.singvar = results_bkmr_singvar_phthalates_alpha, window = "t2", plot_title = "2nd trim. exposure") + 
#   plot_risks.singvar_shannon(risks.singvar = results_bkmr_singvar_phthalates_alpha, window = "t2") + 
#   plot_risks.singvar_faith(risks.singvar = results_bkmr_singvar_phthalates_alpha, window = "t2") + 
#   plot_layout(ncol = 3, nrow = 1)
# 
# plot_risks.singvar_phthalates_t3 <-
#   plot_risks.singvar_sperich(risks.singvar = results_bkmr_singvar_phthalates_alpha, window = "t3", plot_title = "3rd trim. exposure") + 
#   plot_risks.singvar_shannon(risks.singvar = results_bkmr_singvar_phthalates_alpha, window = "t3") + 
#   plot_risks.singvar_faith(risks.singvar = results_bkmr_singvar_phthalates_alpha, window = "t3") + 
#   plot_layout(ncol = 3, nrow = 1)
# 
# plot_risks.singvar_phthalates_M2 <-
#   plot_risks.singvar_sperich(risks.singvar = results_bkmr_singvar_phthalates_alpha, window = "M2", plot_title = "2-month exposure") + 
#   plot_risks.singvar_shannon(risks.singvar = results_bkmr_singvar_phthalates_alpha, window = "M2") + 
#   plot_risks.singvar_faith(risks.singvar = results_bkmr_singvar_phthalates_alpha, window = "M2") +  
#   plot_layout(ncol = 3, nrow = 1)
# 
# plot_risks.singvar_phthalates_Y1 <-
#   plot_risks.singvar_sperich(risks.singvar = results_bkmr_singvar_phthalates_alpha, window = "Y1", plot_title = "12-month exposure") + 
#   plot_risks.singvar_shannon(risks.singvar = results_bkmr_singvar_phthalates_alpha, window = "Y1") + 
#   plot_risks.singvar_faith(risks.singvar = results_bkmr_singvar_phthalates_alpha, window = "Y1") +  
#   plot_layout(ncol = 3, nrow = 1) 


plot_risks.singvar_phthalates <- 
  plot_risks.singvar_sperich(risks.singvar = results_bkmr_singvar_phthalates_alpha, window = "t2", plot_title = bquote("2"^{nd}~trim.~exposure)) + 
  plot_risks.singvar_shannon(risks.singvar = results_bkmr_singvar_phthalates_alpha, window = "t2") + 
  plot_risks.singvar_faith(risks.singvar = results_bkmr_singvar_phthalates_alpha, window = "t2") + 
  
  plot_risks.singvar_sperich(risks.singvar = results_bkmr_singvar_phthalates_alpha, window = "t3", plot_title = bquote("3"^{rd}~trim.~exposure)) + 
  plot_risks.singvar_shannon(risks.singvar = results_bkmr_singvar_phthalates_alpha, window = "t3") + 
  plot_risks.singvar_faith(risks.singvar = results_bkmr_singvar_phthalates_alpha, window = "t3") + 
  
  # plot_risks.singvar_sperich(risks.singvar = results_bkmr_singvar_phthalates_alpha, window = "M2", plot_title = "2-month exposure") + 
  # plot_risks.singvar_shannon(risks.singvar = results_bkmr_singvar_phthalates_alpha, window = "M2") + 
  # plot_risks.singvar_faith(risks.singvar = results_bkmr_singvar_phthalates_alpha, window = "M2") +  
  
  plot_risks.singvar_sperich(risks.singvar = results_bkmr_singvar_phthalates_alpha, window = "Y1", plot_title = "12-month exposure") + 
  plot_risks.singvar_shannon(risks.singvar = results_bkmr_singvar_phthalates_alpha, window = "Y1") + 
  plot_risks.singvar_faith(risks.singvar = results_bkmr_singvar_phthalates_alpha, window = "Y1") +  
  plot_layout(ncol = 3, nrow = 3) 


# ggsave("4_output/phthalates/bkmr_phthalates/Run3 (ms)/plot_risks.singvar_phthalates_t2_run3ms.tiff", 
#        plot_risks.singvar_phthalates_t2, 
#        device = "tiff",
#        units = "mm",
#        dpi = 300, 
#        width = 180,
#        height = 80)
# ggsave("4_output/phthalates/bkmr_phthalates/Run3 (ms)/plot_risks.singvar_phthalates_t3_run3ms.tiff", 
#        plot_risks.singvar_phthalates_t3, 
#        device = "tiff",
#        units = "mm",
#        dpi = 300, 
#        width = 180,
#        height = 80)
# ggsave("4_output/phthalates/bkmr_phthalates/Run3 (ms)/plot_risks.singvar_phthalates_M2_run3ms.tiff", 
#        plot_risks.singvar_phthalates_M2, 
#        device = "tiff",
#        units = "mm",
#        dpi = 300, 
#        width = 180,
#        height = 80)
# ggsave("4_output/phthalates/bkmr_phthalates/Run3 (ms)/plot_risks.singvar_phthalates_Y1_run3ms.tiff", 
#        plot_risks.singvar_phthalates_Y1, 
#        device = "tiff",
#        units = "mm",
#        dpi = 300, 
#        width = 180,
#        height = 80)

ggsave("4_output/bkmr/Run3 (ms)/plot_risks.singvar_phthalates_run3ms.tiff", 
       plot_risks.singvar_phthalates, 
       device = "tiff",
       units = "mm",
       dpi = 300, 
       width = 180,
       height = 250)

# 
# plot_risks.singvar_phthalates_t2_run4 <- 
#   plot_risks.singvar_sperich(risks.singvar = results_bkmr_singvar_phthalates_alpha_run4, window = "t2", plot_title = "2nd trim. exposure") + 
#   plot_risks.singvar_shannon(risks.singvar = results_bkmr_singvar_phthalates_alpha_run4, window = "t2") + 
#   plot_risks.singvar_faith(risks.singvar = results_bkmr_singvar_phthalates_alpha_run4, window = "t2") + 
#   plot_layout(ncol = 3, nrow = 1)
# 
# plot_risks.singvar_phthalates_t3_run4 <-
#   plot_risks.singvar_sperich(risks.singvar = results_bkmr_singvar_phthalates_alpha_run4, window = "t3", plot_title = "3rd trim. exposure") + 
#   plot_risks.singvar_shannon(risks.singvar = results_bkmr_singvar_phthalates_alpha_run4, window = "t3") + 
#   plot_risks.singvar_faith(risks.singvar = results_bkmr_singvar_phthalates_alpha_run4, window = "t3") + 
#   plot_layout(ncol = 3, nrow = 1)
# 
# plot_risks.singvar_phthalates_M2_run4 <-
#   plot_risks.singvar_sperich(risks.singvar = results_bkmr_singvar_phthalates_alpha_run4, window = "M2", plot_title = "2 months exposure") + 
#   plot_risks.singvar_shannon(risks.singvar = results_bkmr_singvar_phthalates_alpha_run4, window = "M2") + 
#   plot_risks.singvar_faith(risks.singvar = results_bkmr_singvar_phthalates_alpha_run4, window = "M2") +  
#   plot_layout(ncol = 3, nrow = 1)
# 
# plot_risks.singvar_phthalates_Y1_run4 <-
#   plot_risks.singvar_sperich(risks.singvar = results_bkmr_singvar_phthalates_alpha_run4, window = "Y1", plot_title = "12 months exposure") + 
#   plot_risks.singvar_shannon(risks.singvar = results_bkmr_singvar_phthalates_alpha_run4, window = "Y1") + 
#   plot_risks.singvar_faith(risks.singvar = results_bkmr_singvar_phthalates_alpha_run4, window = "Y1") +  
#   plot_layout(ncol = 3, nrow = 1) 
# 
# 
# plot_risks.singvar_phthalates_run4 <- 
#   plot_risks.singvar_sperich(risks.singvar = results_bkmr_singvar_phthalates_alpha_run4, window = "t2", plot_title = "2nd trim. exposure") + 
#   plot_risks.singvar_shannon(risks.singvar = results_bkmr_singvar_phthalates_alpha_run4, window = "t2") + 
#   plot_risks.singvar_faith(risks.singvar = results_bkmr_singvar_phthalates_alpha_run4, window = "t2") + 
#   
#   plot_risks.singvar_sperich(risks.singvar = results_bkmr_singvar_phthalates_alpha_run4, window = "t3", plot_title = "3rd trim. exposure") + 
#   plot_risks.singvar_shannon(risks.singvar = results_bkmr_singvar_phthalates_alpha_run4, window = "t3") + 
#   plot_risks.singvar_faith(risks.singvar = results_bkmr_singvar_phthalates_alpha_run4, window = "t3") + 
#   
#   plot_risks.singvar_sperich(risks.singvar = results_bkmr_singvar_phthalates_alpha_run4, window = "M2", plot_title = "2 months exposure") + 
#   plot_risks.singvar_shannon(risks.singvar = results_bkmr_singvar_phthalates_alpha_run4, window = "M2") + 
#   plot_risks.singvar_faith(risks.singvar = results_bkmr_singvar_phthalates_alpha_run4, window = "M2") +  
#   
#   plot_risks.singvar_sperich(risks.singvar = results_bkmr_singvar_phthalates_alpha_run4, window = "Y1", plot_title = "12 months exposure") + 
#   plot_risks.singvar_shannon(risks.singvar = results_bkmr_singvar_phthalates_alpha_run4, window = "Y1") + 
#   plot_risks.singvar_faith(risks.singvar = results_bkmr_singvar_phthalates_alpha_run4, window = "Y1") +  
#   plot_layout(ncol = 3, nrow = 4) 
# 
# 
# ggsave("4_output/phthalates/bkmr_phthalates/Run4 (ms)/plot_risks.singvar_phthalates_t2_run4ms.tiff", 
#        plot_risks.singvar_phthalates_t2_run4, 
#        device = "tiff",
#        units = "mm",
#        dpi = 300, 
#        width = 180,
#        height = 80)
# ggsave("4_output/phthalates/bkmr_phthalates/Run4 (ms)/plot_risks.singvar_phthalates_t3_run4ms.tiff", 
#        plot_risks.singvar_phthalates_t3_run4, 
#        device = "tiff",
#        units = "mm",
#        dpi = 300, 
#        width = 180,
#        height = 80)
# ggsave("4_output/phthalates/bkmr_phthalates/Run4 (ms)/plot_risks.singvar_phthalates_M2_run4ms.tiff", 
#        plot_risks.singvar_phthalates_M2_run4, 
#        device = "tiff",
#        units = "mm",
#        dpi = 300, 
#        width = 180,
#        height = 80)
# ggsave("4_output/phthalates/bkmr_phthalates/Run4 (ms)/plot_risks.singvar_phthalates_Y1_run4ms.tiff", 
#        plot_risks.singvar_phthalates_Y1_run4, 
#        device = "tiff",
#        units = "mm",
#        dpi = 300, 
#        width = 180,
#        height = 80)
# 
# ggsave("4_output/phthalates/bkmr_phthalates/Run4 (ms)/plot_risks.singvar_phthalates_run4ms.tiff", 
#        plot_risks.singvar_phthalates_run4, 
#        device = "tiff",
#        units = "mm",
#        dpi = 300, 
#        width = 180,
#        height = 250)
# 

