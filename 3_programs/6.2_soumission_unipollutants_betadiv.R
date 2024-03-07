## Aline Davias
## 21/02/2024
## Analyses de la beta diversité 

# Packages, functions and data loading ----
library(tidyverse)
library(vegan)
library(phyloseq)
library(broom)
library(forestplot)
library(bkmr)
library(fields)
library(psych)
library(expss)
library(SRS)
library(see)
library(writexl)
library(forcats)
library(egg)
library(corrplot)
load("2_final_data/metadata.RData")
load("2_final_data/bdd_alpha.RData")
load("2_final_data/bdd_taxa.RData")
source("3_programs/4_functions_AD_gumme.R", encoding = 'UTF-8')
rm(model_covar, model_multi, model_summary, model_univ_multi, table_cor, table_cor_sg, test_sensi_sg)
source("3_programs/4_vectors_AD_gumme.R")
rm(pollutant_vec_t2, pollutant_vec_t3, pollutant_vec_M2, pollutant_vec_Y1, 
   covar_vec, covar_vec_cat, covar_vec_cat_i,
   phthalates_vec_cat, phthalates_vec_ln, phthalates_vec, phthalates_vec_ter, 
   taxa_vec, 
   list = ls()[grep("sg", ls())])
rm(list = ls()[grep("phenols", ls())])
rm(list = ls()[grep("pfas", ls())])
rm(list = ls()[grep("num", ls())])
pollutants <- c("mo_DEHP_ms_i_cor_t2_ter", "mo_DEHP_ms_i_cor_t3_ter", "ch_DEHP_ms_i_cor_Y1_ter",
                "mo_MnBP_i_cor_t2_ter", "mo_MnBP_i_cor_t3_ter", "ch_MnBP_i_cor_Y1_ter", 
                "mo_DiNP_ms_i_cor_t2_ter", "mo_DiNP_ms_i_cor_t3_ter", "ch_DiNP_ms_i_cor_Y1_ter",
                "mo_MiBP_i_cor_t2_ter", "mo_MiBP_i_cor_t3_ter", "ch_MiBP_i_cor_Y1_ter", 
                "mo_MBzP_i_cor_t2_ter", "mo_MBzP_i_cor_t3_ter", "ch_MBzP_i_cor_Y1_ter",
                "mo_MEP_i_cor_t2_ter", "mo_MEP_i_cor_t3_ter", "ch_MEP_i_cor_Y1_ter",
                "mo_ohMPHP_i_cor_t2_ter", "mo_ohMPHP_i_cor_t3_ter", "ch_ohMPHP_i_cor_Y1_ter", 
                "mo_DINCH_ms_i_cor_t2_ter", "mo_DINCH_ms_i_cor_t3_ter", "ch_DINCH_ms_i_cor_Y1_ter")

asv_raw_not_rarefied <- read_labelled_csv("0_source_data/gut_microbiota_ASVbased_Y1_labelled_AD_20220504_7239.csv") 

bdd_final_Y1$ch_DINCH_ms_i_cor_Y1_ter <- bdd_final_Y1$ch_DINCH_ms_i_cor_Y1_ter %>%
  fct_recode(
    "3rd tertile" = "3rd etrtile"
  )

# Beta diversity ----
## Data cleaning ----
### Rarefaction ----
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

### Metric calculation (Bray Curtis) ----
### Calcul de la metric Bray Curtis
set.seed(1996)
metric_bray_curtis <- vegan::vegdist(ASV_rarefied_5000_Y1, method = "bray")
### warning message: Plus d’une classe "dist" est trouvée en cache : Utilisation de la première, depuis l’espace de noms 'BiocGenerics'. Aussi défini par ‘spam’
### ok, phyloseq package is supposed to use BiocGeneric 
nmds <- metaMDS(metric_bray_curtis)
scores(nmds) %>%
  as_tibble(rownames = "ident") %>%
  ggplot(aes(x=NMDS1, y=NMDS2)) +
  geom_point()
# nmds <- metaMDS(ASV_rarefied_5000_Y1, autotransform = FALSE) # ?
# scores(nmds) %>%
#   as_tibble(rownames = "ident") %>%
#   ggplot(aes(x=NMDS1, y=NMDS2)) +
#   geom_point()

metric_bray_curtis <- 
  metric_bray_curtis %>% 
  as.matrix %>% 
  as.data.frame() %>%
  rownames_to_column(var = "ident")

### Metadata ----
metadata <- metadata %>%
  select(ident, 
         all_of(covar_vec_i), 
         all_of(pollutants)) %>%
  mutate(ident = as.character(ident))

#### merge betadiversity data and metadata
bdd_final <- inner_join(metric_bray_curtis, metadata, by = "ident")
bdd_final[pollutants] <- lapply(bdd_final[pollutants], as.factor)
covariables <- bdd_final %>% select(all_of(covar_vec_i)) %>% select(-mo_interpreg_3cat) %>% colnames()

#### filter t2 (because of the NA on the pollutants)
bdd_final %>% select(ident, contains("t2")) %>% filter_all(any_vars(is.na(.)))
bdd_final_t2 <- bdd_final %>% 
  select(ident, contains("t2"), everything()) %>%
  filter(ident != 15804) %>%
  select(-"15804")
all_dist_t2 <- bdd_final_t2 %>%
  select(all_of(.[["ident"]])) %>%
  as.dist()

#### filter t3 (because of the NA on the pollutants)
bdd_final %>% select(ident, contains("t3")) %>% filter_all(any_vars(is.na(.)))
bdd_final_t3 <- bdd_final %>% 
  select(ident, contains("t3"), everything()) %>%
  filter(!ident %in% c(17827, 15929, 26891, 25668, 23330, 28199)) %>%
  select(-"15929", -"17827", -"23330", -"25668", -"26891", -"28199")
all_dist_t3 <- bdd_final_t3 %>%
  select(all_of(.[["ident"]])) %>%
  as.dist()

#### filter Y1 (because of the NA on the pollutants)
bdd_final %>% select(ident, contains("Y1")) %>% filter_all(any_vars(is.na(.)))
bdd_final_Y1 <- bdd_final %>% 
  select(ident, contains("Y1"), everything()) %>%
  filter(!ident %in% c(23994, 25166, 26766, 14668, 26923)) %>%
  select(-"23994", -"25166", -"26766", -"14668", -"26923")
all_dist_Y1 <- bdd_final_Y1 %>%
  select(all_of(.[["ident"]])) %>%
  as.dist()



## Statistical analysis Adonis2 ----
### Trim 2. ----
explanatory_vars_t2 <- 
  bdd_final_t2 %>% 
  select(all_of(pollutants)) %>%
  select(contains("t2")) %>%
  colnames()

results_betadiv_univar_bray_curtis_t2 <- 
  lapply(explanatory_vars_t2, function(x) {
    formula <- reformulate(x, response = "all_dist_t2")
    adonis2(formula, data = bdd_final_t2, permutations = 999)
  })
results_betadiv_univar_bray_curtis_t2 <- 
  do.call(rbind, results_betadiv_univar_bray_curtis_t2) %>%
  rownames_to_column(var = "Explanatory variables") %>%
    mutate(
      Pollutants = c(rep("mo_DEHP_ms_i_cor_t2_ter", times = 3), 
                     rep("mo_MnBP_i_cor_t2_ter", times = 3), 
                     rep("mo_DiNP_ms_i_cor_t2_ter", times = 3),
                     rep("mo_MiBP_i_cor_t2_ter", times = 3),
                     rep("mo_MBzP_i_cor_t2_ter", times = 3),
                     rep("mo_MEP_i_cor_t2_ter", times = 3),
                     rep("mo_ohMPHP_i_cor_t2_ter", times = 3),
                     rep("mo_DINCH_ms_i_cor_t2_ter", times = 3))) %>%
  select(Pollutants, everything())

results_betadiv_multivar_bray_curtis_t2 <-
  lapply(explanatory_vars_t2, function(x) {
    formula <- reformulate(c(x, covariables), response = "all_dist_t2")
    adonis2(formula, data = bdd_final_t2, permutations = 999)
  })
results_betadiv_multivar_bray_curtis_t2 <- 
  do.call(rbind, results_betadiv_multivar_bray_curtis_t2) %>%
  rownames_to_column(var = "Explanatory variables") %>% 
  mutate(
    Pollutants = c(rep("mo_DEHP_ms_i_cor_t2_ter", times = 23), 
                   rep("mo_MnBP_i_cor_t2_ter", times = 23), 
                   rep("mo_DiNP_ms_i_cor_t2_ter", times = 23),
                   rep("mo_MiBP_i_cor_t2_ter", times = 23),
                   rep("mo_MBzP_i_cor_t2_ter", times = 23),
                   rep("mo_MEP_i_cor_t2_ter", times = 23),
                   rep("mo_ohMPHP_i_cor_t2_ter", times = 23),
                   rep("mo_DINCH_ms_i_cor_t2_ter", times = 23))) %>%
  select(Pollutants, everything())

results_betadiv_univar_bray_curtis_t2 %>%                         # Visualisation des résultats significatifs univarié
  filter(`Pr(>F)` < 0.05) %>%
  filter(`Explanatory variables` %in% c(explanatory_vars_t2)) %>%
  View()

results_betadiv_multivar_bray_curtis_t2 %>%                      # Visualisation des résultats significatifs multivarié
  filter(`Pr(>F)` < 0.05) %>%
  filter(`Explanatory variables` %in% c(explanatory_vars_t2)) %>%
  View()

### Trim 3. ----
explanatory_vars_t3 <- 
  bdd_final_t3 %>% 
  select(all_of(pollutants)) %>%
  select(contains("t3")) %>%
  colnames()

results_betadiv_univar_bray_curtis_t3 <- 
  lapply(explanatory_vars_t3, function(x) {
    formula <- reformulate(x, response = "all_dist_t3")
    adonis2(formula, data = bdd_final_t3, permutations = 999)
  })
results_betadiv_univar_bray_curtis_t3 <- 
  do.call(rbind, results_betadiv_univar_bray_curtis_t3) %>%
  rownames_to_column(var = "Explanatory variables") %>%
  mutate(
    Pollutants = c(rep("mo_DEHP_ms_i_cor_t3_ter", times = 3), 
                   rep("mo_MnBP_i_cor_t3_ter", times = 3), 
                   rep("mo_DiNP_ms_i_cor_t3_ter", times = 3),
                   rep("mo_MiBP_i_cor_t3_ter", times = 3),
                   rep("mo_MBzP_i_cor_t3_ter", times = 3),
                   rep("mo_MEP_i_cor_t3_ter", times = 3),
                   rep("mo_ohMPHP_i_cor_t3_ter", times = 3),
                   rep("mo_DINCH_ms_i_cor_t3_ter", times = 3))) %>%
  select(Pollutants, everything())

results_betadiv_multivar_bray_curtis_t3 <-
  lapply(explanatory_vars_t3, function(x) {
    formula <- reformulate(c(x, covariables), response = "all_dist_t3")
    adonis2(formula, data = bdd_final_t3, permutations = 999)
  })
results_betadiv_multivar_bray_curtis_t3 <-
  do.call(rbind, results_betadiv_multivar_bray_curtis_t3) %>%
  rownames_to_column(var = "Explanatory variables") %>%
  mutate(
    Pollutants = c(rep("mo_DEHP_ms_i_cor_t3_ter", times = 23), 
                   rep("mo_MnBP_i_cor_t3_ter", times = 23), 
                   rep("mo_DiNP_ms_i_cor_t3_ter", times = 23),
                   rep("mo_MiBP_i_cor_t3_ter", times = 23),
                   rep("mo_MBzP_i_cor_t3_ter", times = 23),
                   rep("mo_MEP_i_cor_t3_ter", times = 23),
                   rep("mo_ohMPHP_i_cor_t3_ter", times = 23),
                   rep("mo_DINCH_ms_i_cor_t3_ter", times = 23))) %>%
  select(Pollutants, everything())

results_betadiv_univar_bray_curtis_t3 %>%                         # Visualisation des résultats significatifs univarié
  filter(`Pr(>F)` < 0.05) %>%
  filter(`Explanatory variables` %in% c(explanatory_vars_t3)) %>%
  View()

results_betadiv_multivar_bray_curtis_t3 %>%                     # Visualisation des résultats significatifs multivarié
  filter(`Pr(>F)` < 0.05) %>%
  filter(`Explanatory variables` %in% c(explanatory_vars_t3)) %>%
  View()

### Y1 ----
explanatory_vars_Y1 <- 
  bdd_final_Y1 %>% 
  select(all_of(pollutants))  %>%
  select(contains("Y1")) %>%
  colnames()

results_betadiv_univar_bray_curtis_Y1 <- 
  lapply(explanatory_vars_Y1, function(x) {
    formula <- reformulate(x, response = "all_dist_Y1")
    adonis2(formula, data = bdd_final_Y1, permutations = 999)
  })
results_betadiv_univar_bray_curtis_Y1 <- 
  do.call(rbind, results_betadiv_univar_bray_curtis_Y1) %>%
  rownames_to_column(var = "Explanatory variables") %>%
  mutate(
    Pollutants = c(rep("ch_DEHP_ms_i_cor_Y1_ter", times = 3), 
                   rep("ch_MnBP_i_cor_Y1_ter", times = 3), 
                   rep("ch_DiNP_ms_i_cor_Y1_ter", times = 3),
                   rep("ch_MiBP_i_cor_Y1_ter", times = 3),
                   rep("ch_MBzP_i_cor_Y1_ter", times = 3),
                   rep("ch_MEP_i_cor_Y1_ter", times = 3),
                   rep("ch_ohMPHP_i_cor_Y1_ter", times = 3),
                   rep("ch_DINCH_ms_i_cor_Y1_ter", times = 3))) %>%
  select(Pollutants, everything())

results_betadiv_multivar_bray_curtis_Y1 <-
  lapply(explanatory_vars_Y1, function(x) {
    formula <- reformulate(c(x, covariables), response = "all_dist_Y1")
    adonis2(formula, data = bdd_final_Y1, permutations = 999)
  })
results_betadiv_multivar_bray_curtis_Y1 <-
  do.call(rbind, results_betadiv_multivar_bray_curtis_Y1) %>%
  rownames_to_column(var = "Explanatory variables") %>%
  mutate(
    Pollutants = c(rep("ch_DEHP_ms_i_cor_Y1_ter", times = 23), 
                   rep("ch_MnBP_i_cor_Y1_ter", times = 23), 
                   rep("ch_DiNP_ms_i_cor_Y1_ter", times = 23),
                   rep("ch_MiBP_i_cor_Y1_ter", times = 23),
                   rep("ch_MBzP_i_cor_Y1_ter", times = 23),
                   rep("ch_MEP_i_cor_Y1_ter", times = 23),
                   rep("ch_ohMPHP_i_cor_Y1_ter", times = 23),
                   rep("ch_DINCH_ms_i_cor_Y1_ter", times = 23))) %>%
  select(Pollutants, everything())

results_betadiv_univar_bray_curtis_Y1 %>%            # Visualisation des résultats significatifs univarié
  filter(`Pr(>F)` < 0.05) %>%
  filter(`Explanatory variables` %in% c(explanatory_vars_Y1)) %>%
  View()

results_betadiv_multivar_bray_curtis_Y1 %>%         # Visualisation des résultats significatifs multivarié
  filter(`Pr(>F)` < 0.05) %>%
  filter(`Explanatory variables` %in% c(explanatory_vars_Y1)) %>%
  View()


### Sous-groups ----
#### sous groupe MEP t2 ----
##### 1st vs 2nd tertile ----
bdd_mep_t2_1_2 <- bdd_final_t2 %>%
  filter(mo_MEP_i_cor_t2_ter %in% c("1st tertile", "2nd tertile"))
mep_t2_to_keep_1_2 <- bdd_mep_t2_1_2$ident
bdd_mep_t2_1_2 <- bdd_mep_t2_1_2 %>%
  select(ident, all_of(covariables), mo_MEP_i_cor_t2_ter, all_of(mep_t2_to_keep_1_2))
all_dist_mep_t2_1_2 <- bdd_mep_t2_1_2 %>% select(all_of(mep_t2_to_keep_1_2)) %>% as.dist()

results_betadiv_multivar_bray_curtis_mep_t2_1st_2nd_ter <-
  adonis2(all_dist_mep_t2_1_2 ~
            mo_MEP_i_cor_t2_ter +
            ch_feces_RUN_Y1 +
            ch_feces_age_w_Y1_i +
            po_delmod + ch_food_intro_Y1_3cat_i +
            ch_antibio_Y1_2cat_i +
            mo_par_2cat +
            mo_pets_i +
            ch_sex +
            mo_tob_gr_anyt_yn_n2_i +
            Mo_ETS_anyT_yn1_opt_i +
            ch_ETS_12m_opt36m +
            mo_dipl_3cat_i +
            po_w_kg_3cat +
            po_he_3cat_i +
            ch_w_Y1_3cat_i +
            ch_he_Y1_3cat_i +
            po_gd +
            mo_age +
            mo_bmi_bepr_3cat_i +
            bf_duration_till48w_4cat_i,
          data = bdd_mep_t2_1_2,
          permutations = 999)

results_betadiv_multivar_bray_curtis_mep_t2_1st_2nd_ter <-
  results_betadiv_multivar_bray_curtis_mep_t2_1st_2nd_ter %>%
  rownames_to_column(var = "Explanatory variables") %>%
  mutate(Pollutants = c(rep("mo_MEP_i_cor_t2_ter", times = 23)),
         Sous_groups = c(rep("1st_2nd_ter", times = 23))) %>%
  select(Pollutants, Sous_groups, everything())

##### 1st vs 3rd tertile ----
bdd_mep_t2_1_3 <- bdd_final_t2 %>%
  filter(mo_MEP_i_cor_t2_ter %in% c("1st tertile", "3rd tertile"))
mep_t2_to_keep_1_3 <- bdd_mep_t2_1_3$ident
bdd_mep_t2_1_3 <- bdd_mep_t2_1_3 %>%
  select(ident, all_of(covariables), mo_MEP_i_cor_t2_ter, all_of(mep_t2_to_keep_1_3))
all_dist_mep_t2_1_3 <- bdd_mep_t2_1_3 %>% select(all_of(mep_t2_to_keep_1_3)) %>% as.dist()

results_betadiv_multivar_bray_curtis_mep_t2_1st_3rd_ter <-
  adonis2(all_dist_mep_t2_1_3 ~
            mo_MEP_i_cor_t2_ter +
            ch_feces_RUN_Y1 +
            ch_feces_age_w_Y1_i +
            po_delmod + ch_food_intro_Y1_3cat_i +
            ch_antibio_Y1_2cat_i +
            mo_par_2cat +
            mo_pets_i +
            ch_sex +
            mo_tob_gr_anyt_yn_n2_i +
            Mo_ETS_anyT_yn1_opt_i +
            ch_ETS_12m_opt36m +
            mo_dipl_3cat_i +
            po_w_kg_3cat +
            po_he_3cat_i +
            ch_w_Y1_3cat_i +
            ch_he_Y1_3cat_i +
            po_gd +
            mo_age +
            mo_bmi_bepr_3cat_i +
            bf_duration_till48w_4cat_i,
          data = bdd_mep_t2_1_3,
          permutations = 999)

results_betadiv_multivar_bray_curtis_mep_t2_1st_3rd_ter <-
  results_betadiv_multivar_bray_curtis_mep_t2_1st_3rd_ter %>%
  rownames_to_column(var = "Explanatory variables") %>%
  mutate(Pollutants = c(rep("mo_MEP_i_cor_t2_ter", times = 23)),
         Sous_groups = c(rep("1st_3rd_ter", times = 23))) %>%
  select(Pollutants, Sous_groups, everything())



##### 2nd vs 3rd tertile ----
bdd_mep_t2_2_3 <- bdd_final_t2 %>%
  filter(mo_MEP_i_cor_t2_ter %in% c("2nd tertile", "3rd tertile"))
mep_t2_to_keep_2_3 <- bdd_mep_t2_2_3$ident
bdd_mep_t2_2_3 <- bdd_mep_t2_2_3 %>%
  select(ident, all_of(covariables), mo_MEP_i_cor_t2_ter, all_of(mep_t2_to_keep_2_3))
all_dist_mep_t2_2_3 <- bdd_mep_t2_2_3 %>% select(all_of(mep_t2_to_keep_2_3)) %>% as.dist()

results_betadiv_multivar_bray_curtis_mep_t2_2nd_3rd_ter <-
  adonis2(all_dist_mep_t2_2_3 ~
            mo_MEP_i_cor_t2_ter +
            ch_feces_RUN_Y1 +
            ch_feces_age_w_Y1_i +
            po_delmod + ch_food_intro_Y1_3cat_i +
            ch_antibio_Y1_2cat_i +
            mo_par_2cat +
            mo_pets_i +
            ch_sex +
            mo_tob_gr_anyt_yn_n2_i +
            Mo_ETS_anyT_yn1_opt_i +
            ch_ETS_12m_opt36m +
            mo_dipl_3cat_i +
            po_w_kg_3cat +
            po_he_3cat_i +
            ch_w_Y1_3cat_i +
            ch_he_Y1_3cat_i +
            po_gd +
            mo_age +
            mo_bmi_bepr_3cat_i +
            bf_duration_till48w_4cat_i,
          data = bdd_mep_t2_2_3,
          permutations = 999)

results_betadiv_multivar_bray_curtis_mep_t2_2nd_3rd_ter <-
  results_betadiv_multivar_bray_curtis_mep_t2_2nd_3rd_ter %>%
  rownames_to_column(var = "Explanatory variables") %>%
  mutate(Pollutants = c(rep("mo_MEP_i_cor_t2_ter", times = 23)),
         Sous_groups = c(rep("2nd_3rd_ter", times = 23))) %>%
  select(Pollutants, Sous_groups, everything())

rm(bdd_mep_t2_1_2, bdd_mep_t2_1_3, bdd_mep_t2_2_3,
   mep_t2_to_keep_1_2, mep_t2_to_keep_1_3, mep_t2_to_keep_2_3,
   all_dist_mep_t2_1_2, all_dist_mep_t2_1_3, all_dist_mep_t2_2_3)


results_betadiv_multivar_bray_curtis_mep_t2 <-
  bind_rows(results_betadiv_multivar_bray_curtis_mep_t2_1st_2nd_ter,
            results_betadiv_multivar_bray_curtis_mep_t2_1st_3rd_ter,
            results_betadiv_multivar_bray_curtis_mep_t2_2nd_3rd_ter) %>%
  filter(`Explanatory variables` %in% c("mo_MEP_i_cor_t2_ter", "Residual", "Total"))



#### sous groupe DINCH t3 ----
##### 1st vs 2nd tertile ----
bdd_dinch_t3_1_2 <- bdd_final_t3 %>%
  filter(mo_DINCH_ms_i_cor_t3_ter %in% c("1st tertile", "2nd tertile"))
dinch_t3_to_keep_1_2 <- bdd_dinch_t3_1_2$ident
bdd_dinch_t3_1_2 <- bdd_dinch_t3_1_2 %>%
  select(ident, all_of(covariables), mo_DINCH_ms_i_cor_t3_ter, all_of(dinch_t3_to_keep_1_2))
all_dist_dinch_t3_1_2 <- bdd_dinch_t3_1_2 %>% select(all_of(dinch_t3_to_keep_1_2)) %>% as.dist()

results_betadiv_multivar_bray_curtis_dinch_t3_1st_2nd_ter <-
  adonis2(all_dist_dinch_t3_1_2 ~
            mo_DINCH_ms_i_cor_t3_ter +
            ch_feces_RUN_Y1 +
            ch_feces_age_w_Y1_i +
            po_delmod + ch_food_intro_Y1_3cat_i +
            ch_antibio_Y1_2cat_i +
            mo_par_2cat +
            mo_pets_i +
            ch_sex +
            mo_tob_gr_anyt_yn_n2_i +
            Mo_ETS_anyT_yn1_opt_i +
            ch_ETS_12m_opt36m +
            mo_dipl_3cat_i +
            po_w_kg_3cat +
            po_he_3cat_i +
            ch_w_Y1_3cat_i +
            ch_he_Y1_3cat_i +
            po_gd +
            mo_age +
            mo_bmi_bepr_3cat_i +
            bf_duration_till48w_4cat_i,
          data = bdd_dinch_t3_1_2,
          permutations = 999)

results_betadiv_multivar_bray_curtis_dinch_t3_1st_2nd_ter <-
  results_betadiv_multivar_bray_curtis_dinch_t3_1st_2nd_ter %>%
  rownames_to_column(var = "Explanatory variables") %>%
  mutate(Pollutants = c(rep("mo_DINCH_ms_i_cor_t3_ter", times = 23)),
         Sous_groups = c(rep("1st_2nd_ter", times = 23))) %>%
  select(Pollutants, Sous_groups, everything())

##### 1st vs 3rd tertile ----
bdd_dinch_t3_1_3 <- bdd_final_t3 %>%
  filter(mo_DINCH_ms_i_cor_t3_ter %in% c("1st tertile", "3rd tertile"))
dinch_t3_to_keep_1_3 <- bdd_dinch_t3_1_3$ident
bdd_dinch_t3_1_3 <- bdd_dinch_t3_1_3 %>%
  select(ident, all_of(covariables), mo_DINCH_ms_i_cor_t3_ter, all_of(dinch_t3_to_keep_1_3))
all_dist_dinch_t3_1_3 <- bdd_dinch_t3_1_3 %>% select(all_of(dinch_t3_to_keep_1_3)) %>% as.dist()

results_betadiv_multivar_bray_curtis_dinch_t3_1st_3rd_ter <-
  adonis2(all_dist_dinch_t3_1_3 ~
            mo_DINCH_ms_i_cor_t3_ter +
            ch_feces_RUN_Y1 +
            ch_feces_age_w_Y1_i +
            po_delmod + ch_food_intro_Y1_3cat_i +
            ch_antibio_Y1_2cat_i +
            mo_par_2cat +
            mo_pets_i +
            ch_sex +
            mo_tob_gr_anyt_yn_n2_i +
            Mo_ETS_anyT_yn1_opt_i +
            ch_ETS_12m_opt36m +
            mo_dipl_3cat_i +
            po_w_kg_3cat +
            po_he_3cat_i +
            ch_w_Y1_3cat_i +
            ch_he_Y1_3cat_i +
            po_gd +
            mo_age +
            mo_bmi_bepr_3cat_i +
            bf_duration_till48w_4cat_i,
          data = bdd_dinch_t3_1_3,
          permutations = 999)

results_betadiv_multivar_bray_curtis_dinch_t3_1st_3rd_ter <-
  results_betadiv_multivar_bray_curtis_dinch_t3_1st_3rd_ter %>%
  rownames_to_column(var = "Explanatory variables") %>%
  mutate(Pollutants = c(rep("mo_DINCH_ms_i_cor_t3_ter", times = 23)),
         Sous_groups = c(rep("1st_3rd_ter", times = 23))) %>%
  select(Pollutants, Sous_groups, everything())



##### 2nd vs 3rd tertile ----
bdd_dinch_t3_2_3 <- bdd_final_t3 %>%
  filter(mo_DINCH_ms_i_cor_t3_ter %in% c("2nd tertile", "3rd tertile"))
dinch_t3_to_keep_2_3 <- bdd_dinch_t3_2_3$ident
bdd_dinch_t3_2_3 <- bdd_dinch_t3_2_3 %>%
  select(ident, all_of(covariables), mo_DINCH_ms_i_cor_t3_ter, all_of(dinch_t3_to_keep_2_3))
all_dist_dinch_t3_2_3 <- bdd_dinch_t3_2_3 %>% select(all_of(dinch_t3_to_keep_2_3)) %>% as.dist()

results_betadiv_multivar_bray_curtis_dinch_t3_2nd_3rd_ter <-
  adonis2(all_dist_dinch_t3_2_3 ~
            mo_DINCH_ms_i_cor_t3_ter +
            ch_feces_RUN_Y1 +
            ch_feces_age_w_Y1_i +
            po_delmod + ch_food_intro_Y1_3cat_i +
            ch_antibio_Y1_2cat_i +
            mo_par_2cat +
            mo_pets_i +
            ch_sex +
            mo_tob_gr_anyt_yn_n2_i +
            Mo_ETS_anyT_yn1_opt_i +
            ch_ETS_12m_opt36m +
            mo_dipl_3cat_i +
            po_w_kg_3cat +
            po_he_3cat_i +
            ch_w_Y1_3cat_i +
            ch_he_Y1_3cat_i +
            po_gd +
            mo_age +
            mo_bmi_bepr_3cat_i +
            bf_duration_till48w_4cat_i,
          data = bdd_dinch_t3_2_3,
          permutations = 999)

results_betadiv_multivar_bray_curtis_dinch_t3_2nd_3rd_ter <-
  results_betadiv_multivar_bray_curtis_dinch_t3_2nd_3rd_ter %>%
  rownames_to_column(var = "Explanatory variables") %>%
  mutate(Pollutants = c(rep("mo_DINCH_ms_i_cor_t3_ter", times = 23)),
         Sous_groups = c(rep("2nd_3rd_ter", times = 23))) %>%
  select(Pollutants, Sous_groups, everything())

rm(bdd_dinch_t3_1_2, bdd_dinch_t3_1_3, bdd_dinch_t3_2_3,
   dinch_t3_to_keep_1_2, dinch_t3_to_keep_1_3, dinch_t3_to_keep_2_3,
   all_dist_dinch_t3_1_2, all_dist_dinch_t3_1_3, all_dist_dinch_t3_2_3)


results_betadiv_multivar_bray_curtis_dinch_t3 <-
  bind_rows(results_betadiv_multivar_bray_curtis_dinch_t3_1st_2nd_ter,
            results_betadiv_multivar_bray_curtis_dinch_t3_1st_3rd_ter,
            results_betadiv_multivar_bray_curtis_dinch_t3_2nd_3rd_ter) %>%
  filter(`Explanatory variables` %in% c("mo_DINCH_ms_i_cor_t3_ter", "Residual", "Total"))




#### sous groupe DiNP Y1 ----
##### 1st vs 2nd tertile ----
bdd_dinp_Y1_1_2 <- bdd_final_Y1 %>%
  filter(ch_DiNP_ms_i_cor_Y1_ter %in% c("1st tertile", "2nd tertile"))
dinp_Y1_to_keep_1_2 <- bdd_dinp_Y1_1_2$ident
bdd_dinp_Y1_1_2 <- bdd_dinp_Y1_1_2 %>%
  select(ident, all_of(covariables), ch_DiNP_ms_i_cor_Y1_ter, all_of(dinp_Y1_to_keep_1_2))
all_dist_dinp_Y1_1_2 <- bdd_dinp_Y1_1_2 %>% select(all_of(dinp_Y1_to_keep_1_2)) %>% as.dist()

results_betadiv_multivar_bray_curtis_dinp_Y1_1st_2nd_ter <-
  adonis2(all_dist_dinp_Y1_1_2 ~
            ch_DiNP_ms_i_cor_Y1_ter +
            ch_feces_RUN_Y1 +
            ch_feces_age_w_Y1_i +
            po_delmod + ch_food_intro_Y1_3cat_i +
            ch_antibio_Y1_2cat_i +
            mo_par_2cat +
            mo_pets_i +
            ch_sex +
            mo_tob_gr_anyt_yn_n2_i +
            Mo_ETS_anyT_yn1_opt_i +
            ch_ETS_12m_opt36m +
            mo_dipl_3cat_i +
            po_w_kg_3cat +
            po_he_3cat_i +
            ch_w_Y1_3cat_i +
            ch_he_Y1_3cat_i +
            po_gd +
            mo_age +
            mo_bmi_bepr_3cat_i +
            bf_duration_till48w_4cat_i,
          data = bdd_dinp_Y1_1_2,
          permutations = 999)

results_betadiv_multivar_bray_curtis_dinp_Y1_1st_2nd_ter <-
  results_betadiv_multivar_bray_curtis_dinp_Y1_1st_2nd_ter %>%
  rownames_to_column(var = "Explanatory variables") %>%
  mutate(Pollutants = c(rep("ch_DiNP_ms_i_cor_Y1_ter", times = 23)),
         Sous_groups = c(rep("1st_2nd_ter", times = 23))) %>%
  select(Pollutants, Sous_groups, everything())

##### 1st vs 3rd tertile ----
bdd_dinp_Y1_1_3 <- bdd_final_Y1 %>%
  filter(ch_DiNP_ms_i_cor_Y1_ter %in% c("1st tertile", "3rd tertile"))
dinp_Y1_to_keep_1_3 <- bdd_dinp_Y1_1_3$ident
bdd_dinp_Y1_1_3 <- bdd_dinp_Y1_1_3 %>%
  select(ident, all_of(covariables), ch_DiNP_ms_i_cor_Y1_ter, all_of(dinp_Y1_to_keep_1_3))
all_dist_dinp_Y1_1_3 <- bdd_dinp_Y1_1_3 %>% select(all_of(dinp_Y1_to_keep_1_3)) %>% as.dist()

results_betadiv_multivar_bray_curtis_dinp_Y1_1st_3rd_ter <-
  adonis2(all_dist_dinp_Y1_1_3 ~
            ch_DiNP_ms_i_cor_Y1_ter +
            ch_feces_RUN_Y1 +
            ch_feces_age_w_Y1_i +
            po_delmod + ch_food_intro_Y1_3cat_i +
            ch_antibio_Y1_2cat_i +
            mo_par_2cat +
            mo_pets_i +
            ch_sex +
            mo_tob_gr_anyt_yn_n2_i +
            Mo_ETS_anyT_yn1_opt_i +
            ch_ETS_12m_opt36m +
            mo_dipl_3cat_i +
            po_w_kg_3cat +
            po_he_3cat_i +
            ch_w_Y1_3cat_i +
            ch_he_Y1_3cat_i +
            po_gd +
            mo_age +
            mo_bmi_bepr_3cat_i +
            bf_duration_till48w_4cat_i,
          data = bdd_dinp_Y1_1_3,
          permutations = 999)

results_betadiv_multivar_bray_curtis_dinp_Y1_1st_3rd_ter <-
  results_betadiv_multivar_bray_curtis_dinp_Y1_1st_3rd_ter %>%
  rownames_to_column(var = "Explanatory variables") %>%
  mutate(Pollutants = c(rep("ch_DiNP_ms_i_cor_Y1_ter", times = 23)),
         Sous_groups = c(rep("1st_3rd_ter", times = 23))) %>%
  select(Pollutants, Sous_groups, everything())



##### 2nd vs 3rd tertile ----
bdd_dinp_Y1_2_3 <- bdd_final_Y1 %>%
  filter(ch_DiNP_ms_i_cor_Y1_ter %in% c("2nd tertile", "3rd tertile"))
dinp_Y1_to_keep_2_3 <- bdd_dinp_Y1_2_3$ident
bdd_dinp_Y1_2_3 <- bdd_dinp_Y1_2_3 %>%
  select(ident, all_of(covariables), ch_DiNP_ms_i_cor_Y1_ter, all_of(dinp_Y1_to_keep_2_3))
all_dist_dinp_Y1_2_3 <- bdd_dinp_Y1_2_3 %>% select(all_of(dinp_Y1_to_keep_2_3)) %>% as.dist()

results_betadiv_multivar_bray_curtis_dinp_Y1_2nd_3rd_ter <-
  adonis2(all_dist_dinp_Y1_2_3 ~
            ch_DiNP_ms_i_cor_Y1_ter +
            ch_feces_RUN_Y1 +
            ch_feces_age_w_Y1_i +
            po_delmod + ch_food_intro_Y1_3cat_i +
            ch_antibio_Y1_2cat_i +
            mo_par_2cat +
            mo_pets_i +
            ch_sex +
            mo_tob_gr_anyt_yn_n2_i +
            Mo_ETS_anyT_yn1_opt_i +
            ch_ETS_12m_opt36m +
            mo_dipl_3cat_i +
            po_w_kg_3cat +
            po_he_3cat_i +
            ch_w_Y1_3cat_i +
            ch_he_Y1_3cat_i +
            po_gd +
            mo_age +
            mo_bmi_bepr_3cat_i +
            bf_duration_till48w_4cat_i,
          data = bdd_dinp_Y1_2_3,
          permutations = 999)

results_betadiv_multivar_bray_curtis_dinp_Y1_2nd_3rd_ter <-
  results_betadiv_multivar_bray_curtis_dinp_Y1_2nd_3rd_ter %>%
  rownames_to_column(var = "Explanatory variables") %>%
  mutate(Pollutants = c(rep("ch_DiNP_ms_i_cor_Y1_ter", times = 23)),
         Sous_groups = c(rep("2nd_3rd_ter", times = 23))) %>%
  select(Pollutants, Sous_groups, everything())

rm(bdd_dinp_Y1_1_2, bdd_dinp_Y1_1_3, bdd_dinp_Y1_2_3,
   dinp_Y1_to_keep_1_2, dinp_Y1_to_keep_1_3, dinp_Y1_to_keep_2_3,
   all_dist_dinp_Y1_1_2, all_dist_dinp_Y1_1_3, all_dist_dinp_Y1_2_3)


results_betadiv_multivar_bray_curtis_dinp_Y1 <-
  bind_rows(results_betadiv_multivar_bray_curtis_dinp_Y1_1st_2nd_ter,
            results_betadiv_multivar_bray_curtis_dinp_Y1_1st_3rd_ter,
            results_betadiv_multivar_bray_curtis_dinp_Y1_2nd_3rd_ter) %>%
  filter(`Explanatory variables` %in% c("ch_DiNP_ms_i_cor_Y1_ter", "Residual", "Total"))


#### sous groupe MEP Y1 ----
##### 1st vs 2nd tertile ----
bdd_mep_Y1_1_2 <- bdd_final_Y1 %>%
  filter(ch_MEP_i_cor_Y1_ter %in% c("1st tertile", "2nd tertile"))
mep_Y1_to_keep_1_2 <- bdd_mep_Y1_1_2$ident
bdd_mep_Y1_1_2 <- bdd_mep_Y1_1_2 %>%
  select(ident, all_of(covariables), ch_MEP_i_cor_Y1_ter, all_of(mep_Y1_to_keep_1_2))
all_dist_mep_Y1_1_2 <- bdd_mep_Y1_1_2 %>% select(all_of(mep_Y1_to_keep_1_2)) %>% as.dist()

results_betadiv_multivar_bray_curtis_mep_Y1_1st_2nd_ter <-
  adonis2(all_dist_mep_Y1_1_2 ~
            ch_MEP_i_cor_Y1_ter +
            ch_feces_RUN_Y1 +
            ch_feces_age_w_Y1_i +
            po_delmod + ch_food_intro_Y1_3cat_i +
            ch_antibio_Y1_2cat_i +
            mo_par_2cat +
            mo_pets_i +
            ch_sex +
            mo_tob_gr_anyt_yn_n2_i +
            Mo_ETS_anyT_yn1_opt_i +
            ch_ETS_12m_opt36m +
            mo_dipl_3cat_i +
            po_w_kg_3cat +
            po_he_3cat_i +
            ch_w_Y1_3cat_i +
            ch_he_Y1_3cat_i +
            po_gd +
            mo_age +
            mo_bmi_bepr_3cat_i +
            bf_duration_till48w_4cat_i,
          data = bdd_mep_Y1_1_2,
          permutations = 999)

results_betadiv_multivar_bray_curtis_mep_Y1_1st_2nd_ter <-
  results_betadiv_multivar_bray_curtis_mep_Y1_1st_2nd_ter %>%
  rownames_to_column(var = "Explanatory variables") %>%
  mutate(Pollutants = c(rep("ch_MEP_i_cor_Y1_ter", times = 23)),
         Sous_groups = c(rep("1st_2nd_ter", times = 23))) %>%
  select(Pollutants, Sous_groups, everything())

##### 1st vs 3rd tertile ----
bdd_mep_Y1_1_3 <- bdd_final_Y1 %>%
  filter(ch_MEP_i_cor_Y1_ter %in% c("1st tertile", "3rd tertile"))
mep_Y1_to_keep_1_3 <- bdd_mep_Y1_1_3$ident
bdd_mep_Y1_1_3 <- bdd_mep_Y1_1_3 %>%
  select(ident, all_of(covariables), ch_MEP_i_cor_Y1_ter, all_of(mep_Y1_to_keep_1_3))
all_dist_mep_Y1_1_3 <- bdd_mep_Y1_1_3 %>% select(all_of(mep_Y1_to_keep_1_3)) %>% as.dist()

results_betadiv_multivar_bray_curtis_mep_Y1_1st_3rd_ter <-
  adonis2(all_dist_mep_Y1_1_3 ~
            ch_MEP_i_cor_Y1_ter +
            ch_feces_RUN_Y1 +
            ch_feces_age_w_Y1_i +
            po_delmod + ch_food_intro_Y1_3cat_i +
            ch_antibio_Y1_2cat_i +
            mo_par_2cat +
            mo_pets_i +
            ch_sex +
            mo_tob_gr_anyt_yn_n2_i +
            Mo_ETS_anyT_yn1_opt_i +
            ch_ETS_12m_opt36m +
            mo_dipl_3cat_i +
            po_w_kg_3cat +
            po_he_3cat_i +
            ch_w_Y1_3cat_i +
            ch_he_Y1_3cat_i +
            po_gd +
            mo_age +
            mo_bmi_bepr_3cat_i +
            bf_duration_till48w_4cat_i,
          data = bdd_mep_Y1_1_3,
          permutations = 999)

results_betadiv_multivar_bray_curtis_mep_Y1_1st_3rd_ter <-
  results_betadiv_multivar_bray_curtis_mep_Y1_1st_3rd_ter %>%
  rownames_to_column(var = "Explanatory variables") %>%
  mutate(Pollutants = c(rep("ch_MEP_i_cor_Y1_ter", times = 23)),
         Sous_groups = c(rep("1st_3rd_ter", times = 23))) %>%
  select(Pollutants, Sous_groups, everything())



##### 2nd vs 3rd tertile ----
bdd_mep_Y1_2_3 <- bdd_final_Y1 %>%
  filter(ch_MEP_i_cor_Y1_ter %in% c("2nd tertile", "3rd tertile"))
mep_Y1_to_keep_2_3 <- bdd_mep_Y1_2_3$ident
bdd_mep_Y1_2_3 <- bdd_mep_Y1_2_3 %>%
  select(ident, all_of(covariables), ch_MEP_i_cor_Y1_ter, all_of(mep_Y1_to_keep_2_3))
all_dist_mep_Y1_2_3 <- bdd_mep_Y1_2_3 %>% select(all_of(mep_Y1_to_keep_2_3)) %>% as.dist()

results_betadiv_multivar_bray_curtis_mep_Y1_2nd_3rd_ter <-
  adonis2(all_dist_mep_Y1_2_3 ~
            ch_MEP_i_cor_Y1_ter +
            ch_feces_RUN_Y1 +
            ch_feces_age_w_Y1_i +
            po_delmod + ch_food_intro_Y1_3cat_i +
            ch_antibio_Y1_2cat_i +
            mo_par_2cat +
            mo_pets_i +
            ch_sex +
            mo_tob_gr_anyt_yn_n2_i +
            Mo_ETS_anyT_yn1_opt_i +
            ch_ETS_12m_opt36m +
            mo_dipl_3cat_i +
            po_w_kg_3cat +
            po_he_3cat_i +
            ch_w_Y1_3cat_i +
            ch_he_Y1_3cat_i +
            po_gd +
            mo_age +
            mo_bmi_bepr_3cat_i +
            bf_duration_till48w_4cat_i,
          data = bdd_mep_Y1_2_3,
          permutations = 999)

results_betadiv_multivar_bray_curtis_mep_Y1_2nd_3rd_ter <-
  results_betadiv_multivar_bray_curtis_mep_Y1_2nd_3rd_ter %>%
  rownames_to_column(var = "Explanatory variables") %>%
  mutate(Pollutants = c(rep("ch_MEP_i_cor_Y1_ter", times = 23)),
         Sous_groups = c(rep("2nd_3rd_ter", times = 23))) %>%
  select(Pollutants, Sous_groups, everything())

rm(bdd_mep_Y1_1_2, bdd_mep_Y1_1_3, bdd_mep_Y1_2_3,
   mep_Y1_to_keep_1_2, mep_Y1_to_keep_1_3, mep_Y1_to_keep_2_3,
   all_dist_mep_Y1_1_2, all_dist_mep_Y1_1_3, all_dist_mep_Y1_2_3)


results_betadiv_multivar_bray_curtis_mep_Y1 <-
  bind_rows(results_betadiv_multivar_bray_curtis_mep_Y1_1st_2nd_ter,
            results_betadiv_multivar_bray_curtis_mep_Y1_1st_3rd_ter,
            results_betadiv_multivar_bray_curtis_mep_Y1_2nd_3rd_ter) %>%
  filter(`Explanatory variables` %in% c("ch_MEP_i_cor_Y1_ter", "Residual", "Total"))



#### sous groupe ohMPHP Y1 ----
##### 1st vs 2nd tertile ----
bdd_ohMPHP_Y1_1_2 <- bdd_final_Y1 %>%
  filter(ch_ohMPHP_i_cor_Y1_ter %in% c("1st tertile", "2nd tertile"))
ohMPHP_Y1_to_keep_1_2 <- bdd_ohMPHP_Y1_1_2$ident
bdd_ohMPHP_Y1_1_2 <- bdd_ohMPHP_Y1_1_2 %>%
  select(ident, all_of(covariables), ch_ohMPHP_i_cor_Y1_ter, all_of(ohMPHP_Y1_to_keep_1_2))
all_dist_ohMPHP_Y1_1_2 <- bdd_ohMPHP_Y1_1_2 %>% select(all_of(ohMPHP_Y1_to_keep_1_2)) %>% as.dist()

results_betadiv_multivar_bray_curtis_ohMPHP_Y1_1st_2nd_ter <-
  adonis2(all_dist_ohMPHP_Y1_1_2 ~
            ch_ohMPHP_i_cor_Y1_ter +
            ch_feces_RUN_Y1 +
            ch_feces_age_w_Y1_i +
            po_delmod + ch_food_intro_Y1_3cat_i +
            ch_antibio_Y1_2cat_i +
            mo_par_2cat +
            mo_pets_i +
            ch_sex +
            mo_tob_gr_anyt_yn_n2_i +
            Mo_ETS_anyT_yn1_opt_i +
            ch_ETS_12m_opt36m +
            mo_dipl_3cat_i +
            po_w_kg_3cat +
            po_he_3cat_i +
            ch_w_Y1_3cat_i +
            ch_he_Y1_3cat_i +
            po_gd +
            mo_age +
            mo_bmi_bepr_3cat_i +
            bf_duration_till48w_4cat_i,
          data = bdd_ohMPHP_Y1_1_2,
          permutations = 999)

results_betadiv_multivar_bray_curtis_ohMPHP_Y1_1st_2nd_ter <-
  results_betadiv_multivar_bray_curtis_ohMPHP_Y1_1st_2nd_ter %>%
  rownames_to_column(var = "Explanatory variables") %>%
  mutate(Pollutants = c(rep("ch_ohMPHP_i_cor_Y1_ter", times = 23)),
         Sous_groups = c(rep("1st_2nd_ter", times = 23))) %>%
  select(Pollutants, Sous_groups, everything())

##### 1st vs 3rd tertile ----
bdd_ohMPHP_Y1_1_3 <- bdd_final_Y1 %>%
  filter(ch_ohMPHP_i_cor_Y1_ter %in% c("1st tertile", "3rd tertile"))
ohMPHP_Y1_to_keep_1_3 <- bdd_ohMPHP_Y1_1_3$ident
bdd_ohMPHP_Y1_1_3 <- bdd_ohMPHP_Y1_1_3 %>%
  select(ident, all_of(covariables), ch_ohMPHP_i_cor_Y1_ter, all_of(ohMPHP_Y1_to_keep_1_3))
all_dist_ohMPHP_Y1_1_3 <- bdd_ohMPHP_Y1_1_3 %>% select(all_of(ohMPHP_Y1_to_keep_1_3)) %>% as.dist()

results_betadiv_multivar_bray_curtis_ohMPHP_Y1_1st_3rd_ter <-
  adonis2(all_dist_ohMPHP_Y1_1_3 ~
            ch_ohMPHP_i_cor_Y1_ter +
            ch_feces_RUN_Y1 +
            ch_feces_age_w_Y1_i +
            po_delmod + ch_food_intro_Y1_3cat_i +
            ch_antibio_Y1_2cat_i +
            mo_par_2cat +
            mo_pets_i +
            ch_sex +
            mo_tob_gr_anyt_yn_n2_i +
            Mo_ETS_anyT_yn1_opt_i +
            ch_ETS_12m_opt36m +
            mo_dipl_3cat_i +
            po_w_kg_3cat +
            po_he_3cat_i +
            ch_w_Y1_3cat_i +
            ch_he_Y1_3cat_i +
            po_gd +
            mo_age +
            mo_bmi_bepr_3cat_i +
            bf_duration_till48w_4cat_i,
          data = bdd_ohMPHP_Y1_1_3,
          permutations = 999)

results_betadiv_multivar_bray_curtis_ohMPHP_Y1_1st_3rd_ter <-
  results_betadiv_multivar_bray_curtis_ohMPHP_Y1_1st_3rd_ter %>%
  rownames_to_column(var = "Explanatory variables") %>%
  mutate(Pollutants = c(rep("ch_ohMPHP_i_cor_Y1_ter", times = 23)),
         Sous_groups = c(rep("1st_3rd_ter", times = 23))) %>%
  select(Pollutants, Sous_groups, everything())



##### 2nd vs 3rd tertile ----
bdd_ohMPHP_Y1_2_3 <- bdd_final_Y1 %>%
  filter(ch_ohMPHP_i_cor_Y1_ter %in% c("2nd tertile", "3rd tertile"))
ohMPHP_Y1_to_keep_2_3 <- bdd_ohMPHP_Y1_2_3$ident
bdd_ohMPHP_Y1_2_3 <- bdd_ohMPHP_Y1_2_3 %>%
  select(ident, all_of(covariables), ch_ohMPHP_i_cor_Y1_ter, all_of(ohMPHP_Y1_to_keep_2_3))
all_dist_ohMPHP_Y1_2_3 <- bdd_ohMPHP_Y1_2_3 %>% select(all_of(ohMPHP_Y1_to_keep_2_3)) %>% as.dist()

results_betadiv_multivar_bray_curtis_ohMPHP_Y1_2nd_3rd_ter <-
  adonis2(all_dist_ohMPHP_Y1_2_3 ~
            ch_ohMPHP_i_cor_Y1_ter +
            ch_feces_RUN_Y1 +
            ch_feces_age_w_Y1_i +
            po_delmod + ch_food_intro_Y1_3cat_i +
            ch_antibio_Y1_2cat_i +
            mo_par_2cat +
            mo_pets_i +
            ch_sex +
            mo_tob_gr_anyt_yn_n2_i +
            Mo_ETS_anyT_yn1_opt_i +
            ch_ETS_12m_opt36m +
            mo_dipl_3cat_i +
            po_w_kg_3cat +
            po_he_3cat_i +
            ch_w_Y1_3cat_i +
            ch_he_Y1_3cat_i +
            po_gd +
            mo_age +
            mo_bmi_bepr_3cat_i +
            bf_duration_till48w_4cat_i,
          data = bdd_ohMPHP_Y1_2_3,
          permutations = 999)

results_betadiv_multivar_bray_curtis_ohMPHP_Y1_2nd_3rd_ter <-
  results_betadiv_multivar_bray_curtis_ohMPHP_Y1_2nd_3rd_ter %>%
  rownames_to_column(var = "Explanatory variables") %>%
  mutate(Pollutants = c(rep("ch_ohMPHP_i_cor_Y1_ter", times = 23)),
         Sous_groups = c(rep("2nd_3rd_ter", times = 23))) %>%
  select(Pollutants, Sous_groups, everything())

rm(bdd_ohMPHP_Y1_1_2, bdd_ohMPHP_Y1_1_3, bdd_ohMPHP_Y1_2_3,
   ohMPHP_Y1_to_keep_1_2, ohMPHP_Y1_to_keep_1_3, ohMPHP_Y1_to_keep_2_3,
   all_dist_ohMPHP_Y1_1_2, all_dist_ohMPHP_Y1_1_3, all_dist_ohMPHP_Y1_2_3)


results_betadiv_multivar_bray_curtis_ohMPHP_Y1 <-
  bind_rows(results_betadiv_multivar_bray_curtis_ohMPHP_Y1_1st_2nd_ter,
            results_betadiv_multivar_bray_curtis_ohMPHP_Y1_1st_3rd_ter,
            results_betadiv_multivar_bray_curtis_ohMPHP_Y1_2nd_3rd_ter) %>%
  filter(`Explanatory variables` %in% c("ch_ohMPHP_i_cor_Y1_ter", "Residual", "Total"))




#### sous groupe dinch Y1 ----
##### 1st vs 2nd tertile ----
bdd_dinch_Y1_1_2 <- bdd_final_Y1 %>%
  filter(ch_DINCH_ms_i_cor_Y1_ter %in% c("1st tertile", "2nd tertile"))
dinch_Y1_to_keep_1_2 <- bdd_dinch_Y1_1_2$ident
bdd_dinch_Y1_1_2 <- bdd_dinch_Y1_1_2 %>%
  select(ident, all_of(covariables), ch_DINCH_ms_i_cor_Y1_ter, all_of(dinch_Y1_to_keep_1_2))
all_dist_dinch_Y1_1_2 <- bdd_dinch_Y1_1_2 %>% select(all_of(dinch_Y1_to_keep_1_2)) %>% as.dist()

results_betadiv_multivar_bray_curtis_dinch_Y1_1st_2nd_ter <-
  adonis2(all_dist_dinch_Y1_1_2 ~
            ch_DINCH_ms_i_cor_Y1_ter +
            ch_feces_RUN_Y1 +
            ch_feces_age_w_Y1_i +
            po_delmod + ch_food_intro_Y1_3cat_i +
            ch_antibio_Y1_2cat_i +
            mo_par_2cat +
            mo_pets_i +
            ch_sex +
            mo_tob_gr_anyt_yn_n2_i +
            Mo_ETS_anyT_yn1_opt_i +
            ch_ETS_12m_opt36m +
            mo_dipl_3cat_i +
            po_w_kg_3cat +
            po_he_3cat_i +
            ch_w_Y1_3cat_i +
            ch_he_Y1_3cat_i +
            po_gd +
            mo_age +
            mo_bmi_bepr_3cat_i +
            bf_duration_till48w_4cat_i,
          data = bdd_dinch_Y1_1_2,
          permutations = 999)

results_betadiv_multivar_bray_curtis_dinch_Y1_1st_2nd_ter <-
  results_betadiv_multivar_bray_curtis_dinch_Y1_1st_2nd_ter %>%
  rownames_to_column(var = "Explanatory variables") %>%
  mutate(Pollutants = c(rep("ch_DINCH_ms_i_cor_Y1_ter", times = 23)),
         Sous_groups = c(rep("1st_2nd_ter", times = 23))) %>%
  select(Pollutants, Sous_groups, everything())

##### 1st vs 3rd tertile ----
bdd_dinch_Y1_1_3 <- bdd_final_Y1 %>%
  filter(ch_DINCH_ms_i_cor_Y1_ter %in% c("1st tertile", "3rd tertile"))
dinch_Y1_to_keep_1_3 <- bdd_dinch_Y1_1_3$ident
bdd_dinch_Y1_1_3 <- bdd_dinch_Y1_1_3 %>%
  select(ident, all_of(covariables), ch_DINCH_ms_i_cor_Y1_ter, all_of(dinch_Y1_to_keep_1_3))
all_dist_dinch_Y1_1_3 <- bdd_dinch_Y1_1_3 %>% select(all_of(dinch_Y1_to_keep_1_3)) %>% as.dist()

results_betadiv_multivar_bray_curtis_dinch_Y1_1st_3rd_ter <-
  adonis2(all_dist_dinch_Y1_1_3 ~
            ch_DINCH_ms_i_cor_Y1_ter +
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
            mo_dipl_3cat_i +
            po_w_kg_3cat +
            po_he_3cat_i +
            ch_w_Y1_3cat_i +
            ch_he_Y1_3cat_i +
            po_gd +
            mo_age +
            mo_bmi_bepr_3cat_i +
            bf_duration_till48w_4cat_i,
          data = bdd_dinch_Y1_1_3,
          permutations = 999)

results_betadiv_multivar_bray_curtis_dinch_Y1_1st_3rd_ter <-
  results_betadiv_multivar_bray_curtis_dinch_Y1_1st_3rd_ter %>%
  rownames_to_column(var = "Explanatory variables") %>%
  mutate(Pollutants = c(rep("ch_DINCH_ms_i_cor_Y1_ter", times = 23)),
         Sous_groups = c(rep("1st_3rd_ter", times = 23))) %>%
  select(Pollutants, Sous_groups, everything())



##### 2nd vs 3rd tertile ----
bdd_dinch_Y1_2_3 <- bdd_final_Y1 %>%
  filter(ch_DINCH_ms_i_cor_Y1_ter %in% c("2nd tertile", "3rd tertile"))
dinch_Y1_to_keep_2_3 <- bdd_dinch_Y1_2_3$ident
bdd_dinch_Y1_2_3 <- bdd_dinch_Y1_2_3 %>%
  select(ident, all_of(covariables), ch_DINCH_ms_i_cor_Y1_ter, all_of(dinch_Y1_to_keep_2_3))
all_dist_dinch_Y1_2_3 <- bdd_dinch_Y1_2_3 %>% select(all_of(dinch_Y1_to_keep_2_3)) %>% as.dist()

results_betadiv_multivar_bray_curtis_dinch_Y1_2nd_3rd_ter <-
  adonis2(all_dist_dinch_Y1_2_3 ~
            ch_DINCH_ms_i_cor_Y1_ter +
            ch_feces_RUN_Y1 +
            ch_feces_age_w_Y1_i +
            po_delmod + ch_food_intro_Y1_3cat_i +
            ch_antibio_Y1_2cat_i +
            mo_par_2cat +
            mo_pets_i +
            ch_sex +
            mo_tob_gr_anyt_yn_n2_i +
            Mo_ETS_anyT_yn1_opt_i +
            ch_ETS_12m_opt36m +
            mo_dipl_3cat_i +
            po_w_kg_3cat +
            po_he_3cat_i +
            ch_w_Y1_3cat_i +
            ch_he_Y1_3cat_i +
            po_gd +
            mo_age +
            mo_bmi_bepr_3cat_i +
            bf_duration_till48w_4cat_i,
          data = bdd_dinch_Y1_2_3,
          permutations = 999)

results_betadiv_multivar_bray_curtis_dinch_Y1_2nd_3rd_ter <-
  results_betadiv_multivar_bray_curtis_dinch_Y1_2nd_3rd_ter %>%
  rownames_to_column(var = "Explanatory variables") %>%
  mutate(Pollutants = c(rep("ch_DINCH_ms_i_cor_Y1_ter", times = 23)),
         Sous_groups = c(rep("2nd_3rd_ter", times = 23))) %>%
  select(Pollutants, Sous_groups, everything())

rm(bdd_dinch_Y1_1_2, bdd_dinch_Y1_1_3, bdd_dinch_Y1_2_3,
   dinch_Y1_to_keep_1_2, dinch_Y1_to_keep_1_3, dinch_Y1_to_keep_2_3,
   all_dist_dinch_Y1_1_2, all_dist_dinch_Y1_1_3, all_dist_dinch_Y1_2_3)


results_betadiv_multivar_bray_curtis_dinch_Y1 <-
  bind_rows(results_betadiv_multivar_bray_curtis_dinch_Y1_1st_2nd_ter,
            results_betadiv_multivar_bray_curtis_dinch_Y1_1st_3rd_ter,
            results_betadiv_multivar_bray_curtis_dinch_Y1_2nd_3rd_ter) %>%
  filter(`Explanatory variables` %in% c("ch_DINCH_ms_i_cor_Y1_ter", "Residual", "Total"))



### Assemblage ----
results_betadiv_univ <-
  list(
    results_betadiv_univar_bray_curtis_t2,
    results_betadiv_univar_bray_curtis_t3,
    results_betadiv_univar_bray_curtis_Y1)
results_betadiv_univ <- do.call(rbind, results_betadiv_univ, quote = FALSE)
results_betadiv_univ <- results_betadiv_univ %>%
  mutate(`Explanatory variables` = if_else(`Explanatory variables` %in% c(explanatory_vars_t2,
                                                                          explanatory_vars_t3,
                                                                          explanatory_vars_Y1),
                                           `Explanatory variables`,
                                           str_remove(`Explanatory variables`, "\\d+$")))

results_betadiv_multi <-
  list(
    results_betadiv_multivar_bray_curtis_t2,
    results_betadiv_multivar_bray_curtis_t3,
    results_betadiv_multivar_bray_curtis_Y1)
results_betadiv_multi <- do.call(rbind, results_betadiv_multi, quote = FALSE)
results_betadiv_multi <-
  results_betadiv_multi %>%
  mutate(
    `Explanatory variables` = if_else(`Explanatory variables` %in% c(explanatory_vars_t2,
                                                                     explanatory_vars_t3,
                                                                     explanatory_vars_Y1),
                                      `Explanatory variables`,
                                      str_remove(`Explanatory variables`, "\\d+$"))) %>%
  filter(
    `Explanatory variables` %in% c(
      explanatory_vars_t2,
      explanatory_vars_t3,
      explanatory_vars_Y1) |
      str_detect(`Explanatory variables`, "Residual") |
      str_detect(`Explanatory variables`, "Total"))

results_betadiv_detailled <- 
  bind_rows(results_betadiv_multivar_bray_curtis_mep_t2, 
            results_betadiv_multivar_bray_curtis_dinch_t3, 
            results_betadiv_multivar_bray_curtis_dinp_Y1,
            results_betadiv_multivar_bray_curtis_mep_Y1, 
            results_betadiv_multivar_bray_curtis_ohMPHP_Y1,
            results_betadiv_multivar_bray_curtis_dinch_Y1)

results_betadiv_complet <-
  list(
    univar = list(results_betadiv_univar_bray_curtis_t2 = results_betadiv_univar_bray_curtis_t2,
                  results_betadiv_univar_bray_curtis_t3 = results_betadiv_univar_bray_curtis_t3,
                  results_betadiv_univar_bray_curtis_Y1 = results_betadiv_univar_bray_curtis_Y1),
    multivar = list(results_betadiv_multivar_bray_curtis_t2 = results_betadiv_multivar_bray_curtis_t2,
                    results_betadiv_multivar_bray_curtis_t3 = results_betadiv_multivar_bray_curtis_t3,
                    results_betadiv_multivar_bray_curtis_Y1 = results_betadiv_multivar_bray_curtis_Y1), 
    detailled = list(
      MEP_t2 = list(results_betadiv_multivar_bray_curtis_mep_t2_1st_2nd_ter = 
                       results_betadiv_multivar_bray_curtis_mep_t2_1st_2nd_ter, 
                     results_betadiv_multivar_bray_curtis_mep_t2_1st_3rd_ter = 
                       results_betadiv_multivar_bray_curtis_mep_t2_1st_3rd_ter, 
                     results_betadiv_multivar_bray_curtis_mep_t2_2nd_3rd_ter = 
                       results_betadiv_multivar_bray_curtis_mep_t2_2nd_3rd_ter), 
      DINCH_t3 = list(results_betadiv_multivar_bray_curtis_dinch_t3_1st_2nd_ter = 
                       results_betadiv_multivar_bray_curtis_dinch_t3_1st_2nd_ter, 
                     results_betadiv_multivar_bray_curtis_dinch_t3_1st_3rd_ter = 
                       results_betadiv_multivar_bray_curtis_dinch_t3_1st_3rd_ter, 
                     results_betadiv_multivar_bray_curtis_dinch_t3_2nd_3rd_ter =
                       results_betadiv_multivar_bray_curtis_dinch_t3_2nd_3rd_ter), 
      DiNP_Y1 = list(results_betadiv_multivar_bray_curtis_dinp_Y1_1st_2nd_ter = 
                       results_betadiv_multivar_bray_curtis_dinp_Y1_1st_2nd_ter, 
                     results_betadiv_multivar_bray_curtis_dinp_Y1_1st_3rd_ter =
                       results_betadiv_multivar_bray_curtis_dinp_Y1_1st_3rd_ter, 
                     results_betadiv_multivar_bray_curtis_dinp_Y1_2nd_3rd_ter = 
                       results_betadiv_multivar_bray_curtis_dinp_Y1_2nd_3rd_ter), 
      MEP_Y1 = list(results_betadiv_multivar_bray_curtis_mep_Y1_1st_2nd_ter = 
                        results_betadiv_multivar_bray_curtis_mep_Y1_1st_2nd_ter, 
                      results_betadiv_multivar_bray_curtis_mep_Y1_1st_3rd_ter =
                        results_betadiv_multivar_bray_curtis_mep_Y1_1st_3rd_ter, 
                      results_betadiv_multivar_bray_curtis_mep_Y1_2nd_3rd_ter = 
                        results_betadiv_multivar_bray_curtis_mep_Y1_2nd_3rd_ter), 
      ohMPHP_Y1 = list(results_betadiv_multivar_bray_curtis_ohMPHP_Y1_1st_2nd_ter = 
                        results_betadiv_multivar_bray_curtis_ohMPHP_Y1_1st_2nd_ter, 
                      results_betadiv_multivar_bray_curtis_ohMPHP_Y1_1st_3rd_ter =
                        results_betadiv_multivar_bray_curtis_ohMPHP_Y1_1st_3rd_ter, 
                      results_betadiv_multivar_bray_curtis_ohMPHP_Y1_2nd_3rd_ter = 
                        results_betadiv_multivar_bray_curtis_ohMPHP_Y1_2nd_3rd_ter), 
      DINCH_Y1 = list(results_betadiv_multivar_bray_curtis_dinch_Y1_1st_2nd_ter = 
                        results_betadiv_multivar_bray_curtis_dinch_Y1_1st_2nd_ter, 
                      results_betadiv_multivar_bray_curtis_dinch_Y1_1st_3rd_ter =
                        results_betadiv_multivar_bray_curtis_dinch_Y1_1st_3rd_ter, 
                      results_betadiv_multivar_bray_curtis_dinch_Y1_2nd_3rd_ter = 
                        results_betadiv_multivar_bray_curtis_dinch_Y1_2nd_3rd_ter)))

list <- list(results_betadiv_multi, 
             results_betadiv_univ, 
             results_betadiv_multivar_bray_curtis_mep_t2, 
             results_betadiv_multivar_bray_curtis_dinch_t3, 
             results_betadiv_multivar_bray_curtis_dinp_Y1, 
             results_betadiv_multivar_bray_curtis_mep_Y1, 
             results_betadiv_multivar_bray_curtis_ohMPHP_Y1, 
             results_betadiv_multivar_bray_curtis_dinch_Y1)
write_xlsx(list, path = "4_output/betadiv/results_betadiv.xlsx")

rm(list, 
   
   results_betadiv_univar_bray_curtis_t2, 
   results_betadiv_univar_bray_curtis_t3, 
   results_betadiv_univar_bray_curtis_Y1, 
   
   results_betadiv_multivar_bray_curtis_t2, 
   results_betadiv_multivar_bray_curtis_t3, 
   results_betadiv_multivar_bray_curtis_Y1, 
   
   results_betadiv_multivar_bray_curtis_mep_t2_1st_2nd_ter, 
   results_betadiv_multivar_bray_curtis_mep_t2_1st_3rd_ter, 
   results_betadiv_multivar_bray_curtis_mep_t2_2nd_3rd_ter, 
   
   results_betadiv_multivar_bray_curtis_dinch_t3_1st_2nd_ter, 
   results_betadiv_multivar_bray_curtis_dinch_t3_1st_3rd_ter, 
   results_betadiv_multivar_bray_curtis_dinch_t3_2nd_3rd_ter, 
   
   results_betadiv_multivar_bray_curtis_dinp_Y1_1st_2nd_ter, 
   results_betadiv_multivar_bray_curtis_dinp_Y1_1st_3rd_ter, 
   results_betadiv_multivar_bray_curtis_dinp_Y1_2nd_3rd_ter, 
   
   results_betadiv_multivar_bray_curtis_mep_Y1_1st_2nd_ter, 
   results_betadiv_multivar_bray_curtis_mep_Y1_1st_3rd_ter, 
   results_betadiv_multivar_bray_curtis_mep_Y1_2nd_3rd_ter, 
   
   results_betadiv_multivar_bray_curtis_ohMPHP_Y1_1st_2nd_ter, 
   results_betadiv_multivar_bray_curtis_ohMPHP_Y1_1st_3rd_ter, 
   results_betadiv_multivar_bray_curtis_ohMPHP_Y1_2nd_3rd_ter, 
   
   results_betadiv_multivar_bray_curtis_dinch_Y1_1st_2nd_ter, 
   results_betadiv_multivar_bray_curtis_dinch_Y1_1st_3rd_ter, 
   results_betadiv_multivar_bray_curtis_dinch_Y1_2nd_3rd_ter, 
   
   results_betadiv_multivar_bray_curtis_mep_t2, 
   results_betadiv_multivar_bray_curtis_dinch_t3, 
   results_betadiv_multivar_bray_curtis_dinp_Y1, 
   results_betadiv_multivar_bray_curtis_mep_Y1, 
   results_betadiv_multivar_bray_curtis_ohMPHP_Y1, 
   results_betadiv_multivar_bray_curtis_dinch_Y1)


results_betadiv_univ %>%
  filter(`Pr(>F)` < 0.05) %>%
  filter(`Explanatory variables` %in% c(explanatory_vars_t2, explanatory_vars_t3, explanatory_vars_Y1)) %>%
  View()

results_betadiv_multi %>%
  filter(`Pr(>F)` < 0.05) %>%
  filter(`Explanatory variables` %in% c(explanatory_vars_t2, explanatory_vars_t3, explanatory_vars_Y1)) %>%
  View()

table_2_article <- results_betadiv_multi %>%
  filter(!`Explanatory variables` %in% c("Residual", "Total")) %>%
  select(-`Explanatory variables`, -Df, -R2) %>%
  mutate(
    Pollutants = fct_relevel(Pollutants, 
                             "mo_DEHP_ms_i_cor_t2_ter", "mo_DEHP_ms_i_cor_t3_ter", "ch_DEHP_ms_i_cor_Y1_ter",
                              "mo_MnBP_i_cor_t2_ter", "mo_MnBP_i_cor_t3_ter", "ch_MnBP_i_cor_Y1_ter",
                              "mo_DiNP_ms_i_cor_t2_ter", "mo_DiNP_ms_i_cor_t3_ter", "ch_DiNP_ms_i_cor_Y1_ter",
                              "mo_MiBP_i_cor_t2_ter", "mo_MiBP_i_cor_t3_ter", "ch_MiBP_i_cor_Y1_ter",
                              "mo_MBzP_i_cor_t2_ter", "mo_MBzP_i_cor_t3_ter", "ch_MBzP_i_cor_Y1_ter",
                              "mo_MEP_i_cor_t2_ter", "mo_MEP_i_cor_t3_ter", "ch_MEP_i_cor_Y1_ter",
                              "mo_ohMPHP_i_cor_t2_ter", "mo_ohMPHP_i_cor_t3_ter", "ch_ohMPHP_i_cor_Y1_ter",
                              "mo_DINCH_ms_i_cor_t2_ter", "mo_DINCH_ms_i_cor_t3_ter", "ch_DINCH_ms_i_cor_Y1_ter")) %>%
  arrange(Pollutants)
write_xlsx(table_2_article, 
           path = "4_output/betadiv/table_2_article.xlsx")

## Boxplots ----
### Trim2. ----
# on réduit la matrice de dissimilarité pour ne pas avoir les données en double
metric_bray_curtis_red <- metric_bray_curtis %>% column_to_rownames(var = "ident")
metric_bray_curtis_red[lower.tri(metric_bray_curtis_red)]<- NA
metric_bray_curtis_red <- metric_bray_curtis_red %>% rownames_to_column(var = "ident")

# on merge les données d'expo t2 aux données de dissimilarité
bdd_long_t2 <- bdd_final_t2 %>% select(ident, all_of(explanatory_vars_t2))
bdd_long_t2 <- full_join(bdd_long_t2, metric_bray_curtis_red, by = "ident")

# on fait passer la base de donnnées en long
bdd_long_t2[,explanatory_vars_t2] <- lapply(bdd_long_t2[,explanatory_vars_t2], as.character)
bdd_long_t2 <- bdd_long_t2 %>%
  pivot_longer(cols = c(-"ident", -all_of(explanatory_vars_t2)), 
               names_to = "ident_bis", 
               values_to = "Bray_curtis_dissimilarity") %>%
  select(ident, ident_bis, Bray_curtis_dissimilarity, everything()) %>%
  filter(ident != ident_bis) %>%
  filter(!is.na(Bray_curtis_dissimilarity)) 

# on créé une base de données bis pour savoir à quelle groupes d'expos sont comparés les paires 
bdd_long_t2_bis <- bdd_final_t2 %>%
  select(ident, 
         all_of(explanatory_vars_t2)) %>%
  rename(ident_bis = ident, 
         mo_DEHP_ms_i_cor_t2_ter_bis = mo_DEHP_ms_i_cor_t2_ter, 
         mo_MnBP_i_cor_t2_ter_bis = mo_MnBP_i_cor_t2_ter, 
         mo_DiNP_ms_i_cor_t2_ter_bis = mo_DiNP_ms_i_cor_t2_ter, 
         mo_MiBP_i_cor_t2_ter_bis = mo_MiBP_i_cor_t2_ter, 
         mo_MBzP_i_cor_t2_ter_bis = mo_MBzP_i_cor_t2_ter, 
         mo_MEP_i_cor_t2_ter_bis = mo_MEP_i_cor_t2_ter, 
         mo_ohMPHP_i_cor_t2_ter_bis = mo_ohMPHP_i_cor_t2_ter, 
         mo_DINCH_ms_i_cor_t2_ter_bis = mo_DINCH_ms_i_cor_t2_ter)

# on merge les données puis on supprime la base de données bis 
bdd_long_t2 <- full_join(bdd_long_t2, bdd_long_t2_bis, by = "ident_bis") 
bdd_long_t2 <- bdd_long_t2 %>%
  filter(ident != 15804) %>%
  filter(ident_bis != 15804)
rm(bdd_long_t2_bis)

# on créé une variable qui indique de quel intergroupe il s'agit
bdd_long_t2 <- bdd_long_t2 %>%
  mutate(
    DEHP = case_when(mo_DEHP_ms_i_cor_t2_ter == "1st tertile" & mo_DEHP_ms_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile", 
                     mo_DEHP_ms_i_cor_t2_ter == "2nd tertile" & mo_DEHP_ms_i_cor_t2_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     mo_DEHP_ms_i_cor_t2_ter == "3rd tertile" & mo_DEHP_ms_i_cor_t2_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     mo_DEHP_ms_i_cor_t2_ter == "1st tertile" & mo_DEHP_ms_i_cor_t2_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     mo_DEHP_ms_i_cor_t2_ter == "2nd tertile" & mo_DEHP_ms_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     mo_DEHP_ms_i_cor_t2_ter == "1st tertile" & mo_DEHP_ms_i_cor_t2_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     mo_DEHP_ms_i_cor_t2_ter == "3rd tertile" & mo_DEHP_ms_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     mo_DEHP_ms_i_cor_t2_ter == "2nd tertile" & mo_DEHP_ms_i_cor_t2_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     mo_DEHP_ms_i_cor_t2_ter == "3rd tertile" & mo_DEHP_ms_i_cor_t2_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    MnBP = case_when(mo_MnBP_i_cor_t2_ter == "1st tertile" & mo_MnBP_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile", 
                     mo_MnBP_i_cor_t2_ter == "2nd tertile" & mo_MnBP_i_cor_t2_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     mo_MnBP_i_cor_t2_ter == "3rd tertile" & mo_MnBP_i_cor_t2_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     mo_MnBP_i_cor_t2_ter == "1st tertile" & mo_MnBP_i_cor_t2_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     mo_MnBP_i_cor_t2_ter == "2nd tertile" & mo_MnBP_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     mo_MnBP_i_cor_t2_ter == "1st tertile" & mo_MnBP_i_cor_t2_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     mo_MnBP_i_cor_t2_ter == "3rd tertile" & mo_MnBP_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     mo_MnBP_i_cor_t2_ter == "2nd tertile" & mo_MnBP_i_cor_t2_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     mo_MnBP_i_cor_t2_ter == "3rd tertile" & mo_MnBP_i_cor_t2_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    DiNP = case_when(mo_DiNP_ms_i_cor_t2_ter == "1st tertile" & mo_DiNP_ms_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile", 
                     mo_DiNP_ms_i_cor_t2_ter == "2nd tertile" & mo_DiNP_ms_i_cor_t2_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     mo_DiNP_ms_i_cor_t2_ter == "3rd tertile" & mo_DiNP_ms_i_cor_t2_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     mo_DiNP_ms_i_cor_t2_ter == "1st tertile" & mo_DiNP_ms_i_cor_t2_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     mo_DiNP_ms_i_cor_t2_ter == "2nd tertile" & mo_DiNP_ms_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     mo_DiNP_ms_i_cor_t2_ter == "1st tertile" & mo_DiNP_ms_i_cor_t2_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     mo_DiNP_ms_i_cor_t2_ter == "3rd tertile" & mo_DiNP_ms_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     mo_DiNP_ms_i_cor_t2_ter == "2nd tertile" & mo_DiNP_ms_i_cor_t2_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     mo_DiNP_ms_i_cor_t2_ter == "3rd tertile" & mo_DiNP_ms_i_cor_t2_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    MiBP = case_when(mo_MiBP_i_cor_t2_ter == "1st tertile" & mo_MiBP_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile", 
                     mo_MiBP_i_cor_t2_ter == "2nd tertile" & mo_MiBP_i_cor_t2_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     mo_MiBP_i_cor_t2_ter == "3rd tertile" & mo_MiBP_i_cor_t2_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     mo_MiBP_i_cor_t2_ter == "1st tertile" & mo_MiBP_i_cor_t2_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     mo_MiBP_i_cor_t2_ter == "2nd tertile" & mo_MiBP_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     mo_MiBP_i_cor_t2_ter == "1st tertile" & mo_MiBP_i_cor_t2_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     mo_MiBP_i_cor_t2_ter == "3rd tertile" & mo_MiBP_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     mo_MiBP_i_cor_t2_ter == "2nd tertile" & mo_MiBP_i_cor_t2_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     mo_MiBP_i_cor_t2_ter == "3rd tertile" & mo_MiBP_i_cor_t2_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    MBzP = case_when(mo_MBzP_i_cor_t2_ter == "1st tertile" & mo_MBzP_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile", 
                     mo_MBzP_i_cor_t2_ter == "2nd tertile" & mo_MBzP_i_cor_t2_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     mo_MBzP_i_cor_t2_ter == "3rd tertile" & mo_MBzP_i_cor_t2_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     mo_MBzP_i_cor_t2_ter == "1st tertile" & mo_MBzP_i_cor_t2_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     mo_MBzP_i_cor_t2_ter == "2nd tertile" & mo_MBzP_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     mo_MBzP_i_cor_t2_ter == "1st tertile" & mo_MBzP_i_cor_t2_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     mo_MBzP_i_cor_t2_ter == "3rd tertile" & mo_MBzP_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     mo_MBzP_i_cor_t2_ter == "2nd tertile" & mo_MBzP_i_cor_t2_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     mo_MBzP_i_cor_t2_ter == "3rd tertile" & mo_MBzP_i_cor_t2_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    MEP = case_when(mo_MEP_i_cor_t2_ter == "1st tertile" & mo_MEP_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile", 
                    mo_MEP_i_cor_t2_ter == "2nd tertile" & mo_MEP_i_cor_t2_ter_bis == "2nd tertile" ~ "2nd tertile", 
                    mo_MEP_i_cor_t2_ter == "3rd tertile" & mo_MEP_i_cor_t2_ter_bis == "3rd tertile" ~ "3rd tertile", 
                    
                    mo_MEP_i_cor_t2_ter == "1st tertile" & mo_MEP_i_cor_t2_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                    mo_MEP_i_cor_t2_ter == "2nd tertile" & mo_MEP_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                    
                    mo_MEP_i_cor_t2_ter == "1st tertile" & mo_MEP_i_cor_t2_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                    mo_MEP_i_cor_t2_ter == "3rd tertile" & mo_MEP_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                    
                    mo_MEP_i_cor_t2_ter == "2nd tertile" & mo_MEP_i_cor_t2_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                    mo_MEP_i_cor_t2_ter == "3rd tertile" & mo_MEP_i_cor_t2_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    ohMPHP = case_when(mo_ohMPHP_i_cor_t2_ter == "1st tertile" & mo_ohMPHP_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile", 
                       mo_ohMPHP_i_cor_t2_ter == "2nd tertile" & mo_ohMPHP_i_cor_t2_ter_bis == "2nd tertile" ~ "2nd tertile", 
                       mo_ohMPHP_i_cor_t2_ter == "3rd tertile" & mo_ohMPHP_i_cor_t2_ter_bis == "3rd tertile" ~ "3rd tertile", 
                       
                       mo_ohMPHP_i_cor_t2_ter == "1st tertile" & mo_ohMPHP_i_cor_t2_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                       mo_ohMPHP_i_cor_t2_ter == "2nd tertile" & mo_ohMPHP_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                       
                       mo_ohMPHP_i_cor_t2_ter == "1st tertile" & mo_ohMPHP_i_cor_t2_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                       mo_ohMPHP_i_cor_t2_ter == "3rd tertile" & mo_ohMPHP_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                       
                       mo_ohMPHP_i_cor_t2_ter == "2nd tertile" & mo_ohMPHP_i_cor_t2_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                       mo_ohMPHP_i_cor_t2_ter == "3rd tertile" & mo_ohMPHP_i_cor_t2_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    DINCH = case_when(mo_DINCH_ms_i_cor_t2_ter == "1st tertile" & mo_DINCH_ms_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile", 
                     mo_DINCH_ms_i_cor_t2_ter == "2nd tertile" & mo_DINCH_ms_i_cor_t2_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     mo_DINCH_ms_i_cor_t2_ter == "3rd tertile" & mo_DINCH_ms_i_cor_t2_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     mo_DINCH_ms_i_cor_t2_ter == "1st tertile" & mo_DINCH_ms_i_cor_t2_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     mo_DINCH_ms_i_cor_t2_ter == "2nd tertile" & mo_DINCH_ms_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     mo_DINCH_ms_i_cor_t2_ter == "1st tertile" & mo_DINCH_ms_i_cor_t2_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     mo_DINCH_ms_i_cor_t2_ter == "3rd tertile" & mo_DINCH_ms_i_cor_t2_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     mo_DINCH_ms_i_cor_t2_ter == "2nd tertile" & mo_DINCH_ms_i_cor_t2_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     mo_DINCH_ms_i_cor_t2_ter == "3rd tertile" & mo_DINCH_ms_i_cor_t2_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"))%>%
  select(ident, ident_bis, Bray_curtis_dissimilarity, 
         DEHP, MnBP, DiNP, MiBP, MBzP, MEP, ohMPHP, DINCH)



# on fait passer les données en long 
bdd_long_t2 <- bdd_long_t2 %>%
  pivot_longer(cols = -c("ident", "ident_bis", "Bray_curtis_dissimilarity"), 
               names_to = "Pollutant", 
               values_to = "Groups")


# duplication de 1st tertile, 2nd tertile, 3rd tertile
test <- 
  bdd_long_t2 %>% 
  filter(Groups %in% c("1st tertile", "2nd tertile", "3rd tertile")) %>%  
  mutate(Groups = fct_recode(Groups, 
                             " 1st tertile" = "1st tertile", 
                             " 2nd tertile" = "2nd tertile", 
                             " 3rd tertile" = "3rd tertile"))
bdd_long_t2 <- bind_rows(bdd_long_t2, test)
rm(test)

bdd_long_t2 <- bdd_long_t2 %>%
  mutate(
    Pollutant_rec = fct_recode(Pollutant, 
                               "ΣDEHP trim.2" = "DEHP",
                                "ΣDINCH trim.2" = "DINCH",
                                "ΣDiNP trim.2" = "DiNP",
                                "MBzP trim.2" = "MBzP",
                                "MEP trim.2" = "MEP",
                                "MiBP trim.2" = "MiBP",
                                "MnBP trim.2" = "MnBP",
                                "ohMPHP trim.2" = "ohMPHP"),
    Pollutant = fct_relevel(Pollutant,
                            "DEHP", "MnBP", "DiNP", "MiBP", "MBzP", "MEP", "ohMPHP", "DINCH"),
    Pollutant_rec = fct_relevel(Pollutant_rec,
                            "ΣDEHP trim.2", "MnBP trim.2", "ΣDiNP trim.2", "MiBP trim.2", 
                            "MBzP trim.2", "MEP trim.2", "ohMPHP trim.2", "ΣDINCH trim.2"),
    Groups = fct_relevel(Groups,
                         " 3rd tertile",                         # intra groupe forte exposition 
                         "2nd tertile - 3rd tertile",            # inter groupe moyenne et forte exposition
                         " 2nd tertile",                         # intra groupe moyenne exposition 
                         
                         "3rd tertile",                          # intra groupe forte exposition 
                         "1st tertile - 3rd tertile",            # inter groupe faible et forte exposition 
                         " 1st tertile", 
                         
                         "2nd tertile",                          # intra groupe moyenne exposition 
                         "1st tertile - 2nd tertile",            # inter groupe faible et moyenne exposition 
                         "1st tertile"),                         # intra groupe faible exposition
    Groups_rec = fct_recode(Groups, 
                            "Medium vs High" = " 3rd tertile",
                            "Medium vs High" = "2nd tertile - 3rd tertile",
                            "Medium vs High" = " 2nd tertile",
                            "Low vs High" = "3rd tertile",
                            "Low vs High" = "1st tertile - 3rd tertile",
                            "Low vs High" = " 1st tertile",
                            "Low vs Medium" = "2nd tertile",
                            "Low vs Medium" = "1st tertile - 2nd tertile",
                            "Low vs Medium" = "1st tertile"))



#### boxplot version 1 ----
boxplot_t2_phthalates <- bdd_long_t2 %>%                    # intra groupe faible exposition
  ggplot() +
  aes(
    x = Bray_curtis_dissimilarity,
    y = Groups,
    fill = Groups
  ) +
  geom_boxplot() +  
  scale_fill_manual(
    values = c(`1st tertile` = "#FFE0DE",   # les catégories faibles expositions (intragroupes)
               ` 1st tertile` = "#FFE0DE", 
               
               `1st tertile - 2nd tertile` = "gray70",  # les variances intergroupes
               `1st tertile - 3rd tertile` = "gray70",
               `2nd tertile - 3rd tertile` = "gray70",
               
               
               `2nd tertile` = "#FF8F87",   # les catégorires moyennes expositions (intragroupes)
               ` 2nd tertile` = "#FF8F87",  
               
               `3rd tertile` = "#FF4034",   
               ` 3rd tertile` = "#FF4034")
  ) +
  labs(x = "Bray Curtis dissimilarity", 
       y = "") +
  theme_lucid() +
  theme(legend.position = "none", 
        axis.title = element_text(face = "bold")) +
  facet_wrap(~Pollutant_rec, scales = "free", ncol = 2)

### Trim3. ----
# on réduit la matrice de dissimilarité pour ne pas avoir les données en double
metric_bray_curtis_red <- metric_bray_curtis %>% column_to_rownames(var = "ident")
metric_bray_curtis_red[lower.tri(metric_bray_curtis_red)]<- NA
metric_bray_curtis_red <- metric_bray_curtis_red %>% rownames_to_column(var = "ident")

# on merge les données d'expo t3 aux données de dissimilarité
bdd_long_t3 <- bdd_final_t3 %>% select(ident, all_of(explanatory_vars_t3))
bdd_long_t3 <- full_join(bdd_long_t3, metric_bray_curtis_red, by = "ident")

# on fait passer la base de donnnées en long
bdd_long_t3[,explanatory_vars_t3] <- lapply(bdd_long_t3[,explanatory_vars_t3], as.character)
bdd_long_t3 <- bdd_long_t3 %>%
  pivot_longer(cols = c(-"ident", -all_of(explanatory_vars_t3)), 
               names_to = "ident_bis", 
               values_to = "Bray_curtis_dissimilarity") %>%
  select(ident, ident_bis, Bray_curtis_dissimilarity, everything()) %>%
  filter(ident != ident_bis) %>%
  filter(!is.na(Bray_curtis_dissimilarity)) 

# on créé une base de données bis pour savoir à quelle groupes d'expos sont comparés les paires 
bdd_long_t3_bis <- bdd_final_t3 %>%
  select(ident, 
         all_of(explanatory_vars_t3)) %>%
  rename(ident_bis = ident, 
         mo_DEHP_ms_i_cor_t3_ter_bis = mo_DEHP_ms_i_cor_t3_ter, 
         mo_MnBP_i_cor_t3_ter_bis = mo_MnBP_i_cor_t3_ter, 
         mo_DiNP_ms_i_cor_t3_ter_bis = mo_DiNP_ms_i_cor_t3_ter, 
         mo_MiBP_i_cor_t3_ter_bis = mo_MiBP_i_cor_t3_ter, 
         mo_MBzP_i_cor_t3_ter_bis = mo_MBzP_i_cor_t3_ter, 
         mo_MEP_i_cor_t3_ter_bis = mo_MEP_i_cor_t3_ter, 
         mo_ohMPHP_i_cor_t3_ter_bis = mo_ohMPHP_i_cor_t3_ter, 
         mo_DINCH_ms_i_cor_t3_ter_bis = mo_DINCH_ms_i_cor_t3_ter)

# on merge les données puis on supprime la base de données bis 
bdd_long_t3 <- full_join(bdd_long_t3, bdd_long_t3_bis, by = "ident_bis") 
bdd_long_t3 <- bdd_long_t3 %>%
  filter(!ident %in% c(17827, 15929, 26891, 25668, 23330, 28199)) %>%
  filter(!ident_bis %in% c(17827, 15929, 26891, 25668, 23330, 28199)) %>%
  filter(!is.na(ident))
rm(bdd_long_t3_bis)

# on créé une variable qui indique de quel intergroupe il s'agit
bdd_long_t3 <- bdd_long_t3 %>%
  mutate(
    DEHP = case_when(mo_DEHP_ms_i_cor_t3_ter == "1st tertile" & mo_DEHP_ms_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile", 
                     mo_DEHP_ms_i_cor_t3_ter == "2nd tertile" & mo_DEHP_ms_i_cor_t3_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     mo_DEHP_ms_i_cor_t3_ter == "3rd tertile" & mo_DEHP_ms_i_cor_t3_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     mo_DEHP_ms_i_cor_t3_ter == "1st tertile" & mo_DEHP_ms_i_cor_t3_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     mo_DEHP_ms_i_cor_t3_ter == "2nd tertile" & mo_DEHP_ms_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     mo_DEHP_ms_i_cor_t3_ter == "1st tertile" & mo_DEHP_ms_i_cor_t3_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     mo_DEHP_ms_i_cor_t3_ter == "3rd tertile" & mo_DEHP_ms_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     mo_DEHP_ms_i_cor_t3_ter == "2nd tertile" & mo_DEHP_ms_i_cor_t3_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     mo_DEHP_ms_i_cor_t3_ter == "3rd tertile" & mo_DEHP_ms_i_cor_t3_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    MnBP = case_when(mo_MnBP_i_cor_t3_ter == "1st tertile" & mo_MnBP_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile", 
                     mo_MnBP_i_cor_t3_ter == "2nd tertile" & mo_MnBP_i_cor_t3_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     mo_MnBP_i_cor_t3_ter == "3rd tertile" & mo_MnBP_i_cor_t3_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     mo_MnBP_i_cor_t3_ter == "1st tertile" & mo_MnBP_i_cor_t3_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     mo_MnBP_i_cor_t3_ter == "2nd tertile" & mo_MnBP_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     mo_MnBP_i_cor_t3_ter == "1st tertile" & mo_MnBP_i_cor_t3_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     mo_MnBP_i_cor_t3_ter == "3rd tertile" & mo_MnBP_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     mo_MnBP_i_cor_t3_ter == "2nd tertile" & mo_MnBP_i_cor_t3_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     mo_MnBP_i_cor_t3_ter == "3rd tertile" & mo_MnBP_i_cor_t3_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    DiNP = case_when(mo_DiNP_ms_i_cor_t3_ter == "1st tertile" & mo_DiNP_ms_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile", 
                     mo_DiNP_ms_i_cor_t3_ter == "2nd tertile" & mo_DiNP_ms_i_cor_t3_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     mo_DiNP_ms_i_cor_t3_ter == "3rd tertile" & mo_DiNP_ms_i_cor_t3_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     mo_DiNP_ms_i_cor_t3_ter == "1st tertile" & mo_DiNP_ms_i_cor_t3_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     mo_DiNP_ms_i_cor_t3_ter == "2nd tertile" & mo_DiNP_ms_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     mo_DiNP_ms_i_cor_t3_ter == "1st tertile" & mo_DiNP_ms_i_cor_t3_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     mo_DiNP_ms_i_cor_t3_ter == "3rd tertile" & mo_DiNP_ms_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     mo_DiNP_ms_i_cor_t3_ter == "2nd tertile" & mo_DiNP_ms_i_cor_t3_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     mo_DiNP_ms_i_cor_t3_ter == "3rd tertile" & mo_DiNP_ms_i_cor_t3_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    MiBP = case_when(mo_MiBP_i_cor_t3_ter == "1st tertile" & mo_MiBP_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile", 
                     mo_MiBP_i_cor_t3_ter == "2nd tertile" & mo_MiBP_i_cor_t3_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     mo_MiBP_i_cor_t3_ter == "3rd tertile" & mo_MiBP_i_cor_t3_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     mo_MiBP_i_cor_t3_ter == "1st tertile" & mo_MiBP_i_cor_t3_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     mo_MiBP_i_cor_t3_ter == "2nd tertile" & mo_MiBP_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     mo_MiBP_i_cor_t3_ter == "1st tertile" & mo_MiBP_i_cor_t3_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     mo_MiBP_i_cor_t3_ter == "3rd tertile" & mo_MiBP_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     mo_MiBP_i_cor_t3_ter == "2nd tertile" & mo_MiBP_i_cor_t3_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     mo_MiBP_i_cor_t3_ter == "3rd tertile" & mo_MiBP_i_cor_t3_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    MBzP = case_when(mo_MBzP_i_cor_t3_ter == "1st tertile" & mo_MBzP_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile", 
                     mo_MBzP_i_cor_t3_ter == "2nd tertile" & mo_MBzP_i_cor_t3_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     mo_MBzP_i_cor_t3_ter == "3rd tertile" & mo_MBzP_i_cor_t3_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     mo_MBzP_i_cor_t3_ter == "1st tertile" & mo_MBzP_i_cor_t3_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     mo_MBzP_i_cor_t3_ter == "2nd tertile" & mo_MBzP_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     mo_MBzP_i_cor_t3_ter == "1st tertile" & mo_MBzP_i_cor_t3_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     mo_MBzP_i_cor_t3_ter == "3rd tertile" & mo_MBzP_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     mo_MBzP_i_cor_t3_ter == "2nd tertile" & mo_MBzP_i_cor_t3_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     mo_MBzP_i_cor_t3_ter == "3rd tertile" & mo_MBzP_i_cor_t3_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    MEP = case_when(mo_MEP_i_cor_t3_ter == "1st tertile" & mo_MEP_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile", 
                    mo_MEP_i_cor_t3_ter == "2nd tertile" & mo_MEP_i_cor_t3_ter_bis == "2nd tertile" ~ "2nd tertile", 
                    mo_MEP_i_cor_t3_ter == "3rd tertile" & mo_MEP_i_cor_t3_ter_bis == "3rd tertile" ~ "3rd tertile", 
                    
                    mo_MEP_i_cor_t3_ter == "1st tertile" & mo_MEP_i_cor_t3_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                    mo_MEP_i_cor_t3_ter == "2nd tertile" & mo_MEP_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                    
                    mo_MEP_i_cor_t3_ter == "1st tertile" & mo_MEP_i_cor_t3_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                    mo_MEP_i_cor_t3_ter == "3rd tertile" & mo_MEP_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                    
                    mo_MEP_i_cor_t3_ter == "2nd tertile" & mo_MEP_i_cor_t3_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                    mo_MEP_i_cor_t3_ter == "3rd tertile" & mo_MEP_i_cor_t3_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    ohMPHP = case_when(mo_ohMPHP_i_cor_t3_ter == "1st tertile" & mo_ohMPHP_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile", 
                       mo_ohMPHP_i_cor_t3_ter == "2nd tertile" & mo_ohMPHP_i_cor_t3_ter_bis == "2nd tertile" ~ "2nd tertile", 
                       mo_ohMPHP_i_cor_t3_ter == "3rd tertile" & mo_ohMPHP_i_cor_t3_ter_bis == "3rd tertile" ~ "3rd tertile", 
                       
                       mo_ohMPHP_i_cor_t3_ter == "1st tertile" & mo_ohMPHP_i_cor_t3_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                       mo_ohMPHP_i_cor_t3_ter == "2nd tertile" & mo_ohMPHP_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                       
                       mo_ohMPHP_i_cor_t3_ter == "1st tertile" & mo_ohMPHP_i_cor_t3_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                       mo_ohMPHP_i_cor_t3_ter == "3rd tertile" & mo_ohMPHP_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                       
                       mo_ohMPHP_i_cor_t3_ter == "2nd tertile" & mo_ohMPHP_i_cor_t3_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                       mo_ohMPHP_i_cor_t3_ter == "3rd tertile" & mo_ohMPHP_i_cor_t3_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    
    DINCH = case_when(mo_DINCH_ms_i_cor_t3_ter == "1st tertile" & mo_DINCH_ms_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile", 
                      mo_DINCH_ms_i_cor_t3_ter == "2nd tertile" & mo_DINCH_ms_i_cor_t3_ter_bis == "2nd tertile" ~ "2nd tertile", 
                      mo_DINCH_ms_i_cor_t3_ter == "3rd tertile" & mo_DINCH_ms_i_cor_t3_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                      mo_DINCH_ms_i_cor_t3_ter == "1st tertile" & mo_DINCH_ms_i_cor_t3_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                      mo_DINCH_ms_i_cor_t3_ter == "2nd tertile" & mo_DINCH_ms_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                      mo_DINCH_ms_i_cor_t3_ter == "1st tertile" & mo_DINCH_ms_i_cor_t3_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                      mo_DINCH_ms_i_cor_t3_ter == "3rd tertile" & mo_DINCH_ms_i_cor_t3_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                      mo_DINCH_ms_i_cor_t3_ter == "2nd tertile" & mo_DINCH_ms_i_cor_t3_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                      mo_DINCH_ms_i_cor_t3_ter == "3rd tertile" & mo_DINCH_ms_i_cor_t3_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"))%>%
  select(ident, ident_bis, Bray_curtis_dissimilarity, 
         DEHP, MnBP, DiNP, MiBP, MBzP, MEP, ohMPHP, DINCH)



# on fait passer les données en long 
bdd_long_t3 <- bdd_long_t3 %>%
  pivot_longer(cols = -c("ident", "ident_bis", "Bray_curtis_dissimilarity"), 
               names_to = "Pollutant", 
               values_to = "Groups")


# duplication de 1st tertile, 2nd tertile, 3rd tertile
test <- 
  bdd_long_t3 %>% 
  filter(Groups %in% c("1st tertile", "2nd tertile", "3rd tertile")) %>%  
  filter(!Pollutant == "BPS") %>%
  mutate(Groups = fct_recode(Groups, 
                             " 1st tertile" = "1st tertile", 
                             " 2nd tertile" = "2nd tertile", 
                             " 3rd tertile" = "3rd tertile"))
bdd_long_t3 <- bind_rows(bdd_long_t3, test)
rm(test)

bdd_long_t3 <- bdd_long_t3 %>%
  mutate(
    Pollutant_rec = fct_recode(Pollutant, 
                               "ΣDEHP trim.3" = "DEHP",
                               "ΣDINCH trim.3" = "DINCH",
                               "ΣDiNP trim.3" = "DiNP",
                               "MBzP trim.3" = "MBzP",
                               "MEP trim.3" = "MEP",
                               "MiBP trim.3" = "MiBP",
                               "MnBP trim.3" = "MnBP",
                               "ohMPHP trim.3" = "ohMPHP"),
    Pollutant = fct_relevel(Pollutant,
                            "DEHP", "MnBP", "DiNP", "MiBP", "MBzP", "MEP", "ohMPHP", "DINCH"),
    Pollutant_rec = fct_relevel(Pollutant_rec,
                                "ΣDEHP trim.3", "MnBP trim.3", "ΣDiNP trim.3", "MiBP trim.3", 
                                "MBzP trim.3", "MEP trim.3", "ohMPHP trim.3", "ΣDINCH trim.3"),
    Groups = fct_relevel(Groups,
                         " 3rd tertile",                             # intra groupe forte exposition 
                         "2nd tertile - 3rd tertile",             # inter groupe moyenne et forte exposition
                         " 2nd tertile",                                  # intra groupe moyenne exposition 
                         
                         "3rd tertile",                                    # intra groupe forte exposition 
                         "1st tertile - 3rd tertile", # inter groupe faible et forte exposition 
                         " 1st tertile",
                         
                         "2nd tertile",                                  # intra groupe moyenne exposition 
                         "1st tertile - 2nd tertile",             # inter groupe faible et moyenne exposition 
                         "1st tertile"),                               # intra groupe faible exposition
    Groups_rec = fct_recode(Groups, 
                            "Medium vs High" = " 3rd tertile",
                            "Medium vs High" = "2nd tertile - 3rd tertile",
                            "Medium vs High" = " 2nd tertile",
                            "Low vs High" = "3rd tertile",
                            "Low vs High" = "1st tertile - 3rd tertile",
                            "Low vs High" = " 1st tertile",
                            "Low vs Medium" = "2nd tertile",
                            "Low vs Medium" = "1st tertile - 2nd tertile",
                            "Low vs Medium" = "1st tertile"))



#### boxplot version 1 ----
boxplot_t3_phthalates <- bdd_long_t3 %>%                            # intra groupe faible exposition
  ggplot() +
  aes(
    x = Bray_curtis_dissimilarity,
    y = Groups,
    fill = Groups
  ) +
  geom_boxplot() +  
  scale_fill_manual(
    values = c(`1st tertile` = "#FFE0DE",   # les catégories faibles expositions (intragroupes)
               ` 1st tertile` = "#FFE0DE", 
               
               `1st tertile - 2nd tertile` = "gray70",  # les variances intergroupes
               `1st tertile - 3rd tertile` = "gray70",
               `2nd tertile - 3rd tertile` = "gray70",
               
               `2nd tertile` = "#FF8F87",   # les catégorires moyennes expositions (intragroupes)
               ` 2nd tertile` = "#FF8F87",   
               
               `3rd tertile` = "#FF4034",   
               ` 3rd tertile` = "#FF4034")
  ) +
  labs(x = "Bray Curtis dissimilarity", 
       y = "") +
  theme_lucid() +
  theme(legend.position = "none", 
        axis.title = element_text(face = "bold")) +
  facet_wrap(~Pollutant_rec, scales = "free", ncol = 2)




### 1 year ----
# on réduit la matrice de dissimilarité pour ne pas avoir les données en double
metric_bray_curtis_red <- metric_bray_curtis %>% column_to_rownames(var = "ident")
metric_bray_curtis_red[lower.tri(metric_bray_curtis_red)]<- NA
metric_bray_curtis_red <- metric_bray_curtis_red %>% rownames_to_column(var = "ident")

# on merge les données d'expo t3 aux données de dissimilarité
bdd_long_Y1 <- bdd_final_Y1 %>% select(ident, all_of(explanatory_vars_Y1))
bdd_long_Y1 <- full_join(bdd_long_Y1, metric_bray_curtis_red, by = "ident")

# on fait passer la base de donnnées en long
bdd_long_Y1[,explanatory_vars_Y1] <- lapply(bdd_long_Y1[,explanatory_vars_Y1], as.character)
bdd_long_Y1 <- bdd_long_Y1 %>%
  pivot_longer(cols = c(-"ident", -all_of(explanatory_vars_Y1)), 
               names_to = "ident_bis", 
               values_to = "Bray_curtis_dissimilarity") %>%
  select(ident, ident_bis, Bray_curtis_dissimilarity, everything()) %>%
  filter(ident != ident_bis) %>%
  filter(!is.na(Bray_curtis_dissimilarity)) 

# on créé une base de données bis pour savoir à quelle groupes d'expos sont comparés les paires 
bdd_long_Y1_bis <- bdd_final_Y1 %>%
  select(ident, 
         all_of(explanatory_vars_Y1)) %>%
  rename(ident_bis = ident, 
         ch_DEHP_ms_i_cor_Y1_ter_bis = ch_DEHP_ms_i_cor_Y1_ter, 
         ch_MnBP_i_cor_Y1_ter_bis = ch_MnBP_i_cor_Y1_ter, 
         ch_DiNP_ms_i_cor_Y1_ter_bis = ch_DiNP_ms_i_cor_Y1_ter, 
         ch_MiBP_i_cor_Y1_ter_bis = ch_MiBP_i_cor_Y1_ter, 
         ch_MBzP_i_cor_Y1_ter_bis = ch_MBzP_i_cor_Y1_ter, 
         ch_MEP_i_cor_Y1_ter_bis = ch_MEP_i_cor_Y1_ter, 
         ch_ohMPHP_i_cor_Y1_ter_bis = ch_ohMPHP_i_cor_Y1_ter, 
         ch_DINCH_ms_i_cor_Y1_ter_bis = ch_DINCH_ms_i_cor_Y1_ter)

# on merge les données puis on supprime la base de données bis 
bdd_long_Y1 <- full_join(bdd_long_Y1, bdd_long_Y1_bis, by = "ident_bis") 
bdd_long_Y1 <- bdd_long_Y1 %>%
  filter(!ident %in% c(23994, 25166, 26766, 14668, 26923)) %>%
  filter(!ident_bis %in% c(23994, 25166, 26766, 14668, 26923)) %>%
  filter(!is.na(ident))
rm(bdd_long_Y1_bis)

# on créé une variable qui indique de quel intergroupe il s'agit
bdd_long_Y1 <- bdd_long_Y1 %>%
  mutate(
    DEHP = case_when(ch_DEHP_ms_i_cor_Y1_ter == "1st tertile" & ch_DEHP_ms_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile", 
                     ch_DEHP_ms_i_cor_Y1_ter == "2nd tertile" & ch_DEHP_ms_i_cor_Y1_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     ch_DEHP_ms_i_cor_Y1_ter == "3rd tertile" & ch_DEHP_ms_i_cor_Y1_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     ch_DEHP_ms_i_cor_Y1_ter == "1st tertile" & ch_DEHP_ms_i_cor_Y1_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     ch_DEHP_ms_i_cor_Y1_ter == "2nd tertile" & ch_DEHP_ms_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     ch_DEHP_ms_i_cor_Y1_ter == "1st tertile" & ch_DEHP_ms_i_cor_Y1_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     ch_DEHP_ms_i_cor_Y1_ter == "3rd tertile" & ch_DEHP_ms_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     ch_DEHP_ms_i_cor_Y1_ter == "2nd tertile" & ch_DEHP_ms_i_cor_Y1_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     ch_DEHP_ms_i_cor_Y1_ter == "3rd tertile" & ch_DEHP_ms_i_cor_Y1_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    MnBP = case_when(ch_MnBP_i_cor_Y1_ter == "1st tertile" & ch_MnBP_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile", 
                     ch_MnBP_i_cor_Y1_ter == "2nd tertile" & ch_MnBP_i_cor_Y1_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     ch_MnBP_i_cor_Y1_ter == "3rd tertile" & ch_MnBP_i_cor_Y1_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     ch_MnBP_i_cor_Y1_ter == "1st tertile" & ch_MnBP_i_cor_Y1_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     ch_MnBP_i_cor_Y1_ter == "2nd tertile" & ch_MnBP_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     ch_MnBP_i_cor_Y1_ter == "1st tertile" & ch_MnBP_i_cor_Y1_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     ch_MnBP_i_cor_Y1_ter == "3rd tertile" & ch_MnBP_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     ch_MnBP_i_cor_Y1_ter == "2nd tertile" & ch_MnBP_i_cor_Y1_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     ch_MnBP_i_cor_Y1_ter == "3rd tertile" & ch_MnBP_i_cor_Y1_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    DiNP = case_when(ch_DiNP_ms_i_cor_Y1_ter == "1st tertile" & ch_DiNP_ms_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile", 
                     ch_DiNP_ms_i_cor_Y1_ter == "2nd tertile" & ch_DiNP_ms_i_cor_Y1_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     ch_DiNP_ms_i_cor_Y1_ter == "3rd tertile" & ch_DiNP_ms_i_cor_Y1_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     ch_DiNP_ms_i_cor_Y1_ter == "1st tertile" & ch_DiNP_ms_i_cor_Y1_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     ch_DiNP_ms_i_cor_Y1_ter == "2nd tertile" & ch_DiNP_ms_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     ch_DiNP_ms_i_cor_Y1_ter == "1st tertile" & ch_DiNP_ms_i_cor_Y1_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     ch_DiNP_ms_i_cor_Y1_ter == "3rd tertile" & ch_DiNP_ms_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     ch_DiNP_ms_i_cor_Y1_ter == "2nd tertile" & ch_DiNP_ms_i_cor_Y1_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     ch_DiNP_ms_i_cor_Y1_ter == "3rd tertile" & ch_DiNP_ms_i_cor_Y1_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    MiBP = case_when(ch_MiBP_i_cor_Y1_ter == "1st tertile" & ch_MiBP_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile", 
                     ch_MiBP_i_cor_Y1_ter == "2nd tertile" & ch_MiBP_i_cor_Y1_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     ch_MiBP_i_cor_Y1_ter == "3rd tertile" & ch_MiBP_i_cor_Y1_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     ch_MiBP_i_cor_Y1_ter == "1st tertile" & ch_MiBP_i_cor_Y1_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     ch_MiBP_i_cor_Y1_ter == "2nd tertile" & ch_MiBP_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     ch_MiBP_i_cor_Y1_ter == "1st tertile" & ch_MiBP_i_cor_Y1_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     ch_MiBP_i_cor_Y1_ter == "3rd tertile" & ch_MiBP_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     ch_MiBP_i_cor_Y1_ter == "2nd tertile" & ch_MiBP_i_cor_Y1_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     ch_MiBP_i_cor_Y1_ter == "3rd tertile" & ch_MiBP_i_cor_Y1_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    MBzP = case_when(ch_MBzP_i_cor_Y1_ter == "1st tertile" & ch_MBzP_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile", 
                     ch_MBzP_i_cor_Y1_ter == "2nd tertile" & ch_MBzP_i_cor_Y1_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     ch_MBzP_i_cor_Y1_ter == "3rd tertile" & ch_MBzP_i_cor_Y1_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     ch_MBzP_i_cor_Y1_ter == "1st tertile" & ch_MBzP_i_cor_Y1_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     ch_MBzP_i_cor_Y1_ter == "2nd tertile" & ch_MBzP_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     ch_MBzP_i_cor_Y1_ter == "1st tertile" & ch_MBzP_i_cor_Y1_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     ch_MBzP_i_cor_Y1_ter == "3rd tertile" & ch_MBzP_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     ch_MBzP_i_cor_Y1_ter == "2nd tertile" & ch_MBzP_i_cor_Y1_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     ch_MBzP_i_cor_Y1_ter == "3rd tertile" & ch_MBzP_i_cor_Y1_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    MEP = case_when(ch_MEP_i_cor_Y1_ter == "1st tertile" & ch_MEP_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile", 
                    ch_MEP_i_cor_Y1_ter == "2nd tertile" & ch_MEP_i_cor_Y1_ter_bis == "2nd tertile" ~ "2nd tertile", 
                    ch_MEP_i_cor_Y1_ter == "3rd tertile" & ch_MEP_i_cor_Y1_ter_bis == "3rd tertile" ~ "3rd tertile", 
                    
                    ch_MEP_i_cor_Y1_ter == "1st tertile" & ch_MEP_i_cor_Y1_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                    ch_MEP_i_cor_Y1_ter == "2nd tertile" & ch_MEP_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                    
                    ch_MEP_i_cor_Y1_ter == "1st tertile" & ch_MEP_i_cor_Y1_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                    ch_MEP_i_cor_Y1_ter == "3rd tertile" & ch_MEP_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                    
                    ch_MEP_i_cor_Y1_ter == "2nd tertile" & ch_MEP_i_cor_Y1_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                    ch_MEP_i_cor_Y1_ter == "3rd tertile" & ch_MEP_i_cor_Y1_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    ohMPHP = case_when(ch_ohMPHP_i_cor_Y1_ter == "1st tertile" & ch_ohMPHP_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile", 
                       ch_ohMPHP_i_cor_Y1_ter == "2nd tertile" & ch_ohMPHP_i_cor_Y1_ter_bis == "2nd tertile" ~ "2nd tertile", 
                       ch_ohMPHP_i_cor_Y1_ter == "3rd tertile" & ch_ohMPHP_i_cor_Y1_ter_bis == "3rd tertile" ~ "3rd tertile", 
                       
                       ch_ohMPHP_i_cor_Y1_ter == "1st tertile" & ch_ohMPHP_i_cor_Y1_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                       ch_ohMPHP_i_cor_Y1_ter == "2nd tertile" & ch_ohMPHP_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                       
                       ch_ohMPHP_i_cor_Y1_ter == "1st tertile" & ch_ohMPHP_i_cor_Y1_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                       ch_ohMPHP_i_cor_Y1_ter == "3rd tertile" & ch_ohMPHP_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                       
                       ch_ohMPHP_i_cor_Y1_ter == "2nd tertile" & ch_ohMPHP_i_cor_Y1_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                       ch_ohMPHP_i_cor_Y1_ter == "3rd tertile" & ch_ohMPHP_i_cor_Y1_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"), 
    
    DINCH = case_when(ch_DINCH_ms_i_cor_Y1_ter == "1st tertile" & ch_DINCH_ms_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile", 
                     ch_DINCH_ms_i_cor_Y1_ter == "2nd tertile" & ch_DINCH_ms_i_cor_Y1_ter_bis == "2nd tertile" ~ "2nd tertile", 
                     ch_DINCH_ms_i_cor_Y1_ter == "3rd tertile" & ch_DINCH_ms_i_cor_Y1_ter_bis == "3rd tertile" ~ "3rd tertile", 
                     
                     ch_DINCH_ms_i_cor_Y1_ter == "1st tertile" & ch_DINCH_ms_i_cor_Y1_ter_bis == "2nd tertile" ~ "1st tertile - 2nd tertile", 
                     ch_DINCH_ms_i_cor_Y1_ter == "2nd tertile" & ch_DINCH_ms_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile - 2nd tertile", 
                     
                     ch_DINCH_ms_i_cor_Y1_ter == "1st tertile" & ch_DINCH_ms_i_cor_Y1_ter_bis == "3rd tertile" ~ "1st tertile - 3rd tertile", 
                     ch_DINCH_ms_i_cor_Y1_ter == "3rd tertile" & ch_DINCH_ms_i_cor_Y1_ter_bis == "1st tertile" ~ "1st tertile - 3rd tertile", 
                     
                     ch_DINCH_ms_i_cor_Y1_ter == "2nd tertile" & ch_DINCH_ms_i_cor_Y1_ter_bis == "3rd tertile" ~ "2nd tertile - 3rd tertile", 
                     ch_DINCH_ms_i_cor_Y1_ter == "3rd tertile" & ch_DINCH_ms_i_cor_Y1_ter_bis == "2nd tertile" ~ "2nd tertile - 3rd tertile"))%>%
  select(ident, ident_bis, Bray_curtis_dissimilarity, 
         DEHP, MnBP, DiNP, MiBP, MBzP, MEP, ohMPHP, DINCH)



# on fait passer les données en long 
bdd_long_Y1 <- bdd_long_Y1 %>%
  pivot_longer(cols = -c("ident", "ident_bis", "Bray_curtis_dissimilarity"), 
               names_to = "Pollutant", 
               values_to = "Groups")


# duplication de 1st tertile, 2nd tertile, 3rd tertile
test <- 
  bdd_long_Y1 %>% 
  filter(Groups %in% c("1st tertile", "2nd tertile", "3rd tertile")) %>%  
  mutate(Groups = fct_recode(Groups, 
                             " 1st tertile" = "1st tertile", 
                             " 2nd tertile" = "2nd tertile", 
                             " 3rd tertile" = "3rd tertile"))
bdd_long_Y1 <- bind_rows(bdd_long_Y1, test)
rm(test)

bdd_long_Y1 <- bdd_long_Y1 %>%
    mutate(
      Pollutant_rec = fct_recode(Pollutant, 
                                 "ΣDEHP 12 months" = "DEHP",
                                 "ΣDINCH 12 months" = "DINCH",
                                 "ΣDiNP 12 months" = "DiNP",
                                 "MBzP 12 months" = "MBzP",
                                 "MEP 12 months" = "MEP",
                                 "MiBP 12 months" = "MiBP",
                                 "MnBP 12 months" = "MnBP",
                                 "ohMPHP 12 months" = "ohMPHP"),
      Pollutant = fct_relevel(Pollutant,
                              "DEHP", "MnBP", "DiNP", "MiBP", "MBzP", "MEP", "ohMPHP", "DINCH"),
      Pollutant_rec = fct_relevel(Pollutant_rec,
                                  "ΣDEHP 12 months", "MnBP 12 months", "ΣDiNP 12 months", "MiBP 12 months", 
                                  "MBzP 12 months", "MEP 12 months", "ohMPHP 12 months", "ΣDINCH 12 months"),
    Groups = fct_relevel(Groups,
                         " 3rd tertile",                           # intra groupe forte exposition 
                         "2nd tertile - 3rd tertile",             # inter groupe moyenne et forte exposition
                         " 2nd tertile",                                   # intra groupe moyenne exposition 
                         
                         "3rd tertile",                                    # intra groupe forte exposition 
                         "1st tertile - 3rd tertile", # inter groupe faible et forte exposition 
                         " 1st tertile",
                         
                         "2nd tertile",                                  # intra groupe moyenne exposition 
                         "1st tertile - 2nd tertile",             # inter groupe faible et moyenne exposition 
                         "1st tertile"))                               # intra groupe faible exposition



#### boxplot version 1 ----
boxplot_Y1_phthalates <- bdd_long_Y1 %>%                            # intra groupe faible exposition
  ggplot() +
  aes(
    x = Bray_curtis_dissimilarity,
    y = Groups,
    fill = Groups
  ) +
  geom_boxplot() +  
  scale_fill_manual(
    values = c(`1st tertile` = "#FFE0DE",   # les catégories faibles expositions (intragroupes)
               ` 1st tertile` = "#FFE0DE", 
               
               `1st tertile - 2nd tertile` = "gray70",  # les variances intergroupes
               `1st tertile - 3rd tertile` = "gray70",
               `2nd tertile - 3rd tertile` = "gray70",
               
               
               `2nd tertile` = "#FF8F87",   # les catégorires moyennes expositions (intragroupes)
               ` 2nd tertile` = "#FF8F87",   
               
               `3rd tertile` = "#FF4034",  
               ` 3rd tertile` = "#FF4034")
  ) +
  labs(x = "Bray Curtis dissimilarity", 
       y = "") +
  theme_lucid() +
  theme(legend.position = "none", 
        axis.title = element_text(face = "bold")) +
  facet_wrap(~Pollutant_rec, scales = "free", ncol = 2)

ggsave(plot = boxplot_t2_phthalates, 
       filename = "4_output/betadiv/boxplot_t2_phthalates.tiff", 
       device = "tiff", 
       units = "cm", 
       width = 25, 
       height = 25)

ggsave(plot = boxplot_t3_phthalates, 
       filename = "4_output/betadiv/boxplot_t3_phthalates.tiff", 
       device = "tiff", 
       units = "cm", 
       width = 25, 
       height = 25)

ggsave(plot = boxplot_Y1_phthalates, 
       filename = "4_output/betadiv/boxplot_Y1_phthalates.tiff", 
       device = "tiff", 
       units = "cm", 
       width = 25, 
       height = 25)


## Detailled (fig4 pour article) ----
### Version 1 : benzo M2 + mepa Y1 + prpa Y1 ----
#### 1st vs 2nd tertile ----
bdd_boxplot_article <- bdd_long_t2 %>%
  filter(Pollutant_rec == "MEP trim.2") 

bdd_boxplot_article_bis <- bdd_long_t3 %>%
  filter(Pollutant_rec == "ΣDINCH trim.3")

bdd_boxplot_article_bis_bis <- bdd_long_Y1 %>%
  filter(Pollutant_rec %in% c("ΣDiNP 12 months", "MEP 12 months", "ohMPHP 12 months", "ΣDINCH 12 months"))

bdd_boxplot_article <- 
  bind_rows(bdd_boxplot_article, bdd_boxplot_article_bis, bdd_boxplot_article_bis_bis) %>%
  filter(!Groups %in% c(" 1st tertile", " 2nd tertile", " 3rd tertile")) %>%
  mutate(
    Groups = fct_relevel(Groups,       
                         "3rd tertile", "2nd tertile - 3rd tertile", " 3rd tertile", " 2nd tertile",
                         "1st tertile - 3rd tertile", " 1st tertile", "2nd tertile",
                         "1st tertile - 2nd tertile", "1st tertile"))

rm(bdd_boxplot_article_bis, bdd_boxplot_article_bis_bis)

boxplot_1 <- bdd_boxplot_article %>%
  filter(Groups %in% c("1st tertile", "1st tertile - 2nd tertile", "2nd tertile")) %>%
  ggplot() +
  aes(
    x = Bray_curtis_dissimilarity,
    y = Groups,
    fill = Groups
  ) +
  geom_boxplot() +
  #labs(x = "Bray Curtis dissimilarity") +
  theme_lucid() +
  theme(legend.position = "none", 
        axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.line.x = element_blank(), 
        strip.text = element_text(size = 12)) +
  facet_wrap(~Pollutant_rec, ncol = 6) +
  scale_fill_manual(
    values = c(`1st tertile` = "#FFE0DE",   # les catégories faibles expositions (intragroupes)
               `1st tertile - 2nd tertile` = "gray70",  # les variances intergroupes
               `1st tertile - 3rd tertile` = "gray70",
               `2nd tertile - 3rd tertile` = "gray70",
               `2nd tertile` = "#FF8F87",   # les catégorires moyennes expositions (intragroupes) 
               `3rd tertile` = "#FF4034"))+
  scale_y_discrete(labels = function(x) {
    lapply(x, function(lbl) {
      if (grepl("1st tertile - 2nd tertile", lbl)) {
        bquote("1"^st* " tertile - 2"^nd* " tertile")
      } else if (grepl("1st", lbl)) {
        bquote("1"^st* " tertile")
      } else if (grepl("2nd", lbl)) {
        bquote("2"^nd* " tertile")
      } else {
        lbl
      }
    })
  }) 

results_betadiv_detailled %>% 
  filter(Sous_groups == "1st_2nd_ter") %>% 
  filter(!`Explanatory variables` %in% c("Residual", "Total")) %>% 
  View()
my_tag <- c("0.10\n(1.41)", 
            "0.09\n(1.40)", 
            "0.02\n(1.77)", 
            "0.07\n(1.41)", 
            "0.02\n(1.62)", 
            "0.05\n(1.50)") 
boxplot_1 <- tag_facet(boxplot_1, 
          x = 1.05, y = 2, 
          vjust = 0.5, hjust = 0.5,
          open = "", close = "",
          fontface = 1,
          size = 3.5,
          #family = "serif",
          tag_pool = my_tag) +
  theme_lucid() +
  theme(legend.position = "none", 
        axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.line.x = element_blank(), 
        strip.text = element_text(size = 12)) 

my_tag_bis <- c("p-value\n(F statistic)", 
                "p-value\n(F statistic)", 
                "p-value\n(F statistic)", 
                "p-value\n(F statistic)", 
                "p-value\n(F statistic)", 
                "p-value\n(F statistic)")
boxplot_1 <- tag_facet(boxplot_1, 
                       x = 1.07, y = 3.1, 
                       vjust = 0.5, hjust = 0.75,
                       open = "", close = "",
                       fontface = 1,
                       size = 3,
                       #family = "serif",
                       tag_pool = my_tag_bis) +
  theme_lucid() +
  theme(legend.position = "none", 
        axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.line.x = element_blank(), 
        strip.text = element_text(size = 12)) 


boxplot_1

#### 1st vs 3rd tertile ----
boxplot_2 <- bdd_boxplot_article %>%
  filter(Groups %in% c("1st tertile", "1st tertile - 3rd tertile", "3rd tertile")) %>%
  ggplot() +
  aes(
    x = Bray_curtis_dissimilarity,
    y = Groups,
    fill = Groups
  ) +
  geom_boxplot() +
  #labs(x = "Bray Curtis dissimilarity") +
  theme_lucid() +
  theme(legend.position = "none", 
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.text.x = element_blank(), 
        axis.line.x = element_blank(),
        strip.background = element_blank(),  # Supprime le fond des facettes
        strip.text = element_blank()) +     # Supprime le texte des facettes) 
  facet_wrap(~Pollutant_rec, ncol = 6) +
  scale_fill_manual(
    values = c(`1st tertile` = "#FFE0DE",   # les catégories faibles expositions (intragroupes)
               `1st tertile - 2nd tertile` = "gray70",  # les variances intergroupes
               `1st tertile - 3rd tertile` = "gray70",
               `2nd tertile - 3rd tertile` = "gray70",
               `2nd tertile` = "#FF8F87",   # les catégorires moyennes expositions (intragroupes)
               `3rd tertile` = "#FF4034"))+
  scale_y_discrete(labels = function(x) {
    lapply(x, function(lbl) {
      if (grepl("1st tertile - 3rd tertile", lbl)) {
        bquote("1"^st* " tertile - 3"^rd* " tertile")
      } else if (grepl("1st", lbl)) {
        bquote("1"^st* " tertile")
      } else if (grepl("3rd", lbl)) {
        bquote("3"^rd* " tertile")
      } else {
        lbl
      }
    })
  })

results_betadiv_detailled %>% 
  filter(Sous_groups == "1st_3rd_ter") %>% 
  filter(!`Explanatory variables` %in% c("Residual", "Total")) %>% 
  View()
my_tag <- c("0.07\n(1.45)", 
            "0.02\n(1.67)", 
            "0.02\n(1.64)", 
            "0.08\n(1.43)", 
            "0.002\n(2.12)", 
            "0.002\n(2.51)") 

boxplot_2 <- tag_facet(boxplot_2, 
                       x = 1.05, y = 2, 
                       vjust = 0.5, hjust = 0.5,
                       open = "", close = "",
                       fontface = 1,
                       size = 3.5,
                       #family = "serif",
                       tag_pool = my_tag) +
  theme_lucid() +
  theme(legend.position = "none", 
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.text.x = element_blank(), 
        axis.line.x = element_blank(),
        strip.background = element_blank(),  # Supprime le fond des facettes
        strip.text = element_blank())      # Supprime le texte des facettes) 

#my_tag_bis <- c("       \n             ", "       \n             ", "       \n             ")
# boxplot_2 <- tag_facet(boxplot_2, 
#                        x = 1.07, y = 3.1, 
#                        vjust = 0.5, hjust = 0.75,
#                        open = "", close = "",
#                        fontface = 1,
#                        size = 3.5,
#                        #family = "serif",
#                        tag_pool = my_tag_bis) +
#   theme_lucid() +
#   theme(legend.position = "none", 
#         axis.title.y = element_blank(), 
#         axis.title.x = element_blank(),
#         axis.text.y = element_text(size = 10, color = "black"),
#         axis.text.x = element_blank(), 
#         axis.line.x = element_blank(),
#         strip.background = element_blank(),  # Supprime le fond des facettes
#         strip.text = element_blank())     # Supprime le texte des facettes) 


#### 2nd vs 3rd tertile ----
boxplot_3 <- bdd_boxplot_article %>%
  filter(Groups %in% c("2nd tertile", "2nd tertile - 3rd tertile", "3rd tertile")) %>%
  ggplot() +
  aes(
    x = Bray_curtis_dissimilarity,
    y = Groups,
    fill = Groups
  ) +
  geom_boxplot() +
  labs(x = "Bray Curtis dissimilarity") +
  theme_lucid() +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 10, color = "black"),
        #axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 10, color = "black"),
        strip.background = element_blank(),  # Supprime le fond des facettes
        strip.text = element_blank()) +     # Supprime le texte des facettes) 
  facet_wrap(~Pollutant_rec, ncol = 6) +
  scale_fill_manual(
    values = c(`1st tertile` = "#FFE0DE",   # les catégories faibles expositions (intragroupes)
               `1st tertile - 2nd tertile` = "gray70",  # les variances intergroupes
               `1st tertile - 3rd tertile` = "gray70",
               `2nd tertile - 3rd tertile` = "gray70",
               `2nd tertile` = "#FF8F87",   # les catégorires moyennes expositions (intragroupes)
               `3rd tertile` = "#FF4034"))+
  scale_y_discrete(labels = function(x) {
    lapply(x, function(lbl) {
      if (grepl("2nd tertile - 3rd tertile", lbl)) {
        bquote("2"^nd* " tertile - 3"^rd* " tertile")
      } else if (grepl("2nd", lbl)) {
        bquote("2"^nd* " tertile")
      } else if (grepl("3rd", lbl)) {
        bquote("3"^rd* " tertile")
      } else {
        lbl
      }
    })
  })

results_betadiv_detailled %>% 
  filter(Sous_groups == "2nd_3rd_ter") %>% 
  filter(!`Explanatory variables` %in% c("Residual", "Total")) %>% 
  View()
my_tag <- c("0.02\n(1.77)", 
            "0.20\n(1.17)", 
            "0.61\n(0.90)", 
            "0.03\n(1.63)", 
            "0.52\n(0.97)", 
            "0.03\n(1.75)") 
boxplot_3 <- tag_facet(boxplot_3, 
                       x = 1.05, y = 2, 
                       vjust = 0.5, hjust = 0.5,
                       open = "", close = "",
                       fontface = 1,
                       size = 3.5,
                       #family = "serif",
                       tag_pool = my_tag) +
  theme_lucid() +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 10, color = "black"),
        #axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 10, color = "black"),
        strip.background = element_blank(),  # Supprime le fond des facettes
        strip.text = element_blank())    # Supprime le texte des facettes) 


# my_tag_bis <- c("       \n             ", "       \n             ", "       \n             ")
# boxplot_3 <- tag_facet(boxplot_3, 
#                        x = 1.07, y = 3.1, 
#                        vjust = 0.5, hjust = 0.75,
#                        open = "", close = "",
#                        fontface = 1,
#                        size = 3.5,
#                        #family = "serif",
#                        tag_pool = my_tag_bis) +
#   theme_lucid() +
#   theme(legend.position = "none", 
#         axis.title.x = element_text(size = 12, face = "bold"),
#         axis.title.y = element_blank(), 
#         axis.text.y = element_text(size = 10, color = "black"),
#         #axis.title.x = element_blank(), 
#         axis.text.x = element_text(size = 10, color = "black"),
#         strip.background = element_blank(),  # Supprime le fond des facettes
#         strip.text = element_blank())    # Supprime le texte des facettes) 


boxplot <- boxplot_1 + boxplot_2 + boxplot_3 + plot_layout(nrow = 3)
boxplot
boxplot_2 + boxplot_3 + plot_layout(ncol = 2)

rm(boxplot_1, boxplot_2, boxplot_3, 
   bdd_long_t2, bdd_long_t3, bdd_long_M2, bdd_long_Y1)

ggsave(plot = boxplot, 
       filename = "4_output/betadiv/boxplot_figure4.tiff", 
       device = "tiff", 
       units = "cm", 
       dpi = 400,
       width = 40, 
       height = 12)


# Test 1 seule fois vérif ----
set.seed(1996)
results <- 
  adonis2(
    all_dist_t2 ~ mo_DEHP_ms_i_cor_t2_ter,
    data = bdd_final_t2,
    permutations = 999)

results <- 
  adonis2(
    all_dist_t2 ~ 
      mo_DEHP_ms_i_cor_t2_ter +
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
      # mo_interpreg_3cat +
      mo_dipl_3cat_i +
      po_w_kg_3cat +
      po_he_3cat_i +
      ch_w_Y1_3cat_i +
      ch_he_Y1_3cat_i +
      po_gd +
      mo_age +
      mo_bmi_bepr_3cat_i +
      bf_duration_till48w_4cat_i,
    data = bdd_final_t2,
    permutations = 999
  )





# Figure S7 ----
heatmap_alpha <- bdd_alpha %>% 
  select(ident, 
         "Specific richness" = ch_feces_SpecRich_5000_ASV_Y1, 
         "Shannon diversity" = ch_feces_Shannon_5000_ASV_Y1, 
         "Faith diversity" = ch_feces_Faith_5000_ASV_Y1)

heatmap_taxa <- bdd_taxa %>%
  select(ident, 
         "Firmicutes" = ch_feces_rel_p1_Y1, 
         "Actinobacteria" = ch_feces_rel_p2_Y1, 
         "Bacteroidetes" = ch_feces_rel_p3_Y1, 
         "Proteobacteria" = ch_feces_rel_p4_Y1)
heatmap_alpha <- left_join(heatmap_alpha, heatmap_taxa, by = "ident") %>%
  select(-ident)

cormat <- round(cor(heatmap_alpha, 
                    use = "pairwise.complete.obs", 
                    method = "spearman"), 3)

get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

reorder_cormat <- function(cormat){              # Utiliser la corrélation entre les variables comme mesure de distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}


cormat <- reorder_cormat(cormat)                 # réordonner les coef de cor
upper_tri <- get_upper_tri(cormat)               # obtenir que le triangle sup
melted_cormat <- melt(upper_tri, na.rm = TRUE)   # passer en df long rapidement 

heatmap <-                                       # heatmap
  ggplot(melted_cormat, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1, 1),
    space = "Lab",
    name = "Spearman\nCorrelation"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(
    angle = 45,
    vjust = 1,
    size = 12,
    hjust = 1
  )) +
  coord_fixed() +
  geom_text(aes(Var2, Var1, label = value),
            color = "black",
            size = 4) +
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


