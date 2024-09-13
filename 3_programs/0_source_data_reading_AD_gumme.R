# 0_source_data_reading
# Aline Davias 
# 2024/08/27


# 0_source_data_reading ----
## 1) Chargement des packages ----
library(haven)
library(readr)
library(tidyverse)
library(foreign)
library(questionr)

## 2) Chargement des métadonnées ----
metadata <- read_sas(
  "0_source_data/base_aline_211115.sas7bdat",                                   # base de données SEPAGES
  catalog_file = "0_source_data/formats.sas7bcat")

bdd_sg <- read_sas(
  "0_source_data/base_aline_220909.sas7bdat",                                   # base de données SEPAGES plus récente avec les gravités spécifiques
  catalog_file = "0_source_data/formats.sas7bcat") %>%
  select(ident, mo_pool_sg_T1, mo_pool_sg_T3, ch_pool_sg_Y1)

bdd_crb <- read_sas(
  "0_source_data/date_selle_aline210520.sas7bdat",                              # base avec les dates de prélévement selles IAB (recu séparément)
  NULL)         

alphadiv_Y1 <- read_labelled_csv(
  "0_source_data/alpha_diversity_ASVbased_Y1_labelled_AD_20220504_32.csv")

asv_taxa <- read_labelled_csv(
  "0_source_data/gut_microbiota_ASVbased_Y1_labelled_AD_20220504_7239.csv") %>%
  select(ident, 
         ch_feces_ID_Y1, 
         starts_with("ch_feces_rel")) %>%
  select(!starts_with("ch_feces_rel_ASV"))

metadata_microbiote <-
  read_labelled_csv("0_source_data/gut_microbiota_ASVbased_Y1_labelled_AD_20220504_7239.csv") %>%
  select(!starts_with("ch_feces_rel")) %>%
  select(!starts_with("ch_feces_raw"))

age_feces_Y1 <- read_labelled_csv("0_source_data/age_feces_collection_Y1_labelled_AD_20220504_1.csv")
atb_Y1 <- read_labelled_csv("0_source_data/antibiotics_Y1_labelled_AD_20220314_9.csv")
solidfood_Y1 <- read_labelled_csv("0_source_data/solidfood_Y1_labelled_AD_20220302_7.csv")
weight_length_Y1 <- read_labelled_csv("0_source_data/weight_length_Y1_labelled_AD_20220301_2.csv")

## 3) Fusion des bases de données ----
### data child (n=484) + data bdd_crb (n=360)
metadata <- metadata %>% rename(ch_feces_ID_Y1 = CodeEchantillon_selle_un_an)
bdd_crb <- bdd_crb %>% rename(ch_feces_ID_Y1 = Code_Echantillon)

var_lab(metadata$ident) <- NULL
var_lab(metadata$ch_feces_ID_Y1) <- NULL
var_lab(alphadiv_Y1$ident) <- NULL
var_lab(alphadiv_Y1$ch_feces_ID_Y1) <- NULL
var_lab(asv_taxa$ident) <- NULL
var_lab(asv_taxa$ch_feces_ID_Y1) <- NULL
var_lab(age_feces_Y1$ident) <- NULL
var_lab(atb_Y1$ident) <- NULL
var_lab(solidfood_Y1$ident) <- NULL
var_lab(weight_length_Y1$ident) <- NULL
var_lab(metadata_microbiote$ident) <- NULL
var_lab(metadata_microbiote$ch_feces_ID_Y1) <- NULL
var_lab(bdd_sg$ident)  <- NULL

bdd <- 
  left_join(metadata,
            bdd_crb,
            by = c("ident", "ch_feces_ID_Y1")) %>%
  mutate(ch_feces_ID_Y1 = 
           str_replace_all(ch_feces_ID_Y1,
                           c("SELLE000089.1" = "SELLE0089",
                             "SELLE0711.1" = "SELLE0711",
                             "SELLE1066.1" = "SELLE1066")),
         ch_feces_ID_Y1 = str_sub(ch_feces_ID_Y1, 1, 9), 
         ident = as.integer(ident))

bdd <- left_join(bdd, asv_taxa, by = c("ident", "ch_feces_ID_Y1"))
bdd <- left_join(bdd, alphadiv_Y1, by = c("ident", "ch_feces_ID_Y1"))
bdd <- left_join(bdd, metadata_microbiote, by = c("ident", "ch_feces_ID_Y1"))
bdd <- left_join(bdd, age_feces_Y1, by = "ident")
bdd <- left_join(bdd, atb_Y1, by = "ident")
bdd <- left_join(bdd, solidfood_Y1, by = "ident")
bdd <- left_join(bdd, weight_length_Y1, by = "ident")
bdd <- left_join(bdd, bdd_sg, by = "ident")

rm(
  bdd_crb,
  metadata,
  asv_taxa,
  alphadiv_Y1,
  metadata_microbiote,
  age_feces_Y1,
  atb_Y1,
  solidfood_Y1,
  weight_length_Y1, 
  bdd_sg)

save.image("1_intermediate_data/0_source_data_reading_AD_gumme.RData")
