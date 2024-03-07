# Figures article phthalates
# A.Davias
# 23.02.2024



# Fig.1 ----
# Adjusted associations between concentrations of chemicals and the α-diversity 
# indices of the child gut microbiota at 12-month-old (single-chemical models, 
#n between 344 and 349).

# Fig.2 ----
# Expected changes in α-diversity indices associated with concurrently 
# increasing quantiles of all chemicals (continuous variables), relative to when 
# all concentrations are fixed at their 25th percentile (mixture models, n 
# between 344 and 349).

# Fig.3 ----
# Expected changes in phyla relative abundances associated with concurrently 
# increasing quantiles of all chemicals (continuous variables), relative to when 
# all concentrations are fixed at their 25th percentile (mixture models, n 
# between 344 and 349).

# Fig.4 ----
# Adjusted associations between concentrations of chemicals and the 46 most 
# abundant genera in the child gut microbiota at one year (adjusted uni-pollutant 
# models, n between 350 and 355).


# Table 1 ----
# Adjusted associations between concentrations of chemicals and the α-diversity 
# indices of the child gut microbiota at 12-month-old (single chemical models, n 
# between 344 and 349).

# Table 2 ----
# Adjusted associations between concentrations of chemicals and Bray-Curtis 
# β-diversity in one-year children gut microbiota (n between 344 and 349).

# Table 3 ----
# Adjusted associations between concentrations of chemicals with identified 
# effects on α-diversity indices and the four most abundant phyla of the child gut 
# microbiota at 12-month-old (n = 351).

# Table A.1 ----
# Assessment of exposure to phthalate and DINCH metabolites in SEPAGES cohort.

# Table A.2 ----
# Distribution and comparison of exposures to phthalates and DINCH metabolites 
# in SEPAGES cohort (484  mother-child pairs).

# Table A.3 ----
# Distribution of α-diversity indices and major taxa in 12-month child gut 
# microbiota from SEPAGES  cohort.
load("2_final_data/metadata.RData")
load("2_final_data/bdd_alpha.RData")
load("2_final_data/bdd_taxa.RData")
source("3_programs/4_functions_AD_gumme.R", encoding = 'UTF-8')      # fonctions
source("3_programs/4_vectors_AD_gumme.R", echo=FALSE)
taxa_table <- read_csv("0_source_data/taxa_table_ASVbased_Y1_AD_20220504_8.csv")

bdd_alpha <- bdd_alpha %>%
  select(ident, 
         all_of(alpha_vec))
genera_linear <- bdd_taxa %>%
  select(contains("ch_feces_rel_g")) %>%
  na.omit() %>%
  select_if(~ sum(. != 0, na.rm = TRUE) / length(.) >= 0.3) %>%
  colnames()
bdd_taxa <- bdd_taxa %>%
  select(ident, 
         ch_feces_rel_p1_Y1, ch_feces_rel_p2_Y1, ch_feces_rel_p3_Y1, ch_feces_rel_p4_Y1,
         all_of(genera_linear))

bdd <- left_join(bdd_alpha, bdd_taxa, by = "ident")
vars <- colnames(bdd[,-1])

Table_A_3 <- descrip_num(data = bdd, vars = vars)
Table_A_3 <- Table_A_3 %>%
  rename("Min." = Min,
         "1st quartile" = Q1, 
         "3rd quartile" = Q3, 
         "Max." = Max) %>%
  mutate(
    `Variable labels` = str_replace_all(`Variable labels`, 
                                        "One year child feces relative abundance of genus ", ""),
    `Variable labels` = str_replace_all(`Variable labels`, 
                                        "One year child feces relative abundance of phylum ", ""))
genera_linear <- unique(Table_A_3$`Variable labels`) %>% as.character() 
genera_linear <- genera_linear[8:53]
corres <- 
  taxa_table %>% 
  select(Phyla_corres = ch_feces_phylum_ASVbased_Y1, 
         `Variable labels` = ch_feces_genus_ASVbased_Y1) %>%
  filter(`Variable labels` %in% genera_linear) %>%
  distinct(`Variable labels`, .keep_all = TRUE)

Table_A_3 <- left_join(Table_A_3, corres, by = "Variable labels")
Table_A_3 <- Table_A_3 %>%
  mutate(
    `Variable labels` = fct_recode(`Variable labels`,
                                   "Faith phylogenetic diversity" = "Faith diversity in 1Y child gut microbiota ASV based (seq. depth = 5 000",
                                   "Shannon diversity" = "Shannon diversity in 1Y child gut microbiota ASV based (seq. depth = 5 000",
                                   "Specific richness" = "Specific richness diversity in 1Y child gut microbiota ASV based (seq. depth = 5 000",
                                   "Clostridium IV" = "Clostridium_IV",
                                   "Clostridium sensu stricto" = "Clostridium_sensu_stricto",
                                   "Clostridium XlVa" = "Clostridium_XlVa",
                                   "Clostridium XVIII" = "Clostridium_XVIII",
                                   "Erysipelotrichaceae incertae sedis" = "Erysipelotrichaceae_incertae_sedis",
                                   "Escherichia Shigella" = "Escherichia_Shigella",
                                   "Lachnospiracea incertae sedis" = "Lachnospiracea_incertae_sedis",
                                   "Saccharibacteria genera incertae sedis" = "Saccharibacteria_genera_incertae_sedis"), 
    Phyla_corres = fct_recode(Phyla_corres, "Candidatus Saccharibacteria" = "Candidatus_Saccharibacteria"),
    Phyla_corres = fct_relevel(Phyla_corres, 
                               "Firmicutes", "Actinobacteria", "Bacteroidetes", "Proteobacteria", "Verrucomicrobia", "Candidatus Saccharibacteria"))%>%
  arrange(Phyla_corres, desc(Median))
write.xlsx(Table_A_3, file = "4_output/Table_A_3.xlsx")
rm(list = ls())


# Table A.4 ----
# Adjusted associations between concentrations of chemicals during the second 
# trimester of pregnancy and the α-diversity of the child gut microbiota at 12 
# months (mixture effects analysis, n = 349).

# Table A.5 ----
# Adjusted associations between concentrations of chemicals during the third 
# trimester of pregnancy and the α-diversity of the child gut microbiota at 12 
# months (mixture effects analysis, n=344).

# Table A.6 ----
# Adjusted associations between concentrations of chemicals at twelve months and 
# the α-diversity of the child gut microbiota at 12 months (mixture effects 
# analysis, n=345).

# Table A.7 ----
# Adjusted associations between concentrations of chemicals without identified 
# effects on α-diversity and the most abundant phyla in 12-month children gut 
# microbiota (n between 350 and 355).

# Table A.8 ----
# Adjusted associations between concentrations of chemicals and the 46 most 
# abundant genera in one-year children gut microbiota (n between 340 and 355).

# Table A.9 ----
# Adjusted associations between concentrations of chemicals and the gut 
# microbiota α-diversity at different sequencing depths.

# Table A.10 ----
# Adjusted associations between concentrations of chemicals and the gut 
# microbiota composition before and after adjustment for specific gravity.

# Figure A.1 ----
# Directed acyclic graph of the relation between concentrations of chemicals and 
# child gut microbiota at 12 months of age.

# Figure A.2 ----
# Pearson correlations between concentrations of chemicals assessed during 
# various exposure windows.  
load("2_final_data/metadata.RData")
load("2_final_data/bdd_alpha.RData")
load("2_final_data/bdd_taxa.RData")
source("3_programs/4_functions_AD_gumme.R", encoding = 'UTF-8')      # fonctions
source("3_programs/4_vectors_AD_gumme.R", echo=FALSE)
bdd_phthalates <- bdd_taxa %>% 
  select(all_of(phthalates_vec), 
         ch_feces_rel_p1_Y1) %>%
  filter(!is.na(ch_feces_rel_p1_Y1)) %>%
  select(-ch_feces_rel_p1_Y1) %>%
  select(!contains(c("MEOHP", "MECPP", "MEHHP", "MEHP", "MMCHP",     
                     "ohMiNP", "oxoMiNP", "cxMiNP",                 
                     "ohMINCH", "oxoMINCH"))) %>%  
  select(!contains("cat")) %>%
  select(!contains("M2"))

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

bdd_phthalates <- bdd_phthalates %>% na.omit()
for (nom in names(bdd_phthalates)) {                   # enleve les étiquettes des variables 
  var_lab(bdd_phthalates[[nom]]) <- NULL
}

bdd_phthalates <- bdd_phthalates %>%
  select("ΣDEHP t2","ΣDEHP t3", "ΣDEHP Y1",
         "MnBP t2", "MnBP t3", "MnBP Y1",
         "ΣDiNP t2", "ΣDiNP t3",  "ΣDiNP Y1",
         "MiBP t2", "MiBP t3", "MiBP Y1",
         "MBzP t2", "MBzP t3", "MBzP Y1",
         "MEP t2", "MEP t3", "MEP Y1",
         "ohMPHP t2", "ohMPHP t3", "ohMPHP Y1",
         "ΣDINCH t2", "ΣDINCH t3", "ΣDINCH Y1")
  

library(corrplot)

heatmap_phthalates <- cor(bdd_phthalates, 
                          use = "pairwise.complete.obs", 
                          method = "pearson")

tiff(filename = "4_output/heatmap_cor_phthalates.tiff", units = "mm", width = 250, height = 250, res = 300)
corrplot(heatmap_phthalates, 
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

# Figure A.3 ----
# Spearman correlations between concentrations of chemicals and covariates.
load("2_final_data/metadata.RData")
load("2_final_data/bdd_alpha.RData")
load("2_final_data/bdd_taxa.RData")
source("3_programs/4_functions_AD_gumme.R", encoding = 'UTF-8')      # fonctions
source("3_programs/4_vectors_AD_gumme.R", echo=FALSE)
vars_1 <- metadata %>% 
  select(all_of(phthalates_vec)) %>%
  select(!contains(c("MEOHP", "MECPP", "MEHHP", "MEHP", "MMCHP",     
                     "ohMiNP", "oxoMiNP", "cxMiNP",                 
                     "ohMINCH", "oxoMINCH"))) %>%  
  select(!contains("cat")) %>%
  select(!contains("M2")) %>%
  colnames()
bdd_taxa <- bdd_taxa %>% select(ident, ch_feces_rel_p1_Y1)
cormat <- metadata %>% 
  select(ident, 
         all_of(vars_1),
         all_of(covar_vec_num_i)) 
cormat <- left_join(bdd_taxa, cormat, by ="ident") 
cormat <- cormat %>%
  filter(!is.na(ch_feces_rel_p1_Y1)) %>%
  select(-ch_feces_rel_p1_Y1, 
         -ident) %>%
  na.omit()

cormat <- round(cor(cormat, 
                    use = "pairwise.complete.obs", 
                    method = "spearman"), 1)

cormat <- cormat %>% 
  as.data.frame() %>% 
  select(all_of(covar_vec_num_i)) %>% 
  t() %>%
  as.data.frame() %>%
  select(all_of(vars_1)) %>%
  as.matrix()

rownames(cormat) <- rownames(cormat) %>%
  str_replace_all(
    c("ch_feces_age_w_Y1_i" = "Child age (weeks)",
      "ch_antibio_Y1_i" = "Antibiotics use 0-12 months",
      "mo_par" = "Maternal parity",
      "po_w_kg" = "Birth weight (kg)",
      "po_he_i"= "Birth length (cm)", 
      "ch_w_Y1_i"="Weight at one year (Kg)", 
      "ch_he_Y1_i"="Length at one year (cm)", 
      "po_gd"= "Gestational age (weeks)", 
      "mo_age"="Maternal age before pregnancy", 
      "mo_bmi_bepr_i"="Maternal BMI before pregnancy", 
      "bf_duration_till48w_i"="Breastfeeding duration (weeks)"))

colnames(cormat) <- colnames(cormat) %>%
  str_replace_all(
    c("mo_" = "",
      "ch_" = "",
      "_total_i_cor_" = " ",
      "_i_cor_" = " ", 
      "_ms" = "", 
      "DINCH" = "ΣDINCH",
      "DiNP" = "ΣDiNP",
      "DEHP" = "ΣDEHP", 
      "_ln" =""))
cormat <- cormat %>%
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

tiff(filename = "4_output/heatmap_cor_phthalates_covar.tiff", units = "mm", width = 300, height = 150, res = 300)
corrplot(cormat, 
         method = 'color', 
         tl.col = 'black', 
         tl.srt = 45, 
         addCoef.col = "black",
         number.cex = 0.8,
         tl.cex = 0.8,
         number.digits = 1,
         col = rev(COL2(diverging = "RdYlBu")))
dev.off()
rm(list = ls())

# Figure A.4 ----
# Estimated effect of an increase  from the 25th to 75th percentile in a single 
# chemical concentration on α-diversity indices when all other chemicals are 
# fixed at either the 25th, 50th, or 75th percentiles (mixture models, n between 
# 344 and 349).

# Figure A.5 ----
# Pearson correlations between the α-diversity  indices and the most abundant 
# taxa in the gut microbiota of SEPAGES children.
load("2_final_data/metadata.RData")
load("2_final_data/bdd_alpha.RData")
load("2_final_data/bdd_taxa.RData")
source("3_programs/4_functions_AD_gumme.R", encoding = 'UTF-8')      # fonctions
source("3_programs/4_vectors_AD_gumme.R", echo=FALSE)
library(corrplot)
bdd_alpha <- bdd_alpha %>% select(ident, all_of(alpha_vec))
bdd_taxa <- bdd_taxa %>% select(ident, ch_feces_rel_p1_Y1,  ch_feces_rel_p2_Y1, ch_feces_rel_p3_Y1, ch_feces_rel_p4_Y1)
bdd <- left_join(bdd_alpha, bdd_taxa, by = "ident")
bdd <- bdd %>%
  select(-ident) %>%
  rename("Specific richness" = "ch_feces_SpecRich_5000_ASV_Y1",
         "Shannon diversity" = "ch_feces_Shannon_5000_ASV_Y1", 
         "Faith phylogenetic diversity" = "ch_feces_Faith_5000_ASV_Y1",
         "Firmicutes" = "ch_feces_rel_p1_Y1", 
         "Actinobacteria" = "ch_feces_rel_p2_Y1", 
         "Bacteroidetes" = "ch_feces_rel_p3_Y1", 
         "Proteobacteria" = "ch_feces_rel_p4_Y1")


heatmap_phthalates <- cor(bdd, 
                          use = "pairwise.complete.obs", 
                          method = "pearson")

tiff(filename = "4_output/heatmap_cor_phthalates.tiff", units = "mm", width = 150, height = 150, res = 300)
corrplot(heatmap_phthalates, 
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
rm(list = ls())
