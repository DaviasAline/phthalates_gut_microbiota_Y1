### Aline Davias
### 2022/05/04
### Cleaning of gut microbiota SEPAGES data 1 year
### (excel file with 3 sheets, received from Dr. Lepage INRAE)

library(readxl)
library(haven)
library(tidyverse)
library(gtsummary)
library(magrittr)
library(labelled)
library(expss)

# 0_source_data_reading ----
## Data renaming ----
file.rename(from = "0_source_data/DB_16S_ASV_Gumme1Y_R2R3.xlsx",       
            to = "0_source_data/DB_16S_ASV_Gumme_Y1.xlsx")

## Correspondence key ----
ident_key <-
  read_sas(
    "0_source_data/base_aline_211115.sas7bdat",                # data file
    catalog_file = "0_source_data/formats.sas7bcat") %>%       # label file
  select(ident, CodeEchantillon_selle_un_an)                   # select key variables

## Raw ASV table ----
### (not filtered)
asv_raw <- read_excel(
  "0_source_data/DB_16S_ASV_Gumme_Y1.xlsx", 
  sheet = "RAW_ASV_unfiltered")

## Relative ASV table ----
### (filtered at 0.7% at the family taxonomic level)
asv_rel <- read_excel(
  "0_source_data/DB_16S_ASV_Gumme_Y1.xlsx",
  sheet = "ASV_TAB&TAX", 
  range = cell_cols("I:NA"))

## ASV taxonomic correspondence ----
### (do not correspond to a classic dataframe with variables in columns and samples in rows)
taxa_ASVbased <- read_excel(                                                             
  "0_source_data/DB_16S_ASV_Gumme_Y1.xlsx", 
  sheet = "ASV_TAB&TAX", 
  range = cell_cols("A:H"))

## Sequencing metadata (ASV based)----
metadata_ASVbased <- read_excel(
  "0_source_data/DB_16S_ASV_Gumme_Y1.xlsx", 
  sheet = "metadata")


# 1_data_cleaning ----
## Correspondance key ----
ident_key <- ident_key %>%                                                               
  rename(
    ch_feces_ID_Y1 = CodeEchantillon_selle_un_an) %>%
  mutate(
    ch_feces_ID_Y1 = 
      str_replace_all(ch_feces_ID_Y1, 
                      c("SELLE000089.1" = "SELLE0089", 
                        "SELLE0711.1" = "SELLE0711", 
                        "SELLE1066.1" = "SELLE1066")), 
    ch_feces_ID_Y1 = str_sub(ch_feces_ID_Y1, 1, 9))

## Raw ASV table ----
asv_raw <- asv_raw %>%                                          
  mutate(
    ASV = str_replace(ASV, "OTU", "ASV")) %>% 
  filter(!is.na(ASV))%>%
  column_to_rownames(var="ASV")%>%
  select(-starts_with(c("Failed", "GTSELLE")))%>%
  t() %>%
  as.data.frame() %>% 
  rownames_to_column(var = "ch_feces_ID_Y1") %>%
  mutate(
    ch_feces_ID_Y1 = str_sub(ch_feces_ID_Y1, 1, 9))

var_label(asv_raw) <- c(                                       # set correct variable labels                    
  "Child feces identity at one year", 
  paste("One year child feces raw abundance of ", 
        colnames(asv_raw[,2:3387]), 
        sep = ""))

colnames(asv_raw) <- c(                                        # set correct variable names                         
  "ch_feces_ID_Y1",  
  paste("ch_feces_raw_", 
        colnames(asv_raw[,2:3387]), 
        "_Y1", 
        sep = ""))

## Relative ASV table ----
asv_rel <- asv_rel %>%                         
  select(ASV_ID, everything()) %>% 
  mutate(
    ASV_ID = str_replace(ASV_ID, "OTU", "ASV")) %>%
  column_to_rownames(var = "ASV_ID") %>%
  t() %>%
  as.data.frame() %>% 
  rownames_to_column(var = "ch_feces_ID_Y1") 

var_label(asv_rel) <- c(                                       # set correct variable labels                    
  "Child feces identity at one year", 
  paste("One year child feces relative abundance of ", 
        colnames(asv_rel[,2:3238]), 
        sep = ""))

colnames(asv_rel) <- c(                                        # set correct variable names            
  "ch_feces_ID_Y1", 
  paste("ch_feces_rel_", 
        colnames(asv_rel[,2:3238]), 
        "_Y1", 
        sep = ""))

## Taxonomic correspondance (ASV based) ----
colnames(taxa_ASVbased) <- c(
  "ch_feces_ASV_ID_Y1",
  paste("ch_feces_", 
        colnames(taxa_ASVbased[, 2:8]), 
        "_ASVbased_Y1", 
        sep = ""))
taxa_ASVbased <- taxa_ASVbased %>%
  mutate(
    across(all_of(
      c("ch_feces_domain_ASVbased_Y1",
        "ch_feces_phylum_ASVbased_Y1",
        "ch_feces_class_ASVbased_Y1",
        "ch_feces_order_ASVbased_Y1",
        "ch_feces_family_ASVbased_Y1",
        "ch_feces_genus_ASVbased_Y1"
      )
    ), factor),
    ch_feces_ASV_ID_Y1 = str_replace(ch_feces_ASV_ID_Y1, "OTU", "ASV"))

## Sequencing metadata (ASV based) ----
metadata_ASVbased <- metadata_ASVbased %>%
  rename(                                                      # set correct variable names   
    ch_feces_ID_Y1 = SAMPLE_ID, 
    ch_feces_OrderSeq_ASVbased_Y1 = OrderSeq, 
    ch_feces_Nreads_ASVbased_Y1 = Nreads, 
    ch_feces_NASV_Y1 = NASV, 
    ch_feces_RUN_Y1 = RUN, 
    ch_feces_ForAnalyses_Y1 = ForAnalyses, 
    ch_feces_SEQ_Y1 = SEQ) %>%
  mutate(
    ch_feces_RUN_Y1 = fct_recode(ch_feces_RUN_Y1,
                                 "Run1" = "R2",
                                 "Run2" = "R3"))

var_label(metadata_ASVbased) <- list(                          # set correct variable labels
  ch_feces_ID_Y1 = "Child feces identity at one year",
  ch_feces_OrderSeq_ASVbased_Y1 = "Sequencing order (one year child feces, ASV based)",
  ch_feces_Nreads_ASVbased_Y1 = "Sequencing depth (one year child feces, ASV based)",
  ch_feces_NASV_Y1 = "Specific richness before rarefaction (one year child feces, ASV based)",
  ch_feces_RUN_Y1 = "Run number (one year child feces)",
  ch_feces_ForAnalyses_Y1 = "1:sequencing quality validated for analyses (one year child feces)", 
  ch_feces_SEQ_Y1 = "Sequencing type (One year child feces)")

## Taxa raw abundances (ASV based) ----
taxa_raw_ASVbased <- asv_raw %>%                               # data preparation
  column_to_rownames(var="ch_feces_ID_Y1") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "ch_feces_ASV_ID_Y1")%>%
  mutate(
    ch_feces_ASV_ID_Y1 = str_replace_all(ch_feces_ASV_ID_Y1, c("ch_feces_raw_" = "", "_Y1"=""))) %>%
  left_join(taxa_ASVbased, by = "ch_feces_ASV_ID_Y1") %>% 
  select(ch_feces_ASV_ID_Y1, 
         ch_feces_TAX_ASVbased_Y1,
         ch_feces_domain_ASVbased_Y1, 
         ch_feces_phylum_ASVbased_Y1, 
         ch_feces_class_ASVbased_Y1, 
         ch_feces_order_ASVbased_Y1, 
         ch_feces_family_ASVbased_Y1, 
         ch_feces_genus_ASVbased_Y1,
         everything()) %>%
  na.omit() %>%
  mutate(
    across(all_of(
      c("ch_feces_phylum_ASVbased_Y1",
        "ch_feces_class_ASVbased_Y1",
        "ch_feces_order_ASVbased_Y1",
        "ch_feces_family_ASVbased_Y1",
        "ch_feces_genus_ASVbased_Y1")
    ), factor))

### genus level ----
genus_raw <- taxa_raw_ASVbased %>%                             # select the variables 
  select(ch_feces_ASV_ID_Y1, 
         ch_feces_genus_ASVbased_Y1, 
         starts_with("SELLE"))

genus_raw <-                                                   # grouping of ASVs belonging to the same genus
  apply(
    genus_raw[, 3:358], 2, tapply, genus_raw$ch_feces_genus_ASVbased_Y1, sum) %>%
  as.data.frame()

genus_raw <- genus_raw %>%                                     # order the dataframe depending on genus average abundances 
  arrange(desc(rowMeans(genus_raw))) %>%
  t() %>%
  as.data.frame()%>%
  rownames_to_column("ch_feces_ID_Y1")

var_label(genus_raw) <- c(                                     # prepare new variables labels and names 
  "Child feces identity at one year", 
  paste("One year child feces raw abundance of genus", colnames(genus_raw[, 2:193]), sep = " "))

colnames(genus_raw) <- c(
  "ch_feces_ID_Y1", 
  paste("ch_feces_raw_g", 1:192, "_Y1", sep = ""))    

### family level ----
family_raw <- taxa_raw_ASVbased %>%                            # select the variables 
  select(ch_feces_ASV_ID_Y1, 
         ch_feces_family_ASVbased_Y1, 
         starts_with("SELLE")) 

family_raw <-                                                  # grouping of ASVs belonging to the same family
  apply(
    family_raw[, 3:358], 2, tapply, family_raw$ch_feces_family_ASVbased_Y1, sum) %>%
  as.data.frame()

family_raw <- family_raw %>%                                   # order the dataframe depending on family average abundances 
  arrange(desc(rowMeans(family_raw))) %>%
  t() %>%
  as.data.frame()%>%
  rownames_to_column("ch_feces_ID_Y1")

var_label(family_raw) <- c(                                    # prepare new variables labels and names 
  "Child feces identity at one year", 
  paste("One year child feces raw abundance of family", colnames(family_raw[, 2:60]), sep = " "))

colnames(family_raw) <- c(
  "ch_feces_ID_Y1", 
  paste("ch_feces_raw_f", 1:59, "_Y1", sep = ""))  

### order level ----
order_raw <- taxa_raw_ASVbased %>%                             # select the variables 
  select(ch_feces_ASV_ID_Y1, 
         ch_feces_order_ASVbased_Y1, 
         starts_with("SELLE")) 

order_raw <-                                                   # grouping of ASVs belonging to the same order
  apply(
    order_raw[, 3:358], 2, tapply, order_raw$ch_feces_order_ASVbased_Y1, sum) %>%
  as.data.frame()

order_raw <- order_raw %>%                                     # order the dataframe depending on order average abundances 
  arrange(desc(rowMeans(order_raw))) %>%
  t() %>%
  as.data.frame()%>%
  rownames_to_column("ch_feces_ID_Y1")

var_label(order_raw) <- c(                                     # prepare new variables labels and names 
  "Child feces identity at one year", 
  paste("One year child feces raw abundance of order", colnames(order_raw[, 2:27]), sep = " "))

colnames(order_raw) <- c(
  "ch_feces_ID_Y1", 
  paste("ch_feces_raw_o", 1:26, "_Y1", sep = ""))  

### class level ----
class_raw <- taxa_raw_ASVbased %>%                             # select the variables 
  select(ch_feces_ASV_ID_Y1, 
         ch_feces_class_ASVbased_Y1, 
         starts_with("SELLE"))  

class_raw <-                                                   # grouping of ASVs belonging to the same class
  apply(
    class_raw[, 3:358], 2, tapply, class_raw$ch_feces_class_ASVbased_Y1, sum) %>%
  as.data.frame()

class_raw <- class_raw %>%                                     # class the dataframe depending on class average abundances 
  arrange(desc(rowMeans(class_raw))) %>%
  t() %>%
  as.data.frame()%>%
  rownames_to_column("ch_feces_ID_Y1")

var_label(class_raw) <- c(                                     # prepare new variables labels and names 
  "Child feces identity at one year", 
  paste("One year child feces raw abundance of class", colnames(class_raw[, 2:18]), sep = " "))

colnames(class_raw) <- c(
  "ch_feces_ID_Y1", 
  paste("ch_feces_raw_c", 1:17, "_Y1", sep = ""))  

### phylum level ----
phylum_raw <- taxa_raw_ASVbased %>%                            # select the variables corresponding to phyla
  select(ch_feces_ASV_ID_Y1, 
         ch_feces_phylum_ASVbased_Y1, 
         starts_with("SELLE")) 

phylum_raw <-                                                  # grouping of ASVs belonging to the same phylum
  apply(
    phylum_raw[, 3:358], 2, tapply, phylum_raw$ch_feces_phylum_ASVbased_Y1, sum) %>%
  as.data.frame()

phylum_raw <- phylum_raw %>%                                   # phylum the dataframe depending on phylum average abundances 
  arrange(desc(rowMeans(phylum_raw))) %>%
  t() %>%
  as.data.frame()%>%
  rownames_to_column("ch_feces_ID_Y1")

var_label(phylum_raw) <- c(                                    # prepare new variables labels and names 
  "Child feces identity at one year", 
  paste("One year child feces raw abundance of phylum", colnames(phylum_raw[, 2:11]), sep = " "))

colnames(phylum_raw) <- c(
  "ch_feces_ID_Y1", 
  paste("ch_feces_raw_p", 1:10, "_Y1", sep = ""))  

## Taxa relative abundances (ASV based) ----
taxa_rel_ASVbased <- asv_rel %>%                               # data preparation
  column_to_rownames(var="ch_feces_ID_Y1") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "ch_feces_ASV_ID_Y1")%>%
  mutate(
    ch_feces_ASV_ID_Y1 = str_replace_all(ch_feces_ASV_ID_Y1, c("ch_feces_rel_" = "", "_Y1"=""))) %>%
  left_join(taxa_ASVbased, by = "ch_feces_ASV_ID_Y1") %>% 
  select(ch_feces_ASV_ID_Y1, 
         ch_feces_TAX_ASVbased_Y1,
         ch_feces_domain_ASVbased_Y1, 
         ch_feces_phylum_ASVbased_Y1, 
         ch_feces_class_ASVbased_Y1, 
         ch_feces_order_ASVbased_Y1, 
         ch_feces_family_ASVbased_Y1, 
         ch_feces_genus_ASVbased_Y1,
         everything()) %>%
  na.omit() %>%
  mutate(
    across(all_of(
      c("ch_feces_phylum_ASVbased_Y1",
        "ch_feces_class_ASVbased_Y1",
        "ch_feces_order_ASVbased_Y1",
        "ch_feces_family_ASVbased_Y1",
        "ch_feces_genus_ASVbased_Y1")
    ), factor))

### genus level ----
genus_rel <- taxa_rel_ASVbased %>%                             # select the variables
  select(ch_feces_ASV_ID_Y1, 
         ch_feces_genus_ASVbased_Y1, 
         starts_with("SELLE"))

genus_rel <-                                                   # grouping of ASVs belonging to the same genus
  apply(
    genus_rel[, 3:358], 2, tapply, genus_rel$ch_feces_genus_ASVbased_Y1, sum) %>%
  as.data.frame()

genus_rel <- genus_rel %>%                                     # order the dataframe depending on genus average abundances 
  arrange(desc(rowMeans(genus_rel))) %>%
  t() %>%
  as.data.frame()%>%
  rownames_to_column("ch_feces_ID_Y1")

var_label(genus_rel) <- c(                                     # prepare new variables labels and names 
  "Child feces identity at one year", 
  paste("One year child feces relative abundance of genus", colnames(genus_rel[, 2:193]), sep = " "))

colnames(genus_rel) <- c(
  "ch_feces_ID_Y1", 
  paste("ch_feces_rel_g", 1:192, "_Y1", sep = ""))    

### family level ----
family_rel <- taxa_rel_ASVbased %>%                            # select the variables corresponding to families
  select(ch_feces_ASV_ID_Y1, 
         ch_feces_family_ASVbased_Y1, 
         starts_with("SELLE"))

family_rel <-                                                  # grouping of ASVs belonging to the same family
  apply(
    family_rel[, 3:358], 2, tapply, family_rel$ch_feces_family_ASVbased_Y1, sum) %>%
  as.data.frame()

family_rel <- family_rel %>%                                   # order the dataframe depending on family average abundances 
  arrange(desc(rowMeans(family_rel))) %>%
  t() %>%
  as.data.frame()%>%
  rownames_to_column("ch_feces_ID_Y1")

var_label(family_rel) <- c(                                    # prepare new variables labels and names 
  "Child feces identity at one year", 
  paste("One year child feces relative abundance of family", colnames(family_rel[, 2:60]), sep = " "))

colnames(family_rel) <- c(
  "ch_feces_ID_Y1", 
  paste("ch_feces_rel_f", 1:59, "_Y1", sep = ""))  

### order level ----
order_rel <- taxa_rel_ASVbased %>%                             # select the variables corresponding to orders
  select(ch_feces_ASV_ID_Y1, 
         ch_feces_order_ASVbased_Y1, 
         starts_with("SELLE")) 

order_rel <-                                                   # grouping of ASVs belonging to the same order
  apply(
    order_rel[, 3:358], 2, tapply,order_rel$ch_feces_order_ASVbased_Y1, sum) %>%
  as.data.frame()

order_rel <- order_rel %>%                                     # order the dataframe depending on order average abundances 
  arrange(desc(rowMeans(order_rel))) %>%
  t() %>%
  as.data.frame()%>%
  rownames_to_column("ch_feces_ID_Y1")

var_label(order_rel) <- c(                                     # prepare new variables labels and names 
  "Child feces identity at one year", 
  paste("One year child feces relative abundance of order", colnames(order_rel[, 2:27]), sep = " "))

colnames(order_rel) <- c(
  "ch_feces_ID_Y1", 
  paste("ch_feces_rel_o", 1:26, "_Y1", sep = ""))  

### class level ----
class_rel <- taxa_rel_ASVbased %>%                             # select the variables corresponding to classes
  select(ch_feces_ASV_ID_Y1, 
         ch_feces_class_ASVbased_Y1, 
         starts_with("SELLE"))

class_rel <-                                                   # grouping of ASVs belonging to the same class
  apply(
    class_rel[, 3:358], 2, tapply, class_rel$ch_feces_class_ASVbased_Y1, sum) %>%
  as.data.frame()

class_rel <- class_rel %>%                                     # class the dataframe depending on class average abundances 
  arrange(desc(rowMeans(class_rel))) %>%
  t() %>%
  as.data.frame()%>%
  rownames_to_column("ch_feces_ID_Y1")

var_label(class_rel) <- c(                                     # prepare new variables labels and names 
  "Child feces identity at one year", 
  paste("One year child feces relative abundance of class", colnames(class_rel[, 2:18]), sep = " "))

colnames(class_rel) <- c(
  "ch_feces_ID_Y1", 
  paste("ch_feces_rel_c", 1:17, "_Y1", sep = ""))  

### phylum level ----
phylum_rel <- taxa_rel_ASVbased %>%                            # select the variables corresponding to phyla
  select(ch_feces_ASV_ID_Y1, 
         ch_feces_phylum_ASVbased_Y1, 
         starts_with("SELLE")) 

phylum_rel <-                                                  # grouping of ASVs belonging to the same phylum
  apply(
    phylum_rel[, 3:358], 2, tapply, phylum_rel$ch_feces_phylum_ASVbased_Y1, sum) %>%
  as.data.frame()

phylum_rel <- phylum_rel %>%                                   # phylum the dataframe depending on phylum average abundances 
  arrange(desc(rowMeans(phylum_rel))) %>%
  t() %>%
  as.data.frame()%>%
  rownames_to_column("ch_feces_ID_Y1")

var_label(phylum_rel) <- c(                                    # prepare new variables labels and names 
  "Child feces identity at one year", 
  paste("One year child feces relative abundance of phylum", colnames(phylum_rel[, 2:11]), sep = " "))

colnames(phylum_rel) <- c(
  "ch_feces_ID_Y1", 
  paste("ch_feces_rel_p", 1:10, "_Y1", sep = ""))  


# 2_variables_selection ----
## Merge final data ----
final_data <- list(
    ident_key, 
    metadata_ASVbased, 
    asv_rel, 
    asv_raw,
    phylum_raw, 
    class_raw, 
    order_raw, 
    family_raw, 
    genus_raw, 
    phylum_rel, 
    class_rel, 
    order_rel, 
    family_rel, 
    genus_rel
) %>%
  reduce(left_join, by = "ch_feces_ID_Y1")

var_label(final_data$ch_feces_ID_Y1) <- "Child feces identity at one year"

## Export ----
saveRDS(final_data, 
        "1_intermediate_data/gut_microbiota_ASVbased_Y1_AD_20220504_7239.rds")
write_dta(final_data, 
          "1_intermediate_data/gut_microbiota_ASVbased_Y1_AD_20220504_7239.dta")
write.csv(final_data, 
          "1_intermediate_data/gut_microbiota_ASVbased_Y1_AD_20220504_7239.csv", 
          row.names = FALSE)
write_labelled_csv(final_data, 
                   "1_intermediate_data/gut_microbiota_ASVbased_Y1_labelled_AD_20220504_7239.csv", 
                   row.names = FALSE)  

saveRDS(taxa_ASVbased, 
        "1_intermediate_data/taxa_table_ASVbased_Y1_AD_20220504_8.rds")
write_dta(taxa_ASVbased, 
          "1_intermediate_data/taxa_table_ASVbased_Y1_AD_20220504_8.dta")
write.csv(taxa_ASVbased, 
          "1_intermediate_data/taxa_table_ASVbased_Y1_AD_20220504_8.csv", 
          row.names = FALSE)

## SessionInfo ----
writeLines(capture.output(sessionInfo()), "4_output/0_data_cleaning_ASV_Y1_sessionInfo.txt")


# 3_variables_description ----
## Taxa raw abundances (ASV table) ----
taxa_raw_abundances_Y1_summary <- final_data %>%
  select(!contains("ASV")) %>%
  select(starts_with("ch_feces_raw")) %>%
  tbl_summary(
    missing = "no", 
    type = list(where(is.numeric)~ "continuous"), 
    statistic = all_continuous()~ "{min}/{p25}/{median}/{mean}/{p75}/{max}/{N_nonmiss}"
  ) %>%
  bold_labels() %>% 
  as_gt()  %>% 
  as.data.frame()%>% 
  select(label, stat_0, variable) %>%
  separate(
    col = stat_0, 
    into = c("Min", "Q1", "Median", "Mean", "Q3", "Max", "N"), 
    sep = "/", 
    remove = TRUE) %>%
  rename(
    "Variable names" = variable, 
    "Variable labels" = label
  ) %>%
  select("Variable names", "Variable labels", everything())

write.table(taxa_raw_abundances_Y1_summary, 
            "4_output/taxa_raw_abundances_ASVbased_Y1_summary.csv", 
            row.names = FALSE, 
            sep =";", 
            dec = ".")

## Taxa rel abundances (ASV based) ----
taxa_rel_abundances_Y1_summary <- final_data %>%
  select(!contains("ASV")) %>%
  select(starts_with("ch_feces_rel")) %>%
  tbl_summary(
    missing = "no", 
    type = list(where(is.numeric)~ "continuous"), 
    statistic = all_continuous()~ "{min}/{p25}/{median}/{mean}/{p75}/{max}/{N_nonmiss}"
  ) %>%
  bold_labels() %>% 
  as_gt()  %>% 
  as.data.frame()%>% 
  select(label, variable, stat_0) %>%
  separate(
    col = stat_0, 
    into = c("Min", "Q1", "Median", "Mean", "Q3", "Max", "N"), 
    sep = "/", 
    remove = TRUE) %>%
  rename(
    "Variable names" = variable, 
    "Variable labels" = label
  ) %>%
  select("Variable names", "Variable labels", everything())

write.table(taxa_rel_abundances_Y1_summary, 
            "4_output/taxa_rel_abundances_ASVbased_Y1_summary.csv", 
            row.names = FALSE, 
            sep =";", 
            dec = ".")

## RowSums at each taxonomic levels ----
final_data %>% select(contains("_rel_g")) %>% rowSums()
final_data %>% select(contains("_rel_f")) %>% rowSums()
final_data %>% select(contains("_rel_o")) %>% rowSums()
final_data %>% select(contains("_rel_c")) %>% rowSums()
final_data %>% select(contains("_rel_p")) %>% rowSums()

# Rowsums of the relative abundances at each taxonomic level are not exactly equal to 100.
# because the abundances were initially filtered at 0.7% at the family level by Dr. Lepage who provided us the rel metagenomic data 
