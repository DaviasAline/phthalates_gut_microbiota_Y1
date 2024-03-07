### Aline Davias
### 2022/05/04
### Alpha diversity indices in one year child gut microbiota (ASV based)

library(haven)
library(tidyverse)
library(gtsummary)
#install_phyloseq(branch = "github")
library(phyloseq)   # specnumber(), rarefy(), rerecurve(), prune_samples(), estimate_richness()
library(ape)
library(SRS)        # SRS()
library(picante)    # pd()
library(labelled)
library(expss)

# specnumber(): obtain the specific richness in ecological samples
# rarefy(): obtain the richness or diversity after rarefaction in ecological samples
# prune_sample(): keep only the samples with a depth of sequencing > a chosen threshold
# SRS(): randomly select sequences for each sample, recreate a ASV table with the same number of sequences for each sample
# otu_table() : indicate to the package phyloseq that this object is an ASV table
# tax_table(): indicate to the phyloseq package that this object is a taxonomic correspondence table
# samples_data(): indicates to the phyloseq package that this object is a metadata 
# phyloseq(): links phyloseq objects (ASV table, taxonomic correspondence, metadata and phylogenetic tree) 
# pd(): calcul of the Faith's phylogenetic diversity
# estimate_richness(): calculates the most common diversity alpha indices


# 0_source_data_reading ----
## Raw abundance ASV table of the gut microbiota at one year
asv_raw <- 
  read_labelled_csv(
    "1_intermediate_data/gut_microbiota_ASVbased_Y1_labelled_AD_20220504_7239.csv")%>%
  select(
    ident, 
    ch_feces_ID_Y1, 
    starts_with("ch_feces_raw_ASV"))

## taxonomic correspondence
taxa <- read_csv("1_intermediate_data/taxa_table_ASVbased_Y1_AD_20220504_8.csv")

## Metadata of gut microbiota at one year
metadata <- 
  read_labelled_csv("1_intermediate_data/gut_microbiota_ASVbased_Y1_labelled_AD_20220504_7239.csv")%>%
  select(!starts_with("ch_feces_raw"))%>%
  select(!starts_with("ch_feces_rel"))


# 1_data_cleaning ---- 
ASV_not_rarefied <- asv_raw %>%                                      # Raw ASVs table not rarefied
  select(-ident) %>%
  na.omit()
row.names(ASV_not_rarefied) <- NULL 
ASV_not_rarefied <- ASV_not_rarefied %>%
  column_to_rownames("ch_feces_ID_Y1") %>%
  t() %>%
  as.data.frame()%>% 
  rownames_to_column("ch_feces_ASV_ID_Y1") %>%
  mutate(
    ch_feces_ASV_ID_Y1 = str_replace_all(ch_feces_ASV_ID_Y1, c("ch_feces_raw_" = "", "_Y1"=""))) %>%
  column_to_rownames("ch_feces_ASV_ID_Y1") %>%
  otu_table(taxa_are_rows = TRUE)                                    # convert as phyloseq objet


# 2_variables_creation ----
## Rarefaction ----
### from the raw ASV table, choose thresholds for the sequencing depth 
### define a subset of samples to keep = samples with a sequencing depth > the chosen threshold  
### samples with a sequencing depth < the chosen threshold become missing data 
keep_Cmin <-
  names(which(sample_sums(ASV_not_rarefied)>= min(colSums(ASV_not_rarefied)))) %>%
  prune_samples(ASV_not_rarefied) %>% as.data.frame()                # Loss of 0 samples 

keep_5000 <- 
  names(which(sample_sums(ASV_not_rarefied)>= 5000))%>%
  prune_samples(ASV_not_rarefied) %>% as.data.frame()                # Loss of 6 samples 

keep_10000 <- names(which(sample_sums(ASV_not_rarefied)>= 10000))%>%
  prune_samples(ASV_not_rarefied) %>% as.data.frame()                # Loss of 17 samples 

### reduce the sequencing depth of the samples to the chosen threshold
### the sequences kept within each sample are randomly selected 
ASV_rarefied_Cmin_Y1 <- keep_Cmin %>%
  SRS(min(colSums(ASV_not_rarefied)), set_seed = TRUE, seed = 1)
rownames(ASV_rarefied_Cmin_Y1)<- rownames(keep_Cmin)

ASV_rarefied_5000_Y1 <- keep_5000 %>%
  SRS(5000, set_seed = TRUE, seed = 1)
rownames(ASV_rarefied_5000_Y1)<- rownames(keep_5000)

ASV_rarefied_10000_Y1 <- keep_10000 %>%
  SRS(10000, set_seed = TRUE, seed = 1)
rownames(ASV_rarefied_10000_Y1) <- rownames(keep_10000)

## Settings phyloseq package (taxonomic correspondence + metadata + ASV table + phylogenetic tree) ----
### taxonomic correspondence 
TAXA <- taxa %>% 
  column_to_rownames("ch_feces_ASV_ID_Y1") %>%
  select(-ch_feces_TAX_ASVbased_Y1) %>%
  as.matrix()%>%
  tax_table()   

### metadata of gut microbiota at one year
METADATA <- metadata %>% na.omit()
row.names(METADATA) <- NULL
METADATA <- METADATA %>%
  column_to_rownames("ch_feces_ID_Y1") %>%
  sample_data()

### Assembling in phyloseq objet
gumme_rarefied_Cmin_Y1 <- phyloseq(
  otu_table(ASV_rarefied_Cmin_Y1, taxa_are_rows = TRUE), 
  TAXA, 
  METADATA)

gumme_rarefied_5000_Y1 <- phyloseq(
  otu_table(ASV_rarefied_5000_Y1, taxa_are_rows = TRUE), 
  TAXA, 
  METADATA)

gumme_rarefied_10000_Y1 <- phyloseq(
  otu_table(ASV_rarefied_10000_Y1, taxa_are_rows = TRUE), 
  TAXA, 
  METADATA)

### phylogenetic trees 
TREE_rarefied_Cmin <- rtree(
  ntaxa(gumme_rarefied_Cmin_Y1),
  rooted=TRUE, 
  tip.label = taxa_names(gumme_rarefied_Cmin_Y1))

TREE_rarefied_5000 <- rtree(
  ntaxa(gumme_rarefied_5000_Y1), 
  rooted=TRUE, 
  tip.label = taxa_names(gumme_rarefied_5000_Y1))
                                                     
TREE_rarefied_10000 <- rtree(
  ntaxa(gumme_rarefied_10000_Y1), 
  rooted=TRUE, 
  tip.label = taxa_names(gumme_rarefied_10000_Y1))

## alpha diversity indices ----
### Faith's phylogenetic diversity calculated on the ASV table rarefied at the minimal threshold (3580)
faith_Cmin <-                                             
  ASV_rarefied_Cmin_Y1 %>%
  t()%>%
  pd(TREE_rarefied_Cmin, include.root = FALSE) %>%
  as.data.frame()%>%
  rownames_to_column(var = "ch_feces_ID_Y1")

### Most common alpha diversity indices calculated on the ASV table rarefied at the minimal threshold (3580)
### and merge with the Faith's diversity
rich_Cmin <- 
  estimate_richness(gumme_rarefied_Cmin_Y1) %>% 
  as.data.frame()%>%
  rownames_to_column(var = "ch_feces_ID_Y1")%>%
  merge(faith_Cmin, by="ch_feces_ID_Y1") %>%         
  select(-SR, -se.chao1, -se.ACE)%>%
  rename("Faith" = "PD")

### Faith's phylogenetic diversity calculated on the ASV table rarefied at the threshold 5000
faith_5000 <- 
  ASV_rarefied_5000_Y1 %>%                         
  t() %>% 
  pd(TREE_rarefied_5000, include.root = FALSE) %>%
  as.data.frame()%>%
  rownames_to_column(var = "ch_feces_ID_Y1")

### Most common alpha diversity indices calculated on the ASV table rarefied at the threshold 5000
### and merge with the Faith's diversity
rich_5000 <- 
  estimate_richness(gumme_rarefied_5000_Y1)%>%       
  as.data.frame() %>%
  rownames_to_column( var = "ch_feces_ID_Y1") %>%
  merge(faith_5000, by="ch_feces_ID_Y1")%>%                
  select(-SR, -se.chao1, -se.ACE) %>% 
  rename("Faith" = "PD")

### Faith's phylogenetic diversity calculated on the ASV table rarefied at the threshold 10 000
faith_10000 <- ASV_rarefied_10000_Y1 %>%    
  t()%>%
  pd(TREE_rarefied_10000, include.root = FALSE)%>%    
  as.data.frame()%>%
  rownames_to_column(var = "ch_feces_ID_Y1")

### Most common alpha diversity indices calculated on the ASV table rarefied at the threshold 10 000
### and merge with the Faith's diversity
rich_10000 <- estimate_richness(gumme_rarefied_10000_Y1) %>%                   
  as.data.frame() %>%                                                           
  rownames_to_column(var = "ch_feces_ID_Y1") %>%
  merge(faith_10000, by="ch_feces_ID_Y1")%>%                                     
  select(-SR, -se.chao1, -se.ACE)%>% 
  rename("Faith" = "PD")


# 3_variables_selection ----
## selection and labellisation ----
rich_Cmin <- rich_Cmin %>%
  rename(
    ch_feces_SpecRich_cmin_ASV_Y1 = Observed, 
    ch_feces_Chao1_cmin_ASV_Y1 = Chao1, 
    ch_feces_ACE_cmin_ASV_Y1 = ACE, 
    ch_feces_Shannon_cmin_ASV_Y1 = Shannon, 
    ch_feces_Simpson_cmin_ASV_Y1 = Simpson, 
    ch_feces_invSimpson_cmin_ASV_Y1 = InvSimpson, 
    ch_feces_Fisher_cmin_ASV_Y1 = Fisher, 
    ch_feces_Faith_cmin_ASV_Y1 = Faith)

rich_5000 <- rich_5000 %>%
  rename(
    ch_feces_SpecRich_5000_ASV_Y1 = Observed, 
    ch_feces_Chao1_5000_ASV_Y1 = Chao1, 
    ch_feces_ACE_5000_ASV_Y1 = ACE, 
    ch_feces_Shannon_5000_ASV_Y1 = Shannon, 
    ch_feces_Simpson_5000_ASV_Y1 = Simpson, 
    ch_feces_invSimpson_5000_ASV_Y1 = InvSimpson, 
    ch_feces_Fisher_5000_ASV_Y1 = Fisher, 
    ch_feces_Faith_5000_ASV_Y1 = Faith)

rich_10000 <- rich_10000 %>%
  rename(
    ch_feces_SpecRich_10000_ASV_Y1 = Observed, 
    ch_feces_Chao1_10000_ASV_Y1 = Chao1, 
    ch_feces_ACE_10000_ASV_Y1 = ACE, 
    ch_feces_Shannon_10000_ASV_Y1 = Shannon, 
    ch_feces_Simpson_10000_ASV_Y1 = Simpson, 
    ch_feces_invSimpson_10000_ASV_Y1 = InvSimpson, 
    ch_feces_Fisher_10000_ASV_Y1 = Fisher, 
    ch_feces_Faith_10000_ASV_Y1 = Faith)

ident_key <- metadata %>% 
  select(ident, ch_feces_ID_Y1) %>%
  mutate(
    ch_feces_ID_Y1 = as.character(ch_feces_ID_Y1))

final_data <- list(
    ident_key,
    rich_Cmin, 
    rich_5000, 
    rich_10000
    ) %>%
  reduce(left_join, by = "ch_feces_ID_Y1")

var_label(final_data) <- colnames(final_data) %>%
  str_replace_all(
    c("ch_feces_ID_Y1" = "Child feces identity at one year", 
      "ch_feces_" ="",
      "SpecRich" = "Specific richness",
      "_cmin_ASV_Y1" = " diversity in 1Y child gut microbiota ASV based (seq. depth = 3 580)", 
      "_5000_ASV_Y1" = " diversity in 1Y child gut microbiota ASV based (seq. depth = 5 000)", 
      "_10000_ASV_Y1" = " diversity in 1Y child gut microbiota ASV based (seq. depth = 10 000)"))

## export ----
saveRDS(final_data, 
          "2_final_data/alpha_diversity_ASVbased_Y1_AD_20220504_26.rds")

write_dta(final_data, 
          version = 13,
          "2_final_data/alpha_diversity_ASVbased_Y1_AD_20220504_26.dta")

write.csv(final_data, 
          "2_final_data/alpha_diversity_ASVbased_Y1_AD_20220504_26.csv", 
          row.names = FALSE) 

write_labelled_csv(final_data, 
                   "2_final_data/alpha_diversity_ASVbased_Y1_labelled_AD_20220504_26.csv", 
                   row.names = FALSE)  

## sessionInfo ----
writeLines(
  capture.output(sessionInfo()), 
  "4_output/1_alpha_diversity_ASVbased_Y1_sessionInfo.txt")


# 4_variables_description ----
## Descriptive table ----
alpha_diversity_Y1_summary <- final_data %>%
  select(-ch_feces_ID_Y1, -ident) %>%
  tbl_summary(
    missing = "no", 
    type = list(where(is.numeric)~ "continuous"), 
    statistic = all_continuous()~ "{min}/{p25}/{median}/{mean}/{p75}/{max}/{N_nonmiss}"
  ) %>%
  bold_labels() %>% 
  as_gt()  %>% 
  as.data.frame()%>% 
  select(variable, label, stat_0) %>%
  separate(
    col = stat_0, 
    into = c("Min", "Q1", "Median", "Mean", "Q3", "Max", "N"), 
    sep = "/", 
    remove = TRUE) %>%
  rename(
    "Variable names" = variable, 
    "Variable labels" = label)

write.table(alpha_diversity_Y1_summary, 
            "4_output/alpha_diversity_ASVbased_Y1_summary.csv", 
            row.names = FALSE, 
            sep =";", 
            dec = ".")

## Rarecurves ----
subsample <- asv_raw %>% 
  na.omit() %>%
  sample_n(100) %>%
  select(-ch_feces_ID_Y1)

dev.off()
alpha_diversity_Y1_rarecurve5000 <- 
  rarecurve(subsample, 
            step = 100, 
            sample = 5000, 
            col = "blue", 
            ylab = "Observed ASVs", 
            xlab = "Depth of sequencing")
dev.print(device = png, 
          file = "4_output/alpha_diversity_ASVbased_Y1_rarecurve5000.png", 
          width = 900, 
          height = 600)

## Scatterplot ----
data_plot_long <- 
  merge(final_data, 
        metadata, 
        by = "ch_feces_ID_Y1") %>%
  select(ch_feces_ID_Y1, 
         ch_feces_NASV_Y1,                                           # observed number of ASV before rarefaction
         ch_feces_SpecRich_cmin_ASV_Y1 ,                             # observed number of ASV after rarefaction at the minimal threshold (3580)
         ch_feces_SpecRich_5000_ASV_Y1,                              # observed number of ASV after rarefaction at the threshold 5000
         ch_feces_SpecRich_10000_ASV_Y1) %>%                         # observed number of ASV after rarefaction at the threshold 10000
  pivot_longer(
    cols = c("ch_feces_SpecRich_cmin_ASV_Y1", "ch_feces_SpecRich_5000_ASV_Y1", "ch_feces_SpecRich_10000_ASV_Y1"), 
    names_to = "rich_after_rar", 
    values_to = "value") %>%
  rename(
    rich_before_rar = ch_feces_NASV_Y1) %>%
  mutate(
    rich_after_rar = fct_recode(rich_after_rar,
                                "Threshold 10 000" = "ch_feces_SpecRich_10000_ASV_Y1",
                                "Threshold 5 000" = "ch_feces_SpecRich_5000_ASV_Y1",
                                "Minimal threshold 3580" = "ch_feces_SpecRich_cmin_ASV_Y1"), 
    rich_after_rar = fct_relevel(rich_after_rar,
                                 "Minimal threshold 3580", 
                                 "Threshold 5 000", 
                                 "Threshold 10 000"))

dev.off()
alpha_diversity_Y1_scatterplot <- data_plot_long %>%
  ggplot() +
  aes(x = value, y = rich_before_rar, colour = rich_after_rar) +
  geom_point(shape = "circle",
             size = 1.5) +
  geom_abline(intercept = 0, slope = 1) +
  scale_color_hue(direction = 1) +
  labs(x = "Observed ASVs after rarefaction", 
       y = "Observed ASVs before rarefaction",
       color = "Threshold of rarefation") +
  theme_classic()

alpha_diversity_Y1_scatterplot
dev.print(device = png, 
          file = "4_output/alpha_diversity_ASVbased_Y1_scatterplot.png", 
          width = 900, 
          height = 600)