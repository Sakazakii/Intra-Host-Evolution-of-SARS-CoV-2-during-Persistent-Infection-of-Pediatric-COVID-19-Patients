# **Intra-Host Evolution of SARS-CoV-2 during Persistent Infection of Pediatric COVID-19 Patients**

Charlie R. Boyle, Tien Doan, Estefany Rios-Guzman, Jessica Maciuch, Lacy M. Simons, Dulce Garcia, 
David B. Williams, Arghavan Alisoltani, Egon A. Ozer, Ramon Lorenzo-Redondo, and Judd F. Hultquist
  
## **Table of Contents**
### [1. Introduction](#-introduction)
### [2. Significance and Findings](#significance-and-findings)
### [3. Study Design](#study-design)
### [4. Dependencies](#dependencies)
### [5. Epidemiological Data Sets](#epidemiological-data-sets)
### [6. Phylogenetic Data Sets](#phylogenetic-data-sets)
### [7. Phylogenetic Analysis](#phylogenetic-analysis)
### [8. Shannon Entropy Data Sets](#shannon-entropy-data-sets)
### [9. Shannon Entropy Analysis](#shannon-entropy-analysis)
### [10. Supplemental Table 2](supplemental-table-2)

## **Introduction**

This repository was curated to store scripts needed to generate both main and supplementary figures in 
Boyle et. al. 2025 (unpublished). The following data are provided:

1. NCBI Accession IDs for sequences used for phylogenetic and phylodynamic analysis
   
2. Codes (.ipynb & .r) used for the analysis of the aforementioned data sets

## **Significance and Findings**
Pediatric patients, especially immunocompromised, exhibit sustained PCR positivity with SARS-CoV-2.

Persistent SARS-CoV-2 positivity generates immune espacape and antiviral resistance non-synonymous mutations.

## **Study Design**
Through a single center retrospective cohort study, we examined the SARS-CoV-2 genetic diversity over a 3-year period in Chicago, Illinois, USA. Between 2020 to 2023, clinical metadata of 20 patient encounters were captured through Ann and Robert H. Lurie Children's Hospital. From March 22nd 2020 to March 22nd 2023, 11,783 residual diagnostic swabs were collected from SARS-CoV-2-positive patients, of which 50 were successfully utilized for SARS-CoV-2 whole-genome sequencing.

## **Dependencies**

### **Python Dependencies**
- Pandas (v2.2.3)
- Numpy (v.2.2.2)
- matplotlib (v.3.10.0)
- biopython (v.1.85)
- scipy (v.1.15.2)

### **R Dependencies**
- readxl (v.1.4.3)
- dplyr (v.4.4.0)
- lubridate (v.1.9.2)
- ggplot2 (v3.5.2)
- tidyr (v.1.3.0)
- ggtree (v3.14.0)
- ape (v.5.8)
- ggtreeExtra (v.1.16.0)
- ggplot2 (v.3.5.1)
- ggc (v9.2.0)
- TreeTools (v1.10.0)
- phylotools (v1.9-16)
- MAFFT (v7.471)
- scorpio (v0.3.12)
- IQ-Tree (v2.0.5)
- fasttree (v2.1.10)
- treetime (v.0.11.4)
- tidyverse (v.2.0.0)
- lmerTest (v3.1-3)
- emmeans (v1.11.2)
- patchwork (v1.3.1)
- glmmTMB (v1.1.11)
- magrittr (v2.0.3)
- gggenes (v0.5.1)
- ggeffects (v2.3.0)

## **Epidemiological Data Sets**
### **Figures 1 & 2**
<a href="https://github.com/Sakazakii/Intra-Host-Evolution-of-SARS-CoV-2-during-Persistent-Infection-of-Pediatric-COVID-19-Patients/blob/main/Immunocompetent/Immunocompetent%20CSV"> Immunocompetent Host CSV </a>

<a href="https://github.com/Sakazakii/Intra-Host-Evolution-of-SARS-CoV-2-during-Persistent-Infection-of-Pediatric-COVID-19-Patients/blob/main/Immunocompromised/Immunocompromised.csv"> Immunocompromised Host CSV </a>
## **Phylogenetic Data Sets**
### **Supplemental Table 1**
<a href="https://github.com/Sakazakii/Intra-Host-Evolution-of-SARS-CoV-2-during-Persistent-Infection-of-Pediatric-COVID-19-Patients/blob/main/NCBI/SARS-CoV-2%20NCBI"> SARS-CoV-2 NCBI Information

## **Phylogenetic Analysis**
### **Installing Packages and Activating Environment**
```
module load anaconda3
module load gcc/9.2.0
conda create -c biconda -n dsbr_covid ivar mafft=7.471 fasttree=2.1.10 treetime
conda activate dsbr_covid
mkdir dsbr_covid
cd dsbr_covid
```

### **Concatenate the consensus sequences**
```
cat xxx.fasta > xxx.fa > xxx.fasta
```

### **Alignment**
```
mafft -auto xxx.fasta > xxx.mafft.fasta
```

### **Maximum Likelihood (ML) Phylogenetic Construction (IQTree2)**
```
fasttree -nt -gtr -gamma xxx.mafft.fasta > xxx.tre
```

### **Temporal ML Phylogenetic Construction**
```
treetime --tree xxx.tre --aln xxx.mafft.fasta --dates seqs/dates.csv --outdir treetime
```

## **Shannon Entropy Data Sets**
### **Figures 3 & 4 and Supplemental Figure 1 and Supplemental Table 2**
<a href="https://github.com/Sakazakii/Intra-Host-Evolution-of-SARS-CoV-2-during-Persistent-Infection-of-Pediatric-COVID-19-Patients/blob/main/Sh%20Variant%20Files"> Shannon Entropy Variant Files </a>

<a href="https://github.com/Sakazakii/Intra-Host-Evolution-of-SARS-CoV-2-during-Persistent-Infection-of-Pediatric-COVID-19-Patients/blob/main/Sh%20samtools_depth%20files"> Shannon Entropy samtools_depth files </a>

## **Shannon Entropy Analysis**
### **Load functions**
<a href="https://github.com/Sakazakii/Intra-Host-Evolution-of-SARS-CoV-2-during-Persistent-Infection-of-Pediatric-COVID-19-Patients/blob/main/figure_generation.Rmd"> Functions </a>

### **Load Libraries in R Studio**
```{r setup, include=TRUE, echo=FALSE}
# Set working directory
  setwd("xxx")

# Load libraries
  libraries <- c("tidyverse",
                 "magrittr",
                 "data.table",
                 "lme4",
                 "emmeans",
                 "patchwork",
                 "glmmTMB")
    
  suppressPackageStartupMessages(
    lapply(libraries, 
           require, 
           character.only = TRUE))
```

### **Aggregate and Format iVar Variants Data**
1) Pull variants.tsv file generated by iVar into one data table
2) Filter data points that pass F-test and with non-zero nucleotide frequency
3) Select only substitution mutations (drop insertion/deletion data)

```{r}
  dataDir <- "xxx/xxxx"

# Specify pattern to search for iVar output files
  pattern <- "variants.tsv"

# Grab file paths
  files <- list.files(path = dataDir,
                      pattern = pattern, 
                      recursive = T,
                      full.names = T)

# Grab sample IDs from folder name
  names <- list.files(path = dataDir,
                      pattern = pattern, 
                      recursive = T) |>
    str_split_i(pattern = "/", i = 1) |>
    str_split_i(pattern = "\\.", i = 1)
    

# Iteratively read files
  freq_data <- files |> purrr::map(\(x) fread(x))
  
# Name list objects with corresponding sample ID
  names(freq_data) <- names

# Combine files into one data table
  freq_data %<>% rbindlist(idcol = "sample_ID")
  
# Filter out entries that did not pass statistical test
  freq_data %<>% .[PASS == "TRUE",]

# Filter out insertion/deletion entries
  freq_data %<>% .[ALT %in% c("A", "C", "G", "T"),]

```

### **Aggregate & Format samtools Depth Data**
```{r}
  dataDir <- "xxx"

# Specify pattern to search for iVar output files
  pattern <- "depth.txt"

# Grab file paths
  files <- list.files(path = dataDir,
                      pattern = pattern, 
                      recursive = T,
                      full.names = T)

# Grab sample IDs from folder name
  names <- list.files(path = dataDir,
                      pattern = pattern, 
                      recursive = T) |>
    str_split_i(pattern = "\\.", i = 1)
    

# Iteratively read files
  depth_data <- files |> purrr::map(\(x) fread(x))
  
# Name list objects with corresponding sample ID
  names(depth_data) <- names

# Combine files into one data table
  depth_data %<>% rbindlist(idcol = "sample_ID")
  
# Change names  
  setnames(depth_data, c("sample_ID", "ref", "POS", "samtools_depth"))

  
# Merge samtools depth with freq_data
  freq_data <- merge(freq_data, depth_data, by = c("sample_ID", "POS"), all = TRUE)
```

### **Merge With Patient Metadata**

1) Load patient metadata
2) Remove empty columns
3) Recode immune status variable
4) Merge metadata with nucleotide frequency data
5) Format "Days" variable
```{r}

# Load patient data
  metadata <- fread("xxx.csv")

# Remove empty columns
  metadata %<>% .[, Filter(function(x) all(!is.na(x)), .SD)]

# Recode immune status variable
  metadata %<>% .[, immune_status := if_else(`Immune status` == "ICH",
                                             "Immunocompromised",
                                             "Immunocompetent")]
  
# Merge metadata with data
  freq_data %<>% merge(., 
                       metadata, 
                       by.x = "sample_ID", 
                       by.y = "LCH ID")
  
# Make sure "Days" variable is an integer for proper plotting and modeling
  freq_data %<>% .[, Days := Days %>% as.character() %>% as.integer()]
```

### **Merge With Remdesivir Administration Info**

1) Load file containing information on dates of treatment administration
2) Calculate days since end of last RDV administration
3) Create binary variable for whether date of last RDV administration <= 7

```
# Load treatment metadata
  treatment <- fread("xxx.csv")

# Calculate days since RDV
  freq_data$days_since_RDV <- apply(freq_data, 
                                    MARGIN = 1, 
                                    FUN = days_since_RDV, 
                                    treatment = treatment)
  
  freq_data %<>%
    .[, RDV := ifelse(!is.na(days_since_RDV) & days_since_RDV <= 7, 
                      "Yes", 
                      "No") %>% 
        as.factor()]
```
### **Determine genome region for each data entry**

1) Apply assign_feature() function to each POS value
2) Order "feature" factor to ensure proper plotting order

```{r}

# Load manually edited anno_db 
  anno_db <- fread("reference/SARS-CoV-2_genome_annotation_NC_045512(2).csv")

# Assign genome region
  freq_data %<>% 
    .[, feature := sapply(POS, assign_feature, anno_db = anno_db)]
  
# Order features as factor for proper plotting
  freq_data %<>% 
    .[, feature := factor(feature, levels = anno_db$feature)]
```

### **Total_Depth Filter**

1) Create a helper column called total_depth
2) Copy the value from TOTAL_DP (iVar) over
3) NAs value are POS that matches with the Reference base. Replace NAs from total_depth column with depth value from samtools_depth
```{r}
# Create helper column
  freq_data[, total_depth := TOTAL_DP]  # Copy values from TOTAL_DP

# If a position does not show up in variants file (TOTAL_DP is NA), replace NA by depth from samtools_depth column
  freq_data[is.na(total_depth), total_depth := samtools_depth]  

# Save data
  write_csv(freq_data,
            file = "xxx.csv")
  
  saveRDS(freq_data, file = "xxx.Rds")
```

### **Depth fiter for comparing across time points and between immune status**

1) Grouping by POS
2) Find the smallest total_depth > 200
3) Compute a lower bound as 90% of that value
4) Check if all POS in that group greater than or equal to this lower bound. If not discard the group
```{r}

freq_data_filtered <- freq_data[, {
  # Find all total_depth > 200
  valid_depths <- total_depth[total_depth > 200]
  
  if (length(valid_depths) == 0) {
    .SD[0]  # Discard group (no values > 200)
  } else {
    # Calculate lower bound from the minimum valid value
    lower_bound <- min(valid_depths) * 0.9
    
    # Keep group only if all total_depth ≥ lower_bound
    if (all(total_depth >= lower_bound, na.rm = TRUE)) {
      .SD
    } else {
      .SD[0]  # Discard group
    }
  }
}, by = POS]

# Save data
  write_csv(freq_data_filtered,
            file = "xxx.csv")
  
  saveRDS(freq_data_filtered, file = "xxx.Rds")

```
### **Calculate Shannon entropy for comparing across time points and between immune status**

1) Apply function calculate_shannon() per position and sample

```{r}

# After all the filtering, if a POS is retained and ALT_FREQ is NA, that POS matches the reference -> ALT_FREQ = 0
  freq_data_filtered[is.na(ALT_FREQ), ALT_FREQ := 0]
  
# Retain only POS where at least one time point at any patient has ALT_FREQ >= 0.03 
  freq_data_filtered <-  freq_data_filtered[, if (any(ALT_FREQ >= 0.03, na.rm = TRUE)) .SD, by = POS]

# Calculate shannon entropy
  shannon_data <- freq_data_filtered[, 
                            shannon_entropy := calculate_shannon(ALT_FREQ), 
                            by = c("POS",
                                   "sample_ID"
                                   )]
  
# Replace shannon entropy that is NA to 0
  shannon_data[is.na(shannon_entropy), shannon_entropy := 0]
  
# If a POS has >=2 ALT, it will show as duplicate. Remove one row 
  shannon_data <- unique(shannon_data, by = c("sample_ID", "POS"))
  
# Save data
  write_csv(shannon_data,
            file = "filtered_data/COVID_LCH_shannon_data.csv")
  
  saveRDS(shannon_data, file = "filtered_data/COVID_LCH_shannon_data.Rds")
```

### **Depth fiter for only comparing across time points within each participants**

1) Grouping by POS and MRN
2) Find the smallest total_depth > 400
3) Compute a lower bound as 90% of that value
4) Check if all POS in that group greater than or equal to this lower bound. If not discard the group

```{r}

freq_data_filtered_timepoint <- freq_data[, {
  # Find all total_depth > 400
  valid_depths <- total_depth[total_depth > 400]
  
  if (length(valid_depths) == 0) {
    .SD[0]  # Discard group (no values > 400)
  } else {
    # Calculate lower bound from the minimum valid value
    lower_bound <- min(valid_depths) * 0.9
    
    # Keep group only if all total_depth ≥ lower_bound
    if (all(total_depth >= lower_bound, na.rm = TRUE)) {
      .SD
    } else {
      .SD[0]  # Discard group
    }
  }
}, by = .(POS, MRN)]

# Save data
  write_csv(freq_data_filtered_timepoint,
            file = "xxx.csv")
  
  saveRDS(freq_data_filtered_timepoint, file = "xxx.Rds")

```

### **Calculate Shannon entropy for only comparing across time points within each participants**

1) Apply function calculate_shannon() per position and sample

```{r}

# After all the filtering, if a POS is retained and ALT_FREQ is NA, that POS matches the reference -> ALT_FREQ = 0
  freq_data_filtered_timepoint[is.na(ALT_FREQ), ALT_FREQ := 0]
  
# Retain only POS where at least one time point at any patient has ALT_FREQ >= 0.03
  freq_data_filtered_timepoint <-  freq_data_filtered_timepoint[, if (any(ALT_FREQ >= 0.03, na.rm = TRUE)) .SD, by = .(MRN,POS)]

# Calculate shannon entropy
  shannon_data_timepoint <- freq_data_filtered_timepoint[, 
                            shannon_entropy := calculate_shannon(ALT_FREQ), 
                            by = c("POS",
                                   "sample_ID"
                                   )]
# Replace shannon entropy that is NA to 0
  shannon_data_timepoint[is.na(shannon_entropy), shannon_entropy := 0]

# If a POS has >2 ALT, it will show as duplicate. Remove one row  
  shannon_data_timepoint <- unique(shannon_data_timepoint, by = c("sample_ID", "POS"))
  
# Save data
  write_csv(shannon_data_timepoint,
            file = "xxx.csv")
  
  saveRDS(shannon_data_timepoint, file = "xxx.Rds")
```

### **Whole genome Shannon entropy model**

1) Sum Shannon entropy across whole genome per sample
2) Check linearity of associations between predictors and outcome and transform if necessary
3) Run model predicting WG Shannon entropy
4) Evaluate model for normality of residuals

```{r}

# Sum shannon entropy across whole genome
  keep <- c("sample_ID",
            "MRN",
            "Ct",
            "RDV",
            "immune_status",
            "Days")

  wg_df <- shannon_data[, .(total_entropy = sum(shannon_entropy)), by = keep]
  
# Save data set
  write_csv(wg_df,
            file = "filtered_data/COVID_LCH_whole_genome_Shannon_entropy.csv")
  
  saveRDS(wg_df, file = "filtered_data/COVID_LCH_whole_genome_Shannon_entropy.Rds")
  
  
# Check linearity of associations
  ggplot(wg_df, 
         aes(x = Days,
             y = log(total_entropy))) + 
    geom_point() -> p1
  
# Run model 
  model <- lmer(log(total_entropy) ~ Ct + RDV + immune_status * Days + (1|MRN), data = wg_df)
  
# QQ plot
  qqnorm(resid(model))
  qqline(resid(model), col = "dodgerblue", lwd = 2)
  
# Save model
  saveRDS(model, file = "filtered_models/xxxx.Rds")

```

### **Gene Shannon entropy model**

```{r}

# Summing per gene 
  keep <- c("sample_ID",
            "feature",
            "MRN",
            "Ct",
            "RDV",
            "immune_status",
            "Days")

  gene_df <- shannon_data[, .(total_entropy = sum(shannon_entropy)), by = keep]
  
# Filter data points not assigned to gene
  gene_df %<>% .[!is.na(feature),]

# Pull sample metadata
  metadata <- gene_df[, sample_ID:Days] %>% 
    .[, feature := NULL] %>%
    unique()
  
# Set order of features 
  anno_db <- fread("reference/SARS-CoV-2_genome_annotation_NC_045512(2).csv")
  gene_df$feature %<>% factor(., levels = anno_db$feature) %>% droplevels()
  
# Save data set
  write_csv(gene_df,
            file = "filtered_data/xxxx.csv")
  
  saveRDS(gene_df, file = "filtered_data/xxxx.Rds")
 
```

### **Run Model**

1) Check linearity of associations between predictors and outcome
2) Split data set by gene
3) Run separate model total Shannon entropy for each gene
4) Evaluate model for normality of residuals
5) Generate pairwise contrasts to get immune status p-values per gene

```{r}

# Check linearity of associations
  ggplot(gene_df, 
         aes(x = Days,
             y = total_entropy,
             color = immune_status)) + 
    geom_point() + 
    geom_jitter() +
    facet_wrap(vars(feature), nrow = 4)


  ggplot(gene_df, 
         aes(x = Ct,
             y = total_entropy,
             color = immune_status)) + 
    geom_point() +
    geom_jitter() +
    facet_wrap(vars(feature), nrow = 4)

# Remove genes that all gene-entropy values is 0
  nonzero_gene_df <- gene_df[, if (sum(total_entropy, na.rm = TRUE) != 0) .SD, by = feature]
  
# drop empty groups
  nonzero_gene_df[, feature := droplevels(feature)]
  
# Split data set by gene
  df_list <- nonzero_gene_df %>% split(f = .[["feature"]])
 
# Generate models
  models <- df_list |> purrr::map(\(x) generate_model(df = x,
                                                      y = "total_entropy",
                                                      predictors = c("immune_status",
                                                                     "Ct",
                                                                     "Days"),
                                                      interaction = c("immune_status:Days"),
                                                      random = "MRN",
                                                      add_RDV = T))
  
# Evaluate normality of residuals
  list(models, names(models)) |> 
    purrr::pmap(\(model, title) assess_QQ_plot(model, title))
  
# Extract p-values
  pvals <- lapply(models, function(m) {
  coef(summary(m))[5, 5]
})
  
# Correct p-values for multiple hypothesis testing
  pvals %<>% p.adjust(method = "fdr")
  
# Save results
  saveRDS(pvals, file = "filtered_models/xxxx.Rds")

```
### **Figure generation**
<a href="https://github.com/Sakazakii/Intra-Host-Evolution-of-SARS-CoV-2-during-Persistent-Infection-of-Pediatric-COVID-19-Patients/blob/main/figure_generation.Rmd"> Figure generation </a>

## **Supplemental Table 2**
<a href="https://github.com/Sakazakii/Intra-Host-Evolution-of-SARS-CoV-2-during-Persistent-Infection-of-Pediatric-COVID-19-Patients/tree/main/Supplemental%20Table%202"> Supplemental Table 2 </a>
