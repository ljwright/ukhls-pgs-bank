library(tidyverse)
library(glue)
library(vroom)
library(fs)
library(curl)
library(tictoc)

rm(list = ls())

# 1. Set Up Environment and Create Helpers ----
# NOTE: Genetic data are build GRCh37 with rsid for SNP IDs.
# Some of the SNP IDs do not follow the rsid format.
# GRCh37 > rsid > GRCh38

# --snp 0 --A1 1 --A2 2 --stat 3 --pvalue 4 --index 
tibble(line = read_lines("Code/02_download_sumstats.R")) %>%
  filter(str_detect(line, "^\\#"), 
         str_detect(line, "\\-\\-\\-\\-")) %>%
  pull(line) %>%
  glue_collapse("\n")

# dir_delete("Data/step_03")
# dir_create("Data/step_03")

check_dbsnp <- function(search_term){
  glue("https://www.ncbi.nlm.nih.gov/snp/?term={curl_escape(search_term)}") %>%
    browseURL()
}
check_dbsnp("rs188996809")

clean_sumstats <- function(df){
  df %>%
    select(id, effect_allele, other_allele, beta, p) %>%
    drop_na() %>%
    add_count(id) %>%
    filter(n == 1) %>%
    select(-n)
}

save_sumstats <- function(df, file_name){
  if (!str_detect(file_name, "\\.txt\\.gz$")) file_name <- glue("{file_name}.txt.gz")
  
  df %>%
    clean_sumstats() %>%
    write_tsv(glue("Data/step_03/{file_name}"), col_names = FALSE)
  
  return(file_name)
}

check_head <- function(file){
  if (!str_detect(file, "\\.txt\\.gz$")) file <- glue("{file}.txt.gz")
  
  glue("Source Data/{file}") %>%
    vroom_lines(n_max = 50) 
}

check_head_df <- function(file){
  if (!str_detect(file, "\\.txt\\.gz$")) file <- glue("{file}.txt.gz")
  
  glue("Source Data/{file}") %>%
    vroom(n_max = 50) 
}

# 2. Download Summary Statistics ----
## a. Health Behaviour ----
### i. Addictive Behaviour (Hatoum et al., 2023; 10.1038/s44220-023-00034-y) ----
vroom("Source Data/hatoum_etal_2023_addictive_behaviour_grch38.txt.gz",
      skip = 19, delim = "\t", 
      col_select = c(2:3, 4:6),
      col_types = cols(.default = col_character())) %>%
  rename(chr = 1, pos = 2, effect_allele = 3, other_allele = 4, beta = 5) %>%
  mutate(id = glue("chr{chr}:{pos}"),
         beta = as.double(beta),
         p = 0) %>%
  save_sumstats("hatoum_etal_2023_addictive_behaviour_grch37")

### ii. Alcohol and Cigarette Consumption (Liu et al., 2019; 10.1038/s41588-018-0307-5) ----
clean_liu_etal_2019 <- function(file){
  df <- glue("Source Data/{file}") %>%
    vroom(delim = "\t", # n_max = 100,
        col_select = c(1:2, 4:5, 8:9),
        col_types = cols(.default = col_character())) %>%
    rename_with(str_to_lower) %>%
    rename(effect_allele = alt, other_allele = ref, p = pvalue) %>%
    mutate(id = glue("chr{chrom}:{pos}"),
           across(c(beta, p), as.double))
  
  save_sumstats(df, file)
  rm(df)
  return(file)
}

list.files("Source Data", "liu_etal") %>%
  walk(clean_liu_etal_2019)


### iii. Diet (Cole et al., 2020; 10.1038/s41467-020-15193-0) ----
### NOTE: NOT YET DONE.

## b. Anthropometrics ----
### i. Age at Menarche (Day et al., 2017; 10.1038/ng.3841) ----
vroom("Source Data/day_etal_2017_age_menarche_grch37.txt.gz",
      delim = "\t", 
      col_select = c(1:5),
      col_types = cols(.default = col_character())) %>%
  rename(id = 1, effect_allele = 2, other_allele = 3, beta = 4, p = 5) %>%
  mutate(across(c(effect_allele, other_allele), str_to_upper),
         across(c(beta, p), as.double)) %>%
  filter(effect_allele != "I", other_allele != "I") %>%
  save_sumstats("day_etal_2017_age_menarche_grch37")


### ii. Age at Menopause (Ruth et al., 2021; 10.1038/s41586-021-03779-7) ----
vroom("Source Data/ruth_etal_2021_age_menopause_grch37.txt.gz",
      delim = "\t",
      col_select = c(1, 4:5, 7, 9),
      col_types = cols(.default = col_character())) %>%
  rename(id = 1, effect_allele = 2, other_allele = 3, beta = 4, p = 5) %>%
  mutate(across(c(effect_allele, other_allele), str_to_upper),
         across(c(beta, p), as.double)) %>%
  filter(effect_allele != "I", other_allele != "I") %>%
  save_sumstats("ruth_etal_2021_age_menopause_grch37")


### iii Birth Weight (Warrington et al., 2019; 10.1038/s41588-019-0403-1) ----
clean_warrington_etal_2019 <- function(file){
  df <- glue("Source Data/{file}") %>%
    vroom_lines() %>%
    str_squish() %>%
    read_table(col_types = cols(.default = col_character())) %>%
    select(chr, pos, effect_allele = ea, other_allele = nea, beta, p) %>%
    mutate(id = glue("chr{chr}:{pos}"),
           across(c(effect_allele, other_allele), str_to_upper),
           across(c(beta, p), as.double))
  
  save_sumstats(df, file)
  rm(df)
  return(file)
}

list.files("Source Data", "warrington_etal_2019") %>%
  walk(clean_warrington_etal_2019)


### iv. Body Fat Percentage (Roshandel et al., 2023; https://doi.org/10.3389/fendo.2023.1274791) ----
check_head_df("roshandel_etal_2023_body_fat_female_grch37")


### v. (Adult) Body Mass Index (Yengo et al., 2018; 10.1093/hmg/ddy271) ----


### vi. (Childhood) Body Mass Index (Vogelezang et al. 2020; 10.1371/journal.pgen.1008718) ----


### vii. Height (Yengo et al., 2018; 10.1093/hmg/ddy271) ----
### NOTE: SEE ABOVE.

### viii. Height (Yengo et al., 2022; 10.1038/s41586-022-05275-y) ----
### ix. Grip Strength (Jones et al., 2021; 10.1038/s41467-021-20918-w) ----
### x. Grip Strength (Tikkanen et al., 2018; 10.1038/s41598-018-24735-y) ----
## c. Physical Health and Biomarkers ----
### i Asthma (Han et al., 2020; org/10.1038/s41467-020-15649-3) ----
### ii. Blood Pressure (Keaton et al., 2024; 10.1038/s41588-024-01714-w) ----
### iii. C-Reactive Protein (Said et al., 2022; 10.1038/s41467-022-29650-5) ----
### iv. Coronary Artery Disease (Aragam et al., 2022; 10.1038/s41588-022-01233-6) ----
### v. Fasting Blood Glucose and Hba1c (Downie et al., 2022; 10.1007/s00125-021-05635-9) ----
### vi. Hba1c (Sinnott-Armstrong et al., 2021; 10.1038/s41588-020-00757-z) ----
### vii. Parental Lifespan (Timmers et al., 2019; 10.7554/eLife.39856) ----
### ix. Hypertension (Bi et al., 2020; 10.1016/j.ajhg.2020.06.003) ----
### x. Rheumatoid Arthritis (Ishigaki et al., 2022; 10.1038/s41588-022-01213-w) ----
### xi. Type I Diabetes (Chiou et al., 2021; 10.1038/s41586-021-03552) ----
### xi. Type II Diabetes (Suzuki et al., 2024; 10.1038/s41586-024-07019-6) ----
### xii. Waist-Hip Index (Christakoudi et al. 2021; 10.1038/s41598-021-89176-6) ----
## d. Mental Health and Cognition ----
### i. ADHD (Demontis et al., 2023; 10.1038/s41588-022-01285-8) ----
### ii. Alzheimer's Disease (Bellenguez et al., 2022; 10.1038/s41588-022-01024-z) ----
### iii. Panic Disorder (Forstner et al., 2021; 10.1038/s41380-019-0590-2) ----
### iv. Autism Spectrum Disorder (Grove et al., 2019; 10.1038/s41588-019-0344-8) ----
### v. Bipolar Disorder (Mullins et al., 2021; 10.1038/s41588-021-00857-4) ----
### vi. Intelligence (Savage et al., 2018; 10.1038/s41588-018-0152-6) ----
### v. Hippocampal Volume (Liu et al., 2023; 10.1038/s41588-023-01425-8) ----
### vi. Major Depressive Disorder (Howard et al., 2019; 10.1038/s41593-018-0326-7) ----
### vii. Parkinson's Disease (Nalls et al., 2019; 10.1016/S1474-4422(19)30320-5) ----
### vii. Schizophrenia (Trubetskoy et al., 2022; 10.1038/s41586-022-04434-5) ----
## e. Personality ----
### i. Big 5 Personality Traits (Gupta et al., 2024; 10.1038/s41562-024-01951-3) ----
## f. Social Outcomes ----
### i. Educational Attainmment (EA4) (Okbay et al., 2022; 10.1038/s41588-022-01016-z) ----
### ii. Educational Attainmment (EA2) (Okbay et al., 2016; 10.1038/nature17671) ----
### iii. Household Income (Hill et al., 2019; 10.1038/s41467-019-13585-5) ----
### iv. Occupational Status (Akimova et al., 2024; 10.1038/s41562-024-02076-3) ----

# 2. 