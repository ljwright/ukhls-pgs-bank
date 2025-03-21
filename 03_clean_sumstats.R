library(tidyverse)
library(glue)
library(vroom)
library(fs)
library(curl)
library(tictoc)
library(summarytools)

rm(list = ls())

# 1. Set Up Environment and Create Helpers ----
# NOTE: Genetic data are build GRCh37 with rsid for SNP IDs.
# Some of the SNP IDs do not follow the rsid format.
# GRCh37 > GRCh38 > rsid

# --snp 0 --A1 1 --A2 2 --stat 3 --pvalue 4 --index 
tibble(line = read_lines("Code/02_download_sumstats.R")) %>%
  filter(str_detect(line, "^\\#"), 
         str_detect(line, "\\-\\-\\-\\-")) %>%
  pull(line) %>%
  glue_collapse("\n")

check_remaining <- function(){
  files_s <- list.files("Source Data", "\\.txt\\.gz$")
  file_c <- list.files("Data/step_03", "\\.txt\\.gz$")
  setdiff(files_s, file_c)
}

# dir_delete("Data/step_03")
# dir_create("Data/step_03")

check_dbsnp <- function(search_term){
  glue("https://www.ncbi.nlm.nih.gov/snp/?term={curl_escape(search_term)}") %>%
    browseURL()
}
check_dbsnp("rs188996809")

add_gz <- function(file){
  file <- case_when(str_detect(file, "\\.txt$") ~ glue("{file}.gz"),
                    str_detect(file, "\\.txt.gz$") ~ file,
                    TRUE ~ glue("{file}.txt.gz"))
  # if (!str_detect(file, "\\.txt\\.gz$")) file <- glue("{file}.txt.gz")
  return(file)
}

open_source <- function(file, col_nums = matches(".*"), n_max = Inf, skip = 0, delim = "\t"){
  vroom(glue("Source Data/{add_gz(file)}"),
        delim = delim, n_max = n_max, skip = skip, 
        col_select = col_nums,
        col_types = cols(.default = col_character()))
}

check_head_df <- function(file_name, col_nums = matches(".*")){
  open_source(file_name, col_nums, n_max = 50)
}

check_head <- function(file){
  glue("Source Data/{add_gz(file)}") %>%
    vroom_lines(n_max = 50) 
}


clean_sumstats <- function(df){
  df_clean <- df %>%
    select(id, effect_allele, other_allele, beta, p) %>%
    drop_na() %>%
    add_count(id) %>%
    filter(n == 1) %>%
    select(-n)
  
  return(df_clean)
}

save_sumstats <- function(df, file_name){
  write_tsv(clean_sumstats(df),
            glue("Data/step_03/{add_gz(file_name)}"), 
            col_names = FALSE)
  
  return(file_name)
}


# 2. Download Summary Statistics ----
## a. Health Behaviour ----
### i. Addictive Behaviour (Hatoum et al., 2023; 10.1038/s44220-023-00034-y) ----
open_source("hatoum_etal_2023_addictive_behaviour_grch38",
            c(2:3, 4:6), skip = 19) %>%
  rename(effect_allele = 3, other_allele = 4, beta = 5) %>%
  mutate(id = glue("chr{chr_name}:{chr_position}"),
         beta = as.double(beta),
         p = 1) %>%
  save_sumstats("hatoum_etal_2023_addictive_behaviour_grch37")

### ii. Alcohol and Cigarette Consumption (Liu et al., 2019; 10.1038/s41588-018-0307-5) ----
clean_liu_etal_2019 <- function(file){
  df <- open_source(file, c(1:2, 4:5, 8:9)) %>%
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
open_source("day_etal_2017_age_menarche_grch37", 1:5) %>%
  rename(id = 1, effect_allele = 2, other_allele = 3, beta = 4, p = 5) %>%
  mutate(across(c(effect_allele, other_allele), str_to_upper),
         across(c(beta, p), as.double)) %>%
  filter(effect_allele != "I", other_allele != "I") %>%
  save_sumstats("day_etal_2017_age_menarche_grch37")


### ii. Age at Menopause (Ruth et al., 2021; 10.1038/s41586-021-03779-7) ----
open_source("ruth_etal_2021_age_menopause_grch37", c(1, 4:5, 7, 9)) %>%
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
clean_roshandel_etal_2023 <- function(file){
  df <- open_source(file, c(1:5, 8)) %>%
    mutate(id = glue("chr{chromosome}:{base_pair_location}"),
           across(c(beta, neg_log_10_p_value), as.double),
           p = 10^(-neg_log_10_p_value))
  
  save_sumstats(df, file)
  rm(df)
  return(file)
}

list.files("Source Data", "roshandel_etal_2023") %>%
  walk(clean_roshandel_etal_2023)


### v. (Adult) Body Mass Index (Yengo et al., 2018; 10.1093/hmg/ddy271) ----
clean_yengo_etal_2018 <- function(file){
  df <- open_source(file, c(1:2, 4:5, 7, 9)) %>%
    rename_with(str_to_lower) %>%
    rename(effect_allele = tested_allele) %>%
    mutate(id = glue("chr{chr}:{pos}"),
           across(c(beta, p), as.double))
  
  save_sumstats(df, file)
  rm(df)
  return(file)
}

list.files("Source Data", "yengo_etal_2018") %>%
  walk(clean_yengo_etal_2018)


### vi. (Childhood) Body Mass Index (Vogelezang et al. 2020; 10.1371/journal.pgen.1008718) ----
open_source("vogelezang_etal_2020_bmi_grch37", c(2:6, 8)) %>%
  rename(p = p_value) %>%
  mutate(id = glue("chr{chromosome}:{base_pair_location}"),
         across(c(beta, p), as.double)) %>%
  save_sumstats("vogelezang_etal_2020_bmi_grch37")


### vii. Height (Yengo et al., 2018; 10.1093/hmg/ddy271) ----
### NOTE: SEE ABOVE.

### viii. Height (Yengo et al., 2022; 10.1038/s41586-022-05275-y) ----
clean_yengo_etal_2022 <- function(file){
  df <- open_source(file, c(3:6, 8, 10)) %>%
    rename_with(str_to_lower) %>%
    mutate(id = glue("chr{chr}:{pos}"),
           across(c(beta, p), as.double))
  
  save_sumstats(df, file)
  rm(df)
  return(file)
}

list.files("Source Data", "yengo_etal_2022") %>%
  walk(clean_yengo_etal_2022)

### ix. Grip Strength (Jones et al., 2021; 10.1038/s41467-021-20918-w) ----
clean_jones_etal_2021 <- function(file){
  df <- open_source(file, c(2:5, 10, 12)) %>%
    rename(p = p_value) %>%
    mutate(id = glue("chr{chromosome}:{base_pair_location}"),
           across(c(beta, p), as.double))
  
  save_sumstats(df, file)
  rm(df)
  return(file)
}

list.files("Source Data", "jones_etal_2021") %>%
  walk(clean_jones_etal_2021)


### x. Grip Strength (Tikkanen et al., 2018; 10.1038/s41598-018-24735-y) ----
## NOTE: NOT YET DONE.

## c. Physical Health and Biomarkers ----
### i Asthma (Han et al., 2020; org/10.1038/s41467-020-15649-3) ----
open_source("han_etal_2020_asthma_grch37", c(2:7)) %>%
  rename_with(str_to_lower) %>%
  rename(effect_allele = ea, other_allele = nea, beta = z) %>%
  mutate(id = glue("chr{chr}:{bp}"),
         across(c(beta, p), as.double)) %>%
  save_sumstats("han_etal_2020_asthma_grch37")


### ii. Blood Pressure (Keaton et al., 2024; 10.1038/s41588-024-01714-w) ----
clean_keaton_etal_2024 <- function(file){
  df <- open_source(file, c(1:5, 8)) %>%
    rename(p = p_value) %>%
    mutate(id = glue("chr{chromosome}:{base_pair_location}"),
           across(c(beta, p), as.double))
  
  save_sumstats(df, file)
  rm(df)
  return(file)
}

list.files("Source Data", "keaton_etal_2024") %>%
  walk(clean_keaton_etal_2024)


### iii. C-Reactive Protein (Said et al., 2022; 10.1038/s41467-022-29650-5) ----
open_source("said_etal_2022_crp_grch37", c(2:6, 8)) %>%
  rename(p = p_value) %>%
  mutate(id = glue("chr{chromosome}:{base_pair_location}"),
         across(c(beta, p), as.double)) %>%
  save_sumstats("said_etal_2022_crp_grch37")


### iv. Coronary Artery Disease (Aragam et al., 2022; 10.1038/s41588-022-01233-6) ----
clean_aragam_etal_2022 <- function(file){
  df <- open_source(file, c(1:5, 8)) %>% # Do I use OR or Beta?
    rename(p = p_value) %>%
    mutate(across(c(effect_allele, other_allele), str_to_upper),
           across(c(beta, p), as.double),
           id = glue("chr{chromosome}:{base_pair_location}"))
  
  save_sumstats(df, file)
  rm(df)
  return(file)
}

list.files("Source Data", "aragam_etal_2022") %>%
  walk(clean_aragam_etal_2022)


### v. Fasting Blood Glucose and Hba1c (Downie et al., 2022; 10.1007/s00125-021-05635-9) ----
# TODO: EMAILED AUTHOR AS SUMMARY STATISTICS DON'T CONTAIN EFFECT SIZES ?!


### vi. Hba1c (Sinnott-Armstrong et al., 2021; 10.1038/s41588-020-00757-z) ----
open_source("sinnotarmstrong_etal_2021_hba1c_grch37", c(1:2, 4:6, 8)) %>%
  rename(chr = 1, pos = 2, other_allele = 3, effect_allele  = 4, beta = 5, p = 6) %>% # Assumed which is the effect_allele
  mutate(id = glue("chr{chr}:{pos}"),
         across(c(beta, p), as.double)) %>%
  save_sumstats("sinnotarmstrong_etal_2021_hba1c_grch37")


### vii. Parental Lifespan (Timmers et al., 2019; 10.7554/eLife.39856) ----
open_source("timmers_etal_2019_parental_lifespan_grch37", c(3:6, 9, 11)) %>%
  rename(effect_allele = a1, other_allele = a0, beta = beta1) %>% 
  mutate(id = glue("chr{chr}:{pos}"),
         across(c(beta, p), as.double)) %>%
  save_sumstats("timmers_etal_2019_parental_lifespan_grch37")


### ix. Hypertension (Bi et al., 2020; 10.1016/j.ajhg.2020.06.003) ----
vroom("Source Data/bi_etal_2020_hypertension_rsid.txt.gz", 
      delim = " ",
      col_names = FALSE,
      col_select = c(3:7, 9)) %>% # SPACox Score is Effect Size, I think
  rename(chr = 1, pos = 2, other_allele = 3, effect_allele = 4, p = 5, beta = 6) %>% # Assumed which is the effect_allele
  mutate(id = glue("chr{chr}:{pos}"),
         across(c(beta, p), as.double)) %>%
  save_sumstats("bi_etal_2020_hypertension_grch37")


### x. Rheumatoid Arthritis (Ishigaki et al., 2022; 10.1038/s41588-022-01213-w) ----
clean_ishigaki_etal_2022 <- function(file){
  df <- open_source(file, c(2:6, 8)) %>%
    rename(p = p_value)  %>%
    mutate(id = glue("chr{chromosome}:{base_pair_location}"),
           across(c(beta, p), as.double))
  
  save_sumstats(df, file)
  rm(df)
  return(file)
}

list.files("Source Data", "ishigaki_etal_2022") %>%
  map(clean_ishigaki_etal_2022)


### xi. Type I Diabetes (Chiou et al., 2021; 10.1038/s41586-021-03552) ----
open_source("chiou_etal_2021_t1d_grch38", c(2:6, 8)) %>%
  rename(p = p_value) %>%
  mutate(id = glue("chr{chromosome}:{base_pair_location}"),
         across(c(beta, p), as.double)) %>%
  save_sumstats("chiou_etal_2021_t1d_grch38")


### xi. Type II Diabetes (Suzuki et al., 2024; 10.1038/s41586-024-07019-6) ----
clean_suzuki_etal_2024 <- function(file){
  df <- open_source(file, c("Chromsome", "Position", "EffectAllele",
                      "NonEffectAllele", "Pval", "Beta")) %>%
    rename_with(str_to_lower) %>%
    rename(chr = chromsome, pos = position,
           effect_allele = effectallele, other_allele = noneffectallele,
           p = pval, beta = beta) %>%
    mutate(id = glue("chr{chr}:{pos}"),
           across(c(beta, p), as.double))
  
  save_sumstats(df, file)
  rm(df)
  return(file)
}

list.files("Source Data", "suzuki_etal_2024") %>%
  map(clean_suzuki_etal_2024)


### xii. Waist-Hip Index (Christakoudi et al. 2021; 10.1038/s41598-021-89176-6) ----
clean_christakoudi_etal_2021 <- function(file){
  df <- open_source(file, c(1:5, 8)) %>%
    rename(p = p_value) %>%
    mutate(id = glue("chr{chromosome}:{base_pair_location}"),
           across(c(beta, p), as.double))
  
  save_sumstats(df, file)
  rm(df)
  return(file)
}

list.files("Source Data", "christakoudi_etal_2021") %>%
  map(clean_christakoudi_etal_2021)


## d. Mental Health and Cognition ----
### i. ADHD (Demontis et al., 2023; 10.1038/s41588-022-01285-8) ----
open_source("demontis_etal_2023_adhd_grch37",
            c(1, 3:5, 9, 11), delim = " ") %>%
  rename_with(str_to_lower) %>%
  rename(effect_allele = a1, other_allele = a2) %>%
  mutate(id = glue("chr{chr}:{bp}"),
         across(c(or, p), as.double),
         beta = log(or)) %>% # Available in OR, which can go in PRSice-2
  save_sumstats("demontis_etal_2023_adhd_grch37")


### ii. Alzheimer's Disease (Bellenguez et al., 2022; 10.1038/s41588-022-01024-z) ----
open_source("bellenguez_etal_2022_alzheimers_grch38", c(2:6, 11)) %>% # Also available in OR, which can go in PRSice-2
  rename(p = p_value) %>%
  mutate(id = glue("chr{chromosome}:{base_pair_location}"),
         across(c(beta, p), as.double)) %>%
  save_sumstats("bellenguez_etal_2022_alzheimers_grch38")


### iii. Panic Disorder (Forstner et al., 2021; 10.1038/s41380-019-0590-2) ----
open_source("forstner_etal_2021_panic_disorder_grch37", 
            c(1:2, 4:6, 8), skip = 72) %>%
  rename_with(str_to_lower) %>%
  rename(chr = 1, pos = 2, effect_allele = 3, other_allele = 4, p = 6) %>%
  mutate(id = glue("chr{chr}:{pos}"),
         across(c(beta, p), as.double)) %>%
  save_sumstats("forstner_etal_2021_panic_disorder_grch37")


### iv. Autism Spectrum Disorder (Grove et al., 2019; 10.1038/s41588-019-0344-8) ----
open_source("grove_etal_2019_autism_grch37", c(1, 3:5, 7, 9)) %>%
  rename_with(str_to_lower) %>%
  rename(effect_allele = a1, other_allele = a2) %>%
  mutate(id = glue("chr{chr}:{bp}"),
         across(c(or, p), as.double),
         beta = log(or)) %>%
  save_sumstats("grove_etal_2019_autism_grch37")


### v. Bipolar Disorder (Mullins et al., 2021; 10.1038/s41588-021-00857-4) ----
clean_mullins_etal_2021 <- function(file){
  df <- open_source(file, c(1:2, 4:6, 8), skip = 72) %>%
    rename_with(str_to_lower) %>%
    rename(chr = 1, effect_allele = 3, other_allele = 4, p = 6) %>%
    mutate(id = glue("chr{chr}:{pos}"),
           across(c(beta, p), as.double))
  
  save_sumstats(df, file)
  rm(df)
  return(file)
}

list.files("Source Data", "mullins_etal_2021") %>%
  map(clean_mullins_etal_2021)


### vi. Intelligence (Savage et al., 2018; 10.1038/s41588-018-0152-6) ----
open_source("savage_etal_2018_intelligence_grch37", c(2:5, 8, 10)) %>%
  rename(p = p_value) %>%
  mutate(across(c(effect_allele, other_allele), str_to_upper),
         id = glue("chr{chromosome}:{base_pair_location}"),
         across(c(beta, p), as.double)) %>%
  save_sumstats("savage_etal_2018_intelligence_grch37")
         

### v. Hippocampal Volume (Liu et al., 2023; 10.1038/s41588-023-01425-8) ----
# TODO: To decide.

### vi. Major Depressive Disorder (Howard et al., 2019; 10.1038/s41593-018-0326-7) ----
clean_howard_etal_2019 <- function(file){
  df <- open_source(file, c(1:3, 5, 7), delim = " ") %>%
    rename(id = 1, effect_allele = 2, other_allele = 3, beta = 4, p = 5) %>%
    mutate(across(c(effect_allele, other_allele), str_to_upper))
  
  save_sumstats(df, file)
  rm(df)
  return(file)
}

list.files("Source Data", "howard_etal_2019") %>%
  map(clean_howard_etal_2019)


### vii. Parkinson's Disease (Nalls et al., 2019; 10.1016/S1474-4422(19)30320-5) ----
open_source("nalls_etal_2019_parkinsons_grch37", c(1:3, 5, 7)) %>%
  rename_with(str_to_lower) %>%
  rename(id = snp, effect_allele = a1, other_allele = a2, beta = b) %>%
  mutate(across(c(beta, p), as.double)) %>%
  save_sumstats("nalls_etal_2019_parkinsons_grch37")


### viii. Schizophrenia (Trubetskoy et al., 2022; 10.1038/s41586-022-04434-5) ----
clean_trubetskov_etal_2022 <- function(file){
  df <- open_source(file, c(1, 3:5, 9, 11), skip = 73) %>%
    rename_with(str_to_lower) %>%
    rename(effect_allele = a1, other_allele = a2, p = pval) %>%
    mutate(id = glue("chr{chrom}:{pos}"),
           across(c(beta, p), as.double))
  
  save_sumstats(df, file)
  rm(df)
  return(file)
}

list.files("Source Data", "trubetskov_etal_2022") %>%
  map(clean_trubetskov_etal_2022)

## e. Personality ----
### i. Big 5 Personality Traits (Gupta et al., 2024; 10.1038/s41562-024-01951-3) ----
clean_gupta_etal_2024 <- function(file){
  df <- open_source(file, c(2:5, 7, 9)) %>%
    rename_with(str_to_lower) %>%
    rename(effect_allele = a1, other_allele = a2) %>%
    mutate(id = glue("chr{chr}:{bp}"),
           across(c(beta, p), as.double))
  
  save_sumstats(df, file)
  rm(df)
  return(file)
}

list.files("Source Data", "gupta_etal_2024") %>%
  map(clean_gupta_etal_2024)


## f. Social Outcomes ----
### i. Educational Attainmment (EA4) (Okbay et al., 2022; 10.1038/s41588-022-01016-z) ----
clean_okbay_etal_2022 <- function(file){
  df <- open_source(file, c(2:5, 7, 10)) %>%
      rename_with(str_to_lower) %>%
      mutate(id = glue("chr{chr}:{bp}"),
             across(c(beta, p), as.double))
  
  save_sumstats(df, file)
  rm(df)
  return(file)
}

c("okbay_etal_2022_eduyears_top_inc23andme_grch37",
  "okbay_etal_2022_eduyears_all_excl3andme_grch37") %>%
  map(clean_okbay_etal_2022)

open_source("okbay_etal_2022_eduyears_top_exclukhls_grch37",
            c(1:2, 5:6, 16:17)) %>%
  rename_with(str_to_lower) %>%
  rename(effect_allele = ea, other_allele = oa) %>%
  mutate(id = glue("chr{chr}:{bp}"),
         across(c(beta, p), as.double)) %>%
  save_sumstats("okbay_etal_2022_eduyears_top_exclukhls_grch37")


### ii. Educational Attainmment (EA2) (Okbay et al., 2016; 10.1038/nature17671) ----
clean_okbay_etal_2016 <- function(file){
  df <- open_source(file, c(2:5, 7, 9)) %>%
    rename_with(str_to_lower) %>%
    rename(effect_allele = a1, other_allele = a2, p = pval) %>%
    mutate(id = glue("chr{chr}:{pos}"),
           across(c(beta, p), as.double))
  
  save_sumstats(df, file)
  rm(df)
  return(file)
}

list.files("Source Data", "okbay_etal_2016") %>%
  map(clean_okbay_etal_2016)


### iii. Household Income (Hill et al., 2019; 10.1038/s41467-019-13585-5) ----
clean_hill_etal_2019 <- function(file){
  df <- open_source(file, c(1, 3:6, 8), delim = " ") %>%
    rename_with(str_to_lower) %>%
    rename(other_allele = non_effect_allele) %>%
    mutate(id = glue("chr{chr}:{bpos}"),
           across(c(beta, p), as.double))

  save_sumstats(df, file)
  rm(df)
  return(file)
}

list.files("Source Data", "hill_etal_2019") %>%
  map(clean_hill_etal_2019)


### iv. Occupational Status (Akimova et al., 2024; 10.1038/s41562-024-02076-3) ----
clean_akimova_etal_2024 <- function(file){
  df <- open_source(file, c(1:5, 8))  %>%
    rename(p = p_value) %>%
    mutate(id = glue("chr{chromosome}:{base_pair_location}"),
           across(c(beta, p), as.double))
  
  save_sumstats(df, file)
  rm(df)
  return(file)
}

list.files("Source Data", "akimova_etal_2024") %>%
  map(clean_akimova_etal_2024)
