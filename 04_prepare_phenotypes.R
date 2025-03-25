library(tidyverse)
library(labelled)
library(glue)
library(haven)
library(summarytools)

# 1. Load Data and Create Helpers ----
dac_v1_fld <- Sys.getenv("ukhls_dac_v1_fld")
dac_v2_fld <- Sys.getenv("ukhls_dac_v2_fld")

df_source <- list(v1 = dac_v1_fld, v2 = dac_v2_fld) %>% 
  map(~ glue("{.x}/OMIC-54-WRIGHT_sendout.dta") %>%
        read_dta())

df_raw <- df_source$v2 %>%
  full_join(df_source$v1, 
            suffix = c("", "_v1"),
            by = "id") %>%
  select(-matches("_v1$"))

lookup <- lookfor(df_raw) %>%
  as_tibble() %>%
  mutate(variable_low = str_to_lower(variable),
         label_low = str_to_lower(variable))

negative_to_na <- function(x){
  na_range(x) <- c(-Inf, -1)
  user_na_to_na(x)
}

values_to_na <- function(x, values){
  na_values(x) <- values
  user_na_to_na(x)
}

extract_pca <- function(df){
  df %>%
    drop_na() %>%
    PCA() %>%
    predict(newdata = df) %>%
    pluck("coord") %>%
    as_tibble() %>%
    pull(1) %>%
    scale() %>%
    as.double()
}

clean_condition <- function(ever, current){
  case_when(current == 1 ~ "Current",
            current == 2 ~ "Past",
            ever == 0 ~ "Never") %>%
    factor(levels = c("Never", "Past", "Current"))
}

find_files <- function(gwas) list.files("Source Data", gwas)

# tibble(line = read_lines("Code/03_prepare_sumstats.R")) %>%
#   filter(str_detect(line, "^\\#"),
#          str_detect(line, "\\-\\-\\-\\-$")) %>%
#   pull(line) %>%
#   glue_collapse("\n")

# 2. Download Summary Statistics ----
clean <- list()

## a. Health Behaviour ----
### i. Addictive Behaviour (Hatoum et al., 2023; 10.1038/s44220-023-00034-y) ----
### ii. Alcohol and Cigarette Consumption (Liu et al., 2019; 10.1038/s41588-018-0307-5) ----
clean$smoke_alc <- df_raw %>%
  mutate(smoke_status_w2 = case_when(b_smnow == 1 ~ "Current",
                                     b_smcigs == 2 ~ "Occasional",
                                     b_smcigs == 1 ~ "Former",
                                     b_smcigs == 3 ~ "Tried",
                                     b_smever == 2 ~ "Never") %>%
           factor(c("Current", "Occasional", "Former", "Tried", "Never")),
         
         smoke_ncigs_w2 = case_when(smoke_status_w2 == "Current" & between(b_ncigs, 0, 80) ~ b_ncigs,
                                    smoke_status_w2 == "Former" & between(b_smncigs, 0, 80) ~ b_smncigs,
                                    smoke_status_w2 %in% c("Occasional", "Tried", "Never") ~ 0),
         
         smoke_startage_w2 = ifelse(smoke_status_w2 %in% c("Current", "Former") & between(b_smagbg, 10, 35), b_smagbg, NA),
         
         alcohol_freq_w2 = case_when(b_sceverdrnk == 2 ~ 9,
                                     b_scfalcdrnk %in% 1:8 ~ b_scfalcdrnk,
                                     b_scalcl7d == 9 ~ 8) %>%
           ordered(1:9, 
                   c("~ 7 times per week", "5-6 times per week",
                     "3-4 times per week", "1-2 times per week",
                     "1-2 times per month", "Every couple months",
                     "1-2 times per year", "0 times past year",
                     "Never")) %>%
           fct_rev()
  ) %>%
  select(id, matches("_w2$"))


### iii. Diet (Cole et al., 2020; 10.1038/s41467-020-15193-0) ----
## b. Anthropometrics ----
clean$anthro <- df_raw %>%
  mutate(
    ### i. Age at Menarche (Day et al., 2017; 10.1038/ng.3841) ----
    ### ii. Age at Menopause (Ruth et al., 2021; 10.1038/s41586-021-03779-7) ----
    ### iii Birth Weight (Warrington et al., 2019; 10.1038/s41588-019-0403-1) ----
    # TODO: Make data request.
    
    ### iv. Body Fat Percentage (Roshandel et al., 2023; https://doi.org/10.3389/fendo.2023.1274791) ----
    across(matches("(b|c)_bfpcval"), negative_to_na),
    body_fat_nv = coalesce(b_bfpcval, c_bfpcval),
    
    ### v. (Adult) Body Mass Index (Yengo et al., 2018; 10.1093/hmg/ddy271) ----
    # TODO: Add b_bmival to data request
    # across(c(b_bmi_val, c_bmival, a_bmi_dv), negative_to_na),
    # bmi_nv = coalesce(b_bmival, c_bmival, a_bmi_dv),
    across(c(c_bmival, a_bmi_dv), negative_to_na),
    bmi_xwave = coalesce(c_bmival, a_bmi_dv),
    
    ### vi. (Childhood) Body Mass Index (Vogelezang et al. 2020; 10.1371/journal.pgen.1008718) ----
    ### vii. Height (Yengo et al., 2018; 10.1093/hmg/ddy271) ----
    # TODO: Make data request. 
    
    ### viii. Height (Yengo et al., 2022; 10.1038/s41586-022-05275-y) ----
    # TODO: Make data request. 
    
    ### ix. Grip Strength (Jones et al., 2021; 10.1038/s41467-021-20918-w) ----
    across(matches("(b|c)_mmgs(d|n)val"), negative_to_na),
    grip_dominant_nv = coalesce(b_mmgsdval, c_mmgsdval),
    grip_nondominant_nv = coalesce(b_mmgsnval, c_mmgsnval),
    
    ### x. Grip Strength (Tikkanen et al., 2018; 10.1038/s41598-018-24735-y) ----
    # SEE ABOVE
  ) %>%
  select(id, matches("_(nv|xwave)$"))


## c. Physical Health and Biomarkers ----
clean$ph <- x

df_raw %>%
  mutate(
    ### i Asthma (Han et al., 2020; org/10.1038/s41467-020-15649-3) ----
    asthma_w1 = clean_condition(a_hcond1, a_hconds01),
    
    ### ii. Blood Pressure (Keaton et al., 2024; 10.1038/s41588-024-01714-w) ----
    across(matches("(b|c)_(omsysval|omdiaval)"), negative_to_na),
    sbp_nv = coalesce(b_omsysval, c_omsysval),
    dbp_nv = coalesce(b_omdiaval, c_omdiaval),
    
    ### iii. C-Reactive Protein (Said et al., 2022; 10.1038/s41467-022-29650-5) ----
    # TODO: Make data request. 
    
    
    ### iv. Coronary Artery Disease (Aragam et al., 2022; 10.1038/s41588-022-01233-6) ----
    cvd_coronary_disease_w1 = clean_condition(a_hcond4, a_hconds04),
    
    ### v. Fasting Blood Glucose and Hba1c (Downie et al., 2022; 10.1007/s00125-021-05635-9) ----
    # TODO: Make data request. 
    
    
    ### vi. Hba1c (Sinnott-Armstrong et al., 2021; 10.1038/s41588-020-00757-z) ----
    # TODO: Make data request. 
    
    
    ### vii. Parental Lifespan (Timmers et al., 2019; 10.7554/eLife.39856) ----
    mother_alive_w1 = ifelse(a_lvrel1 %in% 0:1, a_lvrel1, NA),
    mother_age_w1 = ifelse(mother_alive_w1 == 0, a_maage, NA),
    father_alive_w1 = ifelse(a_lvrel2 %in% 0:1, a_lvrel2, NA),
    father_age_w1 = ifelse(father_alive_w1 == 0, a_paage, NA),
    
    ### ix. Hypertension (Bi et al., 2020; 10.1016/j.ajhg.2020.06.003) ----
    cvd_hypertension_w1 = clean_condition(a_hcond16, a_hconds16),
    cvd_heart_attack_w1 = clean_condition(a_hcond6, a_hconds06),
    cvd_stroke_w1 = clean_condition(a_hcond7, a_hconds07),
    cvd_heart_failure_w1 = clean_condition(a_hcond3, a_hconds03),
    cvd_any_w1 = pick(matches("^cvd_.*_w1$")) %>%
      map(as.integer) %>%
      reduce(~ pmax(.x, .y)) %>%
      factor(1:3, c("Never", "Past", "Current")),
    
    ### x. Rheumatoid Arthritis (Ishigaki et al., 2022; 10.1038/s41588-022-01213-w) ----
    cvd_hypertension_w1 = clean_condition(a_hcond2, a_hconds02),
    
    ### xi. Type I Diabetes (Chiou et al., 2021; 10.1038/s41586-021-03552) ----
    # TODO: Data request
    
    ### xi. Type II Diabetes (Suzuki et al., 2024; 10.1038/s41586-024-07019-6) ----
    diabetes_w1 = clean_condition(a_hcond14, a_hconds14),
    
    ### xii. Waist-Hip Index (Christakoudi et al. 2021; 10.1038/s41598-021-89176-6) ----
    across(matches("(b|c)_wstval"), negative_to_na),
    waist_nv = coalesce(b_wstval, c_wstval),
    # TODO: Height ratio 
    
  ) %>%
  select(id, matches("_(w1|nv)$")) %>%
  count(cvd_any_w1)


## d. Mental Health and Cognition ----
clean$psych <- x

df_raw %>%
  mutate(
    ### i. ADHD (Demontis et al., 2023; 10.1038/s41588-022-01285-8) ----
    ### ii. Alzheimer's Disease (Bellenguez et al., 2022; 10.1038/s41588-022-01024-z) ----
    ### iii. Panic Disorder (Forstner et al., 2021; 10.1038/s41380-019-0590-2) ----
    across(matches("^(b|c)_scghq[a-l]$"), negative_to_na),
    ghq_likert_w2 = rowSums(pick(matches("^b_scghq[a-l]$"))),
    ghq_likert_w3 = rowSums(pick(matches("^c_scghq[a-l]$"))),
    ghq_likert_xwave = coalesce(ghq_likert_w2, ghq_likert_w3),
    
    ### iv. Autism Spectrum Disorder (Grove et al., 2019; 10.1038/s41588-019-0344-8) ----
    ### v. Bipolar Disorder (Mullins et al., 2021; 10.1038/s41588-021-00857-4) ----
    ### vi. Intelligence (Savage et al., 2018; 10.1038/s41588-018-0152-6) ----
    across(matches("c_cg"), negative_to_na),
    cog_immed_w3 = c_cgwri_dv,
    cog_delay_w3 = c_cgwrd_dv,
    cog_subtract_w3 = c_cgs7cs_dv,
    cog_verbal_w3 = c_cgvfc_dv,
    cog_numeric_w3 = c_cgna_dv,
    cog_number_w3 = coalesce(c_cgns1sc10_dv, c_cgns2sc10_dv),
    across(matches("^cog_.*_w3$"), ~ scale(.x) %>% as.double()),
    cog_pca_w3 = extract_pca(pick(matches("^cog_.*_w3$"))),
    
    ### v. Hippocampal Volume (Liu et al., 2023; 10.1038/s41588-023-01425-8) ----
    ### vi. Major Depressive Disorder (Howard et al., 2019; 10.1038/s41593-018-0326-7) ----
    clinical_depression_w1 = clean_condition(a_hcond17, a_hconds17),
    
    ### vii. Parkinson's Disease (Nalls et al., 2019; 10.1016/S1474-4422(19)30320-5) ----
    ### viii. Schizophrenia (Trubetskoy et al., 2022; 10.1038/s41586-022-04434-5) ----
    
    ## e. Personality ----
    ### i. Big 5 Personality Traits (Gupta et al., 2024; 10.1038/s41562-024-01951-3) ----
    across(matches("c_big5(o|c|e|a|n)_dv"), negative_to_na),
    big5_openness_w3 = c_big5o_dv,
    big5_conscientiousness_w3 = c_big5c_dv,
    big5_extraversion_w3 = c_big5e_dv,
    big5_agreeableness_w3 = c_big5a_dv,
    big5_neuroticism_w3 = c_big5n_dv,
  ) %>%
  select(id, matches("_(w[1-3]|nv|xwave)$"))

## f. Social Outcomes ----
### i. Educational Attainmment (EA4) (Okbay et al., 2022; 10.1038/s41588-022-01016-z) ----
### ii. Educational Attainmment (EA2) (Okbay et al., 2016; 10.1038/nature17671) ----
clean$edu <- df_raw %>%
  select(id, matches("hiqual_dv")) %>%
  mutate(across(matches("hiqual_dv"), ~ negative_to_na(.x) %>% as.integer())) %>%
  pivot_longer(-id) %>%
  mutate(wave = case_when(str_detect(name, "^b[a-z]_") ~ match(str_sub(name, 2, 2), letters),
                          str_detect(name, "^[a-z]_") ~ match(str_sub(name, 1, 1), letters) + 18)) %>%
  drop_na() %>%
  group_by(id) %>%
  slice_max(wave) %>%
  ungroup() %>%
  mutate(hiqual_xwave = factor(value, c(1:5, 9), c("Degree", "Other HE", "FE", "GCSE etc.", "Other", "None"))) %>%
  select(id, hiqual_xwave)

clean$sep <- df_raw %>%
  mutate(
    ### iii. Household Income (Hill et al., 2019; 10.1038/s41467-019-13585-5) ----
    across(matches("^(a|c)_fi.*net_dv"), negative_to_na),
    household_income_xwave = coalesce(a_fihhmnnet1_dv, b_fihhmnnet1_dv, c_fihhmnnet1_dv),
    personal_income_xwave = coalesce(a_fimnnet_dv, b_fimnnet_dv, c_fimnnet_dv),
    
    ### iv. Occupational Status (Akimova et al., 2024; 10.1038/s41562-024-02076-3) ----
    # TODO: Data request
  ) %>%
  select(id, matches("_xwave$"))

## g. Extra Information ----
clean$cov <- df_raw %>%
  mutate(sex_xwave = negative_to_na(sex_dv) %>% as_factor(),
         dob_y_xwave = as.integet(birthy),
         # TODO: Data request for age.
         #age_nv = 
           ) %>%
  select(id, matches("_(xwave|nw)^"))
