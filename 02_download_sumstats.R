library(tidyverse)
library(vroom)
library(glue)
library(curl)
library(tools)
library(R.utils)
library(fs)
library(magrittr)
library(googledrive)

rm(list = ls())

# 1. Set Up Environment ----
check_md5 <- function(file_stub, md5, file_ext = "txt.gz"){
  file_path <- glue("Source Data/{file_stub}.{file_ext}")
  
  md5_check <- ifelse(md5sum(file_path) == md5, "MD5 Check Passed", "MD5 Check Failed") %>%
    set_names(NULL)
  return(md5_check)
}

download_curl <- function(url, file_stub, overwrite = FALSE, file_ext = "txt.gz"){
  file_path <- glue("Source Data/{file_stub}.{file_ext}")
  
  if (overwrite == TRUE | !file.exists(file_path)){
    curl_download(url, destfile = file_path)
    ret <- "Downloaded"
  } else{
    ret <- "File already exists"
  }
  return(ret)
}

download_pgs_catalog <- function(pgs_id, gwas, build = 38, overwrite = FALSE){
  url <- glue("https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/{pgs_id}/ScoringFiles/Harmonized/{pgs_id}_hmPOS_GRCh{build}.txt.gz")
  file_stub <- glue("{gwas}_grch{build}")
  
  download_curl(url, file_stub)
  
  md5 <- glue("{url}.md5") %>%
    read_lines() %>%
    str_extract("^[^ ]+")
  
  md5_check <- check_md5(file_stub, md5)
  return(md5_check)
}


download_michigan <- function(file_id, file_stub, overwrite = FALSE){
  url <- glue("https://conservancy.umn.edu/bitstreams/{file_id}/download")
  
  download_curl(url, file_stub, overwrite = overwrite)
}

download_figshare <- function(file_id, file_stub, overwrite = FALSE){
  url <- glue("https://figshare.com/ndownloader/files/{file_id}")
  download_curl(url, file_stub, overwrite = overwrite)
}

download_zip <- function(url, file_stub, overwrite = FALSE){
  download_curl(url, file_stub, overwrite = overwrite, file_ext = "zip")
  
  file_path <- glue("Source Data/{file_stub}.zip")
  unzip(file_path, exdir = glue("Source Data/{file_stub}"))
}

download_gwas_catalog <- function(id, file_stub, build = 37, md5sum = NULL, url_ext = "tsv.gz"){
  ret <- "No MD5 Check"
  
  id_floor <- plyr::round_any(id, 1000, floor) + 1
  id_ceiling <- plyr::round_any(id, 1000, ceiling)
  url <- glue("http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST{id_floor}-GCST{id_ceiling}/GCST{id}/GCST{id}_buildGRCh{build}.{url_ext}")
  
  if (!str_detect(url_ext, "\\.gz$")){
    download_curl(url, file_stub, file_ext = "txt")
    if (!is.null(md5sum)){
      ret <- check_md5(file_stub, md5sum, "txt")
    }
    g_zip(glue("{file_stub}.txt"), file_stub)
    remove_excess_files(file_stub)
  } else{
    download_curl(url, file_stub)
    if (!is.null(md5sum)){
      ret <- check_md5(file_stub, md5sum)
    }
  }
  
  return(ret)
  
}

remove_excess_files <- function(file_stub){
  dir_ls("Source Data", regexp = file_stub) %>%
    str_subset("\\.gz$", negate = TRUE) %>%
    file_delete()
}

download_gwas_catalog <- function(file_stub, gzipped = TRUE){
  
}

g_zip <- function(file_path, file_stub, overwrite = FALSE, remove = FALSE){
  gzip(glue("Source Data/{file_path}"),
       glue("Source Data/{file_stub}.txt.gz"),
       overwrite = overwrite, remove = remove)
}

vroom_head <- function(file_stub){
  file_path <- glue("Source Data/{file_stub}.txt.gz")
  vroom(file_path, n_max = 50) %>%
    print(n = 50)
}

add_grch <- function(file_stub, build){
  replace_string <- glue("{file_stub}_grch{build}")
  
  tibble(path = dir_ls("Source Data", regexp = file_stub)) %>%
    mutate(new_path = str_replace(path, file_stub, replace_string)) %$%
    walk2(path, new_path, ~ file_move(.x, .y))
}


# 2. Download Summary Statistics ----
# TODO: Add every latest one for a domain from this: https://pgc.unc.edu/for-researchers/download-results/
# TODO: Check whether I should use harmonised sub-folder in GWAS catalogue
# TODO: Add Neale Lab GWAS alernative where overlap

## a. Health Behaviour ----
### i. Addictive Behaviour (Hatoum et al., 2023; 10.1038/s44220-023-00034-y) ----
download_pgs_catalog("PGS003849", "hatoum_etal_2023_addictive_behaviour")
# TODO: Add Levy et al. (2023) Cannabis Use Disorder; download separate Hatoum PGS

### ii. Alcohol and Cigarette Consumption (Liu et al., 2019; 10.1038/s41588-018-0307-5) ----
download_michigan("2829267c-7663-41b1-b52e-c8c1cc1c8641", "liu_etal_2019_age_smoking_initiation_grch37")
download_michigan("c8ea1cce-f266-4f7a-93cd-87ed9e1afe68/", "liu_etal_2019_cigs_perday_grch37")
download_michigan("6b37ee0d-6eb9-4d66-ac12-02092ce2731b", "liu_etal_2019_drinks_perday_grch37")
download_michigan("0ccb7fe4-6092-445f-b44d-ff33cd94c7ee", "liu_etal_2019_smoking_cessation_grch37")
download_michigan("1359186c-fd6a-4a39-a31f-640ed52c5efe", "liu_etal_2019_smoking_initiation_grch37")

### iii. Diet (Cole et al., 2020; 10.1038/s41467-020-15193-0) ----
# TODO: Decide which to download


## b. Anthropometrics ----
### i. Age at Menarche (Day et al., 2017; 10.1038/ng.3841) ----
download_zip("https://www.reprogen.org/Menarche_1KG_NatGen2017_WebsiteUpload.zip",
             "day_etal_2017_age_menarche_grch37")
gzip(glue("Source Data/day_etal_2017_age_menarche_grch37/Menarche_1KG_NatGen2017_WebsiteUpload.txt"),
     glue("Source Data/day_etal_2017_age_menarche_grch37.txt.gz"),
     overwrite = FALSE, remove = FALSE)
remove_excess_files("day_etal_2017_age_menarche_grch37")

### ii. Age at Menopause (Ruth et al., 2021; 10.1038/s41586-021-03779-7) ----
download_curl("https://www.reprogen.org/reprogen_ANM_201K_170621.txt.gz",
              "ruth_etal_2021_age_menopause_grch37")

### iii Birth Weight (Warrington et al., 2019; 10.1038/s41588-019-0403-1) ----
# "1. European-only meta-analysis of own birth weight in up to 298,142 individuals; the data as a gzipped text file can be downloaded here."
# TODO: 8 options - could also be trans-ancestry or offspring not own birthweight and partitioning out maternal effect
download_curl("http://egg-consortium.org/BW5/Fetal_BW_European_meta.NG2019.txt.gz",
              "warrington_etal_2019_birthweight_net_grch37")
download_curl("http://egg-consortium.org/BW5/Fetal_Effect_European_meta_NG2019.txt.gz",
              "warrington_etal_2019_birthweight_direct_grch37")
# download_curl("http://egg-consortium.org/BW5/Maternal_Effect_European_meta_NG2019.txt.gz",
#               "warrington_etal_2019_birthweight_indirect")

### iv. Body Fat Percentage (Roshandel et al., 2023; https://doi.org/10.3389/fendo.2023.1274791) ----
# TODO: Add Neale Lab because this has controls.
download_curl("http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90432001-GCST90433000/GCST90432180/GCST90432180.tsv",
              "roshandel_etal_2023_body_fat_female_grch37", 
              file_ext = "txt")
check_md5("roshandel_etal_2023_body_fat_female_grch37",
          "6a3a9e6e73b732bce387726fb347c2b3",
          "txt")
g_zip("roshandel_etal_2023_body_fat_female_grch37.txt",
      "roshandel_etal_2023_body_fat_female_grch37")
remove_excess_files("roshandel_etal_2023_body_fat_female_grch37")

download_curl("http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90432001-GCST90433000/GCST90432179/GCST90432179.tsv",
              "roshandel_etal_2023_body_fat_male_grch37", 
              file_ext = "txt")
check_md5("roshandel_etal_2023_body_fat_male_grch37",
          "8bccd275913fe2b70d916cdc6d8adc00",
          "txt")
g_zip("roshandel_etal_2023_body_fat_male_grch37.txt",
      "roshandel_etal_2023_body_fat_male_grch37")
remove_excess_files("roshandel_etal_2023_body_fat_male_grch37")

### v. (Adult) Body Mass Index (Yengo et al., 2018; 10.1093/hmg/ddy271) ----
download_curl("http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006900/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt.gz",
              "yengo_etal_2018_bmi_grch37")

### vi. (Childhood) Body Mass Index (Vogelezang et al. 2020; 10.1371/journal.pgen.1008718) ----
download_curl("http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90002001-GCST90003000/GCST90002409/GCST90002409_buildGRCh37.tsv.gz",
              "vogelezang_etal_2020_bmi_grch37")
check_md5("vogelezang_etal_2020_bmi_grch37", "09a140e7e1e62a81e7eda694e04ee296")

### vii. Height (Yengo et al., 2018; 10.1093/hmg/ddy271) ----
## TODO: Add Yengo et al., 2022.
download_curl("http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006901/Meta-analysis_Wood_et_al+UKBiobank_2018.txt.gz",
              "yengo_etal_2018_height_grch37")

### viii. Height (Yengo et al., 2022; 10.1038/s41586-022-05275-y) ----
### TODO: Check (as in this case) what to do when there is a distinction between GWAS and PGS weights
download_curl("https://portals.broadinstitute.org/collaboration/giant/images/4/4e/GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_ALL.gz",
              "yengo_etal_2022_height_all_grch37")

download_curl("https://portals.broadinstitute.org/collaboration/giant/images/f/f7/GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EUR.gz",
              "yengo_etal_2022_height_european_grch37")

### ix. Grip Strength (Jones et al., 2021; 10.1038/s41467-021-20918-w) ----
download_gwas_catalog(90007526, 
                      "jones_etal_2021_grip_strength_ewgsop_grch37", 
                      md5sum = "10ed84976ef644b1ac7ceb06cd34e2dc")

download_gwas_catalog(90007529, 
                      "jones_etal_2021_grip_strength_fnih_grch37", 
                      md5sum = "8a11d6233a6531c9e7d4bc92f3cffc3f")

### x. Grip Strength (Tikkanen et al., 2018; 10.1038/s41598-018-24735-y) ----
# TODO: This, but LD Hub website where sumstats are stored isn't loading.

## c. Physical Health and Biomarkers ----
### i Asthma (Han et al., 2020; org/10.1038/s41467-020-15649-3) ----
download_curl("http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010043/HanY_prePMID_asthma_Meta-analysis_UKBB_TAGC.txt.gz",
              "han_etal_2020_asthma_grch37")
check_md5("han_etal_2020_asthma_grch37", "e2da63d044ad6dd5479e1618a7ea5872")

### ii. Blood Pressure (Keaton et al., 2024; 10.1038/s41588-024-01714-w) ----
download_curl("http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90310001-GCST90311000/GCST90310294/GCST90310294.tsv.gz",
              "keaton_etal_2024_sbp_grch37")
check_md5("keaton_etal_2024_sbp_grch37", "517099adeddf65fbf0a8a2f820d37bf6")

download_curl("http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90310001-GCST90311000/GCST90310295/GCST90310295.tsv.gz",
              "keaton_etal_2024_dbp_grch37")
check_md5("keaton_etal_2024_dbp_grch37", "7285a468f22e912fb4e13f1fdc43e834")

### iii. C-Reactive Protein (Said et al., 2022; 10.1038/s41467-022-29650-5) ----
download_curl("http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90029001-GCST90030000/GCST90029070/GCST90029070_buildGRCh37.tsv.gz",
              "said_etal_2022_crp_grch37")
check_md5("said_etal_2022_crp_grch37", "d07a10dc5d032313de677c74b4476eeb")

### iv. Coronary Artery Disease (Aragam et al., 2022; 10.1038/s41588-022-01233-6) ----
download_curl("http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90132001-GCST90133000/GCST90132314/GCST90132314_buildGRCh37.tsv",
              "aragam_etal_2022_coronary_artery_disease_european_grch37", 
              file_ext = "txt")
check_md5("aragam_etal_2022_coronary_artery_disease_european_grch37",
          "92bc58de063ebbf34c4ee39adbdaf798",
          "txt")
g_zip("aragam_etal_2022_coronary_artery_disease_european_grch37.txt",
      "aragam_etal_2022_coronary_artery_disease_european_grch37")
remove_excess_files("aragam_etal_2022_coronary_artery_disease_european_grch37")

download_curl("http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90132001-GCST90133000/GCST90132315/GCST90132315_buildGRCh37.tsv",
              "aragam_etal_2022_coronary_artery_disease_combined_grch37", 
              file_ext = "txt")
check_md5("aragam_etal_2022_coronary_artery_disease_combined_grch37",
          "dba8af649686ec1632b12ad5b1ffc4ea",
          "txt")
g_zip("aragam_etal_2022_coronary_artery_disease_combined_grch37.txt",
      "aragam_etal_2022_coronary_artery_disease_combined_grch37")
remove_excess_files("aragam_etal_2022_coronary_artery_disease_combined_grch37")

### v. Fasting Blood Glucose and Hba1c (Downie et al., 2022; 10.1007/s00125-021-05635-9) ----
download_gwas_catalog(90094959, 
                      "downie_etal_2022_fasting_glucose_combined_grch37", 
                      md5sum = "ec6f2ad3403c4669b7ca00f1d35918b5", 
                      url_ext = "tsv")

download_gwas_catalog(90094960, 
                      "downie_etal_2022_fasting_glucose_european_grch37", 
                      md5sum = "b9f689302cf0d54bc023c24f3509d3e3")

download_gwas_catalog(90094966,
                      "downie_etal_2022_hba1c_combined_grch37", 
                      md5sum = "5d3f679e884eede7834f3580bead533f")

download_gwas_catalog(90094967, 
                      "downie_etal_2022_hba1c_european_grch37", 
                      md5sum = "550193ff1885707766a521ee83a75dc0")

### vi. Hba1c (Sinnott-Armstrong et al., 2021; 10.1038/s41588-020-00757-z) ----
download_figshare(22770197, "sinnotarmstrong_etal_2021_hba1c_grch37")

### vii. Parental Lifespan (Timmers et al., 2019; 10.7554/eLife.39856) ----
download_curl("https://datashare.ed.ac.uk/bitstream/handle/10283/3209/lifegen_phase2_bothpl_alldr_2017_09_18.tsv.gz",
              "timmers_etal_2019_parental_lifespan_grch37")

### ix. Hypertension (Bi et al., 2020; 10.1016/j.ajhg.2020.06.003) ----
download_curl("https://share.sph.umich.edu/UKBB_SPACox_HRC/Summary-Statistics-UKBB/X401.1.txt",
              "bi_etal_2020_hypertension", 
              file_ext = "txt")
g_zip("bi_etal_2020_hypertension.txt",
      "bi_etal_2020_hypertension")
remove_excess_files("bi_etal_2020_hypertension")

### x. Rheumatoid Arthritis (Ishigaki et al., 2022; 10.1038/s41588-022-01213-w) ----
download_gwas_catalog(90132222,
                      "ishigaki_etal_2022_arthritis_all_combined_grch37",
                      md5sum = "dffee6d90b36d03f129d1d49d305875b")

download_gwas_catalog(90132223,
                      "ishigaki_etal_2022_arthritis_all_european_grch37",
                      md5sum = "6150b0bc3448467e9112a39cfea16616")

download_gwas_catalog(90132225,
                      "ishigaki_etal_2022_arthritis_seropositive_combined_grch37",
                      md5sum = "350b9002dd015e7b5fa9b0a5dd7de300")

download_gwas_catalog(90132226,
                      "ishigaki_etal_2022_arthritis_seropositive_european_grch37",
                      md5sum = "95a30fae97a56389413027e16fbd4c4e")

### xi. Type I Diabetes (Chiou et al., 2021; 10.1038/s41586-021-03552) ----
download_gwas_catalog(90014023,
                      "chiou_etal_2021_t1d_grch38",
                      38,
                      "ea8d72e8261f10d70e62242c9dc8ae3c",
                      "tsv")

### xi. Type II Diabetes (Suzuki et al., 2024; 10.1038/s41586-024-07019-6) ----
# Downloaded from http://www.diagram-consortium.org/downloads.html, 2025-02-10
# Suzuki.Nature2024.T2DGGI.EUR.sumstats.zip; # Suzuki.Nature2024.T2DGGI.MultiAncestry.sumstats.zip
check_md5("Downloaded/Suzuki.Nature2024.T2DGGI.EUR.sumstats", "bc7402c2c68b94be1ac8b42fa7e05431", "zip") # Failed md5!
unzip(glue("Source Data/Downloaded/Suzuki.Nature2024.T2DGGI.EUR.sumstats.zip"), 
      exdir = glue("Source Data/suzuki_etal_2024_t2d_european_grch37"))
gzip(glue("Source Data/suzuki_etal_2024_t2d_european_grch37/EUR_Metal_LDSC-CORR_Neff.v2.txt"),
     glue("Source Data/suzuki_etal_2024_t2d_european_grch37.txt.gz"),
     overwrite = FALSE, remove = FALSE)
remove_excess_files("suzuki_etal_2024_t2d_european_grch37")

check_md5("Downloaded/Suzuki.Nature2024.T2DGGI.MultiAncestry.sumstats", 
          "21503cf0e292e5c2c14c51d74b0aeeea", "zip") # Failed md5!
unzip(glue("Source Data/Downloaded/Suzuki.Nature2024.T2DGGI.MultiAncestry.sumstats.zip"), 
      exdir = glue("Source Data/suzuki_etal_2024_t2d_combined_grch37"))
gzip(glue("Source Data/suzuki_etal_2024_t2d_combined_grch37/All_Metal_LDSC-CORR_Neff.v2.txt"),
     glue("Source Data/suzuki_etal_2024_t2d_combined_grch37.txt.gz"),
     overwrite = FALSE, remove = FALSE)
remove_excess_files("suzuki_etal_2024_t2d_combined_grch37")

### xii. Waist-Hip Index (Christakoudi et al. 2021; 10.1038/s41598-021-89176-6) ----
download_curl("http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST012001-GCST013000/GCST012230/GCST012230.tsv.gz",
              "christakoudi_etal_2021_waist_hip_index_bmiadj_grch37")
check_md5("christakoudi_etal_2021_waist_hip_index_bmiadj_grch37",
          "711313d3cfea303aef5aff5595e8ae01")

download_curl("http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST012001-GCST013000/GCST012228/GCST012228.tsv.gz",
              "christakoudi_etal_2021_waist_hip_index_noadj_grch37")
check_md5("christakoudi_etal_2021_waist_hip_index_noadj_grch37",
          "0066820f13efb4dde54715ccba559acd")

## d. Mental Health and Cognition ----
### i. ADHD (Demontis et al., 2023; 10.1038/s41588-022-01285-8) ----
download_figshare(40036684, "demontis_etal_2023_adhd_grch37")

### ii. Alzheimer's Disease (Bellenguez et al., 2022; 10.1038/s41588-022-01024-z) ----
download_curl("http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90027001-GCST90028000/GCST90027158/GCST90027158_buildGRCh38.tsv.gz",
              "bellenguez_etal_2022_alzheimers_grch38")
check_md5("bellenguez_etal_2022_alzheimers_grch38", "9d23b9ba23532da38ab83fb061bab18f")

### iii. Panic Disorder (Forstner et al., 2021; 10.1038/s41380-019-0590-2) ----
download_curl(30731276, "forstner_etal_2021_panic_disorder_grch37")

### iv. Autism Spectrum Disorder (Grove et al., 2019; 10.1038/s41588-019-0344-8) ----
download_curl(28169292, "grove_etal_2019_autism_grch37")

### v. Bipolar Disorder (Mullins et al., 2021; 10.1038/s41588-021-00857-4) ----
# TODO: Add medRxiv 2024 paper from same consortium?
download_figshare(26603681, "mullins_etal_2021_bipolar_all_grch37")
download_figshare(26603690, "mullins_etal_2021_bipolar_i_grch37")
download_figshare(26603702, "mullins_etal_2021_bipolar_ii_grch37")

### vi. Intelligence (Savage et al., 2018; 10.1038/s41588-018-0152-6) ----
download_curl("http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006250/harmonised/29942086-GCST006250-EFO_0004337-Build37.f.tsv.gz",
              "savage_etal_2018_intelligence_grch37")

### v. Hippocampal Volume (Liu et al., 2023; 10.1038/s41588-023-01425-8) ----
# TODO: Decide which of the 132 scores to include

### vi. Major Depressive Disorder (Howard et al., 2019; 10.1038/s41593-018-0326-7) ----
# TODO: Add Cell 2025 paper from same consortium
download_curl("https://datashare.ed.ac.uk/bitstream/handle/10283/3203/PGC_UKB_23andMe_depression_10000.txt",
              "howard_etal_2019_depression_10K_rsid", file_ext = "txt")
gzip(glue("Source Data/howard_etal_2019_depression_10K_rsid.txt"),
     glue("Source Data/howard_etal_2019_depression_10K_rsid.txt.gz"),
     overwrite = FALSE, remove = FALSE)
remove_excess_files("howard_etal_2019_depression_10K_rsid")

download_curl("https://datashare.ed.ac.uk/bitstream/handle/10283/3203/PGC_UKB_depression_genome-wide.txt",
              "howard_etal_2019_depression_all_rsid", file_ext = "txt")
gzip(glue("Source Data/howard_etal_2019_depression_all_rsid.txt"),
     glue("Source Data/howard_etal_2019_depression_all_rsid.txt.gz"),
     overwrite = FALSE, remove = FALSE)
remove_excess_files("howard_etal_2019_depression_all_rsid")

### vii. Parkinson's Disease (Nalls et al., 2019; 10.1016/S1474-4422(19)30320-5) ----
drive_download(as_id("1FZ9UL99LAqyWnyNBxxlx6qOUlfAnublN"),
               glue("Source Data/nalls_etal_2019_parkinsons.zip"))
unzip(glue("Source Data/nalls_etal_2019_parkinsons.zip"), 
      exdir = glue("Source Data/nalls_etal_2019_parkinsons"))
gzip(glue("Source Data/nalls_etal_2019_parkinsons/nallsEtAl2019_excluding23andMe_allVariants.tab"),
     glue("Source Data/nalls_etal_2019_parkinsons.txt.gz"),
     overwrite = FALSE, remove = FALSE)
remove_excess_files("nalls_etal_2019_parkinsons")

### vii. Schizophrenia (Trubetskoy et al., 2022; 10.1038/s41586-022-04434-5) ----
download_figshare(34517807, "trubetskov_etal_2022_schizophrenia_core_grch37")
download_figshare(34517861, "trubetskov_etal_2022_schizophrenia_primary_grch37")
download_figshare(34517828, "trubetskov_etal_2022_schizophrenia_european_grch37")

## e. Personality ----
### i. Big 5 Personality Traits (Gupta et al., 2024; 10.1038/s41562-024-01951-3) ----
download_curl("https://hugesumstats.yale.edu/dl/Neuro_MVP_UKB_sumstat_file",
              "gupta_etal_2024_neuroticism_grch37")
download_curl("https://hugesumstats.yale.edu/dl/agree_MVP_GPC1_sumstat_file",
              "gupta_etal_2024_agreeableness_grch37")
download_curl("https://hugesumstats.yale.edu/dl/consc_MVP_GPC1_sumstat_file",
              "gupta_etal_2024_conscientiousness_grch37")
download_curl("https://hugesumstats.yale.edu/dl/open_MVP_GPC1_sumstat_file",
              "gupta_etal_2024_openness_grch37")
download_curl("https://hugesumstats.yale.edu/dl/extra_MVP_GPC1_sumstat_file",
              "gupta_etal_2024_extraversion_grch37")

## f. Social Outcomes ----
### i. Educational Attainmment (EA4) (Okbay et al., 2022; 10.1038/s41588-022-01016-z) ----
# Downloaded from SSGAC website 2025-02-10
# TODO: Find ones sent to me.
g_zip("Downloaded/EA4_additive_p1e-5_clumped.txt", 
      "okbay_etal_2022_eduyears_top_inc23andme_grch37")

file_copy(glue("Source Data/Downloaded/EA4_additive_excl_23andMe.txt.gz"),
          glue("Source Data/okbay_etal_2022_eduyears_all_excl3andme_grch37.txt"))

g_zip("Downloaded/EA4_excl_UKHLS_2024_12_09_publicSNPs.txt", 
      "okbay_etal_2022_eduyears_top_exclukhls_grch37")

### ii. Educational Attainmment (EA2) (Okbay et al., 2016; 10.1038/nature17671) ----
# Downloaded from SSGAC website 2025-02-10
g_zip("Downloaded/EduYears_Main.txt", 
      "okbay_etal_2016_eduyears_all_excl23andme_grch37")

g_zip("Downloaded/EduYears_Pooled_5000.txt", 
      "okbay_etal_2016_eduyears_top5k_inc23andme_grch37")


### iii. Household Income (Hill et al., 2019; 10.1038/s41467-019-13585-5) ----
download_curl("https://datashare.ed.ac.uk/bitstream/handle/10283/3441/MTAG_household_Income.txt.gz",
              "hill_etal_2019_income_mtag_grch37")

download_curl("https://datashare.ed.ac.uk/bitstream/handle/10283/3441/Household_Income_UKBiobank.txt.gz",
              "hill_etal_2019_income_ukb_grch37")

### iv. Occupational Status (Akimova et al., 2024; 10.1038/s41562-024-02076-3) ----
download_akimova_2024 <- function(id, measure, md5sum){
  file_stub <- glue("akimova_etal_2024_ses_{measure}_grch37")
  url <- glue("http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90446001-GCST90447000/GCST{id}/GCST{id}.tsv")
  
  download_curl(url, file_stub, file_ext = "tsv")
  ret <- check_md5(file_stub, md5sum, "tsv")
  g_zip(glue("{file_stub}.tsv"), file_stub)
  remove_excess_files(file_stub)
  
  return(ret)
}

download_akimova_2024(90446160, "isei", "39aebeecdb040965bc0059a94762a7c8")
download_akimova_2024(90446162, "camsis", "af0315dda9842cacc03b383f37c85d2c")
download_akimova_2024(90446163, "siops", "d5f6dae95ae15e34dfb504c67147870c")

