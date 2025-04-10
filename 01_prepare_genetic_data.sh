#!/bin/bash -l
#$ -m be
#$ -N prepare_bed
#$ -l h_rt=24:00:00
#$ -l mem=8G
#$ -pe smp 8
#$ -wd "/myriadfs/home/rmjllwr/Projects/UKHLS Fellowship/Polygenic Scores Bank"

# qrsh -pe smp 8 -l mem=8G,h_rt=2:00:00 -now n
# qsub "Code/01_prepare_genetic_data.sh"

# 0. Set up Environment -----
cd "/myriadfs/home/rmjllwr/Projects/UKHLS Fellowship/Polygenic Scores Bank"
source config.cfg

module -f unload compilers mpi gcc-libs
module load r/recommended
module load lapack/3.8.0/gnu-4.9.2

mkdir -p Data/step_01

# 1. Amend Sample File ----
cp "${sample_file}" Data/step_01/original.sample

cat > Code/01a_amend_sample_file.R << EOM
library(tidyverse)

old_sample <- read_tsv("Data/step_01/original.sample") %>%
  rename_with(str_to_lower) %>%
  mutate(id1_int = as.integer(id1),
          new_id = ifelse(!is.na(id1_int), id1_int, max(id1_int, na.rm = TRUE) + sample_order),
          sex = as.character(sex)) 
          
old_sample %>%
  select(ID_1 = new_id, ID_2 = new_id, missing, sex) %>%
  add_row(ID_1 = 0, ID_2 = 0, missing = 0, sex = "D",
          .before = 1) %>%
  write_tsv("Data/step_01/amended.sample")

old_sample %>%
  filter(is.na(id1_int)) %>%
  select(fid = new_id, iid = new_id) %>%
  write_tsv("Data/step_01/extraid.id", col_names = FALSE)
EOM

Rscript  --no-save --no-restore Code/01a_amend_sample_file.R

# 2. Drop High Missingness Variants and Extra IDs and Resave to Remove Trailing Space ----
for i in {1..22}; do
    "${plink2}"   --gen "${geno_dir}/chr${i}.UKHLS.UK10K+1KG-imputed.phwe_1e-4.info_0.4.filtered.gen.gz" ref-unknown \
                  --sample "Data/step_01/amended.sample" \
                  --oxford-single-chr ${i} \
                  --remove "Data/step_01/extraid.id" \
                  --geno 0.03 \
                  --memory 8000 --threads 8 \
                  --recode oxford bgz \
                  --out Data/step_01/chr${i}_geno
done

# 3. Filter on INFO Scores ----
## a. Calculate INFO Scores ----
for i in {1..22}; do
    "${qctool}" -g Data/step_01/chr${i}_geno.gen.gz \
                -snp-stats \
                -osnp Data/step_01/chr${i}_snp_stats.txt
done

## b. Find Variants to Keep ----
cat > Code/01b_filter_info.R << EOM
library(tidyverse)
library(glue)
library(vroom)

rm(list = ls())

clean_snp_info <- function(i){
  print(i)
  df <- glue("Data/step_01/chr{i}_snp_stats.txt") %>%
    vroom(delim = " ", 
          skip = 8,
          col_select = c(rsid, info, impute_info)) %>%
    filter(info >= 0.8) %>%
    select(rsid)
  write_delim(df, glue("Data/step_01/chr{i}_info.txt"), col_names = FALSE)
  return(i)
}

walk(1:22, clean_snp_info)
EOM

Rscript  --no-save --no-restore Code/01b_filter_info.R

## c. Filter on INFO Scores and Convert to .bed ----
rm Data/step_01/merge_info.txt
for i in {1..22}; do
    "${plink2}" --gen "Data/step_01/chr${i}_geno.gen.gz" ref-unknown \
                --sample "Data/step_01/chr${i}_geno.sample" \
                --oxford-single-chr ${i} \
                --extract Data/step_01/chr${i}_info.txt \
                --rm-dup exclude-mismatch \
                --max-alleles 2 \
                --memory 8000 --threads 8 \
                --make-bed \
                --out Data/step_01/chr${i}_info

  echo Data/step_01/chr${i}_info >> Data/step_01/merge_info.txt
  # rm Data/step_01/chr${i}_geno.{gen.gz,sample}
done


## 4. Update .fam files and find unrelateds using KING ----
## a. Combine .bed Files and Temporarily Filter on MAF ----
"${plink2}"   --pmerge-list Data/step_01/merge_info.txt bfile \
              --rm-dup exclude-mismatch \
              --maf 0.01 \
              --maj-ref \
              --sort-vars \
              --make-bed \
              --out Data/step_01/ukhls_rsid_grch37_maf_temp
              
rm Data/step_01/ukhls_rsid_grch37_maf_temp.{pvar,pgen,psam}

md5sum  chr{1..22}_info.fam ukhls_rsid_grch37_maf_temp.fam \
        | awk '{print $1}' \
        | uniq -c

## b. Create Relationship Table using KING ----          
"${king}" -b "Data/step_01/ukhls_rsid_grch37_maf_temp.bed" \
          --cpus 8 \
          --related --degree 3 \
          --prefix "Data/step_01/king_relatedness"

## c. Create new .fam Files ----
cat > Code/01c_update_fam.R << EOM
library(tidyverse)
library(haven)
library(igraph)

rm(list = ls())

args <- commandArgs(trailingOnly = TRUE)
pheno_file <- args[1]

# 1. Load Data ----
pheno <- read_dta(pheno_file) %>%
  mutate(sex = as_factor(sex) %>% fct_drop(),
         dob_y = as.integer(doby_dv)) %>%
  select(iid = id, sex, dob_y)

king <- read_delim("Data/step_01/king_relatedness.kin0",
                   delim = "\t") %>%
  select(iid.x = ID1, iid.y = ID2, rel = InfType) %>%
  mutate(row = row_number()) %>%
  uncount(2, .id = "id") %>%
  group_by(row) %>%
  mutate(iid.x = ifelse(id == 1, iid.x, lag(iid.y)),
         iid.y = ifelse(id == 1, iid.y, lag(iid.x))) %>%
  ungroup() %>%
  select(-id)

# 2. Find Parents (Based on KING Relatedness) ----
fam <- read_delim("Data/step_01/ukhls_rsid_grch37_maf_temp.fam", 
                  delim = "\t", 
                  col_names = c("fid", "iid", "pid", "mid", "sex", "pheno"))

king_chr <- king %>%
  left_join(pheno %>% select(iid, dob_y), by = c("iid.x" = "iid")) %>%
  left_join(pheno %>% select(iid, dob_y), by = c("iid.y" = "iid")) %>%
  left_join(fam %>% select(iid, sex), by = c("iid.y" = "iid")) %>%
  left_join(fam %>% select(iid, sex), by = c("iid.y" = "iid"))

df_parents <- king_chr %>%
  filter(rel == "PO",
         dob_y.x > dob_y.y) %>%
  mutate(type = case_when(sex.y == 1 ~ "pid", 
                          sex.y == 2 ~ "mid")) %>%
  drop_na(type) %>%
  add_count(iid.x, type) %>%
  filter(n == 1) %>%
  pivot_wider(id_cols = iid.x,
              names_from = type,
              values_from = iid.y) %>%
  rename(iid = iid.x)

## b. Unit Tests ----
check_rel <- function(rel){
  king %>%  # UNIT TEST: Siblings have same parents.
    filter(rel == !!rel) %>%
    left_join(df_parents, by = c("iid.x" = "iid")) %>%
    left_join(df_parents, by = c("iid.y" = "iid")) %>%
    replace_na(list(mid.x = 0, mid.y = 0, pid.x = 0, pid.y = 0)) %>%
    filter(mid.x != mid.y | pid.x != pid.y)
}
ret <- check_rel("FS") %>% nrow() # UNIT TEST: Siblings have same parents.

check_rel("2nd") %>% 
  filter(iid.x < iid.y) %>%
  select(matches("^(iid|mid|pid)")) %>%
  print(n = Inf)

## 3. Amend .fam ----
king_fam <- king_chr %>%
  filter(rel %in% c("FS", "PO")) %>%
  select(from = iid.x, to = iid.y) %>%
  graph_from_data_frame(directed = FALSE) %>%
  components() %>%
  pluck("membership") %>%
  enframe(name = "iid", value = "fid") %>%
  mutate(iid = as.integer(iid)) %>%
  left_join(fam %>% select(iid, sex), ., by = "iid") %>%
  mutate(fid = ifelse(is.na(fid), max(fid, na.rm = TRUE) + row_number(), fid)) %>%
  left_join(df_parents, by = "iid") %>%
  replace_na(list(mid = 0, pid = 0)) %>%
  mutate(pheno = -9) %>%
  select(fid, iid, pid, mid, sex, pheno)

write_delim(king_fam, "Data/step_01/king.fam",
            delim = " ", col_names = FALSE)

# 4. Find Unrelated Individuals ----
no_relation <- fam %>%
  filter(!(iid %in% king[["iid.x"]]),
         !(iid %in% king[["iid.y"]])) %>%
  pull(iid)

king_graph <- king %>%
  select(from = iid.x, to = iid.y) %>%
  graph_from_data_frame(directed = FALSE) %>%
  add_vertices(nv = length(no_relation), names = no_relation)

set.seed(1)
unrelated <- largest_ivs(king_graph) %>%
  sample(1) %>%
  pluck(1) %>%
  names() %>%
  as.integer()

king_fam %>%
  filter(iid %in% !!unrelated) %>%
  write_tsv("Data/step_01/king_unrelated.id", col_names = FALSE)
EOM

# R --args "${pheno_file}"
Rscript  --no-save --no-restore Code/01c_update_fam.R "${pheno_file}"

# 5. Filter chromosones on MAF ----
for i in {1..22}; do
    "${plink2}" --bed Data/step_01/chr${i}_info.bed \
                --bim Data/step_01/chr${i}_info.bim \
                --fam Data/step_01/king.fam \
                --keep Data/step_01/king_unrelated.id \
                --maf 0.01 \
                --memory 8000 --threads 8 \
                --write-snplist \
                --out Data/step_01/chr${i}_maf

    "${plink2}" --bed Data/step_01/chr${i}_info.bed \
                --bim Data/step_01/chr${i}_info.bed \
                --fam Data/step_01/king.fam \
                --extract Data/step_01/chr${i}_maf.snplist \
                --memory 8000 --threads 8 \
                --make-bed \
                --out Data/step_01/chr${i}_maf
done

# 6. Clean .bim Files ----
## a. Find SNPs to re-ID ----
for i in {1..22}; do
  awk -F'\t' '$2 !~ /^rs[0-9]+$/ || index($2, ";") > 0' OFS='\t' \
    Data/step_01/chr${i}_maf.bim \
    > Data/step_01/chr${i}_toid.txt
done

## b. Use tabix to find matches ----
# mkdir -p Data/step_01/tabix
for i in {1..22}; do
  # wget -nc -P "Data/step_01/tabix" https://ftp.ensembl.org/pub/grch37/current/variation/vcf/homo_sapiens/homo_sapiens-chr${i}.vcf.gz # Downloaded 2025-03-24
  # "${tabix}" -p vcf "Data/step_01/tabix/homo_sapiens-chr${i}.vcf.gz"

  while read chr id dist pos ref alt; do
    "${tabix}" "Data/step_01/tabix/homo_sapiens-chr${i}.vcf.gz" ${chr}:${pos}-${pos}
  done < "Data/step_01/chr${i}_toid.txt" > "Data/step_01/chr${i}_tabix.txt"
done

## c. Re-ID .bim File ----
cat > Code/01d_reid_bim.R << EOM
library(tidyverse)
library(glue)
library(tictoc)

rm(list = ls())

reid_bim <- function(i){
  # 1. Load Data ----
  bim <- glue("Data/step_01/chr{i}_maf.bim") %>%
    read_delim(col_types = "iciicc", delim = "\t",
               col_names = c("chr", "id", "dist", "pos", "ref", "alt"))
  
  to_id <- glue("Data/step_01/chr{i}_toid.txt") %>%
    read_delim(col_types = "iciicc", delim = "\t",
               col_names = c("chr", "id", "dist", "pos", "ref", "alt"))
  
  tabix <- glue("Data/step_01/chr{i}_tabix.txt") %>%
    read_delim(col_types = "iicccccc", delim = "\t", 
               col_names = c("chr", "pos", "rsid", "ref", "alt", "x1", "x2", "info")) %>%
    separate(info, into = c("database", "info"), sep = ";", extra = "merge") %>%
    separate_rows(alt, sep = ",")
  
  # 2. Merge Files ----
  matched <- to_id %>%
    left_join(tabix, by = c("chr", "pos")) %>%
    filter((ref.x == ref.y & alt.x == alt.y) | 
             (ref.x == alt.y & alt.x == ref.y)) %>%
    add_count(id) %>%
    filter(n == 1) %>%
    select(chr, id, pos, rsid)
  prop <- nrow(matched) / nrow(to_id)
  
  new_bim <- bim %>%
    left_join(matched, by = c("chr", "id", "pos")) %>%
    mutate(id = case_when(str_sub(id, 1, 2) == "rs" & !str_detect(id, "\\\;") ~ id,
                          str_sub(id, 1, 2) == "rs" & str_detect(id, "\\\;") & is.na(rsid) ~ id,
                          str_sub(rsid, 1, 2) == "rs" ~ rsid,
                          TRUE ~ glue("chr{chr}:{pos}"))) %>%
    select(-rsid)
  
  # 3. Save New .bim ----
  write_tsv(new_bim, 
            glue("Data/step_01/chr{i}_reid.bim"), 
            col_names = FALSE)

  glue("Chromosome {i}: {nrow(matched)} of {nrow(to_id)} ({round(100*prop, 1)}%)")
}

tic()
map_chr(1:22, reid_bim) %>% glue_collapse("\n")
toc()
EOM

Rscript  --no-save --no-restore Code/01d_reid_bim.R

## d. Identify AGCT SNPs ----
for i in {1..22}; do
  awk '$5 ~ /^[AGCT]+$/ && $6 ~ /^[AGCT]+$/ { print $2 }' \
    Data/step_01/chr${i}_reid.bim \
    > Data/step_01/chr${i}_agct.snplist
done

## e. Extract SNPs and Combine ----
rm Data/step_01/merge_dedup.txt
for i in {1..22}; do
    "${plink2}" --bed Data/step_01/chr${i}_maf.bed \
                --bim Data/step_01/chr${i}_reid.bim \
                --fam Data/step_01/chr${i}_maf.fam \
                --extract Data/step_01/chr${i}_agct.snplist \
                --rm-dup exclude-mismatch \
                --memory 8000 --threads 8 \
                --make-bed \
                --out Data/step_01/chr${i}_dedup

  echo Data/step_01/chr${i}_dedup >> Data/step_01/merge_dedup.txt
done

"${plink2}"   --pmerge-list Data/step_01/merge_dedup.txt bfile \
              --rm-dup exclude-mismatch \
              --maj-ref \
              --sort-vars \
              --make-bed \
              --out Data/step_01/ukhls_rsid_grch37_dedup
              
rm Data/step_01/ukhls_rsid_grch37_dedup.{pvar,pgen,psam}

# 7. Remove High Missingness ----
"${plink2}"   --bfile Data/step_01/ukhls_rsid_grch37_dedup \
              --geno 0.03 \
              --write-snplist \
              --out Data/step_01/ukhls_rsid_grch37_geno

"${plink2}"   --bfile Data/step_01/ukhls_rsid_grch37_dedup \
              --extract Data/step_01/ukhls_rsid_grch37_geno.snplist \
              --mind 0.02 \
              --write-samples \
              --write-snplist \
              --out Data/step_01/ukhls_rsid_grch37_mind

# 8. Remove Upon Harvey-Weinberg Equilibrium ----
# Q: Should this be based on unrelated samples?
awk 'NR==FNR { ids[$2]; next } $2 in ids' \
  Data/step_01/ukhls_rsid_grch37_mind.id \
  Data/step_01/king_unrelated.id \
  > Data/step_01/unrelated_mind.id

"${plink2}"   --bfile Data/step_01/ukhls_rsid_grch37_dedup \
              --extract Data/step_01/ukhls_rsid_grch37_mind.snplist \
              --keep Data/step_01/unrelated_mind.id \
              --hwe 1e-6 \
              --write-snplist \
              --write-samples \
              --out Data/step_01/ukhls_rsid_grch37_hwe


# 9. Filter on Het Scores ----
## Q: Should this be done with unrelated samples?
## TODO: Calculate Mean, SD of Het using unrelated individuals
## a. Calculate LD and Het Scores ----
"${plink2}"   --bfile Data/step_01/ukhls_rsid_grch37_dedup \
              --extract Data/step_01/ukhls_rsid_grch37_hwe.snplist \
              --keep Data/step_01/ukhls_rsid_grch37_hwe.id \
              --indep-pairwise 200 50 0.25 \
              --out Data/step_01/ukhls_rsid_grch37_ld

wc --lines  Data/step_01/ukhls_rsid_grch37_ld.prune.in \
            Data/step_01/ukhls_rsid_grch37_ld.prune.out

"${plink2}"   --bfile Data/step_01/ukhls_rsid_grch37_dedup \
              --extract Data/step_01/ukhls_rsid_grch37_ld.prune.in \
              --keep Data/step_01/ukhls_rsid_grch37_mind.id \
              --het \
              --out Data/step_01/ukhls_rsid_grch37_pruned

## b. Find Samples with |Het| > 3SD ----
# TODO: Make Het from unrelated samples
cat > Code/01e_munge_het.R << EOM
library(tidyverse)

rm(list = ls())

unrelated <- read_delim("Data/step_01/king_unrelated.id", 
                        col_names = c("fid", "iid")) %>%
  mutate(unrelated = 1)

het <- read_delim("Data/step_01/ukhls_rsid_grch37_pruned.het",
                  col_select = c("#FID", "IID", "F")) %>%
  rename(fid = 1, iid = 2, f = 3) %>%
  left_join(unrelated, by = c("fid", "iid")) %>%
  mutate(f_unrelated = ifelse(unrelated == 1, f, NA),
         f_mean = mean(f_unrelated, na.rm = TRUE),
         f_sd = sd(f_unrelated, na.rm = TRUE),
         f_scaled = (f - f_mean) / f_sd)
  
het_ids <- het %>%
  filter(between(f_scaled, -3, 3)) %>%
  select(fid, iid)

write_tsv(het_ids, "Data/step_01/het.id", col_names = FALSE)
EOM

Rscript  --no-save --no-restore Code/01e_munge_het.R

## c. Create .bed File and Convert to chr:bp ----
"${plink2}"   --bfile Data/step_01/ukhls_rsid_grch37_dedup \
              --extract Data/step_01/ukhls_rsid_grch37_hwe.snplist \
              --keep Data/step_01/het.id \
              --set-all-var-ids chr@:# \
              --rm-dup exclude-mismatch \
              --make-bed \
              --write-samples \
              --out Data/step_01/ukhls_chrbp_grch37_ld

# 10. Create Final Set of Unrelated Individuals ----
## a. Find New Set of Unrelated IIDs ---
# awk 'NR==FNR { ids[$2]; next } $2 in ids' \
#   Data/step_01/ukhls_chrbp_grch37_ld.id \
#   Data/step_01/king_unrelated.id \
#   > Data/step_01/unrelated_ld.id

cat > Code/01f_get_final_unrelated.R << EOM
library(tidyverse)
library(igraph)

rm(list = ls())

# 1. Load Data ----
fam <- read_delim("Data/step_01/ukhls_chrbp_grch37_ld.fam", 
                  delim = "\t", 
                  col_names = c("fid", "iid", "pid", "mid", "sex", "pheno"))

king <- read_delim("Data/step_01/king_relatedness.kin0",
                   delim = "\t") %>%
  select(from = ID1, to = ID2) %>%
  filter(from %in% fam[["iid"]],
         to %in% fam[["iid"]])

no_relation <- fam %>%
  filter(!(iid %in% king[["iid.x"]]),
         !(iid %in% king[["iid.y"]])) %>%
  pull(iid)

king_graph <- king %>%
  select(from = iid.x, to = iid.y) %>%
  graph_from_data_frame(directed = FALSE) %>%
  add_vertices(nv = length(no_relation), names = no_relation)

set.seed(1)
unrelated <- largest_ivs(king_graph) %>%
  sample(1) %>%
  pluck(1) %>%
  names() %>%
  as.integer()

fam %>%
  filter(iid %in% !!unrelated) %>%
  select(fid, iid) %>%
  write_tsv("Data/step_01/final_king_unrelated.id", col_names = FALSE)
EOM

Rscript  --no-save --no-restore Code/01f_get_final_unrelated.R

## b. Filter .bed file ----
"${plink2}"   --bfile "Data/step_01/ukhls_chrbp_grch37_ld" \
              --keep "Data/step_01/final_king_unrelated.id" \
              --make-bed \
              --out Data/step_01/ukhls_chrbp_grch37_unrelated


# 11. Delete Files ----
rm Data/step_01/ukhls_rsid_grch37_dedup.{bed,bim,fam}
rm Data/step_01/chr{1..22}_{info,dedup}.{bed,bim,fam}


# # 7. Find Unrelateds Using GCTA ----
# "${gcta64}" --bfile "Data/step_01/ukhls_rsid_grch37_maf_temp" \
#   --thread-num 8 \
#   --make-grm \
#   --out "Data/step_01/gcta"

# "${gcta64}" --grm "Data/step_01/gcta" \
#   --thread-num 8 \
#   --grm-cutoff 0.05 \
#   --make-grm \
#   --out "Data/step_01/gcta_unrelated"

# rm  Data/step_01/ukhls_rsid_grch37_maf_temp.{bed,bim,fam}

# # 9. Create chr:bp and hg38 versions ----
# ## a. chr:bp grch37 ----
# "${plink2}"   --bfile Data/step_01/ukhls_rsid_grch37_ld \
#               --rm-dup exclude-mismatch \
#               --set-all-var-ids chr@:# \
#               --make-bed \
#               --out Data/step_01/ukhls_chrbp_grch37_ld

# head Data/step_01/ukhls_rsid_grch37_ld.bim
# head Data/step_01/ukhls_chrbp_grch37_ld.bim

# ## b. LiftOver rsid to hg38 ----
# ### Get new coordinates ----
# awk '{ print "chr"$1, $4, $4, $2 }' OFS='\t' \
#   Data/step_01/ukhls_rsid_grch37_ld.bim \
#   > Data/step_01/ukhls_rsid_grch37_tolift.txt 

# "${executable_dir}/liftOver" -bedPlus=3 \
#   Data/step_01/ukhls_rsid_grch37_tolift.txt \
#   "${executable_dir}/hg19ToHg38.over.chain.gz"  \
#   Data/step_01/ukhls_rsid_grch37_lifted.txt \
#   Data/step_01/ukhls_rsid_grch37_unlifted.txt

# ### Extract lifted SNPs ----
# awk '{ print $4 }' \
#   Data/step_01/ukhls_rsid_grch37_lifted.txt \
#   > Data/step_01/ukhls_rsid_grch37_lifted.snplist 

# "${plink2}"   --bfile Data/step_01/ukhls_rsid_grch37_ld \
#               --extract Data/step_01/ukhls_rsid_grch37_lifted.snplist  \
#               --make-bed \
#               --out Data/step_01/ukhls_rsid_grch37_tolift

# ### Clean .bim with new coordinates ----
# cat > Code/01f_lift_bim.R << EOM
# library(tidyverse)
# library(glue)

# rm(list = ls())

# # 1. Load Data ----
# bim <- read_tsv("Data/step_01/ukhls_rsid_grch37_tolift.bim",
#                 col_names = c("chr", "snpid", "dist", "pos_grch37", "a1", "a2"))

# lifted <- read_tsv("Data/step_01/ukhls_rsid_grch37_lifted.txt",
#                    col_select = c("pos_hg38", "snpid"),
#                    col_names = c("chr_chr", "pos_hg38", "end", "snpid"))

# # 2. Merge and Create New SNP IDs ----
# new_bim <- bim %>%
#   left_join(lifted, by = "snpid") %>%
#   drop_na() %>%
#   mutate(new_id = ifelse(str_detect(snpid, "^rs"), snpid, glue("chr{chr}:{pos_hg38}"))) %>%
#   select(chr, new_id, dist, pos_hg38, a1, a2)

# write_tsv(new_bim, "Data/step_01/ukhls_rsid_grch38_lifted.bim",
#           col_names = FALSE)
# EOM

# Rscript  --no-save --no-restore Code/01f_lift_bim.R

# ### Merge new coordinates ----
# "${plink2}"   --bed Data/step_01/ukhls_rsid_grch37_tolift.bed \
#               --bim Data/step_01/ukhls_rsid_grch38_lifted.bim \
#               --fam Data/step_01/ukhls_rsid_grch37_tolift.fam \
#               --sort-vars \
#               --rm-dup exclude-mismatch \
#               --max-alleles 2 \
#               --make-pgen \
#               --out Data/step_01/ukhls_rsid_grch38_ld

# "${plink2}"   --pfile Data/step_01/ukhls_rsid_grch38_ld \
#               --make-bed \
#               --out Data/step_01/ukhls_rsid_grch38_ld

# rm Data/step_01/ukhls_rsid_grch37_tolift.{bed,bim,fam}
# rm Data/step_01/ukhls_rsid_grch38_ld.{psam,pvar,pgen}

# ## c. chr:bp hg38 ----
# "${plink2}"   --bfile Data/step_01/ukhls_rsid_grch38_ld \
#               --rm-dup exclude-mismatch \
#               --set-all-var-ids chr@:# \
#               --make-bed \
#               --out Data/step_01/ukhls_chrbp_grch38_ld

# head Data/step_01/ukhls_rsid_grch38_ld.bim
# head Data/step_01/ukhls_chrbp_grch38_ld.bim

# # 10. Create .bed Files of Unrelated Individuals ----
# ## Q. Is ths in the wrong place? Should I have had a unrelated .bed file to calculate LD?
# ## a. Creat GRM and KING Unrelated ----
# "${gcta64}" --bfile "Data/step_01/ukhls_rsid_grch37_ld" \
#   --thread-num 8 \
#   --make-grm \
#   --out "Data/step_01/gcta"

# "${gcta64}" --grm "Data/step_01/gcta" \
#   --thread-num 8 \
#   --grm-cutoff 0.05 \
#   --make-grm \
#   --out "Data/step_01/gcta_unrelated"

# "${king}" -b "Data/step_01/ukhls_rsid_grch37_ld.bed" \
#       --cpus 8 \
#       --unrelated \
#       --degree 3 \
#       --prefix "Data/step_01/king_unrelated"

# "${king}" -b "Data/step_01/ukhls_rsid_grch37_ld.bed" \
#       --cpus 8 \
#       --related \
#       --degree 3 \
#       --prefix "Data/step_01/king_unrelated"

# ## b. Filter .bed Files ----
# file_stubs=("rsid_grch37" "rsid_grch38" "chrbp_grch37" "chrbp_grch38" ) 
# for stub in "${file_stubs[@]}"; do
#   "${plink2}"   --bfile "Data/step_01/ukhls_${stub}_ld" \
#                 --keep "Data/step_01/gcta_unrelated.grm.id" \
#                 --make-bed \
#                 --out Data/step_01/ukhls_${stub}_unrelated
# done
