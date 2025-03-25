#!/bin/bash -l
#$ -m be
#$ -N prepare_bed
#$ -l h_rt=5:00:00
#$ -l mem=8G
#$ -pe smp 8
#$ -wd "/home/rmjllwr/Scratch/Projects/UKHLS Fellowship/Polygenic Scores Bank"

# 0. Set up Environment -----
# qrsh -pe smp 8 -l mem=8G,h_rt=4:00:00 -now n
# qsub "Code/01_prepare_bed.sh"
## a. Environment ----
cd "/home/rmjllwr/Scratch/Projects/UKHLS Fellowship/Polygenic Scores Bank"
source config.cfg

mkdir -p Data/step_01

module -f unload compilers mpi gcc-libs
module load r/recommended
module load lapack/3.8.0/gnu-4.9.2

## b. Download Software ----
wget -nc -P "${executable_dir}" http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
chmod +x "${executable_dir}/liftOver"
wget -nc -P "${executable_dir}" http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
wget -nc -P "${executable_dir}" http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

wget -nc -P "${executable_dir}" https://www.well.ox.ac.uk/~gav/resources/qctool_v2.2.0-CentOS_Linux7.8.2003-x86_64.tgz
tar -xvzf qctool_v2.2.0-CentOS_Linux7.8.2003-x86_64.tgz -C "${executable_dir}"
cp "${executable_dir}/qctool_v2.2.0-CentOS Linux7.8.2003-x86_64/qctool" "${executable_dir}"
chmod +x "${executable_dir}/qctool" 

wget -nc -P "${executable_dir}" https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2
tar -xvjf "${executable_dir}/htslib-1.21.tar.bz2"  -C "${executable_dir}"
make -C "${executable_dir}/htslib-1.21"
"${tabix}" --version

# mkdir -p Data/step_01/tabix
# for i in {1..22}; do
#   wget -nc -P "Data/step_01/tabix" https://ftp.ensembl.org/pub/grch37/current/variation/vcf/homo_sapiens/homo_sapiens-chr${i}.vcf.gz # Downloaded 2025-03-24
#   "${tabix}" -p vcf "Data/step_01/tabix/homo_sapiens-chr${i}.vcf.gz"
# done


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
  write_tsv("Data/step_01/extraid.tsv", col_names = FALSE)
EOM

Rscript  --no-save --no-restore Code/01a_amend_sample_file.R

# 2. Drop MAF < 0.01 and Extra IDs and Resave to Remove Trailing Space ----
# TODO: IS --snps-only an error? Can look at sumstats at frequency of INDELs.
# TODO: Split out the three steps within here --maf --max-alleles --snps-only
for i in {1..22}; do
    ${plink2} \
      --gen "${geno_dir}/chr${i}.UKHLS.UK10K+1KG-imputed.phwe_1e-4.info_0.4.filtered.gen.gz" ref-unknown \
      --sample "Data/step_01/amended.sample" \
      --oxford-single-chr ${i} \
      --remove "Data/step_01/extraid.tsv" \
      --maf 0.01 \
      --memory 12000 --threads 8 \
      --recode oxford bgz --out Data/step_01/chr${i}_maf
done

# 3. Filter on INFO Scores ----
## a. Calculate INFO Scores ----
for i in {1..22}; do
    "${qctool}" \
        -g Data/step_01/chr${i}_maf.gen.gz \
        -snp-stats -osnp Data/step_01/chr${i}_snp_stats.txt
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
for i in {1..22}; do
    ${plink2} \
      --gen "Data/step_01/chr${i}_maf.gen.gz" ref-unknown \
      --sample "Data/step_01/chr${i}_maf.sample" \
      --oxford-single-chr ${i} \
      --extract Data/step_01/chr${i}_info.txt \
      --memory 8000 --threads 8 \
      --make-bed --out Data/step_01/chr${i}_info

  rm Data/step_01/chr${i}_maf.{gen.gz,sample}
done

# 4. Clean rsids ----
## a. Find SNPs to re-ID ----
for i in {1..22}; do
  awk -F'\t' '$2 !~ /^rs/' OFS='\t' \
    Data/step_01/chr${i}_info.bim \
    > Data/step_01/chr${i}_toid.txt
done

## b. Use tabix to find matches ----
for i in {1..22}; do
  # wget -nc -P "Data/step_01/tabix" https://ftp.ensembl.org/pub/grch37/current/variation/vcf/homo_sapiens/homo_sapiens-chr${i}.vcf.gz # Downloaded 2025-03-24
  # "${tabix}" -p vcf "Data/step_01/tabix/homo_sapiens-chr${i}.vcf.gz"

  while read chr id dist pos ref alt; do
    "${tabix}" "Data/step_01/tabix/homo_sapiens-chr${i}.vcf.gz" ${chr}:${pos}-${pos}
  done < "Data/step_01/chr${i}_toid.txt" > "Data/step_01/chr${i}_tabix.txt"
done

## c. Re-ID .bim File ----
cat > Code/01c_reid_bim.R << EOM
library(tidyverse)
library(glue)
library(tictoc)

rm(list = ls())

reid_bim <- function(i){
  # 1. Load Data ----
  bim <- glue("Data/step_01/chr{i}_info.bim") %>%
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
    mutate(id = case_when(str_sub(id, 1, 2) == "rs" ~ id,
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

Rscript  --no-save --no-restore Code/01c_reid_bim.R


# 5. Combine .bed Files and Remove Duplicates ----
# TODO: Not sure this will work with 22 chromosomes
# TODO: Think whether --max-alleles is a suitable option
rm Data/step_01/merge_list.txt
for i in {1..22}; do
    ${plink2}   --bed Data/step_01/chr${i}_info.bed \
                --bim Data/step_01/chr${i}_reid.bim \
                --fam Data/step_01/chr${i}_info.fam \
                --rm-dup exclude-mismatch \
                --max-alleles 2 \
                --make-bed --out Data/step_01/chr${i}_dedup
    
    echo Data/step_01/chr${i}_dedup >> Data/step_01/merge_list.txt
    # echo "Data/step_01/chr${i}_info.bed Data/step_01/chr${i}_reid.bim Data/step_01/chr${i}_info.fam" >> Data/step_01/merge_list.txt
done

${plink2}   --pmerge-list Data/step_01/merge_list.txt bfile \
            --rm-dup exclude-mismatch \
            --make-bed --out Data/step_01/ukhls_rsid_hg37_dedup

for i in {1..22}; do
    rm Data/step_01/chr${i}_{info,reid,dedup}.{bed, fam}
done


# 6. Remove High Missingness ----
## a. Variant-Level ----
${plink2}   --bfile Data/step_01/ukhls_rsid_hg37_dedup \
            --geno 0.03 \
            --make-bed --out Data/step_01/ukhls_rsid_hg37_geno

## b. Sample-Level ----
${plink2}   --bfile Data/step_01/ukhls_rsid_hg37_geno \
            --mind 0.02 \
            --make-bed --out Data/step_01/ukhls_rsid_hg37_mind

rm Data/step_01/ukhls_rsid_hg37_{dedup,geno}.{bed,fam}


# 7. Remove SNPs with non-AGCT variants ----
awk '$5 ~ /^[AGCT]+$/ && $6 ~ /^[AGCT]+$/ { print $2 }' \
  Data/step_01/ukhls_rsid_hg37_mind \
  > Data/step_01/agct_ids.txt

plink2  --bfile Data/step_01/ukhls_rsid_hg37_mind \
        --extract Data/step_01/agct_ids.txt \
        --make-bed \
        --out Data/step_01/ukhls_rsid_hg37_agct
rm Data/step_01/ukhls_rsid_hg37_mind.{bed,fam,bim}


# 8. Remove High LD Regions ----

# 9. Remove IBD Outliers ----

# 10. Create GRM matrix ----

# 11. Create chr:bp ID versions ----
## a. hg37 ----

## b. LiftOver to hg38 ----

# Q: Do in two steps?
# Q: Split out MAF from info to hard call?
# TODO: Remove incorrect samples.
# TODO: Use sample order instead as IID. Save sample_file so when the bank of PGS are made, correct IIDs can be used.
# Remove the ";" if it creates a data management issue
  # Tim would remove and note how many (and which ones) are being stripped
  # Then for each sumstat file, calculate how many are not in analysis - see if any in the file I dropped are in the sumstats
# Check my EA build for the rsid. 
# info score? TODO: Ask for info scores to check I have calculated them correctly.
# TODO: Remove duplicates on chr:bp and rsid. 
  # Tim assumes you won't lose much information by dropping a small amount of stuff but at least you know what you have now.
# TODO: Get Mandy's PGS and compare against the ones I've made.
# TODO: Send my code and a summary to Gemma.
# If I have to port .bim and sumstats on to the same build - which should go on the same built - should hg38 be preferred?
