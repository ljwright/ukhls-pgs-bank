#!/bin/bash -l
#$ -m be
#$ -N prepare_bed
#$ -l h_rt=5:00:00
#$ -l mem=8G
#$ -pe smp 8
#$ -wd "/myriadfs/home/rmjllwr/Projects/UKHLS Fellowship/Polygenic Scores Bank"

# 0. Set up Environment -----
# qrsh -pe smp 8 -l mem=8G,h_rt=2:00:00 -now n
# qsub "Code/00_download_software.sh"

## a. Environment ----
cd "/myriadfs/home/rmjllwr/Projects/UKHLS Fellowship/Polygenic Scores Bank"
source config.cfg

## b. Download Software ----
wget -nc -P "${executable_dir}" "https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_avx2_20250129.zip"
unzip -o -j "${executable_dir}/plink2_linux_avx2_20250129.zip" -d "${executable_dir}/plink2_linux_avx2_20250129"
chmod +x "${executable_dir}/plink2_linux_avx2_20250129/plink2"

wget -nc -P "${executable_dir}" https://www.well.ox.ac.uk/~gav/resources/qctool_v2.2.0-CentOS_Linux7.8.2003-x86_64.tgz
tar -xvzf qctool_v2.2.0-CentOS_Linux7.8.2003-x86_64.tgz -C "${executable_dir}"
cp "${executable_dir}/qctool_v2.2.0-CentOS Linux7.8.2003-x86_64/qctool" "${executable_dir}"
chmod +x "${executable_dir}/qctool_v2.2.0-CentOS Linux7.8.2003-x86_64/qctool"

wget -nc -P "${executable_dir}" https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2
tar -xvjf "${executable_dir}/htslib-1.21.tar.bz2"  -C "${executable_dir}"
make -C "${executable_dir}/htslib-1.21"

mkdir -p Data/step_01/tabix
for i in {1..22}; do
  wget -nc -P "Data/step_01/tabix" https://ftp.ensembl.org/pub/grch37/current/variation/vcf/homo_sapiens/homo_sapiens-chr${i}.vcf.gz # Downloaded 2025-03-24
  "${tabix}" -p vcf "Data/step_01/tabix/homo_sapiens-chr${i}.vcf.gz"
done

wget -nc -P "${executable_dir}" https://github.com/samtools/bcftools/releases/download/1.21/bcftools-1.21.tar.bz2
tar -xvjf "${executable_dir}/bcftools-1.21.tar.bz2"  -C "${executable_dir}"
make -C "${executable_dir}/bcftools-1.21"

wget -nc -P "${executable_dir}" https://www.kingrelatedness.com/Linux-king.tar.gz
tar -xzvf  "${executable_dir}/Linux-king.tar.gz"  -C "${executable_dir}"

wget -nc -P "${executable_dir}" "https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.94.4-linux-kernel-3-x86_64.zip"
unzip -o -j "${executable_dir}/gcta-1.94.4-linux-kernel-3-x86_64.zip" -d "${executable_dir}/gcta-1.94.4-linux-kernel-3-x86_64"
chmod +x "${executable_dir}/gcta-1.94.4-linux-kernel-3-x86_64/gcta64"

wget -nc -P "${executable_dir}" http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
chmod +x "${executable_dir}/liftOver"
wget -nc -P "${executable_dir}" http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
wget -nc -P "${executable_dir}" http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

wget -nc -P "${executable_dir}" "https://github.com/choishingwan/PRSice/releases/download/2.3.5/PRSice_linux.zip"
unzip -o -j "${executable_dir}/PRSice_linux.zip" -d "${executable_dir}/PRSice_linux"
chmod +x "${executable_dir}/PRSice_linux/PRSice_linux"