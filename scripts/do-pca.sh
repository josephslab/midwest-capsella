#!/bin/bash --login
#SBATCH --job-name=midwest-pca
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --time=5:00:00
#SBATCH --partition=josephsnodes
#SBATCH --account=josephsnodes
#SBATCH --mem=50G

#filter out invariant sites
#module purge
#module load BCFtools/1.22-GCC-13.3.0
#bcftools view -v snps -e 'ALT="." || AC==0 || AC==AN' /mnt/research/josephslab/Maya/capsella/vcf/filtered/CBP_msu225_invariantWithfiltered_snps_rehead_pm.vcf.gz -Oz -o ../data/midwestCapsella-variant.vcf


#convert to plink format and filter out based on maf and ld
module purge
module load PLINK/2.00a3.7-gfbf-2023a
#plink2 --vcf ../data/midwestCapsella-variant.vcf --set-all-var-ids @:# --indep-pairwise 50 10 0.2 --geno 0.05 --maf 0.05 --max-alleles 2 --allow-extra-chr --freq --make-bed --out ../data/midwestCapsella




#do the pca only on the midwest individuals
plink2 -bfile ../data/midwestCapsella --keep ../data/midwest-pops1.txt --extract ../data/midwestCapsella.prune.in --read-freq ../data/midwestCapsella.afreq -pca --make-rel --allow-extra-chr --out ../data/midwestCapsella-pca



