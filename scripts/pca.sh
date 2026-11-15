#!/bin/bash --login
#SBATCH --job-name=midwest-pca
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --time=5:00:00
#SBATCH --partition=josephsnodes
#SBATCH --account=josephsnodes
#SBATCH --mem=50G


module purge
module load BCFtools/1.22-GCC-13.3.0
bcftools view -v snps -e 'ALT="." || AC==0 || AC==AN' /mnt/research/josephslab/Maya/capsella/vcf/filtered/CBP_msu225_invariantWithfiltered_snps_rehead_pm.vcf.gz -Oz -o ../data/midwestCapsella-variant.vcf



module purge
module load PLINK/2.00a3.7-gfbf-2023a
plink2 --vcf ../data/midwestCapsella-variant.vcf --set-all-var-ids @:# --indep-pairwise 100 5 0.2 --geno 0.05 --max-alleles 2 --maf 0.05 --allow-extra-chr --make-bed --out ../data/midwestCapsella

plink2 -bfile ../data/midwestCapsella --extract ../data/midwestCapsella.prune.in -pca --allow-extra-chr --out ../data/midwestCapsella-pca

