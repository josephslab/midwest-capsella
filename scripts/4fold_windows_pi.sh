#!/bin/bash --login
#SBATCH --job-name=MidwestPiWindows
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --time=5:00:00
#SBATCH --partition=josephsnodes
#SBATCH --account=josephsnodes
#SBATCH --mem=50G
#SBATCH -o /mnt/scratch/wils1582/pixycbp.out
#SBATCH -e /mnt/scratch/wils1582/pixycbp.error
#
# Running PIXY on 4-fold degenerate sites for Chrissy in windows
# October 8, 2026
#
# M. Wilson Brown
#

# load modules
module purge
module load Miniforge3/25.11.0-1
module load tabixpp/1.1.2-GCC-12.3.0

# activate conda environment
conda activate pixy

# move to working directory
cd /mnt/research/josephslab/Maya/forChrissy

#### SET THE PATH FOR VARIABLES
VCF=/mnt/research/josephslab/Maya/capsella/vcf/filtered/CBP_msu225_invariantWithfiltered_snps_rehead_pm.vcf.gz
POPS=/mnt/research/josephslab/Maya/forChrissy/midwest_pixy_pops.txt
OUTDIR=/mnt/research/josephslab/Maya/forChrissy/
PREFIX=Midwest_4-fold_10kb

# Optional VARS
# This 4-fold degeneracy sites file is generated in cbp_nyc_variation/hpcc_scripts/withinNYC_cbp_pixy_popgen.sh
# Crucially, it uses the paramasked VCF to call 4-fold sites
SITES=/mnt/research/josephslab/Maya/forChrissy/1-based_4-fold_sites_CBP_allSites_msu_maf.txt

## Then, you can run PIXY
 pixy --stats pi  \
	--sites_file $SITES \
	--window_size 10000 \
 	--vcf $VCF \
	--populations $POPS \
	--n_cores 20 \
	--output_folder "$OUTDIR" \
	--output_prefix "$PREFIX"
# deactivate environment
conda deactivate
