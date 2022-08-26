#!/bin/bash
#SBATCH --job-name=open_targets
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mmartorella@nygenome.org
#SBATCH --mem=100G
##SBATCH --time=8:00:00
#SBATCH --output=logs/200914_novaseq/dz_targets/stdout/%A_%a.out
#SBATCH --error=logs/200914_novaseq/dz_targets/stderr/%A_%a.error
#SBATCH --array=1-9 # change to match length of files list

module purge
module load R/3.6.0

IDX=$((${SLURM_ARRAY_TASK_ID} - 1))

# Inputs
SOURCES=(clingen eva gene2phenotype genomics_england orphanet ot_genetics_portal phewas_catalog uniprot_literature uniprot_variants)
ASSOC_PATH=$(find data/open_targets/evidence/ -type d -name sourceId=${SOURCES[${IDX}]})
OUT=analysis/dz_association/220209_20dz_${SOURCES[${IDX}]}.txt


Rscript analysis/dz_association/open_targets_file_processing.R ${ASSOC_PATH} ${OUT}



