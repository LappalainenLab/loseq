# eQTL mapping in Loseq samples (4 tissue types, 4 collections)
# Molly Martorella
# Non-loseq specific scripts provided by Silva Kasela


### Gene quantification/normalization ###

#SKIPPED STEP 1:
find analysis/qtl/counts_files/ -name "*_maxdepthsamps_*" -exec cp {} analysis/qtl/loseq18/counts_files/ \;

1. Create counts matrices (one for each tissue type) with 1 sample/donor. Select sample with highest protein coding/lncRNA read depth. Do for both raw counts and tpms.

# used script: scripts/qtl/expression.split_counts.R
# file output directory: analysis/qtl/counts_files/


2. Get normalized counts in BED format, update column names to match vcf

# Processing gene expression data according to the [GTEx pipeline](https://github.com/broadinstitute/gtex-pipeline/tree/master/qtl):

# The expression data are normalized as follows:
# * Read counts are normalized between samples using TMM (Robinson & Oshlack, Genome Biology, 2010)
# * Genes are selected based on the following expression thresholds:
# 	* ≥0.1 TPM in ≥20% samples AND
# 	* ≥6 reads (unnormalized) in ≥20% samples
# * Each gene is inverse normal transformed across samples.


TISSUES=('buccal' 'hair' 'saliva' 'urine')
TPM_END='_maxdepthsamps_tpms.txt'
CTS_END='_maxdepthsamps_raw.txt'
OUT_DIR='analysis/qtl/loseq18/counts_files/'
OUT_END='.normalized_expression.bed'

### TMM

for tis in ${TISSUES[@]};
do
NORMALIZATION_METHOD="tmm"
TPM_FILE=${OUT_DIR}${tis}${TPM_END}
COUNTS_FILE=${OUT_DIR}${tis}${CTS_END}
LINKING_FILE='data/samplekey_485.txt'
ANNOT_FILE='references/hg38/gencode26/gencode.v26.annotation.gtf'
VCF_SAMPLES='genotype/loseq18/samples_list_18.txt'
OUT_FILE=${OUT_DIR}${tis}${OUT_END}
sbatch --job-name=expr --mem=10G --mail-type=END,FAIL --mail-user=mmartorella@nygenome.org --output=logs/200914_novaseq/qtl/stdout/process_expression_%j.out --error=logs/200914_novaseq/qtl/stderr/process_expression_%j.error scripts/qtl/expression.process_data.R ${NORMALIZATION_METHOD} ${TPM_FILE} ${COUNTS_FILE} ${LINKING_FILE} ${ANNOT_FILE} ${VCF_SAMPLES} ${OUT_FILE}
done;

# Index bed files

FILES=($(find analysis/qtl/loseq18/counts_files/ -name "*.bed"))

for f in ${FILES[@]};
do
sbatch --job-name=expr_index --mem=10G --mail-type=END,FAIL --mail-user=mmartorella@nygenome.org --output=logs/200914_novaseq/qtl/stdout/index.expression_%j.out --error=logs/200914_novaseq/qtl/stderr/index.expression_%j.error scripts/qtl/index_bed.sh ${f}
done;


### DESeq2
# DESeq2 size factors normalized files for aFC calculations:

OUT_END='.deseq_log2_expression.bed'

for tis in ${TISSUES[@]};
do
NORMALIZATION_METHOD="deseq2"
TPM_FILE=${OUT_DIR}${tis}${TPM_END}
COUNTS_FILE=${OUT_DIR}${tis}${CTS_END}
LINKING_FILE='data/samplekey_485.txt'
ANNOT_FILE='references/hg38/gencode26/gencode.v26.annotation.gtf'
VCF_SAMPLES='genotype/loseq18/samples_list_18.txt'
OUT_FILE=${OUT_DIR}${tis}${OUT_END}
sbatch --job-name=expr --mem=10G --mail-type=END,FAIL --mail-user=mmartorella@nygenome.org --output=logs/200914_novaseq/qtl/stdout/process_expression_%j.out --error=logs/200914_novaseq/qtl/stderr/process_expression_%j.error scripts/qtl/expression.process_data.R ${NORMALIZATION_METHOD} ${TPM_FILE} ${COUNTS_FILE} ${LINKING_FILE} ${ANNOT_FILE} ${VCF_SAMPLES} ${OUT_FILE}
done;

# Index bed files

FILES=($(find analysis/qtl/loseq18/counts_files/ -name "*.deseq_log2_expression.bed"))

for f in ${FILES[@]};
do
sbatch --job-name=expr_index --mem=10G --mail-type=END,FAIL --mail-user=mmartorella@nygenome.org --output=logs/200914_novaseq/qtl/stdout/index.expression_%j.out --error=logs/200914_novaseq/qtl/stderr/index.expression_%j.error scripts/qtl/index_bed.sh ${f}
done;


3. Run PCA on per tissue expression matrices to get covariate PCs

TISSUES=('buccal' 'hair' 'saliva' 'urine')
TPM_END='_maxdepthsamps_tpms.txt'
CTS_END='_maxdepthsamps_raw.txt'
CTS_DIR='analysis/qtl/loseq18/counts_files/'
OUT_DIR='analysis/qtl/loseq18/covs/'

for tis in ${TISSUES[@]};
do
CTS_FILE=${CTS_DIR}${tis}${CTS_END}
TPM_FILE=${CTS_DIR}${tis}${TPM_END}
KEY='data/samplekey_485.txt'
GENO_COV='analysis/qtl/loseq18/covs/loseq.genotype.covariates.txt'
sbatch --job-name=expr_cov --mem=10G --mail-type=END,FAIL --mail-user=mmartorella@nygenome.org --output=logs/200914_novaseq/qtl/stdout/%x_%j.out --error=logs/200914_novaseq/qtl/stderr/%x_%j.error scripts/qtl/expression.pca.covars.R ${CTS_FILE} ${TPM_FILE} ${OUT_DIR} ${KEY} ${GENO_COV} ${tis}
done;











### DID NOT DO THIS ###


2. Calculate PEER factors

Probabilistic Estimation of Expression Residuals (PEER, [Stegle et al., 2010](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000770)) is a Bayesian framework that uses factor analysis methods to infer hidden determinants and their effects from gene expression profiles.

Thus, PEER factors are latent factors that capture and correct for technical variation (e.g., batch effects) and unwanted biological variation (e.g., cellular heterogeneity) in molecular datasets.

### __Estimate PEER factors__

* without covariates
* with genotype PCs and gender (imputed sex using Plink)

```bash
mkdir -p peer/factor_relevance
INPUT='expression/normalized/spiromics.normalized_expression.bed'
COV=('no_cov' 'cov_genopc_sex')

for i in ${COV[@]}; do
    sbatch --job-name=peer --mem=25gb --mail-type=FAIL --mail-user=skasela@nygenome.org --output=log/x1c.peers_${i}_%j.out --error=log/x1c.peers_${i}_%j.error scripts/peer.estimate_factors.R ${INPUT} ${i}
done
```

Partition the variance attributable to PEERs in gene expression data

```bash
mkdir peer/fig
sbatch --job-name=peer --mem=25gb --cpus-per-task=9 --mail-type=FAIL --mail-user=skasela@nygenome.org --output=log/x1c.peers_variance_%j.out --error=log/x1c.peers_variance_%j.error scripts/peer.variance_partition.R
```

### __Optimal number of PEERs to choose__

The number of PEER factors to adjust for was chosen to maximize cis-eGene discovery. The cis-eQTL discovery pipeline was run with increments of 5 PEERs with reduced number of permutations (100 instead of 10,000).

-**_Prepare covariates file_**

```bash
COV=("no_cov" "cov_genopc_sex")
K='6' # how many test files to make (6*5 = 30 - include up to 30 PEERs)

for i in ${COV[@]}; do
    mkdir -p peer/optimal_number/${i}/cov
    sbatch --job-name=peer_files --mem=15gb --mail-type=FAIL --mail-user=skasela@nygenome.org --output=log/x1c.prepare_cov_${i}_%j.out --error=log/x1c.prepare_cov_${i}_%j.error scripts/optimization.prepare_cov.R ${i} ${K}
done
```

-**_Run cis-eQTL mapping with increments of 5 PEER factors and with fewer numbers of permutations_**

```bash
COV=("no_cov" "cov_genopc_sex")

for i in ${COV[@]}; do
    sbatch scripts/optimization.tensorqtl_cis_permutation.sh ${i}
done
```

-**_Results of the optimization of PEER factors for cis-mQTL mapping_**

```bash
sbatch --job-name=peers_opt --mem=5gb --mail-type=FAIL --mail-user=skasela@nygenome.org --output=log/x1c.peers_opt_%j.out --error=log/x1c.peers_opt_%j.error --wrap="Rscript -e \"rmarkdown::render('scripts/optimization.overview_of_the_results.Rmd', output_dir = 'peer/optimal_number/', output_file = 'optimal_number_of_peers.html')\""
```

Based on these results and to avoid potential overfitting, optimal number of covariates for mQTL mapping:

* 15 PEER factors + genotype PCs and sex

-**_Select 15 PEERs factors, top 4 genotype PCs and sex for eQTL mapping_**

```bash
mkdir -p cov
cp peer/optimal_number/no_cov/cov/spiromics.covariates.set_3.txt cov/spiromics.covariates.txt
```
