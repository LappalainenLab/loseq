# eQTL replication in 19 samples of hair, buccal, saliva, and urine tissues

## __1. Standard pipeline for cis-eQTLs__

## d) Map eQTLs with tensorQTL

### Cis-eQTL mapping

-**_Mapping of cis-eQTLs with tensorQTL_**

* COVARIATES: Genotype PCs 1 and 2, expression PCs explaining >= 15% of variance
* variants with minor allele frequency ≥ 0.05
* window-size 1Mb of TSS
* 10,000 permutations to control for multiple testing
* gene-level q-values with a fixed P-value interval for the estimation of pi0 (the ‘lambda’ parameter was set to 0.85)
* False discovery rate (FDR) threshold of ≤ 0.05 for significant eGenes


#######  Permutation mode

mkdir analysis/qtl/cis
FILES=($(find analysis/qtl/counts_files/ -name "*.normalized_expression.bed.gz"))
COV_ENDING='.covariates.txt'
COV_PATH='analysis/qtl/covs/'

for f in ${FILES[@]};
do
TISSUE=$(basename ${f} ".normalized_expression.bed.gz")
COV=${COV_PATH}${TISSUE}${COV_ENDING}
echo ${f} ${TISSUE}
sbatch scripts/qtl/tensorqtl.cis_permutation.sh ${f} ${COV} ${TISSUE}
done;

# Saliva no results - only 17 samples, need to fix later if want to include


####### Allpairs (cis-nominal)
mkdir -p analysis/qtl/cis/allpairs
FILES=($(find analysis/qtl/counts_files/ -name "*.normalized_expression.bed.gz"))

for f in ${FILES[@]};
do
TISSUE=$(basename ${f} ".normalized_expression.bed.gz")
COV=${COV_PATH}${TISSUE}${COV_ENDING}
echo ${f} ${TISSUE}
sbatch scripts/qtl/tensorqtl.cis_nominal.sh ${f} ${COV} ${TISSUE}
done;


####### Merge all_pairs parquet files into txt

sbatch scripts/qtl/run_parquet_to_txt.sh


####### Annotate eGenes and eVariants

sbatch scripts/qtl/tensorqtl.annotate_cis_eqtl.sh


####### Calculate pi1 replication (old was gtex overlap - do not do)

sbatch scripts/qtl/run_calculate_replication_metrics.sh


####### Plot results

sbatch analysis/qtl/run_plot_eqtl_replication.sh







-**_Mapping of conditionally independent cis-eQTLs with tensorQTL_**

Independent eQTLs at FDR 5% were identified by forward stepwise regression followed by a backwards selection step. The probe-level significance threshold was set to be the maximum beta-adjusted P-value (correcting for multiple-testing across the variants) over all eGenes.

* Forward stage: performing a scan for cis-eQTLs, correcting for all previously discovered variants and all covariates used in regular cis-eQTL mapping.
   * If the beta adjusted P-value for the lead variant was significant at the probe-level threshold, the lead variant was added to the list of discovered cis-eQTLs.
* Backwards stage: testing each variant separately, controlling for all other discovered variants.
   * If no variant was significant at the probe-level threshold the variant in question was dropped, otherwise the lead variant from this scan, which controls for all other signals (except one) found in the forward stage, was chosen as the variant that represents the signal best in the full model.

```bash
# conditionally independent eQTLs
mkdir -p eqtl/cis/independent
sbatch scripts/tensorqtl.cis_independent.sh
```

### Summary of cis-eQTL mapping and figures

* Locuszoom plots

```bash
mkdir eqtl/cis/fig_locuszoom
module load tabix
module unload python && module load python/3.7.1

GENE_NAME='ACE2'
VARIANT_ID=($(zcat eqtl/cis/spiromics.cis_eqtl.egenes.txt.gz | grep -w ${GENE_NAME} | cut -f 13))

# locuszoom plots
sbatch --job-name=lz --mail-type=FAIL --mail-user=skasela@nygenome.org --mem=20gb --output=lz_%j.out --error=lz_%j.error scripts/plot_eqtl_locuszoom.py --gene_name ${GENE_NAME} --variant_id ${VARIANT_ID} --window 250000 --output eqtl/cis/fig_locuszoom/fig_eqtl_${GENE_NAME}.pdf

# Using GTEx data
#sbatch --job-name=lz --mail-type=FAIL --mail-user=skasela@nygenome.org --mem=20gb --output=lz_%j.out --error=lz_%j.error scripts/plot_eqtl_locuszoom.gtex_v8.py --tissue_id Lung --gene_name ${GENE_NAME} --variant_id ${VARIANT_ID} --window 250000 --output fig_eqtl_${GENE_NAME}.pdf

```

### aFC estimates

aFC and 95% confidence intervals



* Raw gene count data that was normalized with DESeq size factors and log2 transformed, with the aFC arguments `--min_samps 2` and `--min_alleles 1` and including the same covariates that were used in the cis-eQTL mapping

Note: aFC estimates are capped at log2(100) = 6.64 (_it's suggested to remove eQTLs where the aFC has it the cap value_)

```bash
# Run line by line
scripts/afc_submit.sh
```

## e) Explore eQTLs for genes related to COVID-19 response

* Pull out data for ACE2 (chrX), TMPRSS2, ADAM10, and ADAM17 (and possibly others). This should be done for SPIROMICS and GTEx.

### PhenoScanner

Using PhenoScanner v2 to search if the COVID-19 related eVariants have phenotype associations.

```bash
eqtl/cis/phenoscanner/eqtl.lookup_from_phenoscanner.Rmd
```

### Respiratory infection GWAS from [UKBB GWAS v2](https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit?ts=5b5f17db#gid=227859291)

* Diagnoses - main ICD10: J06 Acute upper respiratory infections of multiple and unspecified sites
* Diagnoses - main ICD10: J22 Unspecified acute lower respiratory infection
* Diagnoses - main ICD10: J39 Other diseases of upper respiratory tract
* Diseases of the respiratory system

```bash
sbatch --job-name=lookup_gwas --mem=30gb --mail-type=FAIL --mail-user=skasela@nygenome.org --output=log/x1e.gwas_lookup_ukk_%j.out --error=log/x1e.gwas_lookup_ukk_%j.error scripts/gwas_lookup.respiratory_infections_ukbb_v2.R
```

No associations with P < 1e-05

### Colocalization

* Using both the standard method and coloc-cond/mask
  * When using coloc-cond/mask, then using `method = mask` in the GWAS dataset and LD from the corresponding 1000G population to match the ancestry of the discovery population (_CEU if the discovery population is of European ancestry_)
  * When using coloc-cond/mask, then using `method = single` in the eQTL dataset if the eGene does not have multiple independent eVaraints, otherwise using conditional p-values (_all the eGenes tested in coloc had only one independent eVariant_)
* Coloc in 500-kb region centered on lead eQTL (+/-250 kb from the lead variant)
* Priors p1 = 1e-4, p2 = 1e-4, p3 = 5e-6

Get summary stats for COVID-19 related genes (v11) that have PheWAS hits.

```bash
mkdir -p coloc/input/{qtl,gwas}
# eQTL summary stats for COVID-19 related genes
module unload python && module load python/3.6.4
cat eqtl/cis/phenoscanner/phewas_phenoscanner.covid19_related_egenes_v11.txt | cut -f 25 | grep -v "gene_id" | sort -u > coloc/input/qtl/covid19_related_egenes.txt

for CHR in {1..22} X; do
   echo ${CHR}
   sbatch --job-name=lookup_genes --mem=50gb --mail-type=FAIL --mail-user=skasela@nygenome.org --output=log/x1e.lookup_qtl_%j.out --error=log/x1e.lookup_qtl_%j.error scripts/lookup_phenotype_from_parquet.py --parquet eqtl/cis/allpairs/spiromics.cis_qtl_pairs.chr${CHR}.parquet --str1-list coloc/input/qtl/covid19_related_egenes.txt --output coloc/input/qtl/spiromics.cis_eqtl.covid19_related_egenes.allpairs.${CHR}.txt.gz
done
```

GWAS data from UKBB and GWAS catalog

* UKBB data is in GRCh37, using liftover to lift the positions over to GRCh38

```bash
# https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit?ts=5b5f17db#gid=227859291
# Download variants file
wget https://www.dropbox.com/s/puxks683vb0omeg/variants.tsv.bgz?dl=0 -O coloc/input/gwas/variants.tsv.bgz

# Lift over positions from hg19 to GRCh38
sbatch --job-name=liftover --mem=10gb --mail-type=FAIL --mail-user=skasela@nygenome.org --output=log/x1e.liftover_%j.out --error=log/x1e.liftover_%j.error scripts/ukbb.liftover_hg19_to_grch38.sh
sbatch --job-name=liftover --mem=30gb --mail-type=FAIL --mail-user=skasela@nygenome.org --output=log/x1e.liftover_%j.out --error=log/x1e.liftover_%j.error scripts/ukbb.add_grch38pos.R

# Download phenotypes file to get information about #cases
wget https://www.dropbox.com/s/d4mlq9ly93yhjyt/phenotypes.both_sexes.tsv.bgz -O coloc/input/gwas/phenotypes.both_sexes.tsv.gz
```

Run coloc for COVID-19 related genes that have associations with blood cell traits (_EFO parent category = Hematological measurement_), pulmonary function traits (_EFO parent category = Pulmonary function measurement_), or respiratory system disease (_EFO parent category = Respiratory disease_)

* Coloc-standard

```bash
mkdir -p coloc/result/eqtl/single/{fig_locuscompare,fig_locuszoom,fig_sensitivity}

# Get sdY for expression
pheno_file="expression/normalized/spiromics.normalized_expression.bed.gz"
cov_file="cov/spiromics.covariates.txt"
interaction_file="NA"
cores="5"
out_file="coloc/input/qtl/spiromics.sdY.txt"
sbatch --job-name=get_sdY --cpus-per-task=${cores} --mem=50gb --mail-type=FAIL --mail-user=skasela@nygenome.org --output=log/x1e.get_sdY_%j.out --error=log/x1e.get_sdY_%j.error scripts/coloc.get_sdY.R ${pheno_file} ${cov_file} ${interaction_file} ${cores} ${out_file}

# Run coloc-standard
METHOD='standard'
sbatch scripts/coloc.run_for_phenoscanner_hits.submit.sh ${METHOD}
## Run one
## G="FOXO1"
## sbatch --job-name=coloc --mem=40gb --mail-type=FAIL --mail-user=skasela@nygenome.org --output=log/x1e.coloc_${G}_${METHOD}_%j.out --error=log/x1e.coloc_${G}_${METHOD}_%j.error scripts/coloc.run_for_phenoscanner_hits.R ${G} ${METHOD}
```

* Coloc-cond/mask: masking in GWAS dataset, LD data based on the corresponding 1000G population

```bash
mkdir -p coloc/result/eqtl/mask/{fig_locuscompare,fig_locuszoom,fig_sensitivity}

# Run coloc-cond/mask: masking in GWAS dataset
METHOD='mask'
sbatch scripts/coloc.run_for_phenoscanner_hits.submit.sh ${METHOD}
## Run one
## G="ERMP1"
## sbatch --job-name=coloc --mem=40gb --mail-type=FAIL --mail-user=skasela@nygenome.org --output=log/x1e.coloc_${G}_${METHOD}_%j.out --error=log/x1e.coloc_${G}_${METHOD}_%j.error scripts/coloc.run_for_phenoscanner_hits.R ${G} ${METHOD}
```

## __2. Cell type interacting eQTLs__

## a) Cell type composition analysis by xCell

* Run xCell

## b) Map epithelial interacting eQTLs (ieQTLs)

* Potential epithelial-specific effects


## __3. Maybe: Run ANEVA-DOT (& check GTEx) for potential rare variants__
