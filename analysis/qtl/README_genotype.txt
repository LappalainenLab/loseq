# loseq eQTL pipeline - GENOTYPE

## 1. VCF processing ##

# see genotype/README.txt for vcf processing steps.

# Check vcf stats using:

bcftools stats -s- genotype/vcfs_post_imp_phasing/loseq19.FINAL.vcf.gz > genotype/vcfs_post_imp_phasing/loseq19.FINAL_STATS.vchk

module load python/3.5.1
plot-vcfstats -p genotype/vcfs_post_imp_phasing/loseq19.FINAL_STATS_plots genotype/vcfs_post_imp_phasing/loseq19.FINAL_STATS.vchk



## 2. Convert to Plink ##

module load plink/1.90-b3.29
PLINK_BINARY='genotype/vcfs_post_imp_phasing/loseq19.FINAL'
sbatch --job-name=plink_vcf --mem=15G --mail-type=END,FAIL --mail-user=mmartorella@nygenome.org --output=logs/200914_novaseq/vcfs/stdout/%j.out --error=logs/200914_novaseq/vcfs/stderr/%j.error --wrap="plink --make-bed --output-chr chrM --vcf ${PLINK_BINARY}.vcf.gz --out ${PLINK_BINARY}"



## 3. Create Look-up table ##

module load bcftools/1.9

bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' genotype/vcfs_post_imp_phasing/loseq19.FINAL.vcf.gz > genotype/vcfs_post_imp_phasing/loseq19.FINAL_lookuptable.txt

gzip genotype/vcfs_post_imp_phasing/loseq19.FINAL_lookuptable.txt

# Annotate lookup table:

sbatch --job-name=lookup --mem=40G --mail-type=END,FAIL --mail-user=mmartorella@nygenome.org --output=logs/200914_novaseq/vcfs/stdout/%j.out --error=logs/200914_novaseq/vcfs/stderr/%j.error scripts/qtl/genotype.create_lookup_table.R



## 4. Principle Component Analysis ##

# Plink linking file should be family ID and individual ID columns, but here they are the same (just copied plink .nosex file because already formatted as 2 columns):

cp genotype/vcfs_post_imp_phasing/loseq19.FINAL.nosex genotype/vcfs_post_imp_phasing/loseq19.FINAL.txt

mkdir -p genotype/pca/fig

module load plink/1.90-b3.29
PLINK_BINARY='genotype/vcfs_post_imp_phasing/loseq19.FINAL'
INDIV='genotype/vcfs_post_imp_phasing/loseq19.FINAL.txt'
OUT_DIR='genotype/pca'


# Run LD pruning of autosomal common SNPs with plink,
## --geno filters out all variants with missing call rates exceeding the provided value (default 0.1) to be removed, while --mind does the same for samples.
## SNPs with a 99% genotyping rate (1% missing), MAF >= 0.01, mind 0.1, hwe 1e-5
## --min 0.0625 and --indep-pairwise 50 5 0.2

sbatch --job-name=ld_prune --mem=25G --mail-type=END,FAIL --mail-user=mmartorella@nygenome.org --output=logs/200914_novaseq/vcfs/stdout/%j.out --error=logs/200914_novaseq/vcfs/stderr/%j.error scripts/qtl/pca.ld_pruning.sh ${PLINK_BINARY} ${INDIV} ${OUT_DIR}



# Check identity-by-decent (kinship >= 1/16 = 0.0625 --> 1st cousin once removed)

OUT_DIR='genotype/pca'

sbatch --job-name=ibd --mem=25G --mail-type=END,FAIL --mail-user=mmartorella@nygenome.org --output=logs/200914_novaseq/vcfs/stdout/%j.out --error=logs/200914_novaseq/vcfs/stderr/%j.error --wrap="plink --file ${OUT_DIR}/ld_pruned --genome --min 0.0625 --out ${OUT_DIR}/ld_pruned"



# Run SmartPCA to calculate PCs for a given race/ethnic group

module load eigensoft/6.1.3
OUT_DIR='genotype/pca'

## Generate the ‘.pedind’ file for smartpca

cat ${OUT_DIR}/ld_pruned.ped | cut -d ' ' -f 1-6 | awk -v group='loseq' 'BEGIN{OFS=" ";FS=" "} {print($1,$2,$3,$4,$5,group)}' > ${OUT_DIR}/ld_pruned.pedind

## Run smartpca using eigenstrat (same options as in GTEx V8)

sbatch --job-name=smartpca --mem=25G --mail-type=END,FAIL --mail-user=mmartorella@nygenome.org --output=logs/200914_novaseq/vcfs/stdout/%j.out --error=logs/200914_novaseq/vcfs/stderr/%j.error --wrap="smartpca.perl -i ${OUT_DIR}/ld_pruned.ped -a ${OUT_DIR}/ld_pruned.map -b ${OUT_DIR}/ld_pruned.pedind -k 20 -m 0 -o ${OUT_DIR}/smartpca.pca -e ${OUT_DIR}/smartpca.eval -p ${OUT_DIR}/smartpca.plot -l ${OUT_DIR}/smartpca.log"



## 5. PLOT PCA ##

# ended up plotting in R, ancestry labeled plot doesn't output correctly

module load R/3.6.1

PREFIX='LOSEQ'
HIGHLIGHT_SAMPLES='data/samplekey_508.txt'
EVEC_FILE='genotype/pca/smartpca.pca.evec'
EVAL_FILE='genotype/pca/smartpca.eval'
FIG_OUT_DIR='genotype/pca/fig/'
ANCESTRY_FILE='genotype/pca/1kg/loseq.knn_populations_1kg.txt'

sbatch --job-name=plot_pca --mem=10G --mail-type=END,FAIL --mail-user=mmartorella@nygenome.org --output=logs/200914_novaseq/vcfs/stdout/%j.out --error=logs/200914_novaseq/vcfs/stderr/%j.error scripts/qtl/pca_plot.R ${PREFIX} ${HIGHLIGHT_SAMPLES} ${EVEC_FILE} ${EVAL_FILE} ${FIG_OUT_DIR} ${ANCESTRY_FILE}




### __Project Loseq samples onto 1KG populations__


Processed 1KG VCF file here: `/gpfs/commons/groups/lappalainen_lab/data/1kg/phase3_GRCh38/merged/merged_genotypes.1kg_phase3_grch38.vcf.gz`. Obtained from Silva, her script for generating this file found here: `/gpfs/commons/groups/lappalainen_lab/data/1kg/phase3_GRCh38/merged/1kg_process_vcf.sh`


Steps for projecting SPIROMICS individuals onto 1KG populations:

* Merge 1000G and SPIROMICS VCF files__
* Merged VCF to binary PLINK format and additional filtering before LD-pruning - only autosomal SNPs, CR 99% and MAF 5% across all the samples,exclude regions of long-range LD (Table 1 of Price et al. 2008 AJHG)
* LD-pruning to create list of SNPs to be used in smartpca, output in ped format for smartpca
* Estimate PCs from 1000G samples, and project SPIROMICS samples onto those eigenvectors

1. Get overlapping sites

# ran scripts/vcf_processing/vcf_intersection.sh:

SAMPLES=('genotype/vcfs_post_imp_phasing/loseq19.FINAL.vcf.gz' '../data/1kg/phase3_GRCh38/merged/merged_genotypes.1kg_phase3_grch38.vcf.gz')
OUT_DIR=genotype/overlapping_variants_1kg/

# compressed using: bcftools view genotype/overlapping_variants_1kg/0000.vcf.gz -Oz -o genotype/overlapping_variants_1kg/loseq19.overlap.vcf.gz
# do same for 1kg, output results to: loseq19.overlap.vcf.gz and 1kg.overlap.vcf.gz
#  sites overlap


2. Merge loseq and 1kg

mkdir -p genotype/pca/1kg

# ran scripts/vcf_processing/vcf_merge.sh:

SAMPS=('genotype/overlapping_variants_1kg/loseq19.overlap.vcf.gz' 'genotype/overlapping_variants_1kg/1kg.overlap.vcf.gz')
OUTDIR='genotype/pca/1kg'

# output results: genotype/pca/1kg/loseq19.1kg.merge.vcf.gz



3.  Process merged VCF for smartpca

module load plink/1.90-b3.29
VCF='genotype/pca/1kg/loseq19.1kg.merge.vcf.gz'

sbatch --job-name=ld_prune --mem=25G --mail-type=END,FAIL --mail-user=mmartorella@nygenome.org --output=logs/200914_novaseq/vcfs/stdout/x1a.ldpruning_merged_%j.out --error=logs/200914_novaseq/vcfs/stderr/x1a.ldpruning_merged_%j.error scripts/qtl/pca.ld_pruning_merged_vcf.sh ${VCF}

# check IBD - just to see

OUT_DIR='genotype/pca/1kg'

sbatch --job-name=ibd --mem=25G --mail-type=END,FAIL --mail-user=mmartorella@nygenome.org --output=logs/200914_novaseq/vcfs/stdout/%j.out --error=logs/200914_novaseq/vcfs/stderr/%j.error --wrap="plink --file ${OUT_DIR}/loseq19.1kg.merge.vcf.gz.filtered.ld_pruned --genome --min 0.0625 --out ${OUT_DIR}/loseq19.1kg.merge.vcf.gz.filtered.ld_pruned"

# ex: cat genotype/pca/1kg/loseq19.1kg.merge.vcf.gz.filtered.ld_pruned.genome | awk '$10 > 0.1' | wc -l
# ex: cat genotype/pca/1kg/loseq19.1kg.merge.vcf.gz.filtered.ld_pruned.genome | awk '{print $1}' | sort -u | wc -l
# 127595 related pairs; 935 individuals out of 2288 individuals
# 121476 >= 0.1; 17143 >= 0.15; 312 >= 0.2; 5 >= 0.25; 0 greater than 0.3



4.  Run smartpca from eigenstrat (same options as in GTEx v8)

module load eigensoft/6.1.3

# Using option `-w poplistname` to infer eigenvectors using only individuals from a subset of populations, and then project individuals from all populations onto those eigenvectors.
# cp poplist.txt file from spiromics folder

FILE='genotype/pca/1kg/loseq19.1kg.merge.vcf.gz.filtered.ld_pruned'
POPLIST='data/1kg/poplist.txt'

sbatch --job-name=smartpca --mem=25G --mail-type=END,FAIL --mail-user=mmartorella@nygenome.org --output=logs/200914_novaseq/vcfs/stdout/x1a.smartpca_merged_%j.out --error=logs/200914_novaseq/vcfs/stderr/x1a.smartpca_merged_%j.error smartpca.perl -i ${FILE}.ped -a ${FILE}.map -b ${FILE}.pedind -w ${POPLIST} -k 20 -m 0 -o ${FILE}.pca -e ${FILE}.eval -p ${FILE}.plot -l ${FILE}.log


5.  Use KNN to inferre the population of the LOSEQ samples

# FIX INPUTS IN RSTUDIO
sbatch --job-name=knn_1kg --mem=10G --mail-type=END,FAIL --mail-user=mmartorella@nygenome.org --output=logs/200914_novaseq/vcfs/stdout/x1a.knn_%j.out --error=logs/200914_novaseq/vcfs/stderr/x1a.knn_%j.error scripts/qtl/pca.predict_1kg_pop_knn.R


# Make figures

# FIX INPUTS IN RSTUDIO
sbatch --job-name=plot_pca --mem=10gb --mail-type=END,FAIL --mail-user=mmartorella@nygenome.org --output=logs/200914_novaseq/vcfs/stdout/x1a.pca_plot_merged_%j.out --error=logs/200914_novaseq/vcfs/stderr/x1a.pca_plot_merged_%j.error scripts/pca_plot_projected_onto_1kg.R



#### SELECT COVARIATES #### - FIX

scripts/qtl/pca.select_pcs_redo.R

# PC1 and PC2 from loseq-only genotyping PCA
# PC1 and PC2 from selected donor sample gene expression PCA (per tissue)
# output per tissue covariate files (after do expression script)
# output just genotyping covariates for now:

### genotype covar file: analysis/qtl/covs/loseq.genotype.covariates.txt


* Considering elbow plot of LOSEQ only samples, PCs 1 and 2 are factors up to and including the bend and will be included (see genotype/pca/fig/LOSEQ.pca_plot.pdf and genotype/pca/fig/LOSEQ.pca_plot_ancestry.pdf). These PCs seem to capture both ancestry and genotyping panel, with PC1 significantly associated with EUR/non-EUR, and PC2 and PC3 significantly associated with genotyping panel. Only PC2 will be included because PC3 is driven by 2 ancestry samples, and PC3 is past the bend and does not contribute as meaningfully to explaining the variance.




