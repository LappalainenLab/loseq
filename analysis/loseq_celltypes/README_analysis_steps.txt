git clone https://github.com/BNadel/GEDIT.git

paper: https://academic.oup.com/gigascience/article/10/2/giab002/6137724#228569927


1. Format cts matrix

# row names should be gene names (hgnc symbols)
# using cts thresholded and downsampled to 5e6
# using unfiltered raw cts --> GEDIT filters and quantile normalizes cts

analysis/loseq_celltypes/gedit_formatting.R


2. Run GEDIT.py

# GEDITv2.0 runs entire pipeline in python2 (vs 1.7, which you need to run python then R), but had issues with pathing to necessary scripts to run
# GEDITv1.7 output is an Rscript command to run (will also tell location of intermediate file paths - usually in GEDITv1.7/scripts/scratch and GEDITv1.7/predictions).

sbatch scripts/gedit/gedit_python_wrapper.sh

# script inputs:
#MIX='analysis/loseq_celltypes/inputs/loseq5mil.txt'
#REF='scripts/gedit/GEDIT/ReferenceMatrices/BlueCodeV1.0.tsv'
#OUT='analysis/loseq_celltypes/results/loseq5mil_bluecode.tsv'

#python scripts/gedit/GEDIT/GEDITv1.7/scripts/GEDIT.py -mix $MIX -ref $REF -outFile $OUT


3. Run GLM_decon.R

# can reference the R command output in the standard out for inputs.

sbatch scripts/gedit/gedit_R_wrapper.sh

# script inputs:
# MIX='scripts/gedit/GEDIT/GEDITv1.7/scripts/scratch/loseq5mil.txt_BlueCodeV1.0.tsv_50_Entropy_0.0_ScaledMix.tsv'
# REF='scripts/gedit/GEDIT/GEDITv1.7/scripts/scratch/loseq5mil.txt_BlueCodeV1.0.tsv_50_Entropy_0.0_ScaledRef.tsv'
# OUT='analysis/loseq_celltypes/results/loseq5mil_bluecode.tsv'

# Rscript scripts/gedit/GEDIT/GEDITv1.7/scripts/GLM_Decon.R $MIX $REF $OUT


4. Plot results

analysis/loseq_celltypes/gedit_plotting.R





