# Run tensorqtl with singularity on slurm

Wiki helppage: https://wiki.nygenome.org/display/rescomp/Submit+jobs+to+slurm+gpu+nodes+using+singularity

```
# start slurm-gpu session
srun --mem=40g -p gpu --pty --gres=gpu:tesla:1 -u bash -i
```

## Setting up tensorqtl to be used with singularity

Create `~/tensorqtl` to your home dir and build sif file from the latest image into ~/tensorqtl (see wiki for instructions)

To get the newest updates, wait just 20 min after a github change, and then pull the latest version:

```bash
cd ~/tensorqtl
module load singularity
singularity cache clean -f
#build sif file from latest image
singularity build tensorqtl.sif docker://gcr.io/nygc-public/tensorqtl
```

## Using tensorqtl with singularity from command line

```bash
# cd /path-to-dir, where you want to submit your jobs

module load singularity

# -e : stop importing your local enviroment inside container
# -B : bind your dir in the specific dir in the container
# bind the current path, and dir where data is stored (if in some other dir)
singularity run --nv -e -B /gpfs/commons/home/skasela/test/tensorqtl_image,/gpfs/commons/home/skasela/test/tensorqtl/geno:/geno,/gpfs/commons/home/skasela/test/tensorqtl/input:/input ~/tensorqtl/tensorqtl.sif

export GENO_PATH='/geno/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered'
export PHENO_FILE='/input/topmed_mesa.exam_1.black.normalized_methylation.bed.gz'
export COV_FILE='/input/cov/GEUVADIS.445_samples.covariates.txt'
export PREFIX='my_results/GEUVADIS.445_samples'

python3 -m tensorqtl ${GENO_PATH} ${PHENO_FILE} ${PREFIX} --mode cis --covariates ${COV_FILE} --permutations 100 --fdr 0.05 --seed 124456677
```

Results appear in /path-to-dir/my_results (this dir needs to already exist)

## Run tensorqtl using slurm scripts - example slurm script

```bash
#!/bin/bash
#SBATCH --job-name=test_tensorqtl                 # Job name
#SBATCH --partition=gpu                    # Partition Name
#SBATCH --mail-type=END,FAIL                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=user@nygenome.org        # Where to send mail
#SBATCH --mem=40gb                            # Job memory request
#SBATCH --time=4:00:00                       # Time limit 4 hours
#SBATCH --output=stdout_%j.log               # Standard output and error log
#SBATCH --error=error_%j.log
#SBATCH --gres=gpu:1

module load singularity
singularity exec --nv -B ~/tensorqtl,~/example-input-dir:/dir-inside-container ~/tensorqtl/tensorqtl.sif ~/tensorqtl/run.sh
```

## Check the version of tensorqtl

Instead of running your script you would add bash in the end which will take you to terminal

```bash
singularity exec -i ~/tensorqtl/tensorqtl.sif bash

# Then go to /opt/tensorqtl directory
head -n 1 /opt/tensorqtl/tensorqtl/__init__.py
```
