#!/bin/bash
#SBATCH --time=02:00:00 # Each file is fast, so less time needed
#SBATCH --mem=50G       # Less RAM needed for single files
#SBATCH -c 16            # Fewer threads per task
#SBATCH -J diann_speclib
#SBATCH -o /lustre07/scratch/jflucier/2021_long_covid/logs/diann-%A_%a.out
#SBATCH -A def-gimap5
#SBATCH -N 1

module load apptainer
DIANN=/lustre07/scratch/jflucier/2021_long_covid/diann-2.1.0.sif
BASE_DIR=/lustre07/scratch/jflucier/2021_long_covid
RAW_DIR=$BASE_DIR/upload/TIMS_top-MS_analysis
OUT_DIR=$BASE_DIR/diann_quant_files

mkdir -p $OUT_DIR
mkdir -p $SLURM_TMPDIR/out

echo "copying speclib"
cp /lustre07/scratch/jflucier/2021_long_covid/UP000005640_9606_combo.fasta $SLURM_TMPDIR/

# Get the specific .d file for this array task
FILES=($RAW_DIR/*/*.d)
CURRENT_FILE=${FILES[$SLURM_ARRAY_TASK_ID]}

export SINGULARITY_BIND=$SLURM_TMPDIR:$SLURM_TMPDIR,/lustre07/scratch/jflucier/2021_long_covid:/lustre07/scratch/jflucier/2021_long_covid
export APPTAINER_BIND=$SLURM_TMPDIR:$SLURM_TMPDIR,/lustre07/scratch/jflucier/2021_long_covid:/lustre07/scratch/jflucier/2021_long_covid

# on narval
cd $BASE_DIR

echo "generating spec lib"
singularity exec  --writable-tmpfs -e \
${DIANN} \
diann-linux --threads 60 --verbose 2 \
--f $SLURM_TMPDIR/data/Sheela_DIA_sample_001b_Slot2-22_1_12799.d \
--fasta "$SLURM_TMPDIR/UP000005640_9606_combo.fasta" \
--predict --fasta-search \
--min-pep-len 7 --max-pep-len 30 \
--min-pr-mz 300 --max-pr-mz 1800 \
--min-fr-mz 200 --max-fr-mz 1800 \
--min-pr-charge 2 --max-pr-charge 4 \
--cut K*,R* --missed-cleavages 1 \
--met-excision --unimod4 --var-mods 2 \
--var-mod UniMod:35,15.994915,M --var-mod UniMod:1,42.010565,*n \
--out-lib "$SLURM_TMPDIR/out/predicted_library.speclib"

if [ ! -f "$SLURM_TMPDIR/out/predicted_library.speclib" ] && [ -f "~/.predicted.speclib" ]; then
    cp ~/.predicted.speclib "$SLURM_TMPDIR/out/predicted_library.speclib"
fi

echo "copying results to $OUT_DIR"
cp $SLURM_TMPDIR/out/* $OUT_DIR/
