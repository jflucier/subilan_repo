#!/bin/bash
#SBATCH --time=08:00:00 # Each file is fast, so less time needed
#SBATCH --mem=50G       # Less RAM needed for single files
#SBATCH -n 16            # Fewer threads per task
#SBATCH -J diann_indiv
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
cp /lustre07/scratch/jflucier/2021_long_covid/.predicted.speclib $SLURM_TMPDIR/predicted_library.speclib
cp /lustre07/scratch/jflucier/2021_long_covid/UP000005640_9606_combo.fasta $SLURM_TMPDIR/

# Get the specific .d file for this array task
FILES=($RAW_DIR/*/*.d)
CURRENT_FILE=${FILES[$SLURM_ARRAY_TASK_ID]}

export SINGULARITY_BIND=$SLURM_TMPDIR:$SLURM_TMPDIR,/lustre07/scratch/jflucier/2021_long_covid:/lustre07/scratch/jflucier/2021_long_covid
export APPTAINER_BIND=$SLURM_TMPDIR:$SLURM_TMPDIR,/lustre07/scratch/jflucier/2021_long_covid:/lustre07/scratch/jflucier/2021_long_covid

echo "running diann"
singularity exec --writable-tmpfs -e $DIANN \
diann-linux --threads 16 --verbose 2 \
--f "$CURRENT_FILE" \
--fasta "$SLURM_TMPDIR/UP000005640_9606_combo.fasta" \
--lib "$SLURM_TMPDIR/predicted_library.speclib" \
--out "$SLURM_TMPDIR/out/$(basename "$CURRENT_FILE").tsv" \
--temp "$SLURM_TMPDIR/out" \
--gen-spec-lib --quant-ori-names \
--qvalue 0.01 --reannotate --peptidoforms \
--mass-acc 20 --mass-acc-ms1 20 --individual-mass-acc --individual-windows \
--reanalyse --rt-profiling

echo "copying results to $OUT_DIR"
cp $SLURM_TMPDIR/out/* $OUT_DIR/
