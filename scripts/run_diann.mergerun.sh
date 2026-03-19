#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mem=1000G      # Merging needs more RAM than single files but less than processing all RAW
#SBATCH -n 60
#SBATCH -J diann_merge
#SBATCH -A def-gimap5
#SBATCH -N 1
#SBATCH -o /lustre07/scratch/jflucier/2021_long_covid/diann-merge.out

module load apptainer
BASE_DIR=/lustre07/scratch/jflucier/2021_long_covid
RAW_DIR=$BASE_DIR/upload/TIMS_top-MS_analysis
QUANT_DIR=$BASE_DIR/diann_quant_files
FINAL_OUT=$BASE_DIR/diann_final_results
FINAL_OUT=$BASE_DIR/diann_final_results_cp0359
OUT_DIR=$BASE_DIR/diann_quant_files

mkdir -p $FINAL_OUT $SLURM_TMPDIR/temp

echo "copying quant to temp"
cp $QUANT_DIR/*.quant $SLURM_TMPDIR/temp/

echo "copying speclib"
cp /lustre07/scratch/jflucier/2021_long_covid/.predicted.speclib $SLURM_TMPDIR/predicted_library.speclib
cp /lustre07/scratch/jflucier/2021_long_covid/UP000005640_9606_combo.fasta $SLURM_TMPDIR/

RAW_FILES=$(ls -d $RAW_DIR/*/*.d | sed 's/^/--f /' | tr '\n' ' ')
#QUANT_FILES=$(ls $QUANT_DIR/*.d.quant | sed 's/^/--f /' | tr '\n' ' ')

export SINGULARITY_BIND=$SLURM_TMPDIR:$SLURM_TMPDIR,/net/nfs-bio/jbod2/def-gimap5/analysis/2021_long_covid:/net/nfs-bio/jbod2/def-gimap5/analysis/2021_long_covid
export APPTAINER_BIND=$SLURM_TMPDIR:$SLURM_TMPDIR,/net/nfs-bio/jbod2/def-gimap5/analysis/2021_long_covid:/net/nfs-bio/jbod2/def-gimap5/analysis/2021_long_covid

singularity exec --writable-tmpfs -e $DIANN \
diann-linux --threads 60 --verbose 2 \
"$RAW_FILES" \
--fasta "$SLURM_TMPDIR/UP000005640_9606_combo.fasta" \
--lib "$SLURM_TMPDIR/predicted_library.speclib" \
--out "$FINAL_OUT/report.tsv" \
--temp "$SLURM_TMPDIR/temp" \
--gen-spec-lib --quant-ori-names \
--qvalue 0.01 --reannotate --peptidoforms --use-quant --matrices --reanalyse \
--mass-acc 20 --mass-acc-ms1 20 --individual-mass-acc --individual-windows \
--reanalyse --rt-profiling
