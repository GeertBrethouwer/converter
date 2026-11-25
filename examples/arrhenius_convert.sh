#!/bin/bash
#
# Conversion of a large restart file on one Arrhenius node (the conversion is I/O bound; it needs only a few
# hundred MB of memory, so a login node also works for small files). Edit the command at the bottom.
# Submit from the directory with the files:  sbatch -A <your-account> $HOME/converter_and_spectra/examples/arrhenius_convert.sh
#
#SBATCH -J rsconv
#SBATCH -p gpu
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=72
#SBATCH --gpus-per-node=1
#SBATCH -t 01:00:00

RS=${RS:-$HOME/converter_and_spectra/rsconv}
cd $SLURM_SUBMIT_DIR
# binaries for the processor of this node (rebuilt by make if they were built on another type, e.g. a login node)
make -s -C $(dirname $RS) rsconv > /dev/null || { echo "building rsconv in $(dirname $RS) FAILED"; exit 1; }
$RS info data/restart

# HiCHAN -> FDCHAN, e.g. with 1.4x more z-points for the compact z-differences of FDCHAN:
srun $RS h2f data/restart restart_fdchan.bin --nz 2688 --threads $SLURM_CPUS_PER_TASK

# FDCHAN -> HiCHAN (step and dt for the HiCHAN header):
# srun $RS f2h restart.bin restart_hichan --step 0 --dt 1e-3 --grid tanh 4.5 8 --threads $SLURM_CPUS_PER_TASK

# FDCHAN single -> double precision:
# srun $RS f2f restart_sp.bin restart_dp.bin --double --threads $SLURM_CPUS_PER_TASK
