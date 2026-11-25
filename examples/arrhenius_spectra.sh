#!/bin/bash
#
# Spectra (streamwise, spanwise and 2D) of u, v, w and the scalars at chosen wall-normal positions from HiCHAN or
# FDCHAN restart files on one Arrhenius node, averaged over all given files (time samples) and over both channel
# halves, and plots of them (if python3 with matplotlib is available).
#
# Build once (login node):   cd $HOME && git clone https://github.com/GeertBrethouwer/converter_and_spectra.git
#                            cd converter_and_spectra && make
# Set the variables below (or pass them on the command line), then submit from any directory:
#   sbatch -A <your-account> $HOME/converter_and_spectra/examples/arrhenius_spectra.sh
#   sbatch -A <your-account> --export=ALL,FILES="/scratch/geert/run/samples/restart_*.bin",RETAU=1000,GRID="tanh 4.5" \
#          $HOME/converter_and_spectra/examples/arrhenius_spectra.sh
#
# Time samples: FDCHAN and HiCHAN (&output) write numbered restart files with tout_samples (or iout_samples), e.g.
# tout_samples = 5 -> data/samples/restart_t<time>.bin (HiCHAN also .hdr); FILES="<run>/data/samples/restart_t*.bin".
#
#SBATCH -J spectra
#SBATCH -p gpu
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=72
#SBATCH --gpus-per-node=1
#SBATCH -t 01:00:00

# ---- settings
FILES=${FILES:-"/scratch/geert/fdchan_run/data/samples/restart_t*.bin"}  # restart files (FDCHAN *.bin or HiCHAN <name>); globs work
RETAU=${RETAU:-1000}            # friction Reynolds number: positions in wall units, wavelengths in wall units in the plots
YPLUS=${YPLUS:-"5,15,30,100,200"}  # positions y+ from the wall (or set Y="y1,y2" for y-values, J="j1,j2" for indices)
Y=${Y:-}
J=${J:-}
GRID=${GRID:-"tanh 4.5"}        # y-grid of the run: type and stretching (taken from a HiCHAN header if empty)
SYM=${SYM:-1}                   # 1: average over both channel halves
FIELDS=${FIELDS:-}              # e.g. u,v,w,s1; empty: all
RE=${RE:-}                      # bulk Reynolds number (U h/nu): velocity spectra in wall units in the plots (u_tau = RETAU/RE)
OUT=${OUT:-$SLURM_SUBMIT_DIR/spectra}
B=${B:-$HOME/converter_and_spectra}

# ---- run
# binaries for the processor of this node (rebuilt by make if they were built on another type, e.g. a login node)
make -s -C $B rsspec > /dev/null || { echo "building rsspec in $B FAILED"; exit 1; }
[ -x $B/rsspec ] || { echo "$B/rsspec not found: build it (make in $B)"; exit 1; }
files=$(ls -d $FILES 2>/dev/null | sed 's/\.hdr$//' | sort -u)
[ -n "$files" ] || { echo "no files match $FILES"; exit 1; }
OPTS="--retau $RETAU --out $OUT --threads $SLURM_CPUS_PER_TASK"
if [ -n "$J" ]; then OPTS="$OPTS --j $J"; elif [ -n "$Y" ]; then OPTS="$OPTS --y $Y"; else OPTS="$OPTS --yplus $YPLUS"; fi
[ -n "$GRID" ] && OPTS="$OPTS --grid $GRID"
[ "$SYM" = 1 ] && OPTS="$OPTS --sym"
[ -n "$FIELDS" ] && OPTS="$OPTS --fields $FIELDS"
echo "== $B/rsspec $OPTS" $files
start=$(date +%s)
srun $B/rsspec $OPTS $files || { echo "rsspec FAILED"; exit 1; }
echo "spectra took $(( $(date +%s) - start )) s"
cat $OUT/summary.dat

# ---- plots
PL="--retau $RETAU"
[ -n "$RE" ] && PL="$PL --utau $(python3 -c "print($RETAU/$RE)")"
[ -n "$FIELDS" ] && PL="$PL --fields $FIELDS"
if python3 -c "import matplotlib" 2>/dev/null; then
  python3 $B/tools/plot_spectra.py $OUT $PL
else
  echo "no python3 with matplotlib here: plot later with  python3 $B/tools/plot_spectra.py $OUT $PL"
fi
