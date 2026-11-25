#!/bin/bash
#
# Spectra (streamwise, spanwise and 2D) of u, v, w and the scalars of an FDCHAN run on one Arrhenius node, from its
# time samples (data/samples/restart_t*.bin, written with tout_samples or iout_samples) or, if there are none, from
# its restart file (data/restart.bin). Single or double precision. Averaged over all samples and both channel
# halves; plots if python3 with matplotlib is available.
#
# Everything is taken from the run directory if not set: the samples, the y-grid and Re (from the input file; FDCHAN
# restart files do not contain the grid), Re_tau (the mean of log.dat over the time of the samples).
#
# Build once (login node):   cd $HOME && git clone https://github.com/GeertBrethouwer/converter_and_spectra.git
#                            cd converter_and_spectra && make
# Submit from any directory:
#   sbatch -A <your-account> --export=ALL,RUN=/scratch/geert/fdchan_run $HOME/converter_and_spectra/examples/arrhenius_spectra_fdchan.sh
#   sbatch -A <your-account> --export=ALL,RUN=/scratch/geert/fdchan_run,YPLUS="10,20,50,100,300,1000",TMIN=500 \
#          $HOME/converter_and_spectra/examples/arrhenius_spectra_fdchan.sh
#
#SBATCH -J spectra_fdchan
#SBATCH -p gpu
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=72
#SBATCH --gpus-per-node=1
#SBATCH -t 01:00:00

# ---- settings
RUN=${RUN:-/scratch/geert/fdchan_run}  # FDCHAN run directory (with the input file and the data directory)
INPUT=${INPUT:-$RUN/input.txt}         # its input file
DATA=${DATA:-}                         # its data_dir (default: data_dir of the input file, else RUN/data)
SAMPLES=${SAMPLES:-}                   # the restart files (globs work); default: $DATA/samples/restart_t*.bin, or
                                       #   $DATA/restart.bin if there are no samples
TMIN=${TMIN:-}                         # use only samples with t >= TMIN (e.g. after the transient)
TMAX=${TMAX:-}                         # and t <= TMAX
YPLUS=${YPLUS:-"5,15,30,100,200"}      # positions y+ (or Y="y1,y2" for y-values, J="j1,j2" for grid indices)
Y=${Y:-}
J=${J:-}
RETAU=${RETAU:-}                       # default: the mean Re_tau of log.dat over the times of the samples
RE=${RE:-}                             # default: re of the input file (2800 if not set there, as FDCHAN)
GRID=${GRID:-}                         # "type alpha"; default: grid_type, grid_alpha of the input (tanh 3.0 if not set)
SYM=${SYM:-1}                          # 1: average over both channel halves
FIELDS=${FIELDS:-}                     # e.g. u,v,w,s1; empty: all
OUT=${OUT:-$RUN/spectra}
B=${B:-$HOME/converter_and_spectra}

[ -x $B/rsspec ] || { echo "$B/rsspec not found: build it (make in $B)"; exit 1; }
# value of an item of the FDCHAN input file ("key = value", ! comments)
item() { [ -f "$INPUT" ] && sed -n "s/^[[:space:]]*$1[[:space:]]*=[[:space:]]*\([^!]*\).*/\1/p" $INPUT | tail -1 | sed 's/[[:space:]]*$//'; }
[ -f "$INPUT" ] || echo "warning: no input file $INPUT: grid and Re from the defaults or GRID, RE"

# ---- the samples
if [ -z "$DATA" ]; then
  dd=$(item data_dir)
  if [ -z "$dd" ]; then DATA=$RUN/data; elif [ "${dd:0:1}" = / ]; then DATA=$dd; else DATA=$RUN/$dd; fi
fi
if [ -z "$SAMPLES" ]; then
  if ls $DATA/samples/restart_t*.bin > /dev/null 2>&1; then SAMPLES="$DATA/samples/restart_t*.bin"
  else SAMPLES="$DATA/restart.bin"; fi
fi
tof() { od -An -t f8 -j 32 -N 8 $1 | tr -d ' '; }  # the time in the FDCHAN header (bytes 32..39)
files=""
for f in $(ls -d $SAMPLES 2>/dev/null | sort -u); do
  [ -f $f ] || continue
  t=$(tof $f)
  if [ -n "$TMIN" ] && awk -v t=$t -v a=$TMIN 'BEGIN{exit !(t < a)}'; then continue; fi
  if [ -n "$TMAX" ] && awk -v t=$t -v a=$TMAX 'BEGIN{exit !(t > a)}'; then continue; fi
  files="$files $f"
done
[ -n "$files" ] || { echo "no FDCHAN restart files found for $SAMPLES (TMIN=$TMIN, TMAX=$TMAX)"; exit 1; }
nfiles=$(echo $files | wc -w)
tfirst=$(tof $(echo $files | awk '{print $1}'))
tlast=$(tof $(echo $files | awk '{print $NF}'))
echo "$nfiles sample(s), t = $tfirst ... $tlast"

# ---- grid, Re and Re_tau from the run
if [ -z "$GRID" ]; then
  gt=$(item grid_type); ga=$(item grid_alpha)
  GRID="${gt:-tanh} ${ga:-3.0}"
fi
[ -z "$RE" ] && RE=$(item re) && RE=${RE:-2800}
if [ -z "$RETAU" ] && [ -f $DATA/log.dat ]; then
  # mean of the Re_tau column (7) of log.dat over the time of the samples (the last line before tfirst if only one)
  RETAU=$(awk -v a=$tfirst -v b=$tlast '!/^#/ && $2 >= a - 1e-9 && $2 <= b + 1e-9 {s += $7; n++} !/^#/ && $2 <= a + 1e-9 {last = $7}
              END {if (n) printf "%.2f", s / n; else if (last) printf "%.2f", last}' $DATA/log.dat)
fi
[ -n "$RETAU" ] || { echo "Re_tau unknown: set RETAU"; exit 1; }
echo "grid $GRID, Re = $RE, Re_tau = $RETAU"

# ---- spectra
OPTS="--retau $RETAU --grid $GRID --out $OUT --threads ${SLURM_CPUS_PER_TASK:-1}"
if [ -n "$J" ]; then OPTS="$OPTS --j $J"; elif [ -n "$Y" ]; then OPTS="$OPTS --y $Y"; else OPTS="$OPTS --yplus $YPLUS"; fi
[ "$SYM" = 1 ] && OPTS="$OPTS --sym"
[ -n "$FIELDS" ] && OPTS="$OPTS --fields $FIELDS"
echo "== $B/rsspec $OPTS <$nfiles files>"
start=$(date +%s)
srun $B/rsspec $OPTS $files || { echo "rsspec FAILED"; exit 1; }
echo "spectra took $(( $(date +%s) - start )) s"
cat $OUT/summary.dat

# ---- plots
PL="--retau $RETAU --utau $(awk -v a=$RETAU -v b=$RE 'BEGIN{printf "%.8g", a / b}')"
[ -n "$FIELDS" ] && PL="$PL --fields $FIELDS"
if python3 -c "import matplotlib" 2>/dev/null; then
  python3 $B/tools/plot_spectra.py $OUT $PL
else
  echo "no python3 with matplotlib here: plot later with  python3 $B/tools/plot_spectra.py $OUT $PL"
fi
