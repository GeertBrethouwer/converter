#!/bin/bash
#
# Spectra (streamwise, spanwise and 2D) of u, v, w and the scalars of a HiCHAN run on one Arrhenius node, from its
# time samples (data/samples/restart_t*, written with tout_samples or iout_samples in &output) or, if there are none,
# from its restart file (data/restart). Averaged over all samples and both channel halves; plots if python3 with
# matplotlib is available.
#
# Everything is taken from the run directory if not set: the samples, the y-grid (from the restart headers), Re
# (from the input file), Re_tau (the mean of log.dat over the time of the samples).
#
# Build once (login node):   cd $HOME && git clone https://github.com/GeertBrethouwer/converter_and_spectra.git
#                            cd converter_and_spectra && make
# Submit from any directory:
#   sbatch -A <your-account> --export=ALL,RUN=/scratch/geert/hichan_run $HOME/converter_and_spectra/examples/arrhenius_spectra_hichan.sh
#   sbatch -A <your-account> --export=ALL,RUN=/scratch/geert/hichan_run,YPLUS="10,20,50,100,300,1000",TMIN=500 \
#          $HOME/converter_and_spectra/examples/arrhenius_spectra_hichan.sh
#
#SBATCH -J spectra_hichan
#SBATCH -p gpu
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=72
#SBATCH --gpus-per-node=1
#SBATCH -t 01:00:00

# ---- settings
RUN=${RUN:-/scratch/geert/hichan_run}  # HiCHAN run directory (with the input file and data/)
DATA=${DATA:-$RUN/data}                # its data_dir
INPUT=${INPUT:-}                       # its input file (default: the *.nml in RUN)
SAMPLES=${SAMPLES:-}                   # the restart files (HiCHAN names, globs work); default: $DATA/samples/restart_t*,
                                       #   or $DATA/restart if there are no samples
TMIN=${TMIN:-}                         # use only samples with t >= TMIN (e.g. after the transient)
TMAX=${TMAX:-}                         # and t <= TMAX
YPLUS=${YPLUS:-"5,15,30,100,200"}      # positions y+ (or Y="y1,y2" for y-values, J="j1,j2" for grid indices)
Y=${Y:-}
J=${J:-}
RETAU=${RETAU:-}                       # default: the mean Re_tau of log.dat over the times of the samples
RE=${RE:-}                             # default: re of the input file (for the velocity spectra in wall units)
SYM=${SYM:-1}                          # 1: average over both channel halves
FIELDS=${FIELDS:-}                     # e.g. u,v,w,s1; empty: all
OUT=${OUT:-$RUN/spectra}
B=${B:-$HOME/converter_and_spectra}

# ---- the samples
[ -x $B/rsspec ] || { echo "$B/rsspec not found: build it (make in $B)"; exit 1; }
if [ -z "$SAMPLES" ]; then
  if ls $DATA/samples/restart_t*.hdr > /dev/null 2>&1; then SAMPLES="$DATA/samples/restart_t*"
  else SAMPLES="$DATA/restart"; fi
fi
files=""
for h in $(for p in $SAMPLES; do ls -d $p $p.hdr; done 2>/dev/null | sed -E 's/\.(bin|hdr)$//' | sort -u); do
  [ -f $h.hdr ] && [ -f $h.bin ] || { echo "skipped $h (no .hdr/.bin pair)"; continue; }
  t=$(awk 'NR==4{print $1}' $h.hdr)
  if [ -n "$TMIN" ] && awk -v t=$t -v a=$TMIN 'BEGIN{exit !(t < a)}'; then continue; fi
  if [ -n "$TMAX" ] && awk -v t=$t -v a=$TMAX 'BEGIN{exit !(t > a)}'; then continue; fi
  files="$files $h"
done
[ -n "$files" ] || { echo "no HiCHAN restart files found for $SAMPLES (TMIN=$TMIN, TMAX=$TMAX)"; exit 1; }
nfiles=$(echo $files | wc -w)
tfirst=$(awk 'NR==4{print $1}' $(echo $files | awk '{print $1}').hdr)
tlast=$(awk 'NR==4{print $1}' $(echo $files | awk '{print $NF}').hdr)
echo "$nfiles sample(s), t = $tfirst ... $tlast"

# ---- Re and Re_tau from the run
[ -z "$INPUT" ] && INPUT=$(ls $RUN/*.nml 2>/dev/null | head -1)
if [ -z "$RE" ] && [ -n "$INPUT" ]; then
  RE=$(sed -n 's/.*[^a-z_]re *= *\([0-9.eEdD+-]*\).*/\1/p' $INPUT | head -1 | tr dD eE)
fi
if [ -z "$RETAU" ] && [ -f $DATA/log.dat ]; then
  # mean of the Re_tau column (7) of log.dat over the time of the samples (the last line before tfirst if only one)
  RETAU=$(awk -v a=$tfirst -v b=$tlast '!/^#/ && $2 >= a - 1e-9 && $2 <= b + 1e-9 {s += $7; n++} !/^#/ && $2 <= a + 1e-9 {last = $7}
              END {if (n) printf "%.2f", s / n; else if (last) printf "%.2f", last}' $DATA/log.dat)
fi
[ -n "$RETAU" ] || { echo "Re_tau unknown: set RETAU"; exit 1; }
echo "Re_tau = $RETAU (from $([ -f $DATA/log.dat ] && echo "$DATA/log.dat" || echo input))${RE:+, Re = $RE}"

# ---- spectra (the y-grid from the HiCHAN headers; a grid from a file: set it with YFILE)
OPTS="--retau $RETAU --out $OUT --threads ${SLURM_CPUS_PER_TASK:-1}"
if [ -n "$J" ]; then OPTS="$OPTS --j $J"; elif [ -n "$Y" ]; then OPTS="$OPTS --y $Y"; else OPTS="$OPTS --yplus $YPLUS"; fi
gt=$(awk 'NR==5{print $1}' $(echo $files | awk '{print $1}').hdr)
if [ "$gt" = file ]; then
  YFILE=${YFILE:-$(sed -n "s/.*grid_file *= *'\([^']*\)'.*/\1/p" $INPUT | head -1)}
  [ -n "$YFILE" ] && [ "${YFILE:0:1}" != / ] && YFILE=$RUN/$YFILE
  [ -f "$YFILE" ] || { echo "grid_type = file: set YFILE to the grid file"; exit 1; }
fi
[ -n "$YFILE" ] && OPTS="$OPTS --yfile $YFILE"
[ "$SYM" = 1 ] && OPTS="$OPTS --sym"
[ -n "$FIELDS" ] && OPTS="$OPTS --fields $FIELDS"
echo "== $B/rsspec $OPTS <$nfiles files>"
start=$(date +%s)
srun $B/rsspec $OPTS $files || { echo "rsspec FAILED"; exit 1; }
echo "spectra took $(( $(date +%s) - start )) s"
cat $OUT/summary.dat

# ---- plots
PL="--retau $RETAU"
[ -n "$RE" ] && PL="$PL --utau $(awk -v a=$RETAU -v b=$RE 'BEGIN{printf "%.8g", a / b}')"
[ -n "$FIELDS" ] && PL="$PL --fields $FIELDS"
if python3 -c "import matplotlib" 2>/dev/null; then
  python3 $B/tools/plot_spectra.py $OUT $PL
else
  echo "no python3 with matplotlib here: plot later with  python3 $B/tools/plot_spectra.py $OUT $PL"
fi
