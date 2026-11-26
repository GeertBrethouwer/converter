#!/bin/bash
#
# 2D slices (xz-planes at chosen wall distances, an xy- and a yz-plane) of u, v, w and the scalars of an FDCHAN run
# on one Arrhenius node, for visualisation: from its latest time sample (data/samples/restart_t*.bin), from all
# samples (e.g. frames of a film) or from its restart file (data/restart.bin); single or double precision. Writes
# RUN/slices/<file>_slices.mat (MATLAB / Octave: tools/plot_slices.m; .npz without scipy) and png pictures
# (tools/plot_slices.py, if matplotlib is there).
#
# Everything is taken from the run directory if not set: the files, data_dir and the y-grid (input file; FDCHAN
# restart files do not contain the grid), Re_tau (log.dat at the time of each file, for YPLUS).
#
# Get the code once (login node):   cd $HOME && git clone https://github.com/GeertBrethouwer/converter_and_spectra.git
# Needs python3 with numpy (scipy for .mat files, matplotlib for the pictures).
# Submit from any directory:
#   sbatch -A <your-account> --export=ALL,RUN=/scratch/geert/fdchan_run $HOME/converter_and_spectra/examples/arrhenius_slices_fdchan.sh
#   sbatch -A <your-account> --export=ALL,RUN=/scratch/geert/fdchan_run,WHICH=all,TMIN=500,YPLUS="15 100",UPSAMPLE=2 \
#          $HOME/converter_and_spectra/examples/arrhenius_slices_fdchan.sh
#
#SBATCH -J slices_fdchan
#SBATCH -p gpu
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gpus-per-node=1
#SBATCH -t 02:00:00

# ---- settings
RUN=${RUN:-/scratch/geert/fdchan_run}  # FDCHAN run directory (with the input file and the data directory)
INPUT=${INPUT:-$RUN/input.txt}         # its input file
DATA=${DATA:-}                         # its data_dir (default: data_dir of the input file, else RUN/data)
GRID=${GRID:-}                         # y-grid "type alpha" (default: grid_type, grid_alpha of the input, else tanh 3.0)
WHICH=${WHICH:-last}                   # last: the latest sample (else data/restart.bin); all: every sample (TMIN,
                                       #   TMAX); restart: data/restart.bin
SAMPLES=${SAMPLES:-}                   # instead: these restart files (globs work)
TMIN=${TMIN:-}                         # WHICH=all: only samples with TMIN <= t <= TMAX
TMAX=${TMAX:-}
YPLUS=${YPLUS:-}                       # xz-planes at these y+ (space separated), e.g. "5 15 100" (Re_tau from log.dat)
Y=${Y:-"-0.9 0"}                       # else at these y (-1 bottom wall, 0 centre)
RETAU=${RETAU:-}                       # for YPLUS; default: from log.dat at the time of each file
X=${X:-0}                              # x of the yz-plane
Z=${Z:-0}                              # z of the xy-plane
UPSAMPLE=${UPSAMPLE:-1}                # finer x-z output grid (exact Fourier interpolation), e.g. 2
FLUCT=${FLUCT:-1}                      # 1: fluctuations (the x-z mean removed)
FIELDS=${FIELDS:-}                     # e.g. "u v s1"; empty: all
CLIM=${CLIM:-}                         # fixed colour range of the pictures, e.g. "-0.3 0.3" (the same for all frames)
FORMAT=${FORMAT:-mat}                  # mat (needs scipy) or npz
OUTDIR=${OUTDIR:-$RUN/slices}
PLOT=${PLOT:-1}                        # 1: png pictures
SCM=${SCM:-}                           # directory of the Scientific Colour Maps (default: $SCIENTIFIC_COLOURMAPS,
                                       #   ~/scientific_colourmaps; tools/install_scientific_colourmaps.sh)
CMAP=${CMAP:-}                         # colour map for fields of one sign (default batlow)
CMAP_DIV=${CMAP_DIV:-}                 # colour map for signed fields (default vik)
B=${B:-$HOME/converter_and_spectra}

# ---- checks
python3 -c "import numpy" 2>/dev/null || { echo "python3 with numpy is needed (module load ..., or pip install --user numpy)"; exit 1; }
if [ "$FORMAT" = mat ] && ! python3 -c "import scipy" 2>/dev/null; then
  echo "no scipy here: writing .npz files (plot with tools/plot_slices.py)"; FORMAT=npz
fi
plot=0
if [ "$PLOT" = 1 ]; then
  python3 -c "import matplotlib" 2>/dev/null && plot=1 || echo "no matplotlib here: no pictures (tools/plot_slices.py later)"
fi
mkdir -p $OUTDIR

# ---- data_dir and the y-grid from the input file
item() { [ -f "$INPUT" ] && sed -n "s/^[[:space:]]*$1[[:space:]]*=[[:space:]]*\([^!]*\).*/\1/p" $INPUT | tail -1 | sed 's/[[:space:]]*$//'; }
[ -f "$INPUT" ] || echo "warning: no input file $INPUT: the grid from GRID or the default (tanh 3.0)"
if [ -z "$DATA" ]; then
  dd=$(item data_dir)
  if [ -z "$dd" ]; then DATA=$RUN/data; elif [ "${dd:0:1}" = / ]; then DATA=$dd; else DATA=$RUN/$dd; fi
fi
if [ -z "$GRID" ]; then gt=$(item grid_type); ga=$(item grid_alpha); GRID="${gt:-tanh} ${ga:-3.0}"; fi
echo "y-grid $GRID"

# ---- the files
tof() { od -An -t f8 -j 32 -N 8 $1 | tr -d ' '; }  # the time in the FDCHAN header (bytes 32..39)
if [ -z "$SAMPLES" ]; then
  case $WHICH in
    last)    SAMPLES=$(ls $DATA/samples/restart_t*.bin 2>/dev/null | sort | tail -1)
             [ -z "$SAMPLES" ] && SAMPLES=$DATA/restart.bin ;;
    all)     SAMPLES="$DATA/samples/restart_t*.bin" ;;
    restart) SAMPLES=$DATA/restart.bin ;;
    *)       echo "WHICH must be last, all or restart"; exit 1 ;;
  esac
fi
files=""
for h in $(ls -d $SAMPLES 2>/dev/null | sort -u); do
  [ -f $h ] || continue
  t=$(tof $h)
  if [ -n "$TMIN" ] && awk -v t=$t -v a=$TMIN 'BEGIN{exit !(t < a)}'; then continue; fi
  if [ -n "$TMAX" ] && awk -v t=$t -v a=$TMAX 'BEGIN{exit !(t > a)}'; then continue; fi
  files="$files $h"
done
[ -n "$files" ] || { echo "no FDCHAN restart files found for $SAMPLES (TMIN=$TMIN, TMAX=$TMAX)"; exit 1; }
echo "$(echo $files | wc -w) file(s)"

GOPT="--grid $GRID"

# ---- slices, one file after the other
start=$(date +%s)
for h in $files; do
  t=$(tof $h)
  OPTS="--x $X --z $Z --upsample $UPSAMPLE $GOPT"
  if [ -n "$YPLUS" ]; then
    rt=$RETAU
    if [ -z "$rt" ]; then  # Re_tau of the last line of log.dat at or before t
      rt=$(awk -v t=$t '!/^#/ && $2 <= t + 1e-9 {r = $7} END {if (r) printf "%.2f", r}' $DATA/log.dat 2>/dev/null)
    fi
    [ -n "$rt" ] || { echo "YPLUS needs RETAU (no log.dat in $DATA)"; exit 1; }
    OPTS="$OPTS --yplus $YPLUS --retau $rt"
  else
    OPTS="$OPTS --y $Y"
  fi
  [ "$FLUCT" = 1 ] && OPTS="$OPTS --fluct"
  [ -n "$FIELDS" ] && OPTS="$OPTS --fields $FIELDS"
  out=$OUTDIR/$(basename $h .bin)_slices.$FORMAT
  echo "== $(basename $h .bin) (t = $t)${rt:+, Re_tau = $rt}"
  srun python3 $B/tools/slices.py $h $OPTS -o $out | sed 's/^/   /' || { echo "slices FAILED for $h"; exit 1; }
  if [ $plot = 1 ]; then
    python3 $B/tools/plot_slices.py $out ${FIELDS:+--fields $FIELDS} ${SCM:+--cmap-dir $SCM} ${CMAP:+--cmap $CMAP} ${CMAP_DIV:+--cmap-div $CMAP_DIV} ${CLIM:+--clim $CLIM} | sed 's/^/   /'
  fi
done
echo "slices of $(echo $files | wc -w) file(s) in $OUTDIR, $(( $(date +%s) - start )) s"
