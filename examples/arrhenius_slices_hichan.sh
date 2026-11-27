#!/bin/bash
#
# 2D slices (xz-planes at chosen wall distances, an xy- and a yz-plane) of u, v, w and the scalars of a HiCHAN run
# on one Arrhenius node, for visualisation: from its latest time sample (data/samples/restart_t*), from all samples
# (e.g. frames of a film) or from its restart file (data/restart). Writes RUN/slices/<file>_slices.mat (MATLAB /
# Octave: tools/plot_slices.m; .npz without scipy) and png pictures (tools/plot_slices.py, if matplotlib is there).
#
# Everything is taken from the run directory if not set: the files, the y-grid (restart headers), Re_tau (log.dat
# at the time of each file, for YPLUS).
#
# Get the code once (login node):   cd $HOME && git clone https://github.com/GeertBrethouwer/converter_and_spectra.git
# Needs python3 with numpy (scipy for .mat files, matplotlib for the pictures).
# Submit from any directory:
#   sbatch -A <your-account> --export=ALL,RUN=/scratch/geert/hichan_run $HOME/converter_and_spectra/examples/arrhenius_slices_hichan.sh
#   sbatch -A <your-account> --export=ALL,RUN=/scratch/geert/hichan_run,WHICH=all,TMIN=500,YPLUS="15 100",UPSAMPLE=2 \
#          $HOME/converter_and_spectra/examples/arrhenius_slices_hichan.sh
#
#SBATCH -J slices_hichan
#SBATCH -p gpu
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gpus-per-node=1
#SBATCH -t 02:00:00

# ---- settings
RUN=${RUN:-/scratch/geert/hichan_run}  # HiCHAN run directory (with the input file and data/)
DATA=${DATA:-$RUN/data}                # its data_dir
INPUT=${INPUT:-}                       # its input file (default: the *.nml in RUN; for a grid_type = 'file' grid)
WHICH=${WHICH:-last}                   # last: the latest sample (else data/restart); all: every sample (TMIN, TMAX);
                                       #   restart: data/restart
SAMPLES=${SAMPLES:-}                   # instead: these restart files (HiCHAN names, globs work)
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

# ---- the files
tof() { awk 'NR==4{print $1}' $1.hdr; }
if [ -z "$SAMPLES" ]; then
  case $WHICH in
    last)    SAMPLES=$(ls $DATA/samples/restart_t*.hdr 2>/dev/null | sort | tail -1)
             [ -z "$SAMPLES" ] && SAMPLES=$DATA/restart.hdr ;;
    all)     SAMPLES="$DATA/samples/restart_t*" ;;
    restart) SAMPLES=$DATA/restart.hdr ;;
    *)       echo "WHICH must be last, all or restart"; exit 1 ;;
  esac
fi
files=""
for h in $(for p in $SAMPLES; do ls -d $p $p.hdr; done 2>/dev/null | sed -E 's/\.(bin|hdr|stats)$//' | sort -u); do
  [ -f $h.hdr ] && [ -f $h.bin ] || continue
  t=$(tof $h)
  if [ -n "$TMIN" ] && awk -v t=$t -v a=$TMIN 'BEGIN{exit !(t < a)}'; then continue; fi
  if [ -n "$TMAX" ] && awk -v t=$t -v a=$TMAX 'BEGIN{exit !(t > a)}'; then continue; fi
  files="$files $h"
done
[ -n "$files" ] || { echo "no HiCHAN restart files found for $SAMPLES (TMIN=$TMIN, TMAX=$TMAX)"; exit 1; }
echo "$(echo $files | wc -w) file(s)"

# a grid from a file
gt=$(awk 'NR==5{print $1}' $(echo $files | awk '{print $1}').hdr)
GOPT=""
if [ "$gt" = file ]; then
  [ -z "$INPUT" ] && INPUT=$(ls $RUN/*.nml 2>/dev/null | head -1)
  YFILE=${YFILE:-$(sed -n "s/.*grid_file *= *'\([^']*\)'.*/\1/p" $INPUT | head -1)}
  [ -n "$YFILE" ] && [ "${YFILE:0:1}" != / ] && YFILE=$RUN/$YFILE
  [ -f "$YFILE" ] || { echo "grid_type = file: set YFILE to the grid file"; exit 1; }
  GOPT="--grid-file $YFILE"
fi

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
  out=$OUTDIR/$(basename $h)_slices.$FORMAT
  echo "== $(basename $h) (t = $t)${rt:+, Re_tau = $rt}"
  srun python3 $B/tools/slices.py $h $OPTS -o $out | sed 's/^/   /' || { echo "slices FAILED for $h"; exit 1; }
  if [ $plot = 1 ]; then
    python3 $B/tools/plot_slices.py $out ${FIELDS:+--fields $FIELDS} ${SCM:+--cmap-dir $SCM} ${CMAP:+--cmap $CMAP} ${CMAP_DIV:+--cmap-div $CMAP_DIV} ${CLIM:+--clim $CLIM} | sed 's/^/   /'
  fi
done
echo "slices of $(echo $files | wc -w) file(s) in $OUTDIR, $(( $(date +%s) - start )) s"
