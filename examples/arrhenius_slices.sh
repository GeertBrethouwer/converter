#!/bin/bash
#
# 2D slices (xz-planes at chosen wall distances, an xy- and a yz-plane) of u, v, w and the scalars from a HiCHAN or
# FDCHAN restart file (single or double precision) on one Arrhenius node, for visualisation: a .mat file (MATLAB /
# Octave: tools/plot_slices.m; needs scipy) or a .npz file, and clean png pictures of every plane with
# tools/plot_slices.py (matplotlib; Scientific Colour Maps, tools/install_scientific_colourmaps.sh).
#
# Get the code once (login node):   cd $HOME && git clone https://github.com/GeertBrethouwer/converter_and_spectra.git
# Needs python3 with numpy (scipy for .mat files, matplotlib for the pictures).
# Submit from any directory:
#   sbatch -A <your-account> --export=ALL,FILE=/scratch/geert/fdchan_run/data/restart.bin,RUN=/scratch/geert/fdchan_run \
#          $HOME/converter_and_spectra/examples/arrhenius_slices.sh
#   sbatch -A <your-account> --export=ALL,FILE=/scratch/geert/hichan_run/data/samples/restart_t000500.0012,YPLUS="15 100",RETAU=1000 \
#          $HOME/converter_and_spectra/examples/arrhenius_slices.sh
#
#SBATCH -J slices
#SBATCH -p gpu
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gpus-per-node=1
#SBATCH -t 01:00:00

# ---- settings
FILE=${FILE:-/scratch/geert/fdchan_run/data/restart.bin}  # FDCHAN restart (.bin) or HiCHAN restart (name, .bin or .hdr)
RUN=${RUN:-}            # optional run directory: the FDCHAN grid (input.txt) and Re_tau (data/log.dat) from it
GRID=${GRID:-}          # FDCHAN y-grid "type alpha" (default: grid_type, grid_alpha of RUN/input.txt, else tanh 3.0);
                        #   HiCHAN files contain it
YPLUS=${YPLUS:-}        # xz-planes at these y+ (space separated, needs RETAU or RUN), e.g. "5 15 100"
Y=${Y:-"-0.9 0"}        # else at these y (-1 bottom wall, 0 centre)
RETAU=${RETAU:-}        # for YPLUS; default: the last Re_tau in RUN/data/log.dat
X=${X:-0}               # x of the yz-plane
Z=${Z:-0}               # z of the xy-plane
UPSAMPLE=${UPSAMPLE:-1} # finer x-z output grid (exact Fourier interpolation), e.g. 2
FLUCT=${FLUCT:-1}       # 1: fluctuations (the x-z mean removed)
FIELDS=${FIELDS:-}      # e.g. "u v s1"; empty: all
OUT=${OUT:-}            # output file .mat or .npz; default: <restart>_slices.mat (.npz without scipy)
PLOT=${PLOT:-1}         # 1: png pictures with tools/plot_slices.py
SCM=${SCM:-}            # directory of the Scientific Colour Maps (default: $SCIENTIFIC_COLOURMAPS, ~/scientific_colourmaps)
CMAP=${CMAP:-}          # colour map for fields of one sign (default batlow), CMAP_DIV for signed fields (default vik)
CMAP_DIV=${CMAP_DIV:-}
CLIM=${CLIM:-}          # fixed colour range, e.g. "-0.3 0.3"
B=${B:-$HOME/converter_and_spectra}

# ---- checks
python3 -c "import numpy" 2>/dev/null || { echo "python3 with numpy is needed (module load ..., or pip install --user numpy)"; exit 1; }
F=$FILE
[ -e $F ] || [ -e $F.hdr ] || [ -e $F.bin ] || { echo "not found: $FILE"; exit 1; }
hichan=0
for h in ${F%.*}.hdr $F.hdr; do [ -f $h ] && hichan=1; done
OPTS="--x $X --z $Z --upsample $UPSAMPLE"
if [ $hichan = 0 ]; then  # FDCHAN: the y-grid from GRID, RUN/input.txt or the default
  item() { [ -f "$RUN/input.txt" ] && sed -n "s/^[[:space:]]*$1[[:space:]]*=[[:space:]]*\([^!]*\).*/\1/p" $RUN/input.txt | tail -1 | sed 's/[[:space:]]*$//'; }
  if [ -z "$GRID" ]; then gt=$(item grid_type); ga=$(item grid_alpha); GRID="${gt:-tanh} ${ga:-3.0}"; fi
  OPTS="$OPTS --grid $GRID"
  echo "FDCHAN file, y-grid $GRID"
else
  echo "HiCHAN file (the y-grid from its header)"
fi
if [ -n "$YPLUS" ]; then
  if [ -z "$RETAU" ] && [ -n "$RUN" ]; then
    RETAU=$(awk '!/^#/ {r = $7} END {if (r) printf "%.2f", r}' $RUN/data/log.dat 2>/dev/null)
  fi
  [ -n "$RETAU" ] || { echo "YPLUS needs RETAU (or RUN with data/log.dat)"; exit 1; }
  OPTS="$OPTS --yplus $YPLUS --retau $RETAU"
else
  OPTS="$OPTS --y $Y"
fi
[ "$FLUCT" = 1 ] && OPTS="$OPTS --fluct"
[ -n "$FIELDS" ] && OPTS="$OPTS --fields $FIELDS"
[ -n "$OUT" ] && OPTS="$OPTS -o $OUT"

# ---- slices
echo "== python3 $B/tools/slices.py $F $OPTS"
start=$(date +%s)
log=$(mktemp)
srun python3 $B/tools/slices.py $F $OPTS | tee $log || { echo "slices FAILED"; exit 1; }
echo "slices took $(( $(date +%s) - start )) s"
outf=$(awk '/^wrote/ {print $2}' $log); rm -f $log

# ---- pictures
if [ "$PLOT" = 1 ] && [ -n "$outf" ]; then
  if python3 -c "import matplotlib" 2>/dev/null; then
    python3 $B/tools/plot_slices.py $outf ${FIELDS:+--fields $FIELDS} ${SCM:+--cmap-dir $SCM} ${CMAP:+--cmap $CMAP} ${CMAP_DIV:+--cmap-div $CMAP_DIV} ${CLIM:+--clim $CLIM}
  else
    echo "no python3 with matplotlib here: pictures later with  python3 $B/tools/plot_slices.py $outf"
  fi
fi
