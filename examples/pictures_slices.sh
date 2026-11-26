#!/bin/bash
#
# Clean pictures of slice files (.mat or .npz from tools/slices.py or the slice job scripts) with the Scientific
# Colour Maps of Fabio Crameri: every plane as a picture without colour bar, axes or margins (tools/plot_slices.py).
# For several files (e.g. the frames of WHICH=all) every field gets the same colour range in all of them, so the
# frames can be compared or put together into a film.
#
# Runs on a desktop:        bash examples/pictures_slices.sh                       (settings as variables, see below)
#   SLICES="/data/run/slices/*_slices.mat" CMAP_DIV=vik bash $HOME/converter_and_spectra/examples/pictures_slices.sh
# or as a job on Arrhenius: sbatch -A <your-account> --export=ALL,SLICES="/scratch/geert/fdchan_run/slices/*.mat" \
#                                  $HOME/converter_and_spectra/examples/pictures_slices.sh
# Needs python3 with numpy and matplotlib (scipy for .mat files). The colour maps are installed into SCM first if
# they are not there (tools/install_scientific_colourmaps.sh; needs internet, i.e. a desktop or a login node).
#
#SBATCH -J pictures
#SBATCH -p gpu
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gpus-per-node=1
#SBATCH -t 01:00:00

# ---- settings
SLICES=${SLICES:-"slices/*_slices.mat slices/*_slices.npz"}  # slice files (globs work) or a directory
SCM=${SCM:-/scratch/geert/scientific_colourmaps}  # the Scientific Colour Maps (installed here if missing)
CMAP_DIV=${CMAP_DIV:-vik}      # map for fields with both signs (fluctuations, v, w); e.g. vik, roma, berlin, broc, cork
CMAP=${CMAP:-batlow}           # map for fields of one sign; e.g. batlow, lajolla, oslo, hawaii (reversed: name_r)
FIELDS=${FIELDS:-}             # e.g. "u v s1"; empty: all
PLANE=${PLANE:-all}            # all, xz, xy or yz
COMMON=${COMMON:-1}            # 1: the same colour range per field in all files (symmetric for signed fields)
CLIM=${CLIM:-}                 # or a fixed range for all fields, e.g. "-0.3 0.3"
SIZE=${SIZE:-2000}             # pixels of the longer side of each picture
SHADING=${SHADING:-gouraud}    # gouraud (smooth) or nearest (one colour per grid point)
OUTDIR=${OUTDIR:-}             # where the pictures go (default: next to each slice file)
B=${B:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." 2>/dev/null && pwd)}  # the repository (default: this script's)
[ -f "$B/tools/plot_slices.py" ] || B=$HOME/converter_and_spectra

# ---- checks
[ -f $B/tools/plot_slices.py ] || { echo "tools/plot_slices.py not found in $B: set B to the repository"; exit 1; }
python3 -c "import numpy, matplotlib" 2>/dev/null || { echo "python3 with numpy and matplotlib is needed"; exit 1; }
[ -d "$SLICES" ] && SLICES="$SLICES/*_slices.mat $SLICES/*_slices.npz"
files=$(ls -d $SLICES 2>/dev/null | sort -u)
[ -n "$files" ] || { echo "no slice files found for $SLICES"; exit 1; }
if echo "$files" | grep -q '\.mat$' && ! python3 -c "import scipy" 2>/dev/null; then
  echo "scipy is needed for .mat files (pip install --user scipy)"; exit 1
fi
echo "$(echo $files | wc -w) slice file(s)"

# ---- the colour maps (install if missing)
have_maps() {  # the map $1 in one of the layouts (cmcrameri, the official collection, plain files)
  local f
  for f in $SCM/cmcrameri/cmaps/$1.txt $SCM/*/$1/$1.txt $SCM/$1/$1.txt $SCM/$1.txt; do [ -f "$f" ] && return 0; done
  return 1
}
m1=${CMAP_DIV%_r}; m2=${CMAP%_r}
if ! have_maps $m1 || ! have_maps $m2; then
  echo "== the Scientific Colour Maps are not in $SCM: installing them"
  bash $B/tools/install_scientific_colourmaps.sh $SCM || { echo "installation failed: install them on a machine with internet"; exit 1; }
fi
have_maps $m1 && have_maps $m2 || { echo "colour map $CMAP_DIV or $CMAP not found in $SCM"; exit 1; }
echo "colour maps: $CMAP_DIV (signed fields), $CMAP (others), from $SCM"

# ---- the colour range of every field over all files
declare -A LO HI
if [ -z "$CLIM" ] && [ "$COMMON" = 1 ] && [ $(echo $files | wc -w) -gt 1 ]; then
  while read f lo hi; do LO[$f]=$lo; HI[$f]=$hi; done < <(python3 - $PLANE "$FIELDS" $files <<'EOF'
import sys
import numpy as np
plane, sel, files = sys.argv[1], sys.argv[2].split(), sys.argv[3:]
rng = {}
for fn in files:
    if fn.endswith(".npz"):
        S = dict(np.load(fn))
    else:
        from scipy.io import loadmat
        S = loadmat(fn, squeeze_me=True)
    for f in [str(x).strip() for x in np.ravel(S["fields"])]:
        if sel and f not in sel:
            continue
        parts = [S[p + f] for p in ("xz_", "xy_", "yz_") if plane in ("all", p[:2])]
        lo, hi = min(float(np.min(p)) for p in parts), max(float(np.max(p)) for p in parts)
        a, b = rng.get(f, (lo, hi))
        rng[f] = (min(a, lo), max(b, hi))
for f, (lo, hi) in rng.items():
    if lo < 0 < hi:
        lo, hi = -max(-lo, hi), max(-lo, hi)
    print(f, lo, hi)
EOF
)
  echo "common colour ranges over all files:"
  for f in "${!LO[@]}"; do printf "  %-4s %12.5e .. %12.5e\n" $f ${LO[$f]} ${HI[$f]}; done
fi

# ---- pictures
start=$(date +%s)
for fn in $files; do
  src=$fn
  if [ -n "$OUTDIR" ]; then  # pictures in OUTDIR: plot a link to the file there
    mkdir -p $OUTDIR; src=$OUTDIR/$(basename $fn); [ -e $src ] || ln -s $(cd $(dirname $fn) && pwd)/$(basename $fn) $src
  fi
  OPTS="--cmap-dir $SCM --cmap $CMAP --cmap-div $CMAP_DIV --plane $PLANE --size $SIZE --shading $SHADING"
  echo "== $(basename $fn)"
  if [ ${#LO[@]} -gt 0 ]; then  # field by field, each with its common range
    for f in "${!LO[@]}"; do
      python3 $B/tools/plot_slices.py $src $OPTS --fields $f --clim ${LO[$f]} ${HI[$f]} | grep -v "colour range" | sed 's/^/   /'
    done
  else
    python3 $B/tools/plot_slices.py $src $OPTS ${FIELDS:+--fields $FIELDS} ${CLIM:+--clim $CLIM} | sed 's/^/   /'
  fi
done
echo "pictures of $(echo $files | wc -w) file(s) in $(( $(date +%s) - start )) s"
if [ ${#LO[@]} -gt 0 ]; then  # the common ranges, for the record
  o=${OUTDIR:-$(dirname $(echo $files | awk '{print $1}'))}/common_clim.txt
  { echo "# common colour ranges of the pictures of $(echo $files | wc -w) slice files: field, min, max"
    for f in "${!LO[@]}"; do printf "%-4s %15.7e %15.7e\n" $f ${LO[$f]} ${HI[$f]}; done; } > $o
  echo "written $o"
fi
