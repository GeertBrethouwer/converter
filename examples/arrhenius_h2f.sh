#!/bin/bash
#
# Convert a HiCHAN restart file to an FDCHAN restart file on one Arrhenius node, and print the lines for the
# FDCHAN input file.
#
# Build the converter once (login node):   cd $HOME && git clone https://github.com/GeertBrethouwer/converter_and_spectra.git
#                                          cd converter_and_spectra && make
# Set the variables below (or pass them on the command line), then submit from any directory:
#   sbatch -A <your-account> $HOME/converter_and_spectra/examples/arrhenius_h2f.sh
#   sbatch -A <your-account> --export=ALL,HICHAN=/scratch/geert/run/data/restart,NZ=2688 $HOME/converter_and_spectra/examples/arrhenius_h2f.sh
#
#SBATCH -J h2f
#SBATCH -p gpu
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=72
#SBATCH --gpus-per-node=1
#SBATCH -t 01:00:00

# ---- settings
HICHAN=${HICHAN:-/scratch/geert/hichan_run/data/restart}   # HiCHAN restart: <name> of <name>.bin + <name>.hdr
OUT=${OUT:-$(dirname $HICHAN)/restart_fdchan.bin}           # FDCHAN restart file to write
NZ=${NZ:-}            # FDCHAN z-points; empty: as HiCHAN's nz (Fourier modes). About 1.4x for the same resolution
                      # with compact 8th-order z (e.g. 1920 -> 2688); nz/#GPUs >= 49 for z_scheme = compact, order 8
NX=${NX:-}            # nx of the output; empty: as HiCHAN
PREC=${PREC:-double}  # double (fdchan) or single (fdchan_sp)
NSCAL=${NSCAL:-}      # keep only the first NSCAL scalars; empty: all
RS=${RS:-$HOME/converter_and_spectra/rsconv}

# ---- checks
H=${HICHAN%.bin}; H=${H%.hdr}
# binaries for the processor of this node (rebuilt by make if they were built on another type, e.g. a login node)
make -s -C $(dirname $RS) rsconv > /dev/null || { echo "building rsconv in $(dirname $RS) FAILED"; exit 1; }
for f in $RS $H.bin $H.hdr; do
  [ -e $f ] || { echo "not found: $f"; exit 1; }
done
if [ -e $OUT ]; then echo "$OUT exists: remove it or set OUT"; exit 1; fi
need=$(stat -c %s $H.bin)
[ -n "$NZ" ] && need=$(( need / $(awk 'NR==2{print $3}' $H.hdr) * NZ ))
[ "$PREC" = single ] && need=$(( need / 2 ))
free=$(df -B1 --output=avail $(dirname $OUT) | tail -1)
echo "output needs about $(( need / 1000000000 )) GB, free on $(dirname $OUT): $(( free / 1000000000 )) GB"
[ $free -gt $need ] || { echo "not enough disk space"; exit 1; }

# ---- convert
OPTS="--$PREC --threads $SLURM_CPUS_PER_TASK"
[ -n "$NZ" ] && OPTS="$OPTS --nz $NZ"
[ -n "$NX" ] && OPTS="$OPTS --nx $NX"
[ -n "$NSCAL" ] && OPTS="$OPTS --nscal $NSCAL"
echo "== $RS h2f $H $OUT $OPTS"
start=$(date +%s)
srun $RS h2f $H $OUT $OPTS || { echo "conversion FAILED"; exit 1; }
echo "conversion took $(( $(date +%s) - start )) s"
$RS info $OUT

# ---- lines for the FDCHAN input file (from the HiCHAN header; check grid_type and grid_alpha against the HiCHAN
# input: they must be the same, the restart files do not contain the y-grid)
read nx ny nz ns < <(awk 'NR==2' $H.hdr)
read lx lz < <(awk 'NR==3' $H.hdr)
read gt ga fo < <(awk 'NR==5' $H.hdr)
[ -n "$NX" ] && nx=$NX
[ -n "$NZ" ] && nz=$NZ
[ -n "$NSCAL" ] && ns=$NSCAL
echo
echo "== for the FDCHAN input file (plus the physics, scalar and output settings of the HiCHAN run):"
echo "nx = $nx"
echo "ny = $ny"
echo "nz = $nz"
echo "lx = $lx"
echo "lz = $lz"
[ -n "$gt" ] && echo "grid_type = $gt" && echo "grid_alpha = $ga"
echo "nscal = $ns"
echo "init_type = restart"
echo "restart_file = $OUT"
