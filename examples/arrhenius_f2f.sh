#!/bin/bash
#
# Convert an FDCHAN restart file between single precision (fdchan_sp) and double precision (fdchan) on one
# Arrhenius node, optionally to another nx or nz, and print the lines for the FDCHAN input file.
#
# Build the converter once (login node):   cd $HOME && git clone https://github.com/GeertBrethouwer/converter.git
#                                          cd converter && make
# Set the variables below (or pass them on the command line), then submit from any directory:
#   sbatch -A <your-account> $HOME/converter/examples/arrhenius_f2f.sh
#   sbatch -A <your-account> --export=ALL,IN=/scratch/geert/fd_sp/data/restart.bin $HOME/converter/examples/arrhenius_f2f.sh
#
#SBATCH -J f2f
#SBATCH -p gpu
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=72
#SBATCH --gpus-per-node=1
#SBATCH -t 01:00:00

# ---- settings
IN=${IN:-/scratch/geert/fdchan_sp_run/data/restart.bin}  # FDCHAN restart file
PREC=${PREC:-double}  # precision of the output: double (fdchan) or single (fdchan_sp)
OUT=${OUT:-$(dirname $IN)/restart_$PREC.bin}             # FDCHAN restart file to write
NZ=${NZ:-}            # z-points of the output; empty: as the input (no change in z)
NX=${NX:-}            # nx of the output; empty: as the input
NSCAL=${NSCAL:-}      # keep only the first NSCAL scalars; empty: all
RS=${RS:-$HOME/converter/rsconv}

# ---- checks
for f in $RS $IN; do
  [ -e $f ] || { echo "not found: $f"; exit 1; }
done
if [ -e $OUT ]; then echo "$OUT exists: remove it or set OUT"; exit 1; fi
read nx ny nz prec < <(od -An -t d4 -N16 $IN)
size=$(stat -c %s $IN)
po=$([ $PREC = single ] && echo 4 || echo 8)
need=$(( (size - 40) / prec * po ))
[ -n "$NZ" ] && need=$(( need / nz * NZ ))
[ -n "$NX" ] && need=$(( need / nx * NX ))
free=$(df -B1 --output=avail $(dirname $OUT) | tail -1)
echo "input: nx = $nx, ny = $ny, nz = $nz, $([ $prec = 8 ] && echo double || echo single) precision; output: $PREC precision"
[ $prec = $po ] && [ -z "$NZ$NX$NSCAL" ] && echo "(the input is already in $PREC precision: the output is a copy)"
echo "output needs about $(( need / 1000000000 )) GB, free on $(dirname $OUT): $(( free / 1000000000 )) GB"
[ $free -gt $need ] || { echo "not enough disk space"; exit 1; }

# ---- convert (same nx and nz: the values are only cast, no transforms)
OPTS="--$PREC --threads $SLURM_CPUS_PER_TASK"
[ -n "$NZ" ] && OPTS="$OPTS --nz $NZ"
[ -n "$NX" ] && OPTS="$OPTS --nx $NX"
[ -n "$NSCAL" ] && OPTS="$OPTS --nscal $NSCAL"
echo "== $RS f2f $IN $OUT $OPTS"
start=$(date +%s)
srun $RS f2f $IN $OUT $OPTS || { echo "conversion FAILED"; exit 1; }
echo "conversion took $(( $(date +%s) - start )) s"
$RS info $OUT

# ---- lines for the FDCHAN input file
echo
echo "== for the FDCHAN input file ($([ $PREC = double ] && echo "run with the double-precision build, \$HOME/fdchan" || echo "run with the single-precision build, \$HOME/fdchan_sp")):"
[ -n "$NX" ] && echo "nx = $NX"
[ -n "$NZ" ] && echo "nz = $NZ"
[ -n "$NSCAL" ] && echo "nscal = $NSCAL"
echo "init_type = restart"
echo "restart_file = $OUT"
