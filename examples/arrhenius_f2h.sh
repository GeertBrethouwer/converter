#!/bin/bash
#
# Convert an FDCHAN restart file (single or double precision) to a HiCHAN restart file (<name>.bin + <name>.hdr) on
# one Arrhenius node, and print the lines for the HiCHAN input file.
#
# Build the converter once (login node):   cd $HOME && git clone https://github.com/GeertBrethouwer/converter.git
#                                          cd converter && make
# Set the variables below (or pass them on the command line), then submit from any directory:
#   sbatch -A <your-account> $HOME/converter/examples/arrhenius_f2h.sh
#   sbatch -A <your-account> --export=ALL,FDCHAN=/scratch/geert/fd_run/data/restart.bin,NZ=1920,DT=1.5e-3 \
#          $HOME/converter/examples/arrhenius_f2h.sh
#
#SBATCH -J f2h
#SBATCH -p gpu
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=72
#SBATCH --gpus-per-node=1
#SBATCH -t 01:00:00

# ---- settings
FDCHAN=${FDCHAN:-/scratch/geert/fdchan_run/data/restart.bin}  # FDCHAN restart file
OUT=${OUT:-$(dirname $FDCHAN)/restart_hichan}                  # HiCHAN restart to write: <OUT>.bin and <OUT>.hdr
NZ=${NZ:-}            # HiCHAN Fourier modes in z (even); empty: as FDCHAN's nz (z-points). FDCHAN with compact z
                      # needs about 1.4x the points of HiCHAN's modes for the same resolution (e.g. 2688 -> 1920)
NX=${NX:-}            # nx of the output; empty: as FDCHAN
NSCAL=${NSCAL:-}      # keep only the first NSCAL scalars; empty: all
STEP=${STEP:-0}       # step number for the HiCHAN header (FDCHAN files do not store it)
DT=${DT:-1e-3}        # time step for the HiCHAN header (HiCHAN's first time step; FDCHAN files do not store it)
GRID=${GRID:-"tanh 4.5 8"}  # grid_type grid_alpha fd_order for the HiCHAN header (information only)
RS=${RS:-$HOME/converter/rsconv}

# ---- checks
for f in $RS $FDCHAN; do
  [ -e $f ] || { echo "not found: $f"; exit 1; }
done
O=${OUT%.bin}; O=${O%.hdr}
if [ -e $O.bin ] || [ -e $O.hdr ]; then echo "$O.bin or $O.hdr exists: remove them or set OUT"; exit 1; fi
read nx ny nz prec < <(od -An -t d4 -N16 $FDCHAN)
size=$(stat -c %s $FDCHAN)
nf=$(( (size - 40) / (ny * nz * nx / 2 * 2 * prec) ))
nzo=${NZ:-$nz}
if [ $(( nzo % 2 )) -ne 0 ]; then echo "HiCHAN needs an even nz (FDCHAN has nz = $nz): set NZ"; exit 1; fi
need=$(( (size - 40) / prec * 8 / nz * nzo ))
[ -n "$NX" ] && need=$(( need / nx * NX ))
free=$(df -B1 --output=avail $(dirname $O) | tail -1)
echo "FDCHAN file: nx = $nx, ny = $ny, nz = $nz, $(( nf - 3 )) scalar(s), $([ $prec = 8 ] && echo double || echo single) precision"
echo "output needs about $(( need / 1000000000 )) GB, free on $(dirname $O): $(( free / 1000000000 )) GB"
[ $free -gt $need ] || { echo "not enough disk space"; exit 1; }

# ---- convert
OPTS="--step $STEP --dt $DT --grid $GRID --threads $SLURM_CPUS_PER_TASK"
[ -n "$NZ" ] && OPTS="$OPTS --nz $NZ"
[ -n "$NX" ] && OPTS="$OPTS --nx $NX"
[ -n "$NSCAL" ] && OPTS="$OPTS --nscal $NSCAL"
echo "== $RS f2h $FDCHAN $O $OPTS"
start=$(date +%s)
srun $RS f2h $FDCHAN $O $OPTS || { echo "conversion FAILED"; exit 1; }
echo "conversion took $(( $(date +%s) - start )) s"
echo "(HiCHAN has no kz Nyquist mode: the 'Nyquist kz' fractions above are dropped; they should be small)"
$RS info $O

# ---- lines for the HiCHAN input file
read hnx hny hnz hns < <(awk 'NR==2' $O.hdr)
read hlx hlz < <(awk 'NR==3' $O.hdr)
read gt ga fo <<< "$GRID"
echo
echo "== for the HiCHAN input file (plus the physics, scalar and output settings of the FDCHAN run; grid_type and"
echo "   grid_alpha must be those of the FDCHAN run, the restart files do not contain the y-grid):"
echo "&grid"
echo "  nx = $hnx, nz = $hnz, ny = $hny, lx = $hlx, lz = $hlz"
echo "  grid_type = '$gt', grid_alpha = $ga, fd_order = $fo"
echo "/"
echo "&init"
echo "  init_type = 'restart', restart_file = '$O'"
echo "/"
echo "&scalars"
echo "  nscal = $hns"
echo "/"
echo "(the first div_max in HiCHAN's log.dat after this restart is larger; from the first time step on it is at round-off)"
