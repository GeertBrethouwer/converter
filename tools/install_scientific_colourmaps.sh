#!/bin/bash
#
# Install the Scientific Colour Maps of Fabio Crameri (https://www.fabiocrameri.ch/colourmaps,
# doi:10.5281/zenodo.1243862) into a directory, for tools/plot_slices.py and other programs:
#
#   bash tools/install_scientific_colourmaps.sh [DIR]          (default DIR: /scratch/geert/scientific_colourmaps)
#
# 1. the Python package cmcrameri from PyPI, without its dependencies (numpy, matplotlib are not touched):
#    DIR/cmcrameri/cmaps/<name>.txt (all maps as 256 RGB lines) and, with PYTHONPATH=DIR, "import cmcrameri.cm";
# 2. the official collection from Zenodo (all maps for MATLAB, Python, GMT, ParaView, ... in DIR/ScientificColourMaps*),
#    if Zenodo can be reached (else only 1).
# Run it on a login node (internet access).
DIR=${1:-/scratch/geert/scientific_colourmaps}
mkdir -p $DIR || exit 1
cd $DIR || exit 1
ok=0

echo "== 1. Python package cmcrameri (PyPI) -> $DIR/cmcrameri"
if python3 -m pip install --no-deps --upgrade --target $DIR cmcrameri > pip.log 2>&1; then
  echo "   $(ls $DIR/cmcrameri/cmaps/*.txt | wc -l) colour maps in $DIR/cmcrameri/cmaps"
  ok=1
else
  echo "   pip failed (see $DIR/pip.log)"
fi

echo "== 2. official collection (Zenodo, latest version)"
url=$(python3 - <<'EOF' 2>/dev/null
import json, urllib.request
# the newest version of the record (concept record 1243862)
q = "https://zenodo.org/api/records?q=conceptrecid:1243862&sort=mostrecent&size=1"
r = json.load(urllib.request.urlopen(q, timeout=60))["hits"]["hits"][0]
for f in r["files"]:
    if f["key"].lower().endswith(".zip"):
        print(f["links"]["self"]); break
EOF
)
if [ -n "$url" ] && curl -fsSL --max-time 600 -o scm.zip "$url"; then
  unzip -q -o scm.zip && rm -f scm.zip
  d=$(ls -d ScientificColourMaps* 2>/dev/null | head -1)
  echo "   $d: $(ls -d $d/*/ 2>/dev/null | wc -l) colour maps (MATLAB: addpath(genpath('$DIR/$d')); see its README)"
  ok=1
else
  echo "   Zenodo not reachable from here (or no zip found): only the maps of 1. (enough for plot_slices.py)"
  echo "   (by hand: download the zip from https://doi.org/10.5281/zenodo.1243862 and unzip it in $DIR)"
fi

[ $ok = 1 ] || { echo "nothing installed"; exit 1; }
echo
echo "tools/plot_slices.py finds the maps in $DIR (or give --cmap-dir $DIR, or export SCIENTIFIC_COLOURMAPS=$DIR)."
echo "Python:  PYTHONPATH=$DIR python3 -c 'import cmcrameri.cm as cmc; print(cmc.vik)'"
