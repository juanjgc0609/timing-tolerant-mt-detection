#!/bin/bash
# Fetch the three public datasets of data/README.md at their registered commits into $TTMT_DATA_ROOT and verify every
# file against data/checksums. Files are written from the git blobs (git cat-file), so they are byte-exact whatever
# core.autocrlf is, and the PyMPA37 file names that Windows cannot check out are never checked out.
#   TTMT_DATA_ROOT=/path/to/realdata bash tools/fetch_data.sh
set -euo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
ROOT=${TTMT_DATA_ROOT:?set TTMT_DATA_ROOT}
SRC=$ROOT/_src
mkdir -p "$SRC"

clone() {   # url commit dir
  [ -d "$SRC/$3/.git" ] || git clone --quiet --filter=blob:none --no-checkout "$1" "$SRC/$3"
  git -C "$SRC/$3" cat-file -e "$2^{commit}" 2>/dev/null || git -C "$SRC/$3" fetch --quiet origin "$2"
  [ "$(git -C "$SRC/$3" rev-parse "$2^{commit}")" = "$2" ]
}
blob() {    # dir commit path dest
  mkdir -p "$(dirname "$4")"
  git -C "$SRC/$1" cat-file blob "$2:$3" > "$4"
}
check() {   # folder list
  (cd "$1" && sha256sum -c --quiet "$REPO/data/checksums/$2") && echo "$2: $(wc -l < "$REPO/data/checksums/$2") files OK"
}

# Ridgecrest 2019 (Task R0)
C=9a1f4dd5f7659f75707435da3ce854b6e48321d1; T=tasks/earthquake-phase-association
clone https://github.com/benchflow-ai/skillsbench.git $C skillsbench
for f in environment/data/wave.mseed environment/data/stations.csv verifier/catalog.csv; do
  blob skillsbench $C "$T/$f" "$ROOT/skillsbench_eqassoc/$(basename $f)"; done
check "$ROOT/skillsbench_eqassoc" ridgecrest_2019_R0.sha256

# Alto Tiberina (C16 and Task 2M)
C=4a14b134afec3c1f527b738ca3bae4cdc7cbbbf6; D=$ROOT/pympa37_altotiberina_100723
clone https://github.com/avuan/PyMPA37.git $C PyMPA37
for f in $(git -C "$SRC/PyMPA37" ls-tree --name-only $C main.pympa.dir/24h/ | grep '/100723\.'); do blob PyMPA37 $C "$f" "$D/$(basename $f)"; done
for f in main.pympa.dir/200.100723.cat main.pympa.dir/templates.zmap LICENSE input.calculate_ttimes.dir/inv.IV \
         input.calculate_ttimes.dir/aquila_kato.tvel input.download_data.dir/stations.list; do
  blob PyMPA37 $C "$f" "$D/$(basename $f)"; done
for f in $(awk '{print $2}' "$REPO/data/checksums/altotiberina_other_days_task2M.sha256"); do
  blob PyMPA37 $C "main.pympa.dir/24h/$f" "$ROOT/pympa37_days/$f"; done
cp "$REPO/data/checksums/altotiberina_other_days_task2M.sha256" "$ROOT/pympa37_days.sha256"
check "$D" altotiberina_20100723_C16.sha256
check "$ROOT/pympa37_days" altotiberina_other_days_task2M.sha256

# Hi-net swarm 2012 (Task R1)
C=bac5ccb00816a4ec0aaab215ea1c342c730d412f
clone https://github.com/Geolab-USTC/MatchLocate2.git $C MatchLocate2
git -C "$SRC/MatchLocate2" ls-tree -r --name-only $C example/ | while read -r f; do
  blob MatchLocate2 $C "$f" "$ROOT/MatchLocate2/$f"; done
check "$ROOT/MatchLocate2/example" hinet_matchlocate2_example.sha256
