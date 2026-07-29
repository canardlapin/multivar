#!/bin/sh
set -eu

if rg -n --glob '*/src/main/**/*.scala' \
  '^(package|import)[[:space:]]+(scalafim\.|multivar\.numerics|breeze\.|org\.apache\.spark\.|org\.neuroimaging\.|neuroimaging\.)|scalafim\.(multivar|linalg)' \
  modules/core modules/ir modules/inference
then
  echo "production sources cross the multivar/resample4s/Gale ownership boundary" >&2
  exit 1
fi

echo "standalone package boundary: clean"
