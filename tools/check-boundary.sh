#!/bin/sh
set -eu

if rg -n --glob '*/src/main/**/*.scala' \
  '^(package|import)[[:space:]]+(scalafim\.|multivar\.numerics)|scalafim\.(multivar|linalg)' \
  modules/core modules/ir
then
  echo "production sources must not depend on scalafim packages or recreate multivar.numerics" >&2
  exit 1
fi

echo "standalone package boundary: clean"
