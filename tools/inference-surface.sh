#!/bin/sh
# Snapshot the intentionally supported multivar-inference API independently
# from the ordinary analysis façade.

set -eu

root="$(CDPATH='' cd -- "$(dirname -- "$0")/.." && pwd)"

MULTIVAR_SURFACE_MODULE=inference \
MULTIVAR_SURFACE_TITLE="multivar-inference API" \
MULTIVAR_SURFACE_PROBE=multivar.inference.InferenceError \
MULTIVAR_SURFACE_CLASSES="$root/tools/public-surface/inference-classes.txt" \
MULTIVAR_SURFACE_SNAPSHOT="$root/tools/public-surface/inference-api.txt" \
  exec "$root/tools/public-surface.sh" "$@"
