#!/bin/sh
# Install the pinned Resample4s revision into the local Ivy/Maven repository.
#
# multivar-inference depends on Resample4s by Maven coordinate so its published
# POM remains resolvable. Until Resample4s is on Maven Central, every clean
# machine (and CI) must run this script before `sbt compileAll`. The installed
# version embeds the revision prefix so a different local SNAPSHOT cannot
# silently substitute.
#
# Usage:
#   tools/publish-resample4s-local.sh
#
# Environment:
#   RESAMPLE4S_CACHE  optional checkout directory
#                      (default: ~/.cache/multivar-resample4s/<rev>)
#   SBT               optional sbt launcher (default: sbt)

set -eu

ROOT=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
BUILD_SBT="$ROOT/build.sbt"

resample4s_revision=$(
  sed -n 's/^lazy val resample4sRevision = "\([^"]*\)".*/\1/p' "$BUILD_SBT" |
    head -n 1
)

if [ -z "$resample4s_revision" ]; then
  echo "could not read resample4sRevision from $BUILD_SBT" >&2
  exit 1
fi

# Mirror build.sbt: 0.1.0-<first 12 hex characters of the revision>.
resample4s_version="0.1.0-$(printf '%s' "$resample4s_revision" | cut -c1-12)"
SBT_BIN=${SBT:-sbt}
CACHE=${RESAMPLE4S_CACHE:-"$HOME/.cache/multivar-resample4s/$resample4s_revision"}

echo "publishing Resample4s $resample4s_revision as $resample4s_version"

if [ ! -d "$CACHE/.git" ]; then
  mkdir -p "$(dirname "$CACHE")"
  git clone --filter=blob:none \
    https://github.com/canardlapin/resample4s.git "$CACHE"
fi

cd "$CACHE"
git fetch --depth 1 origin "$resample4s_revision" 2>/dev/null || git fetch origin
git checkout --force --detach "$resample4s_revision"

"$SBT_BIN" -Dsbt.supershell=false --batch \
  "set ThisBuild / version := \"$resample4s_version\"" \
  "set ThisBuild / scalaVersion := \"3.7.4\"" \
  "coreJVM/publishLocal" \
  "coreJS/publishLocal" \
  "designsJVM/publishLocal" \
  "designsJS/publishLocal" \
  "apiJVM/publishLocal" \
  "apiJS/publishLocal"

echo "Resample4s $resample4s_version installed for local resolution"
