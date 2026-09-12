#!/bin/sh
# Install the pinned Gale revision into the local Ivy/Maven repository.
#
# Multivar depends on Gale by Maven coordinate so published POMs are resolvable.
# Until Gale is on Maven Central, every clean machine (and CI) must run this
# script before `sbt compileAll`. The installed version is Gale's own dynver
# for the pinned full-history revision, so source and artifact consumers use
# identical coordinates.
#
# Usage:
#   tools/publish-gale-local.sh
#
# Environment:
#   GALE_CACHE   optional checkout directory (default: ~/.cache/multivar-gale/<rev>)
#   SBT          optional sbt launcher (default: sbt)

set -eu

ROOT=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
BUILD_SBT="$ROOT/build.sbt"

gale_revision=$(
  sed -n 's/^lazy val galeRevision = "\([^"]*\)".*/\1/p' "$BUILD_SBT" | head -n 1
)

if [ -z "$gale_revision" ]; then
  echo "could not read galeRevision from $BUILD_SBT" >&2
  exit 1
fi

gale_version=$(
  sed -n 's/^lazy val galeVersion[[:space:]]*=[[:space:]]*"\([^"]*\)".*/\1/p' "$BUILD_SBT" | head -n 1
)

if [ -z "$gale_version" ]; then
  echo "could not read galeVersion from $BUILD_SBT" >&2
  exit 1
fi

SBT_BIN=${SBT:-sbt}
CACHE=${GALE_CACHE:-"$HOME/.cache/multivar-gale/$gale_revision"}

echo "publishing Gale $gale_revision as $gale_version"

if [ ! -d "$CACHE/.git" ]; then
  mkdir -p "$(dirname "$CACHE")"
  git clone --filter=blob:none https://github.com/canardlapin/gale.git "$CACHE"
fi

cd "$CACHE"
git fetch --depth 1 origin "$gale_revision" 2>/dev/null || git fetch origin
git checkout --force --detach "$gale_revision"

if [ "$(git rev-parse --is-shallow-repository)" = "true" ]; then
  git fetch --unshallow origin
  git checkout --force --detach "$gale_revision"
fi

actual_version=$(
  "$SBT_BIN" -Dsbt.supershell=false --batch "print coreJVM/version" | tail -n 1
)

if [ "$actual_version" != "$gale_version" ]; then
  echo "Gale version mismatch: build.sbt requires $gale_version but $gale_revision derives $actual_version" >&2
  exit 1
fi

"$SBT_BIN" -Dsbt.supershell=false --batch \
  "set ThisBuild / scalaVersion := \"3.7.4\"" \
  "coreJVM/publishLocal" \
  "coreJS/publishLocal"

echo "Gale $gale_version installed for local resolution"
