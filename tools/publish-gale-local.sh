#!/bin/sh
# Install the pinned Gale revision into the local Ivy/Maven repository.
#
# Multivar depends on Gale by Maven coordinate so published POMs are resolvable.
# Until Gale is on Maven Central, every clean machine (and CI) must run this
# script before `sbt compileAll`. The installed version embeds the revision
# prefix so a different local SNAPSHOT cannot silently substitute.
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

# Mirror build.sbt: 1.0.0-<first 12 hex characters of the revision>.
gale_version="1.0.0-$(printf '%s' "$gale_revision" | cut -c1-12)"
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

"$SBT_BIN" -Dsbt.supershell=false --batch \
  "set ThisBuild / version := \"$gale_version\"" \
  "set ThisBuild / scalaVersion := \"3.7.4\"" \
  "coreJVM/publishLocal" \
  "coreJS/publishLocal"

echo "Gale $gale_version installed for local resolution"
