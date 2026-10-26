#!/usr/bin/env bash
# Install the two library dependencies of GIBD into deps/, at the commits used
# for the paper:
#   deps/sdsl-lite  vgteam fork of sdsl-lite (installed to deps/sdsl)
#   deps/gbwt       GBWT (Siren et al.), built against deps/sdsl-lite
# Requires git, a C++17 compiler, cmake and libzstd (with pkg-config).
# Usage: scripts/install_deps.sh [jobs]
set -euo pipefail

JOBS="${1:-1}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DEPS="$ROOT/deps"
SDSL_COMMIT=349de444ded81547cb55a718abeada41960531b5
GBWT_COMMIT=9e92e4f11bafcb1df92df9adf5991199c0f09f61

mkdir -p "$DEPS"
cd "$DEPS"

if [ ! -d sdsl-lite ]; then
    git clone https://github.com/vgteam/sdsl-lite.git
fi
(cd sdsl-lite && git checkout -q "$SDSL_COMMIT" && ./install.sh "$DEPS/sdsl")

if [ ! -d gbwt ]; then
    git clone https://github.com/jltsiren/gbwt.git
fi
(cd gbwt && git checkout -q "$GBWT_COMMIT" && make -j"$JOBS" SDSL_DIR=../sdsl-lite)

echo "Dependencies installed in $DEPS. Now run: make && make test"
