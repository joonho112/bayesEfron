#!/usr/bin/env bash
# Install bayesEfron 0.4.0 (the public copy at tag v0.4.0, commit 89a0a71) into
# the project library r-library/, leaving the user's R library untouched.
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
src="${1:-$root/../bayesEfron-R-package-public}"
lib="$root/r-library"
mkdir -p "$lib"
cd "$src"
git rev-parse HEAD > "$root/qa/installed-commit.txt"
git describe --tags --exact-match HEAD >> "$root/qa/installed-commit.txt" 2>/dev/null || true
R CMD INSTALL --no-multiarch --library="$lib" . 2>&1 | tee "$root/qa/install.log"
R_LIBS_USER="$lib" Rscript -e 'cat("bayesEfron", as.character(packageVersion("bayesEfron", lib.loc = Sys.getenv("R_LIBS_USER"))), "\n")'
