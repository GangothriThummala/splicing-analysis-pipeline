#!/bin/bash
# Install the PAIRADISE that rMATS --paired-stats needs into the active conda env.
#
# rMATS 4.4.0 (rMATS_R/paired_model.R) passes a data.frame to pairadise(): the
# API of Xinglab/PAIRADISE (GitHub, v1.0), which rMATS's own build_rmats installs.
# Bioconductor's PAIRADISE (1.x on bioconda) is a rewrite that only accepts a
# PDseDataSet; with it the paired model fails, while rMATS still exits 0 and
# writes NA p-values. So it is not in splice.yaml, and this script pins the
# GitHub release by tag and commit.
#
# Usage (after creating the env): conda activate splice && envs/install_pairadise.sh
set -euo pipefail
TAG=v1.0.0
COMMIT=eff57cb3cc3fa23675d07d6c80b63dc6d864278e

[ -n "${CONDA_PREFIX:-}" ] || { echo "Activate the splice env first" >&2; exit 1; }
if Rscript -e 'quit(status = !("PAIRADISE" %in% rownames(installed.packages())))' 2>/dev/null \
   && ! Rscript -e 'quit(status = packageVersion("PAIRADISE") != "1.0")' 2>/dev/null; then
  echo "A different PAIRADISE is installed; remove it first: conda remove --force bioconductor-pairadise" >&2
  exit 1
fi

src=$(mktemp -d)
trap 'rm -rf "$src"' EXIT
git clone -q --branch "$TAG" --depth 1 https://github.com/Xinglab/PAIRADISE.git "$src/PAIRADISE"
got=$(git -C "$src/PAIRADISE" rev-parse HEAD)
[ "$got" = "$COMMIT" ] || { echo "Tag $TAG is $got, expected $COMMIT" >&2; exit 1; }

R CMD INSTALL "$src/PAIRADISE/pairadise/src/pairadise_model"
Rscript -e 'library(PAIRADISE); stopifnot(packageVersion("PAIRADISE") == "1.0", "my.data" %in% names(formals(pairadise)))'
echo "PAIRADISE $TAG ($COMMIT) installed in $CONDA_PREFIX"
