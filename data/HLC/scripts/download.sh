#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
FILENAME="sdata_lung_s1.zarr.zip"

cd "$DATA_DIR"
curl -L -O "https://ftp.ebi.ac.uk/biostudies/fire/S-BIAD/146/S-BIAD2146/Files/STHELAR/sdata_slides/${FILENAME}"
unzip "$FILENAME"
