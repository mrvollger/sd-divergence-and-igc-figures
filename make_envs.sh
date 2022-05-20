#!/usr/bin/env bash
set -euo pipefail
conda env export > env.yml
conda env export --from-history > env.simple.yml

