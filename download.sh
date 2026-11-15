#!/usr/bin/env bash
#
# Thin wrapper kept at the repository root for discoverability.
# The Phase 1 data collection stage lives in data-collection/.
#
set -euo pipefail

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec "$BASE_DIR/data-collection/download.sh" "$@"
