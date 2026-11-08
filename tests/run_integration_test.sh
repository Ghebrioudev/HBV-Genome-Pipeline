#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TEST_DIR="$(mktemp -d)"
trap 'rm -rf "$TEST_DIR"' EXIT

cp -R "$BASE_DIR/." "$TEST_DIR/"

if ! PYTHON=true "$TEST_DIR/run.sh" --clean --skip-download >/dev/null; then
    printf 'Integration test failed: clean pipeline run did not complete.\n' >&2
    exit 1
fi

if "$TEST_DIR/run.sh" --not-a-real-option >/dev/null 2>&1; then
    printf 'Integration test failed: unknown options were accepted.\n' >&2
    exit 1
fi

printf 'Integration test passed.\n'
