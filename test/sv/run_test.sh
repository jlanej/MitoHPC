#!/usr/bin/env bash
#########################################################################################
# Thin wrapper for the SV test harness (test/sv/run_test.py), which evaluates the caller
# against the committed mock BAMs (bams/ + truth.tsv), degenerate inputs, and a VCF-spec
# gate. Needs python3 + pysam (the caller is pure pysam). Point HP_PYTHON at an interpreter
# that has pysam if your default python3 does not.
#
#   bash test/sv/run_test.sh                                       # -> "ALL TESTS PASSED"
#   HP_PYTHON=/path/to/venv/bin/python bash test/sv/run_test.sh    # if pysam is in a venv
#########################################################################################
exec "${HP_PYTHON:-python3}" "$(cd "$(dirname "$0")" && pwd)/run_test.py" "$@"
