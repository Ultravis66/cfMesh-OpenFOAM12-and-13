#!/bin/bash
# cfMesh OpenFOAM 13 - Rotor 37 Regression Test
# Runs cartesianMesh and validates output against known-good baselines.
# Exit 0 = pass, Exit 1 = fail

set -uo pipefail

CASE_DIR="${1:-$HOME/rotor37_mesh}"
LOG="$CASE_DIR/regression.log"

# Baselines
BASELINE_CELLS=2139231
CELL_TOL=0.05

NONORTHO_MAX=15.0

BAD_FACES_MAX=10

KNOWN_BAD_CELLS=13
KNOWN_NEG_VOL=7

PASS=0
FAIL=0

log()  { echo "$@" | tee -a "$LOG"; }
pass() { log "  PASS: $*"; ((PASS++)); }
fail() { log "  FAIL: $*"; ((FAIL++)); }

log ""
log "======================================================"
log " cfMesh Rotor 37 Regression Test"
log " Case: $CASE_DIR"
log " $(date)"
log "======================================================"

# Step 1: Run cartesianMesh
log ""
log "[ 1/3 ] Running cartesianMesh..."
cd "$CASE_DIR"

START=$(date +%s)
cartesianMesh > "$LOG.mesh" 2>&1
EXIT=$?
END=$(date +%s)
RUNTIME=$((END - START))

if [ $EXIT -ne 0 ]; then
    fail "cartesianMesh exited with code $EXIT"
    tail -20 "$LOG.mesh" | tee -a "$LOG"
    exit 1
fi
pass "cartesianMesh completed in ${RUNTIME}s"

BAD=$(grep "Number of bad faces is" "$LOG.mesh" | tail -1 | awk '{print $NF}')
BAD=${BAD:-999}
if [ "$BAD" -le "$BAD_FACES_MAX" ]; then
    pass "Bad faces after untangling: $BAD (limit: $BAD_FACES_MAX)"
else
    fail "Bad faces after untangling: $BAD (limit: $BAD_FACES_MAX)"
fi

# Step 2: Run checkMesh
log ""
log "[ 2/3 ] Running checkMesh..."
CHECKMESH=$(checkMesh -case "$CASE_DIR" 2>&1)

CELLS=$(echo "$CHECKMESH" | grep "    cells:" | awk '{print $2}')
CELL_LOW=$(echo "$BASELINE_CELLS $CELL_TOL" | awk '{printf "%d", $1*(1-$2)}')
CELL_HIGH=$(echo "$BASELINE_CELLS $CELL_TOL" | awk '{printf "%d", $1*(1+$2)}')
if [ -n "$CELLS" ] && [ "$CELLS" -ge "$CELL_LOW" ] && [ "$CELLS" -le "$CELL_HIGH" ]; then
    pass "Cell count: $CELLS (baseline: $BASELINE_CELLS +/-5%)"
else
    fail "Cell count: $CELLS (expected ${CELL_LOW}-${CELL_HIGH})"
fi

NONORTHO=$(echo "$CHECKMESH" | grep "non-orthogonality" | grep "average" | awk '{print $NF}')
if [ -n "$NONORTHO" ]; then
    OK=$(echo "$NONORTHO $NONORTHO_MAX" | awk '{print ($1 <= $2) ? "yes" : "no"}')
    if [ "$OK" = "yes" ]; then
        pass "Average non-orthogonality: ${NONORTHO}deg (limit: ${NONORTHO_MAX}deg)"
    else
        fail "Average non-orthogonality: ${NONORTHO}deg (limit: ${NONORTHO_MAX}deg)"
    fi
else
    fail "Could not parse average non-orthogonality"
fi

BAD_ASPECT=$(echo "$CHECKMESH" | grep "High aspect ratio" | grep -o "number of cells [0-9]*" | awk '{print $NF}')
BAD_ASPECT=${BAD_ASPECT:-0}
if [ "$BAD_ASPECT" -le "$KNOWN_BAD_CELLS" ]; then
    pass "High aspect ratio cells: $BAD_ASPECT (known limit: $KNOWN_BAD_CELLS)"
else
    fail "High aspect ratio cells: $BAD_ASPECT (regression from known $KNOWN_BAD_CELLS)"
fi

NEG_VOL=$(echo "$CHECKMESH" | grep "negative volume" | grep -o "Number of negative volume cells: [0-9]*" | awk '{print $NF}')
NEG_VOL=${NEG_VOL:-0}
if [ "$NEG_VOL" -le "$KNOWN_NEG_VOL" ]; then
    pass "Negative volume cells: $NEG_VOL (known limit: $KNOWN_NEG_VOL)"
else
    fail "Negative volume cells: $NEG_VOL (regression from known $KNOWN_NEG_VOL)"
fi

# Step 3: Summary
log ""
log "[ 3/3 ] Summary"
log "======================================================"
log " PASSED: $PASS"
log " FAILED: $FAIL"
log "======================================================"
log ""

if [ "$FAIL" -eq 0 ]; then
    log " RESULT: PASS"
    exit 0
else
    log " RESULT: FAIL"
    exit 1
fi
