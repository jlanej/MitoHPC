#!/usr/bin/env bash
set -uo pipefail
#########################################################################################
# Evaluate the MitoHPC SV caller (scripts/callSV.sh) against the mock BAMs in bams/ using
# the ground truth in truth.tsv. Real-time, self-contained: no full pipeline needed.
#
#   bash test/sv/run_test.sh
#
# Checks per sample:
#   * a deletion call matches the true breakpoints (within BP_TOL) and size (SVLEN_TOL)
#   * the heteroplasmy estimate AFJ matches truth (within AF_TOL)
#   * del4977 samples carry the REPEAT flag; the non-repeat deletion does not
#   * high-heteroplasmy deletions are FILTER=PASS; the 5% one is detected (split-only tier)
#   * the wild-type negative control yields zero PASS calls
#########################################################################################

HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
export HP_SDIR=$ROOT/scripts
export HP_RDIR=$ROOT/RefSeq
BAMS=$HERE/bams
TRUTH=$HERE/truth.tsv
OUT=$(mktemp -d)

BP_TOL=30; SVLEN_TOL=30; AF_TOL=0.12
fail=0

printf "%-16s %-8s %-22s %-8s %-8s %-12s %-7s %s\n" SAMPLE TRUTH_AF CALL_bp5/bp3 CALL_SVLEN CALL_AFJ FILTER FLAGS RESULT
printf -- "-------------------------------------------------------------------------------------------------------\n"

while IFS=$'\t' read -r s tbp5 tbp3 tsvlen thet _; do
  [ "${s#\#}" = "$s" ] || continue   # skip header
  bam=$BAMS/$s.bam
  [ -s "$bam" ] || { echo "MISSING BAM: $bam"; fail=1; continue; }
  bash "$HP_SDIR/callSV.sh" "$s" "$bam" "$OUT/$s" >/dev/null 2>&1
  tab=$OUT/$s.sv.tab

  if [ "$tbp5" = "." ]; then
    # negative control: no PASS calls allowed
    npass=$(awk -F'\t' '$0!~/^#/ && $13=="PASS"{n++} END{print n+0}' "$tab")
    if [ "$npass" -eq 0 ]; then res="PASS(neg)"; else res="FAIL(neg:$npass PASS calls)"; fail=1; fi
    printf "%-16s %-8s %-22s %-8s %-8s %-12s %-7s %s\n" "$s" "$thet" "-" "-" "-" "-" "-" "$res"
    continue
  fi

  # find the best-matching called row (breakpoints within tolerance)
  read -r cbp5 cbp3 csvlen cjr cafj cfilter cflags < <(awk -F'\t' -v t5="$tbp5" -v t3="$tbp3" -v tol="$BP_TOL" '
    $0~/^#/ {next}
    { d5=($3>t5)?$3-t5:t5-$3; d3=($4>t3)?$4-t3:t3-$4;
      if (d5<=tol && d3<=tol) { score=d5+d3; if (best=="" || score<bestv) {bestv=score; best=$3" "$4" "$5" "$6" "$8" "$13" "$14} } }
    END{ print (best=="" ? ". . . . . . ." : best) }' "$tab")

  ok=1; why=""
  if [ "$cbp5" = "." ]; then ok=0; why="not-detected"; fi
  if [ "$ok" = 1 ]; then
    awk -v a="$csvlen" -v b="$tsvlen" -v tol="$SVLEN_TOL" 'BEGIN{d=(a>b)?a-b:b-a; exit !(d<=tol)}' || { ok=0; why="svlen($csvlen vs $tsvlen)"; }
  fi
  if [ "$ok" = 1 ]; then
    awk -v a="$cafj" -v b="$thet" -v tol="$AF_TOL" 'BEGIN{d=(a>b)?a-b:b-a; exit !(d<=tol)}' || { ok=0; why="AFJ($cafj vs $thet)"; }
  fi
  # repeat-flag expectation
  case "$s" in
    *del4977*) echo "$cflags" | grep -q REPEAT || { ok=0; why="${why} missing-REPEAT"; } ;;
    *del6000*) echo "$cflags" | grep -q REPEAT && { ok=0; why="${why} unexpected-REPEAT"; } ;;
  esac
  # PASS expectation for >=10% heteroplasmy
  hi=$(awk -v h="$thet" 'BEGIN{print (h>=0.10)?1:0}')
  if [ "$hi" = 1 ] && [ "$cfilter" != "PASS" ]; then ok=0; why="${why} expected-PASS(got $cfilter)"; fi

  if [ "$ok" = 1 ]; then res="PASS"; else res="FAIL($why)"; fail=1; fi
  printf "%-16s %-8s %-22s %-8s %-8s %-12s %-7s %s\n" "$s" "$thet" "$cbp5/$cbp3" "$csvlen" "$cafj" "$cfilter" "$cflags" "$res"
done < "$TRUTH"

printf -- "-------------------------------------------------------------------------------------------------------\n"
rm -rf "$OUT"
if [ "$fail" = 0 ]; then echo "ALL TESTS PASSED"; exit 0; else echo "SOME TESTS FAILED"; exit 1; fi
