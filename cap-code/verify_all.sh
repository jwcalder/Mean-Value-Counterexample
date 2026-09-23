#!/bin/bash
#
#   ./verify_all.sh [nproc]
#
# Verifies every certificate file in boxes/ and prints the resulting theorems.
# This is the whole computer-assisted proof.  Exit status is 0 if and only if
# every box of every file passed every check.
#
# The boxes are independent, so the work is split across processes: each
# certificate file is cut into chunks of CHUNK boxes and each chunk is verified
# separately.  Splitting only divides the per-box work -- every chunk re-checks
# (C6) over the whole file, that being a property of the file rather than of any
# box.  After the run the verified boxes are compared against the certificate
# line by line, so that a mistake in the chunking here would be caught rather
# than silently leaving boxes unverified.
set -e

NPROC=${1:-$(nproc)}
CHUNK=64
LOG=verify_all.log
exec > >(tee "$LOG") 2>&1

[ -x ./verify ] || { echo "run 'make' first"; exit 1; }

# Scratch goes outside the source tree on purpose: hundreds of chunk files are
# written concurrently, and a synced folder (Dropbox and the like) resolves the
# races by creating "conflicted copy" files, which the assembly below must never
# splice into a result.  Override with SCRATCH=... if /tmp is small.
SCRATCH=${SCRATCH:-/tmp/capd_chunks_$USER}
W="$SCRATCH"; rm -rf "$W"; mkdir -p "$W"; rm -f certified_*.txt
echo "scratch (not synced): $W"

nboxes () { grep -vc '^#\|^claim' "$1"; }

: > "$W/cmds.txt"
# The dimension-tail corollary needs one extra claim, (C7), for two of the
# certificates: that the spherical mean of the source profile is strictly
# negative.  The --lifting-tail flag checks it in addition to (C1)-(C6), at no
# extra cost, since the mean is already formed inside (C5).
for f in boxes/*.txt; do
  n=$(basename "$f" .txt)
  case "$n" in
    d10_lower|lift_4to10) flag=--lifting-tail ;;
    *)                    flag= ;;
  esac
  N=$(nboxes "$f")
  i=1; k=0
  while [ "$i" -le "$N" ]; do
    j=$((i + CHUNK - 1)); [ "$j" -gt "$N" ] && j=$N
    printf './verify %s %s %d %d > %s/%s.part%06d\n' "$flag" "$f" "$i" "$j" "$W" "$n" "$k" \
      >> "$W/cmds.txt"
    i=$((j + 1)); k=$((k + 1))
  done
done

echo "verifying $(cat boxes/*.txt | grep -vc '^#\|^claim') boxes"
echo "in $(wc -l < "$W/cmds.txt") chunks on $NPROC processes"
date
status=0
time xargs -P "$NPROC" -I CMD -a "$W/cmds.txt" bash -c CMD || status=1
date

if [ "$status" != 0 ]; then
  echo; echo "VERIFICATION FAILED -- see the messages above."; exit 1
fi

# Assemble, and check that what was verified is exactly what the certificate
# says, box for box and in order.  This catches any error in the chunking
# above; it is not a mathematical check, the mathematics being (C1)-(C6).
echo
echo "assembling ..."
for f in boxes/*.txt; do
  n=$(basename "$f" .txt)
  parts=$(ls -v "$W"/$n.part[0-9][0-9][0-9][0-9][0-9][0-9] 2>/dev/null)
  [ -n "$parts" ] || { echo "no chunks for $n"; exit 1; }
  cat $parts > "certified_$n.txt"
  # Compare the exponent boxes as numbers, not as text.  A certificate may
  # write an endpoint as the shortest decimal that round-trips (40.36000000000001)
  # while the verifier prints seventeen significant digits (40.360000000000007);
  # those are the same double.  Reformatting both through %.17g canonicalizes
  # them, and seventeen significant digits determine a double uniquely, so the
  # comparison still distinguishes genuinely different boxes.
  if ! diff -q <(grep '^OK' "certified_$n.txt" | awk '{printf "%.17g %.17g\n", $2, $3}') \
                <(grep -v '^#\|^claim' "$f" | awk '{printf "%.17g %.17g\n", $1, $2}') > /dev/null; then
    echo "MISMATCH: the boxes verified for $n are not the boxes in $f"; exit 1
  fi
  printf "%-18s %6d boxes verified" "$n" "$(grep -c '^OK' "certified_$n.txt")"
  case "$n" in
    d10_lower|lift_4to10)
      awk '$1=="OK"{if(!seen++ || $11>m) m=$11}
           END{if(seen) printf "   (C7): mean <= %.17g", m
               else     printf "   (C7): NO BOXES"}' "certified_$n.txt" ;;
  esac
  echo
done
rm -rf "$W"

# Every box passed (C1)-(C5), every chunk checked (C6) over its whole file, and
# the boxes verified are exactly the boxes of the certificate.  The claims
# stated in the certificate headers are therefore proved.
echo
echo "================ verified statements ================"
awk '/^claim/ {
       d = $2; dt = $3
       sym = (d == dt) ? "C(p)" : "C_{" d "->" dt "}(p)"
       printf "%-18s for every p in [%s, %s]: %s %s 0   (d = %s)\n",
              FILENAME, $4, $5, sym, ($6 == "+" ? ">" : "<"), d
     }' boxes/*.txt

echo
echo "ALL CERTIFICATES VERIFIED.  Log written to $LOG"

# The tables are a convenience, not part of the proof, so a failure here must
# not be read as a failure of the verification reported above.
echo
echo "================ tables ================"
echo "%% Table 1 (direct counterexamples)"
python3 gen_table.py $(for d in 3 4 5 6 7 8 9 10; do \
  echo certified_d${d}_lower.txt certified_d${d}_upper.txt; done) \
  || echo "(table generation failed; the verification above is unaffected)"
echo "%% Table 2 (lifted counterexamples)"
python3 gen_table.py certified_lift_3to4.txt certified_lift_3to5.txt \
  certified_lift_4to6.txt certified_lift_5to7.txt certified_lift_4to8.txt \
  certified_lift_5to9.txt certified_lift_4to10.txt \
  || echo "(table generation failed; the verification above is unaffected)"
