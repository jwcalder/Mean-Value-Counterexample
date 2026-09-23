#!/bin/bash
#
#   ./reproduce.sh [nproc]
#
# Regenerates every certificate file in ../boxes/ from scratch.
#
# This is the search, not the proof.  Its output is a list of boxes; the proof
# is ../verify_all.sh, which re-proves every claim about those boxes and takes a
# small fraction of the time.  Run this only if you want to rediscover the boxes
# rather than check them.
#
# What is searched for, matching Tables 1 and 2 of the paper:
#
#   direct                          lifted
#     d = 3 : [2+e, 15.94257]         3 ->  4 : [24.25, 24.36]
#     d = 3 : [15.94261, 100]         3 ->  5 : [32.2, 32.5]
#     d = 4 : [2+e, 24.25]            4 ->  6 : [40.24, 40.36]
#     d = 4 : [24.36, 100]            5 ->  7 : [48.28, 48.4]
#       ...  and likewise for         4 ->  8 : [56.42, 56.54]
#     d = 5,6,7,8,9                   5 ->  9 : [64.67, 64.79]
#     d = 10: [2+e, 70]               4 -> 10 : [69.5, 76.5]
#     d = 10: [76, 100]
#
# with e = DMIN = 1e-7.  The gaps straddle the crossings p*_d at which C
# vanishes.  The lifted ranges cover them for d >= 4; the d = 3 gap cannot be
# closed, since lifting would need a planar counterexample and none exists.
#
# Cost: about 55000 boxes, of the order of 30 core-hours.
set -e

NPROC=${1:-$(nproc)}
LOG=reproduce.log
exec > >(tee "$LOG") 2>&1

DMIN=0.0000001        # smallest certified p - 2 (a measured limit, not a theorem)
BOXES_PER_CHUNK=40
MINW_DIV=64           # per-band minimum width = startwidth / MINW_DIV

[ -x ./certify_range ] || { echo "binary missing; run 'make' first"; exit 1; }

OUT=../boxes
mkdir -p "$OUT"
W=parts; rm -rf "$W"; mkdir -p "$W"

python3 - <<PYEOF > "$W/cmds.txt"
import math
DMIN=float("$DMIN"); per=$BOXES_PER_CHUNK; W="$W"; MINW_DIV=$MINW_DIV
FLAT=0.02                      # box width away from p = 2

cmds=[]
def chunks(d, a, b, dt, sw, group):
    n=max(1,math.ceil((b-a)/sw/per))
    for i in range(n):
        x1 = a if i==0 else a+(b-a)*i/n
        x2 = b if i==n-1 else a+(b-a)*(i+1)/n
        cmds.append((group, d, x1, x2, sw, dt))

def emit(d, a, b, dt, group, flat=FLAT):
    """Search [a,b].  Below p = 3 the certifiable box width shrinks with p - 2,
    because C(p) -> 0 and lambda*(p) -> 2 there, so that stretch is cut into
    geometric bands [2+t, 2+2t] with width proportional to t.  This keeps the
    box count logarithmic in the distance from 2."""
    if a < 3.0:
        edges=[]; t=a-2.0
        while 2.0+t < min(b,3.0):
            edges.append(2.0+t); t*=2.0
        edges.append(min(b,3.0))
        for k in range(len(edges)-1):
            lo,hi = edges[k], edges[k+1]; tt = lo-2.0
            sw = min(0.01*tt if tt < 1e-2 else 0.05*tt, flat)
            chunks(d, lo, hi, dt, sw, group)
        if b > 3.0:
            chunks(d, 3.0, b, dt, flat, group)
    else:
        chunks(d, a, b, dt, flat, group)

PSTAR3 = 15.94258208563744     # numerically observed crossing for d = 3
LO3    = 15.942570000000002    # 15.94257 rounded outward (up)
HI3    = 15.942609999999998    # 15.94261 rounded outward (down)

def close_on_crossing(d, outer, target, dt, group, below):
    """Bands closing on the d = 3 crossing, with the box width proportional to
    the distance from it.  This is the same device as the geometric bands near
    p = 2: the certifiable width shrinks with the distance because C does, so
    the box count grows only like the logarithm of the distance.

    Bands are emitted in increasing p either way, since the parts of a group are
    concatenated in the order produced and must tile the range."""
    tmin = abs(PSTAR3 - target)
    t = abs(PSTAR3 - outer)
    ts = []
    while t/2 > tmin:
        ts.append(t); t /= 2
    ts.append(t)                       # ts is decreasing: outermost first
    if below:
        seq = [(PSTAR3 - t, PSTAR3 - t/2) for t in ts]
        seq[0] = (outer, seq[0][1]); seq[-1] = (seq[-1][0], target)
    else:
        seq = [(PSTAR3 + t/2, PSTAR3 + t) for t in reversed(ts)]
        seq[0] = (target, seq[0][1]); seq[-1] = (seq[-1][0], outer)
    for (lo, hi), t in zip(seq, ts if below else list(reversed(ts))):
        chunks(d, lo, hi, dt, min(t/64, FLAT), group)

L = 2.0+DMIN
# Some endpoints below are written with their full seventeen digits rather than
# as the short decimal the paper prints.  The certificates record range
# endpoints rounded OUTWARD from the printed decimals -- downward at the left
# end, upward at the right -- because the printed range is a claim about every p
# in it and must therefore be contained in what was certified, and the double
# nearest a decimal such as 32.38 can fall on the wrong side.  Where the nearest
# double already falls on the correct side, as at 24.25, 69.5 or 70, the short
# form is used unchanged.
# direct ranges.  For d = 3 the two ranges close on the crossing from either
# side; for the other dimensions the gap is covered by a lift, so the plain flat
# cut is enough.
emit( 3, L,     15.9,   3, "d3_lower")
close_on_crossing( 3, 15.9,  LO3,  3, "d3_lower", below=True)
close_on_crossing( 3, 16.0,  HI3,  3, "d3_upper", below=False)
emit( 3, 16.0,  100.0,  3, "d3_upper")
emit( 4, L,     24.25,  4, "d4_lower")
emit( 4, 24.36, 100.0,  4, "d4_upper")
emit( 5, L,     32.28,  5, "d5_lower")
emit( 5, 32.379999999999995, 100.0,  5, "d5_upper")
emit( 6, L,     40.24,  6, "d6_lower")
emit( 6, 40.36, 100.0,  6, "d6_upper")
emit( 7, L,     48.28,  7, "d7_lower")
emit( 7, 48.4,  100.0,  7, "d7_upper")
emit( 8, L,     56.42,  8, "d8_lower")
emit( 8, 56.54, 100.0,  8, "d8_upper")
emit( 9, L,     64.67,  9, "d9_lower")
emit( 9, 64.78999999999999, 100.0,  9, "d9_upper")
emit(10, L,     70.0,  10, "d10_lower")
emit(10, 76.0,  100.0, 10, "d10_upper")
# lifted ranges; dtarget > d, and H is rational when dtarget - d is even,
# otherwise computed from a validated log-Gamma
emit( 3, 24.25, 24.360000000000003,  4, "lift_3to4", flat=0.002)
emit( 3, 32.199999999999996, 32.5,   5, "lift_3to5", flat=0.005)
emit( 4, 40.239999999999995, 40.36000000000001,  6, "lift_4to6", flat=0.002)
emit( 5, 48.279999999999994, 48.400000000000006,   7, "lift_5to7", flat=0.002)
emit( 4, 56.419999999999995, 56.540000000000006,  8, "lift_4to8", flat=0.002)
emit( 5, 64.66999999999999, 64.79,  9, "lift_5to9", flat=0.002)
emit( 4, 69.5,  76.5,  10, "lift_4to10")

count={}
for g,d,x1,x2,sw,dt in cmds:
    i=count.get(g,0); count[g]=i+1
    print("./certify_range %d %.17g %.17g %.17g %g %d > %s/%s.part%06d 2>/dev/null"
          % (d,x1,x2,sw/MINW_DIV,sw,dt,W,g,i))
PYEOF

echo "$(wc -l < "$W/cmds.txt") chunks on $NPROC processes"
date
time xargs -P "$NPROC" -I CMD -a "$W/cmds.txt" bash -c CMD || true
date

# Assemble the certificate files: one "claim" line stating what is to be
# proved, then the boxes.  The claim is what ../verify proves; the boxes are
# only the certificate for it.
echo
echo "assembling $OUT ..."
assemble () {   # name  claim-fields  header-comment
  { echo "# $3"; echo "claim $2"; echo "# p1 p2 lambda1 lambda2 tau0"
    cat $(ls -v "$W/$1.part"*); } > "$OUT/$1.txt"
  printf "%-16s %6d boxes\n" "$1.txt" "$(grep -vc '^#' "$OUT/$1.txt")"
}
L=2.0000001
assemble d3_lower   "3 3 $L 15.942570000000002 -"  "Direct counterexample,  d = 3:   C(p) < 0 for every p in [$L, 15.94257]"
assemble d3_upper   "3 3 15.942609999999998 100 +" "Direct counterexample,  d = 3:   C(p) > 0 for every p in [15.94261, 100]"
assemble d4_lower   "4 4 $L 24.25 -"    "Direct counterexample,  d = 4:  C(p) < 0 on [$L, 24.25]"
assemble d4_upper   "4 4 24.36 100 +"   "Direct counterexample,  d = 4:  C(p) > 0 on [24.36, 100]"
assemble d5_lower   "5 5 $L 32.28 -"    "Direct counterexample,  d = 5:  C(p) < 0 on [$L, 32.28]"
assemble d5_upper   "5 5 32.379999999999995 100 +"   "Direct counterexample,  d = 5:  C(p) > 0 on [32.38, 100]"
assemble d6_lower   "6 6 $L 40.24 -"    "Direct counterexample,  d = 6:  C(p) < 0 on [$L, 40.24]"
assemble d6_upper   "6 6 40.36 100 +"   "Direct counterexample,  d = 6:  C(p) > 0 on [40.36, 100]"
assemble d7_lower   "7 7 $L 48.28 -"    "Direct counterexample,  d = 7:  C(p) < 0 on [$L, 48.28]"
assemble d7_upper   "7 7 48.4 100 +"    "Direct counterexample,  d = 7:  C(p) > 0 on [48.4, 100]"
assemble d8_lower   "8 8 $L 56.42 -"    "Direct counterexample,  d = 8:  C(p) < 0 on [$L, 56.42]"
assemble d8_upper   "8 8 56.54 100 +"   "Direct counterexample,  d = 8:  C(p) > 0 on [56.54, 100]"
assemble d9_lower   "9 9 $L 64.67 -"    "Direct counterexample,  d = 9:  C(p) < 0 on [$L, 64.67]"
assemble d9_upper   "9 9 64.78999999999999 100 +"   "Direct counterexample,  d = 9:  C(p) > 0 on [64.79, 100]"
assemble d10_lower  "10 10 $L 70 -"     "Direct counterexample,  d = 10:  C(p) < 0 on [$L, 70]"
assemble d10_upper  "10 10 76 100 +"    "Direct counterexample,  d = 10:  C(p) > 0 on [76, 100]"
assemble lift_3to4  "3 4 24.25 24.360000000000003 +" "Lifted counterexample, 3 -> 4:  C_{3->4}(p) > 0 on [24.25, 24.36]"
assemble lift_3to5  "3 5 32.199999999999996 32.5 -"   "Lifted counterexample, 3 -> 5:  C_{3->5}(p) < 0 on [32.2, 32.5]"
assemble lift_4to6  "4 6 40.239999999999995 40.36000000000001 -" "Lifted counterexample, 4 -> 6:  C_{4->6}(p) < 0 on [40.24, 40.36]"
assemble lift_5to7  "5 7 48.279999999999994 48.400000000000006 -"  "Lifted counterexample, 5 -> 7:  C_{5->7}(p) < 0 on [48.28, 48.4]"
assemble lift_4to8  "4 8 56.419999999999995 56.540000000000006 -" "Lifted counterexample, 4 -> 8:  C_{4->8}(p) < 0 on [56.42, 56.54]"
assemble lift_5to9  "5 9 64.66999999999999 64.79 -" "Lifted counterexample, 5 -> 9:  C_{5->9}(p) < 0 on [64.67, 64.79]"
assemble lift_4to10 "4 10 69.5 76.5 -"  "Lifted counterexample, 4 -> 10: C_{4->10}(p) < 0 on [69.5, 76.5]"

if grep -q "^# FAIL" "$OUT"/*.txt; then
  echo
  echo "WARNING: the search failed on some subintervals:"
  grep -H "^# FAIL" "$OUT"/*.txt | head -20
  echo "Those ranges have no certificate, and ../verify_all.sh will reject the"
  echo "affected files at (C6).  Widen the ladder or lower the minimum width."
  exit 1
fi
rm -rf "$W"

echo
echo "boxes regenerated.  Now run ../verify_all.sh to prove the claims."
