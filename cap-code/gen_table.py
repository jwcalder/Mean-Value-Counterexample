"""Generate the LaTeX results table from the certified range files.

    python3 gen_table.py range_d3_3_12.5.txt range_d3_19_50.txt ...

For each range file it prints one table row:

    d & [a,b] & N & union of lambda enclosures & worst-case C bound & sign

Enclosure columns are rounded OUTWARD from the raw certified values, so the
printed bounds are implied by the computation.  Range endpoints are rounded
INWARD, for the same reason in the opposite direction: the printed range is
claimed for every p in it, so it must lie inside what was certified.  The file name encodes d and
the range; the certified boxes inside are re-checked for exact tiling.
"""
import sys, re, math
from decimal import Decimal, ROUND_FLOOR, ROUND_CEILING, localcontext

def exact(x):
    """The exact value of the binary double x.

    Converting through the shortest round-tripping decimal is not the same
    thing and is not safe here: that decimal can lie on the unsafe side of the
    bound before any directed rounding starts.  For the double nearest -0.3,
    whose exact value is -0.29999999999999998889..., the shortest form is
    "-0.3", so rounding it outward still prints an upper bound stronger than
    the one certified.  Decimal.from_float is exact, so every directed
    rounding below starts from the number the verifier actually produced."""
    if not math.isfinite(x):
        raise ValueError("nonfinite endpoint")
    return Decimal.from_float(x)

def quantized(d, exponent, rounding):
    # An exact double can need ~1080 digits; the default context would raise
    # rather than round it.
    with localcontext() as ctx:
        ctx.prec = 1100
        return d.quantize(Decimal(1).scaleb(exponent), rounding=rounding)

def floor_d(x, n):
    return quantized(exact(x), -n, ROUND_FLOOR)

def ceil_d(x, n):
    return quantized(exact(x), -n, ROUND_CEILING)

def sig_outward(x, nsig=3):
    """Round to nsig significant digits, outward: a negative bound rounds up
    (toward zero) and a positive bound rounds down, so the printed number is
    always implied by the certified one."""
    d = exact(x)
    if d == 0:
        return "0"
    r = quantized(d, d.adjusted() - nsig + 1,
                  ROUND_CEILING if d < 0 else ROUND_FLOOR)
    # The property the captions state, enforced rather than assumed.
    assert (d <= r < 0) if d < 0 else (0 < r <= d), f"{r} not implied by {d}"
    e = r.adjusted()
    if -3 <= e <= 3:
        return ("%g" % r)
    with localcontext() as ctx:
        ctx.prec = 1100
        mant = r.scaleb(-e)
    return "%g\\times 10^{%d}" % (mant, e)

def fmt_endpoint(x, upper):
    """Shortest decimal for a certified range endpoint, rounded INWARD.

    The enclosure columns round outward, because a printed bound must contain
    the computed one.  A range endpoint is the opposite: the printed range is
    a claim about every p in it, so it must be *contained* in what was
    certified.  Printing 15.92 for a range certified up to 15.920000000000002
    is therefore sound; printing it for one certified only up to the nearest
    double below 15.92 would not be, which is why the certificates store
    endpoints rounded outward from the stated decimals in the first place.

    The containment is tested on the exact decimal real number denoted by the
    string, not on the double nearest to it.  Those differ exactly when it
    matters: "15.92" rounds to a double strictly below 15.92, so a float test
    would accept it for a range certified only up to that double -- the very
    error this convention exists to prevent.  Decimal of a float is exact, and
    Decimal of a string is the decimal it denotes.

    Certificate endpoints are the intended decimals moved by at most one unit
    in the last place, so we look for the shortest decimal in the one-ulp
    window on the safe side of x.  p = 2.0000001 must not abbreviate to "2",
    which would assert a counterexample at p = 2, where C vanishes; the
    containment test rules that out automatically.  If no short decimal lies
    in the window -- which happens when a certificate has *not* been rounded
    outward -- we fall back to a safe truncation of x itself, which is ugly
    on purpose: it signals that the certificate needs widening.
    """
    xd = exact(x)
    lo_d = Decimal(math.nextafter(x, -math.inf)) if upper else xd
    hi_d = xd if upper else Decimal(math.nextafter(x, math.inf))
    mid = float(lo_d + (hi_d - lo_d) / 2)
    out = None
    for prec in range(1, 18):
        t = "%.*g" % (prec, mid)
        if not (lo_d <= Decimal(t) <= hi_d):
            continue
        if "e" in t or "E" in t:          # keep a plain decimal if one exists
            out = out or t
            continue
        out = t
        break
    if out is None:
        q = Decimal(1).scaleb(xd.adjusted() - 16)
        out = str(xd.quantize(q, rounding=ROUND_FLOOR if upper else ROUND_CEILING))
    # The property the manuscript states, enforced rather than assumed.
    if upper:
        assert Decimal(out) <= xd, f"printed {out} exceeds certified {xd}"
    else:
        assert Decimal(out) >= xd, f"printed {out} is below certified {xd}"
    return out

def fmt_range(a, b):
    return "[%s,\\, %s]" % (fmt_endpoint(a, False), fmt_endpoint(b, True))

def fmt_lam_lo(x):
    """Round down, keeping enough digits to stay strictly above 1."""
    for n in range(4, 16):
        v = floor_d(x, n)
        if v > 1:
            assert v <= exact(x)
            return str(v)
    return str(exact(x))

def fmt_lam_hi(x):
    """Round up, keeping enough digits to stay strictly below 2.  lambda*(p)
    tends to 2 as p tends to 2, so four decimals would print 2.0000 and
    contradict the requirement that Lambda lie inside (1,2)."""
    for n in range(4, 16):
        v = ceil_d(x, n)
        if v < 2:
            assert v >= exact(x)
            return str(v)
    return str(exact(x))

def main():
    rows = []
    for path in sys.argv[1:]:
        base = path.split("/")[-1]
        # The dimensions come from the file name; the range is taken from the
        # data, so that files whose name does not encode it (final_d3_lower
        # and friends) are handled too.
        mm = (re.match(r"range_d(\d+)_[0-9.]+_[0-9.]+\.txt$", base)
              or re.match(r"certified_d(\d+)_\w+\.txt$", base))
        ml = (re.match(r"lift_(\d+)to(\d+)_[0-9.]+_[0-9.]+\.txt$", base)
              or re.match(r"certified_lift_(\d+)to(\d+)\.txt$", base))
        if ml:
            src, d = int(ml.group(1)), int(ml.group(2))
        elif mm:
            src, d = None, int(mm.group(1))
        else:
            print(f"% skipped {path}", file=sys.stderr)
            continue
        a = b = None
        boxes = []
        for line in open(path):
            f = line.split()
            if f and f[0] == "OK":
                boxes.append(tuple(float(x) for x in f[1:7]) + (f[7],))
        if not boxes:
            print(f"% no boxes in {path}", file=sys.stderr)
            continue
        boxes.sort()
        if a is None:
            a, b = boxes[0][0], boxes[-1][1]
        # exact tiling check
        assert boxes[0][0] == a and boxes[-1][1] == b, f"{path}: endpoints"
        for i in range(len(boxes) - 1):
            assert boxes[i][1] == boxes[i + 1][0], f"{path}: gap at {boxes[i][1]!r}"
        signs = set(x[6] for x in boxes)
        assert len(signs) == 1, f"{path}: mixed signs"
        sign = signs.pop()

        lam_lo = min(x[2] for x in boxes)
        lam_hi = max(x[3] for x in boxes)
        # worst case bound on C: the value closest to zero over all boxes
        if sign == "-":
            worst = max(x[5] for x in boxes)          # least negative upper bound
            cbound = "\\C(p)\\leq %s" % sig_outward(worst)
        else:
            worst = min(x[4] for x in boxes)          # least positive lower bound
            cbound = "\\C(p)\\geq %s" % sig_outward(worst)

        rows.append((d, a, b, len(boxes), fmt_lam_lo(lam_lo), fmt_lam_hi(lam_hi),
                     cbound, sign, src))

    rows.sort(key=lambda r: (str(r[8]), r[0], r[1]))
    for d, a, b, n, l1, l2, cb, sign, src in rows:
        head = ("$%d$ & $%d$ & " % (src, d)) if src else ("$%d$ & " % d)
        if src:
            cb = cb.replace("\\C(p)", "\\C_{%d\\to %d}(p)" % (src, d))
        print("%s$%s$ & $%d$ & $[%s,\\, %s]$ & $%s$ & $%s$ \\\\"
              % (head, fmt_range(a, b), n, l1, l2, cb, sign))

if __name__ == "__main__":
    main()
