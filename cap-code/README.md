# Computer-assisted proof: certified counterexamples to the pointwise asymptotic mean value property

This is the computer-assisted proof behind the two certified theorems of the
paper — the direct counterexamples of Section 3.2 (Theorem 3.2, Table 1) and
the lifted ones of Section 3.3 (Theorem 3.10, Table 2). It is self-contained:
the only external dependency is the CAPD library.

The mathematical content of the computation is the list of claims (C1)–(C6) of
Section 3.2, with (C5) replaced by (C5$'$) for the lifted certificates and one
further claim (C7) for two of them, together with the pole lemma (Lemma 3.5:
the enclosure of the regular solution at the singular endpoint $\theta=0$,
proved in Appendix B). Everything here either checks those claims or explains
how the parameter boxes they refer to were found.

## What is being proved

For each dimension $d$ and each range $[p_-,p_+]$ of Table 1, and for **every**
real $p$ in that range, there is a $p$-harmonic function on $\mathbb R^d$,
homogeneous of degree $\lambda\in(1,2)$ with a critical point at the origin,
whose mean value defect $\mathcal C(p)$ is nonzero with the sign recorded in
the table. Table 2 is the same statement for the cylindrical lifts, with the
lifted defect $\mathcal C_{m\to d}(p)$ in place of $\mathcal C(p)$. Since
$\lambda<2$, a nonzero defect makes $\varepsilon^{-2}\mathcal M_\varepsilon
u(0)$ diverge, so the pointwise expansion fails at the origin. Together the two
tables leave no gap in $[2.0000001,100]$ for $d=4,\dots,10$ (Corollary 3.11).

Two of the certificates carry the extra claim (C7), that the spherical mean of
the source profile is strictly negative. With it, the monotonicity of the lift
in the target dimension (Lemma 3.12, which needs no computation) extends the
conclusion to **every** dimension $d\geq10$ for $p\in[2.0000001,76.5]$
(Corollary 3.13), from the single dimension-10 certificate and the $4\to10$
lift.

The 16 certificates named `d<d>_lower` and `d<d>_upper` are the rows of
Table 1; the 7 named `lift_<m>to<d>` are the rows of Table 2.

The ranges are certified, not sampled. Each is tiled exactly by finitely many
closed exponent boxes, and every real exponent lies in one of them.

## Layout

```
verify.cpp        the verifier: proves (C1)-(C6) for one certificate file
Makefile          builds it
verify_all.sh     runs it over every certificate and prints the theorems
gen_table.py      turns the verifier's output into the LaTeX table rows
boxes/            23 certificate files: the parameter boxes, nothing else
capd/             the pinned CAPD source the computation was built against
search/           how the boxes were found -- NOT part of the proof
SHA256SUMS        checksums of every file above
```

## Running it

CAPD is the one external dependency, and it is pinned: `capd/` holds the exact
source the certified computation was built against, so nothing here depends on
an upstream repository still offering that revision. Build it, then the
verifier:

```
tar xzf capd/capd-7310792.tar.gz
cd capd-7310792 && mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/.local/capd -DCMAKE_BUILD_TYPE=Release
make && make install
cd ../..

make                       # set CAPD_DIR=... if CAPD is not in ~/.local/capd
./verify_all.sh 16         # number of processes
```

That is CAPD commit `7310792`, version 6.1.0, with its bundled **filib**
interval backend and directed rounding; `capd/PROVENANCE.txt` records the
commit, the checksums, how the archive was produced and how to re-derive it
from upstream.

Exit status is `0` if and only if every box of every file passed every check.
The script prints the verified statements and then the two table bodies, which
are the rows appearing in the paper. A single certificate can also be checked
on its own,

```
./verify boxes/d3_lower.txt             # the whole file
./verify boxes/d3_lower.txt 1 100       # boxes 1..100, for splitting the work
```

and a box range still checks (C6) over the whole file, that being a property of
the file rather than of any box.

Verifying all 55472 boxes takes about five and a half hours of processor time,
about a quarter of an hour on twenty-two cores. Finding them in the first place
took of the order of 30 core-hours; see `search/`.

## The certificate format

A certificate file is a header line naming what is to be proved, followed by
one line per exponent box. Blank lines and lines beginning with `#` are
ignored.

```
claim  d  dtarget  p_-  p_+  sign
p1  p2  lambda1  lambda2  tau0
...
```

The claim line says: for every $p\in[p_-,p_+]$ the defect
$\mathcal C_{d\to dtarget}(p)$ has the stated sign, with `dtarget = d` meaning
the direct case. Each box line gives the exponent box $P=[p_1,p_2]$, the
eigenvalue box $\Lambda=[\lambda_1,\lambda_2]$ within which the shooting
eigenvalue is trapped, and the rescaled time $\tau_0$ at which Lemma 3.5 is
applied. Boxes are in increasing order of $p_1$. For a lifted certificate the
eigenvalue problem, and hence the whole integration, is the one of the source
dimension `d`, and only the defect is computed in `dtarget`; this is why the
same verifier and the same box format serve both tables.

Nothing in a certificate is trusted. The verifier reads the boxes and proves
every claim about them from scratch; a box on which any claim fails is
rejected and the run stops with a nonzero status. How the list was produced is
therefore irrelevant to the proof, which is why the search is a separate
program in `search/` and is not part of the trusted base.

Two conventions in the files are worth knowing about. First, the tiling
equalities that (C6) checks are *exact* equalities between double precision
numbers: one unit in the last place between $\sup P_i$ and $\inf P_{i+1}$ would
leave a nonempty set of exponents uncertified, so the verifier compares them
rather than trusting the subdivision. Second, range endpoints are stored
rounded **outward** from the decimals the paper prints — downward at $p_-$,
upward at $p_+$ — because a decimal such as `15.94257` is not a double and the
nearest one can fall inside the claimed range. This is why some endpoints in
`boxes/` and in `search/reproduce.sh` are written out to seventeen digits. The
table generator enforces the matching convention in the other direction: it
rounds a printed range endpoint **inward** and asserts that the decimal it
prints is contained in what was certified.

## Two details the paper leaves here

**Enclosing the extrema of `f`.** Claim (C3) needs bounds `Flo <= f <= Fhi`
valid along the whole tube, and these decide whether (C5) excludes 0. They are
taken from the solution curve on 32 subintervals of each accepted step. On a
subinterval `s` the curve is enclosed twice — directly, and through the mean
value form

    f(t) in f(m) + (df/dt)(s) * (s - m),   m in s,

which is available because `df/dt = H*g` along every trajectory of the tube —
and the two are intersected. Direct evaluation is first order in the length of
`s` and the mean value form is second order; the intersection is contained in
each of them, so the refinement is sound whatever the two enclosures happen to
be, and can only narrow the result.

This is not cosmetic. The effect is one-sided, so it is felt asymmetrically: on
a range where `C < 0` the binding requirement is `sup C < 0`, governed by
`sup Fend`, which is tight, whereas where `C > 0` it is `inf C > 0` that binds,
governed by `Flo`, which is not. With the first-order bound alone the two sides
of the `d = 3` gap closed at very different rates; the refinement removes the
discrepancy and brings that gap down from `4e-2` to `4e-5`. Note that no
monotonicity of `f` is assumed anywhere: the extreme values are enclosed from
the attained values `f(0) = a` and `f(pi/2)` on one side and the tube on the
other.

**Evaluating `H_{m->d}`.** When `d - m` is even it is the rational function of
the paper and interval arithmetic evaluates it directly. For `m = 3, d = 4` it
is not, and we use

    log Gamma(z) = (z - 1/2) log z - z + (1/2) log 2pi
                   + sum_{n=1}^{5} B_{2n} / (2n(2n-1) z^{2n-1}) + R(z),
    |R(z)| <= |B_12| / (132 z^11),

the classical remainder bound for real `z > 0`, after raising the argument
above 20 by `Gamma(z+1) = z Gamma(z)`; at `z >= 20` the remainder is below
1e-18. One point deserves care: `H_{m->d}` has to be enclosed over an interval
of `lambda`, and evaluating the two logarithms of `Gamma` separately there
treats them as independent and inflates the enclosure by orders of magnitude.
We instead use that `H_{m->d}` is strictly decreasing in `lambda`, since

    d/dlambda log H_{m->d} = (1/2) [psi((m+lambda)/2) - psi((d+lambda)/2)] < 0

by the monotonicity of the digamma function, and evaluate at the two endpoints.

## License

This package is distributed under the GNU General Public License, version 3 or
later; the text is in `LICENSE`. CAPD, whose source is pinned in `capd/`, is
itself under GPLv3 and carries its own `COPYING` inside that archive.

To keep the record exact, no license header was added to `verify.cpp` or to the
certificates: `verify.cpp` is byte-identical to the file that produced the
certified computation and was reverified unchanged in the independent rerun of
September 2026, and the `SHA256SUMS` in this package and in that audit refer to
those bytes. The `LICENSE` file covers the whole package.

## What is trusted

* the CAPD library (commit `7310792`, filib backend, directed rounding), whose
  source is pinned in `capd/`;
* `verify.cpp`, this file's few hundred lines;
* the pole lemma of the paper (Lemma 3.5), whose hypotheses (H1)–(H3)
  `verify.cpp` checks in interval arithmetic and whose conclusions it uses as
  the initial enclosure.

Not trusted, and not able to affect the result: everything in `search/`, the
choice of boxes, `gen_table.py`, `verify_all.sh` (a failure there cannot turn a
failed check into a passed one — it can only fail to run a check, which is why
it re-compares the verified boxes against the certificate line by line after
the run), and the non-rigorous numerics of Section 3 of the paper.

Outside the computation, the proof also uses Propositions 3.3 and 3.9 (a
spherical profile of the right shape yields a counterexample, and its lift
yields one in the higher dimension), Proposition 3.8 (the average of a
cylindrical lift), the standard equivalence of viscosity and weak solutions for
the $p$-Laplacian, and interior elliptic regularity. These are proved or cited
in the paper.

## A note on `verify.cpp`

The file is shipped byte-identical to the one that produced the archived
certification, so its header comments use the numbering of an earlier draft.
The correspondence is:

| in `verify.cpp` | in the paper |
| --- | --- |
| Theorem 1.1 (direct) | Theorem 3.2, Section 3.2 |
| Theorem 1.2 (lifted) | Theorem 3.10, Section 3.3 |
| claims (C1)–(C6) of Section 4 | Section 3.2, with (C5$'$) in Section 3.3 |
| Lemma 3.1 (the pole lemma) | Lemma 3.5 |

The `--lifting-tail` flag enables one extra claim, (C7): that the interval
enclosing the spherical mean of the source profile is strictly negative. The
mean is already formed inside (C5), so the check costs nothing beyond a sign
test. It is what the dimension-tail corollary of Section 3.3 rests on --- a
negative defect cannot be undone by lifting to a higher dimension once the
source mean is negative --- and `verify_all.sh` therefore passes the flag for
the two certificates that corollary uses, `boxes/d10_lower.txt` and
`boxes/lift_4to10.txt`. The flag also requires the claimed defect sign to be
negative, and refuses a certificate claiming a positive one.

## Reproducing the tables

`verify_all.sh` prints them at the end of a successful run. They can also be
regenerated from the assembled output of an earlier run:

```
python3 gen_table.py certified_d{3,4,5,6,7,8,9,10}_{lower,upper}.txt
python3 gen_table.py certified_lift_*.txt
```

All interval endpoints are rounded outward and all range endpoints inward, so
every number printed is implied by the computation. Every rounding starts from
the exact value of the binary double the verifier produced, via
`Decimal.from_float`, rather than from its shortest round-tripping decimal,
which can already sit on the unsafe side of the bound; each printed bound is
then checked against the exact one by an assertion. The rows are pasted into
the paper unedited.
