The search -- not part of the proof.
====================================

The proof is ../verify.cpp together with the certificate files in ../boxes/.
That verifier takes the boxes as given and re-proves every claim about them in
interval arithmetic.  It contains no search, no heuristics and no floating
point shooting.

The program here is how the boxes in ../boxes/ were found.  It does ordinary
double precision shooting to propose an eigenvalue interval, sweeps a ladder of
paddings and initial times until a box certifies, and bisects the exponent box
when none of them does.  None of that is trusted.  A bug here can propose a box
that is wrong, but a wrong box cannot pass the verifier: it would fail one of
(C1)-(C6) and the verifier would stop with a nonzero exit status.  This is the
whole reason the search is a separate program.

To regenerate ../boxes/ from scratch:

    make
    ./reproduce.sh 16          # number of processes; about 30 core-hours

The certificate files this produces need not agree line for line with the
archived ones: a different compiler or processor can move a padding by one unit
in the last place and take the subdivision down a different path.  What must
hold, and what the paper claims, is that the archived boxes verify.  That is
checked by ../verify_all.sh and is reproducible exactly.

Usage of the search program directly:

    ./certify_range d a b [minwidth] [startwidth] [dtarget]

searches the exponent range [a,b] in dimension d, starting from boxes of width
startwidth and bisecting down to minwidth, and prints certificate body lines to
standard output.  With dtarget > d the lifted defect C_{d->dtarget} is the one
whose sign is certified.  Diagnostics go to standard error.

shooting.py is the original non-rigorous numerics of Section 3 of the paper,
kept for reference.  It produced the figure and the approximate values of p*_d,
and is part of neither the search nor the proof.
