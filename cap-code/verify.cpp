// verify.cpp
//
// Verification of the certified boxes for Theorems 1.1 (direct) and 1.2
// (lifted) of the accompanying paper.
//
// This program proves nothing new; it re-proves, in interval arithmetic, the
// claims (C1)-(C6) of Section 4 for the boxes listed in a certificate file.
// It performs NO search of any kind.  Every quantity that a search would have
// to choose -- the exponent box P = [p1,p2], the eigenvalue box
// Lambda = [lam1,lam2], and the matching time tau0 -- is read from the file.
// How that list was produced is irrelevant to the proof (see search/), and
// nothing in it is trusted: a wrong line simply fails to verify.
//
// For each line of the file the program verifies, for the box P x Lambda:
//
//   (C1)  kappa_p(lambda) > 0 and the hypotheses (H1)-(H3) of the pole lemma
//         (Lemma 3.1) hold at theta0 = (pi/2) tau0, on P x {lam1},
//         P x {lam2} and P x Lambda;
//   (C2)  f'(pi/2) < 0 on P x {lam1} and f'(pi/2) > 0 on P x {lam2}, so that
//         for each fixed p in P the intermediate value theorem gives an
//         eigenvalue lam*(p) in Lambda;
//   (C3)  the solution exists on [tau0,1] over the whole box, f(pi/2) < 0,
//         and the tube gives bounds Flo <= f <= Fhi and an enclosure Jend of
//         the mean integral;
//   (C4)  Lambda is contained in (1,2);
//   (C5)  the resulting enclosure C of the mean value defect (of the lift
//         C_{d->dt} when dt > d) does not contain 0.
//
// and then, over the whole file,
//
//   (C6)  the exponent boxes tile the claimed range [a,b] exactly, as
//         equalities between double precision numbers, and the sign of C is
//         the same for every box.
//
// Certificate file format.  Blank lines and lines beginning with '#' are
// ignored.  One header line
//
//     claim  d  dtarget  a  b  sign
//
// states the theorem being verified: for every p in [a,b] the lifted defect
// C_{d->dtarget}(p) has the stated sign (dtarget = d is the direct case).
// It is followed by one line per box,
//
//     p1  p2  lam1  lam2  tau0
//
// in increasing order of p1.
//
// Usage:  ./verify [--lifting-tail] certificate.txt [first last]
//
// With --lifting-tail, the claimed sign must be negative and an additional
// check (C7) proves that the spherical mean is strictly negative. The lifting
// monotonicity lemma in lifting/lifting.tex then extends each verified box to
// every integer target dimension n >= dtarget. MEAN columns record the bounds.
//
// With no range, the whole file is verified and the certified statement is
// printed.  With a range (1-based, inclusive) only those boxes are verified,
// so that one file can be split across processes; (C6) is checked in full
// either way, since it costs nothing.
// Output: one line "OK p1 p2 lam1 lam2 C_lo C_hi sign" per verified box,
//         then the verified claim.  Exit status 0 if and only if every check
//         passed; the first failure is reported on stderr and stops the run.
//
// Trusted: the CAPD library, this file, and Lemma 3.1 of the paper.

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <cfenv>
#include "capd/capdlib.h"

using capd::interval;
using capd::IMap;
using capd::IOdeSolver;
using capd::ITimeMap;
using capd::IVector;
using capd::C0Rect2Set;

static const int ORDER = 20;
static const int NSUB  = 32;   // subintervals per step for the tube bounds

// Shortest decimal string that reads back as exactly x.  Used only for the
// human readable summary; the certified numbers are printed in full.
static std::string shortest(double x) {
  char buf[32];
  for (int prec = 6; prec < 17; ++prec) {
    std::snprintf(buf, sizeof buf, "%.*g", prec, x);
    if (std::strtod(buf, nullptr) == x) return buf;
  }
  std::snprintf(buf, sizeof buf, "%.17g", x);
  return buf;
}
      // Taylor order of the rigorous solver

// ---------------------------------------------------------------------------
// Interval quantities entering the defect
// ---------------------------------------------------------------------------

// kappa_p(lambda) = lambda ((p-1) lambda + d - p)
static interval kappaOf(const interval& lam, int d, const interval& p) {
  return lam * ((p - 1) * lam + (interval(d) - p));
}

// sigma_d exactly, from sigma_2 = 1/pi, sigma_3 = 1/2, sigma_{m+2} = m/(m-1) sigma_m
static interval sigmaD(int d) {
  int m = (d % 2 == 0) ? 2 : 3;
  interval s = (d % 2 == 0) ? interval(1) / interval::pi() : interval(1) / interval(2);
  while (m < d) { s *= interval(m) / interval(m - 1); m += 2; }
  return s;
}

// Rigorous enclosure of log Gamma(x) for an interval x with inf(x) > 0.
//
// Uses the recurrence log Gamma(x) = log Gamma(x+1) - log x to push the
// argument above 20, then Stirling's series
//   log Gamma(z) = (z-1/2) log z - z + log(2 pi)/2
//                  + sum_{n=1}^{5} B_{2n} / (2n(2n-1) z^{2n-1}) + R,
// whose remainder for real z > 0 is bounded in absolute value by the first
// omitted term, |R| <= |B_12| / (12 * 11 * z^11) = 691 / (360360 z^11).
// At z >= 20 this is below 1e-18, so the enclosure is essentially sharp.
static interval lgammaEncl(interval x) {
  interval acc(0.0);
  int guard = 0;
  while (x.leftBound() < 20.0) {
    acc -= log(x);
    x = x + 1;
    if (++guard > 100) throw std::runtime_error("lgamma: shift did not terminate");
  }
  interval iz = interval(1.0) / x, iz2 = sqr(iz);
  interval p1 = iz;                    // z^-1
  interval p3 = p1 * iz2, p5 = p3 * iz2, p7 = p5 * iz2, p9 = p7 * iz2;
  interval p11 = p9 * iz2;
  interval ser = p1 / 12 - p3 / 360 + p5 / 1260 - p7 / 1680 + p9 / 1188;
  double rb = (interval(691) / 360360 * p11).rightBound();
  interval rem(-rb, rb);
  return acc + (x - interval(1) / 2) * log(x) - x
       + log(2 * interval::pi()) / 2 + ser + rem;
}

// H_{m->dt}(lam) = Gamma(dt/2) Gamma((m+lam)/2) / (Gamma(m/2) Gamma((dt+lam)/2)).
// For dt-m even this telescopes to an exact rational function of lam and no
// Gamma evaluation is needed; otherwise it is computed from lgammaEncl.
// Equals 1 when dt == m, which is the direct case.
static interval Hlift(int m, int dt, const interval& lam) {
  if ((dt - m) % 2 == 0) {
    interval h(1.0);
    for (int j = 0; j < (dt - m) / 2; ++j)
      h *= interval(m + 2 * j) / (interval(m + 2 * j) + lam);
    return h;
  }
  // H is strictly decreasing in lam: (d/dlam) log H
  //   = (psi((m+lam)/2) - psi((dt+lam)/2)) / 2 < 0,
  // because psi is increasing and (m+lam)/2 < (dt+lam)/2 for m < dt.  We
  // therefore evaluate at the two endpoints of lam and take the hull.  This
  // matters: evaluating lgammaEncl on a wide lam treats the two log-Gamma
  // terms as independent and grossly overestimates.
  interval c = lgammaEncl(interval(dt) / 2) - lgammaEncl(interval(m) / 2);
  interval hiEnd = exp(c + lgammaEncl((interval(m) + lam.leftBound()) / 2)
                         - lgammaEncl((interval(dt) + lam.leftBound()) / 2));
  interval loEnd = exp(c + lgammaEncl((interval(m) + lam.rightBound()) / 2)
                         - lgammaEncl((interval(dt) + lam.rightBound()) / 2));
  return interval(loEnd.leftBound(), hiEnd.rightBound());
}

// ---------------------------------------------------------------------------
// (C1): the pole lemma at theta0
// ---------------------------------------------------------------------------

struct LocalData { interval f0, g0, J0, fLoc; };

// Lemma 3.1 with p and lambda both interval valued: verifies (H1)-(H3) and
// returns the enclosures of f(theta0), f'(theta0), the integral over
// [0,theta0] and f on [0,theta0].  Throws if a hypothesis fails.
static LocalData localLemma(int d, const interval& p, const interval& lam,
                            const interval& theta0) {
  interval a   = interval(d - 1) / interval(d);
  interval eps(theta0.rightBound());
  interval kap = kappaOf(lam, d, p);
  if (!(kap.leftBound() > 0))     throw std::runtime_error("kappa <= 0");
  if (!(eps.rightBound() <= 0.5)) throw std::runtime_error("eps > 1/2");

  interval m = 2 * kap / d;
  interval la2 = sqr(lam), a2 = sqr(a), a3 = a2 * a;
  interval beta = 4 * (p - 2) * sqr(m) / (la2 * a2);
  interval c0 = (d - 2) * (interval(34) / 100 + beta) * m
              + 2 * kap * (p - 2) * sqr(m) / (la2 * a) + 2 * (p - 2) * sqr(m) / a;
  interval gf = 16 * (p - 2) * sqr(m) / (la2 * a3);
  interval gg = 8 * (p - 2) * m / (la2 * a2);
  interval cf = (d - 2) * m * gf + 2 * kap * beta + 4 * (p - 2) * sqr(m) / a2;
  interval Lf = kap + cf * sqr(eps);
  interval Lg = (d - 2) * (interval(34) / 100 + beta) + (d - 2) * m * gg
              + interval(3) / 2 * a * kap * gg + 4 * (p - 2) * m / a;

  if (!((kap * sqr(eps)).rightBound() <= (a * d / 2).leftBound()))
    throw std::runtime_error("(H1)");
  if (!((sqr(kap) * sqr(eps) / (d * (d - 1)) + c0 * sqr(eps) / (d + 1)).rightBound()
        < (kap / d).leftBound()))
    throw std::runtime_error("(H2)");
  interval lip = sqrt(sqr(Lf) + sqr(Lg * eps));
  if (!((2 * eps * lip).rightBound() < interval(d - 1).leftBound()))
    throw std::runtime_error("(H3)");

  interval ef = sqr(kap) * sqr(eps) * sqr(theta0) / (2 * d * (d - 1))
              + c0 * power(theta0, 4) / (4 * (d + 1));
  interval eg = sqr(kap) * sqr(eps) * theta0 / (d * (d - 1))
              + c0 * power(theta0, 3) / (d + 1);
  interval rf(-ef.rightBound(), ef.rightBound());
  interval rg(-eg.rightBound(), eg.rightBound());

  LocalData L;
  L.f0 = a - kap * sqr(theta0) / (2 * d) + rf;
  L.g0 = -kap * theta0 / d + rg;
  interval th2(0, sqr(theta0).rightBound());
  interval fRange = a - kap * th2 / (2 * d) + rf;
  L.fLoc = interval(fRange.leftBound(), a.rightBound());
  interval S = power(theta0, d - 1) / (d - 1)
             * interval(power(1 - sqr(theta0) / 6, d - 2).leftBound(), 1);
  L.J0 = L.fLoc * S;
  return L;
}

// ---------------------------------------------------------------------------
// (C2), (C3): rigorous integration from tau0 to 1
// ---------------------------------------------------------------------------

// State (f,g,J,lam,q) with lam' = q' = 0; q carries the exponent p and the
// time is rescaled by theta = H tau, H = pi/2, so that the equator is tau = 1.
static std::string makeField(int d) {
  std::ostringstream sp;
  for (int i = 0; i < d - 2; ++i) sp << (i ? "*" : "") << "sin(H*t)";
  std::ostringstream o;
  o << "time:t;var:f,g,J,lam,q;par:H;fun:"
    << "H*g,"
    << "-H*(((" << d - 2 << ")*(cos(H*t)/sin(H*t))*g"
    << "+(lam*((q-1)*lam+(" << d << "-q)))*f)"
    << "*(g^2+lam^2*f^2)+(q-2)*lam^2*f*g^2)"
    << "/((q-1)*g^2+lam^2*f^2),"
    << "H*(" << sp.str() << ")*f,"
    << "0,0;";
  return o.str();
}

struct RunOut { interval fEnd, gEnd, JEnd; double fMin, fMax; };

// The bounds fMin, fMax on f over [tau0,1] are read off the solution curve on
// NSUB subintervals of each accepted step; no monotonicity of f is assumed.
// On each subinterval the curve is enclosed twice and the two are intersected:
// directly, and through the mean value form
//     f(t) in f(m) + (df/dt)(sub) * (sub - m),    m in sub,
// which is legitimate because df/dt = H g along every trajectory of the tube,
// so for a fixed trajectory the mean value theorem supplies some xi in sub with
// f(t) = f(m) + H g(xi) (t - m).  Direct evaluation is first order in the
// subinterval length, the mean value form second order, and the intersection is
// contained in each, so this can only tighten the enclosure.  It matters near an
// exponent where C is small: there the width of min f is what decides whether
// the enclosure of C excludes 0.
static RunOut rigorousRun(int d, const interval& p, const interval& lam,
                          const LocalData& L, double tau0) {
  IMap vf(makeField(d));
  vf.setParameter("H", interval::pi() / 2);
  IOdeSolver solver(vf, ORDER);
  ITimeMap tm(solver);

  IVector x0(5);
  x0[0] = L.f0; x0[1] = L.g0; x0[2] = L.J0; x0[3] = lam; x0[4] = p;
  C0Rect2Set s(x0, interval(tau0));

  tm.stopAfterStep(true);
  double fMin = 1e300, fMax = -1e300;
  do {
    tm(1.0, s);
    interval step = solver.getStep();
    const IOdeSolver::SolutionCurve& curve = solver.getCurve();
    interval domain = interval(0, 1) * step;
    for (int i = 0; i < NSUB; ++i) {
      interval sub = interval(i, i + 1) * step / NSUB;
      intersection(domain, sub, sub);
      IVector v = curve(sub);
      interval fe = v[0];
      interval m(sub.mid());
      intersection(m, sub, m);
      interval mvf = curve(m)[0] + (interval::pi() / 2) * v[1] * (sub - m);
      interval both;
      if (intersection(fe, mvf, both)) fe = both;
      fMin = std::min(fMin, fe.leftBound());
      fMax = std::max(fMax, fe.rightBound());
    }
  } while (!tm.completed());

  IVector x = IVector(s);
  RunOut r; r.fEnd = x[0]; r.gEnd = x[1]; r.JEnd = x[2];
  r.fMin = fMin; r.fMax = fMax;
  return r;
}

// ---------------------------------------------------------------------------
// One box: (C1)-(C5)
// ---------------------------------------------------------------------------

// Verifies (C1)-(C5) for P x Lambda at the given tau0 and returns the
// enclosure of the (lifted) mean value defect.  Throws on any failed check.
static interval certifyBox(int d, int dt, const interval& P,
                           double lo, double hi, double tau0,
                           interval* sourceMean = nullptr) {
  // (C4)
  if (!(lo > 1.0 && hi < 2.0))
    throw std::runtime_error("(C4) Lambda not inside (1,2)");
  interval Lam(lo, hi);
  interval theta0 = interval::pi() / 2 * tau0;

  // (C1) at the two endpoints of Lambda and on the whole box
  LocalData La = localLemma(d, P, interval(lo), theta0);
  LocalData Lb = localLemma(d, P, interval(hi), theta0);
  LocalData LB = localLemma(d, P, Lam, theta0);

  // (C2) opposite signs of f'(pi/2) at the two endpoints, for every p in P
  RunOut ra = rigorousRun(d, P, interval(lo), La, tau0);
  if (!(ra.gEnd.rightBound() < 0))
    throw std::runtime_error("(C2) f'(pi/2) at lam1 not < 0");
  RunOut rb = rigorousRun(d, P, interval(hi), Lb, tau0);
  if (!(rb.gEnd.leftBound() > 0))
    throw std::runtime_error("(C2) f'(pi/2) at lam2 not > 0");

  // (C3) the whole box, giving the tube bounds and the mean integral
  RunOut r = rigorousRun(d, P, Lam, LB, tau0);
  if (!(r.fEnd.rightBound() < 0))
    throw std::runtime_error("(C3) f(pi/2) not < 0");

  // (C5) the mean value defect.  max_S v is attained either at theta = 0,
  // where f = a, or inside the tube; min_S v either at theta = pi/2 or
  // inside the tube.
  interval a = interval(d - 1) / interval(d);
  interval maxv(a.leftBound(), std::max(a.rightBound(), r.fMax));
  interval minv(std::min(LB.fLoc.leftBound(), r.fMin), r.fEnd.rightBound());
  interval mean = 2 * sigmaD(d) * r.JEnd;      // average of v over S^{d-1}
  if (sourceMean) *sourceMean = mean;
  interval C = (P - 2) / 2 * (maxv + minv)
             + interval(dt) * (dt + 2) / (Lam + dt) * Hlift(d, dt, Lam) * mean;
  if (!(C.leftBound() > 0 || C.rightBound() < 0))
    throw std::runtime_error("(C5) enclosure of C contains 0");
  return C;
}

// ---------------------------------------------------------------------------
// The certificate file, and (C6)
// ---------------------------------------------------------------------------

struct Box { double p1, p2, lam1, lam2, tau0; };

int main(int argc, char** argv) {
  // Parse certificate decimals in the usual round-to-nearest mode, regardless
  // of the mode left by library static initialization.
  std::fesetround(FE_TONEAREST);
  bool liftingTail = argc > 1 && std::string(argv[1]) == "--lifting-tail";
  if (liftingTail) { --argc; ++argv; }
  if (argc != 2 && argc != 4) {
    std::cerr << "usage: verify [--lifting-tail] certificate.txt [first last]\n";
    return 2;
  }
  std::ifstream in(argv[1]);
  if (!in) { std::cerr << "cannot open " << argv[1] << "\n"; return 2; }

  int d = 0, dt = 0;
  double a = 0, b = 0;
  std::string sign;
  bool haveClaim = false;
  std::vector<Box> boxes;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream ls(line);
    std::string junk;
    // Prefix dispatch, so the keyword is checked for equality below:
    // "claims ..." must not be read as a claim line.
    if (line.compare(0, 5, "claim") == 0) {
      if (haveClaim) { std::cerr << "more than one claim line\n"; return 2; }
      std::string kw;
      if (!(ls >> kw >> d >> dt >> a >> b >> sign) || kw != "claim" || (ls >> junk))
        { std::cerr << "malformed claim line: " << line << "\n"; return 2; }
      haveClaim = true;
      continue;
    }
    Box x;
    if (!(ls >> x.p1 >> x.p2 >> x.lam1 >> x.lam2 >> x.tau0) || (ls >> junk))
      { std::cerr << "malformed box line: " << line << "\n"; return 2; }
    if (!(x.p1 < x.p2) || !(x.lam1 < x.lam2) || !(x.tau0 > 0))
      { std::cerr << "degenerate box: " << line << "\n"; return 2; }
    // Domain of Lemma 3.1: p >= 2.  Its constants carry factors of (p-2),
    // which change sign below 2, so the enclosures would not be valid there.
    // (d >= 3, kappa > 0 and 0 < T <= 1/2 are the lemma's other hypotheses;
    // the last two are checked in localLemma.)
    if (!(x.p1 >= 2.0))
      { std::cerr << "box below p = 2, outside the lemma: " << line << "\n"; return 2; }
    boxes.push_back(x);
  }
  if (!haveClaim) { std::cerr << "no claim line\n"; return 2; }
  if (boxes.empty()) { std::cerr << "no boxes\n"; return 2; }
  if (d < 3) { std::cerr << "d < 3, outside the lemma\n"; return 2; }
  if (dt < d) { std::cerr << "dtarget < d\n"; return 2; }
  if (sign != "+" && sign != "-") { std::cerr << "sign must be + or -\n"; return 2; }
  if (liftingTail && sign != "-") {
    std::cerr << "lifting-tail requires a negative defect claim\n"; return 2;
  }

  // (C6), tiling part.  These are equalities between double precision
  // numbers: a discrepancy of one unit in the last place would leave a
  // nonempty set of exponents uncertified.
  if (boxes.front().p1 != a) {
    std::cerr << "FAILED (C6): first box starts at " << shortest(boxes.front().p1)
              << ", not at " << shortest(a) << "\n";
    return 1;
  }
  if (boxes.back().p2 != b) {
    std::cerr << "FAILED (C6): last box ends at " << shortest(boxes.back().p2)
              << ", not at " << shortest(b) << "\n";
    return 1;
  }
  for (size_t i = 0; i + 1 < boxes.size(); ++i)
    if (boxes[i].p2 != boxes[i + 1].p1) {
      std::cerr << "FAILED (C6): gap or overlap at " << shortest(boxes[i].p2)
                << " -> " << shortest(boxes[i + 1].p1) << "\n";
      return 1;
    }

  // Boxes to verify in this run.  The (C6) check above covers the whole file
  // whichever range is given, so splitting a file across processes only
  // splits the per-box work; every process re-checks the tiling.
  size_t first = 1, last = boxes.size();
  if (argc == 4) {
    long f = std::atol(argv[2]), l = std::atol(argv[3]);
    if (f < 1 || l > (long)boxes.size() || f > l) {
      std::cerr << "box range " << f << ".." << l << " outside 1.."
                << boxes.size() << "\n";
      return 2;
    }
    first = f; last = l;
  }

  std::cout << std::setprecision(17);
  for (size_t i = first - 1; i < last; ++i) {
    const Box& x = boxes[i];
    interval C, mean;
    try {
      C = certifyBox(d, dt, interval(x.p1, x.p2), x.lam1, x.lam2, x.tau0, &mean);
      if (liftingTail && !(mean.rightBound() < 0))
        throw std::runtime_error("(C7) spherical mean not strictly negative");
    } catch (std::exception& e) {
      std::cerr << "FAILED box " << i + 1 << " [" << x.p1 << ", " << x.p2
                << "]: " << e.what() << "\n";
      return 1;
    }
    std::string s = (C.leftBound() > 0) ? "+" : "-";
    if (s != sign) {                                    // (C6), sign part
      std::cerr << "FAILED box " << i + 1 << " [" << x.p1 << ", " << x.p2
                << "]: sign " << s << " does not match the claimed " << sign << "\n";
      return 1;
    }
    // Directed rounding is used by the arithmetic above. Decimal formatting
    // must instead round to nearest so 17 digits round-trip to these doubles.
    int arithmeticRounding = std::fegetround();
    std::fesetround(FE_TONEAREST);
    std::cout << "OK " << x.p1 << " " << x.p2 << " " << x.lam1 << " " << x.lam2
              << " " << C.leftBound() << " " << C.rightBound() << " " << s;
    if (liftingTail)
      std::cout << " MEAN " << mean.leftBound() << " " << mean.rightBound();
    std::cout << "\n" << std::flush;
    std::fesetround(arithmeticRounding);
  }

  std::fesetround(FE_TONEAREST);
  if (first != 1 || last != boxes.size()) {
    std::cout << "PARTIAL boxes " << first << ".." << last << " of "
              << boxes.size() << " verified; (C6) holds for the file.\n";
    return 0;
  }
  std::cout << "VERIFIED " << boxes.size() << " boxes tiling ["
            << shortest(a) << ", " << shortest(b)
            << "]: for every p in that range there is an eigenvalue "
            << "lambda*(p) in (1,2) with C";
  if (dt != d) std::cout << "_{" << d << "->" << dt << "}";
  std::cout << "(p) " << (sign == "+" ? ">" : "<") << " 0.\n";
  if (liftingTail)
    std::cout << "VERIFIED TAIL: the spherical mean is negative and the lifted "
              << "defect is negative for every integer target dimension >= "
              << dt << ".\n";
  return 0;
}
