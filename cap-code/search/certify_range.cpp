// search/certify_range.cpp
//
// SEARCH, NOT PROOF.  This program looks for boxes on which the counterexample
// can be certified, and writes them to a certificate file.  Nothing here is
// part of the proof: the proof is what ../verify checks, and it re-verifies
// every claim from scratch given only the box list.  A bug in this file can
// produce a bad box, but a bad box cannot pass ../verify.  It is included so
// that the certificate files in ../boxes/ can be regenerated from scratch.
//
// It searches for boxes P x Lambda on which the counterexample of the paper
// certifies, for fixed dimension d and p in [a,b].  The exponent p enters the
// spherical ODE and the mean value defect only affinely, so it is carried as a
// state variable with p' = 0, exactly as lambda is.
//
// For a given exponent box P = [p1,p2] the eigenvalue box Lambda is proposed by
// ordinary double precision shooting at the two endpoints (Runge-Kutta 4 plus
// bisection), padded on both sides.  The padding must be large enough that the
// sign change (C2) holds for every exponent in P and small enough that the
// enclosure (C5) excludes 0, so a short ladder of paddings and of
// tau0 in {2^-10, 2^-12, 2^-14} is swept until one of them certifies.  If none
// does, the box is bisected and its halves retried, down to a minimum width.
// Bisection is exact in floating point, as is the initial cut, so the boxes
// tile [a,b]; ../verify checks this a posteriori and does not assume it.
//
// Cylindrical lifting.  With a target dimension dt > d (same parity), the
// same ODE in dimension d is used, but the defect is the lifted one,
//   C_{d->dt}(p) = (p-2)/2 (max v + min v)
//                  + dt(dt+2) H_{d->dt}(lam) / (lam+dt) * avg_{S^{d-1}} v,
//   H_{d->dt}(lam) = prod_{j=0}^{k-1} (d+2j)/(d+lam+2j),   k = (dt-d)/2,
// which is an exact rational function of lam when dt-d is even; for dt-d odd,
// H is computed from a validated log-Gamma (Stirling plus recurrence).
// Taking dt = d gives H = 1 and the direct case.
//
// Usage:  ./certify_range d a b [minwidth] [startwidth] [dtarget]
//
// [a,b] is first cut into pieces of length startwidth (a starting guess at
// the certifiable width; this avoids burning the whole (pad,tau0) ladder on
// boxes that are obviously too wide), and each piece is then subdivided
// adaptively down to minwidth.
//
// Output: the body of a certificate file, one line per box,
//         p1 p2 lam1 lam2 tau0
// and "# FAIL p1 p2" for a subinterval below minwidth that did not certify.
// Prepend a "claim d dtarget a b sign" line and pass the result to ../verify.
//
// NOTE: C(p) vanishes at p = 2 and at the crossing p*_d, so no interval whose
// closure contains such a point can be certified; the run will report FAIL
// there.  This is a genuine feature of the problem, not a numerical artifact.

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <cstdlib>
#include "capd/capdlib.h"

using capd::interval;
using capd::IMap;
using capd::IOdeSolver;
using capd::ITimeMap;
using capd::IVector;
using capd::C0Rect2Set;

static const int ORDER = 20;
static const int NSUB  = 32;   // must match ../verify.cpp

// ---------------------------------------------------------------------------
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
// Equals 1 when dt == m.
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

struct LocalData { interval f0, g0, J0, fLoc; };

// Lemma 3.1 with p and lambda both interval valued; throws if (H1)-(H3) fail.
static LocalData localLemma(int d, const interval& p, const interval& lam,
                            const interval& theta0) {
  interval a   = interval(d - 1) / interval(d);
  interval eps(theta0.rightBound());
  interval kap = kappaOf(lam, d, p);
  if (!(kap.leftBound() > 0))     throw std::runtime_error("kappa<=0");
  if (!(eps.rightBound() <= 0.5)) throw std::runtime_error("eps>1/2");

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

// State (f,g,J,lam,q), lam' = q' = 0; q carries the exponent p.
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
    // Must match ../verify.cpp exactly: direct enclosure intersected with the
    // mean value form f(t) in f(m) + H g(sub) (sub - m).
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

// -------------------- nonrigorous seeding --------------------
static double rhsG(double th, double f, double g, double lam, int d, double p) {
  double kap = lam * ((p - 1) * lam + d - p);
  double q = g * g + lam * lam * f * f;
  double num = ((d - 2) * (std::cos(th) / std::sin(th)) * g + kap * f) * q
             + (p - 2) * lam * lam * f * g * g;
  return -num / ((p - 1) * g * g + lam * lam * f * f);
}
static double shootEnd(double lam, int d, double p) {
  const double th0 = 1e-3;
  double kap = lam * ((p - 1) * lam + d - p), a = 1.0 - 1.0 / d;
  double f = a - kap * th0 * th0 / (2 * d), g = -kap * th0 / d, th = th0;
  const int N = 50000; const double h = (M_PI / 2 - th0) / N;
  for (int i = 0; i < N; ++i) {
    double k1f=g, k1g=rhsG(th,f,g,lam,d,p);
    double k2f=g+h/2*k1g, k2g=rhsG(th+h/2,f+h/2*k1f,g+h/2*k1g,lam,d,p);
    double k3f=g+h/2*k2g, k3g=rhsG(th+h/2,f+h/2*k2f,g+h/2*k2g,lam,d,p);
    double k4f=g+h*k3g,   k4g=rhsG(th+h,f+h*k3f,g+h*k3g,lam,d,p);
    f += h/6*(k1f+2*k2f+2*k3f+k4f); g += h/6*(k1g+2*k2g+2*k3g+k4g); th += h;
  }
  return g;
}
static double seedLambda(int d, double p) {
  double lo = 4.0/3 + 1e-12, hi = 2.0, ghi = shootEnd(hi - 1e-12, d, p);
  for (int i = 0; i < 60; ++i) {
    double mid = (lo+hi)/2, gm = shootEnd(mid, d, p);
    if ((gm > 0) == (ghi > 0)) { hi = mid; ghi = gm; } else lo = mid;
  }
  return (lo+hi)/2;
}

// -------------------- one attempt on a box --------------------
struct Cert { double lamLo, lamHi, cLo, cHi, tau0; };

static bool attempt(int d, int dt, double p1, double p2,
                    double lhat1, double lhat2,
                    double pad, double tau0, Cert& out) {
  static const bool VERBOSE = std::getenv("CERTIFY_VERBOSE") != nullptr;
  interval P(p1, p2);
  double lo = std::min(lhat1, lhat2) - pad, hi = std::max(lhat1, lhat2) + pad;
  if (!(lo > 1.0 && hi < 2.0)) {
    if (VERBOSE) std::cerr << "    pad=" << pad << " tau0=" << tau0
                           << " : Lambda=[" << lo << "," << hi << "] not in (1,2)\n";
    return false;
  }
  interval Lam(lo, hi);
  interval theta0 = interval::pi() / 2 * tau0;
  try {
    LocalData La = localLemma(d, P, interval(lo), theta0);
    LocalData Lb = localLemma(d, P, interval(hi), theta0);
    LocalData LB = localLemma(d, P, Lam, theta0);

    RunOut ra = rigorousRun(d, P, interval(lo), La, tau0);
    if (!(ra.gEnd.rightBound() < 0)) {
      if (VERBOSE) std::cerr << "    pad=" << pad << " tau0=" << tau0
                             << " : g(lo)=" << ra.gEnd << " not < 0\n";
      return false;
    }
    RunOut rb = rigorousRun(d, P, interval(hi), Lb, tau0);
    if (!(rb.gEnd.leftBound() > 0)) {
      if (VERBOSE) std::cerr << "    pad=" << pad << " tau0=" << tau0
                             << " : g(hi)=" << rb.gEnd << " not > 0\n";
      return false;
    }
    RunOut r = rigorousRun(d, P, Lam, LB, tau0);
    if (!(r.fEnd.rightBound() < 0)) {
      if (VERBOSE) std::cerr << "    pad=" << pad << " : f(pi/2) not < 0\n";
      return false;
    }

    interval a = interval(d - 1) / interval(d);
    interval maxv(a.leftBound(), std::max(a.rightBound(), r.fMax));
    interval minv(std::min(LB.fLoc.leftBound(), r.fMin), r.fEnd.rightBound());
    interval mean = 2 * sigmaD(d) * r.JEnd;      // average over S^{d-1}
    interval C = (P - 2) / 2 * (maxv + minv)
               + interval(dt) * (dt + 2) / (Lam + dt) * Hlift(d, dt, Lam) * mean;
    if (!(C.leftBound() > 0 || C.rightBound() < 0)) {
      if (VERBOSE) std::cerr << "    pad=" << pad << " tau0=" << tau0
                             << " : C=" << C << " contains 0\n";
      return false;
    }

    out.lamLo = lo; out.lamHi = hi; out.tau0 = tau0;
    out.cLo = C.leftBound(); out.cHi = C.rightBound();
    return true;
  } catch (std::exception& e) {
    if (VERBOSE) std::cerr << "    pad=" << pad << " tau0=" << tau0
                           << " : " << e.what() << "\n";
    return false;
  }
}

// Ladder of (pad, tau0).  pad must be large enough for the uniform sign change
// but small enough that the enclosure of C excludes 0, so we sweep both ways.
static bool tryBox(int d, int dt, double p1, double p2, Cert& out) {
  double l1 = seedLambda(d, p1), l2 = seedLambda(d, p2);
  // The pad must be small enough that Lambda stays inside (1,2) and that the
  // enclosure of C excludes 0, and large enough that the sign change in (C2)
  // is resolved.  Near p = 2 the eigenvalue tends to 2 and the defect tends
  // to 0, so both upper constraints tighten and very small pads are needed.
  static const double pads[] = {3e-4, 1e-3, 1e-4, 3e-5, 3e-3, 1e-5, 1e-6,
                                1e-7, 1e-8, 3e-9, 1e-9, 1e-10, 1e-11};
  static const double taus[] = {1.0/1024, 1.0/4096, 1.0/16384};
  for (double tau0 : taus)
    for (double pad : pads)
      if (attempt(d, dt, p1, p2, l1, l2, pad, tau0, out)) return true;
  return false;
}

int main(int argc, char** argv) {
  if (argc < 4) {
    std::cerr << "usage: " << argv[0]
              << " d a b [minwidth] [startwidth] [dtarget]\n";
    return 2;
  }
  int d = std::atoi(argv[1]);
  double a = std::atof(argv[2]), b = std::atof(argv[3]);
  double minw = (argc > 4) ? std::atof(argv[4]) : 1e-4;
  double startw = (argc > 5) ? std::atof(argv[5]) : (b - a);
  int dt = (argc > 6) ? std::atoi(argv[6]) : d;
  if (dt < d) { std::cerr << "dtarget must be >= d\n"; return 2; }

  std::cout << std::setprecision(17);
  // depth-first, left to right, so output is ordered and tiles [a,b]
  std::vector<std::pair<double,double>> stack;
  {
    long n = (long)std::ceil((b - a) / startw - 1e-12);
    if (n < 1) n = 1;
    for (long i = n - 1; i >= 0; --i) {          // push right to left
      // The first and last endpoints are pinned to a and b: a + (b-a)*n/n
      // need not equal b in floating point, and a sub-ulp mismatch there
      // would leave a nonempty set of exponents uncertified.
      double x1 = (i == 0)     ? a : a + (b - a) * i / n;
      double x2 = (i == n - 1) ? b : a + (b - a) * (i + 1) / n;
      stack.push_back({x1, x2});
    }
  }
  long nOK = 0, nFail = 0;
  while (!stack.empty()) {
    auto [p1, p2] = stack.back(); stack.pop_back();
    Cert c;
    if (tryBox(d, dt, p1, p2, c)) {
      // one certificate line, in the format read by ../verify
      std::cout << p1 << " " << p2 << " " << c.lamLo << " " << c.lamHi
                << " " << c.tau0 << "\n" << std::flush;
      ++nOK;
    } else if (p2 - p1 <= minw) {
      std::cout << "# FAIL " << p1 << " " << p2 << "\n" << std::flush;
      ++nFail;
    } else {
      double mid = 0.5 * (p1 + p2);
      stack.push_back({mid, p2});     // push right first: pop left first
      stack.push_back({p1, mid});
    }
  }
  std::cerr << "# certified " << nOK << " intervals, " << nFail << " failures\n";
  return nFail ? 1 : 0;
}
