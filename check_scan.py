"""Verification of the threshold outside the retained nodes (condition (i) of Theorem 5.10, Lemmas 4.2 and 5.8).

The following is verified.

    1) the retained nodes contain all odd k <= K2 and all j <= J (in the odd sector 1 <= j <= J);
    2) mu(r_k) >= tau for all odd k in (K2, k_scan] outside the retained ones;
    3) mu(r_{2j}) >= tau_e for all even 2j in (2J, k_scan^e] outside the retained ones;
    4) monotone lower bound: mu_low(r_{k_scan + 1}) >= tau and mu_low(r_{k_scan^e + 1}) >= tau_e.

Items 2 and 3 for large `c` require lengthy computation. The flag `--sample N` verifies `N` random nodes of
each kind, and the flags `--odd-range A B` and `--even-range A B` verify an interval of indices. The result is
then marked as partial.
"""
import argparse
import random
import sys
from multiprocessing import Pool

from flint import arb

from common import Certificate, Symbol

_S = {}


def _sym(c, prec):
    if _S.get(prec, (None,))[0] != c:
        _S[prec] = (c, Symbol(c, prec))
    return _S[prec][1]


def _work(task):
    c, prec_lo, prec_hi, thr, nodes = task
    lo, hi = _sym(c, prec_lo), _sym(c, prec_hi)
    t = arb(thr)
    bad = []
    for k in nodes:
        if not lo.mu(k) >= t and not hi.mu(k) >= t:
            bad.append(k)
    return len(nodes), bad


def scan(pool, c, prec_hi, thr, nodes, label):
    tasks = [(c, 96, prec_hi, thr, nodes[i:i + 2000]) for i in range(0, len(nodes), 2000)]
    done, bad = 0, []
    for n, b in pool.imap_unordered(_work, tasks):
        done += n
        bad += b
        if done % 200000 < 2000:
            print("  " + "%s: nodes %d out of %d" % (label, done, len(nodes)), flush=True)
    print("%s: nodes checked %d, below the threshold: %d%s"
          % (label, done, len(bad), "" if not bad else " " + str(sorted(bad)[:20])))
    return not bad


def run(cert, procs=4, sample=0, odd_range=None, even_range=None, seed=1):
    odd, even = cert.balls()
    K = {k for k, _, _ in odd}
    E = {j for j, _, _ in even}
    ok = True
    j0 = 1 if cert.odd else 0
    miss_k = [k for k in range(1, cert.K2 + 1, 2) if k not in K]
    miss_j = [j for j in range(j0, cert.J + 1) if j not in E]
    print("retained nodes: |K| = %d, |E| = %d; odd k <= K2 = %d, missing: %d; j <= J = %d, missing: %d"
          % (len(K), len(E), cert.K2, len(miss_k), cert.J, len(miss_j)))
    ok &= not miss_k and not miss_j
    t1 = float(cert.t1_s)
    ok &= t1 <= cert.tau
    sym = Symbol(cert.c, cert.bits + 128)
    m1 = sym.mu_low(cert.k_scan + 1)
    m2 = sym.mu_low(cert.kk_scan + 1)
    mono = bool(m1 >= arb(cert.tau)) and bool(m2 >= arb(cert.tau_e))
    print("monotone bound: mu_low(r_%d) = %s (>= %g), mu_low(r_%d) = %s (>= %g): %s"
          % (cert.k_scan + 1, m1.str(8), cert.tau, cert.kk_scan + 1, m2.str(8), cert.tau_e, "yes" if mono else "NO"))
    ok &= mono
    a_odd, b_odd = odd_range or (cert.K2 + 1, cert.k_scan)
    a_even, b_even = even_range or (2 * cert.J + 1, cert.kk_scan)
    o0 = max(a_odd, cert.K2 + 1)
    odd_r = range(o0 + (1 - o0 % 2), b_odd + 1, 2)
    even_r = range(max((a_even + 1) // 2, cert.J + 1), b_even // 2 + 1)
    full = not sample and odd_range is None and even_range is None
    if sample:
        rnd = random.Random(seed)
        odd_nodes = sorted([k for k in rnd.sample(odd_r, min(sample + len(K), len(odd_r))) if k not in K][:sample])
        even_nodes = sorted([2 * j for j in rnd.sample(even_r, min(sample + len(E), len(even_r))) if j not in E][:sample])
    else:
        odd_nodes = [k for k in odd_r if k not in K]
        even_nodes = [2 * j for j in even_r if j not in E]
    with Pool(procs) as pool:
        ok &= scan(pool, cert.c, cert.bits + 128, cert.tau, odd_nodes, "scan of odd nodes")
        ok &= scan(pool, cert.c, cert.bits + 128, cert.tau_e, even_nodes, "scan of even nodes")
    print("scan: %s" % ("complete" if full else "PARTIAL (sample or interval), condition (i) not fully verified"))
    return ok, full


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("data", help="Unpacked folder of the certificate data (it contains certificates/ and balls/).")
    ap.add_argument("certificate", help="Name of the certificate folder, e.g. c10_even_eps.")
    ap.add_argument("--procs", type=int, default=4, help="Number of worker processes.")
    ap.add_argument("--sample", type=int, default=0, help="Check N random nodes of each kind instead of an exhaustive pass.")
    ap.add_argument("--seed", type=int, default=1, help="Seed of the random number generator for the sampling.")
    ap.add_argument("--odd-range", type=int, nargs=2, help="Check odd nodes with indices from A to B.")
    ap.add_argument("--even-range", type=int, nargs=2, help="Check even nodes with indices from A to B.")
    a = ap.parse_args()
    ok, _ = run(Certificate(a.data, a.certificate), a.procs, a.sample, a.odd_range, a.even_range, a.seed)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
