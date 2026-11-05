"""Verification of the symbol balls at the retained nodes (condition (ii) of Theorem 5.10).

For each node from the ball file `balls/c<c>_<sector>.csv.gz`, the symbol `mu(r_k)` is recomputed from the
definition by means of Arb, and it is checked that the Arb ball lies in the published ball `[m - rho, m +
rho]`. Then the lower bounds for the weights, which the matrix verification takes from the published balls,
are valid.
"""
import sys
from multiprocessing import Pool

from flint import arb

from common import Certificate, Symbol, ball

_S = {}


def _work(task):
    c, prec, items = task
    if _S.get("key") != (c, prec):
        _S["key"] = (c, prec)
        _S["sym"] = Symbol(c, prec)
    sym = _S["sym"]
    bad = []
    for node, mid_s, rad_s in items:
        A = sym.mu(node)
        m, rho = ball(mid_s, rad_s)
        if not (A - m).abs_upper() <= rho:
            bad.append(node)
    return len(items), bad


def run(cert, procs=4):
    odd, even = cert.balls()
    nodes = [(k, m, r) for k, m, r in odd] + [(2 * j, m, r) for j, m, r in even]
    prec = cert.bits + 128
    tasks = [(cert.c, prec, nodes[i:i + 100]) for i in range(0, len(nodes), 100)]
    total, bad = 0, []
    with Pool(procs) as pool:
        for n, b in pool.imap_unordered(_work, tasks):
            total += n
            bad += b
    ok = not bad
    print("symbol: retained nodes %d, Arb ball lies in the published ball: %d%s"
          % (total, total - len(bad), "" if ok else "; failed at the nodes " + str(sorted(bad)[:20])))
    return ok


if __name__ == "__main__":
    sys.exit(0 if run(Certificate(sys.argv[1], sys.argv[2])) else 1)
