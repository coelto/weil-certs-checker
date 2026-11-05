"""Verification of the matrix inequality `Id - W_up >= 0` (condition (iii) of Theorem 5.10).

The matrix is assembled anew from the published balls of the symbol by Proposition 5.2, Lemma 5.3 and, in the
odd sector, Proposition 5.9. Nothing other than the balls is taken from the program that constructed the
certificate. The weights are set as follows.

    w_k = (mu(r_k) - t1)/L                      for retained odd k,
    d_j - eps = (mu(omega_j) + t1)/2 - eps      for retained even j,
    d_eff = (tau_e + t1)/2 - eps.

Positive definiteness is proved via the ball Cholesky decomposition `arb_mat_cho`: if all pivots are positive,
then every symmetric matrix whose entries lie in these balls is positive definite (Lemma 5.4).
"""
import sys
import time
from multiprocessing import Pool

from flint import arb, ctx, fmpz

import flintmat as fm
from common import Certificate, ball

_W = {}


def _init(cs, xs, prec):
    ctx.prec = prec
    _W["c"] = [fm.load(s) for s in cs]
    _W["x"] = xs


def _sums(task):
    out = []
    for a, y in task:
        yb = arb(fmpz(y)) if isinstance(y, int) else fm.load(y)
        f = arb(0)
        f2 = arb(0)
        for cj, xj in zip(_W["c"], _W["x"]):
            inv = 1 / (yb - xj)
            t = cj * inv
            f += t
            f2 += t * inv
        out.append((a, fm.dump(f), fm.dump(f2)))
    return out


def assemble(cert, prec, procs):
    ctx.prec = prec
    odd_sector = cert.odd
    Lb = arb(cert.c).log()
    pi = arb.pi()
    eps = arb(cert.eps_s)
    t1 = arb(cert.t1_s)
    d_eff = (arb(cert.tau_e) + t1) / 2 - eps
    assert d_eff > 0
    inv_deff = 1 / d_eff
    odd, even = cert.balls()
    kp, kn, good, bad = [], [], [], []
    for k, ms, rs in odd:
        m, r = ball(ms, rs)
        w = (arb(m, r) - t1) / Lb
        if w > 0:
            kp.append((k, w))
        elif w <= 0:
            kn.append((k, -w))
        else:
            wl = w.lower()
            (kp if wl > 0 else kn).append((k, arb(wl) if wl > 0 else arb(-wl)))
    for j, ms, rs in even:
        if odd_sector and j == 0:
            continue
        m, r = ball(ms, rs)
        dj = (arb(m, r) + t1) / 2 - eps
        if dj > 0:
            good.append((j, dj))
        elif dj <= 0:
            bad.append((j, d_eff - dj))
        else:
            dl = dj.lower()
            (good, bad)[0 if dl > 0 else 1].append((j, arb(dl) if dl > 0 else d_eff - arb(dl)))
    a2 = 2 * Lb / pi ** 2
    a2_0 = Lb / pi ** 2

    def asq(j):
        if odd_sector:
            return a2 * (4 * j) ** 2
        return a2_0 if j == 0 else a2

    cj = [(1 / dj - inv_deff) * asq(j) for j, dj in good]
    xs = [4 * j * j for j, _ in good]
    if odd_sector:
        cols = [("odd", k) for k, _ in kp] + [("odd", k) for k, _ in kn] + [("pole", None)]
        nk = len(kp)
        bv = (Lb / 4).sinh()
        vv = (Lb / 2).sinh() - Lb / 2
    else:
        cols = [("odd", k) for k, _ in kp] + [("pole", None)] + [("odd", k) for k, _ in kn]
        nk = len(kp) + 1
        bv = -(Lb / pi) * (Lb / 4).sinh()
        vv = (Lb / 2).sinh() + Lb / 2
    yv = -(Lb ** 2) / (4 * pi ** 2)
    ch4 = (Lb / 4).cosh()

    def sig(k):
        return 1 if (k // 2) % 2 == 0 else -1

    def ycol(a):
        kind, k = cols[a]
        return arb(k * k) if kind == "odd" else yv

    def bcol(a):
        kind, k = cols[a]
        if kind == "pole":
            return bv
        return arb(sig(k)) if odd_sector else arb(2 * k * sig(k))

    tasks = [(a, cols[a][1] ** 2 if cols[a][0] == "odd" else fm.dump(yv)) for a in range(len(cols))]
    chunks = [tasks[i:i + 16] for i in range(0, len(tasks), 16)]
    F, F2 = [None] * len(cols), [None] * len(cols)
    with Pool(procs, initializer=_init, initargs=([fm.dump(x) for x in cj], xs, prec)) as pool:
        for res in pool.imap_unordered(_sums, chunks):
            for a, f, f2 in res:
                F[a], F2[a] = fm.load(f), fm.load(f2)
    ctx.prec = prec

    def plain(a, b):
        ka, kb = cols[a], cols[b]
        if ka[0] == "odd" and kb[0] == "odd":
            return Lb / 2 if ka[1] == kb[1] else None
        if ka[0] == "odd" or kb[0] == "odd":
            k = ka[1] if ka[0] == "odd" else kb[1]
            r = pi * k / Lb
            q = r * r + arb(1) / 4
            return sig(k) * ch4 / q if odd_sector else 2 * sig(k) * r * ch4 / q
        return vv

    def ip(a, b):
        s = bcol(a) * bcol(b) * (F2[a] if a == b else (F[a] - F[b]) / (ycol(b) - ycol(a)))
        p = plain(a, b)
        return s if p is None else p * inv_deff + s

    def phi_ip(j, a):
        aj = asq(j).sqrt()
        if j % 2 == 1:
            aj = -aj
        return aj * bcol(a) / (ycol(a) - 4 * j * j)

    nbe = len(bad)
    nb = len(cols) - nk + nbe
    sq = [w.nonnegative_part().sqrt() for _, w in bad] + [w.nonnegative_part().sqrt() for _, w in kn]
    if odd_sector:
        sq.append(arb(2).sqrt())
    S = fm.Mat(nk, nk)
    for a in range(nk):
        for b in range(a + 1):
            v = ip(a, b)
            if a == b:
                v += 1 / (kp[a][1] if a < len(kp) else arb(2))
            S.set(a, b, v)
            S.set(b, a, v)
    Y = fm.Mat(nk, nb)
    for a in range(nk):
        for i in range(nb):
            v = phi_ip(bad[i][0], a) * inv_deff if i < nbe else ip(a, nk + i - nbe)
            Y.set(a, i, v * sq[i])
    X = fm.Mat(nb, nb)
    for i in range(nb):
        for j in range(i + 1):
            if i < nbe and j < nbe:
                v = inv_deff if i == j else arb(0)
            elif j < nbe:
                v = phi_ip(bad[j][0], nk + i - nbe) * inv_deff
            else:
                v = ip(nk + i - nbe, nk + j - nbe)
            v = v * sq[i] * sq[j]
            X.set(i, j, v)
            X.set(j, i, v)
    Sm, Ym, Z, Zm = fm.Mat(nk, nk), fm.Mat(nk, nb), fm.Mat(nk, nb), fm.Mat(nk, nb)
    fm.L.arb_mat_get_mid(Sm.buf, S.buf)
    fm.L.arb_mat_get_mid(Ym.buf, Y.buf)
    fm.L.arb_mat_approx_solve(Z.buf, Sm.buf, Ym.buf, min(prec, 512))
    fm.L.arb_mat_get_mid(Zm.buf, Z.buf)
    Sm.free(); Ym.free(); Z.free()
    Zt, T, SZ, Q, Tt = fm.Mat(nb, nk), fm.Mat(nb, nb), fm.Mat(nk, nb), fm.Mat(nb, nb), fm.Mat(nb, nb)
    fm.L.arb_mat_transpose(Zt.buf, Zm.buf)
    fm.L.arb_mat_mul(T.buf, Zt.buf, Y.buf, prec)
    fm.L.arb_mat_mul(SZ.buf, S.buf, Zm.buf, prec)
    fm.L.arb_mat_mul(Q.buf, Zt.buf, SZ.buf, prec)
    fm.L.arb_mat_transpose(Tt.buf, T.buf)
    W = fm.Mat(nb, nb)
    fm.L.arb_mat_sub(W.buf, X.buf, T.buf, prec)
    fm.L.arb_mat_sub(W.buf, W.buf, Tt.buf, prec)
    fm.L.arb_mat_add(W.buf, W.buf, Q.buf, prec)
    for m in (S, Y, X, Zm, Zt, T, SZ, Q, Tt):
        m.free()
    A = fm.Mat(nb, nb)
    for i in range(nb):
        for j in range(i + 1):
            v = (1 if i == j else 0) - W.get(i, j)
            A.set(i, j, v)
            A.set(j, i, v)
    W.free()
    info = "|K+| = %d, |K-| = %d, |G| = %d, bad even = %d, order Id - W_up = %d" % (len(kp), len(kn), len(good),
                                                                                      nbe, nb)
    return A, nb, info


def run(cert, procs=4, threads=4):
    prec = cert.bits + 128
    fm.L.flint_set_num_threads(threads)
    t = time.time()
    A, n, info = assemble(cert, prec, procs)
    print("matrix: %s (%.0f s)" % (info, time.time() - t))
    Lc = fm.Mat(n, n)
    ok = bool(fm.L.arb_mat_cho(Lc.buf, A.buf, prec))
    Lc.free(); A.free()
    print("matrix: ball Cholesky decomposition for Id - W_up: %s" % ("done, Id - W_up > 0" if ok else "FAILED"))
    return ok


if __name__ == "__main__":
    sys.exit(0 if run(Certificate(sys.argv[1], sys.argv[2])) else 1)
