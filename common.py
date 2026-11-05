"""Common part: reading the certificate from the published data and the symbol of the truncated form in Arb ball
arithmetic.

The symbol at the node `r` (Section 3 of the paper and Lemma 4.2) is given by the formulas below.

    mu(r) = Re psi(1/4 + i r/2) - ln pi - T(r) - 2 sum_{n<c} Lambda(n) n^{-1/2} cos(r ln n),
    T(r)  = -2 sum_{n>=0} c^{-a_n} (a_n cos(rL) - r sin(rL)) / (a_n^2 + r^2),   a_n = 2n + 1/2,   L = ln c.

At the nodes `r_k = pi k / L` the relations `cos(r_k L) = (-1)^k` and `sin(r_k L) = 0` hold exactly. The tail
of the series `T` after `N` terms does not exceed `2 c^{-a_N} / (a_N (1 - c^{-2}))` and is added to the
radius.

The monotone lower bound (Lemma 4.2) is as follows.

    mu_low(r) = Re psi(1/4 + i r/2) - ln pi - 2 sigma_c - Tbar(r),   Tbar(r) = 2 sum_{n>=0} c^{-a_n} (1/r + a_n/r^2),
    2 sigma_c = 2 sum_{n<c} Lambda(n)/sqrt(n).
"""
import gzip
import json
import math
import os
import re

from flint import acb, arb, ctx, fmpz


def prime_powers(c):
    out = []
    for n in range(2, c):
        p = next(q for q in range(2, n + 1) if n % q == 0)
        m = n
        while m % p == 0:
            m //= p
        if m == 1:
            out.append((n, p))
    return out


class Symbol:
    def __init__(self, c, prec):
        ctx.prec = prec
        self.c, self.prec = c, prec
        self.L = arb(c).log()
        self.pi = arb.pi()
        self.lnpi = self.pi.log()
        self.terms = [(arb(n).log(), arb(p).log() / arb(n).sqrt()) for n, p in prime_powers(c)]
        self.sigma2 = 2 * sum((w for _, w in self.terms), arb(0))
        self.tail = []
        n = 0
        while (2 * n + 0.5) * math.log2(c) <= prec + 40:
            a = arb(2 * n) + arb(1) / 2
            self.tail.append((a, (-a * self.L).exp()))
            n += 1
        self.aN = arb(2 * n) + arb(1) / 2
        self.caN = (-self.aN * self.L).exp()
        self.rem = 2 * self.caN / (self.aN * (1 - (-2 * self.L).exp()))

    def r(self, k):
        return self.pi * k / self.L

    def mu(self, k):
        ctx.prec = self.prec
        r = self.r(k)
        h = acb(arb(1) / 4, r / 2).digamma().real - self.lnpi
        sgn = 1 if k % 2 == 0 else -1
        r2 = r * r
        T = arb(0)
        for a, ca in self.tail:
            T += ca * a / (a * a + r2)
        T = -2 * sgn * T + arb(0, self.rem.upper())
        ps = arb(0)
        for lnn, w in self.terms:
            ps += w * (r * lnn).cos()
        return h - T - 2 * ps

    def mu_low(self, k):
        ctx.prec = self.prec
        r = self.r(k)
        h = acb(arb(1) / 4, r / 2).digamma().real - self.lnpi
        tb = arb(0)
        for a, ca in self.tail:
            tb += ca * (1 / r + a / (r * r))
        rem = 2 * self.caN * (1 / r + self.aN / (r * r)) / (1 - (-2 * self.L).exp()) * 2
        tb = 2 * tb + arb(0, rem.upper())
        return h - self.sigma2 - tb


def parse_hex(s):
    neg = s.startswith("-")
    s = s[1:] if neg else s
    body, e = s[2:].split("p")
    e = int(e)
    if "." in body:
        a, b = body.split(".")
        m = int(a + b, 16)
        e -= 4 * len(b)
    else:
        m = int(body, 16)
    return (-m if neg else m), e


def exact(m, e):
    return arb(fmpz(m)) * arb(2) ** e


def ball(mid_s, rad_s):
    m, e = parse_hex(mid_s)
    rm, re_ = parse_hex(rad_s)
    return exact(m, e), exact(rm, re_)


class Certificate:
    def __init__(self, data, name):
        self.name = name
        self.dir = os.path.join(data, "certificates", name)
        self.man = json.load(open(os.path.join(self.dir, "manifest.json"), encoding="utf-8"))
        p = self.man["params"]
        num = lambda key: re.search(key + r": ([0-9.eE+-]+)", p).group(1)
        self.c = int(num("c_num"))
        assert int(num("c_den")) == 1
        self.bits = int(num("bits"))
        self.K2 = int(num("k2"))
        self.J = int(num("j_exact"))
        self.tau = float(num("tau_o"))
        self.tau_e = float(num("tau_e"))
        self.sector = self.man["sector"]
        self.odd = self.sector == "odd"
        self.eps_s = self.man["eps"]
        self.t1_s = self.man["t1"]
        self.k_scan = int(self.man["k_scan_odd"])
        self.kk_scan = int(self.man["kk_scan_even"])
        self.balls_path = os.path.join(data, "balls", "c%d_%s.csv.gz" % (self.c, self.sector))

    def balls(self):
        odd, even = [], []
        with gzip.open(self.balls_path, "rt") as f:
            head = f.readline().strip().split(",")
            assert head[:4] == ["parity", "index", "mid_exact", "rad_exact"], head
            for line in f:
                parity, idx, mid_s, rad_s = line.strip().split(",")[:4]
                (odd if parity == "odd" else even).append((int(idx), mid_s, rad_s))
        return odd, even
