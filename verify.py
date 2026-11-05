"""Full verification of one certificate: conditions (i), (ii) and (iii) of Theorem 5.10 of the paper.

    python verify.py <data> <certificate> [--procs N] [--threads N] [--sample N | --odd-range A B --even-range A B]

Here `<data>` is the unpacked certificate data folder (containing `certificates/` and `balls/`), and
`<certificate>` is the name of the certificate folder, for example `c10_even_eps`.

The outcome `VERIFIED` means that all three conditions have been verified in full, and then `Q_c(f) >= eps
||f||^2` in the sector of the certificate. The outcome `VERIFIED (partial scan)` means that the symbol and the
matrix have been verified in full, while the scan of the nodes has been carried out only on a sample or an
interval.
"""
import argparse
import sys
import time

import check_matrix
import check_scan
import check_symbol
from common import Certificate


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("data", help="Unpacked folder of the certificate data (it contains certificates/ and balls/).")
    ap.add_argument("certificate", help="Name of the certificate folder, e.g. c10_even_eps.")
    ap.add_argument("--procs", type=int, default=4, help="Number of worker processes.")
    ap.add_argument("--threads", type=int, default=4, help="Number of FLINT threads for matrix operations.")
    ap.add_argument("--sample", type=int, default=0, help="Check N random nodes of each kind instead of an exhaustive pass.")
    ap.add_argument("--seed", type=int, default=1, help="Seed of the random number generator for the sampling.")
    ap.add_argument("--odd-range", type=int, nargs=2, help="Check odd nodes with indices from A to B.")
    ap.add_argument("--even-range", type=int, nargs=2, help="Check even nodes with indices from A to B.")
    a = ap.parse_args()
    cert = Certificate(a.data, a.certificate)
    print("certificate %s: c = %d, sector %s, eps = %s, bits %d, K2 = %d, J = %d, scan endpoints %d and %d"
          % (a.certificate, cert.c, cert.sector, cert.eps_s, cert.bits, cert.K2, cert.J, cert.k_scan, cert.kk_scan))
    t = time.time()
    ok_symbol = check_symbol.run(cert, a.procs)
    ok_scan, full = check_scan.run(cert, a.procs, a.sample, a.odd_range, a.even_range, a.seed)
    ok_matrix = check_matrix.run(cert, a.procs, a.threads)
    ok = ok_symbol and ok_scan and ok_matrix
    verdict = ("VERIFIED" if full else "VERIFIED (partial scan)") if ok else "NOT VERIFIED"
    print("%s: %s (%.0f s)" % (a.certificate, verdict, time.time() - t))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
