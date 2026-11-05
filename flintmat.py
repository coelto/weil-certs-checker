"""Arb matrices via the FLINT library shipped together with python-flint.

python-flint does not provide access to the ball Cholesky decomposition `arb_mat_cho`, so here documented
FLINT functions are called directly. The ball `arb` is transferred to the matrix exactly: the midpoint and the
radius are written as an integer multiplied by a power of two.
"""
import ctypes
import glob
import os

import flint


def _find_library():
    if os.environ.get("FLINT_LIB"):
        return os.environ["FLINT_LIB"]
    base = os.path.dirname(flint.__file__)
    parent = os.path.dirname(base)
    patterns = [os.path.join(parent, "python_flint.libs", "libflint*"),
                os.path.join(base, ".dylibs", "libflint*"),
                os.path.join(base, "*.libs", "libflint*"),
                os.path.join(parent, "*.libs", "libflint*"),
                os.path.join(base, "libflint*")]
    for p in patterns:
        hits = sorted(glob.glob(p))
        if hits:
            if os.name == "nt":
                os.add_dll_directory(os.path.dirname(hits[0]))
            return hits[0]
    raise OSError("FLINT library from python-flint not found. Specify the path to it in the variable FLINT_LIB.")


L = ctypes.CDLL(_find_library())
slong = ctypes.c_longlong if os.name == "nt" else ctypes.c_long
vp = ctypes.c_void_p
for name, args, res in [
        ("flint_set_num_threads", [ctypes.c_int], None),
        ("arb_mat_init", [vp, slong, slong], None),
        ("arb_mat_clear", [vp], None),
        ("arb_mat_entry_ptr", [vp, slong, slong], vp),
        ("arb_mat_cho", [vp, vp, slong], ctypes.c_int),
        ("arb_mat_mul", [vp, vp, vp, slong], None),
        ("arb_mat_add", [vp, vp, vp, slong], None),
        ("arb_mat_sub", [vp, vp, vp, slong], None),
        ("arb_mat_transpose", [vp, vp], None),
        ("arb_mat_get_mid", [vp, vp], None),
        ("arb_mat_approx_solve", [vp, vp, vp, slong], ctypes.c_int),
        ("arb_load_str", [vp, ctypes.c_char_p], ctypes.c_int),
        ("arb_dump_str", [vp], vp),
        ("flint_free", [vp], None)]:
    f = getattr(L, name)
    f.argtypes = args
    f.restype = res


def _hx(v):
    v = int(v)
    return ("-%x" % -v) if v < 0 else ("%x" % v)


def dump(x):
    m, e = x.mid().man_exp()
    r, f = x.rad().man_exp()
    m, e, r, f = int(m), int(e), int(r), int(f)
    if m == 0:
        e = 0
    if r == 0:
        f = 0
    return "%s %s %s %s" % (_hx(m), _hx(e), _hx(r), _hx(f))


def load(s):
    a, b, c, d = s.split()
    return flint.arb(flint.arb(flint.fmpz(int(a, 16))) * flint.arb(2) ** int(b, 16),
                     flint.arb(flint.fmpz(int(c, 16))) * flint.arb(2) ** int(d, 16))


class Mat:
    def __init__(self, r, c):
        self.r, self.c = r, c
        self.buf = ctypes.create_string_buffer(64)
        L.arb_mat_init(self.buf, r, c)

    def ptr(self, i, j):
        return L.arb_mat_entry_ptr(self.buf, i, j)

    def set(self, i, j, x):
        if L.arb_load_str(self.ptr(i, j), dump(x).encode()) != 0:
            raise ValueError("Failed to write the ball to the matrix (arb_load_str).")

    def get(self, i, j):
        p = L.arb_dump_str(self.ptr(i, j))
        s = ctypes.string_at(p).decode()
        L.flint_free(p)
        return load(s)

    def free(self):
        if self.buf is not None:
            L.arb_mat_clear(self.buf)
            self.buf = None
