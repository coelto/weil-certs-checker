"""Negative controls of the verification program. Each control introduces one defect into a temporary copy of the
certificate data or code and runs the corresponding check in a separate process. The check must fail. If the
control passes, the defect has been missed.

    python controls.py <data> <certificate>

First the unmodified certificate must pass all three checks, otherwise the controls are meaningless.
"""
import gzip
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = HERE
data, name = sys.argv[1], sys.argv[2]
procs = sys.argv[sys.argv.index("--procs") + 1] if "--procs" in sys.argv else "8"
man = json.load(open(os.path.join(data, "certificates", name, "manifest.json"), encoding="utf-8"))
c, sector = None, man["sector"]
import re  # noqa: E402
c = int(re.search(r"c_num: ([0-9]+)", man["params"]).group(1))
BALLS = "c%d_%s.csv.gz" % (c, sector)
results = []


def copy_data(tmp, edit_manifest=None, edit_balls=None):
    """Folders certificates/<name>/ and balls/<file> from DATA in the temporary folder, modified as necessary."""
    d = os.path.join(tmp, "data")
    shutil.copytree(os.path.join(data, "certificates", name), os.path.join(d, "certificates", name))
    os.makedirs(os.path.join(d, "balls"))
    src = os.path.join(data, "balls", BALLS)
    if edit_balls:
        with gzip.open(src, "rt") as f:
            lines = f.readlines()
        with gzip.open(os.path.join(d, "balls", BALLS), "wt") as f:
            f.writelines(edit_balls(lines))
    else:
        shutil.copy(src, os.path.join(d, "balls", BALLS))
    if edit_manifest:
        p = os.path.join(d, "certificates", name, "manifest.json")
        m = json.load(open(p, encoding="utf-8"))
        edit_manifest(m)
        json.dump(m, open(p, "w", encoding="utf-8"))
    return d


def copy_pkg(tmp, old=None, new=None, module=None):
    p = os.path.join(tmp, "package")
    shutil.copytree(PKG, p, ignore=shutil.ignore_patterns("__pycache__", ".git", "*.pdf"))
    if module:
        f = os.path.join(p, module)
        s = open(f, encoding="utf-8").read()
        assert s.count(old) == 1, (module, old)
        open(f, "w", encoding="utf-8").write(s.replace(old, new))
    return p


def run(label, script, pkg, d, extra=()):
    t = time.time()
    r = subprocess.run([sys.executable, os.path.join(pkg, script), d, name, *extra], capture_output=True, text=True,
                       encoding="utf-8", errors="replace", env=dict(os.environ, PYTHONIOENCODING="utf-8"))
    caught = r.returncode != 0
    results.append((label, caught))
    tail = (r.stdout.strip().splitlines() or [""])[-1][:110]
    print("CONTROL %-50s %s (%.0f s) | %s" % (label, "caught" if caught else "MISSED", time.time() - t, tail), flush=True)


def shift_ball(lines):
    """The midpoint of one retained odd ball is shifted by three radii."""
    head, rest = lines[0], lines[1:]
    i = next(j for j, l in enumerate(rest) if l.startswith("odd,"))
    parity, idx, mid, rad = rest[i].strip().split(",")[:4]
    neg = mid.startswith("-")
    body, e = (mid[1:] if neg else mid)[2:].split("p")
    rb, re_ = rad[2:].split("p")
    def val(b):  # noqa: E306
        a, _, f = b.partition(".")
        return int(a + f, 16), 4 * len(f)
    m, fm = val(body)
    r, fr = val(rb)
    # The sum m*2^(e - fm) + 3r*2^(re - fr) written with a common exponent.
    E = min(int(e) - fm, int(re_) - fr)
    M = (-m if neg else m) * 2 ** (int(e) - fm - E) + 3 * r * 2 ** (int(re_) - fr - E)
    newmid = ("-" if M < 0 else "") + "0x%xp%d" % (abs(M), E)
    rest[i] = ",".join([parity, idx, newmid] + rest[i].strip().split(",")[3:]) + "\n"
    return [head] + rest


def params_edit(key, f):
    def g(m):
        m["params"] = re.sub(key + r": ([0-9.eE+-]+)", lambda x: "%s: %s" % (key, f(x.group(1))), m["params"])
    return g


with tempfile.TemporaryDirectory() as tmp0:
    # First the unmodified certificate must pass all three checks, otherwise the controls are meaningless.
    pkg = copy_pkg(os.path.join(tmp0, "p0"))
    d = copy_data(os.path.join(tmp0, "d0"))
    for s in ("check_symbol.py", "check_scan.py", "check_matrix.py"):
        r = subprocess.run([sys.executable, os.path.join(pkg, s), d, name], capture_output=True, text=True,
                           encoding="utf-8", errors="replace", env=dict(os.environ, PYTHONIOENCODING="utf-8"))
        print("UNCHANGED %-16s %s" % (s, "passed" if r.returncode == 0 else "NOT PASSED (the controls are meaningless)"), flush=True)
        if r.returncode != 0:
            sys.exit(2)
    # Corrupting the code only matters when it makes the check reject a valid certificate.
    controls = [
        ("symbol: one published ball shifted by three radii", "check_symbol.py", {}, {"edit_balls": shift_ball}),
        ("matrix: eps doubled", "check_matrix.py", {},
         {"edit_manifest": lambda m: m.__setitem__("eps", repr(2 * float(m["eps"]) if float(m["eps"]) > 0 else 1e-3))}),
        ("scan: odd scan end divided by 4", "check_scan.py", {},
         {"edit_manifest": lambda m: m.__setitem__("k_scan_odd", int(m["k_scan_odd"]) // 4)}),
        ("scan: the threshold tau_o increased by 1", "check_scan.py", {},
         {"edit_manifest": params_edit("tau_o", lambda v: repr(float(v) + 1))}),
        ("scan: the initial segment K2 increased by 40", "check_scan.py", {},
         {"edit_manifest": params_edit("k2", lambda v: str(int(v) + 40))}),
        ("matrix: the sign of the pole constant b_v changed", "check_matrix.py",
         {"module": "check_matrix.py", "old": "bv = -(Lb / pi) * (Lb / 4).sinh()", "new": "bv = (Lb / pi) * (Lb / 4).sinh()"}
         if sector == "even" else
         {"module": "check_matrix.py", "old": "bv = (Lb / 4).sinh()", "new": "bv = -(Lb / 4).sinh()"}, {}),
        ("symbol: the weights of the atoms doubled", "check_symbol.py",
         {"module": "common.py", "old": "return h - T - 2 * ps", "new": "return h - T - 4 * ps"}, {}),
    ]
    for i, (label, script, code, dat) in enumerate(controls):
        tmp = os.path.join(tmp0, "c%d" % i)
        os.makedirs(tmp)
        pkg = copy_pkg(tmp, **code) if code else copy_pkg(tmp)
        d = copy_data(tmp, **dat)
        run(label, script, pkg, d)

missed = [l for l, ok in results if not ok]
print("controls caught: %d out of %d%s" % (len(results) - len(missed), len(results), "" if not missed else "; MISSED " + str(missed)))
sys.exit(1 if missed else 0)
