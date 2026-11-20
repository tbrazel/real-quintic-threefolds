#!/usr/bin/env python3
"""
verify_lines.py -- check a certificate for the lines on a quintic threefold X = {F = 0} in P^4.

    pip install python-flint
    python verify_lines.py quintic_2875.txt lines_2875.txt
    python verify_lines.py quintic_15.txt   lines_15.txt

INPUT
  quintic file : one monomial per row, "e0 e1 e2 e3 e4 num k", meaning the term (num / 2^k) x0^e0 x1^e1 ... x4^e4.
                 The coefficients are read exactly.
  lines file   : one box per row, either
                     i j r R  y1 ... y6                          (a real box), or
                     i j r C  Re y1  Im y1  ...  Re y6  Im y6    (a complex box).
                 The line with chart coordinates y is the column span of the 5x2 matrix Y whose rows i and j
                 (0-based) are (1,0) and (0,1), and whose other three rows, in increasing order, are
                 (y1,y2), (y3,y4), (y5,y6). The box is y +- r in every real (and imaginary) part.

WHAT IS CHECKED  (in Arb ball arithmetic, so every inequality below is proved, not estimated)
  Let Phi_k(y) = F(Y(y) (1, tau_k)^T) for tau_k in {-3/2, -3/4, -1/4, 1/4, 3/4, 3/2}. A binary quintic vanishes iff
  it vanishes at 6 points, so Phi(y) = 0 iff the line lies on X; and D Phi(y) is invertible iff the line is a
  reduced isolated point of the Fano scheme of X (the 6 values are an invertible image of the 6 coefficients of F
  restricted to the line).
  1. Krawczyk test on every box: K(box) is contained in the interior of the box. Then Phi has exactly one zero in
     the box, and D Phi is invertible there.
  2. Reality. An R box is a real box, and the test is run in real arithmetic, so its zero is real: a real line.
     For a C box, some coordinate has imaginary part bounded away from 0 on the whole box, so its zero is not
     real. The chart is defined over R, so a line is real iff its chart coordinates are real: a non-real line.
  3. Distinctness. The Plücker vectors of the boxes are pairwise non-proportional, so no line is counted twice.
  4. Negative control: the same test must fail on boxes moved off the zero.

If 2875 distinct reduced lines are certified, they are ALL the lines of X: the incidence correspondence
{(F, l) : l in X} -> P^125 is proper of degree 2875 and etale at each certified line, so there is no room for any
other line (see the paper). The script prints the number of real ones.
"""
import argparse
import sys
import time
from flint import arb, acb, arb_mat, acb_mat, ctx

TAUS_Q = [(-3, 2), (-3, 4), (-1, 4), (1, 4), (3, 4), (3, 2)]   # exact rationals p/q


# ------------------------------------------------------------------------------------------------ input --
def read_quintic(path):
    """Monomials and their exact dyadic coefficients."""
    mons = []
    for row in open(path):
        row = row.split("#")[0].split()
        if not row:
            continue
        e = tuple(int(v) for v in row[:5])
        num, k = int(row[5]), int(row[6])
        assert sum(e) == 5 and min(e) >= 0, f"not a quintic monomial: {e}"
        c = arb(num) / arb(2) ** k               # exact: a power of 2 only shifts the exponent
        assert c.is_exact(), "coefficient not represented exactly; raise --prec"
        mons.append((e, c))
    return mons


def read_lines(path):
    boxes = []
    for row in open(path):
        row = row.split("#")[0].split()
        if not row:
            continue
        i, j, r, kind = int(row[0]), int(row[1]), row[2], row[3]
        vals = row[4:]
        if kind == "R":
            assert len(vals) == 6
            y = [(v, "0") for v in vals]
        elif kind == "C":
            assert len(vals) == 12
            y = [(vals[2 * a], vals[2 * a + 1]) for a in range(6)]
        else:
            raise ValueError(f"unknown box type {kind!r}")
        boxes.append(dict(i=i, j=j, r=r, real=(kind == "R"), y=y))
    return boxes


# ---------------------------------------------------------------------------------- the line system --
class LineSystem:
    """F and its gradient, and the 6x6 chart system Phi with its Jacobian."""

    def __init__(self, mons):
        self.mons = mons
        # gradient: dF/dx_v = sum over monomials of c * e_v * x^(e - unit_v)
        self.grad = [[] for _ in range(5)]
        for e, c in mons:
            for v in range(5):
                if e[v] > 0:
                    f = list(e); f[v] -= 1
                    self.grad[v].append((tuple(f), c * e[v]))
        self.taus = [arb(p) / q for p, q in TAUS_Q]

    @staticmethod
    def _monomial(pw, e):
        val = None
        for v in range(5):
            if e[v]:
                val = pw[v][e[v]] if val is None else val * pw[v][e[v]]
        return val

    def F_grad(self, x, need_F):
        one = x[0] * 0 + 1
        pw = []
        for v in range(5):
            p = [one, x[v]]
            for _ in range(4):
                p.append(p[-1] * x[v])
            pw.append(p)
        F = sum(c * self._monomial(pw, e) for e, c in self.mons) if need_F else None
        G = [sum(c * self._monomial(pw, f) for f, c in self.grad[v]) for v in range(5)]
        return F, G

    def phi_jac(self, i, j, y, need_F=True):
        """Phi_k(y) = F(Y (1, tau_k)) and its Jacobian in the six chart coordinates y."""
        K = [v for v in range(5) if v not in (i, j)]
        Phi, J = [], []
        for tau in self.taus:
            x = [None] * 5
            x[i] = y[0] * 0 + 1
            x[j] = y[0] * 0 + tau
            for m in range(3):
                x[K[m]] = y[2 * m] + tau * y[2 * m + 1]
            F, G = self.F_grad(x, need_F)
            Phi.append(F)
            row = []
            for m in range(3):
                row += [G[K[m]], tau * G[K[m]]]
            J.append(row)
        return Phi, J


# ------------------------------------------------------------------------------------ Krawczyk test --
def make_box(b, shift=0):
    """Exact centre m (balls of radius 0) and the box X = m +- r, in arb (real box) or acb (complex box)."""
    r = arb(b["r"])
    rad = r.mid()           # the box is exactly m +- rad as stored by Arb (rad rounded up, never down)
    m, X = [], []
    for a, (re, im) in enumerate(b["y"]):
        mre = arb(re).mid() + (shift * r.mid() if a == 0 else 0)
        mim = arb(im).mid()
        if b["real"]:
            m.append(arb(mre)); X.append(arb(mre, rad))
        else:
            m.append(acb(mre, mim)); X.append(acb(arb(mre, rad), arb(mim, rad)))
    return m, X


def krawczyk(sys_, b, shift=0):
    """True iff K(X) lies in the interior of X:  K = m - C Phi(m) + (I - C DPhi(X)) (X - m)."""
    m, X = make_box(b, shift)
    Phi_m, J_m = sys_.phi_jac(b["i"], b["j"], m)
    _, J_X = sys_.phi_jac(b["i"], b["j"], X, need_F=False)
    Mat = arb_mat if b["real"] else acb_mat
    try:   # C: any approximate inverse of DPhi(m); only the success of the test depends on it
        Cb = Mat([[J_m[a][c].mid() for c in range(6)] for a in range(6)]).inv()
    except ZeroDivisionError:
        return False, X
    C = [[Cb[a, c].mid() for c in range(6)] for a in range(6)]
    for a in range(6):
        Ka = m[a] - sum(C[a][c] * Phi_m[c] for c in range(6))
        for c in range(6):
            Mac = (1 if a == c else 0) - sum(C[a][d] * J_X[d][c] for d in range(6))
            Ka += Mac * (X[c] - m[c])
        if not X[a].contains_interior(Ka):
            return False, X
    return True, X


def certainly_nonreal(X):
    """Some coordinate's imaginary part excludes 0 on the whole (complex) box."""
    return any(not x.imag.contains(0) for x in X)


# ----------------------------------------------------------------------------------- distinctness --
def plucker(b, X):
    i, j = b["i"], b["j"]
    K = [v for v in range(5) if v not in (i, j)]
    Y = [[None, None] for _ in range(5)]
    Y[i] = [acb(1), acb(0)]
    Y[j] = [acb(0), acb(1)]
    for m in range(3):
        Y[K[m]] = [acb(X[2 * m]), acb(X[2 * m + 1])]
    return [Y[a][0] * Y[c][1] - Y[a][1] * Y[c][0] for a in range(5) for c in range(a + 1, 5)]


def float_interval(x):
    """Float bounds lo <= x <= hi for a real ball x (padding covers the rounding to double)."""
    mid, rad = float(x.mid()), float(x.rad())
    d = rad * (1 + 2.0 ** -40) + abs(mid) * 2.0 ** -50 + 1e-300
    return mid - d, mid + d


def distinct(boxes, Xs):
    """Prove the lines pairwise distinct. For each line, normalise its Plücker vector p by a coordinate c with
    |p_c| >= max|p|/4. Two lines are distinct if (a) the coordinate where line a is largest is < 1/4 of the max
    for line b (the ratios |p_c|/max|p| differ), or (b) their normalised vectors are disjoint in some coordinate."""
    N = len(boxes)
    amax, norm = [], []
    for b, X in zip(boxes, Xs):
        p = plucker(b, X)
        mags = [float(abs(q).mid()) for q in p]
        M = max(mags)
        assert all(float(q.rad()) < 1e-30 * M for q in p), "Plücker balls too wide"
        amax.append(mags.index(M))
        nv = {}
        for c in range(10):
            if mags[c] >= 0.25 * M:
                v = [q / p[c] for q in p]
                fl = []
                for q in v:
                    fl += [*float_interval(q.real), *float_interval(q.imag)]
                nv[c] = (fl, v)
        norm.append(nv)
    n_ratio = n_float = n_ball = 0
    bad = []
    closest = float("inf")
    for a in range(N):
        c = amax[a]
        fa, va = norm[a][c]
        for b in range(a + 1, N):
            nb = norm[b]
            if c not in nb:
                n_ratio += 1
                continue
            fb, vb = nb[c]
            if any(fa[2 * k + 1] < fb[2 * k] or fb[2 * k + 1] < fa[2 * k] for k in range(20)):
                n_float += 1
                gap = max(abs(fa[2 * k] - fb[2 * k]) for k in range(20))
                closest = min(closest, gap)
            elif any(not va[k].overlaps(vb[k]) for k in range(10)):
                n_ball += 1
            else:
                bad.append((a, b))
    return bad, (n_ratio, n_float, n_ball), closest


# ------------------------------------------------------------------------------------------------ main --
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("quintic")
    ap.add_argument("lines")
    ap.add_argument("--prec", type=int, default=768, help="working precision in bits (default 768)")
    ap.add_argument("--controls", type=int, default=20, help="number of shifted boxes in the negative control")
    args = ap.parse_args()
    ctx.prec = args.prec

    mons = read_quintic(args.quintic)
    boxes = read_lines(args.lines)
    sys_ = LineSystem(mons)
    print(f"{args.quintic}: {len(mons)} monomials, exact coefficients")
    print(f"{args.lines}: {len(boxes)} boxes ({sum(b['real'] for b in boxes)} real, "
          f"{sum(not b['real'] for b in boxes)} complex); precision {args.prec} bits")

    t0 = time.time()
    Xs, failed, n_real, n_nonreal, unclassified = [], [], 0, 0, []
    for n, b in enumerate(boxes):
        ok, X = krawczyk(sys_, b)
        Xs.append(X)
        if not ok:
            failed.append(n)
        elif b["real"]:
            n_real += 1
        elif certainly_nonreal(X):
            n_nonreal += 1
        else:
            unclassified.append(n)
        if (n + 1) % 500 == 0:
            print(f"  ... {n + 1} boxes ({time.time() - t0:.0f} s)", flush=True)
    print(f"1. Krawczyk: {len(boxes) - len(failed)} of {len(boxes)} boxes contain exactly one line of X, "
          f"each reduced  ({time.time() - t0:.0f} s)")
    print(f"2. Reality: {n_real} real lines, {n_nonreal} non-real lines"
          + (f", {len(unclassified)} complex boxes meeting the real locus" if unclassified else ""))

    t0 = time.time()
    bad, (n_ratio, n_float, n_ball), closest = distinct(boxes, Xs)
    npairs = len(boxes) * (len(boxes) - 1) // 2
    print(f"3. Distinctness: {npairs} pairs; separated {npairs - len(bad)} "
          f"(ratio test {n_ratio}, float bounds {n_float}, balls {n_ball}); not separated {len(bad)}  "
          f"({time.time() - t0:.0f} s)")
    if closest < float("inf"):
        print(f"   closest pair of normalised Plücker vectors: {closest:.2e}")

    nc = min(args.controls, len(boxes))
    passed = sum(krawczyk(sys_, boxes[n], shift=10)[0] for n in range(nc))
    print(f"4. Negative control: {passed} of {nc} boxes shifted by 10r pass (should be 0)")

    ok = not failed and not unclassified and not bad and passed == 0
    if ok:
        print(f"RESULT: {len(boxes)} distinct reduced lines on X, of which {n_real} are real.")
        if len(boxes) == 2875:
            print("        2875 = the number of lines on a quintic threefold, counted with multiplicity, so these are")
            print(f"        all the lines of X: X has exactly 2875 lines, exactly {n_real} of them real.")
    else:
        print("RESULT: certificate NOT verified."
              + (f" Failed boxes: {failed[:10]}" if failed else "")
              + (f" Unseparated pairs: {bad[:10]}" if bad else ""))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
