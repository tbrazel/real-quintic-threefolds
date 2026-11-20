# Certificates for two real quintic threefolds

These files let you check two theorems yourself:

| quintic | file | claim |
|---|---|---|
| X₁₅ = {x₀⁵+⋯+x₄⁵ + (1/128) Σ_{i<j} (xᵢ−xⱼ)³ Σ_{k≠i,j} x_k² = 0} | `quintic_15.txt` | exactly 2875 lines, exactly **15** of them real |
| X₂₈₇₅ = {x₀x₁x₂x₃x₄ + t(L₁⋯L₅ + εg) = 0}, t = 2⁻⁴⁷, ε = 2⁻⁶⁰ | `quintic_2875.txt` | exactly 2875 lines, **all real** |

## Running it

```
pip install python-flint
python verify_lines.py quintic_15.txt   lines_15.txt
python verify_lines.py quintic_2875.txt lines_2875.txt
```

`python-flint` is the only dependency. It provides the Arb ball arithmetic. Each check takes about 15 seconds on
a laptop (tested with Python 3.15 and python-flint 0.9.0).

The last lines of the output should read:

```
RESULT: 2875 distinct reduced lines on X, of which 15 are real.
        2875 = the number of lines on a quintic threefold, counted with multiplicity, so these are
        all the lines of X: X has exactly 2875 lines, exactly 15 of them real.
```

and the same with `2875 of them real` for the second quintic. The exit code is 0 if and only if the certificate
is verified.

## Files

- `quintic_*.txt`: the polynomial, one monomial per row: `e0 e1 e2 e3 e4 num k`, meaning the term
  (num/2ᵏ)·x₀^e0⋯x₄^e4. Every coefficient is a dyadic rational and is read exactly. For X₂₈₇₅, the linear
  forms Lₛ and the quintic g have short dyadic coefficients m·2ᵉ, listed in the paper. The file is the exact
  expansion of F, and 2³¹²·F has integer coefficients.
- `lines_*.txt`: one box per line of X: `i j r R y1 … y6` for a real box, or `i j r C Re y1 Im y1 … Re y6 Im y6`
  for a complex one. The line is the span of the columns of the 5×2 matrix Y whose rows i, j (0-based) are
  (1,0), (0,1) and whose other rows, in order, are (y1,y2), (y3,y4), (y5,y6). The box is y ± r, with r = 10⁻⁶⁰.
- `verify_lines.py`: the checker (about 300 lines; the docstring at the top explains the argument).

## What the script proves

For a line in the chart (i, j), let Φₖ(y) = F(Y(y)·(1, τₖ)ᵀ) for the six values τₖ = ±3/2, ±3/4, ±1/4. Then:

- Φ(y) = 0 if and only if the line lies on X, because a binary quintic that vanishes at six points is zero.
- DΦ(y) is invertible if and only if the line is a reduced point of the Fano scheme of X, because the six values
  are an invertible (Vandermonde) image of the six coefficients of F restricted to the line.

The script runs four checks, all in interval arithmetic, so every inequality is proved rather than estimated:

1. **Krawczyk test.** For each box, K = m − CΦ(m) + (I − C·DΦ(box))(box − m) lies in the interior of the box.
   So the box contains exactly one line of X, and that line is reduced.
2. **Reality.**
   - A real box is tested in real arithmetic, so its line is real.
   - In a complex box, the imaginary part of some coordinate is bounded away from 0, so its line is not real.
     The chart is defined over ℝ, so a line is real exactly when its coordinates are real.
3. **Distinctness.** The Plücker vectors of the boxes are pairwise non-proportional.
4. **Negative control.** The same test fails on boxes moved by 10r off the line.

**Why there are no further lines.** The incidence variety {(F, ℓ) : ℓ ⊂ {F=0}} → ℙ¹²⁵ is proper of degree 2875.
It is étale at every reduced line. If X had any further line, or a curve of lines, a general quintic near X
would have more than 2875 lines. So 2875 distinct reduced lines are all the lines of X. See the paper for
details.

**Trust base.** The argument relies on Arb's interval arithmetic and on this script. Nothing depends on how the
approximate lines were found (homotopy continuation, then Newton's method at 768 bits). Any wrong or missing
input would make a check fail.
