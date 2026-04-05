"""
Protonation state calculator for Tris and Citrate buffers.
Ratios are normalised against the highest-charge species at each pH.
"""

import math


def monoprotic_fractions(pH, pKa):
    """Return (f_protonated, f_base) for a monoprotic acid."""
    H  = 10 ** -pH
    Ka = 10 ** -pKa
    D  = H + Ka
    return H / D, Ka / D          # (HA, A-)


def triprotic_fractions(pH, pKa1, pKa2, pKa3):
    """Return (f_H3A, f_H2A-, f_HA2-, f_A3-) for a triprotic acid."""
    H   = 10 ** -pH
    Ka1 = 10 ** -pKa1
    Ka2 = 10 ** -pKa2
    Ka3 = 10 ** -pKa3

    n0 = H**3
    n1 = Ka1 * H**2
    n2 = Ka1 * Ka2 * H
    n3 = Ka1 * Ka2 * Ka3
    D  = n0 + n1 + n2 + n3

    return n0/D, n1/D, n2/D, n3/D


def print_tris(pH):
    pKa = 8.06
    f_TrisH, f_base = monoprotic_fractions(pH, pKa)

    # Highest charge = Tris-H+ (protonated)
    # Normalise both against Tris-H+
    ref = f_TrisH

    print(f"\n{'='*55}")
    print(f"  Tris  |  pH {pH}  |  pKa {pKa}")
    print(f"{'='*55}")
    print(f"  {'Species':<18} {'Charge':>6}  {'Fraction':>10}  {'Ratio':>8}")
    print(f"  {'-'*47}")
    print(f"  {'Tris-H+':<18} {'(+1)':>6}  {f_TrisH:>9.1%}  {'1.00':>8}   ← reference")
    print(f"  {'Tris base':<18} {'(0)':>6}  {f_base:>9.1%}  {f_base/ref:>8.4f}")
    print(f"\n  log([base]/[acid]) = {math.log10(f_base/f_TrisH):.3f}")
    print(f"  [Tris-H+]/[Tris base] = {f_TrisH/f_base:.2f} : 1")


def print_citrate(pH):
    pKa1, pKa2, pKa3 = 3.13, 4.76, 6.40
    f0, f1, f2, f3 = triprotic_fractions(pH, pKa1, pKa2, pKa3)

    # Highest (most negative) charge = Cit3- at this pH range
    ref = f3  # normalise against Cit3-

    print(f"\n{'='*55}")
    print(f"  Citrate  |  pH {pH}  |  pKa1={pKa1}, pKa2={pKa2}, pKa3={pKa3}")
    print(f"{'='*55}")
    print(f"  {'Species':<18} {'Charge':>6}  {'Fraction':>10}  {'Ratio':>8}")
    print(f"  {'-'*47}")

    species = [
        ("H3Cit",   "(0)",  f0),
        ("H2Cit-",  "(-1)", f1),
        ("HCit2-",  "(-2)", f2),
        ("Cit3-",   "(-3)", f3),
    ]

    for name, charge, f in species:
        marker = "   ← reference" if name == "Cit3-" else ""
        print(f"  {name:<18} {charge:>6}  {f:>9.1%}  {f/ref:>8.4f}{marker}")

    avg_charge = -(0*f0 + 1*f1 + 2*f2 + 3*f3)
    print(f"\n  Average charge = {avg_charge:.2f}")


# ── Run calculations ──────────────────────────────────────

print_tris(7.4)
print_tris(7.2)
print_citrate(4.9)