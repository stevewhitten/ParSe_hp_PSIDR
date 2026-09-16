#!/usr/bin/env python3
"""
parse_hp_psidr.py

ParSe v2 plus a high-pathogenicity PS-IDR scan.
Window labels, domain growing (>=20 residues, >=90% one label),
embedded-domain splits, and adjacent-overlap repair follow
parse_hp_psidr.f90.

Labels:
  F  folded
  D  conventional ID
  P  phase-separating ID
Second track (disordered windows only):
  X  pathogenic-site-like ID
  Y  other ID

High-pathogenicity PS-IDRs are X-runs that also sit on the P side
of the helix / nu_model line. X-runs are not sent through the
embedded-domain or adjacent-overlap repair used for F/D/P domains.

Usage:
  python3 parse_hp_psidr.py SEQUENCE
  python3 parse_hp_psidr.py proteome.fasta
"""

from __future__ import annotations

import math
import os
import sys
from typing import List, Sequence, Tuple

MAXN = 10000
WIN = 25
PCT_CUT = 0.90
HYDR_CUT = 0.08280152
PD_M = -0.244078945
PD_B = 0.7885823
HELIX_PS = 0.9327272
NU_PS = 0.5416
HELIX_ID = 1.022552
NU_ID = 0.5582901
XY_M = 20.0788
XY_B = -20.0254

# Scale order: A C D E F G H I K L M N P Q R S T V W Y
AA = "ACDEFGHIKLMNPQRSTVWY"
AA_INDEX = {c: i for i, c in enumerate(AA)}

PPII = [
    0.37, 0.25, 0.30, 0.42, 0.17, 0.13, 0.20, 0.39, 0.56, 0.24,
    0.36, 0.27, 1.00, 0.53, 0.38, 0.24, 0.32, 0.39, 0.25, 0.25,
]
HELIX_SC = [
    1.42, 0.73, 1.01, 1.63, 1.16, 0.50, 1.20, 1.12, 1.24, 1.29,
    1.21, 0.71, 0.65, 1.02, 1.06, 0.71, 0.78, 0.99, 1.05, 0.67,
]
HYDR_SC = [
    0.0728, 0.3557, -0.0552, -0.0295, 0.4201, -0.0589, 0.0874, 0.3805, -0.0053, 0.3819,
    0.1613, -0.0390, -0.0492, 0.0126, 0.0394, -0.0282, 0.0239, 0.2947, 0.4114, 0.3113,
]
SHEET_SC = [
    0.90, 1.24, 0.47, 0.62, 1.23, 0.56, 1.12, 1.54, 0.74, 1.26,
    1.09, 0.62, 0.42, 1.18, 1.02, 0.87, 1.30, 1.53, 1.75, 1.68,
]
# Radzicka-Wolfenden vap-to-octanol; C, D, P missing and set to 0
HYDR2_SC = [
    1.42, 0.00, 0.00, -9.45, -2.85, 2.39, -11.22, 0.11, -9.60, 0.52,
    -2.80, -9.67, 0.00, -9.31, -18.60, -5.10, -5.15, 0.81, -8.39, -7.74,
]


def window_props(cnt: Sequence[int], nres: int) -> Tuple[float, float, float, float, float]:
    hydr = helix = sheet = hydr2 = fppii = 0.0
    q = abs((cnt[2] + cnt[3]) - (cnt[8] + cnt[14]))  # |D+E - K+R|
    for k in range(20):
        n = cnt[k]
        hydr += n * HYDR_SC[k]
        helix += n * HELIX_SC[k]
        sheet += n * SHEET_SC[k]
        hydr2 += n * HYDR2_SC[k]
        fppii += n * PPII[k]
    rn = float(nres)
    hydr /= rn
    helix /= rn
    sheet /= rn
    hydr2 /= rn
    fppii /= rn
    if fppii == 1.0:
        fppii = 0.98
    vexp = 0.503 - 0.11 * math.log(1.0 - fppii)
    rh = (
        2.16 * ((4 * nres) ** vexp)
        + 0.26 * (4 * q)
        - 0.29 * math.sqrt(4 * nres)
    )
    nu = math.log(rh / 2.16) / math.log(4 * nres)
    return hydr, helix, nu, sheet, hydr2


def count_segment(seq: str, start: int, end: int) -> List[int]:
    """Counts for seq[start:end] (Python slice, end exclusive)."""
    cnt = [0] * 20
    for c in seq[start:end]:
        cnt[AA_INDEX[c]] += 1
    return cnt


def grow_runs(lab: Sequence[str], want: str) -> List[Tuple[int, int]]:
    """1-based inclusive runs using the original grow-from-20 / >=90% rule."""
    npep = len(lab)
    runs: List[Tuple[int, int]] = []
    i = 1
    while i + 19 <= npep:
        count_p = 0
        count_w = 0
        region = 0
        # DO j = i, i+19. After a completed Fortran DO, j is i+20 (gfortran).
        for jj in range(i, i + 20):
            count_w += 1
            if lab[jj - 1] == want:
                count_p += 1
        j = i + 20
        while True:
            percent_p = count_p / float(count_w)
            if percent_p < PCT_CUT:
                break
            region = 1
            pstart = i
            pend = j
            if j >= npep:
                break
            j += 1
            count_w += 1
            if lab[j - 1] == want:
                count_p += 1
        if region == 1:
            runs.append((pstart, pend))
            i = j
        else:
            i += 1
    return runs


def split_embedded(
    a_runs: List[Tuple[int, int]],
    b_runs: Sequence[Tuple[int, int]],
) -> List[Tuple[int, int]]:
    """If a B-run lies strictly inside an A-run, trim or split A."""
    a = list(a_runs)
    restart = True
    while restart:
        restart = False
        for i, (a0, a1) in enumerate(a):
            for b0, b1 in b_runs:
                if b0 > a0 and b1 < a1:
                    nterm = b0 - a0
                    cterm = a1 - b1
                    if nterm > 20 and cterm <= 20:
                        a[i] = (a0, b0 - 1)
                    elif nterm <= 20 and cterm > 20:
                        a[i] = (b1 + 1, a1)
                    elif nterm > 20 and cterm > 20:
                        a[i] = (a0, b0 - 1)
                        a.append((b1 + 1, a1))
                        restart = True
                        break
            if restart:
                break
    return a


def build_order(
    p_runs: Sequence[Tuple[int, int]],
    d_runs: Sequence[Tuple[int, int]],
    f_runs: Sequence[Tuple[int, int]],
) -> List[int]:
    starts = [s for s, _ in p_runs] + [s for s, _ in d_runs] + [s for s, _ in f_runs]
    return sorted(starts)


def trim_adjacent(
    p_runs: List[Tuple[int, int]],
    d_runs: List[Tuple[int, int]],
    f_runs: List[Tuple[int, int]],
    order: List[int],
) -> List[int]:
    extra = int(round((1.0 - PCT_CUT) * 20.0 / 2.0))  # Fortran int() of 1.0; Python int(0.999...) is 0
    if len(order) <= 1:
        return order

    def _update(runs: List[Tuple[int, int]], start: int, new_start: int = None, new_end: int = None):
        for k, (s, e) in enumerate(runs):
            if s == start:
                ns = s if new_start is None else new_start
                ne = e if new_end is None else new_end
                runs[k] = (ns, ne)
                return ns
        return start

    for i in range(len(order) - 1):
        cur = order[i]
        nxt = order[i + 1]
        for runs in (p_runs, d_runs, f_runs):
            for s, e in list(runs):
                if s == cur and e >= nxt:
                    new_end = nxt + extra
                    _update(runs, cur, new_end=new_end)
                    for other in (p_runs, d_runs, f_runs):
                        ns = _update(other, nxt, new_start=new_end + 1)
                        if ns != nxt:
                            order[i + 1] = ns
                            nxt = ns
    return order


def classify(seq: str) -> Tuple[str, str, float]:
    npep = len(seq)
    lab = [" "] * npep
    lab2 = [" "] * npep

    m = -1.0 / PD_M
    b = NU_ID - m * HELIX_ID
    x = (b - PD_B) / (PD_M - m)
    y = m * x + b
    id_dist = math.hypot(HELIX_ID - x, NU_ID - y)
    b = NU_PS - m * HELIX_PS
    x = (b - PD_B) / (PD_M - m)
    y = m * x + b
    ps_dist = math.hypot(HELIX_PS - x, NU_PS - y)

    p_sum = 0.0
    for i in range(0, npep - WIN + 1):
        mid = i + WIN // 2
        cnt = count_segment(seq, i, i + WIN)
        hydr, helix, nu, sheet, hydr2 = window_props(cnt, WIN)
        if hydr >= HYDR_CUT:
            lab[mid] = "F"
            lab2[mid] = "F"
            continue
        lab2[mid] = "X" if hydr2 > (XY_M * sheet + XY_B) else "Y"
        m = -1.0 / PD_M
        b = nu - m * helix
        x = (b - PD_B) / (PD_M - m)
        y = m * x + b
        if ((nu - PD_B) / PD_M) <= helix:
            lab[mid] = "D"
        else:
            lab[mid] = "P"
            p_sum += math.hypot(helix - x, nu - y) / ps_dist

    half = WIN // 2
    for j in range(half):
        lab[j] = lab[half]
        lab2[j] = lab2[half]
    for j in range(npep - half, npep):
        lab[j] = lab[npep - half - 1]
        lab2[j] = lab2[npep - half - 1]
    return "".join(lab), "".join(lab2), p_sum


def analyze_one(hdr: str, seq: str) -> None:
    npep = len(seq)
    if npep < WIN:
        print("could not parse sequence")
        print("input sequence is too short")
        return
    if any(c not in AA_INDEX for c in seq):
        print("could not parse sequence")
        print("sequence contains noncommon amino acid type")
        return
    if npep > MAXN:
        print("could not parse sequence")
        print("input sequence is too long")
        return

    lab, lab2, p_sum = classify(seq)
    p_runs = grow_runs(lab, "P")
    d_runs = grow_runs(lab, "D")
    f_runs = grow_runs(lab, "F")
    x_runs = grow_runs(lab2, "X")

    p_runs = split_embedded(p_runs, d_runs)
    p_runs = split_embedded(p_runs, f_runs)
    d_runs = split_embedded(d_runs, p_runs)
    d_runs = split_embedded(d_runs, f_runs)
    f_runs = split_embedded(f_runs, d_runs)
    f_runs = split_embedded(f_runs, p_runs)

    order = build_order(p_runs, d_runs, f_runs)
    order = trim_adjacent(p_runs, d_runs, f_runs, order)

    hp = []
    for x0, x1 in x_runs:
        cnt = count_segment(seq, x0 - 1, x1)
        _hydr, helix, nu, _sheet, _hydr2 = window_props(cnt, x1 - x0 + 1)
        if ((nu - PD_B) / PD_M) > helix:
            hp.append((x0, x1))

    print()
    if hdr.strip():
        print(f"> {hdr}")
    print(f"length {npep}")
    print(f"PS potential (summed P-window classifier distance) {p_sum:10.3f}")
    print("ParSe labels (F/D/P):")
    print(lab)
    print("pathogenic-site ID labels (X) vs other ID (Y) or folded (F):")
    print(lab2)
    print()
    print("domains (>=20 residues, >=90% one label)")
    lookup = {s: ("PS-ID (P)    ", e) for s, e in p_runs}
    lookup.update({s: ("nonPS-ID (D) ", e) for s, e in d_runs})
    lookup.update({s: ("folded (F)   ", e) for s, e in f_runs})
    for start in order:
        tag, end = lookup[start]
        print(f"{tag} first {start}  last {end}  length {end - start + 1}")
    print()
    print(f"high-pathogenicity PS-IDRs = {len(hp)}")
    if hp:
        print("index   first    last   length")
        for i, (a, b) in enumerate(hp, 1):
            print(f"{i:5d}  {a:7d}  {b:7d}  {b - a + 1:7d}")


def load_raw(s: str) -> str:
    out = []
    for c in s:
        if c == " ":
            break
        if c in "\t\n\r":
            continue
        out.append(c.upper())
    return "".join(out)


def read_fasta(path: str):
    hdr = None
    chunks: List[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if hdr is not None:
                    yield hdr, "".join(chunks).upper()
                hdr = line[1:]
                chunks = []
            else:
                if hdr is None:
                    hdr = "unnamed"
                chunks.append("".join(c for c in line if not c.isspace()))
        if hdr is not None:
            yield hdr, "".join(chunks).upper()


def main(argv: List[str]) -> int:
    if len(argv) < 2:
        print("no input argument, exiting program")
        return 1
    arg = argv[1]
    if os.path.isfile(arg):
        nrec = 0
        for hdr, seq in read_fasta(arg):
            nrec += 1
            analyze_one(hdr, seq)
        if nrec == 0:
            print("no sequences found in file")
            return 1
        return 0
    analyze_one("", load_raw(arg))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
