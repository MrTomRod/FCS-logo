"""Build the crest: outlines, union, SVG export, geometric checks, overlay.

Run: python -m fcs.build
"""

from __future__ import annotations

import math
import subprocess
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pathops

from .geometry import Belt, Prim, Seg, arc_cubics, end_of, primitives, reverse, start_of
from .logo import Logo, Params, build, top_y

RED = "#f40500"
ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / (sys.argv[1] if len(sys.argv) > 1 else "out")
ORIGINAL = ROOT / "original" / "old_logo2.jpg"
# original raster: row 0 centre at y = 18.5 px, 4.2 px per unit (row pitch 63 px)
PX_SCALE, PX_ROW0 = 4.2, 18.5


def outline(belt: Belt, half: float) -> list[Prim]:
    """Closed outline of one stripe: left edge, end cap, right edge back, start cap."""
    lft = primitives(belt.offset(half), strict=False)
    rgt = primitives(belt.offset(-half), strict=False)
    return [
        *lft,
        Seg(end_of(lft[-1]), end_of(rgt[-1])),
        *[reverse(p) for p in reversed(rgt)],
        Seg(start_of(rgt[0]), start_of(lft[0])),
    ]


def polygon(pts: list[tuple[float, float]]) -> pathops.Path:
    path = pathops.Path()
    path.moveTo(*pts[0])
    for q in pts[1:]:
        path.lineTo(*q)
    path.close()
    return path


def union(contours: list[list[Prim]], logo: Logo) -> pathops.Path:
    path = pathops.Path(fillType=pathops.FillType.WINDING)
    for contour in contours:
        path.moveTo(*start_of(contour[0]))
        for prim in contour:
            if isinstance(prim, Seg):
                path.lineTo(*prim.b)
            else:
                for c1, c2, p3 in arc_cubics(prim):
                    path.cubicTo(*c1, *c2, *p3)
        path.close()
    out = pathops.op(path, polygon(logo.bar), pathops.PathOp.UNION)
    out.simplify(fix_winding=True)
    return out


def fmt(v: float) -> str:
    s = f"{v:.4f}".rstrip("0").rstrip(".")
    return "0" if s in ("-0", "") else s


class Frame:
    """Maps mathematical coordinates to SVG user units (y down, origin top-left)."""

    def __init__(self, x_min: float, y_top: float) -> None:
        self.x_min, self.y_top = x_min, y_top

    def pt(self, x: float, y: float) -> str:
        return f"{fmt(x - self.x_min)} {fmt(self.y_top - y)}"


def path_d(path: pathops.Path, fr: Frame) -> str:
    out = []
    for verb, pts in path.segments:
        if verb == "moveTo":
            out.append("M" + fr.pt(*pts[0]))
        elif verb == "lineTo":
            out.append("L" + fr.pt(*pts[0]))
        elif verb == "curveTo":
            out.append("C" + " ".join(fr.pt(*q) for q in pts))
        elif verb == "qCurveTo":
            out.append("Q" + " ".join(fr.pt(*q) for q in pts))
        elif verb == "closePath":
            out.append("Z")
        else:
            raise ValueError(verb)
    return "".join(out)


def exact_d(prims: list[Prim], fr: Frame, close: bool) -> str:
    """SVG path with true circular arcs (A commands); y flip turns CCW into sweep 0."""
    out = ["M" + fr.pt(*start_of(prims[0]))]
    for p in prims:
        if isinstance(p, Seg):
            out.append("L" + fr.pt(*p.b))
        else:
            large = 1 if abs(p.sweep) > math.pi else 0
            sweep = 0 if p.sweep > 0 else 1
            out.append(f"A{fmt(p.r)} {fmt(p.r)} 0 {large} {sweep} " + fr.pt(*p.point(1.0)))
    return "".join(out) + ("Z" if close else "")


def boundary(shape: pathops.Path, step: float) -> list[np.ndarray]:
    """Flatten the union outline into dense point lists, one per contour."""
    contours, cur, last = [], [], (0.0, 0.0)
    for verb, pts in shape.segments:
        if verb == "moveTo":
            cur, last = [pts[0]], pts[0]
        elif verb == "lineTo":
            n = max(1, math.ceil(math.dist(last, pts[0]) / step))
            cur += [(last[0] + (pts[0][0] - last[0]) * i / n, last[1] + (pts[0][1] - last[1]) * i / n)
                    for i in range(1, n + 1)]
            last = pts[0]
        elif verb == "curveTo":
            a, b, c = pts
            n = max(2, math.ceil((math.dist(last, a) + math.dist(a, b) + math.dist(b, c)) / step))
            for i in range(1, n + 1):
                t = i / n
                m = 1 - t
                cur.append(tuple(m**3 * last[k] + 3 * m * m * t * a[k] + 3 * m * t * t * b[k] + t**3 * c[k]
                                 for k in range(2)))
            last = c
        elif verb == "closePath":
            contours.append(np.array(cur))
    return contours


def checks(p: Params, belts: list[Belt], shape: pathops.Path) -> list[str]:
    """Every white gap >= white, measured between red boundary points across white."""
    bad = []
    half = p.red / 2
    for b in belts:
        for pul in b.pulleys:
            if pul.lines is None and 1e-9 < abs(pul.r) < half - 1e-9:
                bad.append(f"{b.name}: pulley radius {abs(pul.r):.3f} < red/2")
    step = 0.1
    conts = boundary(shape, step)
    pts = np.concatenate(conts)
    cid = np.concatenate([np.full(len(c), i) for i, c in enumerate(conts)])
    pos = np.concatenate([np.arange(len(c)) for c in conts])
    size = np.array([len(c) for c in conts])
    cells: dict[tuple[int, int], list[int]] = {}
    keys = np.floor(pts / p.white).astype(int)
    for i, k in enumerate(map(tuple, keys)):
        cells.setdefault(k, []).append(i)
    min_gap, where, min_sep = math.inf, None, 3 * p.pitch / step
    for (kx, ky), idx in cells.items():
        near = [j for dx in (-1, 0, 1) for dy in (-1, 0, 1) for j in cells.get((kx + dx, ky + dy), [])]
        a, b = np.array(idx), np.array(near)
        d = np.linalg.norm(pts[a][:, None] - pts[b][None], axis=2)
        sep = np.abs(pos[a][:, None] - pos[b][None])
        sep = np.minimum(sep, size[cid[a]][:, None] - sep)
        far = (cid[a][:, None] != cid[b][None]) | (sep > min_sep)
        cand = np.argwhere(far & (d < min_gap))
        for i, j in cand:
            if d[i, j] >= min_gap:
                continue
            mid = (pts[a[i]] + pts[b[j]]) / 2
            if not shape.contains(tuple(mid)):
                min_gap, where = float(d[i, j]), mid
    report = [f"boundary points: {len(pts)}, smallest white gap: {min_gap:.5f} "
              f"(target {p.white}) near x={where[0]:.1f} y={where[1]:.1f}"]
    if min_gap < p.white - 1e-3:
        bad.append(report[0])
    return report + ([] if not bad else ["", "FAILURES:", *bad])


def write_svg(file: Path, w: float, h: float, body: str) -> None:
    file.write_text(
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {fmt(w)} {fmt(h)}" '
        f'width="{fmt(w * 4)}" height="{fmt(h * 4)}">\n{body}\n</svg>\n'
    )


def main() -> int:
    p = Params()
    logo = build(p)
    belts = logo.belts
    for b in belts:
        primitives(b)  # centerlines must keep every straight on its grid line
    half = p.red / 2
    contours = [outline(b, half) for b in belts]
    shape = union(contours, logo)
    x0, y0, x1, y1 = shape.bounds
    fr = Frame(x0, y1)
    w, h = x1 - x0, y1 - y0
    OUT.mkdir(exist_ok=True)

    d = path_d(shape, fr)
    write_svg(OUT / "fcs.svg", w, h, f'<path fill="{RED}" d="{d}"/>')

    lines = [
        f'<path fill="none" stroke="#000" stroke-width="0.3" d="{exact_d(primitives(b), fr, False)}"/>'
        for b in belts
    ]
    edges = [f'<path fill="{RED}" fill-opacity="0.35" d="{exact_d(c, fr, True)}"/>' for c in contours]
    write_svg(OUT / "construction.svg", w, h, "\n".join(edges + lines))

    # overlay on the original raster (px coordinates of old_logo2.jpg)
    s = PX_SCALE
    ty = PX_ROW0 + 8 * p.pitch * s
    raw = Frame(0.0, 0.0)
    d_raw = path_d(shape, raw)  # x unchanged, y negated
    overlay = (
        f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
        f'width="1193" height="541" viewBox="0 0 1193 541">\n'
        f'<path transform="matrix({s} 0 0 {s} 0 {fmt(ty)})" fill="#0050ff" fill-opacity="0.35" '
        f'stroke="#002080" stroke-width="{fmt(0.8 / s)}" d="{d_raw}"/>\n</svg>\n'
    )
    (OUT / "overlay.svg").write_text(overlay)
    subprocess.run(["magick", "-background", "white", "-density", "96", str(OUT / "fcs.svg"),
                    "-flatten", str(OUT / "fcs.png")], check=True)
    subprocess.run(["magick", str(ORIGINAL), "-fill", "white", "-colorize", "55%",
                    "(", "-background", "none", "-density", "96", str(OUT / "overlay.svg"), ")",
                    "-composite", str(OUT / "overlay.png")], check=True)

    report = checks(p, belts, shape)
    print(f"size {w:.3f} x {h:.3f} units (top {top_y(p):.3f})")
    print(f"derived S bowl inner radius s_bowl = {p.s_bowl:.4f}")
    print("\n".join(report))
    return 1 if "FAILURES:" in report else 0


if __name__ == "__main__":
    sys.exit(main())
