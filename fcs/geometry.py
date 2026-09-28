"""Exact planar geometry for stripe centerlines built from lines and circular arcs.

Coordinates are mathematical (y points up). A stripe centerline is a "belt":
it starts on a line, wraps around a sequence of pulleys (circles), and ends on
a line. Offsetting a belt by a distance d is exact: lines shift by d and every
pulley keeps its centre while its signed radius changes by -d. That property
is what guarantees equal stripe spacing.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

Vec = tuple[float, float]
EPS = 1e-9


def add(a: Vec, b: Vec) -> Vec:
    return (a[0] + b[0], a[1] + b[1])


def sub(a: Vec, b: Vec) -> Vec:
    return (a[0] - b[0], a[1] - b[1])


def mul(a: Vec, k: float) -> Vec:
    return (a[0] * k, a[1] * k)


def dot(a: Vec, b: Vec) -> float:
    return a[0] * b[0] + a[1] * b[1]


def cross(a: Vec, b: Vec) -> float:
    return a[0] * b[1] - a[1] * b[0]


def norm(a: Vec) -> float:
    return math.hypot(a[0], a[1])


def left(u: Vec) -> Vec:
    """Left-hand normal of a direction (counter-clockwise rotation by 90 deg)."""
    return (-u[1], u[0])


def direction(deg: float) -> Vec:
    return (math.cos(math.radians(deg)), math.sin(math.radians(deg)))


@dataclass(frozen=True)
class Line:
    """Directed line through point q with unit direction u."""

    q: Vec
    u: Vec

    def shifted(self, d: float) -> Line:
        return Line(add(self.q, mul(left(self.u), d)), self.u)

    def reversed(self) -> Line:
        return Line(self.q, mul(self.u, -1.0))

    def distance(self, p: Vec) -> float:
        """Signed distance of p from the line, positive on the left side."""
        return dot(sub(p, self.q), left(self.u))


def intersect(a: Line, b: Line) -> Vec:
    den = cross(a.u, b.u)
    if abs(den) < EPS:
        raise ValueError("parallel lines do not intersect")
    t = cross(sub(b.q, a.q), b.u) / den
    return add(a.q, mul(a.u, t))


@dataclass(frozen=True)
class Pulley:
    """Circle wrapped by a belt. Signed radius: > 0 turns left, < 0 turns right.

    Special kinds:
    - fillet (`lines` set, hard False): rounded corner between two lines. Any
      radius >= 0 is allowed; where an offset would make the radius negative
      (inner side of a tight corner) the edge becomes a mitred vertex, which is
      the exact offset of the corner;
    - sharp (`lines` set, hard True): mitred on both sides;
    - knee: the belt arrives along `knee` and enters the circle at the
      intersection, without tangency (hard corner), then wraps normally.
    """

    c: Vec
    r: float
    lines: tuple[Line, Line] | None = None
    knee: Line | None = None
    hard: bool = False

    def offset(self, d: float) -> Pulley:
        if self.lines:
            a, b = self.lines[0].shifted(d), self.lines[1].shifted(d)
            s = 1.0 if cross(a.u, b.u) > 0 else -1.0
            radius = s * self.r - s * d
            if self.hard or radius <= 0:
                return Pulley(intersect(a, b), 0.0, (a, b), None, self.hard)
            return corner(a, b, radius)
        return Pulley(self.c, self.r - d, None, self.knee.shifted(d) if self.knee else None)


def sharp(l_in: Line, l_out: Line) -> Pulley:
    """Hard mitred corner between two directed lines."""
    return Pulley(intersect(l_in, l_out), 0.0, (l_in, l_out), None, True)


def knee(c: Vec, r: float, l_in: Line) -> Pulley:
    """Circle entered with a hard corner from line l_in."""
    return Pulley(c, r, None, l_in)


def circle_entry(line: Line, p: Pulley) -> Vec:
    """First intersection of a directed line with the pulley circle."""
    b = dot(sub(line.q, p.c), line.u)
    cc = dot(sub(line.q, p.c), sub(line.q, p.c)) - p.r * p.r
    disc = b * b - cc
    if disc < 0:
        raise ValueError(f"knee line misses circle {p}")
    return add(line.q, mul(line.u, -b - math.sqrt(disc)))


def corner(l_in: Line, l_out: Line, radius: float) -> Pulley:
    """Fillet pulley of the given radius (>= 0) between two directed lines."""
    s = 1.0 if cross(l_in.u, l_out.u) > 0 else -1.0
    c = intersect(l_in.shifted(s * radius), l_out.shifted(s * radius))
    return Pulley(c, s * radius, (l_in, l_out))


def tangent(p1: Pulley, p2: Pulley) -> tuple[Vec, Vec, Vec]:
    """Belt segment leaving p1 and arriving at p2: (start, end, direction)."""
    d = sub(p2.c, p1.c)
    dist = norm(d)
    k = (p2.r - p1.r) / dist
    if abs(k) > 1.0 + 1e-9:
        raise ValueError(f"no belt tangent between {p1} and {p2}")
    k = max(-1.0, min(1.0, k))
    e = mul(d, 1.0 / dist)
    n = add(mul(e, k), mul(left(e), math.sqrt(max(0.0, 1.0 - k * k))))
    u = (n[1], -n[0])
    return sub(p1.c, mul(n, p1.r)), sub(p2.c, mul(n, p2.r)), u


@dataclass(frozen=True)
class Seg:
    a: Vec
    b: Vec


@dataclass(frozen=True)
class Arc:
    c: Vec
    r: float  # absolute radius
    a0: float  # start angle (rad)
    sweep: float  # signed sweep (rad), > 0 counter-clockwise

    def point(self, t: float) -> Vec:
        ang = self.a0 + self.sweep * t
        return (self.c[0] + self.r * math.cos(ang), self.c[1] + self.r * math.sin(ang))


Prim = Seg | Arc


def start_of(p: Prim) -> Vec:
    return p.a if isinstance(p, Seg) else p.point(0.0)


def end_of(p: Prim) -> Vec:
    return p.b if isinstance(p, Seg) else p.point(1.0)


def reverse(p: Prim) -> Prim:
    if isinstance(p, Seg):
        return Seg(p.b, p.a)
    return Arc(p.c, p.r, p.a0 + p.sweep, -p.sweep)


@dataclass(frozen=True)
class Belt:
    """Open stripe centerline: start line, pulleys, end line, plus the two cut lines."""

    name: str
    start: Line
    pulleys: tuple[Pulley, ...]
    end: Line
    start_cut: Line
    end_cut: Line

    def offset(self, d: float) -> Belt:
        return Belt(
            self.name,
            self.start.shifted(d),
            tuple(p.offset(d) for p in self.pulleys),
            self.end.shifted(d),
            self.start_cut,
            self.end_cut,
        )


def _turn(u_in: Vec, u_out: Vec, r: float) -> float:
    ang = math.atan2(cross(u_in, u_out), dot(u_in, u_out))
    if abs(ang) < 1e-12:
        return 0.0
    if r > 0 and ang < 0:
        ang += 2 * math.pi
    if r < 0 and ang > 0:
        ang -= 2 * math.pi
    return ang


def primitives(belt: Belt) -> list[Prim]:
    """Exact line/arc decomposition of a belt. Raises if the belt folds back."""
    points: list[Vec] = [intersect(belt.start, belt.start_cut)]
    dirs: list[Vec] = [belt.start.u]
    # tangent points: entry/exit per pulley
    entries: list[Vec] = []
    exits: list[Vec] = []
    first = belt.pulleys[0]
    entries.append(sub(first.c, mul(left(belt.start.u), first.r)))
    if abs(belt.start.distance(entries[0])) > 1e-6:
        raise ValueError(f"{belt.name}: first pulley not tangent to start line")
    for p1, p2 in zip(belt.pulleys, belt.pulleys[1:]):
        if p2.knee:
            k = p2.knee
            a = p1.c if abs(p1.r) < EPS else sub(p1.c, mul(left(k.u), p1.r))
            if abs(k.distance(a)) > 1e-6:
                raise ValueError(f"{belt.name}: knee line does not leave previous pulley")
            a, b, u = a, circle_entry(k, p2), k.u
        else:
            a, b, u = tangent(p1, p2)
        exits.append(a)
        entries.append(b)
        dirs.append(u)
    last = belt.pulleys[-1]
    exits.append(sub(last.c, mul(left(belt.end.u), last.r)))
    if abs(belt.end.distance(exits[-1])) > 1e-6:
        raise ValueError(f"{belt.name}: last pulley not tangent to end line")
    dirs.append(belt.end.u)
    end_pt = intersect(belt.end, belt.end_cut)

    prims: list[Prim] = []
    seg_starts = [points[0], *exits]
    seg_ends = [*entries, end_pt]
    for i, (a, b) in enumerate(zip(seg_starts, seg_ends)):
        if dot(sub(b, a), dirs[i]) < -1e-7:
            raise ValueError(f"{belt.name}: segment {i} runs backwards")
        if norm(sub(b, a)) > 1e-9:
            prims.append(Seg(a, b))
        if i < len(belt.pulleys):
            pul = belt.pulleys[i]
            a0 = math.atan2(entries[i][1] - pul.c[1], entries[i][0] - pul.c[0])
            if pul.knee:  # arrival is not tangent: sweep from the two angles
                a1 = math.atan2(exits[i][1] - pul.c[1], exits[i][0] - pul.c[0])
                sweep = (a1 - a0) % (2 * math.pi) if pul.r > 0 else -((a0 - a1) % (2 * math.pi))
            else:
                sweep = _turn(dirs[i], dirs[i + 1], pul.r)
            if abs(pul.r) > 1e-9 and abs(sweep) > 1e-12:
                prims.append(Arc(pul.c, abs(pul.r), a0, sweep))
    return prims


def sample(prims: list[Prim], step: float = 0.25) -> list[Vec]:
    pts: list[Vec] = []
    for p in prims:
        length = norm(sub(p.b, p.a)) if isinstance(p, Seg) else p.r * abs(p.sweep)
        n = max(1, math.ceil(length / step))
        for i in range(n):
            t = i / n
            pts.append(add(p.a, mul(sub(p.b, p.a), t)) if isinstance(p, Seg) else p.point(t))
    pts.append(end_of(prims[-1]))
    return pts


def arc_cubics(arc: Arc, max_step: float = math.pi / 6) -> list[tuple[Vec, Vec, Vec]]:
    """Cubic Bezier approximation of an arc (max error ~4e-7 * r at 30 deg)."""
    n = max(1, math.ceil(abs(arc.sweep) / max_step))
    step = arc.sweep / n
    k = 4.0 / 3.0 * math.tan(step / 4.0)
    out = []
    for i in range(n):
        a0 = arc.a0 + step * i
        a1 = a0 + step
        p0 = (arc.c[0] + arc.r * math.cos(a0), arc.c[1] + arc.r * math.sin(a0))
        p3 = (arc.c[0] + arc.r * math.cos(a1), arc.c[1] + arc.r * math.sin(a1))
        t0 = (-math.sin(a0), math.cos(a0))
        t1 = (-math.sin(a1), math.cos(a1))
        out.append((add(p0, mul(t0, k * arc.r)), sub(p3, mul(t1, k * arc.r)), p3))
    return out
