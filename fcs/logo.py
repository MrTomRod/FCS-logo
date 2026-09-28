"""Construction of the FCS crest from a handful of design parameters.

Grid (mathematical coordinates, y up, units arbitrary):
- 9 horizontal rows, centre of row k at y = (8 - k) * pitch. Row 0 is the top.
- Slanted columns at angle `slant` (deg), identified by their x at row 0.
  Neighbouring columns sit exactly one pitch apart (perpendicular distance).
- Red stripe width `red`, white gap `white`, pitch = red + white.

The logo is 6 continuous stripes plus the F crossbar (a straight bar joining
B and D with hard corners). Concentric corners share a centre and their radii
differ by exactly one pitch; hard corners are mitred, so their offsets stay
exactly parallel too. Both keep every white gap equal.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .geometry import Belt, Line, Pulley, Vec, corner, direction, intersect, knee, sharp


@dataclass(frozen=True)
class Params:
    red: float = 8.0
    white: float = 7.0  # white / red = 0.875 as in the original
    slant: float = 70.0
    # columns: x of the outermost stripe centre at row 0
    col_f: float = 50.0  # F stem, outer stripe
    col_c: float = 126.0  # C outer stripe
    col_s_top: float = 204.0  # S upper-left, outer stripe
    col_s_bottom: float = 252.0  # S lower-right, inner stripe (see s_bowl/s_waist)
    # innermost corner radii of each concentric corner group (centerline)
    # (all corners around the F crossbar are hard, mitred corners)
    f_top: float = 18.0  # F top-left (H stripe; B, A are +1, +2 pitch)
    c_top: float = 32.0  # C top-left
    c_bowl: float = 20.0  # C bottom-left
    s_foot: float = 4.0  # S bottom-right
    # S middle runs horizontally on rows 3 (Main), 4 (D), 5 (E). Its length for
    # D is (col_s_bottom - col_s_top) - (s_bowl + s_waist + 2 * pitch) * 1.428
    # (1.428 = tan 55 deg for the 110 deg turns); it must stay >= 0.
    s_bowl: float = 4.0  # S lower bowl, top-right
    s_waist: float = 1.0  # S middle, left
    s_top: float = 10.0  # S top-left

    @property
    def pitch(self) -> float:
        return self.red + self.white

    @property
    def col_step(self) -> float:
        """Horizontal distance between neighbouring slanted columns."""
        return self.pitch / math.sin(math.radians(self.slant))


@dataclass(frozen=True)
class Logo:
    belts: list[Belt]
    bar: list[Vec]


def build(p: Params) -> Logo:
    P, step = p.pitch, p.col_step
    up = direction(p.slant)
    down = (-up[0], -up[1])
    top = 8 * P

    def row(k: int, east: bool = True) -> Line:
        return Line((0.0, (8 - k) * P), (1.0, 0.0) if east else (-1.0, 0.0))

    def col(x0: float, j: float = 0, going_up: bool = True) -> Line:
        return Line((x0 + j * step, top), up if going_up else down)

    def cut(x0: float) -> Line:
        return Line((x0, top), up)

    half = p.red / 2
    edge = half / math.sin(math.radians(p.slant))  # horizontal half stripe width
    gap_cut = (half + p.white) / math.sin(math.radians(p.slant))
    bottom = Line((0.0, -half), (1.0, 0.0))
    cut_f = cut(p.col_c - gap_cut)  # one white gap before the C outer stripe
    cut_c = cut(p.col_s_top - gap_cut)  # one white gap before the S outer stripe
    cut_s = cut(p.col_s_bottom + 2 * step + edge)  # flush with the S bowl

    f0, f1, f2 = col(p.col_f), col(p.col_f, 1), col(p.col_f, 2)
    c0, c1, c2 = (col(p.col_c, j) for j in range(3))
    su = [col(p.col_s_top, j) for j in range(3)]  # E, D, Main (upper S)
    sl = [col(p.col_s_bottom, j) for j in range(3)]  # E, D, Main (lower S)

    def rev(line: Line) -> Line:
        return line.reversed()

    # C bowl (shared centre for E, D; Main follows D one pitch further out)
    bowl_d = corner(rev(c1), row(7), p.c_bowl + P)
    bowl_e = corner(rev(c2), row(6), p.c_bowl)
    # Main runs along row 5 and drops straight onto the bowl circle one pitch
    # outside D (hard corner), so it stays exactly one pitch from D.
    main_bowl = knee(bowl_d.c, p.c_bowl + 2 * P, row(5))

    def s_bowl(k: int) -> Pulley:  # k = 0 E (inner) .. 2 Main (outer)
        return corner(sl[k], row(5 - k, east=False), p.s_bowl + k * P)

    def s_waist(k: int) -> Pulley:  # k = 0 Main (inner) .. 2 E (outer)
        return corner(row(3 + k, east=False), su[2 - k], p.s_waist + k * P)

    belts = [
        Belt("A", f0, (corner(f0, row(0), p.f_top + 2 * P),), row(0), bottom, cut_f),
        Belt("B", f1, (corner(f1, row(1), p.f_top + P),), row(1), bottom, cut_f),
        Belt(
            "H",
            row(2, east=False),
            (
                corner(row(2, False), rev(f2), p.f_top),
                sharp(rev(f2), row(3)),
                sharp(row(3), c0),
                corner(c0, row(0), p.c_top + 2 * P),
            ),
            row(0),
            cut_f,
            cut_c,
        ),
        Belt(
            "Main",
            f2,
            (
                sharp(f2, row(5)),
                main_bowl,
                corner(row(8), sl[2], p.s_foot + 2 * P),
                s_bowl(2),
                s_waist(0),
                corner(su[2], row(2), p.s_top),
            ),
            row(2),
            bottom,
            cut_s,
        ),
        Belt(
            "D",
            row(1, east=False),
            (
                corner(row(1, False), rev(c1), p.c_top + P),
                bowl_d,
                corner(row(7), sl[1], p.s_foot + P),
                s_bowl(1),
                s_waist(1),
                corner(su[1], row(1), p.s_top + P),
            ),
            row(1),
            cut_c,
            cut_s,
        ),
        Belt(
            "E",
            row(2, east=False),
            (
                corner(row(2, False), rev(c2), p.c_top),
                bowl_e,
                corner(row(6), sl[0], p.s_foot),
                s_bowl(0),
                s_waist(2),
                corner(su[0], row(0), p.s_top + 2 * P),
            ),
            row(0),
            cut_c,
            cut_s,
        ),
    ]

    # crossbar (row 4): straight bar from the B centreline into D, hard corners
    half_bar = [row(4).shifted(-half), row(4).shifted(half)]
    x_d = intersect(c1, row(4))[0]  # vertical bar end inside D, even where D bends
    bar = [intersect(f1, half_bar[0]), (x_d, half_bar[0].q[1]), (x_d, half_bar[1].q[1]), intersect(f1, half_bar[1])]
    return Logo(belts, bar)


def top_y(p: Params) -> float:
    return 8 * p.pitch + p.red / 2


