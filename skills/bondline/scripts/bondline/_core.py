"""bondline engine core: style, fonts, atoms, bonds, molecules, structures from SMILES and molfiles."""
# Part of the bondline engine (see __init__.py): the package's modules share one namespace, so names
# defined in any of them are visible here at call time, as in a single module.
import math, re, os, copy as _copy, itertools, warnings, functools
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import PathPatch, Polygon, Circle
from matplotlib.path import Path as _MPath
from matplotlib.transforms import Bbox as _Bbox
from matplotlib.textpath import TextPath, text_to_path
from matplotlib.font_manager import FontProperties, findfont
matplotlib.rcParams["lines.scale_dashes"] = False

# ------------------------------------------------------------------ geometry
def unit(v):
    n = math.hypot(v[0], v[1]); return v / n if n > 1e-12 else np.zeros(2)
def dirv(deg):
    r = math.radians(deg); return np.array([math.cos(r), math.sin(r)])
def ang(v): return math.degrees(math.atan2(v[1], v[0]))
def rotv(v, deg):
    r = math.radians(deg); c, s = math.cos(r), math.sin(r)
    return np.array([c * v[0] - s * v[1], s * v[0] + c * v[1]])
def perp(v): return np.array([-v[1], v[0]])
def cross(a, b): return a[0] * b[1] - a[1] * b[0]
def adiff(a, b): return (a - b + 180.0) % 360.0 - 180.0
def norm(v): return math.hypot(v[0], v[1])
def box_dist(p, b):
    dx = max(b[0] - p[0], 0.0, p[0] - b[2]); dy = max(b[1] - p[1], 0.0, p[1] - b[3])
    return math.hypot(dx, dy)
def box_gap(a, b):
    """White between two boxes (negative when they overlap)."""
    dx = max(b[0] - a[2], a[0] - b[2]); dy = max(b[1] - a[3], a[1] - b[3])
    return max(dx, dy) if dx < 0 and dy < 0 else math.hypot(max(dx, 0.0), max(dy, 0.0))
def seg_box_gap(p, q, b):
    """Distance from segment pq to box b (0 when the segment enters the box)."""
    p, q = np.asarray(p, float), np.asarray(q, float)
    t0, t1, d = 0.0, 1.0, q - p                                  # clip the segment to the box (Liang-Barsky)
    for k, lo, hi in ((0, b[0], b[2]), (1, b[1], b[3])):
        if abs(d[k]) < 1e-12:
            if p[k] < lo or p[k] > hi: t0, t1 = 1.0, 0.0
            continue
        a0, a1 = (lo - p[k]) / d[k], (hi - p[k]) / d[k]
        t0, t1 = max(t0, min(a0, a1)), min(t1, max(a0, a1))
    if t0 <= t1: return 0.0
    C = np.array([[b[0], b[1]], [b[2], b[1]], [b[2], b[3]], [b[0], b[3]]])
    return min(box_dist(p, b), box_dist(q, b), float(_segs_dist(C, np.array([[*p, *q]])).min()))
def _mix(col, t):
    """A tint of col: t of the colour on white (the fill of a shaded lobe drawn in a colour)."""
    from matplotlib.colors import to_rgb
    return "#%02X%02X%02X" % tuple(round(255 * (1 - t + t * c)) for c in to_rgb(col))
def arc(p0, p1, bend, n=48):
    """Points of a circular arc from p0 to p1 that turns through bend degrees (bulging to the left of p0->p1 for
    bend > 0, as a cubic through the same tangents would): even curvature, no kink in the middle."""
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float); th = math.radians(bend)
    if abs(th) < 1e-6: return np.array([p0 + (p1 - p0) * t for t in np.linspace(0, 1, n)])
    c = p1 - p0; L = norm(c); u = c / L; v = perp(u)
    R = L / (2 * math.sin(abs(th) / 2)); C = (p0 + p1) / 2 - v * math.copysign((L / 2) / math.tan(abs(th) / 2), th)
    a0 = math.atan2(*(p0 - C)[::-1]); sw = -th                   # the tangent turns clockwise for a left bulge
    return np.array([C + R * np.array([math.cos(a0 + sw * t), math.sin(a0 + sw * t)]) for t in np.linspace(0, 1, n)])
def bezier(P, n=64):
    t = np.linspace(0, 1, n)[:, None]
    return (1-t)**3*P[0] + 3*(1-t)**2*t*P[1] + 3*(1-t)*t*t*P[2] + t**3*P[3]
def _segs_dist(pts, S):           # pts (n,2), S (m,4) -> (n,m) distances (component-wise: same arithmetic, half the time)
    ax, ay, bx, by = S[:, 0], S[:, 1], S[:, 2], S[:, 3]; px, py = pts[:, 0:1], pts[:, 1:2]
    dx, dy = bx - ax, by - ay; L2 = np.maximum(dx * dx + dy * dy, 1e-12)
    t = np.clip(((px - ax) * dx + (py - ay) * dy) / L2, 0, 1)
    return np.hypot(ax + t * dx - px, ay + t * dy - py)
def _polys_dist(Q, E, starts):
    """Distances (n, k) from points Q to k convex CCW polygons given as one edge list E (rows x0 y0 x1 y1, grouped
    by polygon; `starts` = first edge of each polygon). 0 inside."""
    if not len(starts): return np.zeros((len(Q), 0))
    D = _segs_dist(Q, E)
    cr = (E[None, :, 2] - E[None, :, 0]) * (Q[:, None, 1] - E[None, :, 1]) - (E[None, :, 3] - E[None, :, 1]) * (Q[:, None, 0] - E[None, :, 0])
    out = np.minimum.reduceat(D, starts, axis=1)
    out[np.add.reduceat((cr < 0).astype(np.int32), starts, axis=1) == 0] = 0.0
    return out
def _seg_x(p, q, r, s):           # proper intersection of segments pq and rs
    d = cross(q - p, s - r)
    if abs(d) < 1e-9: return False
    t = cross(r - p, s - r) / d; u = cross(r - p, q - p) / d
    return 0.02 < t < 0.98 and 0.02 < u < 0.98

def _hull(P):
    P = sorted(set(map(tuple, np.round(P, 4))))
    if len(P) < 3: return np.array(P)
    def half(pts):
        h = []
        for p in pts:
            while len(h) >= 2 and cross(np.subtract(h[-1], h[-2]), np.subtract(p, h[-2])) <= 0: h.pop()
            h.append(p)
        return h
    lo, up = half(P), half(P[::-1]); return np.array(lo[:-1] + up[:-1])
def _thin(H, tol=0.008):
    """Drop hull vertices that lie within tol (em) of the line through their neighbours (curves keep ~16 points)."""
    H = np.asarray(H, float)
    while len(H) > 8:
        a, c = np.roll(H, 1, axis=0), np.roll(H, -1, axis=0); e = c - a
        dev = np.abs(e[:, 0] * (H[:, 1] - a[:, 1]) - e[:, 1] * (H[:, 0] - a[:, 0])) / np.maximum(np.hypot(*e.T), 1e-12)
        i = int(np.argmin(dev))
        if dev[i] >= tol: break
        H = np.delete(H, i, axis=0)
    return H
def _hull_dist(p, H):
    """Distance from p to a convex polygon (0 inside)."""
    e = np.roll(H, -1, axis=0) - H; w = p - H
    if (e[:, 0] * w[:, 1] - e[:, 1] * w[:, 0] >= 0).all(): return 0.0
    t = np.clip((w * e).sum(1) / np.maximum((e * e).sum(1), 1e-12), 0, 1)
    return float(np.hypot(*(H + t[:, None] * e - p).T).min())

# ------------------------------------------------------------------ style & fonts
class Style:
    """One Style per document. Every size derives from `bond` (pt per bond length)."""
    def __init__(self, bond=14.4, font=("Liberation Sans", "Arial", "Helvetica", "DejaVu Sans"),
                 ink="#000000", arrow="#FF2600", hetero=False, charge="circle", preset="default", greek="serif",
                 chemtext=True):
        """hetero=True colours every heteroatom label by element (labels only, never bonds).

        preset="default": ChemDraw's ACS 1996 geometry (14.4 pt bonds, 10 pt Arial, 120 deg) with 1.1 pt lines,
        the weight ChemDraw figures for slides and handouts usually use.
        preset="acs": the strict ACS 1996 sheet (0.6 pt lines) for journal-style output."""
        self.bond = float(bond); k = self.bond / 14.4          # ACS 1996 geometry = 14.4 pt bonds
        self.font, self.ink, self.arrow, self.hetero, self.charge = font, ink, arrow, hetero, charge
        self.greek = greek                # "serif": Greek letters from a Times-style face (ChemDraw's Symbol look); "sans"
        self.chemtext = chemtext          # chemical typography in text: H2O -> H₂O, SN2 -> S_N2, pKa (see _chemtext)
        self.preset = preset
        self.lw = (0.6 if preset == "acs" else 1.1) * k     # bond line width (7.7 % of the bond by default)
        self.bold = 2.0 * k                                 # bold-bond bar
        self.size = 10.0 * k              # atom label size (Arial/Helvetica metrics)
        self.tsize = 10.0 * k             # conditions, captions, notes
        self.dbl = (0.18 if preset == "acs" else 0.19) * self.bond   # double-bond line spacing
        # wedges: envelope over the full bond, narrow end = line width, cut where a label begins
        self.wedge = (2.0 * k) if preset == "acs" else 0.27 * self.bond      # width at the far atom
        self.wedge_tip = self.lw
        self.hash_step = (2.5 * k) if preset == "acs" else 0.206 * self.bond # hash bar spacing
        self.hash_wide = (2.0 * k) if preset == "acs" else 0.30 * self.bond
        self.hash_tip = 0.4 * self.lw
        self.margin = 1.6 * k             # bond-to-label margin
        self.lp_r = 0.0825 * self.size    # lone-pair dot radius (dots 0.165 em)
        self.lp_sep = 0.377 * self.size   # dot centre spacing
        self.lp_gap = 0.15 * self.size    # clear space between the label and the dots
        self.chg_r = 0.33 * self.size     # circled charge radius: the ring is as tall as a capital letter, so it reads as
                                          # a mark on the atom rather than another letter (ChemDraw's own runs ~15 % larger)
        self.alw = self.lw                # curly-arrow line width (same weight as bonds)
        self.head_len = 10.0 * self.lw    # ChemDraw default arrowhead: HeadSize 1000 % of line width
        self.head_w = 2.5 * self.lw       # ArrowheadWidth 250 % (half-width)
        self.head_notch = 0.875           # ArrowheadCenterSize 875 / HeadSize 1000
        self.c_head_len = 0.40 * self.bond   # curly (electron-pushing) arrowhead: ChemDraw proportions (length 4 x half-width)
        self.c_head_w = 0.10 * self.bond     # half-width (full width 0.2 bond)
        self.tip_gap = 1.0 * k               # space between an arrow tip and the label it points at
        self.rxn = 2.6                    # reaction arrow length (BU)
        self.gap = 0.6                    # default item gap in rows (BU)

HETERO = {"O": "#D01C1C", "N": "#1F4FD1", "S": "#B8860B", "P": "#D96A00", "F": "#2E8B57",
          "Cl": "#2E8B57", "Br": "#8B3A1A", "I": "#6A1B9A", "B": "#B07800", "Si": "#6D6D6D"}

LPCH = "\ue000"                   # a lone pair written in text ("B:", ":NEt3"): drawn as the dots used on atoms
VGLYPH = {"→": 1.25, "←": 1.25, "⇌": 1.25, "⊕": 0.80, "⊖": 0.80, LPCH: 0.42}   # drawn as vector shapes at the house weight
                                                                    # (advance in em); ⊕ ⊖ match the circled charges on atoms
_VSPLIT = re.compile("([" + "".join(VGLYPH) + "])")
SIGNS = "+−"                      # drawn as bars no lighter than the charge marks' (Fonts._sign)
_SSPLIT = re.compile("([" + "".join(VGLYPH) + re.escape(SIGNS) + "])")
SUPCH = dict(zip("⁺⁻⁰¹²³⁴⁵⁶⁷⁸⁹ⁿ", "+−0123456789n")); SUBCH = dict(zip("₀₁₂₃₄₅₆₇₈₉₊₋", "0123456789+−"))

NONROUND = set("ILTYJ17ijlrft")    # letters far narrower than their cell: bonds and marks clear their real outline,
                                    # not the cap-height circle (which leaves a gap beside I, L, the t of Et)
GREEK = re.compile("[\u0370-\u03ff\u1f00-\u1fff]+")
GREEK_FONTS = ("Liberation Serif", "Times New Roman", "Times", "Nimbus Roman", "DejaVu Serif")
class Fonts:
    def __init__(self, families, greek="serif"):
        """greek="serif": Greek letters (Δ, δ, ν, α ...) come from a Times-style serif, as ChemDraw's Symbol font and most
        journals set them (a sans hν is hard to tell from hv), scaled to the label font's cap height."""
        self.family = "DejaVu Sans"
        for f in families:
            try:
                findfont(FontProperties(family=f), fallback_to_default=False); self.family = f; break
            except Exception:
                pass
        self.fallback = self.family not in families[:-1] if len(families) > 1 else False
        self._p, self._w = {}, {}
        self.cap = TextPath((0, 0), "H", size=100, prop=self.prop()).get_extents().y1 / 100.0
        self.gfamily = None; self.gscale = 1.0
        if greek == "serif":
            for f in GREEK_FONTS:
                try:
                    findfont(FontProperties(family=f), fallback_to_default=False); self.gfamily = f; break
                except Exception:
                    pass
            if self.gfamily:
                self.gscale = self.cap / (TextPath((0, 0), "H", size=100, prop=self.gprop()).get_extents().y1 / 100.0)
    def gprop(self, it=False, bd=False):
        k = ("g", it, bd)
        if k not in self._p:
            self._p[k] = FontProperties(family=self.gfamily, style="italic" if it else "normal", weight="bold" if bd else "normal", size=100)
        return self._p[k]
    def _gsplit(self, s):
        """s as (text, is_greek) pieces; Greek pieces are drawn from the serif face when there is one."""
        if not self.gfamily or not GREEK.search(s): return [(s, False)]
        out, i = [], 0
        for m in GREEK.finditer(s):
            if m.start() > i: out.append((s[i:m.start()], False))
            out.append((m.group(0), True)); i = m.end()
        if i < len(s): out.append((s[i:], False))
        return out
    def prop(self, it=False, bd=False):
        k = (it, bd)
        if k not in self._p:
            self._p[k] = FontProperties(family=self.family, style="italic" if it else "normal",
                                        weight="bold" if bd else "normal", size=100)
        return self._p[k]
    def adv(self, s, size, it=False, bd=False):
        k = (s, it, bd)
        if k not in self._w:
            parts = self._gsplit(s)
            if len(parts) > 1 or parts[0][1]:
                self._w[k] = sum((text_to_path.get_text_width_height_descent(t, self.gprop(it, bd), ismath=False)[0] / 100.0 * self.gscale)
                                 if g else self.adv(t, 1.0, it, bd) for t, g in parts)
            elif any(ch in VGLYPH for ch in s):
                self._w[k] = sum(VGLYPH[g] if g in VGLYPH else self.adv(g, 1.0, it, bd) for g in _VSPLIT.split(s) if g)
            else:
                with warnings.catch_warnings():               # a missing glyph is reported by check() and raised on saving
                    warnings.simplefilter("ignore"); self._w[k] = text_to_path.get_text_width_height_descent(s, self.prop(it, bd), ismath=False)[0] / 100.0
        return self._w[k] * size
    def missing(self, s):
        """Characters of s that the active font has no glyph for (vector glyphs and spaces excepted)."""
        if not hasattr(self, "_ft"):
            from matplotlib import ft2font
            self._ft = ft2font.FT2Font(findfont(self.prop(), fallback_to_default=True)); self._miss = {}
        if self.gfamily and not hasattr(self, "_gft"):
            from matplotlib import ft2font
            self._gft = ft2font.FT2Font(findfont(self.gprop(), fallback_to_default=True))
        out = []
        for ch in s:
            if ch in VGLYPH or ch == " ": continue
            if ch not in self._miss:
                ft = self._gft if (self.gfamily and GREEK.match(ch)) else self._ft
                self._miss[ch] = ft.get_char_index(ord(ch)) == 0
            if self._miss[ch]: out.append(ch)
        return out
    def _vglyph(self, g, x, y, size):
        """Vector arrow glyphs (→ ← ⇌) in reaction-arrow proportions, centred on the cap height; ⊕ ⊖ as the circled
        charges drawn on atoms (ring 0.33 em, lighter than the sign, centred a little above mid cap height)."""
        from matplotlib.path import Path as _P
        em = size; yc = y + self.cap * size / 2; t2 = 0.0375 * em; hh = 0.34 * em; hw = 0.10 * em; k = 0.875
        x0, x1 = x + 0.125 * em, x + 1.125 * em; polys = []
        if g == LPCH:                            # two lone-pair dots (0.165 em, centres 0.377 em apart) stacked like a
            r = 0.0825 * em; cx = x + 0.21 * em; verts, codes = [], []           # colon, 0.13 em clear of the letter
            th = np.linspace(0, 2 * math.pi, 37)[:-1]
            for cy in (yc - 0.1885 * em, yc + 0.1885 * em):
                P = [(cx + r * math.cos(t), cy + r * math.sin(t)) for t in th]
                verts += P + [P[0]]; codes += [_P.MOVETO] + [_P.LINETO] * (len(P) - 1) + [_P.CLOSEPOLY]
            return _P(np.array(verts, float), codes)
        if g in ("⊕", "⊖"):
            c = np.array([x + 0.40 * em, y + 0.62 * self.cap * size]); r = 0.33 * em; tr, tm = 0.0693 * em, 0.0792 * em
            th = np.linspace(0, 2 * math.pi, 73)[:-1]
            outer = c + np.column_stack([np.cos(th), np.sin(th)]) * (r + tr / 2)
            inner = (c + np.column_stack([np.cos(th), np.sin(th)]) * (r - tr / 2))[::-1]
            s = 0.58 * r + tm / 2
            bars = [[(c[0] - s, c[1] - tm / 2), (c[0] + s, c[1] - tm / 2), (c[0] + s, c[1] + tm / 2), (c[0] - s, c[1] + tm / 2)]]
            if g == "⊕": bars.append([(c[0] - tm / 2, c[1] - s), (c[0] + tm / 2, c[1] - s), (c[0] + tm / 2, c[1] + s), (c[0] - tm / 2, c[1] + s)])
            verts, codes = [], []
            for P in [outer, inner] + bars:
                P = [tuple(p) for p in P]; verts += P + [P[0]]; codes += [_P.MOVETO] + [_P.LINETO] * (len(P) - 1) + [_P.CLOSEPOLY]
            return _P(np.array(verts, float), codes)
        def full(xa, xb, yy):                    # tail xa, tip xb (either direction)
            sg = 1 if xb > xa else -1; xn = xb - sg * k * hh; xs = xn - sg * (t2 / hw) * (1 - k) * hh
            return [(xa, yy - t2), (xs, yy - t2), (xb - sg * hh, yy - hw), (xb, yy), (xb - sg * hh, yy + hw), (xs, yy + t2), (xa, yy + t2)]
        def harpoon(xa, xb, yy, up):             # half head: one barb, above (up) or below the shaft
            sg = 1 if xb > xa else -1; sy = 1 if up else -1
            return [(xa, yy - sy * t2), (xb, yy - sy * t2), (xb - sg * hh, yy + sy * hw), (xb - sg * k * hh, yy + sy * t2), (xa, yy + sy * t2)]
        if g == "→": polys = [full(x0, x1, yc)]
        elif g == "←": polys = [full(x1, x0, yc)]
        else: polys = [harpoon(x0, x1, yc + 0.09 * em, True), harpoon(x1, x0, yc - 0.09 * em, False)]
        verts, codes = [], []
        for P in polys:
            verts += list(P) + [P[0]]; codes += [_P.MOVETO] + [_P.LINETO] * (len(P) - 1) + [_P.CLOSEPOLY]
        return _P(np.array(verts, float), codes)
    sign_w = 0.0                         # pt: the least bar weight of a + or − in text (Scene sets it to its charge marks')
    def _sign(self, g, x, y, size, it=False, bd=False, op=False):
        """A + or − in text as two (one) flat bars where the font draws it, but never lighter than sign_w: at the
        superscript size the font's own sign is half a letter's stem (H⁺, Na⁺ read faint), so every sign in a figure
        carries the weight of the signs in its circled charges, at any size. (An operator, op, keeps the font's weight.)"""
        from matplotlib.path import Path as _P
        k = ("sign", g, it, bd)
        if k not in self._w:
            e = TextPath((0, 0), g, size=100, prop=self.prop(it, bd)).get_extents()
            em = TextPath((0, 0), "−", size=100, prop=self.prop(it, bd)).get_extents()          # the font's bar weight
            self._w[k] = ((e.x0 + e.x1) / 200, (e.y0 + e.y1) / 200, (e.x1 - e.x0) / 200, (em.y1 - em.y0) / 100)
        cx, cy, h, t0 = self._w[k]; c = (x + cx * size, y + cy * size); h *= size; t = max(t0 * size, 0 if op else self.sign_w) / 2
        bars = [[(c[0] - h, c[1] - t), (c[0] + h, c[1] - t), (c[0] + h, c[1] + t), (c[0] - h, c[1] + t)]]
        if g == "+": bars.append([(c[0] - t, c[1] - h), (c[0] + t, c[1] - h), (c[0] + t, c[1] + h), (c[0] - t, c[1] + h)])
        verts, codes = [], []
        for P in bars: verts += P + [P[0]]; codes += [_P.MOVETO] + [_P.LINETO] * 3 + [_P.CLOSEPOLY]
        return _P(np.array(verts, float), codes)
    def path(self, s, x, y, size, it=False, bd=False):
        miss = self.missing(s)
        if miss:
            raise ValueError("the house font has no glyph for " + ", ".join(f"'{c}' (U+{ord(c):04X})" for c in dict.fromkeys(miss))
                             + f" in {s!r}: write it another way (Δ, hν, °C, −, ‡, δ+ and → are available; use S.arrow(kind='<=>') for equilibria)")
        gp = self._gsplit(s)
        if any(ch in VGLYPH or ch in SIGNS for ch in s) or len(gp) > 1 or gp[0][1]:
            from matplotlib.path import Path as _P
            parts, cx = [], x
            for t, g in gp:
                vs = _SSPLIT.split(t) if not g else [t]
                for i, v in enumerate(vs):
                    if not v: continue
                    if v in VGLYPH: parts.append(self._vglyph(v, cx, y, size))
                    elif g: parts.append(TextPath((cx, y), v, size=size * self.gscale, prop=self.gprop(it, bd)))
                    elif v in SIGNS: parts.append(self._sign(v, cx, y, size, it, bd, op="".join(vs[i + 1:])[:1].isspace()))
                    elif v.strip(): parts.append(TextPath((cx, y), v, size=size, prop=self.prop(it, bd)))
                    cx += self.adv(v, size, it, bd)
            return _P.make_compound_path(*parts)
        return TextPath((x, y), s, size=size, prop=self.prop(it, bd))
    def hull(self, ch, it=False, bd=False):
        """Convex hull of one glyph's outline at size 1 (array of vertices, CCW), or None for blank glyphs."""
        k = (ch, it, bd)
        if not hasattr(self, "_h"): self._h = {}
        if ch.isspace(): return None
        if k not in self._h:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                polys = (self._vglyph(ch, 0, 0, 100) if ch in VGLYPH else
                         TextPath((0, 0), ch, size=100 * self.gscale, prop=self.gprop(it, bd)) if (self.gfamily and GREEK.match(ch)) else
                         TextPath((0, 0), ch, size=100, prop=self.prop(it, bd))).to_polygons(closed_only=False)
            v = np.vstack(polys) if polys else np.zeros((0, 2))       # points on the outline itself (curves flattened),
            self._h[k] = _thin(_hull(v / 100.0)) if len(v) > 2 else None   # not the Bezier control points around it
        return self._h[k]

# ------------------------------------------------------------------ rich text and labels
SUBF = 0.70
SUPY = 0.50      # superscript raise (em): its foot clears a bond drawn to the letter's middle (R¹, Pd⁰) by ~0.9 pt
class _CStr(str):
    """The text of a run set in its own colour (<c=..>): the colour rides on the string, which drawing reads; slices and
    joins (measuring, checking) drop it."""
def parse_rich(s, _sz=1.0, _y=0.0, _it=False, _bd=False, _c=None):
    """'<i>..</i>', '<b>..</b>', '<c=#E8261A>..</c>' colour, '_{..}' / '_x' subscript, '^{..}' / '^x' superscript; they
    nest (10^{pK_a(H_3O^+)})."""
    runs, it, bd, cs, i, buf = [], _it, _bd, [_c], 0, ""
    def flush():
        nonlocal buf
        if buf:
            if cs[-1]: buf = _CStr(buf); buf.c = cs[-1]
            runs.append((buf, _sz, _y, it, bd)); buf = ""
    while i < len(s):
        mc = re.match(r"<c=([^<>\s]+)>|</c>", s[i:])
        if mc:
            flush(); i += mc.end()
            if mc.group(1): cs.append(mc.group(1))
            elif len(cs) > 1: cs.pop()
            continue
        for tag, val in (("<i>", ("it", True)), ("</i>", ("it", False)), ("<b>", ("bd", True)), ("</b>", ("bd", False))):
            if s.startswith(tag, i):
                flush(); i += len(tag)
                if val[0] == "it": it = val[1]
                else: bd = val[1]
                break
        else:
            c = s[i]
            if c in "_^" and i + 1 < len(s):
                flush()
                if s[i + 1] == "{":
                    d, j = 0, i + 1
                    while j < len(s):
                        d += (s[j] == "{") - (s[j] == "}")
                        if d == 0: break
                        j += 1
                    txt = s[i + 2:j]; i = j + 1
                else: txt = s[i + 1]; i += 2
                runs += parse_rich(txt, _sz * SUBF, _y + (-0.25 if c == "_" else SUPY) * _sz, it, bd, cs[-1]); continue
            buf += c; i += 1
    flush(); return runs

ELEMENTS = set(("D H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb "
                 "Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W "
                 "Re Os Ir Pt Au Hg Tl Pb Bi Th U").split())
TEXT_UNITS = sorted(set(("TBDMS TBDPS TIPS TBS TMS TES Fmoc Boc Cbz Troc PMP PMB DMAP Mes Tol Piv Ph Me Et Pr Bu Bn Bz Ac "
                         "Ts Ms Tf Ns Nf Cy Py Ad Np Tr Xyl Hex Pent Hept Oct Dec Cp Cp* IMes SIMes IPr SIPr BINAP "
                         "dppf dppe dppp dppb dmpe acac dba cod COD nbd bpy bipy phen ppy dtbbpy dtbpy bpz dap pin cat "
                         "dme glyme py tmeda TMEDA en terpy tpy dF OAc OTf OTs OMs").split()), key=len, reverse=True)
TEXT_GENERIC = ("R", "Ar", "X", "Y", "Z", "L", "M", "E", "Nu", "Hal", "LG", "PG")
def _formula_markup(tok):
    """A word written as a formula (H2O, CH2Cl2, Pd(PPh3)4, Et3N, iPr2NEt, H3O+, PF6-, Fe2+, R1, Pd0, [Ru(bpy)3]2+)
    in rich-text markup: counts subscript, numbered groups (R¹) and a metal's zero oxidation state superscript, charges
    superscript, alkyl prefixes italic superscripts (ᵗBuOK) as in labels. None if the word is not a formula or needs
    no change (THF, PCC, C1, E2, a number)."""
    units = _formula_units(tok)
    if units is None: return None
    if not any(k in ("d", "q", "p", "ph") for k, _ in units): return None
    atoms = [k for k, _ in units if k in ("u", "g")]
    if not atoms: return None
    return _formula_out(units, atoms)
def _formula_units(tok):
    """A word parsed as a formula: [(kind, text)] with kind 'u' unit, 'g' generic group, 'p' alkyl prefix, 'o'/'c'
    bracket, 'd' digits, 'q' closing charge, "'" prime; None if it is not one."""
    if not tok or not tok[0].isalpha() and tok[0] not in "([": return None
    units, i, n = [], 0, len(tok)
    while i < n:
        c = tok[i]
        mo = re.match(r"\((?:0|I{1,3}|IV|VI{0,3}|V)\)", tok[i:])
        if mo and units and units[-1][0] == "u" and units[-1][1] in METALS_T:   # Pd(0), Pd(II): as written
            units.append(("x", mo.group(0))); i += mo.end(); continue
        mo = re.match(r"([omp])-(tolyl|anisyl|xylyl|Tol|An)(?![a-z])", tok[i:])   # P(o-tolyl)3
        if mo: units.append(("x", "<i>" + mo.group(1) + "</i>-" + mo.group(2))); i += mo.end(); continue
        if c == "." and units and tok[i + 1:i + 2].isalnum():   # an adduct or salt written with a full stop: Et3N.HCl
            units.append(("x", ".")); i += 1; continue
        if c in "([": units.append(("o", c)); i += 1; continue
        if c in ")]": units.append(("c", c)); i += 1; continue
        if c in "'": units.append(("'", c)); i += 1; continue
        if c in "0123456789":                                   # (ASCII digits only: H₂O typed with its own
            j = i                                               # subscript is kept as typed)
            while j < n and tok[j] in "0123456789": j += 1
            units.append(("d", tok[i:j])); i = j; continue
        if c in "+-−" and i == n - 1: units.append(("q", "−" if c in "-−" else "+")); i += 1; continue
        hit = None
        for p in ("tert-", "sec-", "t-", "i-", "n-", "s-", "t", "i", "n", "s"):   # an alkyl prefix only before an alkyl
            if tok.startswith(p, i):                                              # group (NaOt-Bu keeps its hyphen)
                for a in ("Bu", "Pr", "Pent", "Hex", "Am"):
                    if tok.startswith(a, i + len(p)): hit = ("ph" if p.endswith("-") else "p", p.rstrip("-")); break
            if hit: break
        if hit: units.append(hit); i += len(hit[1]) + (1 if hit[0] == "ph" else 0); continue
        for a in TEXT_UNITS:
            if tok.startswith(a, i): hit = ("u", a); break
        if not hit:
            for g in sorted(TEXT_GENERIC, key=len, reverse=True):
                if tok.startswith(g, i) and not (g in ("E", "M", "L", "X", "Y", "Z", "R") and tok[i:i + 2] in ELEMENTS and len(g) == 1 and tok[i:i + 2] != g):
                    hit = ("g", g); break
        if not hit and tok[i:i + 2] in ELEMENTS and tok[i + 1:i + 2].islower(): hit = ("u", tok[i:i + 2])
        if not hit and tok[i] in ELEMENTS: hit = ("u", tok[i])
        if not hit and tok[i] == "n" and units and units[-1][0] in ("c", "g") and units[-1][1] in (")", "]", "L"):
            hit = ("d", "n")                                       # (CH2)n, Ln: a count
        if not hit: return None
        units.append(hit); i += len(hit[1])
    return units
def _formula_out(units, atoms):
    if len(atoms) == 1 and units[0][0] == "u" and units[0][1] == "C" and len(units) == 2 and units[1][0] == "d": return None   # C1, C2: locants
    numbered = sum(1 for j, (k, v) in enumerate(units[:-1]) if k == "g" and units[j + 1][0] == "d") > 1
    out, j = [], 0
    while j < len(units):
        k, v = units[j]
        if k == "p": out.append("<i>^{" + v + "}</i>")
        elif k == "ph": out.append("<i>" + v + "</i>-")
        elif k == "x": out.append(v)
        elif k == "d":
            pk, pv = units[j - 1] if j else (None, None)
            if pk is None or pk == "o": return None
            if pk == "x" and pv == ".": out.append(v); j += 1; continue   # MgSO4.7H2O: a coefficient
            last = j == len(units) - 1 or (j == len(units) - 2 and units[j + 1][0] == "q")
            if j == len(units) - 2 and units[j + 1][0] == "q":    # digits before a closing charge: Fe2+, SO42−, [..]2+, NH4+
                if pk == "c" and pv == "]" or len(atoms) == 1 and pk == "u" and len(units) == 3 and \
                        not (pv in ("N", "I", "Br") and v == "3"):         # N3-, I3-: azide, triiodide (N_3^-)
                    out.append("^{" + v + units[j + 1][1] + "}"); j += 2; continue
                if len(v) == 2 and pk in ("u", "c", "g"):
                    out.append("_{" + v[0] + "}^{" + v[1] + units[j + 1][1] + "}"); j += 2; continue
            if pk == "g" and v != "n" and (numbered or len(units) == 2 or v == "1"): out.append("^{" + v + "}")
            elif v == "0" and pk == "u" and pv in METALS_T: out.append("^{0}")
            elif v == "1" and pk != "c": return None                # nobody writes a count of one: not a formula (C1-C2)
            else: out.append("_{" + v + "}")
        elif k == "q": out.append("^{" + v + "}")
        else: out.append(v)
        j += 1
    return "".join(out)
METALS_T = set("Li Na K Cs Mg Ca B Al Sc Ti Zr Hf V Nb Ta Cr Mo W Mn Re Fe Ru Os Co Rh Ir Ni Pd Pt Cu Ag Au Zn Cd Hg Sn Pb Bi La Ce Sm Yb".split())
_SN = re.compile(r"(?<![A-Za-z])S([NE])(1|2|i|Ar)(?=$|[^A-Za-z0-9])")
_WORD = re.compile(r"(?:<[^<>]*>|[^\s,;/:·=<])+")     # a tag (<c=#..>) stays inside its word
def _chemtext(s):
    """Chemical typography for running text (conditions, captions, notes): formulas typeset as in labels (H₂O, CH₂Cl₂,
    Pd(PPh₃)₄, Et₃N, H₃O⁺, Fe²⁺, ᵗBuOK, R¹, Pd⁰), S_N2 / S_EAr, E1_cB, p<i>K</i>_a, italic stereodescriptors ((R)-, (E)-),
    italic locant and isomer prefixes (O-silylation, N-oxide, p-TsOH, m-CPBA, cis-, trans-, syn-, anti-, fac-) and unit
    exponents (mol⁻¹, cm⁻¹). Words already holding markup (_ ^ <) are left alone; a backslash before a word keeps it as
    typed (\\C2)."""
    if not s or not re.search(r"[A-Za-z]", s): return s
    def symbolic(w):                                       # B, Nu, Base, X, NEt3, R3N, HO ...: a species, not a word of prose
        return w in ("B", "Base", "Nu", "E", "X", "Y", "Z", "A", "M", "L", "R", "Ar") or (len(w) <= 6 and _formula_units(w) is not None)
    s = re.sub(r"(?<![\w:\\])([A-Za-z][A-Za-z0-9]{0,5}):(?=$|[⊖⊕\-−+,;)^]|\s+(?![A-Za-z]{3,}))",
               lambda q: q.group(1) + LPCH if symbolic(q.group(1)) else q.group(0), s)
    s = re.sub(r"(?:(?<=^)|(?<=[\s(])):(?=([A-Za-z][A-Za-z0-9]{0,5})(?:$|[\s⊖⊕\-−+,;)^]))",
               lambda q: LPCH if symbolic(q.group(1)) else q.group(0), s)
    def word(m):
        w = m.group(0); W = lambda x: word(re.match(r".*", x))
        if w.startswith("\\"): return w[1:]
        mt = re.fullmatch(r"((?:<[^<>]*>)*)(.*?)((?:<[^<>]*>)*)", w)      # <i>−CO2</i>, <b>Ru(bpy)3Cl2</b>: the
        if mt.group(1) or mt.group(3):                                     # word inside its tags
            return mt.group(1) + (W(mt.group(2)) if mt.group(2) and "<" not in mt.group(2) else mt.group(2)) + mt.group(3)
        if "<" in w or ">" in w: return w
        if _MARK.search(w):                                                # H3O^+: the plain parts typeset, markup kept;
            ps = re.split("(" + _MARK.pattern + ")", w)                    # a part after markup goes on the same word
            return "".join(p if not p or _MARK.fullmatch(p) else (W(p) if k == 0 or p[0] not in "-–" else p[0] + W(p[1:]))
                           for k, p in enumerate(ps))
        if re.search(r"[{}]", w): return w
        if len(w) > 1 and w[0] in (LPCH, "•"): return w[0] + W(w[1:])
        if len(w) > 1 and w[-1] in (LPCH, "•"): return W(w[:-1]) + w[-1]   # B:, R2CH• (the dot as typed)
        if len(w) > 1 and w[-1] in ("⊖", "⊕"): return W(w[:-1]) + "^{" + w[-1] + "}"   # H3O⊕: a charge, raised
        mu = re.fullmatch(r"([A-Za-z0-9()\[\]]*[A-Za-z0-9)\]])([⁰¹²³⁴⁵⁶⁷⁸⁹]*[⁺⁻])", w)   # NMe3⁺, N3⁻: unicode charge signs
        if mu:
            f = _formula_markup(mu.group(1) + "".join(SUPCH[c] for c in mu.group(2)))
            if f is not None: return f
        mx = re.match(r"^(\d+(?:\.\d+)?)[×x](?=[A-Z(\[])", w)     # 2×LDA, 2xMe3SiCl: a multiple of a reagent
        if mx: return mx.group(1) + "×" + word(re.match(r".*", w[mx.end():]))
        sg = re.match(r"^([-−+])(?=[A-Z])", w)                    # -H2O, +H2O, -CO2: a species lost or gained
        if sg: return ("−" if sg.group(1) in "-−" else "+") + word(re.match(r".*", w[1:]))
        if re.fullmatch(r"E[12]'?", w): return w                  # E1 and E2 eliminations
        w = _SN.sub(lambda q: "S_{" + q.group(1) + "}" + q.group(2), w)
        if "_{" in w: return w
        w = re.sub(r"^E1c([Bb])$", r"E1_{c\1}", w)
        if "_{" in w: return w
        w = re.sub(r"^p(K)(a|aH|b)('?)$", lambda q: "p<i>K</i>_{" + q.group(2) + "}" + ("′" if q.group(3) else ""), w)
        if "<" in w: return w
        w = re.sub(r"^(mol|cm|s|h|min|M|L|g|kg|K|nm|m)[-−](\d)$", lambda q: q.group(1) + "^{−" + q.group(2) + "}", w)
        if "^" in w: return w
        w = re.sub(r"\(([RSEZ](?:,[RSEZ])*)\)", lambda q: "(" + ",".join("<i>" + x + "</i>" for x in q.group(1).split(",")) + ")", w)
        w = re.sub(r"^(?:(cis|trans|syn|anti|endo|exo|meso|rac|fac|mer|ortho|meta|para|ipso)-(?=[A-Za-z(0-9])"
                   r"|([ompNOSCP])-(?=[A-Za-z(])|([tnsi])-(?=Bu|Pr|Am|Pent|Hex))",
                   lambda q: "<i>" + (q.group(1) or q.group(2) or q.group(3)) + "</i>-", w)
        if w.startswith("<i>") and re.match(r"<i>[NOSCP]</i>-[A-Z]", w): w = re.sub(r"^<i>([NOSCP])</i>", r"\1", w)   # C-H, N-H: bonds
        head = re.match(r"^((?:<i>[^<]*</i>-)?)(.*)$", w); pre, body = head.group(1), head.group(2)
        tail = ""
        if re.search(r"[A-Za-z0-9)\]][-−]$", body) and _formula_markup(body[:-1] + "−") is not None:
            return pre + _formula_markup(body[:-1] + "−")         # PF6-, BF4-: a trailing hyphen is a charge
        f = _formula_markup(body) if "-" in body else None          # NaOt-Bu, P(o-tolyl)3: a formula holding a
        if f is not None: return pre + f                             # hyphen, tried whole before splitting there
        parts = re.split(r"([-–])", body); out = []
        for p in parts:
            if p in ("-", "–") or not p: out.append(p); continue
            f = _formula_markup(p); out.append(f if f is not None else p)
        return pre + "".join(out)
    return _WORD.sub(word, s)

def _typo(s, F=None, minus=True, chem=False):
    """Typographic clean-up of text: unicode super/subscripts become ^{..}/_{..} runs (the house font lacks several of
    them, and runs keep them all alike), ℃ becomes °C, and a hyphen before a number (−78 °C) becomes a true minus.
    chem=True also applies the chemical typography of _chemtext (running text, not labels)."""
    s = s.replace("℃", "°C")
    if chem: s = _chemtext(re.sub(r"(?:(?<=^)|(?<=[\s(]))[–-](?=\d)", "−", s))
    def conv(tab, mark):
        def f(m):
            g = m.group(0)
            if F is not None and not F.missing(g): return g        # labels: only where the font cannot draw them
            return mark + "{" + "".join(tab[c] for c in g) + "}"
        return f
    s = re.sub("[" + "".join(SUPCH) + "]+", conv(SUPCH, "^"), s)
    s = re.sub("[" + "".join(SUBCH) + "]+", conv(SUBCH, "_"), s)
    if minus: s = re.sub(r"(?:(?<=^)|(?<=[\s(]))-(?=\d)", "−", s)
    return s

ABBR = sorted(("TBDMS TBDPS TIPS TBS TMS TES Fmoc Boc Cbz Troc Alloc Teoc PMP PMB DMB DMAP Base Hal Mes Tol Piv Ph "
               "Me Et Pr Bu Pent Hex Hept Oct Bn Bz Ac Ar Het Ts Ms Tf Cy Py Nu LG PG FG DG EWG EDG THP MOM MEM SEM "
               "BOM DBS PNB POM LA LB Ad Np Tr Xyl dppf dppe dppp acac Cp").split(),
              key=len, reverse=True)
ELEM2 = set(("He Li Be Ne Na Mg Al Si Cl Ca Sc Ti Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Zr Nb "
             "Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te Xe Cs Ba La Ce Nd Sm Eu Gd Yb Lu Hf Ta Re Os Ir Pt Au "
             "Hg Tl Pb Bi Po At Rn Th").split())
PREFIX = ("tert", "sec", "iso", "t", "i", "n", "s", "c")
ALKYL = ("Bu", "Pr", "Am", "Pent", "Hex", "Hept", "Oct")
LOWER_UNITS = ("pin", "cat", "dan", "neop", "mida", "dba", "cod", "nbd", "acac", "bpy", "phen", "py", "tmeda", "en")
_MARK = re.compile(r"[\^_](?:\{(?:[^{}]|\{[^{}]*\})*\}|[^{])")     # one ^{..} / _{..} (one level of nesting) or ^x / _x

def _alkyl_at(lbl, i):
    """The prefixed alkyl group starting at i (tBu, iPr, nBu, cHex, nC6H13), or None: a prefix only before an alkyl."""
    for p in PREFIX:
        if lbl.startswith(p, i):
            k = i + len(p)
            for a in ALKYL:
                if lbl.startswith(a, k): return p + a
            if lbl.startswith("C", k) and lbl[k + 1:k + 2].isdigit(): return p + "C"
    return None

def _group(lbl, i):
    for a in ABBR:
        if lbl.startswith(a, i): return a
    if lbl[i].isupper() and lbl[i:i + 2] in ELEM2:
        m = re.match(r"([NPBOS])([tinsc](?:%s))(\d*)" % "|".join(ALKYL), lbl[i:])   # NiPr2, PtBu3, OiPr: N(iPr)2, P(tBu)3
        if not (m and int(m.group(3) or 1) <= (3 if m.group(1) in "NPB" else 1)):   # and O-iPr, but SnBu3, SiPr3 and
            return lbl[i:i + 2]                                                      # NiCl2 are tin, silicon, nickel
    return lbl[i]

def _label_style(lbl):
    """A label wrapped whole in <b>..</b> and/or <i>..</i> (a bold [M]): (the label inside, italic, bold)."""
    it = bd = False
    while lbl:
        m = re.fullmatch(r"<([bi])>(.*)</\1>", lbl)
        if not m or f"</{m.group(1)}>" in m.group(2): break
        it, bd, lbl = it or m.group(1) == "i", bd or m.group(1) == "b", m.group(2)
    return lbl, it, bd

def _ucore(u):
    """A label unit without its markup, count and primes: Pd^{II} -> Pd, ^{13}C -> C, L_n -> L, Ln -> L, Me2 -> Me."""
    u = _MARK.sub("", u)
    return re.sub(r"(?<=[L)\]])n$", "", re.sub(r"[\d']+$", "", u))

WHOLE = {"C6H5", "C6F5"}                                    # formula labels read as one block either way (never H5C6)
def label_units(lbl): return list(_label_units(lbl))
@functools.lru_cache(maxsize=4096)
def _label_units(lbl):
    """'CO2Me' -> ['C','O2','Me'] (the first unit is the atom the bonds attach to). A unit is an element, an
    abbreviation (OTBS -> O, TBS; PG, FG, EWG are one word), a prefixed alkyl (tBu), a bracketed group ((OMe)2, [M]),
    or a lowercase word: its own unit after an element (Bpin, B2pin2), else part of the word before it (Nuc, Phth).
    Counts and primes stay with their unit (Me2, R'), as does markup after it (Pd^{II}, L_n) or before the first one
    (^{13}C); a count n after L or a bracket is a count (Ln, (CH2)n). A whole-label <b>/<i> is set aside."""
    lbl = _label_style(lbl)[0]
    if lbl in WHOLE or lbl in ACRONYMS or "<" in lbl or lbl.startswith("\\"): return (lbl,)
    fa = re.fullmatch(r"C(\d+)([HF])(\d+)", lbl)                # C6H13, C5H11, C6F13: an alkyl formula reads as written
    if fa and int(fa.group(1)) >= 2 and int(fa.group(3)) == 2 * int(fa.group(1)) + 1: return (lbl,)   # (CH3 flips: H3C)
    units, i, n, pre = [], 0, len(lbl), ""
    while i < n:
        mk = _MARK.match(lbl, i)
        if mk:                                             # markup belongs to the unit before it (an isotope's to the next)
            if units: units[-1] += mk.group(0)
            else: pre += mk.group(0)
            i = mk.end(); continue
        c = lbl[i]
        if c in "([":
            d, j = 0, i
            while j < n:
                d += (lbl[j] == c) - (lbl[j] == {"(": ")", "[": "]"}[c])
                if d == 0: break
                j += 1
            u = lbl[i:j + 1]
        else:
            u = _alkyl_at(lbl, i)
            if u is None and c.islower():                  # a lowercase word: Bpin, Nuc, cat; or the n of Ln
                j = i + 1
                while j < n and lbl[j].islower() and _alkyl_at(lbl, j) is None: j += 1
                run, prev = lbl[i:j], (units[-1] if units else None)
                if prev is not None and not _MARK.search(prev):
                    core = _ucore(prev)
                    if (run == "n" and (core == "L" or prev[-1:] in ")]")) or \
                            not (run in LOWER_UNITS and (core in ELEM2 or core in ELEMENTS)):
                        units[-1] = prev + run; i = j          # Ln, (CH2)n; Nuc, Phth: one word; Bpin, Bcat: their own unit
                        m = re.match(r"\d+'*|'+", lbl[i:])
                        if m: units[-1] += m.group(0); i += len(m.group(0))
                        continue
                u = run
            u = u or _group(lbl, i)
        i += len(u)
        if u == "C" and units and units[-1][:1] == "(" and re.match(r"H\d$", lbl[i:]):   # (CH2)11CH3: its CH3 end is
            u += lbl[i:]; i = n                                                              # one unit
        m = re.match(r"\d+'*|'+", lbl[i:])
        if m: u += m.group(0); i += len(m.group(0))
        elif u in METALS_T:                                # PdII, NiIII, TiIV: an oxidation state (one I or V is an
            m = re.match(r"(?:VIII|VII|VI|IV|III|II)(?![a-z])", lbl[i:])   # element: CuI, AgI)
            if m: u += "^{" + m.group(0) + "}"; i += len(m.group(0))
        units.append(pre + u); pre = ""
    if pre: units.append(pre)
    if any(_ucore(u) in ("G", "J", "Q") for u in units): return (lbl,)  # not a formula: an acronym, read as one word
    return tuple(units)

VALENCE = {"O": 2, "S": 2, "Se": 2, "Te": 2, "N": 3, "P": 3, "As": 3, "Sb": 3, "B": 3, "Al": 3, "Ga": 3,
           "C": 4, "Si": 4, "Ge": 4, "Sn": 4, "Pb": 4}
def r_ident(units, k, bonds=1, onium=False):
    """True if the digits of R-group unit k (R1, Ar2) in a condensed label number the group (OR¹, CO₂R², NR¹R², CH₂Ar¹)
    rather than count it (NR₂, SiR₃, PAr₃, CR₃): a count needs that many free bonds on the atom the group hangs from.
    bonds: bond orders at the label's atom (0 for a lone species). onium: the atom is a cation or bonded to a metal, so
    N, P, O or S there has one bond more (NR₃⁺, Ar₃P→M)."""
    m = re.fullmatch(r"(R|Ar)(\d+)", units[k])
    if not m or k == 0: return False                       # a lone R1, X1, Ar2: unit_runs(first=True) decides
    n = int(m.group(2))
    if n <= 1: return True                                 # nobody counts one group
    if sum(bool(re.fullmatch(r"(?:R|Ar)\d*'*", u)) for u in units) > 1: return True     # NR1R2: two different groups
    used, j = 0, k - 1
    while j > 0:
        h = re.fullmatch(r"([A-Z][a-z]?)(\d*)", units[j])
        if h and h.group(1) != "H" and h.group(1) in VALENCE: break
        c = re.search(r"(\d+)'*$", units[j]); used += int(c.group(1)) if c else 1      # H2, Me2, (OMe): bonds taken
        j -= 1
    host = re.match(r"[A-Z][a-z]?", units[j]); el = host.group(0) if host else ""
    if el not in VALENCE: return False                     # a metal or a group: a count (PdL2, Ph3)
    used += bonds if j == 0 else 1                         # the bonds into the host: the skeleton's, or the chain's
    free = VALENCE[el] - used + (1 if (onium and j == 0 and el in ("N", "P", "As", "Sb", "O", "S", "Se")) else 0)
    return n > free

def unit_runs(u, first=True, ident=False):
    """Runs of one label unit. Digits are subscripts (counts) except: the number of a lone R-group label (R1, X2, Ar1),
    an R-group number inside a condensed label (ident=True, see r_ident), and an oxidation state of zero (Pd0, Ni0):
    those are superscripts. The n of Ln or (CH2)n is a count too; markup in the unit (Pd^{II}, L_n) is set as written."""
    runs, prev, i = [], "", 0
    for seg in _MARK.split(u) if not _MARK.search(u) else re.split("(" + _MARK.pattern + ")", u):
        if not seg: continue
        if seg[0] in "^_" and _MARK.fullmatch(seg): runs += parse_rich(seg); prev = seg; i += 1; continue
        for tok in re.findall(r"\d+|L(?=n(?![a-z]))|(?<=[L)\]])n(?![a-z])|[A-Z][a-z]*|[a-z]+|.", seg):
            if tok.isdigit():
                sup = (i == 1 and (ident or (first and prev in ("R", "X", "Z", "Ar", "L", "E")))) or (tok == "0" and prev[:1].isalpha())
                runs.append((tok, SUBF, SUPY if sup else -0.25, False, False))
            elif tok == "n" and prev[-1:] in ("L", ")", "]"): runs.append((tok, SUBF, -0.25, False, False))   # Lₙ
            elif tok in PREFIX: runs.append((tok, SUBF, SUPY, True, False))
            else: runs.append((tok, 1.0, 0.0, False, False))
            prev = tok; i += 1
    return runs

# ------------------------------------------------------------------ molecules
_UID = itertools.count(1)

class LP:
    """A lone pair on an atom; angle in degrees (None = placed automatically)."""
    def __init__(self, atom, angle=None):
        self.atom, self.angle, self.auto, self._pull = atom, angle, angle is None, None
    def __repr__(self): return f"LP({self.atom!r}, {self.angle})"

class Mid:
    """A point between two atoms (where a new bond forms): Mid(a, b, t=0.5)."""
    def __init__(self, a, b, t=0.5): self.a, self.b, self.t = a, b, t
    def __repr__(self): return f"Mid({self.a!r}, {self.b!r})"

HETMET = set("O N S P Se Te Si Ge Sn B As Sb Bi Mg Li Na K Cs Be Ca Zn Cd Hg Cu Ag Au Al Ga In Tl Ti Zr Pb".split())
METALS2 = set("Mg Li Na K Cs Be Ca Zn Cd Hg Cu Ag Au Al Ga In Tl Ti Zr Pb Sn B Si".split())
METALS = set("Li Na K Cs Mg Ca B Al Ga In Sc Ti Zr Hf V Nb Ta Cr Mo W Mn Re Fe Ru Os Co Rh Ir Ni Pd Pt Cu Ag Au Zn Cd Hg "
             "Sn Pb Bi La Ce Sm Yb M".split())         # acceptors a lone pair can be dative to (Et₃N→BH₃, Ar₃P→Pd)
def anchor_first(lbl):
    """A label typed the way it reads with its bond coming from the right (HO, H2N, H3C, MeO, EtO, AcO, TsO, Me2N,
    BocHN, BrMg, F3C ...) rewritten anchor-first (OH, NH2, CH3, OMe, OEt, OAc, OTs, NMe2, NHBoc, MgBr, CF3), so bonds,
    charges and lone pairs attach to the bonded atom. Returns (label, was_reversed). Rich labels are left alone."""
    if not lbl or re.search(r"[\^<]", lbl): return lbl, False
    plain = re.sub(r"_\{?(\d+)\}?", r"\1", lbl) if re.fullmatch(r"(?:[A-Za-z()']+(?:_\{?\d+\}?)?)+", lbl) else lbl
    if "_" in plain: return lbl, False
    m = re.fullmatch(r"H(\d*)([A-Z][a-z]?)", plain)
    if m and (m.group(2) in HETMET or m.group(2) in ("C", "F", "Cl", "Br", "I", "X", "Y", "Z", "Nu", "A", "E")):
        return m.group(2) + "H" + m.group(1), True       # (HY, HX, HNu: the generic H-Y bonds through Y)
    m = re.fullmatch(r"(F|Cl|Br|I)(\d+)(C|Si)", plain)
    if m: return m.group(3) + m.group(1) + m.group(2), True
    m = re.fullmatch(r"(F|Cl|Br|I)([A-Z][a-z]?)", plain)
    if m and m.group(2) in METALS2: return m.group(2) + m.group(1), True
    m = re.fullmatch(r"((?:tert|sec|iso|t|i|n|s|c)?(?:%s|R)'*)(\d*)(H\d*)?([A-Z][a-z]?)" % "|".join(ABBR), plain)
    if m and m.group(4) in HETMET: return m.group(4) + (m.group(3) or "") + m.group(1) + m.group(2), True
    grp = r"((?:tert|sec|iso|t|i|n|s|c)?(?:%s|R|H)'*)" % "|".join(ABBR)
    m = re.fullmatch(grp + r"(\d*)O2C", plain)                       # EtO2C, MeO2C, tBuO2C, HO2C -> CO2Et, CO2Me, CO2H
    if m: return "CO2" + m.group(1) + m.group(2), True
    m = re.fullmatch(grp + r"(\d*)C(H\d*)?", plain)                  # R2CH, Me2CH, Ph3C, PhCH2, R2C -> CHR2, CHMe2, CPh3, CH2Ph, CR2
    if m and m.group(1) != "H": return "C" + (m.group(3) or "") + m.group(1) + m.group(2), True
    m = re.fullmatch(r"(?:HO(H\d*)C|HOC(H\d*))", plain)              # HOH2C, HOCH2 -> CH2OH
    if m: return "C" + (m.group(1) or m.group(2)) + "OH", True
    m = re.fullmatch(r"((?:%s|R)'*)(\d*)SO2" % "|".join(ABBR), plain)   # PhSO2, MeSO2 -> SO2Ph, SO2Me
    if m: return "SO2" + m.group(1) + m.group(2), True
    REV = {"O2N": "NO2", "NC": "CN", "OHC": "CHO", "HOOC": "COOH", "H3CO": "OCH3", "H3CS": "SCH3", "F3CO": "OCF3",
           "(HO)2B": "B(OH)2", "pinB": "Bpin", "H3CN": "NCH3"}
    if plain in REV: return REV[plain], True
    return lbl, False

def _formula_only(lbl):
    """True for a formula of element symbols and counts after its anchor (PF6, BF4, HSO4, NO3): read only as typed."""
    u = label_units(lbl)
    return len(u) > 1 and all((lambda c: c in ELEMENTS and c not in ABBR and c != "H")(_ucore(x)) for x in u[1:])

def _sub_unit(u):
    """A label unit that can stand before the atom it hangs from: H, an alkyl, aryl or acyl group (Me2, tBu, Ph, Ac, R1,
    Ar'), a ligand (L, Ln, L_2, Cp*) or a bracketed group ((PPh3)4, (Ph))."""
    c = _ucore(u)
    return u[:1] in "([" or c in ("H", "L", "X", "R", "M") or c in ABBR or c in ("Cp*",) or bool(_alkyl_at(c, 0) == c)

def _anchor(lbl):
    """How a label is stored and read: (label stored anchor-first, typed reading leftwards, pivot). Besides the
    reversals of anchor_first (HO, MeO, Et3N, F3C ...): a metal written after its ligands (L2Pd, LnNi^{II}, (PPh3)4Pd^0,
    L2(Ph)Pd^{II}, LnM) is anchored on the metal and reads that way leftwards; a heteroatom or metal written between
    groups (EtOH, MeNH2, iPr2NEt, PhSH, L2Pd(Br)Ph, LnPdX) is anchored where it stands, its left part stored last:
    pivot is the number of units drawn to the anchor's left, in the order typed."""
    lab, rev = anchor_first(lbl)
    if rev or not lbl or "<" in lbl or lbl.startswith("\\"): return lab, rev, 0
    units = label_units(lbl)
    if len(units) < 2: return lbl, False, 0
    cores = [_ucore(u) for u in units]
    k = next((i for i, c in enumerate(cores) if i and (c in METALS_T or c == "M")), None)
    if k is None:
        k = next((i for i, c in enumerate(cores) if i and c in ("O", "N", "S", "P")), None)
        if k is None or k == len(units) - 1 or not all(re.fullmatch(r"H\d*", u) or _sub_unit(u) for u in units[k + 1:]):
            return lbl, False, 0                           # (MeOCH2: the CH2 holds the bond, not the O)
    if not all(_sub_unit(u) for u in units[:k]): return lbl, False, 0
    if k == len(units) - 1: return units[k] + "".join(reversed(units[:k])), True, 0
    return units[k] + "".join(units[k + 1:]) + "".join(units[:k]), False, k

SIDED = {"NC"} | {"H" + g for g in ("A", "X", "Y", "Z", "E", "B", "Nu")}   # read from the bond outwards (Atom._nc_iso)
class Atom:
    def __init__(self, mol, label="C", xy=(0, 0), name=None):
        self.mol, self.label, self.xy, self.name = mol, (label or "C"), np.array(xy, float), name
        self.uid = next(_UID)
        self.charge, self.radical, self.lps = 0, False, []
        self.ldir = None; self.show = False; self.color = None
        self.charge_at = None; self.radical_at = None; self.partial = None; self.number = None
        self.hidden = False; self.clear = None; self.generic = False; self.site = None; self._site = None
        self.ambi = False                      # the + - notation of a carbene (empty p, filled sp2): net charge 0
        self.z = None                          # depth (BU, + toward the viewer) where bonds cross: the back one is broken
        self._g = None; self._avoid = []; self._lin = False    # _lin: an sp carbon between collinear double bonds
    @property
    def label(self):
        if self._typed in SIDED and self._nc_iso(): return self._typed   # tBu-NC: an isocyanide, bonded through its N
        return self._label
    @label.setter
    def label(self, v):                       # HO, MeO, BrMg ... are stored anchor-first; a lone species keeps its typed order
        self._typed = v if v is not None else "C"
        self._label, self._tl, self._pivot = _anchor(self._typed)
    @property
    def _typed_left(self): return self._tl and not (self._typed in SIDED and self._nc_iso())
    def _nc_iso(self):
        """NC typed: a cyano group written leftwards when its bond comes from the right (NC-CH2-), otherwise (tBu-NC, or
        a bond above or below) an isocyanide, N bonded, as read. NC and CN are different groups, never interchanged. A
        generic acid (HA, HX, HY) reads the same way: HY-R bonds through Y, O-HA through its H."""
        if self.ldir: return self.ldir != "left"
        if self.mol is None: return False
        vs = [unit(n.xy - self.xy) for n in self.mol.nbrs(self)]
        return bool(vs) and sum(v[0] for v in vs) <= 0.3
    @property
    def drawn(self): return bool((self.show or self._lin or self.label not in ("C", "")) and not self.hidden)
    def lp(self, angle=None):
        l = LP(self, angle); self.lps.append(l); return l
    def add_lps(self, n): return [self.lp() for _ in range(n)]
    def __repr__(self): return f"Atom({self.label}{'' if self.name is None else ':' + str(self.name)})"

class Bond:
    STYLES = ("plain", "wedge", "hash", "bold", "dash", "dot", "dative", "partial", "wavy")
    def __init__(self, a, b, order=1, style="plain", side=None, free=False, color=None):
        assert style in self.STYLES, style
        self.a, self.b, self.order, self.style = a, b, order, style
        self.side, self.free, self.color = side, free, color
        self.behind = False                    # True: broken wherever it crosses another bond (drawn behind it)

def _turn_fixed(a, f):
    """Turn atom a's fixed directions with its fragment (lone pairs placed by hand, a pinned charge or radical): f maps
    an old angle (deg) to the new one."""
    for l in a.lps:
        if l.angle is not None and not l.auto: l.angle = f(l.angle)
    for k in ("charge_at", "radical_at"):
        if getattr(a, k, None) is not None: setattr(a, k, f(getattr(a, k)))

class Mol:
    def __init__(self): self.atoms, self.bonds = [], []
    # --- building
    def atom(self, label="C", xy=(0, 0), name=None, **kw):
        a = Atom(self, label, xy, name)
        for k, v in kw.items(): setattr(a, k, v)
        self.atoms.append(a); return a
    def bond(self, a, b, order=1, style="plain", **kw):
        old = self.find(a, b)
        if old:
            old.order, old.style = order, style
            for k, v in kw.items(): setattr(old, k, v)
            return old
        bb = Bond(a, b, order, style, **kw); self.bonds.append(bb); return bb
    def find(self, a, b):
        for x in self.bonds:
            if (x.a is a and x.b is b) or (x.a is b and x.b is a): return x
    def unbond(self, a, b):
        """Break the bond a-b. A ring vertex left with one substituent that stood 30 deg off the exterior bisector (one of
        a gem pair, e.g. after its H leaves) swings back onto the bisector, the reverse of attach() making the pair."""
        x = self.find(a, b)
        if not x: return
        self.bonds.remove(x)
        for v in (a, b):
            rb = [n for n in self.nbrs(v) if self.find(v, n) is not None and self.ring_of(self.find(v, n))]
            ex = [n for n in self.nbrs(v) if n not in rb]
            if len(rb) != 2 or len(ex) != 1: continue
            bis = ang(-(unit(rb[0].xy - v.xy) + unit(rb[1].xy - v.xy))); off = adiff(ang(ex[0].xy - v.xy), bis)
            if 25 <= abs(off) <= 35 and v not in self.branch(ex[0], v): self.swing(ex[0], v, -off)
    def remove(self, a):
        self.bonds = [b for b in self.bonds if b.a is not a and b.b is not a]; self.atoms.remove(a)
    def nbrs(self, a):
        return [b.b if b.a is a else b.a for b in self.bonds if a in (b.a, b.b)]
    def color(self, c, atoms=None):
        """Colour the whole molecule, or only the fragment `atoms`: their labels, the bonds between them, and their
        charges, lone pairs, radicals and δ marks take colour c ('#1F4FD1', 'tab:blue' ...; None = back to ink).
        Single bonds or atoms: set b.color / a.color directly. Colours travel with copy(), merge() and the helpers,
        so a coloured fragment stays coloured through every step of a mechanism. Returns the molecule."""
        sel = {id(a) for a in (self.atoms if atoms is None else atoms)}
        for a in self.atoms:
            if id(a) in sel: a.color = c
        for b in self.bonds:
            if id(b.a) in sel and id(b.b) in sel: b.color = c
        return self
    colour = color
    def std_dirs(self, a, k=1, toward=None):
        """ChemDraw-style directions for k new substituents on atom a: 120 deg zigzag from one bond,
        exterior bisector from two (gem pair = bisector +/- 30), largest-gap bisector from three or more.
        For k=1 several candidates may be returned; attach() picks the one pointing outward."""
        ds = sorted([ang(n.xy - a.xy) % 360 for n in self.nbrs(a)] + ([ang(toward.xy - a.xy) % 360] if toward is not None else []))
        if not ds: return [0.0, 180.0, 90.0, -90.0][:max(k, 1)] if k > 1 else [0.0]
        if len(ds) == 1:
            p = ds[0]
            return [p + 120, p - 120] if k <= 2 else [p + 90, p + 180, p - 90]
        gaps = [((ds[(i + 1) % len(ds)] - ds[i]) % 360) or 360 for i in range(len(ds))]
        i = int(np.argmax(gaps)); e = ds[i] + gaps[i] / 2
        if len(ds) == 2 and min(gaps) > 170:                   # linear: perpendicular
            return [ds[0] + 90, ds[0] - 90]
        if k == 1: return [e]
        if k == 2: return [e - 30, e + 30] if len(ds) == 2 else [e - 25, e + 25]
        return [e - 40, e, e + 40]
    def free_angle(self, a, length=1.0, toward=None):
        """Best direction for one new substituent: a standard direction with the most clearance."""
        cands = self.std_dirs(a, 1, toward); others = [x.xy for x in self.atoms if x is not a]
        prev = [n for n in self.nbrs(a)]
        def clear(t):
            p = a.xy + dirv(t) * length
            return min([norm(p - q) for q in others] or [9.0])
        def zig(t):                                            # continue a zigzag: parallel to the bond two back
            if len(prev) != 1: return 0.0
            g = [x for x in self.nbrs(prev[0]) if x is not a]
            return -min([abs(adiff(t, ang(prev[0].xy - x.xy))) for x in g] or [0.0])
        return max(cands, key=lambda t: (round(min(clear(t), 1.6), 2), zig(t), round(dirv(t)[1], 3), round(dirv(t)[0], 3)))
    def attach(self, a, angle=None, label="C", order=1, style="plain", length=1.0, name=None, toward=None, **kw):
        """New atom bonded to `a` along `angle` (deg). angle=None picks the standard direction (120 deg
        zigzag / exterior bisector) that points outward; toward=atom reserves the direction of a partner
        the atom is about to react with (e.g. the H a base will remove)."""
        if angle is None and getattr(a, "_out", None) is not None and toward is None:
            angle = a._out                                 # on a perspective ring: straight out from the ring, in its plane
        if angle is None:
            rb = [n for n in self.nbrs(a) if self.ring_of(self.find(a, n))]
            ex = [n for n in self.nbrs(a) if n not in rb]
            if len(rb) == 2 and len(ex) == 1 and toward is None:
                bis = ang(-(unit(rb[0].xy - a.xy) + unit(rb[1].xy - a.xy)))
                if abs(adiff(ang(ex[0].xy - a.xy), bis)) < 3:     # a gem pair on a ring vertex: the existing branch
                    br = ex[0]; subl = self.branch(br, a); sub = {id(x) for x in subl}    # swings 30 deg to the side
                    others = [x for x in self.atoms if id(x) not in sub and x is not a]  # where it (and the new atom
                    def room(sg):                                                          # on the other side) has room
                        pts = [(a.xy + rotv(x.xy - a.xy, sg * 30), self._lab_w(x)) for x in subl]
                        pts.append((a.xy + dirv(bis - sg * 30) * length, 0.3 if label not in ("C", "") else 0.0))
                        return round(min([norm(u - v.xy) - w - self._lab_w(v) for u, w in pts for v in others] or [9.0]), 2)
                    big = max(rb, key=lambda x: (len([n for n in self.nbrs(x) if n is not a]), -x.uid))
                    def hside(sg):                     # an explicit H sits next to the more substituted ring
                        return -abs(adiff(bis - sg * 30, ang(big.xy - a.xy))) if label == "H" else 0   # neighbour
                    sg = max((1, -1), key=lambda g: (min(room(g), 0.6), hside(g), g))
                    self.swing(br, a, sg * 30); angle = bis - sg * 30
            elif not rb and len(ex) == 3 and toward is None:
                # a chain vertex (two bonds about 120 deg apart, a branch on their exterior bisector) taking a fourth
                # substituent: the branch swings 30 deg and the new one takes the other side (a gem pair at +-30), as on
                # rings, instead of squeezing in 60 deg from both chain bonds; the side with more room wins
                cand = []
                for k in range(3):                     # the branch: the smallest substituent, the latest on a tie (in a
                    br, (p, q) = ex[k], [x for j, x in enumerate(ex) if j != k]       # perfect Y every bond qualifies)
                    d1, d2 = unit(p.xy - a.xy), unit(q.xy - a.xy); gap = abs(adiff(ang(d1), ang(d2)))
                    subl = self.branch(br, a)
                    if 100 <= gap <= 140 and abs(adiff(ang(br.xy - a.xy), ang(-(d1 + d2)))) < 3 and p not in subl and q not in subl:
                        cand.append((len(subl), -br.uid, k, subl, p, q, ang(-(d1 + d2))))
                if cand:
                    _, _, k, subl, p, q, bis = min(cand, key=lambda c: c[:3]); br = ex[k]
                    sub = {id(x) for x in subl}
                    others = [x.xy for x in self.atoms if id(x) not in sub and x is not a]
                    def room(sg):
                        pts = [a.xy + dirv(bis + sg * 30) * length] + \
                              [a.xy + rotv(x.xy - a.xy, -sg * 30) for x in self.atoms if id(x) in sub]
                        return round(min([norm(u - v) for u in pts for v in others] or [9.0]), 2)
                    big = max((p, q), key=lambda x: (len(self.branch(x, a)), -x.uid))
                    def hside(sg):                     # an explicit H sits next to the bigger neighbour (the carbonyl
                        t = bis + sg * 30 if label == "H" else (bis - sg * 30 if br.label == "H" else None)   # side of an
                        return 0 if t is None else -abs(adiff(t, ang(big.xy - a.xy)))                      # alpha-H)
                    sg = max((1, -1), key=lambda g: (room(g), hside(g)))
                    self.swing(br, a, -sg * 30); angle = bis + sg * 30
            if angle is None: angle = self.free_angle(a, length, toward)
        bkw = {k: kw.pop(k) for k in ("side", "free", "color") if k in kw}
        if abs(length - 1.0) > 1e-9: bkw.setdefault("free", True)    # a length asked for is deliberate
        n = self.atom(label, a.xy + dirv(angle) * length, name=name, **kw)
        self.bond(a, n, order, style, **bkw); return n
    def chain(self, a, angles, labels=None, orders=None):
        out = []
        for i, t in enumerate(angles):
            a = self.attach(a, t, (labels or {}).get(i, "C") if isinstance(labels, dict) else
                            (labels[i] if labels and labels[i] else "C"),
                            (orders[i] if orders else 1)); out.append(a)
        return out
    def ring(self, n, center=(0, 0), rot=90.0, labels=None, kekule=False):
        """Regular n-ring, vertex i at angle rot + 360*i/n about center (CCW)."""
        R = 1.0 / (2 * math.sin(math.pi / n)); c = np.array(center, float)
        at = [self.atom((labels or {}).get(i, "C"), c + dirv(rot + 360.0 * i / n) * R) for i in range(n)]
        for i in range(n): self.bond(at[i], at[(i + 1) % n], 2 if kekule and i % 2 == 0 else 1)
        return at
    def ring_out(self, a, angle, n=6, kekule=False, labels=None, length=1.0):
        """Ring hanging off atom a: its vertex 0 is bonded to a along `angle` (bond `length` BU)."""
        R = 1.0 / (2 * math.sin(math.pi / n)); c = a.xy + dirv(angle) * (length + R)
        at = [self.atom((labels or {}).get(i, "C"), c + dirv(angle + 180 + 360.0 * i / n) * R) for i in range(n)]
        for i in range(n): self.bond(at[i], at[(i + 1) % n], 2 if kekule and i % 2 == 0 else 1)
        self.bond(a, at[0], free=abs(length - 1) > 1e-9); return at
    def ring_at(self, a, angle, n=6, labels=None, kekule=False):
        """Ring that contains atom a as vertex 0, its centre along `angle` from a (spiro centres etc.)."""
        R = 1.0 / (2 * math.sin(math.pi / n)); c = a.xy + dirv(angle) * R
        at = [a] + [self.atom((labels or {}).get(i, "C"), c + dirv(angle + 180 + 360.0 * i / n) * R) for i in range(1, n)]
        for i in range(n): self.bond(at[i], at[(i + 1) % n], 2 if kekule and i % 2 == 0 else 1)
        return at
    def ring_fuse(self, a, b, n=6, side=None, labels=None, kekule=False):
        """New ring sharing bond a-b, on the left (side=1) or right (-1) of a->b; by default on the side with more room,
        away from a's and b's other neighbours (outside a ring it is fused to, not across it). Returns new atoms."""
        u = b.xy - a.xy; ap = 1.0 / (2 * math.tan(math.pi / n))
        if side is None:
            oth = [y.xy for x in (a, b) for y in self.nbrs(x) if y is not a and y is not b]
            side = -1 if oth and cross(u, np.mean(oth, axis=0) - a.xy) > 1e-9 else 1
        c = (a.xy + b.xy) / 2 + perp(unit(u)) * side * ap; R = 1.0 / (2 * math.sin(math.pi / n))
        a0 = ang(b.xy - c); step = 360.0 / n * (1 if cross(a.xy - c, b.xy - c) > 0 else -1)
        new = [self.atom((labels or {}).get(i, "C"), c + dirv(a0 + step * (i + 1)) * R) for i in range(n - 2)]
        seq = [b] + new + [a]
        for i in range(len(seq) - 1): self.bond(seq[i], seq[i + 1], 2 if kekule and i % 2 == 1 else 1)
        return new
    def ligate(self, metal, geom="sqpl", length=None, **sites):
        """Ligands on standard sites. sqpl: E N W S | sqplx: NE NW SW SE | oct: up down fl fr bl br | tbp: up down left fr br |
        tet: up left fr br | lin: E W | bent: up mid down (the open side of a bent metallocene from sandwich(bend=...)).
        Values: label string, or (label, order). Returns {site: atom}.
        length=None spaces each one-label ligand so the bond shows the same length between the two labels whatever
        the letters (L, Br, PPh3, the ⁰ of Pd⁰): equal white either side, as a chemist evens them by eye. A number
        fixes the bond length (BU) instead."""
        G = {"sqpl": {"E": (0, "plain"), "N": (90, "plain"), "W": (180, "plain"), "S": (270, "plain")},
             "sqplx": {"NE": (45, "plain"), "NW": (135, "plain"), "SW": (225, "plain"), "SE": (315, "plain")},
             "oct": {"up": (90, "plain"), "down": (270, "plain"), "fl": (210, "wedge"), "fr": (330, "wedge"),
                     "bl": (150, "hash"), "br": (30, "hash")},
             "tbp": {"up": (90, "plain"), "down": (270, "plain"), "left": (180, "plain"),
                     "fr": (325, "wedge"), "br": (35, "hash")},
             "tet": {"up": (90, "plain"), "left": (210, "plain"), "fr": (330, "wedge"), "br": (10, "hash")},
             "lin": {"E": (0, "plain"), "W": (180, "plain")},
             "bent": {}}[geom]
        if geom == "bent":                                 # the open side of a bent metallocene; up is the upper of the pair
            o = getattr(metal, "_open", 0.0); sg = 1 if math.cos(math.radians(o)) > -1e-9 else -1
            G = {"up": (o + 45 * sg, "plain"), "mid": (o, "plain"), "down": (o - 45 * sg, "plain")}
        out = {}
        for s, v in sites.items():
            lab, order = (v, 1) if isinstance(v, str) else v
            L = 1.3 if length is None else length
            t, sty = G[s]; out[s] = self.attach(metal, t, lab, order, sty, length=L, free=abs(L - 1) > 1e-9)
            if length is None: self.find(metal, out[s])._lig = metal
        return out
    def centroid(self, atoms):
        """Hidden atom at the centroid of `atoms` (for eta bonds to rings / alkenes)."""
        return self.atom("", np.mean([a.xy for a in atoms], axis=0), hidden=True)
    def floating(self, a, b, label="R", length=0.9, inside=0.3):
        """A floating substituent (any position on the ring): its bond starts `inside` BU within the ring, crosses the
        ring edge a-b at its middle, square to it, and ends at `label` `length` BU outside, as ChemDraw draws R1 on a
        generic arene. The inner line of a double bond it crosses is broken there. Returns the label's atom."""
        bd = self.find(a, b); ring = self.ring_of(bd) if bd is not None else None
        if not ring: raise ValueError(f"floating(): {a}-{b} is not a ring bond (a floating substituent crosses a ring edge)")
        mid = (a.xy + b.xy) / 2; n = unit(mid - np.mean([x.xy for x in ring], axis=0))
        p = self.atom("", mid - n * inside, hidden=True); p._float = (a, b)
        r = self.atom(label, mid + n * length); self.bond(p, r, free=True); return r
    def eta(self, metal, atoms):
        """eta-n bond from a metal to a pi system (eta-2 alkene, eta-5 Cp, eta-6 arene): it points at the centroid and
        stops where it meets the pi system's nearest line, as in ChemDraw."""
        c = self.centroid(atoms); c._eta = list(atoms); self.bond(metal, c, free=True)
        if all(getattr(x, "_tmpl", None) == "persp" for x in atoms): c._persp3d = True     # drawn with depth, see persp_ring
        return c
    def persp_ring(self, n=5, center=(0, 0), rot=None, squash=0.38, front="bold", tilt=None, size=1.45, toward=None, axis=None):
        """An n-ring seen in perspective, as in sandwich and half-sandwich complexes: a regular ring seen slightly from
        above, its depth squashed to `squash` of its width. front="bold" draws it as ChemDraw's perspective template
        does: the near edges a bold bar, the edges running back from it tapering to the plain width at the back, all
        joined as one mitred outline; "plain" draws every edge alike. toward=metal (an atom or point) sets the ring
        square to the line from the metal to its centre (a bent metallocene's rings, a ring off to one side); axis=deg
        does the same for a direction; tilt=deg turns the ring in the page directly. rot puts a vertex at that angle in
        the ring's own frame (90: at the back); by default the ring has an edge, not a vertex, where the metal's line
        crosses it: at the front when the metal is in front of the ring (the top ring of ferrocene), at the back when it
        is behind (the bottom ring). The atoms are unlabelled carbons (attach() points substituents straight out of the
        ring); the bonds are free (foreshortened). size: the widest edges in BU (1.45 keeps the flattened ring legible
        next to its labels). eta(metal, ring) to a perspective ring draws its line with depth: behind the near edges,
        in front of the far ones. Returns the atoms, counter-clockwise (seen from above) from the one at `rot`: for
        ferrocene's top ring [back, back left, front left, front right, back right]."""
        R = size / (2 * math.sin(math.pi / n)); c = np.array(center, float)
        P = None if toward is None else (np.array(toward.xy, float) if isinstance(toward, Atom) else np.array(toward, float))
        if P is not None and axis is None and norm(c - P) > 1e-9: axis = ang(c - P)
        if tilt is None:                                   # the ring's frame: +y is whichever half of its axis points up the
            tilt = 0.0 if axis is None else (axis - 90.0 if math.sin(math.radians(axis)) > -1e-9 else axis + 90.0)  # page
        tilt = (tilt + 180.0) % 360.0 - 180.0
        if rot is None:
            back = P is not None and rotv(P - c, -tilt)[1] > 1e-9          # the metal sits behind the ring
            rot = (90.0 if back else 270.0) - 180.0 / n
        ths = [rot + 360.0 * i / n for i in range(n)]
        at = [self.atom("C", c + rotv(np.array([math.cos(math.radians(t)), squash * math.sin(math.radians(t))]) * R, tilt)) for t in ths]
        lim = -0.25 if n == 3 else -0.5
        near = [math.sin(math.radians(t)) < lim - 1e-6 for t in ths]     # the lower part of the ring is nearer the eye
        for i in range(n):
            j = (i + 1) % n; a, b = at[i], at[j]
            if front == "bold" and near[i] and near[j]: bb = self.bond(a, b, style="bold", free=True)
            elif front == "bold" and (near[i] or near[j]):
                bb = self.bond(b, a, style="wedge", free=True) if near[i] else self.bond(a, b, style="wedge", free=True)
                bb.wide = "bold"                                           # its wide end meets the bold bar flush
            else: bb = self.bond(a, b, free=True)
        for x, t in zip(at, ths):
            x._tmpl = "persp"; x._out = ang(rotv(np.array([math.cos(math.radians(t)), squash * math.sin(math.radians(t))]), tilt))
        self._persp = getattr(self, "_persp", []) + [{"atoms": at, "c": c, "R": R, "k": squash, "tilt": tilt}]
        return at
    def sandwich(self, metal="Fe", n=5, at=(0, 0), gap=1.75, squash=0.38, stagger=None, rings=("top", "bottom"),
                 aromatic=True, eta=True, size=1.45, bend=None, open=0.0):
        """Sandwich complex in perspective (ferrocene, other metallocenes, bis(benzene)chromium with n=6): flattened
        rings above and below the metal, bold near edges, an aromatic ellipse in each ring, and lines from the metal
        to the ring centres that pass behind a ring's near edge and in front of its far edge, as in a textbook
        drawing. stagger turns one ring half a step against the other (default: staggered for odd n, the D5d drawing of
        ferrocene; eclipsed for even n, as bis(benzene)chromium is); the metal's line always crosses an edge, never a
        corner. rings=("top",) gives a half-sandwich (then ligate() the metal below). bend=deg draws a bent metallocene
        (Cp2TiCl2: bend=130, the ring-metal-ring angle): each ring square to its own line from the metal, the pair opening
        towards `open` (deg, default 0: to the right), where ligate(M, "bent", up=..., mid=..., down=...) puts the other
        ligands. Returns (metal, [ring atoms ...])."""
        M = self.atom(metal, at); out = []; at = np.array(at, float); half = 180.0 / n
        if stagger is None: stagger = n % 2 == 1
        if bend is not None: M._open = float(open)
        for side in rings:
            top = side == "top"
            if bend is None:
                axis = 90.0 if top else 270.0; back = 90.0 - half       # the lower ring: an edge at the back, where the
                rot = (back + half if stagger else back) if top else back   # metal's line passes in front of it
            else:
                axis = open + 180.0 - bend / 2.0 if top else open + 180.0 + bend / 2.0; rot = None
            ring = self.persp_ring(n, at + dirv(axis) * gap, rot, squash, size=size, toward=M)
            if aromatic: self.deloc(ring)
            if eta:
                c = self.eta(M, ring); c._persp3d = True
            out.append(ring)
        return M, out
    def deloc(self, atoms, charge=0):
        """A delocalised pi system drawn as one hybrid: a ring gets a circle inside it (an ellipse in a perspective
        ring), an open chain (allyl, pentadienyl) dashed partial bonds along its inner side, mitred and symmetric;
        a charge sits at the system's centre, shared by all its atoms rather than on any one (Cp⁻, C7H7⁺, the allyl
        cation and anion). Returns the molecule."""
        atoms = list(atoms); n = len(atoms)
        ring = n >= 3 and all(self.find(atoms[i], atoms[(i + 1) % n]) for i in range(n))
        if not ring:
            c = np.mean([a.xy for a in atoms], axis=0)
            for x, y in zip(atoms, atoms[1:]):
                b = self.find(x, y)
                if b is None: raise ValueError(f"deloc: {x} and {y} are not bonded")
                b.order = 1; b.style = "partial"
                b.side = 1 if cross(y.xy - x.xy, c - x.xy) * (1 if b.a is x else -1) > 0 else -1
        self._deloc = getattr(self, "_deloc", []) + [{"atoms": atoms, "charge": charge, "ring": ring}]
        return self
    CHAIR = (66.0, -15.0, 16.0)        # ring-bond directions C1-C2, C2-C3, C3-C4; the other three are parallel to these
    def chair(self, origin=(0, 0), flip=False):
        """Chair cyclohexane as chemists draw it: steep left end, shallow top and bottom edges sloping slightly down to
        the right, the foot (atom 0) at lower left and the head (atom 3) at upper right; opposite ring bonds parallel,
        axial bonds vertical, each equatorial bond parallel to the ring bonds one removed. Every axial and equatorial
        bond is at least 50 deg from the ring bonds at its atom. flip=True draws the ring-flipped chair (the mirror image
        top to bottom: the same carbons, the head and foot swapped), where every axial position becomes equatorial on the
        same face. Returns (atoms, ax, eq): the axial and equatorial directions (deg) at each atom."""
        o = np.array(origin, float); a, b, c = self.CHAIR; dirs = [a, b, c, a + 180, b + 180, c + 180]
        if flip: dirs = [-d for d in dirs]
        p = o.copy(); at = []
        for d in dirs:
            at.append(self.atom("C", p)); p = p + dirv(d)
        for i in range(6): self.bond(at[i], at[(i + 1) % 6])
        for x in at: x._tmpl = "chair"                            # perspective template: 2D angle rules do not apply
        up = [False, True, False, True, False, True]
        if flip: up = [not u for u in up]
        ax = [90.0 if u else 270.0 for u in up]; cen = np.mean([x.xy for x in at], axis=0); eq = []
        for i in range(6):
            d = dirs[(i + 1) % 6]                                  # ring bond one removed is parallel to the equatorial bond
            e = d if np.dot(dirv(d), at[i].xy - cen) > 0 else d + 180
            eq.append(e % 360)
        return at, ax, eq
    def chair_flip(self, origin=(0, 0)):
        """The ring-flipped chair (chair(flip=True)): a substituent axial on atom i in chair() is equatorial on atom i
        here, on the same face of the ring, and vice versa."""
        return self.chair(origin, flip=True)
    def newman(self, front, back, front_rot=90.0, back_rot=30.0, at=(0, 0), r=0.42):
        """Newman projection: front=[3 labels] on bonds from the centre, back=[3 labels] from the circle edge."""
        c = np.array(at, float); f = self.atom("C", c); bk = self.atom("C", c, hidden=True); bk.clear = None
        self._newman = getattr(self, "_newman", []) + [(f, r)]
        L = 1.15                                                   # a little longer, so front and back labels never touch
        fs = [self.attach(f, front_rot + 120 * i, front[i] or "C", length=L) for i in range(3)]
        bs = [self.attach(bk, back_rot + 120 * i, back[i] or "C", length=L) for i in range(3)]
        for b in self.bonds:
            if b.a in (f, bk) or b.b in (f, bk): b.free = True
        f._tmpl = bk._tmpl = "newman"; bk._newman_r = r; f._newman_pair = bk; bk._newman_pair = f; return f, bk, fs, bs
    def stereo(self, center, other, kind="wedge"):
        """Set the bond center-other to 'wedge' (towards viewer), 'hash' (away), 'wavy' (unknown), 'bold' or 'plain',
        always with the narrow end on the stereocentre."""
        b = self.find(center, other)
        if b is None: raise ValueError("stereo: no bond between these atoms")
        if b.a is not center: b.a, b.b = b.b, b.a
        b.style = kind; return b
    def bridge(self, a, b, k, away=None, labels=None):
        """k new atoms joining a and b with unit bonds on a symmetric arc (chelate or bridged rings),
        bulging away from atom `away` (e.g. the metal). labels={i: 'O'} sets new atom labels."""
        c = norm(b.xy - a.xy); n = k + 1
        if not 0 < c < n: raise ValueError("bridge: atoms must be closer than k+1 bond lengths")
        theta = lambda R: n * 2 * math.asin(1 / (2 * R))
        lo, hi = 1 / (2 * math.sin(math.pi / n)) + 1e-12, 1e6
        for _ in range(200):
            R = (lo + hi) / 2
            if 2 * R * math.sin(theta(R) / 2) < c: lo = R
            else: hi = R
        R = (lo + hi) / 2; th = theta(R); beta = th / n
        u = unit(b.xy - a.xy); out = perp(u)
        if away is not None and np.dot(away.xy - (a.xy + b.xy) / 2, out) > 0: out = -out
        C = (a.xy + b.xy) / 2 - out * R * math.cos(th / 2); fa = ang(a.xy - C)
        sg = min((1, -1), key=lambda s: norm(C + R * dirv(fa + s * math.degrees(th)) - b.xy))
        new = [self.atom((labels or {}).get(i, "C"), C + R * dirv(fa + sg * math.degrees(beta) * (i + 1))) for i in range(k)]
        seq = [a] + new + [b]
        for i in range(n): self.bond(seq[i], seq[i + 1])
        return new
    # --- editing between mechanism steps
    def copy(self, keep_lps=False):
        """Deep copy for the next mechanism step (same names/uids, same coordinates). Lone pairs and site dots (S.dot) are
        dropped: they belong to the step they are drawn in (a Nu that has reacted is no longer one)."""
        m = _copy.deepcopy(self); m._moved = set(); m._fresh = set()
        for a in m.atoms:
            a._g = None; a.site = None; a._site = None; a._arrowed = False
            if not keep_lps: a.lps = []
        return m
    def numbered(self, start=1):
        """A reference drawing for a step list: a copy with every atom numbered start, start + 1, ... (in self.atoms
        order, each number beside its atom). The numbers stay with the atoms of this molecule and of every step copied
        from it, so the user can write "lp on N1 -> H5" and num(1), num(5) find those atoms in any step."""
        for k, a in enumerate(self.atoms): a.ref = start + k
        m = self.copy()
        for a in m.atoms: a.number = a.ref
        return m
    def num(self, k):
        """The atom numbered k in the reference drawing (numbered()), in this molecule or any step copied from it."""
        for a in self.atoms:
            if getattr(a, "ref", None) == k: return a
        raise KeyError(f"no atom numbered {k} here (numbered() first; merged reagents keep their own numbers)")
    def __getitem__(self, k):
        if isinstance(k, Atom):
            for a in self.atoms:
                if a.uid == k.uid: return a
            raise KeyError(k)
        for a in self.atoms:
            if a.name == k: return a
        raise KeyError(k)
    get = __getitem__
    def move(self, dx, dy):
        for a in self.atoms: a.xy = a.xy + (dx, dy)
        return self
    def center(self): return np.mean([a.xy for a in self.atoms], axis=0)
    def branch(self, atom, root):
        """Atoms reachable from `atom` without passing through `root` (the substituent hanging off root)."""
        seen, todo = {atom.uid: atom}, [atom]
        while todo:
            x = todo.pop()
            for y in self.nbrs(x):
                if y is not root and y.uid not in seen: seen[y.uid] = y; todo.append(y)
        return list(seen.values())
    def swing(self, atom, pivot, deg):
        """Rotate the branch starting at `atom` about `pivot` by deg (bond lengths unchanged)."""
        br = self.branch(atom, pivot)
        if pivot in br: raise ValueError("swing: branch is part of a ring through the pivot")
        for x in br: x.xy = pivot.xy + rotv(x.xy - pivot.xy, deg)
        for x in br: _turn_fixed(x, lambda v: v + deg)
        self._moved = getattr(self, "_moved", set()) | {x.uid for x in br}       # a deliberate move (continuity)
        return self
    def frag(self, atom):
        """Atoms of the connected fragment containing `atom`. A pi system bonded eta-n to a metal is connected to it
        through its centroid, so a sandwich or half-sandwich complex is one fragment (it moves, docks and centres as one)."""
        eta = {}
        for x in self.atoms:
            for y in getattr(x, "_eta", None) or []:
                if y.mol is self: eta.setdefault(x.uid, []).append(y); eta.setdefault(y.uid, []).append(x)
        seen, todo = {atom.uid: atom}, [atom]
        while todo:
            x = todo.pop()
            for y in self.nbrs(x) + eta.get(x.uid, []):
                if y.uid not in seen: seen[y.uid] = y; todo.append(y)
        return list(seen.values())
    @staticmethod
    def _dock_dist(a, b):
        """Standard distance (BU) between two reacting partners: 2.0 for a proton transfer or a lone-pair donor (a
        labelled atom: O, N, S, halide ...), 1.7 between two carbons (an enolate, alkene or organometallic carbon
        attacking a carbon). These are the smallest distances that leave every curly arrow a shaft of 0.8 BU."""
        if a.radical or b.radical: return 2.2                 # two half-arrows share the gap (radical abstraction)
        if a.label == "H" or b.label == "H": return 2.0
        return 2.0 if (a.drawn or b.drawn) else 1.7
    def _gap(self, mover, anchor, pos):
        """Smallest distance (BU) between the moved fragment (positions pos) and everything else, leaving out the
        reacting pair itself."""
        P = [(pos[id(y)], y) for y in self.frag(mover)]; ids = {id(y) for _, y in P}
        out = [norm(p - q.xy) for p, y in P for q in self.atoms if id(q) not in ids and not q.hidden and not (y is mover and q is anchor)]
        return min(out + [9.0])
    def dock(self, atom, target, angle, dist=None):
        """Place `atom` dist BU from `target` along `angle`: moves atom's fragment if both are in this Mol,
        otherwise the whole molecule. dist=None takes the standard distance (_dock_dist: 2.0 BU for a proton transfer
        or a lone-pair donor, 1.7 BU carbon to carbon). deprotonate(), protonate() and attack() choose the angle too.
        An atom of another Mol docked against this one is merged in (as attack() does), so arrows and validate() see it."""
        if target in self.atoms and atom not in self.atoms: atom = self._adopt(atom)
        if dist is None:                                        # the standard distance, a little more if that would crowd
            d0 = self._dock_dist(atom, target)                  # the reagent's other atoms against the substrate
            if target in self.atoms and atom in self.atoms:
                for k in range(7):
                    dist = d0 + 0.1 * k; dv = target.xy + dirv(angle) * dist - atom.xy
                    if self._gap(atom, target, {id(y): y.xy + dv for y in self.frag(atom)}) >= 1.1: break
            else: dist = d0
        d = target.xy + dirv(angle) * dist - atom.xy
        if target in self.atoms:
            for x in self.frag(atom): x.xy = x.xy + d
            for x, y in ((atom, target), (target, atom)):  # a lone CH3, HO-, Et3N reads with its anchor towards the partner
                if not self.nbrs(x) and x.ldir is None and len(label_units(x.label)) > 1 and not re.search(r"[_^<]", x.label) \
                        and not _formula_only(x.label):     # (PF6, BF4: a formula is never written backwards)
                    cx = float(np.cos(np.radians(ang(y.xy - x.xy))))
                    if cx > 0.3: x._face = "left"
                    elif cx < -0.3: x._face = "right"
            return self
        return self.move(*d)
    # --- intent-level placement of reagents (standard distance and angle chosen by rule)
    def _adopt(self, x, donor=False):
        """x as an atom of this Mol: a stock reagent Mol (hydroxide() ...) or an atom of another Mol is merged in (and so
        counts as the newly added species of this step); returns the atom. A Mol without .key gives its first atom, or
        with donor=True (a base, a nucleophile) its donor atom: an anion, else an atom with lone pairs drawn, else its
        first N, O, S or P (pyridine from SMILES gives its N)."""
        if isinstance(x, Mol):
            k = getattr(x, "key", None)
            if k is None and donor:
                for test in (lambda a: a.charge < 0, lambda a: bool(a.lps), lambda a: self._elem(a) == "N",
                             lambda a: self._elem(a) == "O", lambda a: self._elem(a) in ("S", "P")):
                    k = next((a for a in x.atoms if test(a)), None)
                    if k is not None: break
            k = k or (x.atoms[0] if x.atoms else None)
            if k is None: raise ValueError("empty reagent Mol")
            self.merge(x); return k
        if x not in self.atoms: self.merge(x.mol)
        return x
    def _proton(self, x):
        """The H an acid gives up, adopted like _adopt(): an H atom; of a reagent Mol its .h (water()), its key if that is
        an H or a lone species written with its H (H3O+, HCl), else its one H on a heteroatom (a hand-built acid).
        Anything else is refused: docking another atom as the proton turns the arrows into nonsense."""
        if isinstance(x, Mol):
            k = getattr(x, "h", None) or getattr(x, "key", None)
            if k is None or (k.label != "H" and x.nbrs(k)):
                hs = [a for a in x.atoms if a.label == "H" and len(x.nbrs(a)) == 1 and x.nbrs(a)[0].drawn]
                if len(hs) != 1: raise ValueError("protonate: say which H the acid gives up (pass that H atom, or set the acid's .key to it)")
                k = hs[0]
            self.merge(x); return k
        x = self._adopt(x)
        if x.label != "H" and self.nbrs(x): raise ValueError(f"protonate: {x!r} is not a proton (pass the acid's H atom)")
        return x
    def _settle(self, x, P, t_in, tx=None, commit=True):
        """Put x at P (BU) with its fragment turned so that t_in (deg, from x towards its partner) lands on one of the
        directions tx at x (deg, in the fragment's current frame; by default x's standard directions for a new bond, as
        if it were the bond x is about to make; an H keeps its own bond in line). Of the turns (and mirror images) that do
        this, the one with most room wins. Returns {id(atom): new xy}."""
        fx = self.frag(x); ids = {id(y) for y in fx}; keep_at, keep_b = self.atoms, self.bonds
        nb = self.nbrs(x)
        if tx is None:
            if x.label == "H" and len(nb) == 1: tx = [ang(x.xy - nb[0].xy)]
            elif not nb: tx = [t_in]
            else:
                try:
                    self.atoms = list(fx); self.bonds = [bb for bb in keep_b if id(bb.a) in ids]
                    tx = self.std_dirs(x, 1)
                finally: self.atoms, self.bonds = keep_at, keep_b
        stereo = any(bb.style in ("wedge", "hash") for bb in self.bonds if id(bb.a) in ids)
        rest = [y.xy for y in self.atoms if id(y) not in ids and not y.hidden]; P = np.asarray(P, float); P0 = x.xy.copy()
        def place(t, mir):
            u = dirv(t); out = {}
            for y in fx:
                v = y.xy - P0
                if mir: v = 2 * (v @ u) * u - v
                out[id(y)] = P + rotv(v, t_in - t)
            return out
        opts = [(t, mir) for t in tx for mir in ((False, True) if not stereo and len(fx) > 2 else (False,))]
        def room(o):
            pos = place(*o); d = [norm(pos[id(y)] - q) for y in fx if y is not x for q in rest]
            return round(min(d + [9.0]), 2)
        t, mir = max(opts, key=lambda o: (min(room(o), 1.2), -opts.index(o)))
        pos = place(t, mir)
        if commit:
            rot = t_in - t; self._moved = getattr(self, "_moved", set()) | {y.uid for y in fx}   # deliberate (continuity)
            for y in fx:
                y.xy = pos[id(y)]; _turn_fixed(y, lambda v: (2 * t - v if mir else v) + rot)
            if not nb and x.ldir is None and len(label_units(x.label)) > 1 and not re.search(r"[_^<]", x.label) \
                    and not _formula_only(x.label):
                cx = math.cos(math.radians(t_in))           # a lone HO-, EtO-, Et3N reads with its anchor towards the partner
                if cx > 0.3: x._face = "left"               # (while it stays lone: once bonded, its bonds decide)
                elif cx < -0.3: x._face = "right"
        return pos
    def _fresh_side(self, x, y):
        """True if y's species was added to this step (merged in) and x's was not: then y is the one that moves."""
        fr = getattr(self, "_fresh", set())
        return any(id(a) in fr for a in self.frag(y)) and not any(id(a) in fr for a in self.frag(x))
    def _box_hw(self, pos):
        """(height, area) of the drawing's atom box with some atoms moved (pos: {id: xy}); schemes run left to right, so
        a lower panel is the more compact one."""
        P = np.array([pos.get(id(a), a.xy) for a in self.atoms if not a.hidden])
        if not len(P): return 0.0, 0.0
        w, h = P.max(0) - P.min(0) + 1.0; return float(h), float(w * h)
    def _free_dirs(self, a, min_sep=55):
        """Directions (deg, 5-deg grid) at least min_sep from every bond at a and not pointing into a ring through a."""
        occ = [ang(n.xy - a.xy) for n in self.nbrs(a)]; out = []; rings = []
        for n in self.nbrs(a):
            r = self.ring_of(self.find(a, n))
            if r: rings.append(np.mean([y.xy for y in r], axis=0))
        for t in range(0, 360, 5):
            if occ and min(abs(adiff(t, o)) for o in occ) < min_sep: continue
            if any(abs(adiff(t, ang(c - a.xy))) < 60 for c in rings): continue
            out.append(float(t))
        return out
    def _new_bond_dirs(self, a, other):
        """Standard directions for a new bond at a, judged within a's own fragment (other's fragment left out)."""
        fo = {id(y) for y in self.frag(other)}; keep_at, keep_b = self.atoms, self.bonds
        try:
            self.atoms = [y for y in keep_at if id(y) not in fo]; self.bonds = [b for b in keep_b if id(b.a) not in fo]
            if not self.nbrs(a): return [float(t) for t in range(0, 360, 30)]
            ds = self.std_dirs(a, 1)
            if len(self.nbrs(a)) == 1 and self.find(a, self.nbrs(a)[0]).order >= 2: ds = list(ds)     # =CH2, C=O: two zigzag spots
            return ds
        finally: self.atoms, self.bonds = keep_at, keep_b
    def _sight(self, mover, anchor, pos):
        """Clearance (BU) of the straight run from the mover (at pos) to the anchor: the nearest other atom's distance
        from it. The arrow between them runs there, so a substituent within about 30 deg of that line is in its way."""
        P = pos[id(mover)]; Q = anchor.xy; ids = {id(y) for y in self.frag(mover)}; v = Q - P; L2 = max(float(v @ v), 1e-9)
        d = [norm(P + v * max(0.0, min(1.0, float((y.xy - P) @ v) / L2)) - y.xy)
             for y in self.atoms if id(y) not in ids and y is not anchor and not y.hidden]
        return min(d + [9.0])
    def _place_pair(self, mover, anchor, dist, tiers, mover_dirs, memo, sight=False):
        """Move mover's fragment next to anchor: dist BU along a direction from `tiers` (lists of directions in the
        anchor's frame, the first tier preferred), turned so one of mover_dirs (None: its standard directions) points back
        at the anchor. Choice: the first tier with a spot clear of everything (1.1 BU, atom to atom; with sight=True also
        a clear straight run to the anchor, 0.5 BU from other atoms, for the arrow of an attack), at the standard
        distance if possible (else up to 0.4 BU more), then the lowest and most compact drawing, then the direction
        remembered from an earlier step (attribute `memo` on the anchor, kept by copy())."""
        prev = getattr(anchor, memo, None); best = None
        for gs in ((0.0, 0.2, 0.4), (0.7, 1.0)):        # farther out only when no direction is clear nearer (a big partner)
            for tier in [tr for tr in tiers if tr]:
                scored = []
                for k, t in enumerate(tier):
                    for g in gs:
                        pos = self._settle(mover, anchor.xy + dirv(t) * (dist + g), t + 180.0, tx=mover_dirs, commit=False)
                        cl = self._gap(mover, anchor, pos)
                        if sight and self._sight(mover, anchor, pos) < 0.5: cl = min(cl, 1.0)
                        if cl >= 1.1 or g == gs[-1]: break
                    same = prev is not None and abs(adiff(t, prev)) < 1; hh, ar = self._box_hw(pos)
                    scored.append(((cl >= 1.1), -g, round(min(cl, 1.1), 2), -round(hh, 1), -round(ar, 0), same, -k, t, g))
                top = max(scored)
                if top[0] or best is None or top[:3] > best[:3]: best = top
                if top[0]: break
            if best[0]: break
        t, g = best[-2], best[-1]
        self._settle(mover, anchor.xy + dirv(t) * (dist + g), t + 180.0, tx=mover_dirs)
        setattr(anchor, memo, float(t) % 360); return t
    def deprotonate(self, base, H, dist=None):
        """Place a base to remove proton H: in line with the X-H bond, beyond the H, dist BU from it (2.0 by default),
        turned so the base's own substituents point away on standard angles (a lone HO-, EtO- or Et3N reads with its
        anchor towards the H). If that spot is crowded, the nearest clear direction within 60 deg is used. `base` is an
        atom or a stock reagent such as hydroxide() (merged in). The species added in this step is the one that moves:
        normally the base; if the molecule carrying H was the one just added, it is placed against the base instead.
        Returns the base atom."""
        base = self._adopt(base, donor=True); X = self.nbrs(H)
        if len(X) != 1: raise ValueError("deprotonate: the H must have exactly one bond (to the atom it leaves)")
        t0 = ang(H.xy - X[0].xy); dist = self._dock_dist(base, H) if dist is None else dist
        if self._fresh_side(base, H): self._place_pair(H, base, dist, [self._new_bond_dirs(base, H)], [t0], "_dep_t")
        else: self._place_pair(base, H, dist, [[t0], [t0 + 15, t0 - 15], [t0 + 30, t0 - 30], [t0 + 45, t0 - 45], [t0 + 60, t0 - 60]], None, "_dep_t")
        return base
    def protonate(self, H, acceptor, dist=None, angle=None):
        """Place an acid so its proton H points at `acceptor`, dist BU away (2.0 by default) along the acceptor's free
        standard direction (where the new bond will point), with the acid's own H-A bond in line. `H` is the acid's H
        atom or a stock acid such as hydronium() or acid("Cl") (merged in). The species added in this step moves:
        normally the acid. Returns the H atom."""
        H = self._proton(H); dist = self._dock_dist(H, acceptor) if dist is None else dist
        A = self.nbrs(H); tH = [ang(H.xy - A[0].xy)] if len(A) == 1 else None
        if self._fresh_side(H, acceptor):
            self._place_pair(acceptor, H, dist, [tH or [0.0]], None, "_prot_t")
        else:
            ds = [angle] if angle is not None else self._new_bond_dirs(acceptor, H)
            self._place_pair(H, acceptor, dist, [ds], tH, "_prot_t")
        return H
    def attack(self, nu, E, dist=None, angle=None):
        """Place nucleophile `nu` and electrophile E for attack, dist BU apart (2.0 from a lone-pair donor, 1.7 from a
        carbon): the nucleophile comes in anti to the leaving group of a saturated centre (SN2); to a C=O or C=N carbon in
        the free sector closest to 107 deg from the multiple bond (Burgi-Dunitz, as far as a flat drawing allows);
        otherwise along E's freest standard direction. The partner added in this step moves (normally the nucleophile, or
        the electrophile when it is the one just merged in), turned so its own substituents point away on standard
        angles; a moving partner with bonds of its own keeps them on the 30-degree grid, so the approach is the nearest
        grid direction (120 deg from a C=O). Of several equally good spots the most compact drawing wins, then the one
        used in an earlier step.
        `nu` is an atom or a stock reagent (hydroxide(), halide("I") ...). Returns the nucleophile atom."""
        nu = self._adopt(nu, donor=True); dist = self._dock_dist(nu, E) if dist is None else dist
        fo = {id(y) for y in self.frag(nu)}; keep = self.atoms; sn2 = False
        try:
            self.atoms = [y for y in keep if id(y) not in fo]; free = self._free_dirs(E) or self._free_dirs(E, 40) or [self.free_angle(E)]
            nbs = [(n, self.find(E, n)) for n in self.nbrs(E)]
            lg = [n for n, b in nbs if b.order == 1 and n.drawn and (label_units(n.label)[0] in ("Cl", "Br", "I")
                  or n.label in ("OTs", "OMs", "OTf", "ONs", "OBs", "OH2", "N2", "SMe2", "NMe3") or n.charge > 0)]
            mult = [n for n, b in nbs if b.order >= 2 and n.drawn]
            plan = self._cross_plan(E) if (mult and angle is None) else None
            rc = [np.mean([y.xy for y in r], axis=0) for r in (self.ring_of(self.find(E, n)) for n in self.nbrs(E)) if r]
            inring = lambda t: any(abs(adiff(t, ang(q - E.xy))) < 60 for q in rc)
            if angle is not None: cands = [float(angle)]
            elif plan is not None and self._relay(E, plan[0]) is not None: cands = [plan[1]]
            elif lg:                                       # SN2: from the back, straight anti to the leaving group (a group
                anti = ang(E.xy - lg[0].xy) % 360          # beside that line nudges it below); only a ring in the way, which
                cands = [anti if not inring(anti) else min(free, key=lambda t: abs(adiff(t, anti)))]; sn2 = True   # a flat
                # drawing cannot get behind, falls back to the free direction nearest the back
            elif mult:
                tm = ang(mult[0].xy - E.xy); sc = {t: abs(abs(adiff(t, tm)) - 107.0) for t in free}; lo = min(sc.values())
                cands = sorted([t for t in free if sc[t] <= lo + 0.5], key=lambda t: sc[t])
            else: cands = [self.free_angle(E)]
            pick = (lambda t: sc[t]) if (mult and not sn2 and angle is None and plan is None) else (lambda t: abs(adiff(t, cands[0])))
            alt = sorted([t for t in free if all(abs(adiff(t, u)) > 1e-6 for u in cands)], key=pick)   # if those are blocked
        finally: self.atoms = keep
        fresh = self._fresh_side(nu, E)
        if angle is None and nbs and len(self.frag(E if fresh else nu)) > 1:
            # the partner that moves has bonds of its own: keep them on the 30-degree grid of E's bonds (a structure is
            # never drawn tilted), so the approach is the nearest grid direction to the rule's (120 deg from a C=O for 107)
            ref = ang(nbs[0][0].xy - E.xy); snap = []
            for t in cands:
                u = float((ref + 30.0 * round(adiff(t, ref) / 30.0)) % 360)
                if all(abs(adiff(u, w)) > 1e-6 for w in snap): snap.append(u)
            cands = snap; snap = []
            for t in alt:
                u = float((ref + 30.0 * round(adiff(t, ref) / 30.0)) % 360)
                if all(abs(adiff(u, w)) > 1e-6 for w in snap + cands): snap.append(u)
            alt = snap
        if fresh: self._place_pair(E, nu, dist, [self._new_bond_dirs(nu, E)], cands, "_atk_t", sight=True)
        elif sn2:                                          # backside attack: as near anti to the leaving group as is clear
            c0 = cands[0]; tiers, seen = [], []            # (a group beside the carbon may block the exact line)
            for dv in (0, 15, 30, 45, 60, 75):
                tr = [u for u in ((c0 + dv) % 360, (c0 - dv) % 360) if all(abs(adiff(u, w)) > 1e-6 for w in seen)
                      and (not inring(u) or inring(c0))]
                if nbs and len(self.frag(nu)) > 1:         # a partner with bonds of its own stays on E's 30-degree grid
                    ref = ang(nbs[0][0].xy - E.xy)
                    tr = [u for u in dict.fromkeys(float((ref + 30.0 * round(adiff(t, ref) / 30.0)) % 360) for t in tr)
                          if all(abs(adiff(u, w)) > 1e-6 for w in seen)]
                if tr: tiers.append(tr); seen += tr
            self._place_pair(nu, E, dist, tiers, None, "_atk_t", sight=True)
        else: self._place_pair(nu, E, dist, [cands] + ([alt] if angle is None else []), None, "_atk_t", sight=True)
        return nu
    def tidy(self, keep=(), scene=None, carbon=False):
        """Standard depiction: fold every explicit H on a heteroatom (O-H, N-H, S-H ...) into that atom's label (OH, NH2,
        H2O for a lone water), and on a condensed carbon label (CH2-H becomes CH3), except H atoms in `keep`, marked
        `deliberate`, or touched by a curly arrow of `scene`.
        carbon=True also makes the same H atoms on skeletal carbons implicit (H2C=O drawn with its H atoms becomes a
        plain CH2 vertex once it is part of a chain). Use it on the last structure of a mechanism, and on
        intermediates for any H the next step does not use."""
        used = {a.uid for a in self.atoms if getattr(a, "_arrowed", False)}     # an arrow drawn on this step uses it
        if scene is not None:
            for c in scene.curlies:
                for x in (c.src, c.dst): used.update(y.uid for y in Scene._xatoms(x))
        kid = {id(k) for k in keep}
        def spare(X, h): return (h.label == "H" and len(self.nbrs(h)) == 1 and self.find(X, h).order == 1
                                 and self.find(X, h).style == "plain" and not (h.charge or h.radical or h.lps)
                                 and h.uid not in used and id(h) not in kid and not getattr(h, "deliberate", False))
        if carbon:
            for X in list(self.atoms):
                if X.label == "C" and not X.drawn and len(self.nbrs(X)) > 1:
                    lost = False
                    for h in [h for h in self.nbrs(X) if spare(X, h)]:
                        if len([y for y in self.nbrs(X) if y.label != "H"]) >= 1: self.remove(h); lost = True
                    if lost:                               # a centre left as a T (join() drew it as a cross) gets 120 deg
                        self._trigonal(X, getattr(X, "_cross", False) and len(self.nbrs(X)) == 3); X._cross = False
        for X in list(self.atoms):
            if not X.drawn or X.hidden or re.search(r"[_^<(]", X.label): continue
            u = label_units(X.label); el = re.sub(r"\d+'*$", "", u[0])
            if el not in ("O", "N", "S", "P", "Se", "Te", "B", "Si") and not (el == "C" and len(u) == 2 and re.fullmatch(r"H\d*", u[1])):
                continue                                   # (a condensed CH or CH2 label takes its H in too: CH2-H is CH3)
            hs = [h for h in self.nbrs(X) if h.label == "H" and len(self.nbrs(h)) == 1 and self.find(X, h).order == 1
                  and self.find(X, h).style == "plain" and not (h.charge or h.radical or h.lps) and h.uid not in used
                  and id(h) not in kid and not getattr(h, "deliberate", False)]
            if not hs: continue
            rest = u[1:]; n0 = 0
            if rest and re.fullmatch(r"H\d*", rest[0]): n0 = int(rest[0][1:] or 1); rest = rest[1:]
            if rest: continue
            n = n0 + len(hs)
            for h in hs: self.remove(h)
            hn = "H" + (str(n) if n > 1 else "")
            X.label = (hn + el) if not self.nbrs(X) else (el + hn)
        return self
    def unfold(self, atom, n=None):
        """The reverse of tidy() for one atom: a label holding H atoms (OH, NH2, H2O, OH2) becomes the bare element with
        that many explicit H atoms (n of them, if given) on free standard directions, so a protonation, a proton loss
        or an oxonium can draw its O-H bonds. Returns the new H atoms."""
        u = label_units(atom.label); el = re.sub(r"\d+'*$", "", u[0]); rest = u[1:]; k = 0
        if u and re.fullmatch(r"H\d*", u[0]) and len(u) == 2:           # H2O written H-first
            k = int(u[0][1:] or 1); el = u[1]
        elif rest and re.fullmatch(r"H\d*", rest[0]) and len(rest) == 1: k = int(rest[0][1:] or 1)
        if not k: raise ValueError(f"unfold: {atom.label!r} holds no H atoms")
        atom.label = el; out = []
        for _ in range(k if n is None else min(n, k)): out.append(self.attach(atom, None, "H"))
        if n is not None and n < k: atom.label = el + "H" + (str(k - n) if k - n > 1 else "")
        return out
    def transfer(self, atom, new_parent, angle=None, length=1.0):
        """Move `atom` (e.g. a proton) from its current partner to `new_parent`, placed on a standard direction."""
        for n in self.nbrs(atom): self.unbond(atom, n)
        if angle is None: angle = self.free_angle(new_parent, length)
        atom.xy = new_parent.xy + dirv(angle) * length; self.bond(new_parent, atom); return atom
    def remove_frag(self, atom):
        for x in self.frag(atom): self.remove(x)
    def merge(self, other):
        """Take every atom and bond of `other` (a reagent drawn as its own Mol) into this Mol, positions unchanged, so
        both can be docked, joined and validated as one step. Returns self."""
        if other is self: return self
        self._fresh = getattr(self, "_fresh", set()) | {id(a) for a in other.atoms}
        names = {a.name for a in self.atoms if a.name is not None}   # names stay unique: a merged atom whose name is taken
        for a in other.atoms:                                          # here gets a prime (a second from_smiles Mol's 8 is
            if a.name is not None and a.name in names:                 # 8'), so m[8] is still this molecule's own atom 8
                n = f"{a.name}'"
                while n in names: n += "'"
                a.name = n
            if a.name is not None: names.add(a.name)
        for a in other.atoms: a.mol = self; self.atoms.append(a)
        self.bonds += other.bonds; other.atoms, other.bonds = [], []
        return self
    def _cross_plan(self, x, ignore=()):
        """For a centre about to take a fourth substituent: if x has three bonds at about 120 deg to each other (none in a
        ring, none a wedge or hash), the turns that redraw it as a cross (bonds along 0/90/180/270 deg, each existing bond
        turned as little as possible, at most 60 deg) and the free arm of the cross: ({atom: new angle}, free angle),
        else None. A tetrahedral intermediate is drawn this way, not with the new bond squeezed between two others."""
        ign = {id(y) for y in ignore}; nb = [n for n in self.nbrs(x) if id(n) not in ign]
        if len(nb) != 3: return None
        for n in nb:
            b = self.find(x, n)
            if b.style in ("wedge", "hash", "bold", "wavy") or self.ring_of(b): return None
        th = sorted(ang(n.xy - x.xy) % 360 for n in nb)
        if not all(100 <= (th[(i + 1) % 3] - th[i]) % 360 <= 140 for i in range(3)): return None
        best = None
        for perm in itertools.permutations((0.0, 90.0, 180.0, 270.0), 3):
            rot = [adiff(t, ang(n.xy - x.xy)) for n, t in zip(nb, perm)]
            key = (round(max(abs(r) for r in rot), 1), round(sum(abs(r) for r in rot), 1), perm)
            if best is None or key < best[0]: best = (key, perm)
        if best[0][0] > 60.0: return None
        free = [t for t in (0.0, 90.0, 180.0, 270.0) if t not in best[1]][0]
        return {n: t for n, t in zip(nb, best[1])}, free
    def _relay(self, x, targets, ignore=()):
        """Turn each substituent branch of x about x so its bond to x points along targets[atom] (deg). Returns the new
        positions {id(atom): xy}, or None if a branch runs back into x (a ring) or would come within 0.9 BU of another
        atom (atoms in `ignore` left out)."""
        ign = {id(y) for y in ignore}; pos = {}; owner = {}
        for n, t in targets.items():
            br = self.branch(n, x)
            if any(y is x for y in br): return None
            r = adiff(t, ang(n.xy - x.xy))
            for y in br: pos[id(y)] = x.xy + rotv(y.xy - x.xy, r); owner[id(y)] = id(n)
        others = [y for y in self.atoms if id(y) not in ign and y is not x]
        for y in others:
            if id(y) not in pos: continue
            for z in others:
                if z is y or owner.get(id(z)) == owner[id(y)]: continue
                if norm(pos[id(y)] - pos.get(id(z), z.xy)) < 0.9: return None
        return pos
    def _apply_pos(self, x, pos):
        """Move atoms to pos (from _relay), turning their fixed lone pairs with them; recorded as a deliberate move."""
        for y in self.atoms:
            if id(y) in pos:
                r = ang(pos[id(y)] - x.xy) - ang(y.xy - x.xy) if norm(y.xy - x.xy) > 1e-9 else 0.0
                y.xy = pos[id(y)]; _turn_fixed(y, lambda v: v + r)
        self._moved = getattr(self, "_moved", set()) | {y.uid for y in self.atoms if id(y) in pos}
    @staticmethod
    def _lab_w(x):
        """Rough extra half-width (BU) of a label beyond one letter, for layout choices (Me 0.15, OMe 0.3; none for C, H, O)."""
        if not x.drawn: return 0.0
        return 0.15 * max(0, len(re.sub(r"[^A-Za-z]", "", re.sub(r"<[^>]*>", "", x.label))) - 1)
    def join(self, a, b, angle=None, order=1, length=1.0, style="plain"):
        """Make the bond a-b between two fragments of this Mol and lay the smaller fragment out on it: the new bond is
        `length` BU long along the standard free direction at the larger fragment's atom (or `angle`), and the moved
        fragment is turned so the new bond takes a standard direction at its own atom too. The moved atoms are
        recorded, so continuity() treats them as the deliberate move of this step, not drift. Returns the bond."""
        for x in (a, b):
            if x not in self.atoms: raise ValueError(f"join: {x!r} is not in this Mol (merge its molecule in first)")
        fa, fb = self.frag(a), self.frag(b)
        if any(x is b for x in fa): raise ValueError("join: the atoms are already connected; use bond() to change a bond")
        if len(fb) > len(fa): a, b, fa, fb = b, a, fb, fa                    # the smaller fragment (fb, holding b) moves
        if angle is None and len(fa) == 1:                 # a lone atom: the partner comes in from where it is (H-Br, not Br-H)
            angle = 30.0 * round(ang(b.xy - a.xy) / 30.0)
        if angle is None:                                  # a trigonal centre taking a fourth bond is redrawn as a cross
            plan = self._cross_plan(a, ignore=fb)
            pos = self._relay(a, plan[0], ignore=fb) if plan else None
            if pos is not None: self._apply_pos(a, pos); a._cross = True; angle = plan[1]
        if angle is None:                                  # a linear centre (N=C=N, O=C=O) taking a third bond: the new
            nb = [n for n in self.nbrs(a) if not any(n is y for y in fb)]   # bond square to the line, the old two turned
            if len(nb) == 2 and abs(abs(adiff(ang(nb[0].xy - a.xy), ang(nb[1].xy - a.xy))) - 180.0) < 10 and \
               all(self.find(a, n).style == "plain" and not self.ring_of(self.find(a, n)) for n in nb):   # 30 deg away to 120
                ib0 = {id(x) for x in fb}; keep_at = self.atoms
                try: self.atoms = [x for x in keep_at if id(x) not in ib0]; tn = self.free_angle(a, length)
                finally: self.atoms = keep_at
                t1, t2 = tn + 120.0, tn - 120.0; d0 = ang(nb[0].xy - a.xy)
                tg = {nb[0]: t1, nb[1]: t2} if abs(adiff(t1, d0)) <= abs(adiff(t2, d0)) else {nb[0]: t2, nb[1]: t1}
                for flip in (None, nb[0], nb[1]):            # if the turned branches would meet (both hanging on one side),
                    br = self.branch(flip, a) if flip is not None else []    # one is first mirrored across its own bond
                    if any(bb.style in ("wedge", "hash") for bb in self.bonds if any(bb.a is y for y in br)): continue
                    keep_xy = [y.xy.copy() for y in br]; u = unit(flip.xy - a.xy) if flip is not None else None
                    for y in br: v = y.xy - a.xy; y.xy = a.xy + 2 * (v @ u) * u - v
                    pos = self._relay(a, tg, ignore=fb)
                    if pos is not None: self._apply_pos(a, pos); angle = tn; break
                    for y, q in zip(br, keep_xy): y.xy = q
        ib = {id(x) for x in fb}; keep_at, keep_b = self.atoms, self.bonds
        try:                                                                   # directions judged within each fragment
            self.atoms = [x for x in keep_at if id(x) not in ib]; t_a = self.free_angle(a, length) if angle is None else float(angle)
            self.atoms = list(fb); self.bonds = [bb for bb in keep_b if id(bb.a) in ib]
            tbs = self.std_dirs(b, 1) if self.nbrs(b) else [t_a + 180.0]
        finally: self.atoms, self.bonds = keep_at, keep_b
        stereo = any(bb.style in ("wedge", "hash") for bb in self.bonds if id(bb.a) in ib)
        P0 = b.xy.copy(); target = a.xy + dirv(t_a) * length; rest = [(x.xy, self._lab_w(x)) for x in fa if x is not a]
        def place(t_b, mir):
            u = dirv(t_b); out = {}
            for x in fb:
                v = x.xy - P0
                if mir: v = 2 * (v @ u) * u - v                                 # mirror across b's new-bond axis
                out[id(x)] = target + rotv(v, (t_a + 180.0) - t_b)
            return out
        opts = [(t, mir) for t in tbs for mir in ((False, True) if not stereo and len(fb) > 2 else (False,))]
        def room(o):                                       # clearance between atoms, less what wider labels take up
            pos = place(*o); d = [norm(pos[id(x)] - q) - w - self._lab_w(x) for x in fb if x is not b for q, w in rest]
            return round(min(d + [9.0]), 2)
        t_b, mir = max(opts, key=lambda o: (min(room(o), 2.0), -opts.index(o)))
        pos = place(t_b, mir); rot = (t_a + 180.0) - t_b
        for x in fb:
            x.xy = pos[id(x)]; _turn_fixed(x, lambda v: (2 * t_b - v if mir else v) + rot)
        bd = self.bond(a, b, order, style, free=abs(length - 1) > 1e-9)
        self._moved = getattr(self, "_moved", set()) | {x.uid for x in fb}
        return bd
    def split(self, atom):
        """The fragment holding `atom` (already unbonded from the rest) taken out of this Mol as a Mol of its own,
        returned, so a leaving group or by-product stands in the scheme as its own item behind a '+'."""
        fr = self.frag(atom)
        if len(fr) == len(self.atoms): raise ValueError("split: that fragment is the whole molecule (unbond it first)")
        new = Mol(); ids = {id(x) for x in fr}
        for x in fr: self.atoms.remove(x); new.atoms.append(x); x.mol = new
        mv = [b for b in self.bonds if id(b.a) in ids]
        for b in mv: self.bonds.remove(b)
        new.bonds = mv
        if len(fr) == 1 and fr[0].ldir is None:            # a group now standing alone reads as the free species, its
            u = label_units(fr[0].label); el = self._elem(fr[0])   # anchor last (MeO-, H2O, Et3N: the charge by O, N)
            if len(u) > 1 and (el in ("O", "S") or (el == "N" and any(not re.fullmatch(r"H\d*", x) for x in u[1:]))):
                fr[0]._face = "left"
        if len(fr) == 2:                                     # a lone diatomic shows its carbon (C≡N-, C≡O), never a bare
            for x in fr:                                     # line end
                if x.label == "C": x.show = True
        return new
    TRIG = ("C", "N", "O", "B", "Si", "P", "S", "Se", "Ge", "As")
    @staticmethod
    def _elem(a):
        return re.sub(r"[^A-Za-z].*", "", label_units(a.label)[0]) if a.drawn else "C"
    def _lopsided(self, x):
        """True if x (a main-group centre, its bonds plain and in no ring) has three bonds in a T (two about opposite, the
        third about square to them) or spread even worse. Never right: a three-bond centre has its bonds 120 deg apart."""
        nb = self.nbrs(x)
        if len(nb) != 3 or getattr(x, "_tmpl", None) or self._elem(x) not in self.TRIG: return False
        for n in nb:
            b = self.find(x, n)
            if b.style != "plain" or self.ring_of(b): return False
        th = sorted(ang(n.xy - x.xy) % 360 for n in nb)
        gaps = sorted((th[(i + 1) % 3] - th[i]) % 360 for i in range(3))
        return gaps[0] <= 100 and gaps[2] >= 165
    def _trigonal(self, x, force=False):
        """Redraw a three-bond centre left as a T or spread unevenly (or a cross that lost a bond, force=True) with its
        bonds 120 deg apart: one bond kept, the other two at +-120 from it, the moving branches turned as little as
        possible, weighed by their size (a small =O or Cl turns, not the rest of the molecule). Nothing in a ring moves."""
        nb = self.nbrs(x)
        if len(nb) != 3 or not (force or self._lopsided(x)): return
        size = {id(n): len(self.branch(n, x)) for n in nb}; best = None
        for k in range(3):                                 # keep one bond, put the other two at +-120 deg from it
            t0 = ang(nb[k].xy - x.xy); rest = [n for i, n in enumerate(nb) if i != k]
            for sg in (1, -1):
                tg = {nb[k]: t0, rest[0]: t0 + sg * 120.0, rest[1]: t0 - sg * 120.0}
                rot = [abs(adiff(t, ang(n.xy - x.xy))) for n, t in tg.items()]
                key = (round(sum(r * size[id(n)] for r, n in zip(rot, tg)), 1), round(max(rot), 1), -self.find(x, nb[k]).order, k, sg)
                if best is None or key < best[0]: best = (key, tg)
        pos = self._relay(x, best[1])
        if pos is not None: self._apply_pos(x, pos)
    def release(self, atom, partner):
        """Break the bond atom-partner and return atom's fragment as its own Mol (see split). A centre left with three bonds
        that join() drew as a cross (a tetrahedral intermediate), or that is left as a T or spread unevenly (an SN1
        carbocation, the one left behind or the one departing), is redrawn with them 120 deg apart, the smaller groups
        turned."""
        self.unbond(atom, partner); out = self.split(atom)
        crossed = getattr(partner, "_cross", False); partner._cross = False
        self._trigonal(partner, crossed)
        out._trigonal(atom, getattr(atom, "_cross", False)); atom._cross = False
        return out
    def rotate(self, deg, about=None):
        c = self.center() if about is None else np.array(about, float)
        for a in self.atoms: a.xy = c + rotv(a.xy - c, deg)
        return self
    def flip(self, vertical=False, about=None):
        """Mirror left-right (or top-bottom if vertical). Swaps wedges/hashes."""
        c = self.center() if about is None else np.array(about, float)
        for a in self.atoms:
            d = a.xy - c; a.xy = c + (d * (1, -1) if vertical else d * (-1, 1))
            for l in a.lps:
                if l.angle is not None: l.angle = -l.angle if vertical else 180 - l.angle
            for k in ("charge_at", "radical_at"):             # a pinned charge or radical is mirrored with its atom
                if getattr(a, k, None) is not None: setattr(a, k, -getattr(a, k) if vertical else 180 - getattr(a, k))
        for b in self.bonds:
            if b.style in ("wedge", "hash"): b.style = "hash" if b.style == "wedge" else "wedge"
        for a in self.atoms:                              # a mirror with wedges swapped is the molecule turned over:
            if getattr(a, "z", None) is not None: a.z = -a.z   # what was in front is now behind
        return self
    def ring_of(self, b, maxn=8):
        prev = {b.a.uid: None}; frontier = [b.a]
        for _ in range(maxn - 1):
            nxt = []
            for x in frontier:
                for y in self.nbrs(x):
                    if (x is b.a and y is b.b) or y.uid in prev: continue
                    prev[y.uid] = x; nxt.append(y)
                    if y is b.b:
                        path, z = [], y
                        while z is not None: path.append(z); z = prev[z.uid]
                        return path
            frontier = nxt
        return None

# ------------------------------------------------------------------ stock reagents (always drawn the same way)
def _reagent(atoms, key):
    m = Mol()
    for a in atoms: a.mol = m; m.atoms.append(a)
    m.key = key; return m
def _lone(label, charge=0):
    m = Mol(); a = m.atom(label, (0, 0)); a.charge = charge; m.key = a; return m
def hydroxide(explicit=False):
    """Hydroxide: the single label HO⁻ (it reads ⁻OH when its partner is on the left). explicit=True draws an O-H bond,
    only for when that H itself takes part. The key atom (O) is .key; pass the Mol straight to deprotonate()/attack()."""
    if not explicit: return _lone("HO", -1)
    m = Mol(); o = m.atom("O", (0, 0)); o.charge = -1; m.attach(o, 210, "H"); m.key = o; return m
def water(explicit=True):
    """Water: H-O-H drawn with its two O-H bonds at 120 deg (key atom: O; .h = the first H, the one protonate() uses),
    or the label H2O with explicit=False (a by-product or solvent)."""
    if not explicit: return _lone("H2O")
    m = Mol(); o = m.atom("O", (0, 0)); h1 = m.attach(o, 330, "H"); m.attach(o, 210, "H"); m.key = o; m.h = h1; return m
def hydronium(explicit=True):
    """Hydronium: O⁺ with three O-H bonds (key atom: the H that is transferred, so protonate(hydronium(), acceptor)
    works directly), or the label H3O⁺ with explicit=False."""
    if not explicit: return _lone("H3O", 1)
    m = Mol(); o = m.atom("O", (0, 0)); o.charge = 1; h = m.attach(o, 90, "H"); m.attach(o, 210, "H"); m.attach(o, 330, "H")
    m.key = h; m.o = o; return m
def alkoxide(R="Et", explicit=False):
    """Alkoxide RO⁻: the single label (EtO⁻, MeO⁻, ᵗBuO⁻), or O⁻ bonded to an R label with explicit=True."""
    if not explicit: return _lone(R + "O", -1)
    m = Mol(); o = m.atom("O", (0, 0)); o.charge = -1; m.attach(o, 210, R); m.key = o; return m
def halide(X="Br"):
    """Halide ion X⁻ (Cl⁻, Br⁻, I⁻, F⁻)."""
    return _lone(X, -1)
def amine(label="Et3N"):
    """An amine base or nucleophile as one label: Et3N, iPr2NEt, Me2NH, NH3 (bonds and lone pairs attach to N)."""
    return _lone(label)
def base(label="B", charge=-1, lone_pairs=1):
    """A generic base, B:⁻ (or B: with charge=0, or "Base"): one atom drawn with its lone pair and charge, read as a
    generic symbol (no valence check; formulas count it as itself), so its conjugate acid B-H reads B-H, not boron.
    Key atom: .key. Pass it straight to deprotonate(), or draw it over an arrow (S.arrow(above=base()))."""
    m = _lone(label, charge); m.key.generic = True
    for _ in range(lone_pairs): m.key.lp()
    return m
def acid(A="Cl", charge=0):
    """A protic acid H-A drawn with its H-A bond (key atom: the H), e.g. acid("Cl"), acid("OTs")."""
    m = Mol(); a = m.atom(A, (0, 0)); a.charge = charge; h = m.attach(a, 0, "H", length=1.3); m.find(a, h).free = True
    m.key = h; m.a = a; return m                           # (1.3 BU: a bond between two labels stays readable at any angle)

# Groups contracted into one label by from_smiles: (name, SMARTS with the label's atom first, label). "auto" holds the
# groups schemes always write as labels; abbrev=True adds the ones authors usually (not always) write that way.
CONTRACT_AUTO = [
    ("NO2", "[N+](=O)[O-]", "NO2"), ("NO2", "N(=O)=O", "NO2"),
    ("OTBDPS", "[OX2][Si](c1ccccc1)(c1ccccc1)C([CH3])([CH3])[CH3]", "OTBDPS"),
    ("OTBS", "[OX2][Si]([CH3])([CH3])C([CH3])([CH3])[CH3]", "OTBS"),
    ("OTIPS", "[OX2][Si]([CH1]([CH3])[CH3])([CH1]([CH3])[CH3])[CH1]([CH3])[CH3]", "OTIPS"),
    ("OSiEt3", "[OX2][Si]([CH2][CH3])([CH2][CH3])[CH2][CH3]", "OSiEt3"),
    ("OSiMe3", "[OX2][Si]([CH3])([CH3])[CH3]", "OSiMe3"),
    ("TBDPS", "[Si](c1ccccc1)(c1ccccc1)C([CH3])([CH3])[CH3]", "TBDPS"),
    ("TBS", "[Si]([CH3])([CH3])C([CH3])([CH3])[CH3]", "TBS"),
    ("TIPS", "[Si]([CH1]([CH3])[CH3])([CH1]([CH3])[CH3])[CH1]([CH3])[CH3]", "TIPS"),
    ("SiEt3", "[Si]([CH2][CH3])([CH2][CH3])[CH2][CH3]", "SiEt3"), ("SiMe3", "[Si]([CH3])([CH3])[CH3]", "SiMe3"),
    ("SnBu3", "[Sn]([CH2][CH2][CH2][CH3])([CH2][CH2][CH2][CH3])[CH2][CH2][CH2][CH3]", "SnBu3"),
    ("SnMe3", "[Sn]([CH3])([CH3])[CH3]", "SnMe3"),
    ("Bpin", "B1OC([CH3])([CH3])C([CH3])([CH3])O1", "Bpin"), ("B(OH)2", "[BX3]([OH1])[OH1]", "B(OH)2"),
    ("OTf", "[OX2]S(=O)(=O)C(F)(F)F", "OTf"), ("OTs", "[OX2]S(=O)(=O)c1ccc([CH3])cc1", "OTs"),
    ("OMs", "[OX2]S(=O)(=O)[CH3]", "OMs"), ("ONf", "[OX2]S(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F", "ONf"),
    ("Ts", "S(=O)(=O)c1ccc([CH3])cc1", "Ts"), ("Tf", "S(=O)(=O)C(F)(F)F", "Tf"),
    ("CF3", "C(F)(F)F", "CF3"), ("CCl3", "C(Cl)(Cl)Cl", "CCl3"),
    ("OLi", "[OX2][Li]", "OLi"), ("ONa", "[OX2][Na]", "ONa"), ("OK", "[OX2][K]", "OK"), ("OMgBr", "[OX2][Mg]Br", "OMgBr"),
    ("MgBr", "[Mg]Br", "MgBr"), ("MgCl", "[Mg]Cl", "MgCl"), ("MgI", "[Mg]I", "MgI"), ("ZnCl", "[Zn]Cl", "ZnCl"),
    ("ZnBr", "[Zn]Br", "ZnBr"), ("ZnI", "[Zn]I", "ZnI")]
CONTRACT_MORE = [
    ("OMOM", "[OX2][CH2][OX2][CH3]", "OMOM"), ("OTHP", "[OX2][CH1]1[CH2][CH2][CH2][CH2][OX2]1", "OTHP"),
    ("OPMB", "[OX2][CH2]c1[cH1][cH1]c(O[CH3])[cH1][cH1]1", "OPMB"), ("OPiv", "[OX2]C(=O)C([CH3])([CH3])[CH3]", "OPiv"),
    ("CO2tBu", "C(=O)[OX2]C([CH3])([CH3])[CH3]", "CO2tBu"), ("CO2Bn", "C(=O)[OX2][CH2]c1ccccc1", "CO2Bn"),
    ("CO2Et", "C(=O)[OX2][CH2][CH3]", "CO2Et"), ("CO2Me", "C(=O)[OX2][CH3]", "CO2Me"), ("CO2H", "C(=O)[OH1]", "CO2H"),
    ("NHBoc", "[NH1]C(=O)OC([CH3])([CH3])[CH3]", "NHBoc"), ("Boc", "C(=O)OC([CH3])([CH3])[CH3]", "Boc"),
    ("NHCbz", "[NH1]C(=O)O[CH2]c1ccccc1", "NHCbz"), ("Cbz", "C(=O)O[CH2]c1ccccc1", "Cbz"),
    ("NHAc", "[NH1]C(=O)[CH3]", "NHAc"), ("OAc", "[OX2]C(=O)[CH3]", "OAc"), ("OBz", "[OX2]C(=O)c1ccccc1", "OBz"),
    ("OBn", "[OX2][CH2]c1ccccc1", "OBn"), ("OtBu", "[OX2]C([CH3])([CH3])[CH3]", "OtBu"),
    ("OEt", "[OX2][CH2][CH3]", "OEt"), ("OMe", "[OX2][CH3]", "OMe"), ("NMe2", "[NX3]([CH3])[CH3]", "NMe2"),
    ("SMe", "[SX2][CH3]", "SMe"), ("SO2Ph", "S(=O)(=O)c1ccccc1", "SO2Ph"), ("SO2Me", "S(=O)(=O)[CH3]", "SO2Me"),
    ("CN", "C#N", "CN")]
CONTRACT_OPT = [                                            # only when named in abbrev=[...]: a phenyl as text, an
    ("Ph", "c1[cH1][cH1][cH1][cH1][cH1]1", "Ph"), ("C6H5", "c1[cH1][cH1][cH1][cH1][cH1]1", "C6H5"),   # alkyl or acyl as
    ("PMB", "[CH2]c1[cH1][cH1]c(O[CH3])[cH1][cH1]1", "PMB"), ("PMP", "c1[cH1][cH1]c(O[CH3])[cH1][cH1]1", "PMP"),   # its name
    ("Bn", "[CH2]c1[cH1][cH1][cH1][cH1][cH1]1", "Bn"), ("Bz", "C(=O)c1[cH1][cH1][cH1][cH1][cH1]1", "Bz"),
    ("Piv", "C(=O)C([CH3])([CH3])[CH3]", "Piv"), ("Ac", "C(=O)[CH3]", "Ac"), ("tBu", "C([CH3])([CH3])[CH3]", "tBu"),
    ("Cy", "[CH1]1[CH2][CH2][CH2][CH2][CH2]1", "Cy"), ("iPr", "[CH1]([CH3])[CH3]", "iPr"),
    ("Bu", "[CH2;!$([CH2]([CH2,CH3])[CH2,CH3])][CH2][CH2][CH3]", "Bu"), ("nBu", "[CH2;!$([CH2]([CH2,CH3])[CH2,CH3])][CH2][CH2][CH3]", "nBu"),
    ("Pr", "[CH2;!$([CH2]([CH2,CH3])[CH2,CH3])][CH2][CH3]", "Pr"), ("nPr", "[CH2;!$([CH2]([CH2,CH3])[CH2,CH3])][CH2][CH3]", "nPr"),
    ("Et", "[CH2;!$([CH2]([CH2,CH3])[CH2,CH3])][CH3]", "Et"),   # (never the end of a longer chain) ("MOM", "[CH2][OX2][CH3]", "MOM"), ("THP", "[CH1]1[CH2][CH2][CH2][CH2][OX2]1", "THP"),
    ("NMe", "[#7X3;R][CH3]", "NMe")]                        # (a ring N-methyl: the N keeps its ring bonds)
KEEP_BONDS = {"NMe"}                                         # labels whose atom keeps its other bonds (ring bonds)
def _contract_plan(rm, abbrev, keep=()):
    """{label atom index: label} and the indices to drop, for the groups of `abbrev` ('auto', True, False/None, or a list
    of names) found in rm. A group contracts only if it hangs from the rest by one single bond from its label atom and is
    not the whole molecule; groups are tried largest first, and never overlap. A group holding an atom of `keep` (one
    the user labels) stays drawn out."""
    from rdkit import Chem
    if not abbrev: return {}, set()
    table = CONTRACT_AUTO + (CONTRACT_MORE if abbrev is True or not isinstance(abbrev, str) else [])
    if not isinstance(abbrev, (str, bool)): table = [g for g in table + CONTRACT_OPT if g[0] in set(abbrev)]
    lab, drop, used = {}, set(), set()
    for name, sma, text in table:
        q = Chem.MolFromSmarts(sma)
        for mt in rm.GetSubstructMatches(q, useChirality=False):
            g = set(mt)
            if g & used or g & set(keep) or len(g) >= rm.GetNumAtoms(): continue
            ext = [(i, n.GetIdx()) for i in mt for n in rm.GetAtomWithIdx(i).GetNeighbors() if n.GetIdx() not in g]
            if name in KEEP_BONDS:                         # the label atom keeps every bond it has outside the group
                if not ext or any(e[0] != mt[0] for e in ext): continue
            else:
                if len(ext) != 1 or ext[0][0] != mt[0]: continue
                if rm.GetBondBetweenAtoms(*ext[0]).GetBondType() != Chem.BondType.SINGLE: continue
            lab[mt[0]] = text; drop |= g - {mt[0]}; used |= g
    return lab, drop
def _layout_bad(rm, skip=()):
    """True if an RDKit 2D layout is out of shape: bonds more than 1.35 times as long as the shortest, two bonds at an atom
    under 40 deg apart, or two crossing. Atoms in skip (contracted into a label) do not count."""
    cf = rm.GetConformer(); P = cf.GetPositions()[:, :2]; skip = set(skip)
    E = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in rm.GetBonds() if not {b.GetBeginAtomIdx(), b.GetEndAtomIdx()} & skip]
    if not E: return False
    L = np.array([norm(P[i] - P[j]) for i, j in E])
    if L.min() <= 1e-6 or L.max() / L.min() > 1.35: return True
    for A in rm.GetAtoms():                                 # two bonds at an atom under 40 deg apart
        if A.GetIdx() in skip: continue
        ts = sorted(ang(P[n.GetIdx()] - P[A.GetIdx()]) % 360 for n in A.GetNeighbors() if n.GetIdx() not in skip)
        if len(ts) > 1 and min((ts[(k + 1) % len(ts)] - ts[k]) % 360 for k in range(len(ts))) < 40: return True
    return any(_seg_x(P[i], P[j], P[k], P[l]) for (i, j), (k, l) in itertools.combinations(E, 2) if not {i, j} & {k, l})
def _fit_layout(rm, alt, target):
    """rm's conformer set to alt's layout turned (and mirrored, if that fits better) onto target {index: (x, y)}."""
    Q = alt.GetConformer().GetPositions()[:, :2]; ks = sorted(target)
    Lq = [norm(Q[b.GetBeginAtomIdx()] - Q[b.GetEndAtomIdx()]) for b in alt.GetBonds()]
    if Lq: Q = Q * (1.5 / max(float(np.median(Lq)), 1e-6))   # (RDKit's own units: 1.5 a bond, as the target's)
    T = np.array([target[k] for k in ks], float); best = None
    for refl in (False, True):
        Pk = Q[ks] * ([-1, 1] if refl else [1, 1]); pc, tc = Pk.mean(0), T.mean(0)
        U, _, Vt = np.linalg.svd((Pk - pc).T @ (T - tc)); R = U @ Vt
        if np.linalg.det(R) < 0: U[:, -1] *= -1; R = U @ Vt
        err = float(np.sqrt((((Pk - pc) @ R + tc - T) ** 2).sum(1).mean()))
        if best is None or err < best[0]: best = (err, refl, R, pc, tc)
    _, refl, R, pc, tc = best; Q = (Q * ([-1, 1] if refl else [1, 1]) - pc) @ R + tc; cf = rm.GetConformer()
    for k in range(rm.GetNumAtoms()): cf.SetAtomPosition(k, (float(Q[k][0]), float(Q[k][1]), 0.0))
def _place_free(rm, pins):
    """Around a pinned drawing, RDKit can lay an unpinned substituent onto a ring bond or onto its neighbour's (the Me and
    OH of a pinned cage side by side, a C=O on a ring bond): each such branch hanging from a pinned atom is turned, as a
    whole, into the free sector round that atom with the most room (a gem pair 30 deg either side of its middle)."""
    cf = rm.GetConformer(); P = cf.GetPositions()[:, :2]; n = rm.GetNumAtoms()
    E = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in rm.GetBonds()]
    L0 = float(np.median([norm(P[i] - P[j]) for i, j in E])) if E else 1.5
    def branch(x, c):                                       # atoms reached from x without passing c
        out, todo = {x}, [x]
        while todo:
            y = todo.pop()
            for z in rm.GetAtomWithIdx(y).GetNeighbors():
                k = z.GetIdx()
                if k != c and k not in out: out.add(k); todo.append(k)
        return out
    def room(S_, c):                                        # clearance of a branch from every other atom and bond (a
        d = 9.0                                             # terminal OH, NH2 counts its H's letter too, written away
        for k in S_:                                        # from its bond)
            pts = [P[k]]
            A_ = rm.GetAtomWithIdx(k)
            if A_.GetSymbol() != "C" and A_.GetDegree() == 1 and A_.GetTotalNumHs():
                pts.append(P[k] + np.array([0.7 if P[k][0] >= P[c][0] - 1e-6 else -0.7, 0.0]) * L0 / 1.5)
            for q in pts:
                for x in range(n):
                    if x not in S_ and x != c: d = min(d, norm(q - P[x]))
                for i, j in E:
                    if i in S_ or j in S_ or c in (i, j): continue
                    a, b = P[i], P[j]; t = np.clip((q - a) @ (b - a) / max((b - a) @ (b - a), 1e-9), 0, 1)
                    d = min(d, norm(q - (a + t * (b - a))))
        return d
    for c in sorted(pins):
        A = rm.GetAtomWithIdx(c)
        fr = [nb.GetIdx() for nb in A.GetNeighbors() if nb.GetIdx() not in pins]
        fr = [x for x in fr if not (branch(x, c) & pins)]
        if not fr or len(fr) > 2: continue
        Bs = [branch(x, c) for x in fr]
        fixed = sorted(ang(P[nb.GetIdx()] - P[c]) % 360 for nb in A.GetNeighbors() if nb.GetIdx() not in fr)
        def crowded():
            for x, S_ in zip(fr, Bs):
                t = ang(P[x] - P[c])
                if room(S_, c) < 0.6 * L0 or any(abs(adiff(t, f)) < 45 for f in fixed): return True
            return len(fr) == 2 and abs(adiff(ang(P[fr[0]] - P[c]), ang(P[fr[1]] - P[c]))) < 45
        if not fixed or not crowded(): continue
        gaps = [((fixed[(k + 1) % len(fixed)] - fixed[k]) % 360 or 360, fixed[k]) for k in range(len(fixed))]
        old = P.copy(); best = None
        for g, t0 in gaps:
            if g < (70 if len(fr) == 1 else 120): continue
            mid = t0 + g / 2
            if len(fr) == 1:                                # a lone group: anywhere 45 deg clear of the bonds
                cands = [[mid]] + [[t0 + d] for d in range(45, int(g) - 44, 10)]
            else: cands = [[mid - min(30, g / 4), mid + min(30, g / 4)]]
            for ts in cands:
                order = sorted(range(len(fr)), key=lambda k: adiff(ang(old[fr[k]] - old[c]), mid))   # keep their order
                for k, t in zip(order, ts):
                    rot = math.radians(adiff(t, ang(old[fr[k]] - old[c]))); cs, sn = math.cos(rot), math.sin(rot)
                    for y in Bs[k]:
                        v = old[y] - old[c]; P[y] = old[c] + np.array([cs * v[0] - sn * v[1], sn * v[0] + cs * v[1]])
                sc = (round(min(min(room(S_, c) for S_ in Bs), 0.9 * L0), 2), -abs(adiff(ts[0], mid)), g)
                if best is None or sc > best[0]: best = (sc, P.copy())
                P[:] = old
        if best is not None and best[0][0] > min(room(S_, c) for S_ in Bs): P[:] = best[1]
    for k in range(n): cf.SetAtomPosition(k, (float(P[k][0]), float(P[k][1]), 0.0))
def _like_marks(m, ref, tol=0.15):
    """At each stereocentre a like= product shares with its drawing (same place), move the wedge or hash to the bond the
    drawing marks, so one scheme keeps one convention (an alcohol and its TBS ether both wedged on C-O); wedge or hash as
    keeps the centre's configuration (the sign of its neighbours' signed volume)."""
    def vol(c, nb, z):
        V = np.array([[*(x.xy - c.xy), z.get(x.uid, 0.0)] for x in nb])
        return np.sign(np.linalg.det(V if len(V) == 3 else np.array([V[1] - V[0], V[2] - V[0], V[3] - V[0]])))
    at = lambda p, ms: next((x for x in ms if norm(x.xy - p) < tol), None)
    for c in m.atoms:
        sb = [b for b in m.bonds if b.a is c and b.style in ("wedge", "hash")]
        nb = m.nbrs(c)
        if len(sb) != 1 or len(nb) not in (3, 4): continue
        cr = at(c.xy, ref.atoms)
        if cr is None: continue
        rb = [b for b in ref.bonds if b.a is cr and b.style in ("wedge", "hash")]
        if len(rb) != 1: continue
        y = at(rb[0].b.xy, nb)
        if y is None or y is sb[0].b or m.find(c, y).order != 1: continue
        v0 = vol(c, nb, {sb[0].b.uid: 1.0 if sb[0].style == "wedge" else -1.0}); v1 = vol(c, nb, {y.uid: 1.0})
        if not v0 or not v1: continue
        old = sb[0]; old.style = "plain"
        m.unbond(c, y); m.bond(c, y, 1, "wedge" if v1 == v0 else "hash")
def _place_ions(m):
    """The counter-ion of a salt (a charged fragment of up to three atoms) beside the opposite charge, as a chemist writes
    it: 2 BU out from that atom along its free direction (further if anything else is nearer than 1.6 BU), not where
    RDKit sets disconnected fragments side by side, 1 BU off."""
    seen, frs = set(), []
    for a in m.atoms:
        if a.uid not in seen: f = m.frag(a); seen.update(x.uid for x in f); frs.append(f)
    if len(frs) < 2: return
    main = max(frs, key=len)
    for f in frs:
        q = sum(x.charge for x in f)
        if f is main or len(f) > 3 or not q: continue
        cen = np.mean([x.xy for x in f], axis=0)
        part = [x for x in main if x.charge * q < 0] or main
        c = min(part, key=lambda x: norm(x.xy - cen))
        v = -sum((unit(n.xy - c.xy) for n in m.nbrs(c)), np.zeros(2)); d = unit(v) if norm(v) > 1e-6 else unit(cen - c.xy)
        ins = {x.uid for x in f}
        for r in (2.0, 2.5, 3.0, 3.5, 4.0):
            sh = c.xy + d * r - cen
            if all(norm(x.xy + sh - y.xy) >= 1.6 for x in f for y in m.atoms if y.uid not in ins): break
        for x in f: x.xy = x.xy + sh
def _keep_ez(rm, soft, hard=()):
    """A layout built round fixed atoms (like=, xy=) can put a free substituent of a double bond on the wrong side (the
    methyl that takes a vinyl bromide's Br trans to the chain): reflect one side of each such bond across it, so every
    E/Z drawn is the SMILES's. The side reflected holds fewer atoms laid out like= another drawing (soft) and none
    pinned by xy= (hard); a bond in a ring keeps the ring's own geometry."""
    from rdkit import Chem
    BS = Chem.BondStereo; cf = rm.GetConformer(); n = rm.GetNumAtoms()
    P = np.array([[cf.GetAtomPosition(i).x, cf.GetAtomPosition(i).y] for i in range(n)])
    def side_of(a, other):                                   # atoms on a's side of the bond a=other, a itself left out
        out, todo = set(), [x.GetIdx() for x in rm.GetAtomWithIdx(a).GetNeighbors() if x.GetIdx() != other]
        while todo:
            x = todo.pop()
            if x in out or x == a: continue
            out.add(x); todo += [y.GetIdx() for y in rm.GetAtomWithIdx(x).GetNeighbors() if y.GetIdx() != a]
        return out
    moved = False
    for b in rm.GetBonds():
        st = b.GetStereo()
        if st not in (BS.STEREOE, BS.STEREOZ, BS.STEREOCIS, BS.STEREOTRANS) or b.IsInRing(): continue
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx(); s0, s1 = list(b.GetStereoAtoms())
        u = P[j] - P[i]; c0, c1 = cross(u, P[s0] - P[i]), cross(u, P[s1] - P[i])
        if abs(c0) < 1e-6 or abs(c1) < 1e-6 or ((c0 > 0) == (c1 > 0)) == (st in (BS.STEREOZ, BS.STEREOCIS)): continue
        sides = [S_ for S_ in (side_of(j, i), side_of(i, j)) if not S_ & set(hard) and not S_ & {i, j}]
        if not sides: continue
        S_ = min(sides, key=lambda S_: (len(S_ & soft), len(S_)))
        e = unit(u)
        for k in S_:                                         # mirrored in the line through i and j
            v = P[k] - P[i]; P[k] = P[i] + 2 * (v @ e) * e - v
        moved = True
    if moved:
        for k in range(n): cf.SetAtomPosition(k, (float(P[k][0]), float(P[k][1]), 0.0))
def _members(rm, c, drop):
    """Atom indices a contracted label at c stands for (c and the dropped atoms reached from it)."""
    seen, todo = {c}, [c]
    while todo:
        x = todo.pop()
        for n in rm.GetAtomWithIdx(x).GetNeighbors():
            if n.GetIdx() in drop and n.GetIdx() not in seen: seen.add(n.GetIdx()); todo.append(n.GetIdx())
    return seen
def _bridged_systems(rm, maxring=8):
    """Atom index sets of the bridged or caged ring systems (rings of up to maxring atoms, two of them sharing three or
    more atoms, or an atom in three rings of a system: norbornane, adamantane, cubane, morphinan, a taxane core)."""
    rings = [set(r) for r in rm.GetRingInfo().AtomRings() if len(r) <= maxring]; n = len(rings)
    adj = {i: set() for i in range(n)}; br = set()
    for i, j in itertools.combinations(range(n), 2):
        k = len(rings[i] & rings[j])
        if k >= 2: adj[i].add(j); adj[j].add(i)
        if k >= 3: br |= {i, j}
    seen, out = set(), []
    for i in range(n):
        if i in seen: continue
        st, comp = [i], set()
        while st:
            x = st.pop()
            if x not in comp: comp.add(x); st += list(adj[x] - comp)
        seen |= comp; atoms = set().union(*(rings[c] for c in comp))
        if comp & br or (len(comp) >= 3 and any(sum(a in rings[c] for c in comp) >= 3 for a in atoms)): out.append(atoms)
    return out
def _depth(rm, P, seed=0xC0DE):
    """Depth of each atom of a 2D layout P (RDKit indices): a 3D conformer (same stereo) turned to fit the drawing
    best, by a proper rotation, so its z says which of two crossing bonds is in front. None if it cannot be embedded."""
    from rdkit import Chem, rdBase
    from rdkit.Chem import AllChem
    quiet = rdBase.BlockLogs()                              # (a dummy atom, R, has no force-field type: no warnings)
    try:
        rw = Chem.RWMol(rm)                                 # a dummy atom (R, NBz written as one label) embeds as carbon
        for A in rw.GetAtoms():
            if A.GetAtomicNum() == 0: A.SetAtomicNum(6); A.SetIsotope(0); A.SetNoImplicit(False)
        mh = rw.GetMol(); Chem.SanitizeMol(mh)
        mh = Chem.AddHs(mh); ps = AllChem.ETKDGv3(); ps.randomSeed = seed
        if AllChem.EmbedMolecule(mh, ps) < 0: return None
        AllChem.MMFFOptimizeMolecule(mh, maxIters=500)
    except Exception: return None
    n = rm.GetNumAtoms(); X = mh.GetConformer().GetPositions()[:n]; Y = np.column_stack([P[:n], np.zeros(n)])
    X = X - X.mean(0); Y = Y - Y.mean(0)
    U, S_, Vt = np.linalg.svd(X.T @ Y); d = np.sign(np.linalg.det(Vt.T @ U.T)) or 1.0
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T; Z = (X @ R.T)[:, 2]
    k = float(np.sqrt((Y[:, :2] ** 2).sum() / max((X ** 2).sum(), 1e-9)))     # into the drawing's units
    return Z * k
def _linear_centres(m):
    """Mark each unlabelled carbon between two collinear double bonds (allene, ketene, CO2, isocyanate, DCC) to be drawn
    as C: the two double bonds would otherwise run into one long line, and N=C=N would read N=N."""
    inc = {}
    for b in m.bonds:
        inc.setdefault(b.a.uid, []).append(b); inc.setdefault(b.b.uid, []).append(b)
    for a in m.atoms:
        a._lin = False; bs = inc.get(a.uid, [])
        if a.label not in ("C", "") or a.hidden or len(bs) != 2: continue
        if all(b.order == 2 and b.style not in ("dash", "dot") for b in bs):
            u, v = (unit((b.b if b.a is a else b.a).xy - a.xy) for b in bs)
            a._lin = bool(float(u @ v) < -0.97)                # within 14 deg of straight
def _label_room(m, a):
    """Rough box (BU) a label takes around its atom: about 0.45 BU a letter (subscript digits less), 0.62 BU tall, the
    first letter on the atom and the rest running away from its bond."""
    if not a.drawn or a.hidden: return None
    s = re.sub(r"[^A-Za-z0-9+−-]", "", a.label or ""); n = sum(0.3 if ch.isdigit() else 0.45 for ch in s)
    nb = [x for x in m.nbrs(a) if not x.hidden]
    sx = sum(unit(x.xy - a.xy)[0] for x in nb); sy = sum(unit(x.xy - a.xy)[1] for x in nb)   # the drawing's own rule
    ld = getattr(a, "ldir", None)                             # (Scene._ldir): written leftwards when the bonds come from
    if ld is None and len(nb) > 1 and abs(sx) <= 0.3 and abs(sy) > 0.5: ld = "down" if sy > 0.5 else "up"   # the right,
    if ld in ("up", "down"):                                  # the rest stacked over or under the atom between two bonds
        hw = max(0.23, (n - 0.45) / 2); up = ld == "up"       # (the H of an NH in a zigzag)
        return (a.xy[0] - hw, a.xy[1] - 0.31 - (0 if up else 0.62), a.xy[0] + hw, a.xy[1] + 0.31 + (0.62 if up else 0))
    left = (ld or ("left" if sx > 0.3 else "right" if sx < -0.3 else
            ("left" if getattr(a, "_typed_left", False) else "right"))) == "left"
    x0, x1 = (a.xy[0] - n + 0.22, a.xy[0] + 0.23) if left else (a.xy[0] - 0.23, a.xy[0] + n - 0.22)
    return (x0, a.xy[1] - 0.31, x1, a.xy[1] + 0.31)
def _relieve_labels(m):
    """Swing a terminal labelled group (OH, OMe, OAc ...) about its atom, by as little as clears it, where its label
    would land on another label or on a bond (a layout that is fine as a skeleton can still put HO onto OMe)."""
    def boxes_hit(b1, b2, pad=0.1):
        return not (b1[2] + pad < b2[0] or b2[2] + pad < b1[0] or b1[3] + pad < b2[1] or b2[3] + pad < b1[1])
    def seg_hit(bx, p, q):
        return seg_box_gap(p, q, bx) < 0.06
    def clashes(a):
        bx = _label_room(m, a)
        if bx is None: return 0
        n = 0
        for x in m.atoms:
            if x is a: continue
            bx2 = _label_room(m, x)
            if bx2 is not None and boxes_hit(bx, bx2): n += 1
        for b in m.bonds:
            if a in (b.a, b.b): continue
            if seg_hit(bx, b.a.xy, b.b.xy): n += 1
        return n
    for a in m.atoms:
        nb = m.nbrs(a)
        if len(nb) != 1 or not a.drawn or not clashes(a): continue
        if abs(unit(nb[0].xy - a.xy)[0]) <= 0.5 and len(re.sub(r"[^A-Za-z]", "", a.label or "")) > 1:
            k0 = clashes(a); a.ldir = "right" if _label_room(m, a)[0] < a.xy[0] - 0.3 else "left"  # an upright bond
            if clashes(a) < k0:                               # leaves the text free to run either way: the other way
                if not clashes(a): continue                   # (MeO for OMe) may clear it
            else: a.ldir = None
        c = nb[0]; r = norm(a.xy - c.xy); t0 = ang(a.xy - c.xy); others = [ang(x.xy - c.xy) for x in m.nbrs(c) if x is not a]
        if len(m.nbrs(c)) == 2:                               # a chain centre keeps its 120 deg: the group may only go to
            g_ = [x for x in m.nbrs(c) if x is not a][0]      # the zigzag's other side (mirrored in the centre's other bond)
            e = unit(c.xy - g_.xy); v = a.xy - g_.xy; old = a.xy.copy(); a.xy = g_.xy + 2 * (v @ e) * e - v
            if clashes(a): a.xy = old
            continue
        if len(m.nbrs(c)) < 4: continue                       # only a four-bond centre (a cross) swings a group: at a
                                                              # trigonal centre any swing breaks its 120 deg (and a ring
                                                              # atom's lone substituent keeps its bisector)
        best = None
        for d in (15, -15, 30, -30, 45, -45, 60, -60):
            t = t0 + d
            if others and min(abs(adiff(t, o)) for o in others) < 55: continue
            old = a.xy.copy(); a.xy = c.xy + dirv(t) * r
            k = clashes(a); a.xy = old
            if k == 0: best = t; break
        if best is not None: a.xy = c.xy + dirv(best) * r
# Curated skeleton layouts, for ring systems RDKit draws badly or not the way chemists do: a core SMILES (matched by
# element, any bond orders) and its atoms' positions in BU (y up), laid out as textbooks draw that skeleton. from_smiles()
# lays a molecule that contains one out on it, the rest round it; xy= and like= take precedence.
TEMPLATES = {
    "penam": ("S1CCN2C1CC2=O",                                    # beta-lactam square upright, S on top, C=O down-left
              "0.95,1.31 1.54,0.5 0.95,-0.31 0,0 0,1 -1,1 -1,0 -1.71,-0.71"),
    "epoxymorphinan": ("N1CCC23C4C1CC5=C2C(=CC=C5)OC3CC=C4",   # morphine, codeine...: A/B/C flat, the 4,5-ether and the
                       "2.10,-0.05 2.30,-1.40 1.00,-1.42 0.00,-1.00 0.87,-0.50 0.87,0.50 0.00,1.00 -0.87,0.50 -0.87,-0.50 "
                       "-1.73,-1.00 -2.60,-0.50 -2.60,0.50 -1.73,1.00 -1.12,-1.95 0.00,-2.00 0.87,-2.50 1.73,-2.00 1.73,-1.00"),
    "morphinan": ("N1CCC23C4C1CC5=C2C=CC=C5.C3CCC4",               # N-bridge in perspective over ring C (free bonds)
                  "2.10,-0.05 2.30,-1.40 1.00,-1.42 0.00,-1.00 0.87,-0.50 0.87,0.50 0.00,1.00 -0.87,0.50 -0.87,-0.50 "
                  "-1.73,-1.00 -2.60,-0.50 -2.60,0.50 -1.73,1.00 0.00,-2.00 0.87,-2.50 1.73,-2.00 1.73,-1.00"),
    "strychnos": ("C1CN2CC3=CCOC4CC(=O)N5C6C4C3CC2C61C7=CC=CC=C75",   # strychnine's ring system, flat, as Woodward drew it
                  "-1.22,1.48 -0.61,2.27 0.37,2.10 1.36,1.89 2.14,1.25 3.08,0.92 3.09,-0.08 2.47,-0.87 1.49,-1.07 0.97,-1.92 "
                  "-0.03,-1.87 -0.55,-2.73 -0.53,-1.01 -0.08,-0.12 0.91,-0.26 1.71,0.35 0.95,0.98 -0.02,1.18 -0.81,0.57 -1.71,0.12 "
                  "-2.63,0.50 -3.38,-0.16 -3.20,-1.15 -2.27,-1.52 -1.52,-0.86"),
}
_TMPL_Q = {}
def _clear_wedges(rm, prefer=()):
    """Move a stereocentre's wedge where RDKit's choice leaves its other three bonds within 180 deg (a T, which reads as
    ambiguous) to a bond that leaves them spread round it: one to a non-stereocentre first, an exocyclic one before a ring
    bond. On the atoms in prefer (a template's) a wedge joining two stereocentres moves to a bond to a non-stereocentre
    too, as IUPAC prefers. Returns the indices of bonds now wedged from their end atom (their narrow end goes there)."""
    from rdkit import Chem
    BD = Chem.BondDir; cf = rm.GetConformer(); rev = set()
    stereo = {A.GetIdx() for A in rm.GetAtoms() if A.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED}
    def gap(c, bs):
        p = cf.GetAtomPosition(c)
        ds = sorted(math.degrees(math.atan2(cf.GetAtomPosition(b.GetOtherAtomIdx(c)).y - p.y,
                                            cf.GetAtomPosition(b.GetOtherAtomIdx(c)).x - p.x)) % 360 for b in bs)
        return max((ds[(k + 1) % len(ds)] - ds[k]) % 360 for k in range(len(ds)))
    for c in sorted(stereo):
        bs = list(rm.GetAtomWithIdx(c).GetBonds())
        own = [b for b in bs if b.GetBeginAtomIdx() == c and b.GetBondDir() in (BD.BEGINWEDGE, BD.BEGINDASH)]
        if len(own) != 1: continue
        amb = len(bs) == 4 and gap(c, [b for b in bs if b.GetIdx() != own[0].GetIdx()]) >= 180
        two = c in prefer and own[0].GetOtherAtomIdx(c) in stereo
        if not (amb or two) or (amb and sum(b.IsInRing() for b in bs) == 2): continue   # (a ring atom's gem pair takes a
        cand = sorted((b for b in bs if b.GetBondDir() == BD.NONE and b.GetBondType() == Chem.BondType.SINGLE   # wedge and
                       and (len(bs) < 4 or gap(c, [x for x in bs if x.GetIdx() != b.GetIdx()]) < 179)          # a hash: below)
                       and (amb or b.GetOtherAtomIdx(c) not in stereo)),
                      key=lambda b: (b.GetOtherAtomIdx(c) in stereo, b.IsInRing(), b.GetIdx()))
        if not cand: continue
        own[0].SetBondDir(BD.NONE); b = cand[0]
        Chem.WedgeBond(b, c, cf)
        if b.GetBeginAtomIdx() != c: rev.add(b.GetIdx())
    return rev
def _template_map(rm):
    """({atom index: (x, y)} (BU), {atom indices drawn in perspective}) for the atoms of rm that make up the first
    TEMPLATES skeleton it contains, or None. A template bond not 1 BU long is a perspective bond (a bridge drawn standing
    out of the page): its atoms are drawn as a cage's are (free bonds, the depth read from a 3D conformer)."""
    from rdkit import Chem
    for name, (core, pts) in TEMPLATES.items():
        if name not in _TMPL_Q:
            p = Chem.AdjustQueryParameters(); p.makeBondsGeneric = True; p.aromatizeIfPossible = False
            p.adjustDegree = False; p.adjustRingCount = False
            cm = Chem.MolFromSmiles(core); P = [tuple(float(v) for v in c.split(",")) for c in pts.split()]
            per = {k for b in cm.GetBonds() for k in (b.GetBeginAtomIdx(), b.GetEndAtomIdx())
                   if abs(norm(np.subtract(P[b.GetBeginAtomIdx()], P[b.GetEndAtomIdx()])) - 1) > 0.03}
            _TMPL_Q[name] = (Chem.AdjustQueryProperties(cm, p), P, per)
        q, P, per = _TMPL_Q[name]
        hit = rm.GetSubstructMatch(q)
        if hit: return {i: P[k] for k, i in enumerate(hit)}, {hit[k] for k in per}
    return None
def from_smiles(smi, rotate=0.0, flip=False, align=None, at=(0, 0), labels=None, explicit_h=False, abbrev="auto", like=None,
                xy=None, wedges=None, _fh=None):
    """Mol with RDKit 2D coordinates scaled to bond length 1. align=(i, j, deg) rotates so atom i->j
    points along deg; align=(i, j, deg, k) also puts atom k on the left-hand side of the bond i->j (above it when the
    bond points right), mirroring the drawing about the bond if needed. flip=True then mirrors the drawing left-right (a
    bond aligned at deg ends up at 180 - deg). labels={atom_index:
    'Ar'} relabels atoms (e.g. dummy '*'). Atom names = RDKit index.
    abbrev: groups drawn as one label. "auto" (default) contracts those a scheme always writes as a label: NO2, silyl
    (OSiMe3, OTBS, TIPS ...), SnBu3, B(OH)2, Bpin, OTf/OTs/OMs, Ts, CF3, CCl3, OLi, MgBr ...; True also esters and
    ethers and other usual labels (CO2Me, CO2Et, CO2H, OMe, OEt, OAc, NHBoc, Boc, NHAc, NMe2, SMe, SO2Ph, CN); a list
    picks names; False draws every atom. like=other (a Mol already drawn): lay the part this molecule shares with
    `other` out exactly as it is drawn there (same orientation, same zigzag), so a product matches its starting
    material; other atoms are placed around it. xy={atom_index: (x, y)} pins atoms where a drawing puts them (BU, y up),
    for a perspective cage read off a slide or a sketch; the rest is laid out around them. Bonds between two pinned
    atoms are drawn plain (the picture's perspective carries the depth); wedges="exo" keeps the wedges of pinned
    substituents (OH, OAc on a pinned cage) and wedges=True keeps every wedge. A bridged or caged system
    takes a curated layout (TEMPLATES: penam, morphinan, strychnos) or RDKit's ring-system template where one exists
    (norbornane, adamantane, cubane); where
    bonds of the result cross, the one behind is broken, read from a 3D conformer fitted to the drawing."""
    from rdkit import Chem
    from rdkit.Chem import rdDepictor
    pins = {int(k): v for k, v in xy.items()} if xy else None
    rdDepictor.SetPreferCoordGen(True)
    ps = Chem.SmilesParserParams(); ps.removeHs = not explicit_h
    rm = Chem.MolFromSmiles(smi, ps)
    if rm is None: raise ValueError(f"from_smiles: RDKit cannot parse {smi!r}")
    if labels:
        bad = [i for i in labels if not 0 <= int(i) < rm.GetNumAtoms()]
        if bad: raise ValueError(f"from_smiles labels: no atom {bad[0]} (atom indices follow the SMILES string from 0; "
                                 f"this one has {rm.GetNumAtoms()} atoms, 0-{rm.GetNumAtoms() - 1})")
    clab, drop = _contract_plan(rm, abbrev, keep=set(labels or ()))   # (aromatic patterns match before kekulisation)
    cmap = _like_map(rm, like) if like is not None else None
    cage = set(pins) if pins else set()                     # atoms drawn as a picture has them (a template, pins): their
    if pins:                                                # bonds keep the lengths and angles of that perspective
        bad = [i for i in pins if not 0 <= i < rm.GetNumAtoms()]
        if bad: raise ValueError(f"from_smiles xy: no atom {bad[0]} (atom indices follow the SMILES string from 0)")
        cmap = dict(cmap or {}); cmap.update({i: (float(v[0]), float(v[1])) for i, v in pins.items()})
    tpl = _template_map(rm) if cmap is None else None      # a curated skeleton (TEMPLATES), drawn as textbooks draw it
    if tpl: cmap = tpl[0]; cage |= tpl[1]                   # (its perspective bonds free, as a pinned cage's)
    if cmap:
        from rdkit.Geometry import Point2D
        rdDepictor.SetPreferCoordGen(False)                  # coordMap is honoured by RDKit's own depiction
        from rdkit import rdBase
        quiet = rdBase.BlockLogs()                           # (RDKit's own banner for a failure: the error below says it)
        try: rdDepictor.Compute2DCoords(rm, coordMap={i: Point2D(float(x) * 1.5, float(y) * 1.5) for i, (x, y) in cmap.items()})
        except RuntimeError:
            if not pins: raise
            raise ValueError("from_smiles xy: RDKit cannot lay the molecule out around these pins; check the indices (they "
                             "follow the SMILES string from 0): " + ", ".join(f"{i} is {rm.GetAtomWithIdx(i).GetSymbol()}"
                             + (" (bonded to " + ", ".join(str(n.GetIdx()) for n in rm.GetAtomWithIdx(i).GetNeighbors()) + ")")
                             for i in sorted(pins))) from None
        finally: rdDepictor.SetPreferCoordGen(True); del quiet
        if not pins and _layout_bad(rm, drop):               # like= alone: where fitting the shared part to the drawing
            alt = Chem.Mol(rm); rdDepictor.Compute2DCoords(alt)   # bends the rest out of shape (a cage, a ring the drawing
            if not _layout_bad(alt, drop): _fit_layout(rm, alt, {i: (x * 1.5, y * 1.5) for i, (x, y) in cmap.items()})   # lacks),
        if pins: _place_free(rm, set(pins))                 # the molecule's own layout turned onto the shared atoms
        _keep_ez(rm, set(cmap), set(pins or ()))
    else:
        rdDepictor.Compute2DCoords(rm)                       # CoordGen first (it lays out most systems best) ...
        def quality(m):                                      # (bonds crossing, longest / shortest bond)
            c = m.GetConformer(); P = np.array([[c.GetAtomPosition(i).x, c.GetAtomPosition(i).y] for i in range(m.GetNumAtoms())])
            E = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in m.GetBonds()]
            L = np.array([norm(P[i] - P[j]) for i, j in E])
            nx = sum(_seg_x(P[i], P[j], P[k], P[l]) for (i, j), (k, l) in itertools.combinations(E, 2) if not {i, j} & {k, l})
            med = float(np.median(L)) if len(L) else 1.0     # labels that would land on one another (HO onto OMe) or
            bx = {}                                          # on a bond count as crossings too
            for A in m.GetAtoms():
                i = A.GetIdx()
                if i in drop or (A.GetSymbol() == "C" and i not in clab): continue
                lab = clab.get(i) or (A.GetSymbol() + "H" * A.GetTotalNumHs())
                ext = 0.45 * (len(re.sub(r"[^A-Za-z]", "", lab)) - 1) * med; h = 0.225 * med   # first letter on the atom,
                nbx = [P[x.GetIdx()][0] for x in A.GetNeighbors() if x.GetIdx() not in drop]  # the rest away from its bond
                left = len(nbx) == 1 and nbx[0] > P[i][0] + 0.2 * med
                bx[i] = (P[i][0] - h - (ext if left else 0), P[i][1] - 0.31 * med, P[i][0] + h + (0 if left else ext), P[i][1] + 0.31 * med)
            ks = sorted(bx)
            for x, y in itertools.combinations(ks, 2):
                a_, b_ = bx[x], bx[y]
                if not (a_[2] < b_[0] or b_[2] < a_[0] or a_[3] < b_[1] or b_[3] < a_[1]): nx += 1
            for i in ks:
                for j, k in E:
                    if i not in (j, k) and j not in drop and k not in drop and seg_box_gap(P[j], P[k], bx[i]) < 0.02 * med: nx += 1
            return (nx, round(float(L.max() / L.min()), 2) if len(L) and L.min() > 0 else 1.0)
        if rm.GetRingInfo().NumRings() > 1:                  # ... but a bridged system may come out uneven or crossed:
            q = quality(rm)                                  # RDKit's own depiction is the second opinion
            if q != (0, 1.0) and q[1] > 1.05 or q[0]:
                alt = Chem.Mol(rm); rdDepictor.SetPreferCoordGen(False)
                rdDepictor.Compute2DCoords(alt)
                qa = quality(alt)
                if qa < q: rm, q = alt, qa
                for sd in (1, 2, 3, 4, 5, 6):                # still crossed or a label on a label: RDKit's sampled
                    if not q[0]: break                       # layouts (seeded, so still deterministic)
                    alt = Chem.Mol(rm); rdDepictor.Compute2DCoords(alt, nFlipsPerSample=1, nSample=12, sampleSeed=sd)
                    qa = quality(alt)
                    if qa < q: rm, q = alt, qa
                rdDepictor.SetPreferCoordGen(True)
            sysb = _bridged_systems(rm)
            if sysb:                                         # a curated template for the cage, when RDKit has one,
                rdDepictor.SetPreferCoordGen(False)          # is the drawing chemists know (morphinan, cubane...)
                plain = Chem.Mol(rm); rdDepictor.Compute2DCoords(plain)
                tm = Chem.Mol(rm); rdDepictor.Compute2DCoords(tm, useRingTemplates=True); rdDepictor.SetPreferCoordGen(True)
                pa, pb = plain.GetConformer().GetPositions(), tm.GetConformer().GetPositions()
                hit = [S_ for S_ in sysb if np.abs(pa[list(S_)] - pb[list(S_)]).max() > 1e-3]
                if hit: rm = tm
                P_ = rm.GetConformer().GetPositions()[:, :2]
                for S_ in sysb:                              # a cage drawn in perspective (a template, or RDKit's own
                    Sl = sorted(S_); E_ = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in rm.GetBonds()   # adamantane):
                                           if b.GetBeginAtomIdx() in S_ and b.GetEndAtomIdx() in S_]          # foreshortened
                    L_ = np.array([norm(P_[i] - P_[j]) for i, j in E_]); L_ = L_ / np.median(L_)   # bonds and sharp
                    nb_ = {frozenset(e) for e in E_}                                               # angles by design,
                    far = min([norm(P_[i] - P_[j]) / np.median([norm(P_[x] - P_[y]) for x, y in E_])   # but no atom on
                               for i, j in itertools.combinations(Sl, 2) if frozenset((i, j)) not in nb_] + [9])   # another
                    if (S_ in hit or L_.max() / L_.min() > 1.15) and L_.min() >= 0.6 and far >= 0.45: cage |= set(S_)
    # a stereocentre at a ring fusion (three ring bonds and an H) shows that H on the wedge or hash, pointing out of the
    # notch opposite the fusion bond on a shorter bond (the steroid convention), rather than wedging a ring bond
    fus0 = [A.GetIdx() for A in rm.GetAtoms() if A.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED and A.GetTotalNumHs() == 1
            and A.GetDegree() == 3 and all(b.IsInRing() for b in A.GetBonds())]
    fus = [i for i in fus0 if _fh and i in _fh]               # (decided on a first build, in the final orientation)
    if fus:
        ri = rm.GetRingInfo(); cf0 = rm.GetConformer()
        P0 = np.array([[cf0.GetAtomPosition(i).x, cf0.GetAtomPosition(i).y] for i in range(rm.GetNumAtoms())])
        med = float(np.median([norm(P0[b.GetBeginAtomIdx()] - P0[b.GetEndAtomIdx()]) for b in rm.GetBonds()]))
        rh = Chem.AddHs(rm, onlyOnAtoms=fus); from rdkit.Geometry import Point3D; cfh = rh.GetConformer()
        for H in rh.GetAtoms():
            if H.GetIdx() < rm.GetNumAtoms(): continue
            c = H.GetNeighbors()[0].GetIdx(); A = rm.GetAtomWithIdx(c)
            fb = [b for b in A.GetBonds() if ri.NumBondRings(b.GetIdx()) >= 2]
            if fb: u = -unit(P0[fb[0].GetOtherAtomIdx(c)] - P0[c])
            else:
                v = -sum(unit(P0[n.GetIdx()] - P0[c]) for n in A.GetNeighbors()); u = unit(v) if norm(v) > 1e-6 else np.array([0.0, 1.0])
            q = P0[c] + u * 0.8 * med; cfh.SetAtomPosition(H.GetIdx(), Point3D(float(q[0]), float(q[1]), 0.0))
        fus_h = list(range(rm.GetNumAtoms(), rh.GetNumAtoms())); rm = rh
    else: fus_h = []
    Chem.Kekulize(rm, clearAromaticFlags=True)
    Chem.WedgeMolBonds(rm, rm.GetConformer())
    rev = _clear_wedges(rm, set(tpl[0]) if tpl else ())     # (a wedge that left its centre reading as a T, moved)
    cf = rm.GetConformer()
    xy = np.array([[cf.GetAtomPosition(i).x, cf.GetAtomPosition(i).y] for i in range(rm.GetNumAtoms())])
    if rm.GetNumBonds():
        if cmap: xy /= 1.5
        else:                                               # (a template's cage keeps its own proportions: the scale comes
            bl = [norm(xy[b.GetBeginAtomIdx()] - xy[b.GetEndAtomIdx()]) for b in rm.GetBonds()   # from the bonds outside it)
                  if not (b.GetBeginAtomIdx() in cage and b.GetEndAtomIdx() in cage)]
            xy /= np.median(bl or [norm(xy[b.GetBeginAtomIdx()] - xy[b.GetEndAtomIdx()]) for b in rm.GetBonds()])
    m = Mol(); amap = {}
    for A in rm.GetAtoms():
        if A.GetIdx() in drop: continue
        s, h = A.GetSymbol(), A.GetTotalNumHs()
        lab = "R" if s == "*" else ("C" if s == "C" else s + ("H" if h == 1 else f"H{h}" if h else ""))
        if A.GetIdx() in clab: lab = clab[A.GetIdx()]
        if labels and A.GetIdx() in labels:
            lab = labels[A.GetIdx()]; deg = sum(1 for x in A.GetNeighbors() if x.GetIdx() not in drop)
            if deg > 1 + (A.GetFormalCharge() > 0 and A.GetSymbol() in ("N", "O", "S", "P", "Se")) and lab in ONE_BOND:   # e.g.
                                                           # "Me" on a ring carbon: the index is off (indices follow
                raise ValueError(f"labels={{{A.GetIdx()}: {lab!r}}}: atom {A.GetIdx()} ({A.GetSymbol()}) has {deg} bonds, "
                                 f"but {lab} is a one-bond group; atom indices follow the SMILES string from 0")   # the SMILES)
        a = m.atom(lab, xy[A.GetIdx()], name=A.GetIdx()); amap[A.GetIdx()] = a
        a.charge = 0 if A.GetIdx() in clab else A.GetFormalCharge(); a.radical = A.GetNumRadicalElectrons() > 0
    BT = Chem.BondType
    for B in rm.GetBonds():
        i, j = B.GetBeginAtomIdx(), B.GetEndAtomIdx()
        if i in drop or j in drop: continue
        if B.GetIdx() in rev: i, j = j, i
        t = B.GetBondType(); d = B.GetBondDir()
        sty = ("wedge" if d == Chem.BondDir.BEGINWEDGE else "hash" if d == Chem.BondDir.BEGINDASH
               else "dative" if t == BT.DATIVE else "plain")
        order = {BT.DOUBLE: 2, BT.TRIPLE: 3}.get(t, 1)
        m.bond(amap[i], amap[j], order, sty)
    for k, a in amap.items():                               # a skeletal carbon's H count, while its bonds and charge are
        if a.label == "C" and rm.GetAtomWithIdx(k).GetSymbol() == "C":   # as the SMILES has them: a lone pair or radical
            a._h0 = (rm.GetAtomWithIdx(k).GetTotalNumHs(), _bsum(m, a), a.charge)   # added there that leaves no room is a slip
    for k in cage:
        if k in amap: amap[k]._tmpl = "cage"
    for b in list(m.bonds):                                 # a substituent on a cage stands at the usual 1 BU (a template
        c, x = (b.a, b.b) if getattr(b.a, "_tmpl", None) == "cage" else (b.b, b.a)   # can leave it long)
        if getattr(c, "_tmpl", None) != "cage" or getattr(x, "_tmpl", None) == "cage" or (pins and x.name in pins): continue
        seen, todo = {x.uid}, [x]
        while todo:
            y = todo.pop()
            for z in m.nbrs(y):
                if z is not c and z.uid not in seen: seen.add(z.uid); todo.append(z)
        if any(getattr(z, "_tmpl", None) == "cage" for z in m.atoms if z.uid in seen): continue
        d = unit(x.xy - c.xy) * (1.0 - norm(x.xy - c.xy))
        for z in m.atoms:
            if z.uid in seen: z.xy = z.xy + d
    for b in m.bonds:
        if getattr(b.a, "_tmpl", None) == "cage" and getattr(b.b, "_tmpl", None) == "cage":
            b.free = True
            if pins and wedges is not True and b.style in ("wedge", "hash"):          # a perspective drawing's own
                rb = rm.GetBondBetweenAtoms(int(b.a.name), int(b.b.name))                # depth; wedges="exo" keeps
                if wedges != "exo" or (rb is not None and rb.IsInRing()): b.style = "plain"   # the substituents' wedges
    if pins and wedges is not True:                         # a pinned bridgehead drawn in perspective (its ring bonds
        def _head(x):                                       # all on one side, as camphor's): the cage carries its
            if x is None or int(x.name) not in pins: return False   # configuration and a wedge on its substituent says
            ra = rm.GetAtomWithIdx(int(x.name))                    # nothing more; a flat drawing's (taxol's C1) stays
            nb = [rb.GetOtherAtomIdx(ra.GetIdx()) for rb in ra.GetBonds() if rb.IsInRing() and rb.GetOtherAtomIdx(ra.GetIdx()) in pins]
            if len(nb) < 3: return False
            ds = sorted(ang(np.asarray(pins[j], float) - np.asarray(pins[int(x.name)], float)) % 360 for j in nb)
            return max(((ds[(i + 1) % len(ds)] - ds[i]) % 360) for i in range(len(ds))) >= 180
        for b in m.bonds:
            if b.style in ("wedge", "hash") and (_head(b.a) or _head(b.b)): b.style = "plain"
    if like is not None: _like_marks(m, like)               # a like= product's wedges on the bonds its drawing has them
    for c in m.atoms:                                       # a ring atom whose two exocyclic groups sit as a gem pair, one
        ex = [b for b in m.bonds if c in (b.a, b.b) and not m.ring_of(b)]   # wedged and one plain, can leave its three
        rb = [b for b in m.bonds if c in (b.a, b.b) and m.ring_of(b)]       # plain bonds within 180 deg (a small ring):
        sb = [b for b in ex if b.style in ("wedge", "hash") and b.a is c]   # the other one takes the opposite mark, as
        pl = [b for b in ex if b.style == "plain" and b.order == 1]         # IUPAC draws two exocyclic groups
        if len(ex) != 2 or len(rb) != 2 or len(sb) != 1 or len(pl) != 1: continue
        ds = sorted(ang((b.b if b.a is c else b.a).xy - c.xy) % 360 for b in rb + pl)
        if max((ds[(k + 1) % 3] - ds[k]) % 360 for k in range(3)) < 180 - 1e-6: continue
        x = pl[0].b if pl[0].a is c else pl[0].a
        m.unbond(c, x); m.bond(c, x, 1, "hash" if sb[0].style == "wedge" else "wedge")
    E = [(b.a.xy, b.b.xy, b) for b in m.bonds]
    if any(_seg_x(p0, p1, q0, q1) for (p0, p1, b1), (q0, q1, b2) in itertools.combinations(E, 2)
           if not {b1.a, b1.b} & {b2.a, b2.b}):            # bonds cross (a cage, a pinned perspective drawing): their
        Z = _depth(rm, xy)                                  # depths say which is in front, and that one stays whole
        if Z is not None:
            for k, a in amap.items(): a.z = float(Z[k])
    if pins and not (align or rotate or flip): return m      # (like= goes on: its fusion H's as the drawing it copies)
    if align:
        for k in align[:2] + tuple(align[3:4]):
            if k not in amap:
                raise ValueError(f"from_smiles align: atom {k} " + (f"is part of an abbreviation here ({clab[next(c for c in clab if k in _members(rm, c, drop))]})"
                                 if k in drop else "is not in the molecule") + "; name a drawn atom, or pass abbrev=False")
        i, j, deg = align[:3]; m.rotate(deg - ang(amap[j].xy - amap[i].xy))
        if len(align) > 3 and cross(dirv(deg), amap[align[3]].xy - amap[i].xy) < 0:   # atom k to the left of i->j
            m.rotate(-deg); m.flip(vertical=True); m.rotate(deg)
    if rotate: m.rotate(rotate)
    if flip: m.flip()
    if _fh is None and fus0:                                  # ring-fusion stereocentres drawn as in the textbooks (one bond
        ok = []                                               # near vertical, the others at 120 deg): rebuilt with their H
        for i in fus0:
            if i not in amap: continue
            c = amap[i]; ts = sorted(ang(n.xy - c.xy) % 360 for n in m.nbrs(c))
            if len(ts) != 3: continue
            if min(min(abs(adiff(t, 90)), abs(adiff(t, 270))) for t in ts) > 12: continue
            tv = min(ts, key=lambda t: min(abs(adiff(t, 90)), abs(adiff(t, 270))))
            o1, o2 = [t for t in ts if t != tv][:2]; a_ = (o2 - o1) % 360          # the arc between the other two bonds,
            g0, w0 = (o2, 360 - a_) if (tv - o1) % 360 < a_ else (o1, a_)          # across from the vertical one: wide
            if w0 < 110: continue                                                  # enough for the H (6/6, 5/7, 4/7)
            q = c.xy + dirv(g0 + w0 / 2) * FUSION_H           # where the H would stand: clear of every other atom
            if all(norm(q - x.xy) >= 0.8 for x in m.atoms if x is not c and x not in m.nbrs(c)): ok.append(i)
        if ok:
            return from_smiles(smi, rotate, flip, align, at, labels, explicit_h, abbrev, like, pins, wedges, _fh=set(ok))
    for k in fus_h:                                           # the ring-fusion H stands straight up or down, opposite the
        if k not in amap: continue                            # centre's most nearly vertical bond (the steroid drawing)
        h = amap[k]; c = m.nbrs(h)[0]; ns = [n for n in m.nbrs(c) if n is not h]
        ts = sorted(ang(n.xy - c.xy) % 360 for n in ns); gaps = [((ts[(i + 1) % 3] - ts[i]) % 360, ts[i]) for i in range(3)]
        v = min(ns, key=lambda n: min(abs(adiff(ang(n.xy - c.xy), 90)), abs(adiff(ang(n.xy - c.xy), 270))))
        tv = ang(v.xy - c.xy) % 360                           # H across from the most vertical bond, in the middle of the
        g, t0 = next(((g, t0) for g, t0 in gaps if abs(adiff(t0, tv)) > 1 and abs(adiff((t0 + g) % 360, tv)) > 1),
                     max(gaps))                               # gap between the other two (straight up on a 6/6 fusion)
        t = (t0 + g / 2) % 360                                # the free sector's middle; within 15 deg of vertical it
        for v_ in (90.0, 270.0):                              # is drawn vertical (a 6/5 fusion's sector leans a few deg)
            if abs(adiff(t, v_)) <= 15: t = v_
        h.xy = c.xy + dirv(t) * FUSION_H; m.find(c, h).free = True   # shorter than a bond, so it fits inside a ring
    _place_ions(m); _relieve_labels(m)
    for b in m.bonds: b._smiles = True                      # (check() then knows a stretched bond is RDKit's layout)
    return m.move(*(np.array(at, float) - m.center()))
FUSION_H = 0.9    # BU from a ring-fusion centre to its H: about half a bond shows past the H's letter
def _like_map(rm, ref):
    """{atom index in rm: (x, y) in BU} for the atoms rm shares with the drawn Mol `ref` (maximum common substructure,
    elements matched, any bond orders, rings to rings), at their coordinates in `ref`; None if under 3 atoms match."""
    from rdkit import Chem
    from rdkit.Chem import rdFMCS
    rr, idx = to_rdkit(ref); inv = {v: k for k, v in idx.items()}
    res = rdFMCS.FindMCS([rr, rm], timeout=20, atomCompare=rdFMCS.AtomCompare.CompareElements,
                         bondCompare=rdFMCS.BondCompare.CompareAny, ringMatchesRingOnly=True, completeRingsOnly=True)
    if not res.smartsString or res.numAtoms < 3: return None
    q = Chem.MolFromSmarts(res.smartsString)
    A_ = rr.GetSubstructMatches(q, uniquify=True, maxMatches=20); B_ = rm.GetSubstructMatches(q, uniquify=False, maxMatches=200)
    if not A_ or not B_: return None
    qb = [(bd.GetBeginAtomIdx(), bd.GetEndAtomIdx()) for bd in q.GetBonds()]
    def agree(ab):                                          # of the ways the shared part sits in each, the pair whose bond
        a, b = ab; n = 0                                    # orders agree most: I-CH2... onto I-CH2..., not onto the
        for x, y in qb:                                     # vinyl iodide at the other end
            b1, b2 = rr.GetBondBetweenAtoms(a[x], a[y]), rm.GetBondBetweenAtoms(b[x], b[y])
            n += b1 is not None and b2 is not None and b1.GetBondType() == b2.GetBondType()
        return n
    a, b = max(((a, b) for a in A_ for b in B_), key=agree)
    byuid = {x.uid: x for x in ref.atoms}
    out = {j: tuple(byuid[inv[i]].xy) for i, j in zip(a, b) if i in inv}
    # an atom that gains substituents up to four bonds (a ring carbon becoming a gem pair, a chain carbon a cross) must
    # not keep its old substituent pinned where the lone one stood: that branch is laid out afresh with the new one
    for i, j in zip(a, b):
        if j not in out or i not in inv: continue
        A = rm.GetAtomWithIdx(j); deg = A.GetDegree()       # (bonds counted on the reference read out: a CF3 label is
        if deg < 4 or deg <= rr.GetAtomWithIdx(i).GetDegree(): continue   # a carbon with four, not one)
        ring = A.IsInRing()
        brs = []
        for nb in A.GetNeighbors():
            k = nb.GetIdx(); bd = rm.GetBondBetweenAtoms(j, k)
            if k not in out or bd.IsInRing(): continue
            seen, todo = {k}, [k]
            while todo:
                x = todo.pop()
                for y in rm.GetAtomWithIdx(x).GetNeighbors():
                    if y.GetIdx() != j and y.GetIdx() not in seen: seen.add(y.GetIdx()); todo.append(y.GetIdx())
            brs.append(seen)
        if not ring and brs: brs = [min(brs, key=len)]
        for br in brs:
            for x in br: out.pop(x, None)
    return out if len(out) >= 3 else None

GROUPS = {"Me": "C", "Et": "CC", "Pr": "CCC", "nPr": "CCC", "iPr": "C(C)C", "Bu": "CCCC", "nBu": "CCCC", "iBu": "CC(C)C",
          "sBu": "C(C)CC", "tBu": "C(C)(C)C", "Ph": "c1ccccc1", "C6H5": "c1ccccc1", "Bn": "Cc1ccccc1", "Ar": "c1ccccc1", "PMP": "c1ccc(OC)cc1",
          "Ac": "C(C)=O", "Bz": "C(=O)c1ccccc1", "Boc": "C(=O)OC(C)(C)C", "Cbz": "C(=O)OCc1ccccc1", "Ts": "S(=O)(=O)c1ccc(C)cc1",
          "Ms": "S(=O)(=O)C", "Tf": "S(=O)(=O)C(F)(F)F", "TMS": "[Si](C)(C)C", "TBS": "[Si](C)(C)C(C)(C)C",
          "OMe": "OC", "OEt": "OCC", "OAc": "OC(C)=O", "OTf": "OS(=O)(=O)C(F)(F)F", "OTs": "OS(=O)(=O)c1ccc(C)cc1",
          "OBn": "OCc1ccccc1", "OtBu": "OC(C)(C)C", "CO2H": "C(=O)O", "CO2Me": "C(=O)OC", "CO2Et": "C(=O)OCC",
          "CO2tBu": "C(=O)OC(C)(C)C", "CO2R": "C(=O)O*", "CHO": "C=O", "CN": "C#N", "NO2": "[N+](=O)[O-]",
          "NHMe": "NC", "NMe2": "N(C)C", "NHAc": "NC(C)=O", "NHBoc": "NC(=O)OC(C)(C)C", "SMe": "SC", "SiMe3": "[Si](C)(C)C",
          "CF3": "C(F)(F)F", "CCl3": "C(Cl)(Cl)Cl", "OH": "O", "NH2": "N", "NH": "N", "SH": "S",
          "TBDMS": "[Si](C)(C)C(C)(C)C", "TBDPS": "[Si](c1ccccc1)(c1ccccc1)C(C)(C)C", "TIPS": "[Si](C(C)C)(C(C)C)C(C)C",
          "TES": "[Si](CC)(CC)CC", "PMB": "Cc1ccc(OC)cc1", "MOM": "COC", "Cy": "C1CCCCC1", "THP": "C1CCCCO1", "Piv": "C(=O)C(C)(C)C",
          "Bpin": "B1OC(C)(C)C(C)(C)O1", "ONf": "OS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F", "NC": "[N+]#[C-]"}

_SUBST = {}
ONE_BOND = (set(GROUPS) - {"NH"}) | {"R", "R'", "R''", "Ar", "Hal"}   # labels that take exactly one bond (X is a halide
                                                                          # or a generic heteroatom in a ring or C=X)
GENERIC = {"R", "Ar", "X", "Y", "Z", "Nu", "E", "L", "LG", "M", "Hal", "PG", "Base", "A"}
def _subst(g):
    """Formula (Counter) of an abbreviation as a substituent: Et = C2H5, Ac = C2H3O, OTs = C7H7O3S."""
    if g not in _SUBST:
        from rdkit import Chem
        from collections import Counter
        rm = Chem.AddHs(Chem.MolFromSmiles("*" + GROUPS[g])); c = Counter()
        for A in rm.GetAtoms():
            if A.GetSymbol() != "*": c[A.GetSymbol()] += 1
        _SUBST[g] = c
    return _SUBST[g]
ELEM1 = set("H B C N O F P S K V Y I W U".split())
ACRONYMS = set("LDA DBU DBN DMAP NBS NCS NIS mCPBA MCPBA DIBAL DIBAL-H TBAF PCC PDC DMP DCC EDC HOBt DIAD DEAD TEMPO AIBN "
               "DDQ CSA PPTS TFA TFAA THF DMF DMSO DCM HMPA TMEDA KHMDS NaHMDS LiHMDS LHMDS BINAP dppf dppe dba acac cod "
               "COD NMO TPAP IBX".split())
_ORGANIC = {"H", "B", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I", "Si"}
def label_formula(label):
    """Atoms a condensed label stands for, as a Counter: 'Et3N' -> C6H15N, 'CO2Me' -> C2H3O2, 'OH' -> HO, 'BF3' -> BF3,
    'PPh3' -> C18H15P, 'LiAlH4' -> AlH4Li, 'R1' -> R. Element symbols take their written count; abbreviations (Me, Ph,
    Ts, Boc...) count as the group; generic symbols (R, X, Nu, L...) count as themselves. At each point the longest
    match wins, an abbreviation or generic symbol before an element of the same length (Ac acetyl, Ar aryl, Pr propyl)."""
    from collections import Counter
    lab = re.sub(r"<[^>]*>|[_^{}/]", "", label or "")
    if lab in ACRONYMS: return Counter({lab: 1})          # a reagent named by its acronym counts as itself (LDA, NBS)
    lab = lab.translate(str.maketrans("₀₁₂₃₄₅₆₇₈₉⁰¹²³⁴⁵⁶⁷⁸⁹ᵗⁱⁿˢ", "01234567890123456789tins"))
    lab = re.sub(r"['′″]+", "", lab); lab = re.sub(r"^(tert|sec|iso|n|i|t|s)-", lambda m_: {"tert": "t", "sec": "s", "iso": "i"}.get(m_.group(1), m_.group(1)), lab)
    literal = lambda k: all(t in _ORGANIC for t in re.findall(r"[A-Z][a-z]?", k)) and re.fullmatch(r"(?:[A-Z][a-z]?\d*)+", k)
    abbr = {k for k in GROUPS if not literal(k)}
    def at(t, i):
        best = None
        for nm in abbr | GENERIC:
            if t.startswith(nm, i) and (best is None or len(nm) > len(best[0])): best = (nm, "abbr" if nm in abbr else "gen")
        two = t[i:i + 2]
        if len(two) == 2 and two in ELEM2 and (best is None or len(best[0]) < 2): best = (two, "el")
        if best is None and t[i] in ELEM1: best = (t[i], "el")
        return best
    def parse(t, i=0, depth=0):
        c = Counter()
        while i < len(t):
            if t[i] == "(":
                sub, i = parse(t, i + 1, depth + 1); m_ = re.match(r"\d+", t[i:]); k = int(m_.group()) if m_ else 1
                i += len(m_.group()) if m_ else 0
                for e, v in sub.items(): c[e] += v * k
                continue
            if t[i] == ")" and depth: return c, i + 1
            hit = at(t, i)
            if hit is None: i += 1; continue
            nm, kind = hit; i += len(nm); m_ = re.match(r"\d+", t[i:]); k = int(m_.group()) if m_ else 1
            i += len(m_.group()) if m_ else 0
            if kind == "abbr":
                if GROUPS[nm] == "c1ccccc1" and nm == "Ar": c["Ar"] += k
                else:
                    for e, v in _subst(nm).items(): c[e] += v * k
            elif kind == "gen": c["X" if nm == "Hal" else nm] += k
            else: c[nm] += k
        return c, i
    return parse(lab)[0]
def formulas_in(text):
    """Every formula written in a piece of text (arrow conditions such as '-H2O', 'Et3N, CH2Cl2'), as Counters."""
    t = re.sub(r"<[^>]*>|[_^{}]", "", str(text)).translate(str.maketrans("₀₁₂₃₄₅₆₇₈₉", "0123456789"))
    return [label_formula(w) for w in re.findall(r"[A-Za-z0-9()]+", t) if re.search(r"[A-Z]", w)]
def hill(c):
    """A Counter of atoms as a Hill formula: C first, then H, then the rest alphabetically (generic R, X... last). Small
    hydrides are written as chemists write them: HBr, HCl, H2O, H2S, OH, NH3, PH3, BH3."""
    c = {k: v for k, v in c.items() if v > 0}
    if not c.get("C") and c.get("H") and len(c) == 2:
        (x, n), = [(k, v) for k, v in c.items() if k != "H"]; h = c["H"]
        if x in ("F", "Cl", "Br", "I", "O", "S", "Se", "Te") and n == 1:
            return ("OH" if (x == "O" and h == 1) else "H" + (str(h) if h > 1 else "") + x)
    gen = sorted(k for k in c if k in GENERIC); el = sorted(k for k in c if k not in GENERIC and k not in ("C", "H"))
    order = (["C", "H"] if c.get("C") else []) + (el if c.get("C") else sorted(el + ["H"] * bool(c.get("H")))) + gen
    if not c.get("C") and c.get("H") and len(c) == 2 and order[0] == "H":
        order = order[1:] + ["H"]                         # NH3, PH3, BH3, SiH4: the element first
    return "".join(k + (str(c[k]) if c[k] != 1 else "") for k in order if c.get(k, 0) > 0)
def _bsum(m, a):
    """Bond orders at atom a as the valence checks count them: dative, dashed and dotted bonds count nothing."""
    return sum({"dative": 0, "dash": 0, "dot": 0}.get(b.style, b.order) for b in m.bonds if a in (b.a, b.b))
def formula(ms):
    """Atoms of a Mol (or list of Mols), as (Counter, net charge): skeletal carbons take the hydrogens their valence
    leaves (one fewer for a charge, a radical or a share in a delocalised system, two fewer for a lone pair drawn on a
    neutral carbon, a carbene); labels count as written (label_formula); hidden centroids and partial bonds count
    nothing."""
    from collections import Counter
    c, q = Counter(), 0
    for m in (ms if isinstance(ms, (list, tuple)) else [ms]):
        deloc = {id(a) for d in getattr(m, "_deloc", []) for a in d["atoms"]}
        q += sum(d.get("charge", 0) for d in getattr(m, "_deloc", []))
        for a in m.atoms:
            if a.hidden:
                if getattr(a, "_float", None): c["H"] -= 1         # a floating substituent stands in for a ring H
                continue
            q += a.charge
            lab = a.label if a.drawn else "C"
            if lab == "C":
                bo = sum(b.order for b in m.bonds if a in (b.a, b.b) and b.style not in ("dash", "dot", "dative")
                         and not (getattr(b.a, "hidden", False) or getattr(b.b, "hidden", False)))
                lx = max(0, len(a.lps) + bool(getattr(a, "ambi", False)) - (1 if a.charge < 0 else 0))   # pairs beyond a
                                                                                    # carbanion's own (a carbene's, its + -)
                c["C"] += 1; c["H"] += max(0, 4 - bo - abs(a.charge) - (1 if a.radical else 0) - (1 if id(a) in deloc else 0) - 2 * lx)
            elif getattr(a, "generic", False): c[lab] += 1                 # a generic base B: is not boron
            else:
                for e, v in label_formula(lab).items(): c[e] += v
    return c, q
_VAL = dict(VALENCE, H=1, F=1, Cl=1, Br=1, I=1, Li=1, Na=1, K=1, Cs=1, Cu=1, Ag=1, Mg=2, Zn=2, Hg=2, Ca=2, Cd=2)
_OXO = {"C": 1, "N": 1, "S": 2, "Se": 1, "P": 1}            # oxo O's written straight after the atom: COMe, SO2Ph, PO(OEt)2
_GENFAM = {"R": "R", "R'": "R", "R''": "R", "X": "X", "Hal": "X"}   # generic labels stereo_report reads as a family
_LMOL = {}
def label_mol(lab, bonds=1, charge=0):
    """What a label stands for, as an RDKit Mol with the label's own (bonded) atom first, or None if it cannot be read
    (R, Nu, a name): an abbreviation (GROUPS), or a condensed label read left to right, each group on the last atom
    with a free valence and an O straight after C, S or P an oxo (OMs, OTBS, OSiMe3, OCH2Ph, CH2OH, CO2Bn, COMe, SO2Ph,
    B(OH)2, PO(OEt)2). bonds: bond orders from the drawing into the label's atom; charge: its formal charge."""
    key = (lab, bonds, charge)
    if key not in _LMOL:
        try: _LMOL[key] = _read_label(lab, bonds, charge)
        except Exception: _LMOL[key] = None                # (RDKit's sanitisation errors included)
    return _LMOL[key]
def _read_label(lab, bonds, charge):
    from rdkit import Chem
    pp = Chem.SmilesParserParams(); pp.removeHs = False
    smi = GROUPS.get(lab) or {"N3": "N=[N+]=[N-]", "NCO": "N=C=O", "NCS": "N=C=S", "OCN": "OC#N", "SCN": "SC#N"}.get(lab)
    if smi: return Chem.MolFromSmiles(smi, pp)
    try: return _read_units(label_units(lab), lab, bonds, charge)
    except Exception: pass
    hd = re.fullmatch(r"([A-Z][a-z]?)(H\d*)?([A-Z][A-Za-z]+)", lab)   # one atom and one group the units do not split
    if hd and hd.group(3) in GROUPS:                                  # out: OMOM, OTHP, NHPMB
        return _read_units([hd.group(1)] + ([hd.group(2)] if hd.group(2) else []) + [hd.group(3)], lab, bonds, charge)
    return None
def _read_units(units, lab, bonds, charge):
    from rdkit import Chem
    pp = Chem.SmilesParserParams(); pp.removeHs = False
    rw = Chem.RWMol(); st = []                               # [index, element, free valence] of atoms that take more
    def use(n):
        if not st or st[-1][2] < n: raise ValueError(lab)
        st[-1][2] -= n
        if st[-1][2] == 0: st.pop()
    def put(sub, n=1):                                       # n copies of a terminal group on the last atom with room
        for _ in range(n):
            if not st: raise ValueError(lab)
            k = rw.GetNumAtoms(); rw.InsertMol(sub); rw.AddBond(st[-1][0], k, Chem.BondType.SINGLE); use(1)
    def new(el, order=1):
        if not st: raise ValueError(lab)
        k = rw.AddAtom(Chem.Atom(el)); rw.AddBond(st[-1][0], k, Chem.BondType.DOUBLE if order == 2 else Chem.BondType.SINGLE)
        if order == 1: use(1)
        elif st[-1][1] == "C" or st[-1][1] == "N": use(2)   # (S, Se and P take their oxo O's on top of their valence)
        if order == 1 and _VAL.get(el, 1) > 1: st.append([k, el, _VAL[el] - 1])
    prev = None
    for k, u in enumerate(units):
        mm = re.fullmatch(r"(.*?)(\d*)'*", u); core, n = mm.group(1), int(mm.group(2) or 1)
        nxt = units[k + 1] if k + 1 < len(units) else ""
        if k == 0:
            if core not in ELEM1 | ELEM2 or mm.group(2) or core not in _VAL: raise ValueError(lab)
            free = _VAL[core] - bonds + (charge if core in ("N", "P", "O", "S", "Se") else -abs(charge))
            if free < 0: raise ValueError(lab)
            rw.AddAtom(Chem.Atom(core))
            if free: st.append([0, core, free])
        elif core == "H": use(n)                             # H's written on the last atom
        elif core in ("O", "(O)", "(=O)") and prev in _OXO and st and st[-1][1] == prev and not (n == 1 and nxt.startswith("H")):
            dbl = min(n, _OXO[prev]) if core == "O" else 1
            for _ in range(dbl): new("O", 2)
            for _ in range(n - dbl): new("O")
        elif core.startswith("(") and core.endswith(")"):
            sub = label_mol(core[1:-1])
            if sub is None: raise ValueError(lab)
            put(sub, n)
        elif core in GROUPS: put(Chem.MolFromSmiles(GROUPS[core], pp), n)
        elif core in GENERIC or core in ("R'", "R''"):
            d = Chem.MolFromSmiles("*"); d.GetAtomWithIdx(0).SetProp("generic", _GENFAM.get(core, core))
            d.GetAtomWithIdx(0).SetProp("gname", u.rstrip("0123456789") if not r_ident(units, k, bonds) else u)
            put(d, 1 if r_ident(units, k, bonds) else n)
        elif core in ELEM1 | ELEM2 and core in _VAL:
            if n > 1 and _VAL[core] > 1: raise ValueError(lab)
            for _ in range(n): new(core)
        else: raise ValueError(lab)
        prev = core if core in _OXO and n == 1 else None
    mol = rw.GetMol(); Chem.SanitizeMol(mol)
    return mol
def to_rdkit(m):
    """RDKit molecule (with 2D conformer and wedge directions) for a drawn Mol; returns (rdmol, {uid: rd index}). Labels
    become the atoms they stand for (label_mol); a label that cannot be read is a dummy atom with property 'unread'."""
    from rdkit import Chem
    from rdkit.Geometry import Point3D
    rw = Chem.RWMol(); idx = {}; pos = []
    for a in m.atoms:
        if a.hidden: continue                              # an eta centroid is not an atom
        lab = a.label if a.drawn else "C"
        lab = re.sub(r"</?[ib]>|[_^<>{}/]", "", lab)
        frag = fm = None
        u = label_units(lab) or ["C"]; el = re.sub(r"\d+'*$", "", u[0]); gen = el in GENERIC or el in ("R'", "R''")
        if lab in GROUPS: frag = GROUPS[lab]
        elif len(u) == 1 and el in ELEM2 | set("H B C N O F P S K V Y I W U".split()): frag = "[" + el + "]" if el not in ("C", "N", "O", "S", "P", "F", "Cl", "Br", "I", "B") else el
        elif len(u) == 2 and u[1].startswith("H") and el in ("C", "N", "O", "S", "P", "B", "Si"): frag = el
        elif gen: frag = "*"
        else:                                              # a condensed label: OMs, OTBS, CH2OH, CO2Bn ...
            bo = sum(b.order for b in m.bonds if a in (b.a, b.b) and b.style not in ("dash", "dot"))
            fm = label_mol(lab, bo, a.charge); frag = "*" if fm is None else None
        pp = Chem.SmilesParserParams(); pp.removeHs = False     # an explicit H is kept (and RDKit stays quiet about it)
        if fm is None: fm = Chem.MolFromSmiles(frag, pp) if frag != "*" else Chem.MolFromSmiles("*")
        if fm is None: fm = Chem.MolFromSmiles("*")
        base = rw.GetNumAtoms(); rw.InsertMol(fm)
        first = rw.GetAtomWithIdx(base)
        if a.charge: first.SetFormalCharge(first.GetFormalCharge() + a.charge)
        if a.radical: first.SetNumRadicalElectrons(1); first.SetNoImplicit(True)
        if first.GetAtomicNum() == 0 and a.drawn:          # read as a dummy: R, X, Nu as meant, anything else unread
            first.SetProp("molFileAlias", lab)
            if not gen: first.SetProp("unread", lab)
            elif el in _GENFAM: first.SetProp("generic", _GENFAM[el]); first.SetProp("gname", lab)
        for k in range(1, fm.GetNumAtoms()):               # the R of CO2R
            x = rw.GetAtomWithIdx(base + k)
            if x.GetAtomicNum() == 0 and not x.HasProp("generic"): x.SetProp("generic", "R"); x.SetProp("gname", "R")
        idx[a.uid] = base
        for k in range(fm.GetNumAtoms()):
            pos.append(a.xy + (dirv(37 * k) * 0.3 * k if k else 0))
            if fm.GetNumAtoms() > 1:                       # an abbreviation drawn out: remembered, so a molfile can
                rw.GetAtomWithIdx(base + k).SetProp("cd_group", f"{base}:{lab}")   # store it as a contracted group
    BT = Chem.BondType
    for b in m.bonds:
        if b.style in ("dash", "dot") or b.a.uid not in idx or b.b.uid not in idx: continue
        t = {1: BT.SINGLE, 2: BT.DOUBLE, 3: BT.TRIPLE}.get(b.order, BT.SINGLE)
        if b.style == "dative": t = BT.DATIVE
        i, j = idx[b.a.uid], idx[b.b.uid]; rw.AddBond(i, j, t); bd = rw.GetBondBetweenAtoms(i, j)
        if bd.GetBeginAtomIdx() != i: bd = rw.GetBondBetweenAtoms(i, j)
        if b.style == "wedge": bd.SetBondDir(Chem.BondDir.BEGINWEDGE)
        elif b.style == "hash": bd.SetBondDir(Chem.BondDir.BEGINDASH)
        elif b.style == "wavy": bd.SetBondDir(Chem.BondDir.UNKNOWN)
    mol = rw.GetMol(); conf = Chem.Conformer(mol.GetNumAtoms())
    for k, p in enumerate(pos): conf.SetAtomPosition(k, Point3D(float(p[0]), float(p[1]), 0.0))
    mol.AddConformer(conf, assignId=True)
    try: Chem.SanitizeMol(mol)
    except Exception: mol.UpdatePropertyCache(strict=False)
    return mol, idx

def from_molfile(src, labels=None, rotate=0.0, flip=False, align=None, at=(0, 0), contract=True):
    """Mol from a MDL molfile (a path, or the text itself; V2000 or V3000), as exported by ChemDraw: the drawing's own
    layout kept (scaled so the median bond is 1 BU), wedges, hashes and wavy bonds (narrow end at the bond's first atom),
    charges, radicals, heteroatom labels with their hydrogens (OH, NH2), explicit H atoms kept as drawn, atom aliases
    (OMe, Boc...) as labels, and abbreviations stored as contracted superatoms (SUP groups) collapsed back into one label
    (contract=False draws them out). Atom names are the molfile's atom indices from 0; labels={i: 'R'} relabels."""
    from rdkit import Chem
    text = src
    if "\n" not in str(src) and os.path.exists(str(src)): text = open(src).read()
    rm = Chem.MolFromMolBlock(text, sanitize=False, removeHs=False, strictParsing=False)
    if rm is None: raise ValueError("not a molfile that RDKit can read")
    BT = Chem.BondType
    if any(b.GetBondType() == BT.AROMATIC for b in rm.GetBonds()):
        Chem.FastFindRings(rm)
        try: Chem.Kekulize(rm, clearAromaticFlags=True)
        except Exception: pass
    try: Chem.SanitizeMol(rm, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE)
    except Exception: rm.UpdatePropertyCache(strict=False)
    cf = rm.GetConformer(); n = rm.GetNumAtoms()
    xy = np.array([[cf.GetAtomPosition(i).x, cf.GetAtomPosition(i).y] for i in range(n)], float)
    if rm.GetNumBonds():
        xy /= np.median([norm(xy[b.GetBeginAtomIdx()] - xy[b.GetEndAtomIdx()]) for b in rm.GetBonds()]) or 1.0
    drop, relab, place = set(), {}, {}
    if contract:
        for sg in Chem.GetMolSubstanceGroups(rm):
            if not sg.HasProp("TYPE") or sg.GetProp("TYPE") != "SUP" or not sg.HasProp("LABEL"): continue
            inside = set(sg.GetAtoms())
            cross = [b for b in rm.GetBonds() if (b.GetBeginAtomIdx() in inside) != (b.GetEndAtomIdx() in inside)]
            if len(cross) > 1: continue                    # a group bonded twice stays drawn out
            keep = (cross[0].GetBeginAtomIdx() if cross[0].GetBeginAtomIdx() in inside else cross[0].GetEndAtomIdx()) if cross \
                else min(inside)
            relab[keep] = sg.GetProp("LABEL"); drop |= inside - {keep}
            if not cross: place[keep] = xy[list(inside)].mean(0)
    m = Mol(); A = {}
    for at_ in rm.GetAtoms():
        i = at_.GetIdx()
        if i in drop: continue
        s = at_.GetSymbol(); h = at_.GetTotalNumHs()
        if i in relab: lab = relab[i]
        elif at_.HasProp("molFileAlias"): lab = at_.GetProp("molFileAlias")
        elif s == "C" and not at_.GetIsotope(): lab = "C"
        elif s in ("*", "R", "R#"): lab = "R"
        else: lab = (f"{at_.GetIsotope()}" if at_.GetIsotope() else "") + s + ("H" if h == 1 else f"H{h}" if h else "")
        if labels and i in labels: lab = labels[i]
        a = m.atom(lab, place.get(i, xy[i]), name=i)
        if i not in relab and not at_.HasProp("molFileAlias"):
            a.charge = at_.GetFormalCharge(); a.radical = at_.GetNumRadicalElectrons() > 0
        A[i] = a
    for b in rm.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if i not in A or j not in A: continue
        t = b.GetBondType(); st_ = b.GetPropsAsDict().get("_MolFileBondStereo", 0) or b.GetPropsAsDict().get("_MolFileBondCfg", 0)
        sty = {1: "wedge", 6: "hash", 4: "wavy"}.get(int(st_), "plain")
        if t == BT.DATIVE: sty = "dative"
        order = {BT.DOUBLE: 2, BT.TRIPLE: 3}.get(t, 1)
        m.bond(A[i], A[j], order, sty)
    if align: i, j, deg = align; m.rotate(deg - ang(A[j].xy - A[i].xy))
    if rotate: m.rotate(rotate)
    if flip: m.flip()
    return m.move(*(np.array(at, float) - m.center()))
def to_molfile(m, path=None, expand=True):
    """The drawing as an MDL molfile (V2000; returns the text, and writes it to `path` if given), for ChemDraw or any
    other program: the layout as drawn (1 BU = 1.5 molfile units, RDKit's bond length), wedges and hashes, charges and
    radicals. Abbreviations (OMe, Ph, Boc...) are drawn out as real atoms, laid out around the atom they hang from
    while every drawn atom keeps its place; expand=False writes them as labelled pseudo-atoms (aliases) instead.
    Partial bonds and eta centroids are left out."""
    from rdkit import Chem
    from rdkit.Chem import rdDepictor
    from rdkit.Geometry import Point2D
    if expand:
        rm, idx = to_rdkit(m)
        fixed = {idx[a.uid]: Point2D(*(a.xy * 1.5)) for a in m.atoms if a.uid in idx}
        rdDepictor.Compute2DCoords(rm, coordMap=fixed)
        # RDKit may not honour the fixed coordinates exactly in a bridged ring system (it can even return the mirror
        # image, which with the wedges written as drawn would be the enantiomer): fit its layout onto the drawing
        # (rotation, shift, reflection if that fits better), then put every drawn atom back exactly where it is drawn
        cf = rm.GetConformer(); ks = sorted(fixed)
        Q = np.array([[cf.GetAtomPosition(i).x, cf.GetAtomPosition(i).y] for i in range(rm.GetNumAtoms())])
        T = np.array([[fixed[i].x, fixed[i].y] for i in ks]); Pk = Q[ks]
        best = None
        for refl in (False, True):
            P = Pk * ([-1, 1] if refl else [1, 1]); pc, tc = P.mean(0), T.mean(0)
            U, _, Vt = np.linalg.svd((P - pc).T @ (T - tc)); R = U @ Vt
            if np.linalg.det(R) < 0: U[:, -1] *= -1; R = U @ Vt
            err = float(np.sqrt((((P - pc) @ R + tc - T) ** 2).sum(1).mean()))
            if best is None or err < best[0]: best = (err, refl, R, pc, tc)
        _, refl, R, pc, tc = best
        Q = (Q * ([-1, 1] if refl else [1, 1]) - pc) @ R + tc
        for i in range(rm.GetNumAtoms()):
            x, y = (fixed[i].x, fixed[i].y) if i in fixed else Q[i]
            cf.SetAtomPosition(i, (float(x), float(y), 0.0))
        rm = Chem.RWMol(rm); Chem.Kekulize(rm, clearAromaticFlags=True); BT = Chem.BondType
        for b in m.bonds:                                   # the drawn double bonds stay where they were drawn
            if b.style in ("dash", "dot") or b.a.uid not in idx or b.b.uid not in idx: continue
            bd = rm.GetBondBetweenAtoms(idx[b.a.uid], idx[b.b.uid])
            if bd is not None and b.style != "dative": bd.SetBondType({1: BT.SINGLE, 2: BT.DOUBLE, 3: BT.TRIPLE}.get(b.order, BT.SINGLE))
        rm = rm.GetMol()
    else:
        rw = Chem.RWMol(); idx = {}; BT = Chem.BondType
        for a in m.atoms:
            if a.hidden: continue
            lab = a.label if a.drawn else "C"; u = label_units(lab)
            el = re.sub(r"\d+'*$", "", u[0]) if u else "C"
            if len(u) <= 2 and el in ELEM2 | ELEM1 and (len(u) == 1 or u[1].startswith("H")): A_ = Chem.Atom(el)
            else: A_ = Chem.Atom(0); A_.SetProp("molFileAlias", lab)
            A_.SetFormalCharge(a.charge)
            if a.radical: A_.SetNumRadicalElectrons(1)
            idx[a.uid] = rw.AddAtom(A_)
        for b in m.bonds:
            if b.style in ("dash", "dot") or b.a.uid not in idx or b.b.uid not in idx: continue
            t = BT.DATIVE if b.style == "dative" else {1: BT.SINGLE, 2: BT.DOUBLE, 3: BT.TRIPLE}.get(b.order, BT.SINGLE)
            k = rw.AddBond(idx[b.a.uid], idx[b.b.uid], t) - 1; bd = rw.GetBondWithIdx(k)
            if b.style == "wedge": bd.SetBondDir(Chem.BondDir.BEGINWEDGE)
            elif b.style == "hash": bd.SetBondDir(Chem.BondDir.BEGINDASH)
            elif b.style == "wavy": bd.SetBondDir(Chem.BondDir.UNKNOWN)
        rm = rw.GetMol(); rm.UpdatePropertyCache(strict=False)
        cf = Chem.Conformer(rm.GetNumAtoms())
        for a in m.atoms:
            if a.uid in idx: cf.SetAtomPosition(idx[a.uid], (float(a.xy[0] * 1.5), float(a.xy[1] * 1.5), 0.0))
        rm.AddConformer(cf, assignId=True)
    if expand:                                              # abbreviations drawn out become contracted SUP groups,
        groups = {}                                         # as ChemDraw writes them (from_molfile collapses them)
        for at_ in rm.GetAtoms():
            if at_.HasProp("cd_group"): groups.setdefault(at_.GetProp("cd_group"), []).append(at_.GetIdx())
        for key, ids in groups.items():
            sg = Chem.CreateMolSubstanceGroup(rm, "SUP")
            for i in ids: sg.AddAtomWithIdx(i)
            for bd in rm.GetBonds():
                if (bd.GetBeginAtomIdx() in ids) != (bd.GetEndAtomIdx() in ids): sg.AddBondWithIdx(bd.GetIdx())
            sg.SetProp("LABEL", key.split(":", 1)[1])
    dirs = {}
    for b in m.bonds:                                       # the wedges as drawn, narrow end first
        if b.style in ("wedge", "hash", "wavy") and b.a.uid in idx and b.b.uid in idx:
            dirs[(idx[b.a.uid], idx[b.b.uid])] = {"wedge": 1, "hash": 6, "wavy": 4}[b.style]
    txt = Chem.MolToMolBlock(rm, kekulize=False)
    lines = txt.split("\n"); na, nb = int(lines[3][:3]), int(lines[3][3:6])
    for k in range(4 + na, 4 + na + nb):                    # write the stereo flags ourselves, in the drawn direction
        i, j = int(lines[k][:3]) - 1, int(lines[k][3:6]) - 1
        if (i, j) in dirs: lines[k] = lines[k][:9] + f"{dirs[(i, j)]:3d}" + lines[k][12:]
        elif (j, i) in dirs: lines[k] = f"{j + 1:3d}{i + 1:3d}" + lines[k][6:9] + f"{dirs[(j, i)]:3d}" + lines[k][12:]
        else: lines[k] = lines[k][:9] + "  0" + lines[k][12:]
    txt = "\n".join(lines)
    if path: open(path, "w").write(txt)
    return txt
# what a generic group is read as, at both ends of its range: R the lowest-ranking alkyl other than methyl (ethyl) and
# a very high one, X fluoride and iodide, and two different names (R, R') told apart by isotopes, first one way round,
# then the other; a descriptor that differs between the two readings depends on what the groups stand for
_FILL = ({"R": "CC", "X": "F"}, {"R": "C(C(C)(C)C)(C(C)(C)C)C(C)(C)C", "X": "I"})
def _fill_generic(mol, fill, rev=False):
    """A copy of a to_rdkit molecule with each generic dummy (R, X) turned into the group `fill` gives for its family."""
    from rdkit import Chem
    rw = Chem.RWMol(mol); cf = rw.GetConformer()
    names = sorted({a.GetProp("gname") for a in rw.GetAtoms() if a.HasProp("gname")}, reverse=rev)
    for at in list(rw.GetAtoms()):
        fam = at.GetProp("generic") if at.HasProp("generic") else None
        if at.GetAtomicNum() or fam not in fill: continue
        f = Chem.MolFromSmiles(fill[fam]); i = at.GetIdx(); p = cf.GetAtomPosition(i); mp = {0: i}
        at.SetAtomicNum(f.GetAtomWithIdx(0).GetAtomicNum()); at.SetNoImplicit(False); at.ClearProp("molFileAlias")
        g = at.GetProp("gname") if at.HasProp("gname") else fam
        if len(names) > 1 and g in names: at.SetIsotope(100 + names.index(g))
        for x in f.GetAtoms():
            if x.GetIdx(): mp[x.GetIdx()] = rw.AddAtom(Chem.Atom(x.GetAtomicNum())); cf.SetAtomPosition(mp[x.GetIdx()], p)
        for b in f.GetBonds(): rw.AddBond(mp[b.GetBeginAtomIdx()], mp[b.GetEndAtomIdx()], b.GetBondType())
    out = rw.GetMol()
    try: Chem.SanitizeMol(out)
    except Exception: out.UpdatePropertyCache(strict=False)
    return out
def stereo_report(m, quiet=False):
    """CIP labels read back from the drawing: {'atoms': {Atom: 'R'/'S'}, 'bonds': {(a, b): 'E'/'Z'}}. Labels count as the
    atoms they stand for (OMs, OTBS, CH2OH...); a generic R counts as an alkyl group and X as a halogen, and a descriptor
    that depends on which one they are is '?'. A label that cannot be read is named in the printed line."""
    from rdkit import Chem
    from rdkit.Chem import rdCIPLabeler
    mol, idx = to_rdkit(m); inv = {v: k for k, v in idx.items()}; ua = {a.uid: a for a in m.atoms}
    def cip(mol):
        Chem.AssignChiralTypesFromBondDirs(mol); Chem.DetectBondStereochemistry(mol)
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        P = mol.GetConformer().GetPositions()[:, :2]          # a double bond whose only neighbour bond on one side carries
        for bd in mol.GetBonds():                             # a stereocentre's wedge (RDKit keeps one direction a bond)
            if bd.GetBondType() != Chem.BondType.DOUBLE or bd.GetStereo() != Chem.BondStereo.STEREONONE: continue
            if any(bd.IsInRingSize(k) for k in range(3, 8)): continue   # is read straight off the drawing
            i, j = bd.GetBeginAtomIdx(), bd.GetEndAtomIdx()
            ni = [x.GetIdx() for x in mol.GetAtomWithIdx(i).GetNeighbors() if x.GetIdx() != j]
            nj = [x.GetIdx() for x in mol.GetAtomWithIdx(j).GetNeighbors() if x.GetIdx() != i]
            if not ni or not nj: continue
            u = P[j] - P[i]; c0, c1 = cross(u, P[ni[0]] - P[i]), cross(u, P[nj[0]] - P[i])
            if abs(c0) < 1e-6 or abs(c1) < 1e-6: continue
            bd.SetStereoAtoms(ni[0], nj[0])
            bd.SetStereo(Chem.BondStereo.STEREOCIS if (c0 > 0) == (c1 > 0) else Chem.BondStereo.STEREOTRANS)
        try: rdCIPLabeler.AssignCIPLabels(mol)
        except Exception: pass
        A = {at.GetIdx(): at.GetProp("_CIPCode") for at in mol.GetAtoms() if at.HasProp("_CIPCode") and at.GetIdx() in inv}
        B = {(bd.GetBeginAtomIdx(), bd.GetEndAtomIdx()): bd.GetProp("_CIPCode") for bd in mol.GetBonds()
             if bd.HasProp("_CIPCode") and bd.GetBeginAtomIdx() in inv and bd.GetEndAtomIdx() in inv}
        return A, B
    gen = any(at.HasProp("generic") and at.GetProp("generic") in _FILL[0] for at in mol.GetAtoms())
    runs = [cip(_fill_generic(mol, f, rev=k > 0)) for k, f in enumerate(_FILL)] if gen else [cip(mol)]
    out = {"atoms": {}, "bonds": {}}
    for kind, key in ((0, "atoms"), (1, "bonds")):
        for k in sorted(set().union(*(r[kind] for r in runs))):
            vs = {r[kind].get(k) for r in runs}
            lab = vs.pop() if len(vs) == 1 else "?"
            if kind == 0: out["atoms"][ua[inv[k]]] = lab
            else: out["bonds"][(ua[inv[k[0]]], ua[inv[k[1]]])] = lab
    if not quiet:
        s = ", ".join(f"{a!r}={c}" for a, c in out["atoms"].items()) or "no stereocentres"
        e = ", ".join(f"{a!r}={b!r}:{c}" for (a, b), c in out["bonds"].items())
        notes = []
        if "?" in list(out["atoms"].values()) + list(out["bonds"].values()):
            notes.append("? depends on what the R or X group stands for")
        un = sorted({at.GetProp("unread") for at in mol.GetAtoms() if at.HasProp("unread")})
        if un and (out["atoms"] or out["bonds"]):
            notes.append("label" + "s" * (len(un) > 1) + " " + ", ".join(repr(u) for u in un) + " not read as atoms: check "
                         "descriptors near " + ("them" if len(un) > 1 else "it") + " by hand")
        print("stereo:", s + (("; " + e) if e else "") + "".join(f" ({n})" for n in notes))
    return out

# ------------------------------------------------------------------ scene items
class Item:
    def move(self, dx, dy): self.x += dx; self.y += dy   # pt

class Text(Item):
    def __init__(self, s, x=0.0, y=0.0, size=None, ha="center", va="center", color=None):
        self.s, self.x, self.y, self.size, self.ha, self.va, self.color = s, x, y, size, ha, va, color

class RArrow(Item):
    """Reaction arrow. kind: '->', '<-', '<=>' (equilibrium), '<->' (resonance), '=>' and '<=' (retrosynthesis),
    '-/->' and '<-/-' (a crossed arrow: the reaction does not happen), '->>' (several steps: two short arrows, stacked
    and staggered, as Clayden draws them), '<=>>' and '<<=>' (an unequal equilibrium: the favoured half full length,
    the other half short). dashed=True draws the shaft dashed (a step that does not happen, or one only proposed).
    above/below hold text and/or reagent structures (Mol)."""
    KINDS = ("->", "<-", "<=>", "<->", "=>", "<=", "-/->", "<-/-", "->>", "<=>>", "<<=>")
    def __init__(self, above=None, below=None, kind="->", length=None, step=None, vertical=False, dashed=False, split=False):
        if kind not in self.KINDS: raise ValueError(f"arrow kind {kind!r}: use one of {', '.join(self.KINDS)}")
        self.above, self.below, self.kind, self.length, self.step = above, below, kind, length, step
        self.vertical, self.x, self.y, self.dashed, self.split = vertical, 0.0, 0.0, dashed, split
        self.branches = []; self._bmols = []         # curved arrows joining or leaving it (Scene.branch), their structures
        self._mols = [p for side in (above, below) for p in (side if isinstance(side, (list, tuple)) else [side])
                      if isinstance(p, Mol)]
        self._B = 14.4
    def move(self, dx, dy):
        self.x += dx; self.y += dy
        for m in self._mols + getattr(self, "_bmols", []): m.move(dx / self._B, dy / self._B)   # its structures travel with it

class Group:
    def __init__(self, *mols): self.mols = list(mols)

class Pointer:
    def __init__(self, targets, text, side, color, dashed, length, bend):
        self.targets, self.text, self.side, self.color, self.dashed, self.length, self.bend = targets, text, side, color, dashed, length, bend
    def mols(self):
        out = []
        for t in self.targets:
            for a in ((t.a, t.b) if isinstance(t, Mid) else t if isinstance(t, tuple) else (t,)):
                if a.mol not in out: out.append(a.mol)
        return out

class Curly:
    def __init__(self, src, dst, kind="full", side=None, bend=None, color=None, dashed=False):
        self.src, self.dst, self.kind, self.side, self.bend, self.color = src, dst, kind, side, bend, color
        self.dashed = dashed

# ------------------------------------------------------------------ scene

__all__ = [n for n in list(globals()) if not n.startswith("__")]   # every name, for the other parts
