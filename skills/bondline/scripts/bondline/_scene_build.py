"""bondline engine: Scene methods, part 1 of 3."""
# Part of the bondline engine (see __init__.py): the package's modules share one namespace, so names
# defined in any of them are visible here at call time, as in a single module.
from ._core import *  # noqa: F401,F403


class _SceneBuild:
    def __init__(self, style=None):
        self.st = style or Style(); self.F = Fonts(self.st.font, getattr(self.st, "greek", "serif"))
        if self.F.fallback and not getattr(Scene, "_warned_font", False):
            Scene._warned_font = True                          # figures from machines with different fonts will not match
            print(f"bondline: none of the fonts {', '.join(self.st.font[:-1])} is installed; using {self.F.family} "
                  "(metrics differ, so figures will not match ones made elsewhere)")
        self._cell = (self.F.cap / 2 + 0.015) * self.st.size      # clearance circle of a one-letter label (pt)
        self.mols, self.curlies, self.items, self.captions, self.extras = [], [], [], [], []; self.pointers = []
        self._rowid = 0; self._ctx = None

    # ---------- public API
    def add(self, *ms):
        for m in ms:
            if m not in self.mols: self.mols.append(m)
        return ms[0] if len(ms) == 1 else ms
    def curly(self, src, dst, kind="full", side=None, bend=None, color=None, dashed=False):
        """Electron-pushing arrow. src: LP | Atom (radical/anion) | (a, b) bond | Mid. dst: Atom | (a, b) bond | Mid.
        kind 'full' (2 e-) or 'fish' (1 e-). side 'left'/'right' of the tail->head chord forces the bulge; bend (deg)
        asks for a total turn other than the family's usual one. dashed=True draws the shaft dashed (a flow that does
        not happen, or a minor pathway)."""
        if isinstance(src, Atom) and src is dst:
            raise ValueError(f"curly: the arrow from {src!r} has it as its target too; for a free end (an arrow drawn "
                             "wrong on purpose) give a point (x, y) in BU as dst")
        c = Curly(src, dst, kind, side, bend, color, dashed); self.curlies.append(c)
        for x in (src, dst):                               # tidy() leaves the atoms an arrow uses (an H being moved)
            for a in self._xatoms(x): a._arrowed = True
        return c
    def fish(self, src, dst, **kw): return self.curly(src, dst, kind="fish", **kw)
    def arrow(self, above=None, below=None, kind="->", length=None, step=None, at=None, vertical=False, dashed=False,
              split=False):
        """Reaction arrow (kinds in RArrow.KINDS; dashed=True for a dashed shaft). above/below: a string (\\n between
        lines), a list of lines, a Mol (a reagent drawn as a structure), or a list mixing them, top to bottom:
        above=[catalyst, "(30 mol%)"]. The arrow grows to fit what it carries. A vertical arrow sets both to its right
        (text: above= over its middle and below= under it; one of them alone, or a structure, as one stack centred on
        the shaft); split=True sets above= to its left and below= to its right, each centred on the shaft (a reagent
        structure on one side, the conditions on the other)."""
        a = RArrow(above, below, kind, length, step, vertical, dashed, split); a._B = self.st.bond; self.items.append(a)
        for m in a._mols: self.add(m)
        if at is not None: a.x, a.y = at[0] * self.st.bond, at[1] * self.st.bond
        return a
    def text(self, s, at=(0, 0), size=None, ha="center", va="center", color=None, highlight=None):
        """Free text at `at` (BU); highlight= a fill behind it (a yellow tag)."""
        t = Text(s, at[0] * self.st.bond, at[1] * self.st.bond, size, ha, va, color); t.highlight = highlight
        self.items.append(t); return t
    def frame(self, items=None, pad=0.4, label=None, size=None, at=(0, 0)):
        """A rectangle round laid-out items (structures, arrows, texts: an equation boxed on an exam page), pad BU clear
        of their ink, or, with size=(w, h) in BU and no items, an empty box centred at `at` (an answer box); label= a
        bold letter in its top left corner ("A")."""
        its = [] if items is None else (list(items) if isinstance(items, (list, tuple)) else [items])
        if not its and size is None: raise ValueError("frame(): give the items to frame, or size=(w, h) for an empty box")
        self.extras.append(("frame", its, pad, label, size, at)); return its
    def plus(self): return self.text("+")
    def ring_text(self, atoms, s, color=None, size=None):
        """A label in the middle of a ring (the A, B, C, D of a steroid), given by its atoms; it moves with the ring."""
        self.add(atoms[0].mol); self.extras.append(("ringtext", list(atoms), s, color, size)); return s
    def caption(self, mol, s, gap=0.35, color=None, side="below", highlight=None):
        """Caption under a structure, or centred under several (a list: "Crossed aldol adducts" under a pair).
        color= sets its colour (a red "2 C-C bonds formed" under the bonds drawn red); side="above" sets it over the
        structure instead (a heading: "Aryl halide (electrophile)"); highlight= a fill behind it (a highlighter yellow
        "#FFFC79" behind a by-product or a step name)."""
        if side not in ("below", "above"): raise ValueError(f"caption(): side {side!r}: use 'below' or 'above'")
        ms = tuple(mol) if isinstance(mol, (list, tuple)) else mol
        for m in (ms if isinstance(ms, tuple) else (ms,)): self.add(m)
        self.captions.append([ms, s, gap, None, None, color, side, highlight]); return mol
    def orbital(self, atom, kind="p", angle=90.0, shade=None, electrons=0, size=1.0, color=None, label=None):
        """ChemDraw-style orbital on an atom. kind 'p' (two lobes), 's' (sphere) or 'lobe' (a hybrid: its lobe and a
        small back lobe of the other phase). shade: None, 'up' (lobe along angle filled; an s sphere shaded), 'down'
        (the opposite lobe; an s sphere white: the other phase). electrons: 0, 1 (radical) or 2, drawn in the lobe
        along angle. size scales the lobes (1.6 for the larger coefficient of a polarised pi bond); color sets the
        outline, and a shaded lobe is filled with a light tint of it. label names the orbital ("n_{N(sp3)}", "empty p"),
        set just past the lobe's tip in the orbital's colour. Lobes are soft teardrops (Scene.LOBE_*) meeting at a bare
        atom's nucleus; a labelled atom's come out from behind its letter, their outline open there. The pair of p
        orbitals on two atoms is their pi bond: check() counts it, so draw the sigma bond alone."""
        self.extras.append(("orbital", atom, kind, angle, shade, electrons, size, color, label)); self.add(atom.mol)
    def dot(self, atom, color="#E8261A", label=None, label_color=None):
        """A site dot on a vertex (the nucleophilic carbon in red, the electrophilic one in blue: color="#1F4FD1"), with
        an optional label beside it in the same colour ("Nu", "E"), placed clear of bonds, labels, marks and arrows like
        a charge. On a labelled atom (an H marking a C-H position) it is a halo: a disc behind the label, in a light
        colour ("#D9D9D9"). A site dot belongs to its step: copy() leaves it behind, as it does lone pairs."""
        atom.site = (color, label, label_color); self.add(atom.mol); return atom
    def _halo_r(self, a):
        """Radius (pt) of the halo disc behind a labelled atom's label (S.dot on a labelled atom)."""
        return self._cell + 2.2 * self.st.bond / 14.4
    def dipole(self, bond, color=None, side=None):
        """The crossed dipole arrow beside bond (a, b): parallel to it on its roomier side (side="left" | "right" of
        a->b to choose), pointing to the negative end b, its crossbar at the positive end."""
        a, b = bond; self.add(a.mol); self.extras.append(("dipole", a, b, color, side)); return bond
    def wiggle(self, bond, at=0.5, color=None, dashed=False, length=None):
        """The wiggly line across a bond (a, b): at its middle (at=0.5) it marks the disconnection in a retrosynthesis;
        at=1 it crosses b's free end, the attachment point of a fragment or substituent (draw the bond to a bare vertex
        b). dashed=True draws a straight dashed line instead (a bond broken). color as for other marks; length in BU
        (0.75 across a bond, 0.5 at an end)."""
        a, b = bond; self.add(a.mol)
        self.extras.append(("wiggle", a, b, float(at), color, bool(dashed), length)); return bond
    def pointer(self, targets, text=None, side="below", color=None, dashed=False, length=1.2, bend=0.0):
        """A note with arrows pointing at one or more atoms or bonds: targets an Atom, (a, b), Mid, or a list of them.
        side "below" | "above" | "left" | "right" puts the note there with straight arrows (all tails level, one note for
        several targets: "2 possible sites for enolate formation"); a number (deg) puts the note in that direction from a
        single target and, with bend=deg, draws a curved callout ("C–C bond formed"). length: arrow length (BU).
        dashed=True for the dashed pointers of teaching figures; color for arrows and text (default ink)."""
        tg = targets if isinstance(targets, list) else [targets]
        for t in tg:
            if not (isinstance(t, (Atom, Mid)) or (isinstance(t, tuple) and len(t) == 2 and all(isinstance(x, Atom) for x in t))):
                raise ValueError(f"pointer(): {t!r} is not an atom, a bond (a, b) or Mid(a, b); a note about a reaction "
                                 "arrow goes in its above= or below= text")
        p = Pointer(tg, text, side, color, dashed, length, bend); self.pointers.append(p)
        for t in tg:
            for a in ((t.a, t.b) if isinstance(t, Mid) else t if isinstance(t, tuple) else (t,)): self.add(a.mol)
        return p
    def brackets(self, mols, charge=None, pad=0.25, count=None):
        """Square brackets round one or more structures: a transition state (charge="‡"), an ion (charge="−"), or a
        repeat with count= set as a subscript at the lower right ([PhC≡CCu]₂: count=2; a polymer: count="n"). Captions
        of the bracketed structures go under the brackets, and rows space other items from them."""
        self.extras.append(("brackets", list(mols) if isinstance(mols, (list, tuple)) else [mols], charge, pad, count))

    # ---------- coordinates
    def P(self, a): return a.xy * self.st.bond

    # ---------- labels
    def _runs_w(self, runs, size):
        lay = self._lay(runs, 0.0, 0.0, size)
        return max([x + self.F.adv(t, sz, it, bd) for t, x, y, sz, it, bd in lay] or [0.0])
    def _lay(self, runs, x, y, size):
        out, prev = [], None
        for t, r, dy, it, bd in runs:
            w = self.F.adv(t, size * r, it, bd)
            if prev and r < 1 and abs(r - prev[1]) < 1e-9 and (dy > 0) != (prev[2] > 0) and not prev[5]:   # H_2^+ :
                x = prev[3]; w = max(w, prev[4]); stk = True    # stack (scripts of one base only: not the ^3 inside a
            else: stk = False                                   # subscript, as in n_{N(sp^3)}; and only two high)
            out.append((t, x, y + dy * size, size * r, it, bd)); prev = (t, r, dy, x, w, stk); x += w
        return out
    def _ubox(self, runs, x, y, w, size):
        capH = self.F.cap * size; ins = 0.06 * size
        lo = min([dy for _, _, dy, _, _ in runs] + [0]); hi = max([dy for _, _, dy, _, _ in runs] + [0])
        y0 = y + (lo * size - 0.03 * size if lo < 0 else 0)
        if any(re.search("[gjpqy]", t) for t, *_ in runs): y0 = min(y0, y - 0.2 * size)
        y1 = y + capH + (hi * size + (SUBF - 1) * capH if hi > 0 else 0)
        return (x + ins, y0, x + w - ins, y1)
    def _untangle_labels(self, m):
        """A terminal label of several letters (OH, NH2, OMe) whose letters come within 1 pt of another atom's label is
        turned to read another way (the other side, or stacked above or below its atom, never back along its own bond)
        when that leaves 1.5 pt of white to every other label and keeps its letters off other bonds; else left as it is
        (check() still reports it)."""
        st = self.st
        labs = [a for a in m.atoms if a._g and a.drawn and not a.hidden]
        def white(a, others):
            return min([box_gap(g1["box"], g2["box"]) for x in others for g1 in a._g for g2 in x._g] + [99.0])
        def on_bond(a):
            for b in m.bonds:
                if a in (b.a, b.b): continue
                if any(seg_box_gap(self.P(b.a), self.P(b.b), g["box"]) < st.lw for g in a._g): return True
            return False
        for a in labs:
            nb = m.nbrs(a)
            if len(nb) != 1 or a.ldir or len(re.sub(r"[^A-Za-z]", "", a.label or "")) < 2: continue
            others = [x for x in labs if x is not a and x is not nb[0]]
            if white(a, others) >= 1.0: continue
            v = unit(self.P(nb[0]) - self.P(a)); cur = self._ldir(a)
            into = {"right": v[0] > 0.3, "left": v[0] < -0.3, "up": v[1] > 0.5, "down": v[1] < -0.5}
            for d in ("left", "right", "up", "down"):
                if d == cur or into[d]: continue
                a.ldir = d; a._g = self._label(a); a._hE = a._hEc = a._hEw = None
                if white(a, others) >= 1.5 and not on_bond(a): break
            else:
                a.ldir = None; a._g = self._label(a); a._hE = a._hEc = a._hEw = None
    def _ldir(self, a):
        if a.ldir: return a.ldir
        typed = "left" if getattr(a, "_typed_left", False) else "right"    # HO, EtO typed: reads that way when free
        vs = [unit(n.xy - a.xy) for n in a.mol.nbrs(a)]
        if not vs: return getattr(a, "_face", None) or typed       # a lone reagent faces the partner it was placed by
        sx = sum(v[0] for v in vs); sy = sum(v[1] for v in vs)
        if sx > 0.3: return "left"
        if sx < -0.3: return "right"
        if len(vs) > 1 and sy > 0.5: return "down"
        if len(vs) > 1 and sy < -0.5: return "up"
        return typed
    def _label(self, a):
        size = self.st.size; lbl = a.label
        if "\n" in lbl: raise ValueError(f"atom label {lbl!r} contains a line break: use a caption or S.text for notes")
        lbl, it0, bd0 = _label_style(lbl)                  # <b>[M]</b>: a label set bold (or italic) as a whole
        if self.F.missing(lbl): lbl = _typo(lbl, F=self.F, minus=False)
        lit = lbl.startswith("\\")                         # \Nuc: printed as written, no subscripts, never turned round
        if lit: lbl = lbl[1:]
        elif re.fullmatch(r"(?:[A-Za-z()]+(?:_\{?\d+\}?)?)+", lbl):   # NR_2, CO_2Me, CH_3: plain formula, auto subscripts
            lbl = re.sub(r"_\{?(\d+)\}?", r"\1", lbl)
        rich = lit or "<" in lbl                           # tags inside a label: set as written (Pd^{II}, L_n are units)
        units = [lbl] if rich else label_units(lbl)
        d = self._ldir(a); pv = 0 if rich else min(getattr(a, "_pivot", 0), len(units) - 1)
        if pv: d = "pivot"                                 # EtOH, iPr2NEt: anchored in the middle, read as typed
        ident = [False] * len(units)
        if not rich and len(units) > 1 and any(re.fullmatch(r"(?:R|Ar)\d+", u) for u in units):
            m = a.mol; bo = sum(b.order for b in m.bonds if a in (b.a, b.b)) if m else 0
            onium = (a.charge or 0) > 0 or any(re.sub(r"[^A-Za-z].*", "", n.label) in METALS for n in (m.nbrs(a) if m else []))
            ident = [r_ident(units, k, bo, onium) for k in range(len(units))]
        runs = [[(lbl, 1.0, 0.0, False, False)]] if lit else \
               [parse_rich(u) if rich else unit_runs(u, first=(k == 0 and len(units) == 1), ident=ident[k]) for k, u in enumerate(units)]
        tu = label_units(_label_style(getattr(a, "_typed", ""))[0]) if d == "left" and not rich and a._typed_left else []
        if len(tu) == len(units) > 1 and _ucore(tu[-1]) == _ucore(units[0]) and not any(ident):
            units = [units[0]] + tu[-2::-1]                # typed reading leftwards ((PPh3)4Pd^0, (HO)2B, Et3N): kept as typed
            runs = [runs[0]] + [unit_runs(u, first=False) for u in units[1:]]
        elif d == "left" and not rich:                     # written leftwards, a bracketed group reads leftwards too:
            for k, u in enumerate(units):                  # B(OH)2 -> (HO)2B, P(OMe)3 -> (MeO)3P
                mm = re.fullmatch(r"\(([^()]+)\)(\d*'*)", u)
                if k and mm and len(label_units(mm.group(1))) > 1:
                    units[k] = "(" + "".join(reversed(label_units(mm.group(1)))) + ")" + mm.group(2)
                    runs[k] = unit_runs(units[k], first=False, ident=ident[k])
        if it0 or bd0: runs = [[(t, r, dy, it or it0, bd or bd0) for t, r, dy, it, bd in rr] for rr in runs]
        ws = [self._runs_w(r, size) for r in runs]
        a._dir = d if len(units) > 1 else "single"
        c = self.P(a); yb = c[1] - self.F.cap * size / 2; geo = []
        def place(r, w, x, y):
            lay = self._lay(r, x, y, size); hulls = []; sup = []
            for t, gx, gy, sz, it, bd in lay:
                box = None
                for ch in t:
                    H = self.F.hull(ch, it, bd); aw = self.F.adv(ch, sz, it, bd)
                    if H is not None and len(H) >= 3:      # a bond may pass under a trailing superscript (R¹, Pd⁰), never
                        sp = sz < size and gy > y + 0.01 and not it                                   # a prefix (ᵗBu)
                        G = H * sz + (gx, gy)
                        if sp:                             # a superscript keeps clear by its box, as in ChemDraw: a bond that
                            (bx0, by0), (bx1, by1) = G.min(0), G.max(0)     # cannot pass under it starts beyond it,
                            box = (bx0, by0, bx1, by1) if box is None else \
                                  (min(box[0], bx0), min(box[1], by0), max(box[2], bx1), max(box[3], by1))   # one box for
                        else: hulls.append(G); sup.append(False)          # the whole run (II, 2+): nothing slips between
                    gx += aw                                               # its glyphs
                if box is not None:
                    hulls.append(np.array([[box[0], box[1]], [box[2], box[1]], [box[2], box[3]], [box[0], box[3]]])); sup.append(True)
            geo.append({"runs": lay, "box": self._ubox(r, x, y, w, size), "hulls": hulls, "sup": sup, "y": y})
        x0 = c[0] - ws[0] / 2; ch = ""; core = "" if rich else _ucore(units[0])
        lay0 = [r for r in self._lay(runs[0], 0.0, 0.0, size) if abs(r[3] - size) < 1e-6 and r[0].strip()]
        if rich and not lit:                       # tags inside: anchor on the leading element symbol
            m0 = re.match(r"([A-Z][a-z]?)(?![a-z])", lbl)
            if m0 and d != "left": x0 = c[0] - self.F.adv(m0.group(1), size, it0, bd0) / 2
            el = m0.group(1) if m0 else lbl
        elif units[0][:1] in "([" and abs(sum(unit(n.xy - a.xy)[0] for n in (a.mol.nbrs(a) if a.mol else []))) <= 0.3:
            el = units[0]                          # a bracketed unit ([M]) on a bond above or below: centred
        elif not lit and (len(core) == 1 or core in ELEM2):     # an element (or R, X, L, M): centred on its letters, not on
            fr = next((r for r in lay0 if r[0].startswith(core)), None)   # a count or oxidation state (Cl2, Pd^{II}, Lₙ)
            if fr is not None: x0 = c[0] - (fr[1] + self.F.adv(core, size, fr[4], fr[5]) / 2)
            el = core
        else:                                      # a group (Me, Nu, ᵗBu, PG): its first letter, or its last when it
            t, x, _, sz, it, bd = lay0[-1] if d == "left" else lay0[0]     # reads leftwards (the e of Me)
            end = [r for r in self._lay(runs[0], 0.0, 0.0, size) if r[0].strip()][-1:]
            if d == "left" and len(units) == 1 and end and end[0][3] < size - 1e-6 and end[0][2] < 0:
                t, x, _, sz, it, bd = end[0]       # C6H13, C6H5 read leftwards: the bond meets the label's end, past its count
            ch = t[-1] if d == "left" else t[0]
            pre = self.F.adv(t[:-1], sz, it, bd) if (d == "left" and len(t) > 1) else 0.0
            x0 = c[0] - (x + pre + self.F.adv(ch, sz, it, bd) / 2); el = "x"
        x0 += self._stem_shift(a, el if el != "x" else ch, size)
        place(runs[0], ws[0], x0, yb)
        c0 = self.P(a); hs0 = geo[0]["hulls"]
        a._anchor_hull = min(hs0, key=lambda h: norm(h.mean(0) - c0)) if hs0 else None
        ach = el if el != "x" else ch        # the letter the bonds point at (the t of Et written leftwards)
        a._round = not (len(ach) == 1 and ach in NONROUND)      # I, L, T, Y, the t of Et ...: their own outline counts
        a._anchor1 = len(el) == 1 and a._round   # one-letter anchor (O, N, H, the M or e of Me): its circle stands in for it
        if d == "pivot":                           # right part as typed, then the left part ending at the anchor
            x = x0 + ws[0]
            for r, w in zip(runs[1:len(units) - pv], ws[1:len(units) - pv]): place(r, w, x, yb); x += w
            x = x0
            for r, w in reversed(list(zip(runs[len(units) - pv:], ws[len(units) - pv:]))): x -= w; place(r, w, x, yb)
        elif d == "right":
            x = x0 + ws[0]
            for r, w in zip(runs[1:], ws[1:]): place(r, w, x, yb); x += w
        elif d == "left":
            x = x0; seq, blk = [], []                  # written leftwards (R₂N, R²O₂C), but numbered groups keep
            for item in zip(runs[1:], ws[1:], ident[1:]):   # their reading order: NR1R2 becomes R¹R²N, not R²R¹N
                if item[2]: blk.append(item); continue
                seq += blk[::-1] + [item]; blk = []
            for r, w, _ in seq + blk[::-1]: x -= w; place(r, w, x, yb)
        elif len(units) > 1:
            rr = sum(runs[1:], []); w = self._runs_w(rr, size)
            dy = (self.F.cap * size + 0.12 * size) * (-1 if d == "down" else 1)
            place(rr, w, c[0] - w / 2, yb + dy)
        self._sup_nudge(a, geo)
        return geo
    OPTICAL = "LFPEK"          # letters whose ink sits well left of their box (L, F, P, E, K: 0.45-0.85 pt at 10 pt)
    def _stem_shift(self, a, ch, size):
        """On a vertical bond, centring the letter's box on the bond (as ChemDraw does) puts a lopsided letter visibly off
        the bond: the bond to an L below an atom lands in the empty space beside its stem. Centre the letter's ink
        instead (its area centroid), so L, F, P, E and K sit on a vertical bond the way O and N do, the same above and
        below. Returns the shift (pt, along x)."""
        if ch not in self.OPTICAL or a.mol is None: return 0.0
        A = self.P(a)
        if not any(abs(v[0]) < 0.5 for v in (unit(self.P(n) - A) for n in a.mol.nbrs(a))): return 0.0
        if not hasattr(self, "_inkc"): self._inkc = {}
        if ch not in self._inkc:
            A_ = Cx = 0.0
            for Q in TextPath((0, 0), ch, size=100, prop=self.F.prop()).to_polygons():
                x, y = Q[:, 0], Q[:, 1]; x1, y1 = np.roll(x, -1), np.roll(y, -1); cr = x * y1 - x1 * y
                A_ += cr.sum() / 2; Cx += ((x + x1) * cr).sum() / 6          # holes wind the other way and subtract
            self._inkc[ch] = self.F.adv(ch, 1.0) / 2 - (Cx / A_ / 100.0 if A_ else self.F.adv(ch, 1.0) / 2)
        return self._inkc[ch] * size
    WAVE = (0.14, 0.055, 4)  # wavy bond: half-wave length and amplitude (BU), fewest half-waves
    SITE_R = 1.9                  # a site dot (nucleophilic/electrophilic carbon) is 4.2 pt (0.29 bond) across, as in teaching figures
    MARK_W, RING_W = 0.72, 0.63   # stroked marks (the signs in ⊕ ⊖, the ‡) carry a letter's stem weight at their size
    DAGGER_W = 0.72        # (0.79 pt: a 10 pt stem is 0.93 pt, a superscript's 0.79), the charge ring a touch lighter
    SEAM = 0.12             # pt: hairline round each polygon of a perspective ring's outline, so neighbours meet without a seam
    CAP_SEP = 1.0           # BU: captions side by side in a row keep this much white, so they never read as one line
    TEXT_CLEAR = 5.0        # pt (at 14.4 pt bonds): a note, caption or free text keeps this much white from the ink of a
                            # structure it does not belong to (a third of a bond), so it never reads as part of it
    SUP_GAP = 1.05          # pt of white between a bond and a superscript it passes: the same as for any other letter
    SUP_NUDGE = 2.0         # pt: the furthest a superscript moves aside for a bond before the bond goes round it instead
    def _sup_nudge(self, a, geo):
        """A bond leaving a label steeply up and to the right (R³ with its bond at about 58-70 deg) runs past the
        corner of the superscript. Move the superscript right, as little as possible and by at most SUP_NUDGE, until
        every bond at the atom passes it with SUP_GAP of white. A shallower bond, level ones included, cannot get past
        that way: the superscript stays put and the bond starts beyond it, clear of the label like any other letter
        (the label reads as one block, as in ChemDraw)."""
        m = a.mol
        if m is None or not any(any(g["sup"]) for g in geo): return
        nb = m.nbrs(a)
        if not nb: return
        k = self.st.bond / 14.4; A = self.P(a); R0 = self._cell + self.st.margin
        need = self.st.lw / 2 + self.SUP_GAP * k
        segs = np.array([[*(A + unit(self.P(n) - A) * R0), *self.P(n)] for n in nb if norm(self.P(n) - A) > R0])
        if not len(segs): return
        size = self.st.size
        for g in geo:
            if not any(g["sup"]): continue
            ks = [i for i, (t, gx, gy, sz, it, bd) in enumerate(g["runs"]) if sz < size and gy > g["y"] + 0.01 and not it]
            if not ks: continue
            r0 = ks[0]; xs = g["runs"][r0][1]
            if any(abs(r[1] - xs) < 1e-6 for r in g["runs"][:r0]): continue      # stacked over a subscript (H₂⁺): leave it
            moving = [i for i, h in enumerate(g["hulls"]) if h[:, 0].min() >= xs - 0.05 * size]   # the superscript and all after it
            sh = [g["hulls"][i] for i, s in enumerate(g["sup"]) if s and i in moving]
            if not sh: continue
            gap = lambda dx: min(float(_segs_dist(h + (dx, 0.0), segs).min()) for h in sh)
            if gap(0.0) >= need: continue
            for dx in np.arange(0.25, self.SUP_NUDGE + 1e-9, 0.25) * k:
                if gap(dx) >= need:
                    g["runs"] = [(t, gx + dx, gy, sz, it, bd) if i >= r0 else (t, gx, gy, sz, it, bd)
                                 for i, (t, gx, gy, sz, it, bd) in enumerate(g["runs"])]
                    g["hulls"] = [h + (dx, 0.0) if i in moving else h for i, h in enumerate(g["hulls"])]
                    b = g["box"]; g["box"] = (b[0], b[1], b[2] + dx, b[3])
                    break
    def _exit(self, a, u, origin=None, margin=None, cell=False):
        """Distance along u (from the atom, or from `origin`) at which a line or mark clears the label by `margin`.
        The bonded letter counts as a circle one cap-height across, centred on the atom (like ChemDraw's character
        cell), so bond ends, lone pairs, charges and arrow tips sit at the same distance in every direction whatever
        the letter's shape (H, O, the small e of Me); the label's other letters count with their real outlines.
        Bonds and arrows meet a narrow letter (I, L, the t of Et: NONROUND) at its own outline instead, which the circle
        would leave far behind; marks (cell=True: lone pairs, charges, radicals) keep the circle round every letter,
        so the dots of I⁻ sit in an even ring as they do round O⁻."""
        if a.clear is not None: return a.clear
        g = a._g
        if not g: return 0.0
        m = self.st.margin if margin is None else margin
        o = self.P(a) if origin is None else origin
        w = self.P(a) - o; along = float(w @ u); off2 = float(w @ w) - along * along; R = self._cell + m
        rnd = cell or getattr(a, "_round", True)
        disc = along + math.sqrt(R * R - off2) if (off2 < R * R and rnd) else 0.0
        hx = self._exit_hull(a, u, o, m, disc, cell=cell and not getattr(a, "_round", True))
        if hx > disc: self._limited = True                      # this end is set by another letter of the label
        return max(disc, hx)
    def _exit_hull(self, a, u, o, m, disc=0.0, cell=False, whole=False):
        """Where a ray from o along u, starting inside the label, leaves it for good: clear of the clearance circle
        (up to `disc`) and `m` away from the outlines of the label's other letters (a bond may pass under a
        superscript); pt. whole=True measures every letter, the anchor too (an orbital lobe leaves the O itself)."""
        key = "_hEw" if whole else "_hEc" if cell else "_hE"   # cell: a narrow anchor letter counts as its circle, as for marks
        if getattr(a, key, None) is None:
            one = (getattr(a, "_anchor1", False) or cell) and not whole
            hs = [(h, s) for x in a._g for h, s in zip(x["hulls"], x["sup"]) if not (one and h is getattr(a, "_anchor_hull", None))]
            setattr(a, key, (np.vstack([np.hstack([h, np.roll(h, -1, axis=0)]) for h, s in hs]),
                             np.cumsum([0] + [len(h) for h, s in hs[:-1]]).astype(int), np.array([s for h, s in hs], bool)) if hs else None)
        if getattr(a, key) is None: return 0.0
        E, st0, sup = getattr(a, key)
        mv = np.where(sup, min(m, self.st.lw / 2 + self.SUP_GAP * self.st.bond / 14.4), m)
        clear = lambda T: (_polys_dist(o + T[:, None] * u, E, st0) >= mv[None, :]).all(1) & (T >= disc)
        T = np.arange(0.0, 6 * self.st.size, 0.25); out = clear(T)
        if out[0]: return 0.0
        if not out.any(): return float(T[-1])
        k = int(np.argmax(out)); TT = np.linspace(T[k - 1], T[k], 9)
        return float(TT[int(np.argmax(clear(TT)))])

    # ---------- marks (charges, radicals, lone pairs)
    def _inner(self, a):
        """Angular sectors (lo, width) pointing into rings that contain atom a."""
        m = a.mol; out = []
        for b in m.bonds:
            if a not in (b.a, b.b): continue
            ring = m.ring_of(b)
            if not ring: continue
            nb = [x for x in m.nbrs(a) if x in ring]
            if len(nb) < 2: continue
            c = np.mean([x.xy for x in ring], axis=0); t1, t2 = ang(nb[0].xy - a.xy), ang(nb[1].xy - a.xy)
            w = (t2 - t1) % 360; lo = t1
            if abs(adiff(t1 + w / 2, ang(c - a.xy))) > 90: lo, w = t2, (t1 - t2) % 360
            out.append((lo, w))
        return out
    def _in_ring(self, a, t, pad=8):
        return any(((t - lo) % 360) > pad and ((t - lo) % 360) < w - pad for lo, w in self._inner(a))
    def _owner_ok(self, a, c, skip=(), f=1.8, d=None):
        """True if a mark centred at c (pt) is clearly nearer its own atom than any other atom (f times; 1.8 for lone
        pairs, 1.3 for charges, which sit further out), so its owner is never unclear (d: its distance to its own atom's
        ink, when an orbital's lobe drawn on the atom is nearer than the atom)."""
        d0 = f * (norm(c - self.P(a)) if d is None else d)
        return all(norm(c - self.P(x)) >= d0 for m in (self._ctx or [a.mol]) for x in m.atoms
                   if x is not a and not x.hidden and x not in skip)
    def _occ(self, a):
        occ = [ang(n.xy - a.xy) for n in a.mol.nbrs(a)]
        if a._g and len(a._g) > 1:
            for e in {"right": [0], "left": [180], "down": [-90], "up": [90], "pivot": [0, 180]}[a._dir]:
                occ += [e - 30, e, e + 30]                     # (a label read through its middle, EtOH: both sides)
        for ex in getattr(self, "extras", []):                 # orbital lobes occupy their directions too
            if ex[0] == "orbital" and ex[1] is a:
                occ += [ex[3] - 20, ex[3], ex[3] + 20] + ([ex[3] + 160, ex[3] + 180, ex[3] + 200] if ex[2] == "p" else
                                                          [ex[3] + 180] if ex[2] == "lobe" else [])
        return occ
    def _lobe_at(self, a, deg):
        """The empty hybrid lobe (S.orbital "lobe", electrons=0) drawn on atom a at angle deg, or None."""
        if deg is None: return None
        for e in getattr(self, "extras", []):
            if e[0] == "orbital" and e[1] is a and e[2] == "lobe" and not e[5] and abs(adiff(e[3], deg)) < 1: return e
        return None
    def _lp_geom(self, a, deg):
        """Lone pair at angle deg: (centre, dot 1, dot 2). Both dots lie on one circle about the atom, lp_gap clear of
        the label, so they are always the same distance from it; a pair in a hybrid lobe drawn there (S.orbital "lobe",
        no electrons of its own) sits inside the lobe, where the orbital's own electrons would."""
        st = self.st; A = self.P(a); u = dirv(deg)
        lb = self._lobe_at(a, deg)
        if lb is not None:
            c, v = self._lobe_pair(A, u, self._lobe_tip(a, u, A, lb[6]), lb[6], form=self._lobe_form(a, deg))
            return c, c - v, c + v
        R = (self._cell if a.drawn else 0.16 * st.bond) + st.lp_gap + st.lp_r
        if a.drawn and a._g:                  # (a wide lone label, Nu or Base: clear of every letter, not just the cell)
            R = max(R, self._exit_hull(a, u, A, st.lp_gap + st.lp_r, whole=self._lone_wide(a) is not None))
        dl = math.degrees(math.asin(min(0.9, st.lp_sep / 2 / R)))
        d1, d2 = A + dirv(deg - dl) * R, A + dirv(deg + dl) * R
        return (d1 + d2) / 2, d1, d2
    def _lone_wide(self, a):
        """For an atom with no bonds whose label runs on to the right of its anchor letter by more than 0.55 em (Base,
        Nu, Hal): the label's glyph outlines; else None."""
        if a.mol.nbrs(a) or not a.drawn or not a._g: return None
        busy = a.lps or any(isinstance(c.src, Atom) and c.src is a or isinstance(c.src, LP) and c.src.atom is a
                            or c.dst is a for c in self.curlies)
        if len(label_units(a.label)) != 1 and busy: return None   # a reacting ion (:OH with its arrow) keeps the charge on
        H = [h for g in a._g for h in g.get("hulls", [])]        # the O or N it belongs to; a spectator reads NMe4⊕, Base⊖
        if not H: return None
        right = max(float(h[:, 0].max()) for h in H)
        return H if right - self.P(a)[0] > 0.55 * self.st.size else None
    def _lone_mark(self, a, r, H):
        """Centre of a charge or radical after a lone label's last letter: out along 45 deg from that letter's middle until
        the mark keeps 1.2 pt of white from every letter."""
        A = self.P(a); right = max(float(h[:, 0].max()) for h in H)
        last = max(H, key=lambda h: float(h[:, 0].max()))
        L = np.array([(float(last[:, 0].min()) + right) / 2, A[1]])
        E = np.vstack([np.hstack([h, np.roll(h, -1, axis=0)]) for h in H])
        u = dirv(45.0)
        for t in np.arange(0.0, 3.0 * self.st.size, 0.1):
            c = L + u * t
            if float(_segs_dist(c[None, :], E).min()) >= r + 1.2: return c
        return L + u * 3.0 * self.st.size
    def _textmark(self, a, kind):
        """(text, size) of a text mark: a δ, or a carbene's + - (circled as the house charges, or plain)."""
        if kind == "pc": return "δ" + ("+" if a.partial == "+" else "−"), 0.85 * self.st.size
        return ("⊕⊖" if self.st.charge == "circle" else "+−"), (self.st.size if self.st.charge == "circle" else 0.8 * self.st.size)
    def _mark_gap(self, a, c, r, t, kind):
        """(white, most allowed) between a charge or radical at c (radius r, direction t) and its atom: white to the
        label's ink (from a bare vertex), allowed 3.2 pt (4.5) or 1 pt more than the mark's uncrowded place in that
        direction leaves (beside a narrow letter the character cell alone leaves more): placement and check() agree."""
        A = self.P(a)
        if a.drawn and a._g:
            H = np.vstack([h for x in a._g for h in x["hulls"]]); gap = lambda q: float(np.hypot(*(H - q).T).min()) - r; lim = 3.2
        else: gap = lambda q: norm(A - q) - r; lim = 4.5
        c0 = A + dirv(t) * self._mark_dist(a, t, r + {"rad": 0.08, "chg": 0.12}[kind] * self.st.size)
        return gap(c), max(lim, gap(c0) + 1.0)
    def _lobe_dist(self, a, c):
        """Distance (pt) from c to the outline of the orbital lobes drawn on atom a (inf if none)."""
        lob = [np.asarray(q[1], float) for e in getattr(self, "extras", []) if e[0] == "orbital" and e[1] is a
               for q in self._orbital_prims(e) if q[0] == "orb"]
        return float(np.hypot(*(np.vstack(lob) - c).T).min()) if lob else math.inf
    def _by_own_lobe(self, a, c, r):
        """True if a mark at c (radius r) sits beside a lobe of an orbital drawn on its own atom (2.6 pt of white at
        most): the lobe says whose it is (C⊕ by its empty p orbital), however far the lobe made it stand off."""
        return self._lobe_dist(a, c) - r - self.st.lw / 2 <= 2.6
    def _mark_dist(self, a, deg, extra):
        if not a.drawn: return extra + 0.16 * self.st.bond
        return self._exit(a, dirv(deg), margin=0.05 * self.st.size, cell=True) + extra
    def _ctx_geo(self, a):
        """What a mark on atom a must keep clear of (pt): nearby atoms, label boxes, bond centre-lines, curly arrows and
        orbital lobes."""
        mols = self._ctx or [a.mol]; A = self.P(a); R = 2.2 * self.st.bond
        pts = [self.P(x) for m in mols for x in m.atoms if x is not a and norm(self.P(x) - A) < R]
        boxes = [g["box"] for m in mols for x in m.atoms if x is not a for g in (x._g or [])
                 if box_dist(A, g["box"]) < R]
        boxes += [bk for x in a.mol.atoms if x is not a                                   # marks already placed on this
                  for k in ("_chg", "_rad", "_pc", "_num", "_site", "_amb")                         # molecule's other atoms (a site
                  for bk in [self._mark_box(x, k)] if bk and norm(getattr(x, k)[0] - A) < R]   # label, a charge) are in the way
        segs = []
        for m in mols:
            for b in m.bonds:
                if min(norm(self.P(b.a) - A), norm(self.P(b.b) - A)) >= R: continue
                segs.append((self.P(b.a), self.P(b.b), b))
                if b.style in ("plain", "bold", "dash", "partial") and (b.order == 2 or b.style == "partial"):
                    sd = self._dside(b) if b.order == 2 else (self._dside(b) or 1)      # the second line, set off to one side
                    n = perp(unit(self.P(b.b) - self.P(b.a)))
                    for off in ([-self.st.dbl / 2, self.st.dbl / 2] if sd == 0 else [sd * self.st.dbl]):
                        if off: segs.append((self.P(b.a) + n * off, self.P(b.b) + n * off, b))
                elif b.order == 3 and b.style in ("plain", "bold"):            # a triple bond's outer lines (N≡C, C≡C)
                    n = perp(unit(self.P(b.b) - self.P(b.a)))
                    segs += [(self.P(b.a) + n * off, self.P(b.b) + n * off, b) for off in (-self.st.dbl, self.st.dbl)]
        for p_ in getattr(self, "pointers", []):                 # the line a pointer will take, out from its target
            if any(m is a.mol for m in p_.mols()): segs += [(q0, q1, None) for q0, q1 in self._ptr_rays(p_)]
        lobes = [np.asarray(pr[1]) for ex in getattr(self, "extras", []) if ex[0] == "orbital" for pr in self._orbital_prims(ex) if pr[0] == "orb"]
        src = getattr(self, "_arrow_src", []); geo_a = getattr(self, "_arrow_geo", [])
        own = [g for g, sa in zip(geo_a, src) if sa is a]  # an arrow that starts at a's radical dot does not push the dot away
        arr = [g[np.hypot(*(g - A).T) < R + 4] for g in [g for g in geo_a if not any(g is o for o in own)] + lobes]   # arrows and orbital lobes
        return pts, boxes, segs, (np.vstack(arr) if arr else np.zeros((0, 2)))
    def _clear(self, a, c, r, geo):
        """Free distance around a mark circle (centre c, radius r): to every bond line, other labels and curly arrows."""
        pts, boxes, segs, arr = geo; d = 99.0
        for p, q, b in segs:
            ab = q - p; t = max(0.0, min(1.0, ((c - p) @ ab) / max(ab @ ab, 1e-9))); d = min(d, norm(p + t * ab - c) - r - self.st.lw / 2)
        for bx in boxes: d = min(d, box_dist(c, bx) - r)
        if len(arr): d = min(d, float(np.hypot(*(arr - c).T).min()) - r - self.st.alw / 2 - 0.6)
        return d
    def _clear_box(self, bx, geo):
        """Free distance around a text box (x0, y0, x1, y1): to every bond line, other labels and curly arrows."""
        pts, boxes, segs, arr = geo; d = 99.0
        for p, q, b in segs: d = min(d, seg_box_gap(p, q, bx) - self.st.lw / 2)
        for b2 in boxes: d = min(d, box_gap(bx, b2))
        if len(arr):
            c = np.clip(arr, [bx[0], bx[1]], [bx[2], bx[3]]); d = min(d, float(np.hypot(*(arr - c).T).min()) - self.st.alw / 2 - 0.6)
        return d
    def _site_block(self, a, taken, geo, home):
        """A site label (Nu, E, 2nd E, "E" over "more δ+"): placed as its true text box, not a circle round it, beside
        the dot as close as it fits (1.3 pt of white to bonds, labels and marks), in the roomiest direction near home;
        a label of several lines is aligned on the side facing the dot."""
        st = self.st; O = self._site_xy(a); R0 = (self._halo_r(a) if a.drawn else self.SITE_R * st.lw) + 1.6; multi = "\n" in a.site[1]; rel = {}
        best = None
        for p in range(0, 360, 15):
            if self._in_ring(a, p): continue
            u = dirv(p); cx = math.cos(math.radians(p))
            ha = ("left" if cx > 0.3 else "right" if cx < -0.3 else "center") if multi else "center"
            if ha not in rel: rel[ha] = self._runs_box(self._site_geo(a, (0.0, 0.0), ha)["runs"])
            B = rel[ha]; d = 0.0
            while d < 4 * st.bond and box_dist((0.0, 0.0), (B[0] + u[0] * d, B[1] + u[1] * d, B[2] + u[0] * d, B[3] + u[1] * d)) < R0:
                d += 0.25                                      # the box just clear of the dot, then out until clear of the rest
            for k in range(24):
                c = O + u * (d + 0.03 * st.bond * k); bx = (c[0] + B[0], c[1] + B[1], c[0] + B[2], c[1] + B[3])
                if self._clear_box(bx, geo) >= 1.3: break
            else: continue
            own = all(norm(c - self.P(x)) >= 1.05 * norm(c - O) for x in a.mol.atoms if x is not a and x.site)
            key = (own, k <= 6, min([abs(adiff(p, t)) for t in taken] + [90]) >= 45, -abs(adiff(p, home)), u[1] > 0.1, -k)
            # (close enough: within 2.6 pt of touching; then the usual side, and above before below, win over a point or
            # two nearer elsewhere)
            if best is None or key > best[0]: best = (key, c, p)
        if best is None: c, p = O + dirv(home) * (R0 + 0.5 * st.bond), home
        else: _, c, p = best
        a._site_ha = "left" if math.cos(math.radians(p)) > 0.3 else "right" if math.cos(math.radians(p)) < -0.3 else "center"
        bx = self._runs_box(self._site_geo(a, c)["runs"])
        return (c, 0.5 * math.hypot(bx[2] - bx[0], bx[3] - bx[1]), p)
    def _runs_box(self, runs):
        """Box (pt) of laid-out text runs as check() measures it: a quarter em under the baseline to cap height, each run
        (a superscript too) where it sits."""
        pts = [(x, y - 0.25 * sz, x + self.F.adv(s, sz, it, bd), y + self.F.cap * sz) for s, x, y, sz, it, bd in runs if s.strip()]
        return (min(p[0] for p in pts), min(p[1] for p in pts), max(p[2] for p in pts), max(p[3] for p in pts))
    def _site_geo(self, a, c=None, ha=None):
        """The site label's text as drawn (rich-text geometry) at its place c: centred, or aligned to the side facing its
        dot when it has several lines."""
        st = self.st; c = a._site[0] if c is None else c
        ha = (ha or getattr(a, "_site_ha", "center")) if "\n" in a.site[1] else "center"
        gb = self._rich_geo(a.site[1], 0.0, 0.0, st.tsize)["box"]
        x0 = c[0] + {"left": -(gb[2] - gb[0]) / 2, "right": (gb[2] - gb[0]) / 2}.get(ha, 0.0)
        return self._rich_geo(a.site[1], x0, c[1], st.tsize, ha)
    def _mark_box(self, a, k):
        """Box (pt) of mark k ('_chg', '_rad', '_pc', '_num', '_site') of atom a: a site label by its text."""
        mk = getattr(a, k, None)
        if not mk: return None
        if k == "_site": return self._runs_box(self._site_geo(a)["runs"])
        if k == "_amb": t_, sz_ = self._textmark(a, "amb"); return self._runs_box(self._rich_geo(t_, mk[0][0], mk[0][1], sz_)["runs"])
        return (mk[0][0] - mk[1], mk[0][1] - mk[1], mk[0][0] + mk[1], mk[0][1] + mk[1])
    def _lp_set(self, a, taken, auto, geo, pfix=()):
        """Angles for the automatic lone pairs `auto`, clear of the directions `taken`. Two or more pairs on one atom
        (pfix: the angles of its other, fixed pairs) stand at least 90 deg apart where the room allows, 120 for two
        (the textbook water, ether and carbonyl O), so the gap between pairs is twice the gap inside one and each
        pair reads as a pair, never one even arc of dots."""
        st = self.st; A = self.P(a); step = 5 if len(auto) < 3 else 10
        pulls = [l._pull for l in auto if l._pull is not None]
        npair = len(auto) + len(pfix); S = 120 if npair == 2 else 90
        def gp(cb):                            # the closest two pairs on the atom
            return min([abs(adiff(x, y)) for x, y in itertools.combinations(list(cb) + list(pfix), 2)] or [180])
        oc = sorted({round(o % 360, 3) for o in self._occ(a)}); secs = []   # the sectors between bonds (and lobes,
        if len(oc) >= 2 and npair >= 2:                                     # a label's other letters)
            secs = [(oc[i], (oc[(i + 1) % len(oc)] - oc[i]) % 360 or 360) for i in range(len(oc))]
            if max(w for _, w in secs) < 180: secs = []
        def cramp(cb):                         # pairs inside a bond angle under 180 deg while a wider open side is free
            return sum(1 for x in cb for lo, w in secs if (x - lo) % 360 < w and w < 180)
        bd = [ang(n.xy - a.xy) for n in a.mol.nbrs(a)]   # (a label's own letters only need the 50 deg kept off them)
        def lp_clear(t):
            c = self._lp_geom(a, t)[0]
            return self._clear(a, c, st.lp_sep / 2 + st.lp_r, geo)
        memo = getattr(self, "_lpmemo", None); mk = (id(a), step, len(auto), id(geo), tuple(pfix))
        if memo is not None and mk in memo: free, cl, combos = memo[mk]
        else:
            free = [t for t in range(0, 360, step) if not self._in_ring(a, t) and self._owner_ok(a, self._lp_geom(a, t)[0])] \
                or [t for t in range(0, 360, step) if not self._in_ring(a, t)] or list(range(0, 360, step))
            cl = {t: lp_clear(t) for t in free}
            combos = []
            if npair >= 2:                     # several pairs: only the sets 88 deg apart, if there are any (a few hundred
                def rec(i0, cur):              # instead of every set of directions)
                    if len(cur) == len(auto): combos.append(tuple(cur)); return
                    for i in range(i0, len(free)):
                        if all(abs(adiff(free[i], y)) >= 88 for y in cur + list(pfix)): rec(i + 1, cur + [free[i]])
                rec(0, [])
            combos = combos or list(itertools.combinations(free, len(auto)))
            if memo is not None: memo[mk] = (free, cl, combos)
        def mn(combo):
            ds = [abs(adiff(x, y)) for x in combo for y in taken] + [abs(adiff(x, y)) for x, y in itertools.combinations(combo, 2)]
            return round(min(ds)) if ds else 180
        mns = {cb: mn(cb) for cb in combos}
        thr = min(50, max(mns.values()))
        lone = not a.mol.nbrs(a) and not pulls # a lone atom or ion (B:, Cl-) keeps spare pairs on the axes: ":B⊖", not a
        def score(cb):                         # pair slanted down-left opposite the charge
            m = mns[cb]; pull = sum(min(abs(adiff(x, p)) for x in cb) for p in pulls)
            crowd = sum(min(0.0, cl[x] - 1.0) for x in cb)
            ax = sum(min(abs(adiff(x, k)) for k in (0, 90, 180, 270)) < 1 for x in cb) if lone else 0
            if npair >= 2:                     # several pairs: apart, out of any narrow bond angle, spread, then each
                g = gp(cb)                     # pair's nearest bonds as far off as can be; an arrow's pair is the one
                lex = tuple(sorted(round(abs(adiff(x, y)) / 5) for x in cb for y in bd)[:2 * len(cb)])   # nearest its
                return (m >= thr, g >= 88, -cramp(cb), round(crowd), min(g, S) // 15, lex, -round(pull / 10), ax, m,   # pull
                        -sum(abs(adiff(x, 90)) for x in cb))
            if lone:
                return (m >= thr, round(crowd), -round(pull / 10), ax, m, -sum(abs(adiff(x, 90)) for x in cb))
            return (m >= thr, round(crowd), -round(pull / 10), m, -sum(abs(adiff(x, 90)) for x in cb))
        best = max(combos, key=score); return best, mns[best]
    def _marks(self, a):
        self._lpmemo = {}                     # free directions and their clearances do not change within one call
        try: return self._marks_(a)
        finally: self._lpmemo = None
    def _marks_(self, a):
        st = self.st; A = self.P(a); occ = self._occ(a) + list(a._avoid)
        isauto = lambda l: l.auto and not getattr(l, "_locked", False)
        taken = occ + [l.angle for l in a.lps if not isauto(l)]
        auto = [l for l in a.lps if isauto(l)]; geo = self._ctx_geo(a)
        for t in [e[3] for e in getattr(self, "extras", []) if e[0] == "orbital" and e[1] is a and e[2] == "lobe" and not e[5]]:
            if auto: auto.pop(0).angle = float(t)   # the pair of an empty hybrid lobe drawn on the atom (a carbene's sp2)
        a._chg = a._rad = None
        def mindist(t, ts): return min([abs(adiff(t, o)) for o in ts] or [180])
        a._pc = a._num = a._site = a._amb = None
        orb_e = any(e[0] == "orbital" and e[1] is a and e[5] for e in getattr(self, "extras", []))   # an orbital drawn
        orb_l = any(e[0] == "orbital" and e[1] is a and e[2] != "s" for e in getattr(self, "extras", []))   # lobes on it
        lp_home = None                       # a bare carbanion's one lone pair (no arrow from it): straight out along the
        if not a.drawn and len(auto) == 1 and auto[0]._pull is None and a.charge and a.charge_at is None and occ:
            bs = sorted(o % 360 for o in occ)               # middle of the widest gap outside any ring, the charge beside
            gaps = [((bs[(i + 1) % len(bs)] - bs[i]) % 360 or 360, bs[i]) for i in range(len(bs))]   # it (of equal gaps,
            mids = [(round(g), (b0 + g / 2) % 360) for g, b0 in gaps if not self._in_ring(a, (b0 + g / 2) % 360)] or \
                   [(round(g), (b0 + g / 2) % 360) for g, b0 in gaps]   # the one furthest from the ring, then lowest)
            inn = [lo + w / 2 for lo, w in self._inner(a)]
            lp_home = float(max(mids, key=lambda gm: (gm[0], min([abs(adiff(gm[1], x)) for x in inn] or [0]) // 5,
                                                        -round(math.sin(math.radians(gm[1])), 2)))[1]); taken = taken + [lp_home]
            _, d1, d2 = self._lp_geom(a, lp_home); e = 0.5    # the charge keeps 1 pt clear of the pair's dots (their cores
            geo = (geo[0], geo[1] + [(d[0] - e, d[1] - e, d[0] + e, d[1] + e) for d in (d1, d2)], geo[2], geo[3])   # here)
        def assign(angs):                    # an arrow's pair takes the place nearest its pull
            angs = list(angs)
            for l in sorted(auto, key=lambda l: l._pull is None):
                t = min(angs, key=lambda x: abs(adiff(x, l._pull)) if l._pull is not None else 0)
                angs.remove(t); l.angle = float(t)
        ptaken = []                          # two or more pairs go first, as the atom's electron domains, so they stay
        if len(a.lps) >= 2:                  # clearly apart; a charge then sits in a gap between two of them, just clear
            pfix = [l.angle for l in a.lps if not isauto(l)]   # of their dots (⊖ beside :Ö–H, as textbooks draw it)
            if auto:
                pins = [t for t, f in ((a.charge_at, a.charge), (a.radical_at, a.radical)) if f and t is not None]
                assign(self._lp_set(a, occ + pfix + pins, auto, geo, pfix)[0]); auto = []
            ptaken = [l.angle for l in a.lps]; taken = list(occ)
        pdots = [d for l in a.lps if l.angle is not None for d in self._lp_geom(a, l.angle)[1:]] if ptaken else \
                (list(self._lp_geom(a, lp_home)[1:]) if lp_home is not None else [])
        def md(p): return min(mindist(p, taken), mindist(p, ptaken) + 20)   # a mark may come nearer a pair than a bond
        ownb = []                            # boxes of this atom's charge and radical, once placed
        for kind, flag, at, prefs in (("chg", a.charge, a.charge_at, [45, 135, -45, -135, 90, -90, 0, 180]),   # with its
                                      ("rad", a.radical and not orb_e, a.radical_at, [90, 45, 135, -90, 0, 180, -45, -135]),  # electron
                                      ("pc", a.partial, None if a.charge else a.charge_at, [90, -90, 45, 135, -45, -135, 0, 180]),   # δ
                                      ("amb", getattr(a, "ambi", False), None, [-90, 90, -45, -135, 45, 135, 0, 180]),   # + -
                                      ("num", a.number, None, [-45, -135, 45, 135, -90, 90, 0, 180]),
                                      ("site", a.site and a.site[1], None, [45, 135, 30, 150, -45, -135, -30, -150, 90, -90, 0, 180])):
            if not flag: continue
            if kind == "num" and a.drawn and a._g:         # a locant on a labelled atom (the 1 of a Me) keeps clear of it
                geo = (geo[0], geo[1] + [g_["box"] for g_ in a._g], geo[2], geo[3])
            if kind == "pc": tw = self._runs_w(parse_rich("δ+"), 0.85 * st.size); r = 0.5 * max(tw, 0.85 * st.size)   # its ink,
                                                                     # not the corners of its box: δ sits as close as a charge
            elif kind == "num": tw = self._runs_w(parse_rich(str(a.number)), 0.7 * st.size); r = 0.5 * math.hypot(tw, 0.7 * st.size)
            elif kind == "amb":
                t_, sz_ = self._textmark(a, "amb"); bb_ = self._runs_box(self._rich_geo(t_, 0.0, 0.0, sz_)["runs"])
                r = 0.5 * max(bb_[2] - bb_[0], bb_[3] - bb_[1])
            elif kind == "site":                            # a site label (Nu, E) beside its dot: clear of the dot as well
                tw = self._runs_w(parse_rich(_typo(a.site[1], chem=True)), st.tsize); r = 0.5 * math.hypot(tw, self.F.cap * st.tsize) + 0.4
            else: r = {"chg": st.chg_r if st.charge == "circle" else 0.3 * st.size, "rad": 1.2 * st.lp_r}[kind]
            ex = r + {"rad": 0.08, "chg": 0.12, "pc": 0.12, "amb": 0.12, "site": 0.1}.get(kind, 0.04) * st.size
            if kind == "pc" and a.drawn and a._g:           # δ on a label: its ink sits 1.5 pt from the label's ink, as a
                dg = self._rich_geo("δ+", 0.0, 0.0, 0.85 * st.size)["box"]          # superscript would, not a circle apart
                dh = ((dg[2] - dg[0]) / 2, 0.5 * self.F.cap * 0.85 * st.size + 0.3)
                HP = np.vstack([h for x in a._g for h in x["hulls"]]) - A
            if kind == "site":                              # the dot names the site's atom; its label need only be nearer
                def own_ok(c):                              # its own dot than any other dot
                    return all(norm(c - self.P(x)) >= 1.05 * norm(c - A) for x in a.mol.atoms if x is not a and x.site)
            else:                                           # (beside its own orbital's lobe, measured from the lobe:
                def own_ok(c):                              # the lobe shows whose it is)
                    d = self._lobe_dist(a, c) if orb_l and kind in ("chg", "rad") and self._by_own_lobe(a, c, r) else None
                    return self._owner_ok(a, c, f=1.3, d=None if d is None else min(d, norm(c - A)))
            pcb = self._runs_box(self._rich_geo(self._textmark(a, kind)[0], 0.0, 0.0, self._textmark(a, kind)[1])["runs"]) \
                if kind in ("pc", "amb") else None    # a δ (a carbene's + -) is kept clear by its text box, as check() measures it
            def room(c, r=r, pcb=pcb):
                if pcb is None: return self._clear(a, c, r, geo)
                return self._clear_box((c[0] + pcb[0], c[1] + pcb[1], c[0] + pcb[2], c[1] + pcb[3]), geo) + 0.3
            def pos(p, ex=ex, r=r, fine=False):
                """Nearest radius (stepping outward; fine: in 0.14 pt steps) at which the mark clears other bonds, labels
                and atoms."""
                if kind == "pc" and a.drawn and a._g:
                    u = dirv(p); base = float((HP @ u).max()) + 1.5 + abs(u[0]) * dh[0] + abs(u[1]) * dh[1]
                elif kind == "amb" and not a.drawn:        # a carbene's + -: its box just clear of the vertex
                    u_ = dirv(p); base = 0.0
                    while base < 2 * st.bond and box_dist((0.0, 0.0), (pcb[0] + u_[0] * base, pcb[1] + u_[1] * base,
                                                                       pcb[2] + u_[0] * base, pcb[3] + u_[1] * base)) < 0.16 * st.bond:
                        base += 0.25
                else: base = self._mark_dist(a, p, ex)
                O = self._site_xy(a) if kind == "site" else A       # a site label stands off its dot, wherever the dot is
                u = dirv(p)
                for d in pdots:                                     # between the atom's own pairs: 1 pt of white to their
                    w = d - O; t = float(u @ w); q = (r + st.lp_r + 1.0) ** 2 - (float(w @ w) - t * t)   # dots, no more
                    if q > 0: base = max(base, t + math.sqrt(q))
                sbx = [h for x in (a._g or []) for h, sp in zip(x["hulls"], x["sup"]) if sp] if kind != "site" else []
                while sbx and base < 4 * st.size and min(_hull_dist(O + u * base, h) for h in sbx) < r + 0.6:
                    base += 0.25                                    # clear of the label's own superscript (the II of Pd^II)
                for k in (np.arange(0.0, 0.37, 0.01) if pdots or fine else np.arange(0.0, 0.73, 0.03) if orb_l else (0.0, 0.12, 0.24, 0.36)):
                    c = O + u * (base + k * st.bond)
                    if room(c) >= 1.3: return c
                return O + u * base
            def far(c, t, r=r):                                # beyond what check() allows a charge or radical
                g, lim = self._mark_gap(a, c, r, t, kind); return g > lim and not (orb_l and self._by_own_lobe(a, c, r))
            if kind == "site":                                 # placed as its text box, from the middle of the widest gap
                hm = 45                                        # (of equal gaps, the one nearest the upper right)
                if occ:
                    bs = sorted(o % 360 for o in occ); gaps = [((bs[(i + 1) % len(bs)] - bs[i]) % 360 or 360, bs[i]) for i in range(len(bs))]
                    gm = max(g for g, _ in gaps)
                    hm = min((round((b0 + g / 2) % 360, 1) for g, b0 in gaps if g >= gm - 1), key=lambda h: abs(adiff(h, 45)))
                a._site = self._site_block(a, taken, (geo[0], geo[1] + ownb, geo[2], geo[3]), hm); taken = taken + [a._site[2]]; continue
            lw_ = self._lone_wide(a) if kind == "chg" and at is None else None   # (a radical stays by its own atom)
            if lw_ is not None:                                # a lone ion with a wide label (Base⊖, Nu⊖): the charge goes
                c = self._lone_mark(a, r, lw_)                 # after the last letter, as a superscript charge would, not over
                taken = taken + [45]                           # the middle of the word
                if kind == "chg": a._chg = (c, r, 45)
                else: a._rad = (c, r, 45)
                continue
            if at is None:
                home = prefs[0]
                if kind == "chg" and lp_home is not None:      # beside the pair, on the side nearer the usual upper right,
                    cn = [ang(n.xy - a.xy) for n in a.mol.nbrs(a) if n.charge]   # not towards a charged neighbour (an
                    home = min(((lp_home + 60) % 360, (lp_home - 60) % 360),      # ylide's C⊖ by its N⊕ reads N⊕⊖)
                               key=lambda h: (min([abs(adiff(h, x)) for x in cn] + [180]) < 90, abs(adiff(h, 45))))
                if kind == "site" and occ:                      # a site label (Nu, E) takes the middle of the widest gap left
                    bs = sorted(o % 360 for o in occ)           # by the bonds and any pointer coming in; of equal gaps, the one
                    gaps = [((bs[(i + 1) % len(bs)] - bs[i]) % 360 or 360, bs[i]) for i in range(len(bs))]   # nearest the
                    gm = max(g for g, _ in gaps)                # upper right
                    home = min((round((b0 + g / 2) % 360, 1) for g, b0 in gaps if g >= gm - 1), key=lambda h: abs(adiff(h, 45)))
                if kind == "amb" and occ:                   # a carbene's + -: in the widest gap, below first
                    bs = sorted(o % 360 for o in occ); gaps = [((bs[(i + 1) % len(bs)] - bs[i]) % 360 or 360, bs[i]) for i in range(len(bs))]
                    gm = max(g for g, _ in gaps)
                    home = min((round((b0 + g / 2) % 360, 1) for g, b0 in gaps if g >= gm - 1), key=lambda h: abs(adiff(h, 270)))
                if (kind == "num" or (kind == "pc" and not a.drawn)) and occ:   # a locant, or a δ on a bare vertex, sits
                    bs = sorted(o % 360 for o in occ)       # out along the middle of the widest gap between the atom's bonds:
                    gaps = [((bs[(i + 1) % len(bs)] - bs[i]) % 360 or 360, bs[i]) for i in range(len(bs))]   # radially
                    gm = max(g for g, _ in gaps)            # outside a ring, on the outer side of a chain bend; of equal
                    mids = [round((b0 + g / 2) % 360, 1) for g, b0 in gaps if g >= gm - 1]   # gaps, above, then below
                    home = min(mids, key=lambda h: (abs(adiff(h, 90)) > 1, abs(adiff(h, 270)) > 1, abs(adiff(h, 90))))
                def roomy(p, ang_min, clr):
                    return mindist(p, taken) >= ang_min - 0.5 and mindist(p, ptaken) >= 40 and room(pos(p)) >= clr and \
                        (not auto or lp_home is not None or self._lp_set(a, taken + [p], auto, geo)[1] >= 60)
                pull = getattr(a, "_radpull", None) if kind == "rad" else None
                if pull is not None:
                    ok = [p for p in range(0, 360, 15) if not self._in_ring(a, p) and roomy(p, 50, 1.5)]
                    if ok: at = min(ok, key=lambda p: abs(adiff(p, pull)))
                if at is not None: pass
                elif roomy(home, 60, 1.3 if kind == "pc" or ptaken else 2.5) and not self._in_ring(a, home) and own_ok(pos(home)): at = home
                else:
                    cands = [p for p in range(0, 360, 15) if not self._in_ring(a, p)] or list(range(0, 360, 15))
                    # room counts up to 2 pt of white; past that the usual place (upper right) wins over a roomier one
                    def good(p):                        # first, places where a charge stays as near as check() wants
                        c = pos(p); cl = room(c); own = own_ok(c)
                        near = kind not in ("chg", "rad") or not far(c, p) or not far(pos(p, fine=True), p)
                        return (cl >= 1.2 and own and near, cl >= 1.2 and own, cl >= 1.2, own, cl >= 0.6)
                    if kind == "num" and not any(good(p)[0] for p in cands):   # a locant with no room outside its rings
                        cands = list(range(0, 360, 15))                        # sits in a ring's corner beside its atom
                    at = max(cands, key=lambda p: (*good(p), min(md(p), 90) // 15, min(room(pos(p)), 2.0) // 0.5,
                                                   not auto or lp_home is not None or self._lp_set(a, taken + [p], auto, geo)[1] >= 60, -abs(adiff(p, home))))
            taken = taken + [at]
            c_ = pos(at)
            if kind in ("chg", "rad") and far(c_, at): c_ = pos(at, fine=True)   # the coarse steps overshot: as close as fits
            if kind in ("chg", "rad"): ownb.append((c_[0] - r, c_[1] - r, c_[0] + r, c_[1] + r))   # a site label keeps
                                                                                                # clear of these
            if kind == "chg": a._chg = (c_, r, at)
            elif kind == "rad": a._rad = (c_, r, at)
            elif kind == "pc": a._pc = (pos(at), r, at)
            elif kind == "amb": a._amb = (pos(at), r, at)
            elif kind == "site": a._site = (pos(at), r, at)
            else: a._num = (pos(at), r, at)
        if auto and lp_home is not None: auto[0].angle = lp_home
        elif auto: assign(self._lp_set(a, taken, auto, geo)[0])
        a._lp = [(l, *self._lp_geom(a, l.angle)) for l in a.lps]

    # ---------- bonds
    @staticmethod
    def _dsubs(m, x, other):
        """Substituents that count for double-bond sidedness: drawn neighbours other than the partner.
        An H drawn on a labelled atom (C=O+-H, C=N-H) counts as part of that label."""
        out = []
        for y in m.nbrs(x):
            if y is other or getattr(y, "hidden", False): continue
            if x.drawn and y.label == "H" and len(m.nbrs(y)) == 1: continue
            out.append(y)
        return out
    def _dside_rule(self, b, force_offset=False):
        """IUPAC GR-1.10 (ChemDraw's convention): +1 second line to the left of a->b, -1 right, 0 centred."""
        m = b.a.mol; u = b.b.xy - b.a.xy; ring = m.ring_of(b)
        if ring:
            c = np.mean([x.xy for x in ring], axis=0); return 1 if cross(u, c - b.a.xy) > 0 else -1
        na, nb = self._dsubs(m, b.a, b.b), self._dsubs(m, b.b, b.a)
        ul = unit(u)                                                   # a group in line with the bond (the far end of
        side = lambda xs, y: [x for x in xs if abs(cross(ul, unit(x.xy - y.xy))) > 0.17]   # an allene) is on neither side
        sa, sb = side(na, b.a), side(nb, b.b)
        ea, eb = (sa if b.a.drawn else na), (sb if b.b.drawn else nb)  # (past a label it does not meet the pair: C=C=O)
        if not force_offset:
            if b.a.drawn and b.b.drawn: return 0                       # between two labels: N=N, S=O
            if not ea and not eb: return 0                             # O=O, isolated ethene
            wide = any(bb.style == "wedge" and bb.b in (b.a, b.b) and not bb.b.drawn for bb in m.bonds if bb is not b)
            if ((not ea and len(nb) >= 2) or (not eb and len(na) >= 2)) and not wide:
                return 0                                               # GR-1.10.2 carbonyl, =CH2 (offset beside the wide
                                                                       # end of a wedge, which a centred pair cannot meet)
        if not sa and not sb: return 1
        s = sum(np.sign(cross(u, x.xy - b.a.xy)) for x in sa + sb)
        if s != 0: return int(np.sign(s))                              # GR-1.10.1 side with more substituents
        inring = lambda x: len(m.nbrs(x)) > 1 and any(m.ring_of(bb) for bb in m.bonds if x in (bb.a, bb.b) and bb is not b)
        heavy = max(sa + sb, key=lambda x: (inring(x), len(m.nbrs(x)), x.drawn, -x.uid))   # tie: ring, then bigger group
        return 1 if cross(u, heavy.xy - b.a.xy) > 0 else -1
    _trim = False                 # single bonds stop at a centred C=O's near line only in the drawn lines (_build_ sets
                                  # it after the arrows are placed, so their solutions never depend on it)
    def _centred_trim(self, X, b, w):
        """How far short of bare vertex X (pt, along w, which points from X into bond b) a single bond stops when X also
        carries a centred double bond (an acid's or ketone's C=O): at that bond's near line, where the line meets it, so
        nothing runs on into the gap between its two lines and the two meet in a clean corner. 0 if there is none."""
        if X.drawn or X.hidden or X.clear is not None: return 0.0
        st = self.st; m = X.mol; P = self.P(X)
        for b2 in m.bonds:
            if b2 is b or X not in (b2.a, b2.b) or b2.order != 2 or b2.style not in ("plain", "bold"): continue
            Z = b2.b if b2.a is X else b2.a                                       # a C=O, C=S or C=N to a label
            if self._dside(b2) != 0 or not (Z.drawn and not Z.hidden): continue
            w2 = unit(self.P(Z) - P); n2 = perp(w2)
            if abs(cross(w2, w)) < math.sin(math.radians(40)): continue            # (nearly in line: leave it)
            o = n2 * (st.dbl / 2) * (1 if float(w @ n2) > 0 else -1)
            M = np.array([[w[0], -w2[0]], [w[1], -w2[1]]])
            try: t, s_ = np.linalg.solve(M, o)
            except np.linalg.LinAlgError: continue
            if 0 < t < 0.4 * st.bond and abs(s_) < 0.4 * st.bond: return float(t)
        return 0.0
    def _dside(self, b):
        if b.side is not None: return {"left": 1, "right": -1, "center": 0}.get(b.side, b.side)
        us = getattr(b, "_uside", None)
        if us is not None: return us
        return self._dside_rule(b)
    def _unify_sides(self):
        """One depiction per double bond across a scheme: if the same bond (same atom uids, via Mol.copy) in the same
        surroundings (in a ring or not, the same number of substituents at each end) is offset in any step, it is offset
        in every step, never jumping between the two. A bond whose surroundings change keeps its own natural form (an
        exocyclic enolate C=CH2 is centred even if the same carbons later form a ring C=C, which is offset)."""
        groups = {}
        for m in self.mols:
            for b in m.bonds:
                b._uside = None
                if b.order == 2 and b.style in ("plain", "bold") and b.side is None:
                    env = (bool(m.ring_of(b)), tuple(sorted((len(self._dsubs(m, b.a, b.b)), len(self._dsubs(m, b.b, b.a))))))
                    groups.setdefault((frozenset((b.a.uid, b.b.uid)), env), []).append(b)
        for key, bs in groups.items():
            if len(bs) < 2: continue
            nat = [self._dside_rule(b) for b in bs]
            if any(v != 0 for v in nat) and any(v == 0 for v in nat):
                for b, v in zip(bs, nat):
                    if v == 0: b._uside = self._dside_rule(b, force_offset=True)
    def _bond(self, b):
        st = self.st; A, B = self.P(b.a), self.P(b.b); u = unit(B - A); n = perp(u)
        col = b.color or st.ink; prims = []; m = b.a.mol
        def junction(X, other, w, ov, mode):
            """Distance along w (into the bond) where a line offset by ov should end at unlabelled atom X."""
            if mode == "axis" and b.order == 1 and b.style in ("plain", "bold", "dash") and self._trim:
                return self._centred_trim(X, b, w)       # a centred C=O at that vertex: stop at its near line
            if mode == "axis" or norm(ov) < 1e-9: return 0.0
            sg = np.sign(cross(w, ov)); best = None
            for Y in m.nbrs(X):
                if Y is other: continue
                nb = m.find(X, Y)                            # every drawn bond counts, dashed and wavy ones too: the
                v = unit(Y.xy - X.xy)                        # inner line keeps clear of whatever line is drawn there
                if np.sign(cross(w, v)) != sg: continue
                th = math.acos(max(-1.0, min(1.0, float(w @ v))))
                if best is None or th < best[0]: best = (th, v, nb)
            if best is None: return 0.0                      # unsubstituted end or trans substituent: full length (IUPAC GR-1.10)
            th, v, nb = best
            if mode == "inner":
                t0 = norm(ov) / math.tan(th / 2)
                if nb is not None and nb.style == "wedge" and nb.b is X:     # beside a wedge's wide end: keep the same
                    wq = (st.bold if getattr(nb, "wide", None) == "bold" else st.wedge) / 2   # clearance as from a line
                    c0, c1 = cross(v, ov), cross(v, w)
                    if abs(c1) > 1e-6 and np.sign(c0) == np.sign(c1):
                        t0 = max(t0, (wq + norm(ov) - st.lw / 2 - abs(c0)) / abs(c1))
                return t0
            M = np.array([[w[0], -v[0]], [w[1], -v[1]]])
            if abs(np.linalg.det(M)) < 1e-6: return 0.0
            t, s = np.linalg.solve(M, -ov)
            if th > math.radians(150):                     # nearly in line with the neighbour (a C=O off a cage's sharp
                e = perp(v) * np.sign(cross(v, ov)) * st.lw / 2       # angle): stop where the line enters the
                try: t = float(np.linalg.solve(M, e - ov)[0])         # neighbour's edge, not its centre line,
                except np.linalg.LinAlgError: pass                    # or a sliver of it runs on outside
            return float(max(-2 * norm(ov), min(2 * norm(ov), t)))
        def seg(off, sh=(0.0, 0.0), mode="axis"):
            oa, ob = A + n * off, B + n * off
            # at a label every line of a multiple bond ends level with the centre line, as in ChemDraw
            ta = max(self._exit(b.a, u, oa), self._exit(b.a, u)) if (b.a.drawn or b.a.clear is not None) else junction(b.a, b.b, u, n * off, mode)
            tb = max(self._exit(b.b, -u, ob), self._exit(b.b, -u)) if (b.b.drawn or b.b.clear is not None) else junction(b.b, b.a, -u, n * off, mode)
            return oa + u * max(ta, sh[0]) if mode != "center" else oa + u * ta, \
                   ob - u * max(tb, sh[1]) if mode != "center" else ob - u * tb
        lw = st.bold if b.style == "bold" else st.lw
        dash = {"dash": (3.2 * st.lw, 2.4 * st.lw), "dot": (0.1, 2.2 * st.lw)}.get(b.style)
        offs = [0.0]
        if b.style in ("plain", "bold", "dash", "dot", "partial"):
            if b.order == 2 or b.style == "partial":
                s = self._dside(b) if b.order == 2 else (self._dside(b) or 1)
                offs = [-st.dbl / 2, st.dbl / 2] if s == 0 else [0.0, s * st.dbl]
                inner = s != 0
            elif b.order == 3: offs, inner = [-st.dbl, 0.0, st.dbl], False
            else: inner = False
            for i, off in enumerate(offs):
                mode = "inner" if (inner and i == 1) else ("center" if (b.order == 2 and not inner and off != 0) else "axis")
                p, q = seg(off, (0.0, 0.0), mode)
                if b.style == "dot" and norm(q - p) > 2 * lw:   # dots clear of the letters, centred, and as many as on
                    p, q, k_, off_ = self._dot_layout(b); dash = (dash[0], off_)   # H-bonds alike beside them
                d = (3.2 * st.lw, 2.4 * st.lw) if (b.style == "partial" and i == 1) else dash
                if d is not None: d = self._fit_dashes(norm(q - p), d, least=2)   # (a lone dash reads as a plain stub)
                prims.append(("line", [p, q], lw, col, d, "round" if b.style == "dot" else "butt"))
        elif b.style == "wedge":
            p, q = seg(0.0); Lf = norm(B - A)
            wmax = st.bold if getattr(b, "wide", None) == "bold" else st.wedge      # into a bold bar: as wide as the bar
            s0 = max(0.0, float((p - A) @ u))                 # where the wedge comes out of its narrow end's label
            def wid(X):                                       # cut by a label, it starts narrow there and widens at a
                if s0 > 0.5: t = max(0.0, min(1.0, (float((X - A) @ u) - s0) / st.bond))   # full wedge's rate: slim
                else: t = max(0.0, min(1.0, float((X - A) @ u) / Lf))                     # between two labels
                return (st.wedge_tip + (wmax - st.wedge_tip) * t) / 2
            wp, wq = wid(p), wid(q)
            chain = {1: [q + n * wq], -1: [q - n * wq]}; edge = {}; convex = set()
            X = b.b
            if not X.drawn and not X.hidden and X.clear is None:
                # the wide end at a bare vertex joins the bonds there: on a side where the next bond round is less than
                # 180 deg away, the wedge's edge runs on to that bond's near edge; on the open side of a two-bond vertex
                # it meets the other bond's outer edge in a mitre (a bevel if that would be long), as two lines would
                nbs = []                                       # (direction, half-width at X, at the far end, length)
                for Y in m.nbrs(X):
                    if Y is b.a: continue
                    nb = m.find(X, Y)
                    if nb is None: continue
                    if nb.style == "wedge" and nb.a is X:      # a wedge going on from its narrow end (two wedges in a
                        h = st.wedge_tip / 2                   # row): its edges, widening, take the joins, so the two
                        h2 = (st.bold if getattr(nb, "wide", None) == "bold" else st.wedge) / 2   # fill as one shape
                    elif nb.style in ("plain", "bold", "partial", "dative"):
                        h = (st.bold if nb.style == "bold" else st.lw) / 2
                        if nb.order == 2 and nb.style in ("plain", "bold") and self._dside(nb) == 0: h += st.dbl / 2   # centred pair
                        h2 = h
                    else: continue
                    nbs.append((unit(self.P(Y) - B), h, h2, norm(self.P(Y) - B)))
                wdir = -u; ang_w = lambda t: math.acos(max(-1.0, min(1.0, float(t[0] @ wdir))))
                for sg in (1, -1):
                    side = [t for t in nbs if float(t[0] @ (n * sg)) > 1e-6]
                    far = side == []
                    if far:                                    # nothing this side: the bond met going the long way round
                        if not nbs: continue
                        v, h, h2, Ln = max(nbs, key=ang_w)
                    else:
                        v, h, h2, Ln = min(side, key=ang_w)
                        if ang_w((v,)) > math.radians(170): continue
                    face = perp(v); face = face if float(face @ wdir) > 0 else -face   # the bond's edge towards the wedge
                    eg = -face if far else face                # on the open side: its other edge, facing the gap
                    e0 = B + eg * (h if far else h - 0.3)      # (0.3 pt inside on the closed side: no hairline seam)
                    ve = unit(v * Ln + eg * (h2 - h))          # that edge's own direction (a wedge's widens)
                    M = np.array([[u[0], -ve[0]], [u[1], -ve[1]]])
                    if abs(np.linalg.det(M)) < 1e-6: continue
                    s_, t_ = np.linalg.solve(M, e0 - chain[sg][0]); C = chain[sg][0] + u * s_
                    if far:                                    # mitre (bevel if long), then along that edge into the bond
                        pts = ([C] if norm(C - B) <= 3 * wq + h else [chain[sg][0]]) + [e0 + ve * h, B + v * h]
                        chain[sg] = pts if sg == -1 else pts[::-1]; convex.add(sg)
                    elif t_ >= -h and norm(C - B) <= 3 * wq + h: chain[sg] = [C]; edge[sg] = (e0, ve)
            J = [q]                                        # between the two sides: where the bonds' near edges meet
            if convex: J = []
            elif len(edge) == 2:
                (P1, v1), (P2, v2) = edge[1], edge[-1]; M = np.array([[v1[0], -v2[0]], [v1[1], -v2[1]]])
                if abs(np.linalg.det(M)) > 1e-6:
                    a_, _ = np.linalg.solve(M, P2 - P1); Jc = P1 + v1 * a_
                    if norm(Jc - B) <= 3 * wq: J = [Jc]
            elif edge: J = [B]
            prims.append(("poly", [p + n * wp, p - n * wp] + chain[-1] + J + chain[1], col))
        elif b.style == "hash":
            p, q = seg(0.0); Lf = norm(B - A)
            d0 = max(0.0, float((p - A) @ u)); d1 = float((q - A) @ u)
            if not b.b.drawn and self._trim: d1 -= self._centred_trim(b.b, b, -u)   # (a centred C=O there: from its near line)
            if not b.b.drawn and len(m.nbrs(b.b)) > 1: d1 -= 0.9 * st.hash_step   # clear the far vertex's other bonds
            else: d1 -= 0.5 * st.lw
            span = d1 - d0; k = max(int(round(span / st.hash_step)), 2)
            for i in range(1, k + 1):
                d = d0 + span * i / k; c = A + u * d
                t = min(1.0, (d - d0) / st.bond) if d0 > 0.5 else min(1.0, d / Lf)    # as for wedges: slim from a label
                w = (st.hash_tip + (st.hash_wide - st.hash_tip) * t) / 2
                prims.append(("line", [c - n * w, c + n * w], st.lw, col, None, "butt"))
        elif b.style == "wavy":                                # waves from end to end: it comes off each atom already
            p, q = seg(0.0); L = norm(q - p); half, amp, nmin = self.WAVE              # waving, no straight lead, and ends on
            nh = max(int(round(L / (half * st.bond))), nmin)                            # the bond's line (whole half-waves);
            t = np.linspace(0, 1, 40 * nh)                                              # round ends meet the vertex cleanly
            pts = p + (q - p) * t[:, None] + n[None, :] * (amp * st.bond * np.sin(math.pi * nh * t))[:, None]
            prims.append(("line", list(pts), st.lw, col, None, "round"))
        elif b.style == "dative":
            p, q = seg(0.0); h, w = 0.65 * st.head_len, 0.65 * st.head_w
            prims.append(("line", [p, q - u * h * 0.6], lw, col, None, "butt"))
            prims.append(("poly", [q, q - u * h + n * w, q - u * h * 0.72, q - u * h - n * w], col))
        e = b.b if getattr(b.b, "_eta", None) else b.a if getattr(b.a, "_eta", None) else None
        if e is not None and prims and prims[0][0] == "line" and not getattr(e, "_persp3d", False):   # eta: stop at the first line
            P0, P1 = np.array(prims[0][1][0]), np.array(prims[0][1][-1])
            if e is b.a: P0, P1 = P1, P0
            best = None
            for bb in m.bonds:
                if bb is b or not (bb.a in e._eta and bb.b in e._eta): continue
                for pr in getattr(bb, "_prims", []):
                    if pr[0] != "line": continue
                    q0, q1 = np.array(pr[1][0]), np.array(pr[1][-1]); den = cross(P1 - P0, q1 - q0)
                    if abs(den) < 1e-9: continue
                    t = cross(q0 - P0, q1 - q0) / den; s = cross(q0 - P0, P1 - P0) / den
                    if 0 < t <= 1.0001 and -0.05 <= s <= 1.05 and (best is None or t < best): best = t
            if best is not None:
                Q = P0 + (P1 - P0) * best; seg2 = [P0, Q] if e is b.b else [Q, P0]
                prims[0] = (prims[0][0], seg2) + tuple(prims[0][2:])
        b._offs = offs if b.style not in ("wedge", "hash") else [-st.wedge / 4, st.wedge / 4]
        self._limited = False; p, q = seg(0.0); b._limited = self._limited
        b._mid = (p + q) / 2; b._vis = (p, q); b._u = u
        b._gaps = self._behind(b)
        if b._gaps and b.style in ("plain", "bold", "dash"):
            prims = [pr for x in prims for pr in self._cut(x, b._gaps)]
        b._prims = prims
        return prims
    MARGIN = 2.0                  # white (pt) kept either side of a bond where another passes behind it (ChemDraw's margin)
    STUB = 3.0                    # what a broken bond keeps showing either side of its gap (pt), next to an atom or a label
    def _gap_fit(self, b, b2, t, stub=None):
        """The gap in b where b2 crosses in front of it at fraction t: half-widths (towards b.a, towards b.b) in pt, or None
        if b cannot take one. The gap keeps MARGIN of white beside b2 (0.8 pt at the least, where room is short) and leaves
        `stub` pt of b showing on each side of it, between the gap and the label or atom at that end."""
        if b.style not in ("plain", "bold", "dash"): return None
        st = self.st; stub = self.STUB if stub is None else stub
        A, B = self.P(b.a), self.P(b.b); C, D = self.P(b2.a), self.P(b2.b); L = norm(B - A); u = (B - A) / max(L, 1e-9)
        w = {2: st.dbl, 3: st.dbl}.get(b2.order, 0.0) + (st.wedge / 2 if b2.style in ("wedge", "hash") else st.lw / 2)
        sin = max(abs(cross(u, unit(D - C))), 0.35)
        s0 = self._exit(b.a, u) if (b.a.drawn or b.a.clear is not None) else 0.0           # where b shows: from one
        s1 = L - (self._exit(b.b, -u) if (b.b.drawn or b.b.clear is not None) else 0.0)    # label's edge to the other's
        c = t * L; gf, gm = (w + self.MARGIN) / sin, (w + 0.8) / sin
        ga, gb = min(gf, c - s0 - stub), min(gf, s1 - c - stub)
        return (ga, gb) if min(ga, gb) >= gm else None
    def _gap_owner(self, b, b2, t, s_):
        """Which of two crossing bonds is broken, and its gap: (bond, (towards a, towards b)) or None (drawn whole). The one
        marked behind=True, else the deeper one (atom depths z), else, with no depth known, a dashed partial bond under a
        solid one; that bond is broken if it keeps STUB pt showing either side, else the other one is (as an illustrator
        does near an atom), else the crossing stays whole. A bond marked behind is broken in any case, its gap narrowed
        to keep a stub where it can."""
        ex, ex2 = getattr(b, "behind", False), getattr(b2, "behind", False)
        if ex or ex2:
            x, y, tx = (b, b2, t) if ex else (b2, b, s_)
            if x.style not in ("plain", "bold", "dash"): return None
            f = self._gap_fit(x, y, tx) or self._gap_fit(x, y, tx, stub=1.0)
            if f is None:                                            # no room for a stub: the full gap, as asked
                A, B = self.P(x.a), self.P(x.b); C, D = self.P(y.a), self.P(y.b); st = self.st
                w = {2: st.dbl, 3: st.dbl}.get(y.order, 0.0) + (st.wedge / 2 if y.style in ("wedge", "hash") else st.lw / 2)
                g = (w + self.MARGIN) / max(abs(cross(unit(B - A), unit(D - C))), 0.35); f = (g, g)
            return x, f
        zs = [getattr(x, "z", None) for x in (b.a, b.b, b2.a, b2.b)]
        if all(z is not None for z in zs):
            behind_b = zs[0] + t * (zs[1] - zs[0]) < zs[2] + s_ * (zs[3] - zs[2])
        elif (b.style == "dash") != (b2.style == "dash"): behind_b = b.style == "dash"
        else: return None                                            # (check() reports it)
        for x, y, tx in (((b, b2, t), (b2, b, s_)) if behind_b else ((b2, b, s_), (b, b2, t))):
            f = self._gap_fit(x, y, tx)
            if f is not None: return x, f
        return None
    def _behind(self, b):
        """Where bond b passes behind another bond of its molecule: [(crossing point (pt), half-width of the gap along b,
        the other bond, half-widths towards b.a and towards b.b, b's direction)] (_gap_owner decides which is broken)."""
        m = b.a.mol; A, B = self.P(b.a), self.P(b.b); out = []
        for b2 in m.bonds:
            if b2 is b or {b.a, b.b} & {b2.a, b2.b} or b2.order == 0: continue
            C, D = self.P(b2.a), self.P(b2.b); den = cross(B - A, D - C)
            if abs(den) < 1e-9: continue
            t = cross(C - A, D - C) / den; s_ = cross(C - A, B - A) / den
            if not (0.02 < t < 0.98 and 0.02 < s_ < 0.98): continue
            own = self._gap_owner(b, b2, t, s_)
            if own is None or own[0] is not b: continue
            ga, gb = own[1]; out.append((A + (B - A) * t, max(ga, gb), b2, ga, gb, unit(B - A)))
        return out
    def _cut(self, pr, gaps):
        """A line prim with the stretches under the given gaps taken out."""
        if pr[0] != "line" or len(pr[1]) != 2: return [pr]
        P, Q = np.array(pr[1][0], float), np.array(pr[1][1], float); L = norm(Q - P)
        if L < 1e-9: return [pr]
        u = (Q - P) / L; cuts = []
        for g in gaps:                                     # (a gap's two half-widths run along its bond, a to b)
            x = float((g[0] - P) @ u); ga, gb = (g[3], g[4]) if len(g) > 5 else (g[1], g[1])
            cuts.append((x - ga, x + gb) if len(g) <= 5 or float(u @ g[5]) >= 0 else (x - gb, x + ga))
        def piece(p0, p1):                                 # a dashed line's pieces each hold whole dashes again
            d = pr[4] if len(pr) > 4 else None
            if d is not None and d[0] > 0.2: d = self._fit_dashes(norm(p1 - p0), d)
            return (pr[0], [p0, p1]) + tuple(pr[2:4]) + ((d,) + tuple(pr[5:]) if len(pr) > 4 else ())
        cuts.sort(); out, t0 = [], 0.0
        for c0, c1 in cuts:
            if c0 > t0 + 0.3: out.append(piece(P + u * t0, P + u * min(c0, L)))
            t0 = max(t0, c1)
        if t0 < L - 0.3: out.append(piece(P + u * t0, Q))
        return out

    def _dot_span(self, b):
        """Where dotted bond b's dots may run (pt): from where a solid bond would end at each label, lw/2 further in for
        the dots' round ends, and a bond end's margin from the label's letters themselves too (a line leaving an H by its
        corner passes nearer the letter than the H's circle allows for)."""
        st = self.st; A, B = self.P(b.a), self.P(b.b); u = unit(B - A)
        def cut(X, w):                                   # at a bare vertex where other bonds meet, the dots keep their
            if not (X.drawn or X.clear is not None):     # rhythm from it, the bonds' ends standing in for a dot
                return 2.2 * st.lw - st.lw / 2 if any(x is not b and X in (x.a, x.b) for x in X.mol.bonds) else 0.0
            t = self._exit(X, w)
            if X.clear is None and X._g: t = max(t, self._exit_hull(X, w, self.P(X), st.margin, whole=True))
            return t + st.lw / 2
        return A + u * cut(b.a, u), B - u * cut(b.b, -u)
    def _dot_layout(self, b):
        """(p, q, n, gap) of dotted bond b: n dots centred between p and q, their round ends `gap` apart. A dotted bond
        holds the number of dots its length holds at the usual spacing, rounded, the spacing eased (to no less than 0.7
        of the usual) when rounded up: a short H-bond keeps three dots rather than two, which could read as a lone pair.
        Dotted bonds meeting at one atom and about as long (within 40%) hold the same number: the two H-bonds of
        N-H...Cl...H-N look alike even where one ends at the l of Cl and the other at its C."""
        on, off = 0.1, 2.2 * self.st.lw
        p, q = self._dot_span(b); L = norm(q - p); u = unit(q - p) if L > 0 else np.array([1.0, 0.0])
        Ls = [L]
        for X in (b.a, b.b):
            for b2 in X.mol.bonds:
                if b2 is b or b2.style != "dot" or X not in (b2.a, b2.b): continue
                p2, q2 = self._dot_span(b2); L2 = norm(q2 - p2)
                if max(L, L2) <= 1.4 * min(L, L2): Ls.append(L2)
        most = lambda v: int((v + 0.7 * off) // (on + 0.7 * off))          # dots that fit, eased to the limit
        each = [min(int((v + off) / (on + off) + 0.5), most(v)) for v in Ls]
        n = max(2, min(max(each), min(most(v) for v in Ls)))
        g = min(off, (L - n * on) / (n - 1)); e = (L - n * on - (n - 1) * g) / 2
        return p + u * e, q - u * e, n, g
    @staticmethod
    def _fit_dashes(L, d, least=1):
        """Dash pattern (on, off) stretched or shrunk a little so a line of length L holds whole dashes and begins and ends
        with one: dashed bonds look the same at both ends, and dashes meet cleanly where two dashed lines join. Dots (on
        tiny) are spaced evenly from end to end the same way, so the white at both ends is the same."""
        on, off = d
        if L <= 0: return d
        n = max(least, int(round((L + off) / (on + off)))); k = L / (n * on + (n - 1) * off)
        return (on * k, off * k)
    def _joints(self, m):
        """Fill the open outer corner wherever bond lines end at one point beside an unlabelled atom: at the atom itself
        (bonds meeting) or inside the angle, where the inner lines of two double or partial bonds meet (mitre join,
        bevel if sharp). Solid and dashed lines and dative bonds alike: a dashed line starts and ends with a whole dash, so
        its end joins like a solid one (dots and wavy lines have round ends and need nothing). Only a gap wider than 180 deg has a fill; elsewhere
        the line ends already overlap, and a fill would poke out."""
        st = self.st; out = []
        for a in m.atoms:
            if a.drawn or a.hidden: continue
            V = self.P(a); ends = []
            for b in m.bonds:
                if a not in (b.a, b.b) or b.style not in ("plain", "bold", "partial", "dash", "dative"): continue
                for pr in b._prims:
                    if pr[0] != "line" or (pr[4] and pr[4][0] <= 0.2): continue
                    P = np.asarray(pr[1], float)
                    if len(P) < 2: continue
                    e, f = (P[0], P[1]) if norm(P[0] - V) <= norm(P[-1] - V) else (P[-1], P[-2])
                    if norm(e - V) < 0.5 * st.bond and norm(f - e) > 1e-6: ends.append((e, unit(f - e), pr[2] / 2, pr[3]))
            for b in m.bonds:                              # a wedge's narrow end is a line end too (as thin as a bond)
                if b.style == "wedge" and b.a is a:
                    ends.append((V, unit(self.P(b.b) - V), st.wedge_tip / 2, b.color or st.ink))
            groups = []                                    # line ends that meet at one point
            for t in ends:
                for g in groups:
                    if norm(g[0][0] - t[0]) < 0.05: g.append(t); break
                else: groups.append([t])
            for g in groups:
                if len(g) < 2: continue
                G = np.mean([t[0] for t in g], axis=0); lines = sorted([t[1:] for t in g], key=lambda t: ang(t[0]) % 360)
                for i in range(len(lines)):
                    (u2, h2, c2), (u1, h1, c1) = lines[i], lines[(i + 1) % len(lines)]
                    if not 181 < (ang(u1) - ang(u2)) % 360 < 359: continue  # only a gap wider than 180 deg leaves a corner
                    n1 = -perp(u1); n2 = perp(u2)                   # the two edges facing that gap
                    A1, A2 = G + n1 * h1, G + n2 * h2
                    M = np.array([[u1[0], -u2[0]], [u1[1], -u2[1]]])
                    try: s, t = np.linalg.solve(M, A2 - A1)
                    except np.linalg.LinAlgError: continue
                    X = A1 + u1 * s
                    d1, d2 = u1 * h1, u2 * h2                      # reach into each line body: no seams, no stroke
                    gw = (ang(u1) - ang(u2)) % 360                 # another bond of the atom inside the corner (a C=O off a
                    bis = ang(u2) + gw / 2                         # a centred double bond (a C=O off a cage's sharp angle)
                    busy = gw > 270 and any(abs(adiff(ang(self.P(bb.b if bb.a is a else bb.a) - V), bis)) < 35   # on the
                                            for bb in m.bonds if a in (bb.a, bb.b) and bb.order == 2 and self._dside(bb) == 0)
                                                                             # corner: its two lines straddle where the
                                                                             # mitre would spike, so the corner is bevelled
                    if norm(X - G) < 3 * max(h1, h2) and not busy: pts = [G + d1, A1 + d1, X, A2 + d2, G + d2]
                    else: pts = [G + d1, A1 + d1, A2 + d2, G + d2]; X = (A1 + d1 + A2 + d2) / 2
                    if c1 == c2: out.append(("poly", pts, c1))
                    else:                                          # two colours meet: each line fills its own half
                        out.append(("poly", [G + d1, A1 + d1, X, G], c1)); out.append(("poly", [G, X, A2 + d2, G + d2], c2))
        return out
    # ---------- preparation
    def _ptbu(self, x):
        if isinstance(x, LP): return x.atom.xy
        if isinstance(x, Atom): return x.xy
        if isinstance(x, Mid): return x.a.xy * (1 - x.t) + x.b.xy * x.t
        if isinstance(x, (tuple, list)) and len(x) == 2 and isinstance(x[0], Atom): return (x[0].xy + x[1].xy) / 2
        return np.array(x, float)
    def _lig_space(self, m):
        """Ligands from ligate(): slide each one-label ligand along its bond until the bond shows LIG_SHOW between the
        two labels (the white a 1.3 BU bond leaves between two round one-letter labels). Moves the atom itself, so
        layout, arrows and checks all see the spaced position; the same result every time."""
        st = self.st; show0 = 1.3 * st.bond - 2 * (self._cell + st.margin)
        movable = {}
        for b in m.bonds:
            M = getattr(b, "_lig", None)
            if M is None or M not in (b.a, b.b): continue
            X = b.b if M is b.a else b.a
            if not (M.drawn and X.drawn) or X.hidden or len(m.nbrs(X)) != 1 or not M._g or not X._g: continue
            movable.setdefault(id(M), (M, []))[1].append(X)
        for M, Xs in movable.values():
            fixed = []                                     # the metal's other label-to-label bonds (a chelate ring's) set the
            for Y in m.nbrs(M):                            # length, so every bond round one metal looks the same
                if Y in Xs or not Y.drawn or Y.hidden or not Y._g or m.find(M, Y).free: continue   # (not one drawn at a
                                                                                              # length of its own)
                u = unit(self.P(Y) - self.P(M)); fixed.append(norm(self.P(Y) - self.P(M)) - self._exit(M, u) - self._exit(Y, -u))
            show = max(float(np.median(fixed)), 0.3 * st.bond) if fixed else show0
            for X in Xs:
                u = unit(self.P(X) - self.P(M)); D = norm(self.P(X) - self.P(M))
                want = show + self._exit(M, u) + self._exit(X, -u)
                if abs(want - D) > 0.02:
                    X.xy = M.xy + u * want / st.bond
                    X._g = self._label(X); X._hE = X._hEc = X._hEw = None
    LBL_SHOW = 6.0          # pt of line a bond between two labels shows at least (ChemDraw's H-Cl at the standard length)
    LOBE_L = 0.78           # an orbital lobe's length (BU at size 1): a soft teardrop, LOBE_Q as wide as it is long,
    LOBE_Q = 0.6            # its sides leaving the nucleus LOBE_TIP deg off its axis and bowing out into a round far end
    LOBE_TIP = 30.0         # (met LOBE_CAP deg off the axis), as textbooks draw a p orbital; on a labelled atom it comes
    LOBE_CAP = 12.0         # out from behind the letter, already LOBE_OPEN line widths across there, its outline open
    LOBE_OPEN = 2.5
    LOBE_GAP = 0.35         # pt: a lobe from behind a letter shows from just off its ink, so it reads as passing behind it
    LBL_SHOW_STEREO = 8.0   # pt a wedge or hash between two labels shows at least (three hash lines, a readable wedge)
    def _lbl_space(self, m):
        """A bond between two labels at the standard length shows a stub that reads as a hyphen (O-H, H-Cl, Pd-H, Pd-Ar²,
        C-R¹): the label at its far end, if that has no other bond, slides outward along the bond until the bond shows
        LBL_SHOW, as ligate() spaces its ligands (the bond becomes free, longer by design; a bond already drawn longer
        is lengthened only if a superscript still leaves a stub: Ni^III-Ar). A ring bond keeps its place."""
        st = self.st; need = self.LBL_SHOW * st.bond / 14.4
        arrowed = {y.uid for c in self.curlies for x in (c.src, c.dst) for y in self._xatoms(x)}   # a reacting atom stays
        others = [x for mm in (self._ctx or [m]) for x in mm.atoms if not x.hidden]               # where its helper put it
        moves = {}                                         # centre -> [(label, bond, unit vector, extra length)]
        for b in m.bonds:
            if getattr(b, "_lig", None) is not None or b.style in ("dative", "partial") or getattr(b.a, "_tmpl", None): continue
            if not (b.a.drawn and b.b.drawn) or b.a.hidden or b.b.hidden or not (b.a._g and b.b._g): continue
            for X, Y in ((b.b, b.a), (b.a, b.b)):          # the end with no other bond moves (of H-Cl alone, the second)
                if len(m.nbrs(X)) != 1 or X.clear is not None or X.uid in arrowed: continue
                u = unit(self.P(X) - self.P(Y)); D = norm(self.P(X) - self.P(Y))
                vis = D - self._exit(Y, u) - self._exit(X, -u)
                nd = self.LBL_SHOW_STEREO * st.bond / 14.4 if b.style in ("wedge", "hash") else need
                if vis < nd - 0.02: moves.setdefault(id(Y), (Y, []))[1].append((X, b, u, (nd - vis) / st.bond))
                break
        for Y, xs in moves.values():                       # one centre's labels move together (both O of SO2, all H of
            near = lambda X, p: min([norm(p - z.xy) for z in others if z is not X and z is not Y] or [9.0])   # BH4), as far
            f = min(next((f for f in (1.0, 0.75, 0.5, 0.25, 0.0)                                              # as the most
                          if near(X, X.xy + u * f * e) >= min(near(X, X.xy), 0.9) - 1e-9), 0.0) for X, b, u, e in xs)   # hemmed
            if f <= 0: continue                                                                               # in can go
            for X, b, u, e in xs:
                X.xy = X.xy + u * f * e; X._g = self._label(X); X._hE = X._hEc = X._hEw = None; b.free = True
    def _dative_need(self):
        return 1.6 * 0.65 * self.st.head_len                # the head and a shaft 0.6 of its length
    def _dative_space(self, m):
        """A dative bond between two labels: its arrowhead alone is longer than the white between the labels, so the
        head would sit against the donor with no shaft. A partner with no other bond (BF3, a ligand) slides out along
        the bond until head and shaft both show; the bond is then free (longer by design). Otherwise check() says so."""
        st = self.st; need = self._dative_need()
        for b in m.bonds:
            if b.style != "dative" or not (b.a.drawn and b.b.drawn) or b.a.hidden or b.b.hidden or not (b.a._g and b.b._g): continue
            for X, Y in ((b.b, b.a), (b.a, b.b)):          # move the acceptor if it stands alone, else the donor
                if len(m.nbrs(X)) != 1: continue
                u = unit(self.P(X) - self.P(Y)); D = norm(self.P(X) - self.P(Y))
                vis = D - self._exit(Y, u) - self._exit(X, -u)
                if vis < need - 0.02:
                    X.xy = Y.xy + u * (D + need - vis) / st.bond; X._g = self._label(X); X._hE = X._hEc = X._hEw = None; b.free = True
                break
    def _prep_mol(self, m):
        _linear_centres(m)
        for a in m.atoms:
            if getattr(a, "_newman_r", None) is not None: a.clear = a._newman_r * self.st.bond
        if self._ctx is None or m not in self._ctx:
            for a in m.atoms: a._g = self._label(a) if (a.drawn and not a.hidden) else None; a._hE = a._hEc = a._hEw = None
        self._lig_space(m); self._lbl_space(m); self._dative_space(m)
        for a in m.atoms: a._avoid = []
        for a in m.atoms: a._radpull = None
        for c in self.curlies:
            cd = self._edst(c)
            if isinstance(c.src, Atom) and c.src.mol is m and c.src.radical:
                c.src._radpull = ang(self._ptbu(cd) - c.src.xy)
            if isinstance(c.src, LP) and c.src.atom.mol is m and c.src.auto:
                a = c.src.atom; t = ang(self._ptbu(cd) - a.xy); c.src._pull = t
                if norm(self._ptbu(cd) - a.xy) > 1.3 and self._aclass(c) == "inter":
                    # an arrow to another molecule leaves a pair set 30 deg off the line to its target, so it curves like a
                    # hand-drawn arrow instead of running straight out of the pair: on the side with more room
                    occ = [ang(n.xy - a.xy) for n in m.nbrs(a)]
                    if len(label_units(a.label)) > 1:                       # the side the label's other letters take
                        occ += [0.0, 180.0] if getattr(a, "_pivot", 0) else [180.0 if self._ldir(a) == "left" else 0.0]
                                                                    # (the H of an OH, the Et of EtO-; both for EtOH)
                    occ.append(a.charge_at if a.charge_at is not None else 45.0)   # the charge's usual spot, charge or not, so
                                                                                    # a neutral base (Et3N) curls like HO- and EtO-
                    def room(d): return min([abs(adiff(t + d, o)) for o in occ] or [180.0])
                    c.src._pull = t + max((30.0, -30.0), key=lambda d: (round(room(d) / 10), round(dirv(t + d)[1], 3), round(dirv(t + d)[0], 3)))
            if isinstance(cd, Atom) and cd.mol is m:
                cd._avoid.append(ang(self._ptbu(c.src) - cd.xy))
            if isinstance(c.dst, Mid):
                if c.dst.a.mol is m: c.dst.a._avoid.append(ang(c.dst.b.xy - c.dst.a.xy))
                if c.dst.b.mol is m: c.dst.b._avoid.append(ang(c.dst.a.xy - c.dst.b.xy))
        for pt in getattr(self, "pointers", []):              # marks keep out of the way a pointer comes in
            sd = {"below": 270.0, "above": 90.0, "left": 180.0, "right": 0.0}.get(pt.side, pt.side)
            for t in pt.targets:
                if isinstance(t, Atom) and t.mol is m and sd is not None: t._avoid.append(float(sd))
        for a in m.atoms: a._chg = a._rad = a._pc = a._num = a._site = a._amb = None   # placed afresh, in atom order
        for a in m.atoms: self._marks(a)
        for b in m.bonds: self._bond(b)
        self._float_gaps(m)
        m._xprims = self._ring_extras(m)
    @staticmethod
    def _cut_poly(P, X, g):
        """Polyline P with the stretch within arc length g of its point nearest X removed: a list of 0-2 polylines."""
        P = np.asarray(P, float); seg = np.diff(P, axis=0); L = np.hypot(*seg.T); cum = np.concatenate([[0.0], np.cumsum(L)])
        best = (1e9, 0.0)
        for i in range(len(seg)):
            if L[i] < 1e-9: continue
            t = max(0.0, min(1.0, float((X - P[i]) @ seg[i]) / L[i] ** 2)); d = norm(P[i] + seg[i] * t - X)
            if d < best[0]: best = (d, cum[i] + t * L[i])
        s0, s1 = best[1] - g, best[1] + g
        def at(sv):
            i = int(np.clip(np.searchsorted(cum, sv) - 1, 0, len(seg) - 1)); t = (sv - cum[i]) / max(L[i], 1e-9)
            return P[i] + seg[i] * t
        out = []
        if s0 > 0.3: out.append(np.vstack([P[cum < s0], at(s0)]))
        if s1 < cum[-1] - 0.3: out.append(np.vstack([at(s1), P[cum > s1]]))
        return out
    @staticmethod
    def _clip_half(P, nv, d):
        """Convex polygon P clipped to the half-plane x . nv <= d."""
        out = []; k = len(P)
        for i in range(k):
            a, b = np.asarray(P[i], float), np.asarray(P[(i + 1) % k], float); da, db = float(a @ nv) - d, float(b @ nv) - d
            if da <= 0: out.append(a)
            if (da <= 0) != (db <= 0): out.append(a + (b - a) * (da / (da - db)))
        return out
    def _persp_outline(self, m, r):
        """Outline of a perspective ring as one stroke of varying width: each edge becomes a polygon (bold bar, side edges
        tapering from bold to plain, plain back edges) and neighbouring edges meet in a mitred corner (bevelled if
        sharper than a 3:1 mitre), so no butt end or notch shows at any corner. Replaces the edges' line and wedge
        prims; each polygon carries its spine and end widths. Labelled ring atoms keep plain ends. Returns extra
        prims (bevel triangles)."""
        st = self.st; A = r["atoms"]; n = len(A); E = []
        for i in range(n):
            x, y = A[i], A[(i + 1) % n]; b = m.find(x, y)
            if b is None or b.order != 1 or b.style not in ("plain", "bold", "wedge") or getattr(b, "_vis", None) is None: E.append(None); continue
            if b.style == "bold": w0 = w1 = st.bold
            elif b.style == "wedge":
                wide = st.bold if getattr(b, "wide", None) == "bold" else st.wedge
                w0, w1 = (st.wedge_tip, wide) if b.a is x else (wide, st.wedge_tip)
            else: w0 = w1 = st.lw
            p, q = (b._vis[0], b._vis[1]) if b.a is x else (b._vis[1], b._vis[0])     # the visible stretch, x to y
            E.append({"b": b, "p": np.asarray(p, float), "q": np.asarray(q, float), "w": (w0, w1)})
        def sides(e):                                    # (narrowed by the hairline _draw adds round each polygon)
            u = unit(e["q"] - e["p"]); nl = perp(u); h0, h1 = (e["w"][0] - self.SEAM) / 2, (e["w"][1] - self.SEAM) / 2
            return ((e["p"] + nl * h0, e["q"] + nl * h1), (e["p"] - nl * h0, e["q"] - nl * h1))
        def meet(l1, l2):
            (a0, a1), (b0, b1) = l1, l2; d1, d2 = a1 - a0, b1 - b0; den = cross(d1, d2)
            if abs(den) < 1e-9: return None
            return a0 + d1 * (cross(b0 - a0, d2) / den)
        ends = [None] * n; extra = []                    # ends[i] = [[L, R] at the start, [L, R] at the end]
        for i, e in enumerate(E):
            if e is not None: (L0, L1), (R0, R1) = sides(e); ends[i] = [[L0, R0], [L1, R1]]
        for i in range(n):                               # vertex i: the end of edge i-1 meets the start of edge i
            e0, e1 = E[i - 1], E[i]; V = A[i]
            if e0 is None or e1 is None or V.drawn or V.hidden: continue
            s0, s1 = sides(e0), sides(e1); PV = self.P(V); hmax = max(e0["w"][1], e1["w"][0]) / 2
            for k in (0, 1):                             # the left sides meet, and the right sides
                X = meet(s0[k], s1[k])
                if X is None: continue
                if norm(X - PV) > 3.0 * hmax:            # too sharp a corner: bevel it
                    extra.append(("poly", [PV, ends[i - 1][1][k], ends[i][0][k]], e1["b"].color or st.ink)); continue
                ends[i - 1][1][k] = X; ends[i][0][k] = X
        for i, e in enumerate(E):
            if e is None: continue
            (L0, R0), (L1, R1) = ends[i]
            e["b"]._prims = [("poly", [L0, L1, R1, R0], e["b"].color or st.ink, {"spine": (e["p"], e["q"]), "w": e["w"]})]
        return extra
    def _bond_lines(self, b):
        """A bond's drawn lines for the checks, as (start, end, extra): each straight line prim, and the spine of each
        outline polygon of a perspective ring (extra: how much it reaches beyond a plain line on each side)."""
        out = []
        for pr in b._prims:
            if pr[0] == "line": out.append((np.asarray(pr[1][0], float), np.asarray(pr[1][-1], float), 0.0))
            elif pr[0] == "poly" and len(pr) > 3 and isinstance(pr[3], dict):
                p, q = pr[3]["spine"]; out.append((np.asarray(p, float), np.asarray(q, float), max(0.0, (max(pr[3]["w"]) - self.st.lw) / 2)))
            elif pr[0] == "poly" and b.style == "wedge" and len(pr[1]) >= 4:   # a wedge: its two long edges (its wide end
                P = [np.asarray(x, float) for x in pr[1]]                       # can reach a label its spine stays clear of)
                out += [(P[1], P[2], -self.st.lw / 2), (P[0], P[-1], -self.st.lw / 2)]
        return out
    def _float_gaps(self, m):
        """Where a floating substituent's bond (Mol.floating) crosses its ring edge, the inner line of a double bond
        there is broken, as ChemDraw draws it; the edge's own line runs on under the crossing."""
        st = self.st
        for f in m.atoms:
            e = getattr(f, "_float", None); bd = m.find(*e) if e else None
            if bd is None or not m.nbrs(f): continue
            S0, S1 = self.P(f), self.P(m.nbrs(f)[0]); A = self.P(bd.a); nrm = perp(unit(self.P(bd.b) - A)); out = []
            for pr in bd._prims:
                if pr[0] == "line" and len(pr[1]) == 2:
                    p0, p1 = np.asarray(pr[1][0], float), np.asarray(pr[1][1], float)
                    if abs(float(((p0 + p1) / 2 - A) @ nrm)) > 0.3 * st.dbl and _seg_x(p0, p1, S0, S1):
                        r, sv = p1 - p0, S1 - S0; X = p0 + r * (cross(S0 - p0, sv) / cross(r, sv))
                        d = unit(r); g = st.lw / 2 + 1.2
                        out += [(pr[0], [p0, X - d * g], *pr[2:]), (pr[0], [X + d * g, p1], *pr[2:])]; continue
                out.append(pr)
            bd._prims = out
    def _ring_extras(self, m):
        """Prims for perspective rings and delocalised systems of m: each perspective ring's outline as one mitred stroke
        (_persp_outline), the circle or ellipse of each delocalised ring and the charge at its centre (or, for a chain,
        just clear of its partial bonds), and depth for a metal's lines to a perspective ring: the line goes behind the
        ring's near edges and in front of its far ones, broken where the nearer element passes."""
        st = self.st; out = []; m._xtexts = []
        persp = getattr(m, "_persp", [])
        for r in persp: out += self._persp_outline(m, r)  # 1. the outline
        rings = {}                                         # 2. delocalised systems
        for d in getattr(m, "_deloc", []):
            P = np.array([self.P(a) for a in d["atoms"]]); c = P.mean(0); mc = d["atoms"][0].color or st.ink
            if d["ring"]:
                r = next((r for r in persp if set(map(id, r["atoms"])) == set(map(id, d["atoms"]))), None)
                if r is not None:
                    n = len(r["atoms"]); rho = 0.64 * r["R"] * st.bond; k = 0.7 * r["k"]; tl = r.get("tilt", 0.0)   # flatter than
                                                                           # the ring, so its near arc keeps clear of the bold bar
                else:
                    segs = np.array([[*P[i], *P[(i + 1) % len(P)]] for i in range(len(P))])
                    rho = 0.6 * float(_segs_dist(c[None, :], segs).min()); k = 1.0; tl = 0.0
                t = np.linspace(0, 2 * math.pi, 121)
                E = c + np.array([rotv(v, tl) for v in np.column_stack([rho * np.cos(t), k * rho * np.sin(t)])])
                pr = ("line", list(E), st.lw, mc, None, "round"); out.append(pr); rings[id(d["atoms"][0])] = (pr, c)
                q = c
            else:
                mid = d["atoms"][len(d["atoms"]) // 2]; v = unit(c - self.P(mid)) if norm(c - self.P(mid)) > 1e-6 else np.array([0.0, -1.0])
                lines = [np.asarray(pr[1], float) for x, y in zip(d["atoms"], d["atoms"][1:]) for pr in m.find(x, y)._prims if pr[0] == "line"]
                segs = np.array([[*L_[i], *L_[i + 1]] for L_ in lines for i in range(len(L_) - 1)])
                q = c.copy()
                for _ in range(60):
                    if float(_segs_dist(q[None, :], segs).min()) >= st.chg_r + st.lw / 2 + 1.5: break
                    q = q + v * 0.5
            if d["charge"]:
                sgn = d["charge"] > 0
                if st.charge == "circle" and abs(d["charge"]) == 1:
                    sz = 0.58 * st.chg_r + self.MARK_W * st.lw / 2    # flat-ended strokes, like bonds, as long as before
                    out += [("circle", q, st.chg_r, mc, False, self.RING_W * st.lw), ("line", [q - (sz, 0), q + (sz, 0)], self.MARK_W * st.lw, mc, None, "butt")]
                    if sgn: out.append(("line", [q - (0, sz), q + (0, sz)], self.MARK_W * st.lw, mc, None, "butt"))
                else:
                    txt = (str(abs(d["charge"])) if abs(d["charge"]) > 1 else "") + ("+" if sgn else "−")
                    g = self._rich_geo(txt, q[0], q[1], 0.8 * st.size); m._xtexts.append({"runs": g["runs"], "color": mc, "box": g["box"]})
        for b in m.bonds:                                  # 3. a metal's lines to a perspective ring, with depth
            e = b.b if getattr(b.b, "_persp3d", False) else b.a if getattr(b.a, "_persp3d", False) else None
            if e is None or not b._prims or b._prims[0][0] != "line": continue
            r = next((r for r in persp if set(map(id, r["atoms"])) == set(map(id, e._eta))), None)
            if r is None: continue
            pr0 = b._prims[0]; P0, P1 = np.asarray(pr0[1][0], float), np.asarray(pr0[1][-1], float)
            if e is b.a: P0, P1 = P1, P0                   # P0 at the metal, P1 at the ring's centre
            Ce = self.P(e); tl = r.get("tilt", 0.0); L = norm(P1 - P0); u = (P1 - P0) / L; cuts = []
            elems = [(bb, pr) for bb in m.bonds if bb.a in e._eta and bb.b in e._eta for pr in bb._prims]
            ell = rings.get(id(r["atoms"][0]))
            if ell is not None: elems.append((None, ell[0]))
            for bb, pr in elems:
                outline = pr[0] == "poly" and len(pr) > 3 and isinstance(pr[3], dict)
                if outline: Q = np.array(pr[3]["spine"], float); w0, w1 = pr[3]["w"]
                elif pr[0] == "poly": Q = np.array([np.mean(pr[1][:2], axis=0), np.mean(pr[1][2:], axis=0)]); w0 = w1 = st.bold
                else: Q = np.asarray(pr[1], float); w0 = w1 = pr[2]
                for i in range(len(Q) - 1):
                    q0, q1 = Q[i], Q[i + 1]; den = cross(P1 - P0, q1 - q0)
                    if abs(den) < 1e-9: continue
                    t = cross(q0 - P0, q1 - q0) / den; sv = cross(q0 - P0, P1 - P0) / den
                    if not (0 < t < 1 and 0 <= sv <= 1): continue
                    X = P0 + (P1 - P0) * t; w = (w0 + (w1 - w0) * sv) / 2 if len(Q) == 2 else w0 / 2
                    if rotv(X - Ce, -tl)[1] < -0.01:       # near side (lower in the ring's own frame): the line goes behind
                        cuts.append((t * L, w + (0.9 if bb is not None else 0.6)))
                    else:                                  # far side: the line passes in front, so the element breaks
                        g = pr0[2] / 2 + 0.9
                        if bb is None:
                            parts = self._cut_poly(Q, X, g); out[:] = [x for x in out if x is not pr]
                            out += [("line", list(Pp), pr[2], pr[3], None, "round") for Pp in parts]
                            rings[id(r["atoms"][0])] = (None, ell[1])
                        elif outline:                      # cut the edge's polygon square across, g either side
                            ue = unit(q1 - q0); Lq = norm(q1 - q0); sx = float((X - q0) @ ue); new = []
                            for a0, a1 in ((None, sx - g), (sx + g, None)):
                                lo_, hi_ = (0.0 if a0 is None else a0), (Lq if a1 is None else a1)
                                if hi_ - lo_ < 0.5: continue
                                Pp = list(pr[1])
                                if a0 is not None: Pp = self._clip_half(Pp, -ue, -float((q0 + ue * a0) @ ue))
                                if a1 is not None: Pp = self._clip_half(Pp, ue, float((q0 + ue * a1) @ ue))
                                if len(Pp) < 3: continue
                                wl = lambda s_: w0 + (w1 - w0) * s_ / Lq
                                new.append(("poly", Pp, pr[2], {"spine": (q0 + ue * lo_, q0 + ue * hi_), "w": (wl(lo_), wl(hi_))}))
                            bb._prims = [x for x in bb._prims if x is not pr] + new
                        else:
                            parts = self._cut_poly(Q, X, g); bb._prims = [x for x in bb._prims if x is not pr]
                            bb._prims += [(pr[0], list(Pp)) + tuple(pr[2:]) for Pp in parts]
                    break
            keep, s0 = [], 0.0                             # the metal's line, hidden behind the ring's near side: from
            if cuts:                                       # its near edge to its ellipse, so only the tick at the centre shows
                lo = min(sv - g for sv, g in cuts); hi = max(sv + g for sv, g in cuts)
                if lo > 0.8: keep.append((0.0, lo))
                s0 = hi
            if L > s0 + 0.5: keep.append((s0, L))
            b._prims = [("line", [P0 + u * a0, P0 + u * a1]) + tuple(pr0[2:]) for a0, a1 in keep] + b._prims[1:]
        return out
    def _mol_box(self, m, curlies=True, atoms=None, prep=True):
        # once a build has drawn the marks (charges, lone pairs) round the arrows, a box is measured on them as drawn:
        # preparing the molecule again would re-place them without the arrows, and check() would then measure marks
        # that are not in the picture
        if prep and not getattr(self, "_drawn", False): self._prep_mol(m)
        xs, ys = [], []; ids = None if atoms is None else {id(a) for a in atoms}
        for a in (m.atoms if atoms is None else atoms):
            p = self.P(a); xs.append(p[0]); ys.append(p[1])
            for g in a._g or []: xs += [g["box"][0], g["box"][2]]; ys += [g["box"][1], g["box"][3]]
            if a.site and a.drawn: R = self._halo_r(a); xs += [p[0] - R, p[0] + R]; ys += [p[1] - R, p[1] + R]   # its halo
            for k in ("_chg", "_rad", "_pc", "_num", "_site", "_amb"):
                bk = self._mark_box(a, k)
                if bk: xs += [bk[0], bk[2]]; ys += [bk[1], bk[3]]
            for _, c, d1, d2 in a._lp:
                for d in (d1, d2): xs += [d[0] - self.st.lp_r, d[0] + self.st.lp_r]; ys += [d[1] - self.st.lp_r, d[1] + self.st.lp_r]
        for b in m.bonds:
            if ids is not None and (id(b.a) not in ids or id(b.b) not in ids): continue
            for pr in b._prims:
                for p in pr[1]: xs.append(p[0]); ys.append(p[1])
        if atoms is None:
            for pr in getattr(m, "_xprims", []):
                P = np.array(pr[1], float).reshape(-1, 2) if pr[0] in ("line", "poly") else np.array([pr[1] - pr[2], pr[1] + pr[2]])
                xs += list(P[:, 0]); ys += list(P[:, 1])
            for ex in getattr(self, "extras", []):         # its orbitals are part of it (a caption goes below the lobes)
                if ex[0] in ("wiggle", "dipole") and ex[1].mol is m:   # and its wiggles and dipole arrows
                    for pr in (self._wiggle_prims(ex) if ex[0] == "wiggle" else self._dipole_prims(ex)):
                        P = np.asarray(pr[1], float).reshape(-1, 2); xs += list(P[:, 0]); ys += list(P[:, 1])
                if ex[0] != "orbital" or ex[1].mol is not m: continue
                g_ = self._orbital_label(ex)               # and their names
                if g_: xs += [g_["box"][0], g_["box"][2]]; ys += [g_["box"][1], g_["box"][3]]
                for pr in self._orbital_prims(ex):
                    if pr[0] == "orb": P = np.asarray(pr[1], float)
                    elif pr[0] in ("circle", "dot"): P = np.array([pr[1] - pr[2], pr[1] + pr[2]])
                    else: continue
                    xs += list(P[:, 0]); ys += list(P[:, 1])
        if curlies:
            own = [c for c in self.curlies if self._owners(c) <= {id(m)}]
            if own:
                self._solve(own, [m], quick=True)
                for c in own:
                    pts = np.vstack([c._pts, c._head]); xs += list(pts[:, 0]); ys += list(pts[:, 1])
        return (min(xs), min(ys), max(xs), max(ys))
    def _owners(self, c):
        s = set()
        for x in (c.src, c.dst):
            if isinstance(x, LP): s.add(id(x.atom.mol))
            elif isinstance(x, Atom): s.add(id(x.mol))
            elif isinstance(x, (tuple, list, Mid)):
                for a in ((x.a, x.b) if isinstance(x, Mid) else x):
                    if isinstance(a, Atom): s.add(id(a.mol))
        return s

    # ---------- curly arrows
    # Every arrow is one circular arc (as ChemDraw draws them), or two smoothly joined arcs, running tangent into a
    # straight head. The head is aimed at its target: the centre of an atom label or vertex, a bond, or along a
    # forming bond.
    # Candidates that let the head or the shaft touch a bond, label, charge, lone pair or another arrow are
    # rejected outright; the rest are ranked by clearance and by the shape usual for their family.
    ARROW_PREF = {"hook": 160.0, "lpadj": 130.0, "adj": 100.0, "shift": 150.0, "inter": 55.0}

    def _concerted(self, c):
        """True if another arrow also forms the bond c forms (two arrows into the same Mid: they meet nose to nose)."""
        d = c.dst
        if not isinstance(d, Mid): return False
        pair = {d.a.uid, d.b.uid}
        return any(x is not c and isinstance(x.dst, Mid) and {x.dst.a.uid, x.dst.b.uid} == pair for x in self.curlies)
    def _edst(self, c):
        """Where an arrow's head is drawn. An arrow to Mid(a, b) between two unbonded atoms more than 1.3 BU apart,
        with no partner arrow forming the same bond, points at the far atom (the electrophile) like any intermolecular
        attack; a forming bond across a short gap, or one made by two arrows together, keeps the nose-to-nose midpoint."""
        d = c.dst
        if not isinstance(d, Mid) or self._concerted(c): return d
        a, b = d.a, d.b
        if (a.mol is b.mol and a.mol.find(a, b) is not None) or norm(a.xy - b.xy) <= 1.3: return d
        s = c.src; sa = [s.atom] if isinstance(s, LP) else [s] if isinstance(s, Atom) else [x for x in (s.a, s.b)] if isinstance(s, Mid) \
            else [x for x in s if isinstance(x, Atom)] if isinstance(s, (tuple, list)) else []
        ins = {x.uid for x in sa}
        if a.uid in ins and b.uid not in ins: return b     # the source holds one end (the group of a 1,2-shift, an O closing a
        if b.uid in ins and a.uid not in ins: return a     # ring onto a carbonyl): the head goes to the other
        fa = {x.uid for x in a.mol.frag(a)}; fb = {x.uid for x in b.mol.frag(b)}
        if ins and ins <= fa and not ins & fb: return b
        if ins and ins <= fb and not ins & fa: return a
        return d
    def _aclass(self, c):
        """Arrow family (shapes are kept alike within a family): 'hook' (a bond's electrons onto one of its own
        atoms), 'lpadj' (lone pair into a bond of the same atom), 'adj' (bond into a bond or forming bond that
        shares an atom with it), 'shift' (a bond onto an atom bonded to one of its ends: a 1,2-shift, or the C-H of an
        E1 step drawn onto the cation) and 'inter' (everything else)."""
        s, d = c.src, self._edst(c)
        dd = (d.a, d.b) if isinstance(d, Mid) else tuple(d) if isinstance(d, (tuple, list)) else ()
        if isinstance(s, (tuple, list)) and isinstance(d, Atom) and any(d is x for x in s): return "hook"
        if isinstance(s, LP) and any(s.atom is x for x in dd): return "lpadj"
        if isinstance(s, (tuple, list)) and any(x is y for x in s for y in dd): return "adj"
        if self._shift_end(s, d) is not None: return "shift"
        return "inter"
    @staticmethod
    def _shift_end(s, d):
        """For an arrow from bond s onto atom d bonded to exactly one end of it: (that end, the other end), else None."""
        if not (isinstance(s, (tuple, list)) and len(s) == 2 and all(isinstance(x, Atom) for x in s) and isinstance(d, Atom)):
            return None
        m = s[0].mol
        if d.mol is not m or any(d is x for x in s) or m.find(*s) is None: return None
        adj = [x for x in s if m.find(d, x) is not None]
        return (adj[0], s[1] if adj[0] is s[0] else s[0]) if len(adj) == 1 else None
    @staticmethod
    def _hit(c):
        """The atoms an arrow attacks: its target atom, or the two ends of the bond it forms."""
        return (c.dst.a, c.dst.b) if isinstance(c.dst, Mid) else (c.dst,) if isinstance(c.dst, Atom) else ()
    def _curl(self, c, cd):
        """For an arrow that curls back beside its source bond, a 'hook' (X-Y onto X), a 'shift' (X-Y onto Z bonded to X) or
        a bond into the next bond (X-Y into X-Z): (midpoint of the bond, the target, the far end Y, shift?, Y moves?),
        points in pt, else None. Such an arrow may
        leave its bond leaning towards Y, away from where it turns back to; a shift starts on the side facing Z, and when
        Y moves with the electrons (a 1,2-shift, not an E1 proton taken by a base) the lean is its usual shape."""
        s = c.src
        if not (isinstance(s, (tuple, list)) and len(s) == 2 and all(isinstance(x, Atom) for x in s)): return None
        b = s[0].mol.find(*s)
        if b is None: return None
        if isinstance(cd, (tuple, list)) and len(cd) == 2 and all(isinstance(x, Atom) for x in cd) and cd[0].mol.find(*cd) is not None:
            sh = [x for x in s if any(x is y for y in cd)]      # a bond into the next bond (C-H into C-C=O): may lean towards
            if len(sh) != 1: return None                      # its far end and curl back into the corner, like a hook
            far = s[1] if sh[0] is s[0] else s[0]
            return (np.asarray(b._mid, float), self.P(sh[0]), self.P(far), False, False)
        if not isinstance(cd, Atom): return None
        if any(cd is x for x in s):
            far = s[1] if cd is s[0] else s[0]
            return (np.asarray(b._mid, float), self.P(cd), self.P(far), False, False)
        e = self._shift_end(s, cd)
        if e is None: return None
        taken = {id(y) for x in self.curlies if x is not c for y in self._hit(x)}
        moves = id(e[1]) not in taken and b.order == 1
        return (np.asarray(b._mid, float), self.P(cd), self.P(e[1]), True, moves)
