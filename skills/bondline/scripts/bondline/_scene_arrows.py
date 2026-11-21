"""bondline engine: Scene methods, part 2 of 3."""
# Part of the bondline engine (see __init__.py): the package's modules share one namespace, so names
# defined in any of them are visible here at call time, as in a single module.
from ._core import *  # noqa: F401,F403


class _SceneArrows:
    @classmethod
    def _migrant(cls, s, d, others=()):
        """The group of a 1,2-shift drawn as an arrow from the single bond s onto atom d bonded to one end of it: the atom
        at the other end, which moves to d with the electrons. None if that is not the case, or if another arrow takes
        that atom (its uid in `others`: the base removing the H of an E1 step, whose C-H electrons stay behind)."""
        e = cls._shift_end(s, d)
        if e is None or s[0].mol.find(*s).order != 1 or e[1].uid in others: return None
        return e[1]
    def _obst(self, mols, skip=frozenset()):
        """Obstacles in pt: bond line segments (owner = bond), glyph outlines (owner = atom) and mark circles (owner =
        atom for charges, radicals, partial charges and numbers; the LP for lone-pair dots). Marks of atoms in
        `skip` are left out, apart from lone pairs that start an arrow (those marks are re-placed afterwards)."""
        st = self.st; S, Sid, G, Gid, C, Cid = [], [], [], [], [], []
        for m in mols:
            for b in m.bonds:
                for pr in b._prims:
                    q = pr[1]; k = len(q) if pr[0] == "poly" else len(q) - 1
                    for i in range(k): S.append((*q[i], *q[(i + 1) % len(q)])); Sid.append(id(b))
            for a in m.atoms:
                for g in a._g or []:
                    for H in g["hulls"]: G.append(H); Gid.append(id(a))
        srcs = {id(c.src) for c in self.curlies if isinstance(c.src, LP)}
        for m in mols:
            for a in m.atoms:
                mv = id(a) in skip
                for mk in (a._chg, a._rad, a._pc, a._num, a._site):
                    if mk and not mv: C.append((*mk[0], mk[1])); Cid.append(id(a))
                for l, c0, d1, d2 in a._lp:
                    if mv and id(l) not in srcs: continue
                    C.append((*d1, st.lp_r)); Cid.append(id(l)); C.append((*d2, st.lp_r)); Cid.append(id(l))
        f = lambda L, k: np.array(L, float).reshape(-1, k)
        rings = []
        for m in mols:
            seen = set()
            for b in m.bonds:
                r = m.ring_of(b)
                if r and frozenset(x.uid for x in r) not in seen:
                    seen.add(frozenset(x.uid for x in r)); rings.append(np.array([self.P(x) for x in r]))
        HS = []                                   # bond centre lines hidden under a label's clearance: an arrow may
        for m in mols:                            # not slip between a label and its bond (it would cross the bond)
            for b in m.bonds:
                vis = getattr(b, "_vis", None)
                if vis is None or b.style in ("dash", "dot") or getattr(b.a, "_eta", None) or getattr(b.b, "_eta", None): continue
                for X, V in ((b.a, vis[0]), (b.b, vis[1])):
                    if (X.drawn or X.clear is not None) and not X.hidden and norm(V - self.P(X)) > 0.5: HS.append((*self.P(X), *V))
        return {"S": f(S, 4), "Sid": np.array(Sid, np.int64), "Gp": G, "Gid": np.array(Gid, np.int64),
                "Gbox": f([(*H.min(0), *H.max(0)) for H in G], 4), "C": f(C, 3), "Cid": np.array(Cid, np.int64),
                "rings": rings, "AS": np.zeros((0, 4)), "HS": f(HS, 4)}
    @staticmethod
    def _osub(O, lo, hi):
        """The obstacles that reach into the box lo..hi (pt); glyph outlines become one edge list."""
        def bb(X): return np.column_stack([np.minimum(X[:, 0], X[:, 2]), np.minimum(X[:, 1], X[:, 3]),
                                           np.maximum(X[:, 0], X[:, 2]), np.maximum(X[:, 1], X[:, 3])]) if len(X) else np.zeros((0, 4))
        def ov(B): return (B[:, 2] >= lo[0]) & (B[:, 0] <= hi[0]) & (B[:, 3] >= lo[1]) & (B[:, 1] <= hi[1])
        ks = ov(bb(O["S"])); kg = ov(O["Gbox"]); ka = ov(bb(O["AS"])); C = O["C"]; HS = O.get("HS", np.zeros((0, 4))); kh = ov(bb(HS))
        kc = ov(np.column_stack([C[:, :2] - C[:, 2:], C[:, :2] + C[:, 2:]])) if len(C) else np.zeros(0, bool)
        E = [np.hstack([O["Gp"][i], np.roll(O["Gp"][i], -1, axis=0)]) for i in np.flatnonzero(kg)]
        st = np.cumsum([0] + [len(e) for e in E[:-1]]).astype(int) if E else np.zeros(0, int)
        return {"S": O["S"][ks], "Sid": O["Sid"][ks], "G": np.vstack(E) if E else np.zeros((0, 4)), "Gst": st,
                "Gid": O["Gid"][kg], "C": C[kc], "Cid": O["Cid"][kc], "AS": O["AS"][ka], "rings": O["rings"], "HS": HS[kh]}
    def _bond_side_pt(self, b, s, gap, t=0.5):
        ext = max(s * o for o in b._offs) + self.st.lw / 2
        p, q = getattr(b, "_vis", (b._mid, b._mid))
        return p + (q - p) * t + perp(b._u) * s * (ext + gap)
    def _pt(self, x):
        if isinstance(x, LP): return next(t[1] for t in x.atom._lp if t[0] is x)
        if isinstance(x, Atom): return self.P(x)
        if isinstance(x, Mid): return self.P(x.a) * (1 - x.t) + self.P(x.b) * x.t
        if isinstance(x, (tuple, list)) and len(x) == 2 and isinstance(x[0], Atom):
            b = x[0].mol.find(*x); return b._mid if b is not None else (self.P(x[0]) + self.P(x[1])) / 2
        return np.array(x, float) * self.st.bond
    def _tails(self, x, dpt, tgt=()):
        """Tail options: (start point, preferred departure direction, info). A lone pair may move to any free
        direction where it still clearly belongs to its atom (it may point at the arrow's own target atoms `tgt`)."""
        st = self.st; g0 = st.alw / 2 + 1.0 * st.bond / 14.4; out = []       # tail starts 1 pt clear of its pair or bond
        if isinstance(x, LP):
            a = x.atom; A = self.P(a); angs = [x.angle]
            if x.auto and not getattr(x, "_locked", False) and self._lobe_at(a, x.angle) is None:   # (a lobe's pair stays)
                occ = self._occ(a) + [mk[2] for mk in (a._chg, a._rad) if mk]
                po = [l.angle for l in a.lps if l is not x and l.angle is not None]   # the atom's other pairs: 90 deg off
                sep = 100 if len(a.lps) == 2 else 88                                  # (100 beside a single other one),
                angs += [t for t in range(0, 360, 15) if min([abs(adiff(t, o)) for o in occ] or [180]) >= 50   # so each
                         and min([abs(adiff(t, o)) for o in po] or [180]) >= sep                               # reads as a pair
                         and abs(adiff(t, x.angle)) > 1 and not self._in_ring(a, t)
                         and self._owner_ok(a, self._lp_geom(a, t)[0], skip=tgt)]
            for t in angs:
                u = dirv(t); cen = self._lp_geom(a, t)[0]
                out.append((cen + u * (st.lp_r + g0), u, {"lp": (x, t, abs(adiff(t, x.angle)))}))
        elif isinstance(x, Atom):
            mk = x._rad or x._chg
            for r in (0, 30, -30, 60, -60):
                if mk: u = rotv(unit(dpt - mk[0]), r); out.append((mk[0] + u * (mk[1] + g0), u, {}))
                else:
                    u = rotv(unit(dpt - self.P(x)), r); rr = (self._exit(x, u) + g0) if x.drawn else 0.14 * st.bond
                    out.append((self.P(x) + u * rr, u, {}))
        elif isinstance(x, Mid): out.append((self._pt(x), None, {}))
        elif isinstance(x, (tuple, list)) and len(x) == 2 and isinstance(x[0], Atom):
            b = x[0].mol.find(*x)
            if b is None: return self._tails(Mid(*x), dpt)
            n = perp(b._u)
            for s in (1, -1):
                for t in (0.5, 0.35, 0.65, 0.3, 0.7): out.append((self._bond_side_pt(b, s, g0, t), n * s, {"bond": b, "t": t}))
        else: out.append((np.array(x, float) * st.bond, None, {}))
        return out
    def _tips(self, x, tref, wide=False):
        """Tip options: (tip point, [(head direction, aim error deg)], info). Atom tips sit just outside the label
        (or near the vertex) with the head pointing at the atom's centre, give or take 24 deg."""
        st = self.st; out = []
        if isinstance(x, Atom):
            A = self.P(x); occ = [ang(n.xy - x.xy) for n in x.mol.nbrs(x)]
            ext = {"right": 0, "left": 180, "down": -90, "up": 90}.get(getattr(x, "_dir", None)) if (x._g and len(x._g) > 1) else None
            lim = 24 if x.drawn else 34
            for phi in range(0, 360, 10):
                if any(abs(adiff(phi, o)) < lim for o in occ): continue
                if ext is not None and abs(adiff(phi, ext)) < 40: continue
                u = dirv(phi); r = self._exit(x, u, margin=st.tip_gap) if x.drawn else 0.22 * st.bond
                out.append((A + u * r, [(rotv(-u, dl), abs(dl)) for dl in (0, 12, -12, 24, -24)], {"phi": phi}))
        elif isinstance(x, Mid):
            # a forming bond: the tip stops just on the source's side of the midpoint and the head points along the
            # new bond, from the source's end towards the other atom, so two arrows meeting there end nose to nose
            M = self._pt(x); u = unit(self.P(x.b) - self.P(x.a)); sd = 1.0 if float((tref - M) @ u) >= 0 else -1.0
            out.append((M + u * sd * 0.09 * st.bond, [(rotv(-u * sd, r), abs(r) / 3) for r in (0, 15, -15, 30, -30, 45, -45)], {}))
            if wide:                                  # one arrow alone into a forming bond may also come in from the side
                out.append((M, [(rotv(-u * sd, r), abs(r) / 3) for r in (60, -60, 75, -75, 90, -90)], {}))
        elif isinstance(x, (tuple, list)) and len(x) == 2 and isinstance(x[0], Atom):
            b = x[0].mol.find(*x)
            if b is None: return self._tips(Mid(*x), tref)
            n = perp(b._u)
            for s in (1, -1):
                for t in (0.5, 0.38, 0.62, 0.3, 0.7):
                    out.append((self._bond_side_pt(b, s, st.tip_gap, t), [(rotv(-n * s, dl), abs(dl) / 2) for dl in (0, 20, -20, 40, -40)], {}))
        else: out.append((np.array(x, float) * st.bond, [(dirv(t), 0.0) for t in range(0, 360, 15)], {}))
        return out
    @staticmethod
    def _arc_par(T, B, d):
        """Circle through T and B that arrives at B travelling along unit d: (centre, radius, ccw?, start angle,
        sweep in radians, unit departure at T), centre None for a straight line; None if impossible."""
        wx, wy = T[0] - B[0], T[1] - B[1]; nx, ny = -d[1], d[0]; den = 2.0 * (nx * wx + ny * wy); L2 = wx * wx + wy * wy
        if L2 < 1e-12: return None
        if abs(den) < 1e-9 * math.sqrt(L2):
            if -(wx * d[0] + wy * d[1]) <= 0: return None
            return (None, math.inf, True, 0.0, 0.0, np.array(d, float))
        rho = L2 / den; cx, cy = B[0] + rho * nx, B[1] + rho * ny; ccw = rho > 0
        thT = math.atan2(T[1] - cy, T[0] - cx); thB = math.atan2(B[1] - cy, B[0] - cx)
        sw = (thB - thT) % (2 * math.pi) if ccw else (thT - thB) % (2 * math.pi); s = 1.0 if ccw else -1.0
        return ((cx, cy), abs(rho), ccw, thT, sw, np.array([-math.sin(thT) * s, math.cos(thT) * s]))
    @staticmethod
    def _arc_par_f(Tx, Ty, Bx, By, d0, d1):
        """_arc_par on plain floats (the same arithmetic in the same order, without numpy scalar overhead)."""
        wx, wy = Tx - Bx, Ty - By; nx, ny = -d1, d0; den = 2.0 * (nx * wx + ny * wy); L2 = wx * wx + wy * wy
        if L2 < 1e-12: return None
        if abs(den) < 1e-9 * math.sqrt(L2):
            if -(wx * d0 + wy * d1) <= 0: return None
            return (None, math.inf, True, 0.0, 0.0, np.array((d0, d1), float))
        rho = L2 / den; cx, cy = Bx + rho * nx, By + rho * ny; ccw = rho > 0
        thT = math.atan2(Ty - cy, Tx - cx); thB = math.atan2(By - cy, Bx - cx)
        sw = (thB - thT) % (2 * math.pi) if ccw else (thT - thB) % (2 * math.pi); s = 1.0 if ccw else -1.0
        return ((cx, cy), abs(rho), ccw, thT, sw, np.array([-math.sin(thT) * s, math.cos(thT) * s]))
    @classmethod
    def _biarc_par_f(cls, Tx, Ty, t0x, t0y, Bx, By, t1x, t1y):
        """_biarc_par on plain floats (same arithmetic, same order); the joint comes back as an array, as before."""
        dx, dy = Bx - Tx, By - Ty; vx, vy = t0x + t1x, t0y + t1y
        a = 2.0 * (t0x * t1x + t0y * t1y - 1.0); dv = dx * vx + dy * vy; dd = dx * dx + dy * dy
        if dd < 1e-12: return None
        if abs(a) < 1e-9:
            if dv <= 1e-9: return None
            al = dd / (2.0 * dv)
        else: al = (2.0 * dv - math.sqrt(max(0.0, 4.0 * dv * dv - 4.0 * a * dd))) / (2.0 * a)
        if not al > 1e-9: return None
        Jx = (Tx + al * t0x + Bx - al * t1x) / 2.0; Jy = (Ty + al * t0y + By - al * t1y) / 2.0
        p1 = cls._arc_par_f(Jx, Jy, Tx, Ty, -t0x, -t0y); p2 = cls._arc_par_f(Jx, Jy, Bx, By, t1x, t1y)
        if p1 is None or p2 is None: return None
        s1 = None if p1[0] is None else (not p1[2]); s2 = None if p2[0] is None else p2[2]
        if s1 is not None and s2 is not None and s1 != s2:
            if min(p1[4], p2[4]) > math.radians(4): return None                     # S-bend
            sw = abs(p1[4] - p2[4]); ccw = s1 if p1[4] > p2[4] else s2
        else: sw = p1[4] + p2[4]; ccw = s2 if s2 is not None else (s1 if s1 is not None else True)
        L = sum((p[1] * p[4]) if p[0] is not None else math.hypot(qx - Jx, qy - Jy) for p, qx, qy in ((p1, Tx, Ty), (p2, Bx, By)))
        return (np.array([Jx, Jy]), p1, p2, math.degrees(sw), min(p1[1], p2[1]), ccw, L, max(p1[1], p2[1]))
    @staticmethod
    def _arc_pts(T, B, par, n):
        C, R, ccw, thT, sw, _ = par
        if C is None: return T + (B - T) * np.linspace(0, 1, n)[:, None]
        th = thT + (sw if ccw else -sw) * np.linspace(0, 1, n)
        return np.column_stack([C[0] + R * np.cos(th), C[1] + R * np.sin(th)])
    @classmethod
    def _biarc_par(cls, T, t0, B, t1):
        """Two circular arcs joined smoothly, leaving T along unit t0 and arriving at B along unit t1 (equal tangent
        lengths): (joint, arc 1 params (travel joint -> T), arc 2 params, sweep deg, smallest radius, ccw?, length,
        largest radius), or None for S-bends and impossible cases."""
        dx, dy = B[0] - T[0], B[1] - T[1]; vx, vy = t0[0] + t1[0], t0[1] + t1[1]
        a = 2.0 * (t0[0] * t1[0] + t0[1] * t1[1] - 1.0); dv = dx * vx + dy * vy; dd = dx * dx + dy * dy
        if dd < 1e-12: return None
        if abs(a) < 1e-9:
            if dv <= 1e-9: return None
            al = dd / (2.0 * dv)
        else: al = (2.0 * dv - math.sqrt(max(0.0, 4.0 * dv * dv - 4.0 * a * dd))) / (2.0 * a)
        if not al > 1e-9: return None
        J = (T + al * t0 + B - al * t1) / 2.0
        p1 = cls._arc_par(J, T, -t0); p2 = cls._arc_par(J, B, t1)
        if p1 is None or p2 is None: return None
        s1 = None if p1[0] is None else (not p1[2]); s2 = None if p2[0] is None else p2[2]
        if s1 is not None and s2 is not None and s1 != s2:
            if min(p1[4], p2[4]) > math.radians(4): return None                     # S-bend
            sw = abs(p1[4] - p2[4]); ccw = s1 if p1[4] > p2[4] else s2
        else: sw = p1[4] + p2[4]; ccw = s2 if s2 is not None else (s1 if s1 is not None else True)
        L = sum((p[1] * p[4]) if p[0] is not None else norm(q - J) for p, q in ((p1, T), (p2, B)))
        return (J, p1, p2, math.degrees(sw), min(p1[1], p2[1]), ccw, L, max(p1[1], p2[1]))
    @staticmethod
    def _par_shift(par, o):
        """Arc parameters (single arc or joined pair) moved by o."""
        if len(par) == 6:                                       # single arc: (centre, R, ccw, thT, sweep, departure)
            return par if par[0] is None else ((par[0][0] + o[0], par[0][1] + o[1]),) + tuple(par[1:])
        J, p1, p2 = par[:3]
        return (J + o, Scene._par_shift(p1, o), Scene._par_shift(p2, o)) + tuple(par[3:])
    @staticmethod
    def _arcs_pts(T, B, pars, n):
        """_arc_pts for many arcs at once (same arithmetic, element by element): rows of (n, 2) points."""
        N = len(pars); out = np.empty((N, n, 2)); lin = np.linspace(0, 1, n)
        cur = [i for i, p in enumerate(pars) if p[0] is not None]; st_ = [i for i, p in enumerate(pars) if p[0] is None]
        if cur:
            C = np.array([pars[i][0] for i in cur], float); R = np.array([pars[i][1] for i in cur], float)
            th0 = np.array([pars[i][3] for i in cur], float); ssw = np.array([pars[i][4] if pars[i][2] else -pars[i][4] for i in cur], float)
            th = th0[:, None] + ssw[:, None] * lin[None, :]
            out[cur, :, 0] = C[:, 0:1] + R[:, None] * np.cos(th); out[cur, :, 1] = C[:, 1:2] + R[:, None] * np.sin(th)
        for i in st_: out[i] = T[i] + (B[i] - T[i]) * lin[:, None]
        return out
    def _cand_pts(self, out, h, n):
        """Shaft points of every candidate (n per candidate), computed in batches; identical to _arc_pts/_biarc_pts."""
        if not out: return
        one = [k for k in out if not k["bi"]]; two = [k for k in out if k["bi"]]
        if one:
            T = np.array([k["T"] for k in one], float); B = np.array([k["P"] - k["d"] * h for k in one], float)
            P = self._arcs_pts(T, B, [k["par"] for k in one], n)
            for k, q in zip(one, P): k["pts"] = q
        if two:
            m = max(4, n // 2)
            J = np.array([k["par"][0] for k in two], float); T = np.array([k["T"] for k in two], float)
            B = np.array([k["P"] - k["d"] * h for k in two], float)
            A1 = self._arcs_pts(J, T, [k["par"][1] for k in two], m); A2 = self._arcs_pts(J, B, [k["par"][2] for k in two], n - m + 1)
            for k, a1, a2 in zip(two, A1, A2): k["pts"] = np.vstack([a1[::-1], a2[1:]])
    @classmethod
    def _biarc_pts(cls, T, B, par, n):
        J, p1, p2 = par[:3]; m = max(4, n // 2)
        return np.vstack([cls._arc_pts(J, T, p1, m)[::-1], cls._arc_pts(J, B, p2, n - m + 1)[1:]])
    def _head_poly(self, P, d, kind, sb=1):
        st = self.st; h, w = st.c_head_len, st.c_head_w; nn = perp(d)
        if kind == "full": return np.array([P, P - d * h + nn * w, P - d * st.head_notch * h, P - d * h - nn * w])
        return np.array(self._half_head(P, d, h, w * 1.25, st.head_notch, st.alw / 2, sb))   # fishhook: one barb, sharp
    def _head_samples(self, P, d, kind, sb=1):
        return self._heads(np.asarray(P)[None], np.asarray(d)[None], kind, np.array([sb]))[0]
    def _heads(self, P, d, kind, sb):
        """Outline samples (N, 14, 2) of N arrowheads: corners, points along each edge and the centre."""
        st = self.st; h, w = st.c_head_len, st.c_head_w; n = np.column_stack([-d[:, 1], d[:, 0]])
        if kind == "full": H = np.stack([P, P - d * h + n * w, P - d * st.head_notch * h, P - d * h - n * w], 1)
        else: H = np.stack([P, P - d * h + n * (w * 1.25 * sb)[:, None], P - d * st.head_notch * h], 1)
        E = np.roll(H, -1, axis=1); t = np.array([0.0, 0.33, 0.67])[None, :, None, None]
        pts = (H[:, None] * (1 - t) + E[:, None] * t).reshape(len(P), -1, 2)
        pts = np.concatenate([pts, H.mean(1)[:, None]], 1)
        return np.concatenate([pts, np.repeat(P[:, None], 14 - pts.shape[1], 1)], 1) if pts.shape[1] < 14 else pts[:, :14]
    def _arrow_cands(self, c, relax=False):
        """Candidate arrows: a single arc (ChemDraw's default) or, to leave the source and reach the target in set
        directions, two joined arcs. Tails leave lone pairs along their axis (within 15 deg) and bonds at right angles
        (within 30 deg)."""
        st = self.st; h = st.c_head_len; cd = self._edst(c); dpt = self._pt(cd); out = []
        tgt = (cd,) if isinstance(cd, Atom) else (cd.a, cd.b) if isinstance(cd, Mid) else tuple(cd) if isinstance(cd, (tuple, list)) else ()
        tails = self._tails(c.src, dpt, tgt); tips = self._tips(cd, tails[0][0] if tails else dpt, wide=isinstance(cd, Mid) and not self._concerted(c))
        force = None if c.side is None else (1 if c.side == "left" else -1)
        swmax, rmin, lmin = (300.0, 0.10 * st.bond, 0.12 * st.bond) if relax else (250.0, 0.15 * st.bond, 0.22 * st.bond)
        swmin = 12.0 if relax else 25.0                       # an electron-pushing arrow always reads as curved
        heads = [(P, d, aim, P - d * h) for P, dirs, pi in tips for d, aim in dirs]      # head bases, computed once
        # the same tails and tips again (a second pass, or a rebuild for check()): reuse the candidates
        key = None; lean = self._curl(c, cd) if not relax else None
        if tails and heads:                                     # exact inputs only, so a reuse is bit-identical
            o = np.zeros(2); rk = lambda v: np.asarray(v, float).tobytes()
            key = (relax, force, h, st.bond, None if lean is None else (rk(lean[0]), rk(lean[1]), rk(lean[2]), lean[3], lean[4]),
                   tuple((rk(T), None if u0 is None else rk(u0), "lp" in ti, "bond" in ti) for T, u0, ti in tails),
                   tuple((rk(P), rk(d), aim) for P, d, aim, B in heads))
            cache = self.__dict__.setdefault("_ccache", {}); hit = cache.get(key)
            if hit is not None:                                 # (side= filters a reused set too)
                return self._side_filter([{"T": tails[i][0], "P": heads[j][0], "d": heads[j][1], "par": self._par_shift(par, o), "sw": sw, "R": R,
                                           "L": L, "dep": dep, "aim": heads[j][2], "ccw": ccw, "ti": tails[i][2], "bi": bi, "even": ev, "pts": pts + o}
                                          for i, j, par, sw, R, L, dep, ccw, bi, ev, pts in hit], force)
        idx = []
        for ii, (T, u0, ti) in enumerate(tails):
            # leave a lone pair square to the pair (along its axis, so the arrow visibly starts from its middle) and a
            # bond at right angles to it
            depmax = 150.0 if relax else (15.0 if "lp" in ti else 30.0 if "bond" in ti else 60.0)
            rs = (0, 15, -15) if "lp" in ti else (0, 20, -20) if "bond" in ti else (0, 30, -30)
            tow, lmax = 0.0, depmax                # (a departure leaning towards tow may go as far as lmax)
            if lean is not None and "bond" in ti:
                # an arrow curling back beside its own bond leaves it leaning away from where it turns back to (towards
                # the far end Y), so it can come round without a kink: a hook onto a crowded vertex, a 1,2-shift. A shift
                # starts on the side of the bond facing the atom it goes to (the far side would wrap round the bond); when
                # the group moves with the electrons, leaning towards it is the usual shape, otherwise only if need be
                mid, tg, far, shift, moves = lean
                if shift and float((T - mid) @ (tg - mid)) <= 0: continue
                tow = 1.0 if cross(u0, far - T) > 0 else -1.0
                if moves: u0 = rotv(u0, 40.0 * tow); rs = (0, 20, -20, 40, -40); depmax = lmax = 40.0
                else: rs = (0, 20, -20, 40 * tow, 60 * tow); lmax = 65.0
            deps = [None] + ([(rotv(u0, r), float(abs(r))) for r in rs] if u0 is not None else [])
            Tf = (float(T[0]), float(T[1]))
            for jj, (P, d, aim, B) in enumerate(heads):
                Bf = (float(B[0]), float(B[1])); df = (float(d[0]), float(d[1]))
                if True:
                    for dep0 in deps:
                        if dep0 is None:
                            p = self._arc_par_f(Tf[0], Tf[1], Bf[0], Bf[1], df[0], df[1])
                            if p is None: continue
                            sw = math.degrees(p[4]); R = p[1]; ccw = p[2]; Ls = R * p[4] if p[0] is not None else norm(B - T); even = 1.0
                            dep = 0.0 if u0 is None else math.degrees(math.acos(max(-1.0, min(1.0, float(p[5] @ u0)))))
                            dlim = lmax if tow and cross(u0, p[5]) * tow > 0 else depmax
                        else:
                            p = self._biarc_par_f(Tf[0], Tf[1], float(dep0[0][0]), float(dep0[0][1]), Bf[0], Bf[1], df[0], df[1])
                            if p is None: continue
                            sw, R, ccw, Ls = p[3], p[4], p[5], p[6]; dep = dep0[1]
                            even = R / p[7] if math.isfinite(p[7]) else 0.0      # 1 = one even curve
                            dlim = lmax if tow and cross(u0, dep0[0]) * tow > 0 else depmax
                        if sw > swmax or sw < swmin or R < rmin or dep > dlim or Ls < lmin: continue
                        out.append({"T": T, "P": P, "d": d, "par": p, "sw": sw, "R": R, "L": Ls + h, "dep": dep, "aim": aim,
                                    "ccw": ccw, "ti": ti, "bi": dep0 is not None, "even": even}); idx.append((ii, jj))
        self._cand_pts(out, h, 26)
        if key is not None:
            cache[key] = [(i, j, self._par_shift(k["par"], -o), k["sw"], k["R"], k["L"], k["dep"], k["ccw"], k["bi"], k["even"], k["pts"] - o)
                          for (i, j), k in zip(idx, out)]
        return self._side_filter(out, force)
    @staticmethod
    def _side_filter(out, force):
        """side=: keep the candidates that bulge to that side of their tail-to-head chord."""
        if not force: return out
        return [k for k in out if k["sw"] <= 1 or (1 if cross(k["P"] - k["T"], k["pts"][13] - k["T"]) > 0 else -1) == force]
    def _arrow_rel(self, c):
        """What the arrow may come close to: its own lone pair, source bond and source atom (near the tail only)."""
        s = c.src; sb = -1; sa = None; slp = -1
        if isinstance(s, LP): sa, slp = s.atom, id(s)
        elif isinstance(s, Atom): sa = s
        elif isinstance(s, (tuple, list)) and len(s) == 2 and isinstance(s[0], Atom):
            b = s[0].mol.find(*s); sb = id(b) if b is not None else -1
        sab = {id(b) for b in sa.mol.bonds if sa in (b.a, b.b)} if sa is not None else set()
        return sb, (id(sa) if sa is not None else -1), slp, sab
    def _arrow_req(self, c, O):
        """Required clearances (pt, beyond line edges and mark radii) for shaft points, tail-zone points and head points.
        -1e9 means 'ignore'."""
        st = self.st; lw2 = st.lw / 2; aw2 = st.alw / 2; k = st.bond / 14.4; sb, sa, slp, sab = self._arrow_rel(c); X = -1e9
        Sid, Gid, Cid = O["Sid"], O["Gid"], O["Cid"]; nA = len(O["AS"])
        srcb = Sid == sb; atb = np.isin(Sid, list(sab)) if sab else np.zeros(len(Sid), bool)
        own = Cid == slp; satm = Cid == sa; sgl = Gid == sa
        # what the head is aimed at may come as close as the tip gap; everything else keeps 1 pt of clear white space
        cd = self._edst(c); ta = cd if isinstance(cd, Atom) else None
        tb = cd[0].mol.find(*cd) if isinstance(cd, (tuple, list)) and len(cd) == 2 and isinstance(cd[0], Atom) else None
        tgtb = np.isin(Sid, [id(b) for b in ta.mol.bonds if ta in (b.a, b.b)] if ta is not None else ([id(tb)] if tb is not None else []))
        tgtg = Gid == (id(ta) if ta is not None else -1)
        return [(lambda Q: _segs_dist(Q, O["S"]), np.full(len(Sid), lw2 + aw2 + 1.0 * k),
                 np.where(srcb | atb, X, lw2 + aw2 + 0.8 * k), np.where(tgtb, lw2 + 0.75 * k, lw2 + 1.0 * k)),
                (lambda Q: _polys_dist(Q, O["G"], O["Gst"]), np.full(len(Gid), aw2 + 1.0 * k), np.where(sgl, X, 0.5 * k),
                 np.where(tgtg, 0.45 * k, 1.0 * k)),
                (lambda Q: np.hypot(Q[:, None, 0] - O["C"][None, :, 0], Q[:, None, 1] - O["C"][None, :, 1]) - O["C"][None, :, 2],
                 np.where(own, X, aw2 + 0.75 * k), np.where(own | satm, X, 0.5 * k), np.where(own, X, 0.9 * k)),
                (lambda Q: _segs_dist(Q, O["AS"]), np.full(nA, 2 * aw2 + 1.0 * k), np.full(nA, 2 * aw2 + 0.8 * k), np.full(nA, aw2 + 1.0 * k))]
    def _arrow_eval(self, c, K, O, pref, wshape=1.0, coarse=False):
        """Hard-constraint violation V (0 = nothing touched) and soft cost for each candidate in K."""
        st = self.st; N = len(K)
        SH = np.stack([k["pts"] for k in K]); ns = SH.shape[1]
        HD = self._heads(np.array([k["P"] for k in K]), np.array([k["d"] for k in K]), c.kind,
                         np.array([-1.0 if k["ccw"] else 1.0 for k in K])); nh = HD.shape[1]
        if coarse: SH = SH[:, ::3]; HD = HD[:, [0, 1, 2, 3, 12]]; ns = SH.shape[1]; nh = HD.shape[1]
        seg = np.hypot(*np.diff(SH, axis=1).transpose(2, 0, 1)); TZ = np.concatenate([np.zeros((N, 1)), np.cumsum(seg, 1)], 1) < 2.5
        allp = np.vstack([SH.reshape(-1, 2), HD.reshape(-1, 2)]); O = self._osub(O, allp.min(0) - 6.0, allp.max(0) + 6.0)
        fam = [f for f in self._arrow_req(c, O) if len(f[1])]
        V = np.zeros(N); CS = np.zeros(N); CH = np.zeros(N)
        for i0 in range(0, N, 48):
            n = min(N, i0 + 48) - i0; Qs = SH[i0:i0 + n].reshape(-1, 2); tz = TZ[i0:i0 + n].reshape(-1); Qh = HD[i0:i0 + n].reshape(-1, 2)
            ms = np.full(len(Qs), 99.0); mh = np.full(len(Qh), 99.0); vs = np.zeros(len(Qs)); vh = np.zeros(len(Qh))
            for dist, rm, rt, rh in fam:
                sl = dist(Qs) - np.where(tz[:, None], rt[None, :], rm[None, :])
                vs += (np.minimum(sl, 0) ** 2).sum(1); ms = np.minimum(ms, sl.min(1))
                sl = dist(Qh) - rh[None, :]; vh += (np.minimum(sl, 0) ** 2).sum(1); mh = np.minimum(mh, sl.min(1))
            V[i0:i0 + n] = vs.reshape(n, ns).sum(1) + 2.0 * vh.reshape(n, nh).sum(1) + 4.0 * self._crossings(SH[i0:i0 + n], O["HS"])
            CS[i0:i0 + n] = (np.maximum(0, 2.0 - ms) ** 2).reshape(n, ns).mean(1)
            CH[i0:i0 + n] = (np.maximum(0, 1.5 - mh) ** 2).reshape(n, nh).mean(1)
        ring = np.zeros(N)
        if O["rings"]:
            from matplotlib.path import Path as _Path
            probe = SH[:, [ns // 4, ns // 2, 3 * ns // 4]].reshape(-1, 2)
            inside = np.zeros(len(probe), bool)
            for R in O["rings"]: inside |= _Path(R).contains_points(probe)
            ring = 1.5 * inside.reshape(N, 3).sum(1)
        g = lambda key: np.array([k[key] for k in K], float)
        lpd = np.array([k["ti"]["lp"][2] if "lp" in k["ti"] else 0.0 for k in K])
        tdev = np.array([abs(k["ti"].get("t", 0.5) - 0.5) for k in K])
        R = g("R") / st.bond
        cost = CS + 1.5 * CH + wshape * ((g("sw") - pref) / 40.0) ** 2 + 0.25 * g("L") / st.bond + 0.6 * (g("aim") / 12.0) ** 2 \
            + 0.6 * (g("dep") / 15.0) ** 2 + 0.5 * (lpd / 30.0) ** 2 + 2.0 * tdev + ring \
            + 3.0 * np.maximum(0, 0.22 - R) / 0.07 + 2.0 * (1.0 - g("even"))     # a bend that tightens part-way reads as a kink
        return V, cost
    @staticmethod
    def _crossings(SH, HS):
        """For shaft polylines SH (N, n, 2): how many of the segments HS (M, 4) each one crosses (proper crossings)."""
        if not len(HS) or not len(SH): return np.zeros(len(SH))
        P0, P1 = SH[:, :-1, None, :], SH[:, 1:, None, :]; Q0, Q1 = HS[None, None, :, :2], HS[None, None, :, 2:]
        r, s = P1 - P0, Q1 - Q0; den = r[..., 0] * s[..., 1] - r[..., 1] * s[..., 0]; w = Q0 - P0
        ok = np.abs(den) > 1e-12; den = np.where(ok, den, 1.0)
        t = (w[..., 0] * s[..., 1] - w[..., 1] * s[..., 0]) / den; u = (w[..., 0] * r[..., 1] - w[..., 1] * r[..., 0]) / den
        hit = ok & (t >= 0) & (t <= 1) & (u > 0.02) & (u < 0.98)
        return hit.any(1).sum(1).astype(float)
    def _arrow_segs(self, c):
        P = np.vstack([c._arc, c._P - c._d * 0.8 * self.st.c_head_len]); H = c._head
        return np.vstack([np.hstack([P[:-1], P[1:]]), np.hstack([H, np.roll(H, -1, axis=0)])])
    def _solve_one(self, c, O, pref=None, wshape=1.0, K=None, quick=False):
        cls = self._aclass(c); p = float(c.bend) if c.bend is not None else (pref if pref is not None else self.ARROW_PREF[cls])
        K = K if K is not None else (self._arrow_cands(c) or self._arrow_cands(c, relax=True))
        if not K:                                              # nothing fits at all: a plain arc, flagged by check()
            T = self._tails(c.src, self._pt(self._edst(c)))[0][0]; P = self._tips(self._edst(c), T)[0][0]; d = unit(P - T)
            K = [{"T": T, "P": P, "d": d, "pts": T + (P - d * self.st.c_head_len - T) * np.linspace(0, 1, 26)[:, None], "sw": 0.0,
                  "R": np.inf, "L": norm(P - T), "dep": 0.0, "aim": 0.0, "ccw": True, "ti": {}}]
        extra = lambda K: self._extra(c, cls, K)
        def best(K):
            if len(K) > 150 or quick:
                V, cost = self._arrow_eval(c, K, O, p, wshape, coarse=True); cost = cost + extra(K)
                K = [K[i] for i in np.lexsort((cost, V > 0.05))[:(1 if quick else 150)]]
            V, cost = self._arrow_eval(c, K, O, p, wshape); feas = V <= 0.05; pen = extra(K)
            i = int(np.argmin(np.where(feas, cost + pen, np.inf))) if feas.any() else int(np.argmin(4.0 * V + cost))
            return K[i], bool(feas.any()), float(V[i]), float(cost[i])
        k, ok, v, cs = best(K)
        if not ok:                                             # no arc gets past the ink: a freer, hand-drawn curve
            F = self._free_cands(c)
            if F:
                k2, ok2, v2, cs2 = best(F)
                if ok2 or 4.0 * v2 + cs2 < 4.0 * v + cs: k, ok, v = k2, ok2, v2
        self._finish(c, k); c._class = cls; c._ok = ok; c._V = v
    FREE_DEP = {"lp": (0, 8, -8, 15, -15),    "bond": (0, 25, -25, 50, -50), "atom": (0, 30, -30, 60, -60, 90, -90)}
    FREE_HANDLE = (0.3, 0.55, 0.85)   # handle lengths of a free curve, as fractions of its chord
    def _free_cands(self, c):
        """Candidate arrows beyond arcs and joined arcs, for an arrow none of those can draw without touching or
        crossing ink: smooth cubic curves from the same tails to the same tips, leaving the tail in a wider fan of
        directions and coming into the head with long or short handles, so the shaft can swing out round a bond or a
        label and come back in the way a hand-drawn arrow does. Only curves that bend one way (no S), turn between 25
        and 270 deg and nowhere more tightly than the arcs may are kept. Same fields as _arrow_cands' candidates."""
        st = self.st; h = st.c_head_len; cd = self._edst(c); dpt = self._pt(cd)
        tgt = (cd,) if isinstance(cd, Atom) else (cd.a, cd.b) if isinstance(cd, Mid) else tuple(cd) if isinstance(cd, (tuple, list)) else ()
        tails = self._tails(c.src, dpt, tgt)
        tips = self._tips(cd, tails[0][0] if tails else dpt, wide=isinstance(cd, Mid) and not self._concerted(c))
        force = None if c.side is None else (1 if c.side == "left" else -1)
        rk = lambda v: None if v is None else np.asarray(v, float).tobytes()
        key = (force, h, tuple((rk(T), rk(u0), sorted(ti)) for T, u0, ti in tails),
               tuple((rk(P), tuple((rk(d), a_) for d, a_ in hd)) for P, hd, pi in tips))
        if getattr(c, "_free", None) is not None and c._free[0] == key: return c._free[1]   # (the same inputs: reuse)
        rmin, lmin = 0.15 * st.bond, 0.22 * st.bond; n = 26
        t = np.linspace(0.0, 1.0, n)[:, None]; tf = np.linspace(0.0, 1.0, 90)[:, None]
        bz = lambda C0, C1, C2, C3, t: ((1 - t) ** 3) * C0 + 3 * ((1 - t) ** 2) * t * C1 + 3 * (1 - t) * t * t * C2 + (t ** 3) * C3
        out = []
        for T, u0, ti in tails:
            kind = "lp" if "lp" in ti else "bond" if "bond" in ti else "atom"
            dirs = [(rotv(u0, r), float(abs(r))) for r in self.FREE_DEP[kind]] if u0 is not None else \
                   [(dirv(a_), 0.0) for a_ in range(0, 360, 30)]
            for P, hd, pi in tips:
                for d, aim in hd:
                    B = P - d * h; L = norm(B - T)
                    if L < lmin: continue
                    for t0, dep in dirs:
                        for fa in self.FREE_HANDLE:
                            for fb in self.FREE_HANDLE:
                                C1, C2 = T + t0 * fa * L, B - d * fb * L
                                Q = bz(T, C1, C2, B, t); D = np.diff(Q, axis=0); sl = np.hypot(D[:, 0], D[:, 1])
                                if sl.min() < 1e-6: continue
                                tu = D / sl[:, None]; tu = np.vstack([t0[None], tu, d[None]])
                                cr = tu[:-1, 0] * tu[1:, 1] - tu[:-1, 1] * tu[1:, 0]; dt = tu[:-1, 0] * tu[1:, 0] + tu[:-1, 1] * tu[1:, 1]
                                turn = np.arctan2(cr, dt)
                                if (turn > 0.01).any() and (turn < -0.01).any(): continue          # an S: never
                                sw = float(np.degrees(abs(turn.sum())))
                                if not 25.0 <= sw <= 270.0: continue
                                seg = np.concatenate([[sl[0]], (sl[:-1] + sl[1:]) / 2, [sl[-1]]])
                                curv = np.abs(turn) / np.maximum(seg, 1e-9)
                                R = float(1.0 / curv.max()) if curv.max() > 1e-12 else math.inf
                                if R < rmin: continue
                                Rx = float(1.0 / max(curv.min(), 1e-12))
                                ccw = bool(turn.sum() > 0)
                                out.append({"T": T, "P": P, "d": d, "sw": sw, "R": R, "L": float(sl.sum()) + h, "dep": dep,
                                            "aim": aim, "ccw": ccw, "ti": ti, "bi": False, "even": min(1.0, R / Rx) if math.isfinite(Rx) else 0.0,
                                            "pts": Q, "fine": bz(T, C1, C2, B, tf)})
        out = self._side_filter(out, force)
        c._free = (key, out)
        return out
    def _extra(self, c, cls, K):
        """Shape rules on top of the ranking: arrows that start and end close together (adj, lpadj, hook) want a
        shaft of 0.8 BU (two head lengths), so they never read as all head; the shortfall is weighed against the other
        shape costs rather than forced, so the arrow never loops round a neighbour to gain length (a hook onto an
        unlabelled atom, a radical carbon say, aims for 0.5 BU: its natural curl is tight). The arrows of a cyclic
        array (a pericyclic step) all sit on one side of the ring and turn the same way round it."""
        pen = np.zeros(len(K)); h = self.st.c_head_len
        tgt = self._shaft_target(c, cls)
        lp = getattr(self, "_loopinfo", {}).get(id(c))
        if tgt:
            d = np.maximum(0.0, tgt - (np.array([k["L"] for k in K]) - h) / self.st.bond)
            pen += 120.0 * d ** 2 + (0.0 if lp is not None else self.SHAFT_STEP) * np.minimum(1.0, d / 0.005)
        if lp is not None: pen += self._loop_pen(K, lp, getattr(self, "_loopmode", {}).get(id(c), "in"))
        return pen
    def _frame(self, c):
        """An arrow's own frame for comparing shapes: (origin, x, y) in pt, the origin at the middle of its source bond,
        x along that bond towards the atom it shares with its target, y square to it; None if it has no such atom."""
        s, d = c.src, self._edst(c)
        if not (isinstance(s, (tuple, list)) and len(s) == 2 and all(isinstance(x, Atom) for x in s)): return None
        b = s[0].mol.find(*s)
        if b is None: return None
        sh = [x for x in s if any(x is y for y in self._xatoms(d))]
        if len(sh) != 1: return None
        far = s[1] if sh[0] is s[0] else s[0]; ex = unit(self.P(sh[0]) - self.P(far))
        return np.asarray(b._mid, float), ex, perp(ex)
    def _alike(self, grp, cands, O):
        """The arrows of a cyclic array that are rotated copies of each other (the three of benzene's Kekulé resonance)
        take one shape: of the shapes (in each arrow's own frame) that every one of them can take cleanly, the one with
        the least total cost. Left as placed if they are not alike, no shape suits all, or the copies would touch."""
        cls = {self._aclass(c) for c in grp}; fr = [self._frame(c) for c in grp]
        if len(cls) != 1 or any(f is None for f in fr): return
        cls = cls.pop(); O2 = dict(O); O2["AS"] = np.zeros((0, 4)); feats, costs = [], []
        for c, (o, ex, ey) in zip(grp, fr):
            K = cands[id(c)]
            if not K: return
            V, cost = self._arrow_eval(c, K, O2, self.ARROW_PREF[cls]); cost = np.where(V <= 0.05, cost + self._extra(c, cls, K), np.inf)
            F = np.array([[(k["T"] - o) @ ex, (k["T"] - o) @ ey, (k["P"] - o) @ ex, (k["P"] - o) @ ey, 10 * (k["d"] @ ex),
                           10 * (k["d"] @ ey), k["sw"] / 10, 100.0 * bool(k["ccw"]), 100.0 * bool(k["bi"])] for k in K])
            feats.append(F); costs.append(cost)
        best = None
        for i in np.argsort(costs[0], kind="stable"):
            if not np.isfinite(costs[0][i]): break
            pick, tot = [i], costs[0][i]
            for F, cost in zip(feats[1:], costs[1:]):
                dd = np.abs(F - feats[0][i]).max(1); j = int(np.argmin(dd))
                if dd[j] > 0.05 or not np.isfinite(cost[j]): pick = None; break
                pick.append(j); tot += cost[j]
            if pick is not None and (best is None or tot < best[0] - 1e-9): best = (tot, pick)
        if best is None: return
        keys = ("_pts", "_head", "_ctrl", "_arc", "_T", "_P", "_d", "_sweep", "_R", "_lpc", "_class", "_ok", "_V", "_barb", "_ccw")
        keep = [{k: getattr(c, k) for k in keys} for c in grp]
        for c, j in zip(grp, best[1]): self._finish(c, cands[id(c)][j]); c._ok = True; c._V = 0.0
        segs = [self._arrow_segs(c) for c in grp]
        if any((_segs_dist(np.vstack([S1[:, :2], S1[:, 2:]]), S2) < self.st.alw + 1.0).any()
               for i_, S1 in enumerate(segs) for S2 in segs[i_ + 1:]):
            for c, kv in zip(grp, keep):
                for k, v in kv.items(): setattr(c, k, v)
    def _joint(self, grp, cands, O):
        """A cyclic array with an arrow left touching something (one placed before the other took its room: the two of
        a 4-centre array, each wanting the inside) is placed jointly: of each arrow's best clean candidates, the
        combination with the least total cost whose arrows keep clear of each other."""
        import itertools as _it
        O2 = dict(O); O2["AS"] = np.zeros((0, 4)); top = []; n = 120 if len(grp) == 2 else 30
        for c in grp:
            K = cands[id(c)]
            if not K: return
            cls = self._aclass(c); V, cost = self._arrow_eval(c, K, O2, self.ARROW_PREF[cls]); cost = cost + self._extra(c, cls, K)
            ok = [i for i in np.argsort(cost, kind="stable") if V[i] <= 0.05][:n]
            if not ok: return
            P = [np.vstack([K[i]["pts"], self._head_poly(K[i]["P"], K[i]["d"], c.kind)]) for i in ok]
            top.append((np.array([cost[i] for i in ok]), ok, P))
        lim = self.st.alw + 1.0; clear = {}
        for a, b in _it.combinations(range(len(grp)), 2):           # clearance between every candidate of a and of b
            Sb = np.vstack([np.hstack([q[:-1], q[1:]]) for q in top[b][2]]); nb = len(top[b][2][0]) - 1
            clear[a, b] = np.array([_segs_dist(q, Sb).min(0).reshape(-1, nb).min(1) for q in top[a][2]]) >= lim
        best = None
        for combo in _it.product(*[range(len(t[1])) for t in top]):
            tot = sum(float(top[k][0][i]) for k, i in enumerate(combo))
            if (best is None or tot < best[0]) and all(clear[a, b][combo[a], combo[b]] for a, b in clear): best = (tot, combo)
        if best is None: return
        for k, (c, i) in enumerate(zip(grp, best[1])): self._finish(c, cands[id(c)][top[k][1][i]]); c._ok = True; c._V = 0.0
    SHAFT_STEP = 1.0        # extra cost of falling short of the shaft target at all: a candidate that reaches it wins
                            # unless it costs this much more in shape; the shortfall itself is then weighed (120 d^2).
                            # Not for the arrows of a cyclic array, which share the ring and are kept alike instead.
    def _shaft_target(self, c, cls=None):
        """Shaft length (BU) a short-range arrow aims for: 0.8 BU (two head lengths) for adj, lpadj and hook arrows, 0.5
        BU for a hook onto an unlabelled atom; None for intermolecular arrows (their length is set by the docking)."""
        cls = cls or self._aclass(c)
        if cls not in ("adj", "lpadj", "hook", "shift"): return None
        if cls == "hook":
            d = self._edst(c)
            if isinstance(d, Atom) and not d.drawn: return 0.5
        if cls == "adj" and isinstance(c.dst, (tuple, list)) and id(c) in getattr(self, "_loopinfo", {}):
            b1 = c.src[0].mol.find(*c.src); b2 = c.dst[0].mol.find(*c.dst)    # a cyclic array round a small ring (benzene's
            r1 = c.src[0].mol.ring_of(b1, 6) if b1 is not None else None       # Kekule arrows): adjacent bond middles are
            if r1 and b2 is not None and all(any(x is y for y in r1) for x in (b2.a, b2.b)): return 0.55   # 0.87 BU apart
        return 0.8
    @staticmethod
    def _loop_pen(K, lp, mode):
        """Penalty for an arrow of a cyclic array: inside the ring an arrow that leaves its bond square and points into
        the next bond turns against the flow round the ring; outside, with it (judged on the side it is on). The arrows of
        one array keep to one side (mode) where they can: a 4-centre array has room inside for one only."""
        from matplotlib.path import Path as _Path
        C, poly, flow = lp; mids = np.array([k["pts"][len(k["pts"]) // 2] for k in K]); ins = _Path(poly).contains_points(mids)
        want = np.where(ins, not flow, bool(flow))
        return 80.0 * (np.array([bool(k["ccw"]) for k in K]) != want).astype(float) + 140.0 * (ins != (mode == "in")).astype(float)
    @staticmethod
    def _xatoms(x):
        return [x.atom] if isinstance(x, LP) else [x] if isinstance(x, Atom) else [x.a, x.b] if isinstance(x, Mid) \
            else [y for y in x if isinstance(y, Atom)] if isinstance(x, (tuple, list)) else []
    def _loops(self, curlies):
        """Cyclic arrays of arrows (pericyclic steps): each arrow's destination shares an atom with the next arrow's
        source, closing a loop. Returns {id(arrow): (centre, ring polygon in pt, True if the flow runs anticlockwise)}."""
        cs = [c for c in curlies if not self._concerted(c)]    # (two arrows meeting nose to nose at a forming bond, a
        nxt = {}                                                # radical pair, are not a flow round a ring)
        for i, c in enumerate(cs):
            d = {id(a) for a in self._xatoms(c.dst)}
            nxt[i] = [j for j, e in enumerate(cs) if j != i and d & {id(a) for a in self._xatoms(e.src)}]
        out = {}
        for i0 in range(len(cs)):
            path, i = [i0], i0; seen = {i0}
            while nxt[i]:
                j = nxt[i][0]
                if j == i0: break
                if j in seen: path = None; break
                path.append(j); seen.add(j); i = j
            if path is None or not nxt[i] or nxt[i][0] != i0 or len(path) < 2 or any(id(cs[k]) in out for k in path): continue
            ats = {}
            for k in path:
                for a in self._xatoms(cs[k].src) + self._xatoms(cs[k].dst): ats[id(a)] = a
            P = np.array([self.P(a) for a in ats.values()])
            if len(P) < 3: continue
            C = P.mean(0); H = _hull(P)
            if len(H) < 3: continue
            # the sense of the flow: each arrow's turn about the centre from its tail to its head (the area of the tails'
            # polygon is zero for two arrows, a 4-centre array such as a beta-hydride elimination)
            flow = sum(cross(self._pt(cs[k].src) - C, self._pt(self._edst(cs[k])) - C) for k in path)
            for k in path: out[id(cs[k])] = (C, H, flow > 0)
        return out
    def _geo_sig(self, mols):
        ref = mols[0].atoms[0].xy.copy() if mols and mols[0].atoms else np.zeros(2); s = []
        for m in mols:
            for a in m.atoms:
                s += list(a.xy - ref)
                for mk in (a._chg, a._rad, a._pc, a._num, a._site):
                    if mk: s += list(mk[0] / self.st.bond - ref)
                s += [l.angle if l.angle is not None else -999.0 for l in a.lps]
        return tuple(np.round(s, 6)), ref
    def _solve(self, curlies, mols, skip=frozenset(), quick=False, second=True):
        """Place the curly arrows: one at a time, then once more each with all the others in place (so they avoid
        each other both ways) and with each family pulled towards its common shape. Results are cached by geometry;
        quick=True (layout passes) keeps the best coarse candidate."""
        if not curlies: return
        key = (tuple(id(c) for c in curlies), tuple(id(m) for m in mols), skip, quick, second)
        sig, ref = self._geo_sig(mols); cache = self.__dict__.setdefault("_acache", {}); hit = cache.get(key)
        if hit is not None and hit[0] == sig:
            dv = (ref - hit[1]) * self.st.bond
            for c, r in zip(curlies, hit[2]):
                for k, v in r.items(): setattr(c, k, v + dv if k in ("_pts", "_head", "_ctrl", "_arc", "_T", "_P") else v)
            return
        self._loopinfo = self._loops(curlies); self._loopmode = {}
        O = self._obst(mols, skip); placed = {}; cands = {}
        for c in curlies:
            cands[id(c)] = self._arrow_cands(c) or self._arrow_cands(c, relax=True)
        groups = {}
        for c in curlies:
            lp = self._loopinfo.get(id(c))
            if lp is not None: groups.setdefault((round(lp[0][0], 3), round(lp[0][1], 3)), []).append(c)
        for grp in groups.values():                      # inside or outside the ring: whichever fits the whole array best
            score = {}
            for mode in ("in", "out"):
                tot = 0.0
                for c in grp:
                    K = cands[id(c)]
                    if not K: continue
                    V, cost = self._arrow_eval(c, K, O, self.ARROW_PREF[self._aclass(c)], coarse=True)
                    cost = cost + self._loop_pen(K, self._loopinfo[id(c)], mode); feas = V <= 0.05
                    tot += float(np.min(np.where(feas, cost, np.inf))) if feas.any() else 1e3 + float(np.min(4.0 * V + cost))
                score[mode] = tot
            for c in grp: self._loopmode[id(c)] = min(score, key=lambda md: (score[md], md != "in"))
        for c in curlies:
            O["AS"] = np.vstack([placed[id(x)] for x in curlies if id(x) in placed] or [np.zeros((0, 4))])
            self._solve_one(c, O, K=cands[id(c)], quick=quick); placed[id(c)] = self._arrow_segs(c)
        med = {}
        curlies_b = curlies if (second and not quick) else []
        for c in curlies:
            sws = [x._sweep for x in curlies if x._class == c._class and x._ok]
            if len(sws) >= 2: med[c._class] = float(np.median(sws))
        for c in curlies_b:
            O["AS"] = np.vstack([placed[id(x)] for x in curlies if x is not c] or [np.zeros((0, 4))])
            self._solve_one(c, O, med.get(c._class), 2.0 if c._class in med else 1.0, K=cands[id(c)]); placed[id(c)] = self._arrow_segs(c)
        if second and not quick:
            for grp in groups.values():
                if len(grp) >= 3: O["AS"] = np.zeros((0, 4)); self._alike(grp, cands, O)
                if len(grp) <= 3 and any(not c._ok for c in grp): O["AS"] = np.zeros((0, 4)); self._joint(grp, cands, O)
        keys = ("_pts", "_head", "_ctrl", "_arc", "_T", "_P", "_d", "_sweep", "_R", "_lpc", "_class", "_ok", "_V", "_barb", "_ccw")
        cache[key] = (sig, ref, [{k: getattr(c, k) for k in keys} for c in curlies])
    def _finish(self, c, k):
        st = self.st; h = st.c_head_len; T, P, d = k["T"], k["P"], k["d"]; B = P - d * h
        pts = self._biarc_pts(T, B, k["par"], 90) if k.get("bi") else (self._arc_pts(T, B, k["par"], 90) if "par" in k else k.get("fine", k["pts"]))
        sb = -1 if k["ccw"] else 1                               # fishhook barb on the outside of the curve
        c._head = self._head_poly(P, d, c.kind, sb)
        c._pts = np.vstack([pts, P - d * (0.8 if c.kind == "full" else 0.65) * h])          # the shaft ends inside the head
        c._arc, c._T, c._P, c._d, c._barb = pts, T, P, d, sb
        c._sweep, c._R = k["sw"], k["R"]; c._lpc = k["ti"].get("lp"); c._ctrl = np.vstack([c._pts, c._head]); c._ccw = bool(k["ccw"])
    def _fragment_issues(self):
        """Species within one Mol: two fragments closer than 1.5 BU with no bond or curly arrow between them read as one
        molecule drawn wrongly (a leaving group left where it was attached); a fragment with no arrow lying more than
        3 BU from the rest of its Mol is a stray. Both are cured by split() (its own item behind '+') or by an arrow."""
        out = []
        cat = []
        for c in self.curlies:
            u = set()
            for x in (c.src, c.dst):
                for y in ([x.atom] if isinstance(x, LP) else [x] if isinstance(x, Atom) else [x.a, x.b] if isinstance(x, Mid)
                          else list(x) if isinstance(x, (tuple, list)) else []):
                    if isinstance(y, Atom): u.add(id(y))         # the atom objects themselves: copies share uids
            cat.append(u)
        for m in self.mols:
            frs, seen = [], set()
            for a in m.atoms:
                if a.uid in seen: continue
                fr = m.frag(a); seen.update(x.uid for x in fr)
                if any(not x.hidden for x in fr): frs.append(fr)
            if len(frs) < 2: continue
            ids = [{id(x) for x in f} for f in frs]; where = {u: i for i, f in enumerate(ids) for u in f}
            grp = list(range(len(frs)))                     # fragments drawn as one species: eta bonds, Newman front/back
            def root(i):
                while grp[i] != i: i = grp[i]
                return i
            for f in frs:
                for x in f:
                    for y in list(getattr(x, "_eta", None) or []) + list(getattr(x, "_float", None) or []) + \
                             ([x._newman_pair] if getattr(x, "_newman_pair", None) is not None else []):
                        if id(y) in where: grp[root(where[id(x)])] = root(where[id(y)])
            P = [np.array([x.xy for x in f if not x.hidden]) for f in frs]
            name = lambda f: repr(next((x for x in f if x.drawn), f[0]))
            V = [[x for x in f if not x.hidden] for f in frs]
            for i, j in itertools.combinations(range(len(frs)), 2):
                if root(i) == root(j): continue
                D = np.hypot(*(P[i][:, None, :] - P[j][None, :, :]).transpose(2, 0, 1)); k = np.unravel_index(np.argmin(D), D.shape)
                d = float(D[k])
                if d < 1.5 and not any(u & ids[i] and u & ids[j] for u in cat):
                    out.append(f"{V[i][k[0]]!r} sits {d:.2f} BU from {V[j][k[1]]!r} but no bond or curly arrow joins them "
                               "(split() a departing group off as its own item, or keep it with its arrow)")
            for i, f in enumerate(frs):
                if any(u & ids[k] for u in cat for k in range(len(frs)) if root(k) == root(i)): continue
                others = [j for j in range(len(frs)) if root(j) != root(i)]
                if not others: continue
                d = min(float(np.min(np.hypot(*(P[i][:, None, :] - P[j][None, :, :]).transpose(2, 0, 1)))) for j in others)
                if d > 3.0: out.append(f"stray species: {name(f)} lies {d:.1f} BU from the rest of its molecule with no arrow (split() it off as its own item)")
        return out
    def _head_issues(self):
        """Arrowheads of black arrows (reaction arrows, cycle arcs, branches) and full curly arrows: no part of the shaft
        (with its line width and cap) may show beyond the head's notch outside the head, and on a curved arrow the head's
        axis must follow the curve (within 5 deg of the tangent at the notch)."""
        from matplotlib.path import Path as _Path
        st = self.st; out = []; log = list(getattr(self, "_arrow_log", []))
        for c in self.curlies:
            if c.kind == "full": log.append({"shaft": np.asarray(c._pts), "poly": np.asarray(c._head), "lw": st.alw, "curved": False, "name": "curly"})
        for k, A in enumerate(log):
            H = A.get("poly")
            if H is None:
                tip, d, h, w, kk = A["tip"], A["d"], A["h"], A["w"], A["k"]; n = perp(d)
                H = np.array([tip, tip - d * h + n * w, tip - d * kk * h, tip - d * h - n * w])
            H = np.asarray(H, float); tip = H[0]; notch = H[2] if len(H) == 4 else H[2]
            d = unit(tip - notch) if len(H) == 4 else unit(tip - (H[1] + H[2]) / 2)
            if len(H) != 4: continue                               # half heads (equilibrium, fishhook): the shaft is their spine
            S = A["shaft"]; seg = np.diff(S, axis=0); L = np.hypot(*seg.T)
            Q = np.vstack([S[i] + seg[i] * t for i in range(len(seg)) for t in np.linspace(0, 1, max(2, int(L[i] / 0.25) + 1))] + [S[-1:]])
            hl = norm(tip - notch) / 0.875                         # head length
            back = np.concatenate([np.cumsum(np.hypot(*np.diff(Q, axis=0).T)[::-1])[::-1], [0.0]])
            Q = Q[back <= 0.6 * hl]                                # only the end of the shaft, where it meets its head
            r = A["lw"] / 2; tq = np.gradient(Q, axis=0) if len(Q) > 1 else np.array([d]); tq = tq / np.maximum(np.hypot(*tq.T), 1e-12)[:, None]
            nq = np.column_stack([-tq[:, 1], tq[:, 0]])
            edge = np.vstack([Q + nq * r, Q - nq * r, Q + tq * r])    # the stroke's outline: both sides and the round cap
            ahead = (edge - notch) @ d > 0.02
            if ahead.any() and not _Path(H).contains_points(edge[ahead], radius=-0.04).all():
                out.append(f"arrowhead at ({tip[0]:.0f}, {tip[1]:.0f}) pt: the shaft shows beyond the head")
            if A.get("curved") and len(S) >= 3:
                i0 = int(np.argmin(np.hypot(*(S[:-1] - notch).T)))      # where the curve meets the head
                if i0 >= 1:
                    dev = abs(adiff(ang(d), ang(S[i0] - S[i0 - 1])))
                    back = np.concatenate([[0.0], np.cumsum(np.hypot(*np.diff(S[:i0 + 1][::-1], axis=0).T))])
                    j = i0 - int(np.searchsorted(back, hl))                # one head length further back along the curve
                    turn = abs(adiff(ang(S[i0] - S[i0 - 1]), ang(S[j + 1] - S[j]))) if 0 <= j < i0 - 1 else 0.0
                    # on a curve the head lies along its chord: half the curve's turn across the head, plus 2 deg
                    if dev > 2 + 0.6 * turn: out.append(f"arrowhead at ({tip[0]:.0f}, {tip[1]:.0f}) pt sits {dev:.0f} deg off its curve")
        return out
    def _arrow_issues(self):
        """Arrow QA for check(): heads that miss their atom or touch anything, shafts that touch or cross, loose tails."""
        st = self.st; out = []; O0 = self._obst(self.mols); lw2 = st.lw / 2; self._arrow_notes = []
        segs = [self._arrow_segs(c) for c in self.curlies]
        loops = {}
        for c in self.curlies:
            lp = self._loops(self.curlies).get(id(c))
            if lp is not None: loops.setdefault((round(lp[0][0], 3), round(lp[0][1], 3)), []).append(c)
        from matplotlib.path import Path as _Path
        for grp in loops.values():                   # each turns with the flow round the ring, as seen from its own side
            C, poly, flow = self._loops(self.curlies)[id(grp[0])]
            bad = [c for c in grp if bool(getattr(c, "_ccw", True)) != ((not flow) if _Path(poly).contains_point(c._arc[len(c._arc) // 2]) else bool(flow))]
            if bad:
                out.append("curly arrows " + ", ".join(f"#{self.curlies.index(c)}" for c in grp) + " form a cyclic array but do not all turn the same way round it")
        byid = {id(b): b for m in self.mols for b in m.bonds}; aid = {id(a): a for m in self.mols for a in m.atoms}
        lpid = {id(l): l.atom for m in self.mols for a in m.atoms for l in a.lps}
        for i, c in enumerate(self.curlies):
            nm = f"curly #{i} ({c.src!r} -> {c.dst!r})"; sb, sa, slp, sab = self._arrow_rel(c)
            A = c._arc; s = np.concatenate([[0], np.cumsum(np.hypot(*np.diff(A, axis=0).T))]); body = A[s >= 2.5]
            H = self._head_samples(c._P, c._d, c.kind, c._barb); probs = []
            allp = np.vstack([A, H]); O = self._osub(O0, allp.min(0) - 6.0, allp.max(0) + 6.0)
            if len(O["S"]):
                for Q, lim, where in ((H, 0.4, "head"), (body, 0.5, "shaft")):
                    if not len(Q): continue
                    D = _segs_dist(Q, O["S"]) - lw2 - np.where(O["Sid"] == sb, 0.25 if where == "shaft" else lim, lim)[None, :]
                    j = np.unravel_index(np.argmin(D), D.shape)[1]
                    if D.min() < 0: b = byid[int(O["Sid"][j])]; probs.append(f"{where} touches bond {b.a}-{b.b}")
            if len(O["Gid"]):
                for Q, lim, where in ((H, 0.25, "head"), (body, 0.4, "shaft")):
                    if not len(Q): continue
                    D = _polys_dist(Q, O["G"], O["Gst"]) - lim
                    if (D < 0).any(): j = int(np.argmin(D.min(0))); probs.append(f"{where} touches the label of {aid[int(O['Gid'][j])]}")
            if len(O["C"]):
                ok = O["Cid"] != slp
                for Q, lim, where in ((H, 0.4, "head"), (body, 0.4, "shaft")):
                    if not len(Q) or not ok.any(): continue
                    D = np.hypot(Q[:, None, 0] - O["C"][None, ok, 0], Q[:, None, 1] - O["C"][None, ok, 1]) - O["C"][None, ok, 2] - lim
                    if (D < 0).any():
                        j = int(O["Cid"][ok][np.argmin(D.min(0))]); probs.append(f"{where} touches a charge or electrons on {aid.get(j) or lpid.get(j)}")
            for j, S2 in enumerate(segs):
                if j == i: continue
                Q = np.vstack([body, H]) if len(body) else H
                if (_segs_dist(Q, S2) - lw2 < 0.6).any(): probs.append(f"touches curly #{j}")
            # clear white space: the shaft and head keep at least 1 pt (0.8 tolerated) from every bond and label they do not
            # start or end at; less than that reads as touching at normal size
            if not any("touches" in p for p in probs):
                k = st.bond / 14.4; aw2 = st.alw / 2; cdt = self._edst(c); ta = cdt if isinstance(cdt, Atom) else None
                tbd = cdt[0].mol.find(*cdt) if isinstance(cdt, (tuple, list)) and len(cdt) == 2 and isinstance(cdt[0], Atom) else None
                near_t = {id(b) for b in ta.mol.bonds if ta in (b.a, b.b)} if ta is not None else ({id(tbd)} if tbd is not None else set())
                worst = None
                if len(O["S"]):
                    for Q, extra, skip in ((body, aw2, set()), (H, 0.0, near_t)):
                        if not len(Q): continue
                        D = _segs_dist(Q, O["S"]).min(0) - lw2 - extra
                        for g, sid in zip(D, O["Sid"]):
                            if int(sid) not in skip and g < 0.8 * k and (worst is None or g < worst[0]):
                                b = byid[int(sid)]; worst = (g, f"bond {b.a}-{b.b}")
                if len(O["Gid"]):
                    for Q, extra, skipa in ((body, aw2, -1), (H, 0.0, id(ta) if ta is not None else -1)):
                        if not len(Q): continue
                        D = _polys_dist(Q, O["G"], O["Gst"]).min(0) - extra
                        for g, gid in zip(D, O["Gid"]):
                            if int(gid) != skipa and g < 0.8 * k and (worst is None or g < worst[0]):
                                worst = (g, f"the label of {aid[int(gid)]}")
                if worst is not None: probs.append(f"passes only {max(worst[0], 0):.1f} pt from {worst[1]} (keep at least 1 pt of clear space)")
            cd = self._edst(c)
            if isinstance(cd, Atom):
                X = self.P(cd); err = abs(adiff(ang(c._d), ang(X - c._P)))
                if err > 26: probs.append(f"head points {err:.0f} deg away from {cd!r} (misses the atom)")
                if cd._g:
                    gap = min(_hull_dist(c._P, hh) for g in cd._g for hh in g["hulls"])
                    if gap > 2.6 * st.bond / 14.4: probs.append(f"tip stops {gap:.1f} pt short of {cd!r}")
                elif norm(c._P - X) > 0.32 * st.bond: probs.append(f"tip stops {norm(c._P - X) / st.bond:.2f} BU from {cd!r}")
            short_mid = isinstance(cd, Mid) and norm(cd.a.xy - cd.b.xy) <= 1.3
            # (an atom target has its own tests above: a tip at the edge of a wide label is far from the atom's centre)
            if not isinstance(cd, Atom) and not (isinstance(cd, Mid) and (short_mid or self._concerted(c))):
                da = min([norm(c._P - self.P(x)) for m in self.mols for x in m.atoms if not x.hidden] or [99.0])
                db = min([float(_segs_dist(c._P[None], np.array([[*b._vis[0], *b._vis[1]]]))[0, 0]) for m in self.mols for b in m.bonds
                          if getattr(b, "_vis", None) is not None] or [99.0])
                if min(da, db) > 0.45 * st.bond:
                    probs.append(f"ends {min(da, db) / st.bond:.2f} BU from the nearest atom or bond (points at nothing)")
            elif isinstance(c.dst, (tuple, list)) and c.dst[0].mol.find(*c.dst) is not None:
                b = c.dst[0].mol.find(*c.dst); p, q = b._vis; e = q - p; den = cross(c._d, e)
                hitb = abs(den) > 1e-9 and 0 <= cross(p - c._P, e) / den <= 4.0 and 0 <= cross(p - c._P, c._d) / den <= 1
                if not hitb: probs.append(f"head does not point at bond {b.a}-{b.b}")
            hsx = self._crossings(c._arc[None], O["HS"])[0] if len(O.get("HS", [])) else 0
            if hsx: probs.append("shaft crosses a bond where it runs into a label (it slips between the bond and the label)")
            if np.isfinite(c._R) and c._R < 0.15 * st.bond: probs.append(f"loop radius only {c._R:.1f} pt")
            if c._sweep < 20: probs.append(f"is nearly straight ({c._sweep:.0f} deg of turn): an electron-pushing arrow must read as curved")
            if c._class == "inter" and c._sweep > 150:
                self._arrow_notes.append(f"{nm}: takes the long way round ({c._sweep:.0f} deg of turn): turn or move the partner so the "
                                         "arrow can run direct (attack(), deprotonate() and protonate() place partners for this)")
            shaft = float(np.hypot(*np.diff(c._arc, axis=0).T).sum()) / st.bond; tg = self._shaft_target(c, c._class)
            if tg and shaft < tg - 0.08:
                self._arrow_notes.append(f"{nm}: shaft only {shaft:.2f} BU, little more than its head (no clean placement reaches {tg:.1f} BU here)")
            if isinstance(c.src, Atom) and c.src.label == "H" and not c.src.radical and (c.src.charge or 0) >= 0:
                probs.append(f"starts on {c.src!r}, which has no electrons to give (no lone pair, charge or radical): start "
                             "it at a bond to the H, or make it a hydride (charge -1)")
            t0 = unit(A[min(3, len(A) - 1)] - A[0])
            if isinstance(c.src, LP) and c.src.angle is not None and abs(adiff(ang(t0), c.src.angle)) > 20:
                probs.append(f"leaves its lone pair {abs(adiff(ang(t0), c.src.angle)):.0f} deg off the pair's axis (not from its middle)")
            if isinstance(c.src, (tuple, list)) and c.src[0].mol.find(*c.src) is not None:
                bb = c.src[0].mol.find(*c.src); dev = 90 - math.degrees(math.acos(min(1.0, abs(float(t0 @ bb._u)))))
                ln = self._curl(c, self._edst(c))       # an arrow curling back beside its bond may lean away from where it
                lim = 35                                # turns back to: a 1,2-shift with its group up to 85 deg, else 65
                if ln is not None and float(t0 @ (ln[2] - ln[0])) > 0: lim = 85 if ln[4] else 65
                if dev > lim: probs.append(f"leaves bond {bb.a}-{bb.b} {dev:.0f} deg off the perpendicular")
            out += [f"{nm}: {p}" for p in dict.fromkeys(probs)]
        return out

    # ---------- layout
    def _wrap(self, it):
        if isinstance(it, str): return self.text(it) if it != "+" else self.plus()
        if isinstance(it, Mol): self.add(it); return Group(it)
        if isinstance(it, Group): self.add(*it.mols)
        return it
    @staticmethod
    def _cmem(cp):
        """The structures a caption stands under, as a tuple."""
        return cp[0] if isinstance(cp[0], tuple) else (cp[0],)
    def _sbox(self, m, prep=True, boxes=None):
        """Box (pt) of structure m as its captions see it: its ink, reaching down (and up) over any brackets round it with
        their charge or count and down over its pointers, so a caption goes under them, centred on the structure. boxes:
        structure boxes already measured, by id."""
        bx = dict(boxes or {})
        if id(m) not in bx: bx[id(m)] = self._mol_box(m, curlies=False, prep=prep)
        out = list(bx[id(m)])
        for ex in self.extras:
            if ex[0] == "brackets" and any(x is m for x in ex[1]):
                b = self._bracket_geo(ex, bx)["box"]; out[1] = min(out[1], b[1]); out[3] = max(out[3], b[3])
        for p in self.pointers:                            # and its pointers' arrows and notes below it
            if p.mols() and all(x is m for x in p.mols()): out[1] = min(out[1], self._ptr_geo(p)["box"][1])
        return tuple(out)
    def _own_brackets(self, mols, boxes=None):
        """Boxes (pt) of the brackets that enclose only structures among mols (part of that row item's structure box)."""
        ids = {id(m) for m in mols}
        return [self._bracket_geo(ex, boxes)["box"] for ex in self.extras
                if ex[0] == "brackets" and ex[1] and all(id(x) in ids for x in ex[1])]
    def _caps_geo(self, caps, box_of, pts_of=None, align=False):
        """Captions laid out (pt): {id(caption): rich-text geometry}. Each is centred under the box of its structure(s)
        (box_of(m): ink and brackets) and below their curly arrows (pts_of(caption): point arrays). Captions of one
        structure stack in the order written, as the lines of one block; a caption under several structures goes below
        their own captions, a caption gap lower. align=True (as drawn) sets the first captions of a row's structures on
        one line, the lowest of their natural ones (raised 0.1 em, as row captions always have been); the captions
        stacked under them follow on."""
        st = self.st; B = st.bond; sz = st.tsize; capH = self.F.cap * sz; chem = getattr(st, "chemtext", True)
        up = [cp for cp in caps if len(cp) > 6 and cp[6] == "above"]
        if up:                                             # headings over structures: stacked upwards, the first written
            caps = [cp for cp in caps if cp not in up]; out = {}                    # nearest; a row's on one line
            tops = {}
            nat = lambda c: max(box_of(m)[3] for m in self._cmem(c)) + c[2] * B + self._ink_shift(c[1], "above", sz)
            for cp in up:
                bs = np.array([box_of(m) for m in self._cmem(cp)]); x0, x1 = bs[:, 0].min(), bs[:, 2].max()
                key = frozenset(id(m) for m in self._cmem(cp))
                yb = tops[key] if key in tops else max(nat(c) for c in up if c[3] == cp[3]) if (align and cp[3]) else nat(cp)
                g = self._rich_geo(cp[1], (x0 + x1) / 2, yb, sz, "center", "bottom"); out[id(cp)] = g
                tops[key] = g["box"][3] + 0.2 * sz
            out.update(self._caps_geo(caps, box_of, pts_of, align) if caps else {}); return out
        order = sorted(range(len(caps)), key=lambda k: (len(self._cmem(caps[k])), k))
        mem = [frozenset(id(m) for m in self._cmem(cp)) for cp in caps]
        lev, above = {}, {}
        for n, k in enumerate(order):
            above[k] = [j for j in order[:n] if mem[j] & mem[k]]
            lev[k] = 1 + max(lev[j] for j in above[k]) if above[k] else 0
        out = {}
        for L in sorted(set(lev.values())):
            at = {}
            for k in (k for k in order if lev[k] == L):
                cp = caps[k]; bs = np.array([box_of(m) for m in self._cmem(cp)])
                x0, y0, x1, y1 = (*bs[:, :2].min(0), *bs[:, 2:].max(0))
                for P in (pts_of(cp) if pts_of else []): y0 = min(y0, P[:, 1].min())
                t = y0 - cp[2] * B - self._ink_shift(cp[1], "below", sz)
                for j in above[k]:
                    b = out[id(caps[j])]["box"]
                    if mem[j] == mem[k]:                   # the next line of the same structure's block
                        last = _typo((self._lines(caps[j][1]) or [""])[-1], chem=chem)
                        first = _typo((self._lines(cp[1]) or [""])[0], chem=chem)
                        pitch = max(1.2 * sz, self._line_ink(last, sz)[1] + self._line_ink(first, sz)[0] + 0.2 * sz)
                        t = min(t, b[1] + 0.22 * sz - pitch + capH)
                    else: t = min(t, b[1] - cp[2] * B - self._ink_shift(cp[1], "below", sz))
                at[k] = ((x0 + x1) / 2, t)
            if align and L == 0:                           # one line for a row's first captions
                low = {}
                for k, (cx, t) in at.items():
                    key = caps[k][3] if caps[k][3] else id(caps[k]); low[key] = min(low.get(key, t + 0.1 * sz), t + 0.1 * sz)
                at = {k: (cx, low[caps[k][3] if caps[k][3] else id(caps[k])]) for k, (cx, t) in at.items()}
            for k, (cx, t) in at.items(): out[id(caps[k])] = self._rich_geo(caps[k][1], cx, t, sz, "center", "top")
        return out
    def _cap_geo(self, cp, curly_pts=None, box=None):
        """One caption alone, under its structure(s) (box: given instead of measured)."""
        bx = None if box is None else {id(m): box for m in self._cmem(cp)}
        return self._caps_geo([cp], lambda m: self._sbox(m, boxes=bx), lambda c: curly_pts or [])[id(cp)]
    def _ibox(self, it):
        st = self.st
        if isinstance(it, Group):
            bs = [self._mol_box(m) for m in it.mols]; ids = {id(m) for m in it.mols}
            own0 = [c for c in self.curlies if self._owners(c) and self._owners(c) <= {id(m) for m in it.mols}]
            if own0: self._solve(own0, it.mols, quick=True)
            sbx = {}
            def sbox(m):
                if id(m) not in sbx: sbx[id(m)] = self._sbox(m)
                return sbx[id(m)]
            caps = [cp for cp in self.captions if all(id(m) in ids for m in self._cmem(cp))]
            if caps:
                G = self._caps_geo(caps, sbox, lambda cp: [np.vstack([c._pts, c._head]) for c in own0
                                                            if {id(m) for m in self._cmem(cp)} & self._owners(c)])
                bs += [G[id(cp)]["box"] for cp in caps]
            bs += self._own_brackets(it.mols)
            own = [c for c in self.curlies if self._owners(c) and self._owners(c) <= {id(m) for m in it.mols}]
            if own:
                self._solve(own, it.mols, quick=True)
                for c in own:
                    pts = np.vstack([c._pts, c._head]); bs.append((*pts.min(0), *pts.max(0)))
            for p in self.pointers:
                if p.mols() and all(any(m is x for x in it.mols) for m in p.mols()): bs.append(self._ptr_geo(p)["box"])
            bs = np.array(bs); return (bs[:, 0].min(), bs[:, 1].min(), bs[:, 2].max(), bs[:, 3].max())
        g = self._item_geo(it); return g["box"]
    def _imove(self, it, dx, dy):
        if isinstance(it, Group):
            for m in it.mols: m.move(dx / self.st.bond, dy / self.st.bond)
        else: it.move(dx, dy)
    def _icy(self, it):
        if isinstance(it, RArrow) and not it.vertical: return it.y
        b = self._ibox(it); return (b[1] + b[3]) / 2
    def _lboxes(self, it):
        """Layout boxes of a row item (pt): (full box with curly arrows and captions; structure box: atoms, labels, marks
        and brackets; y of the line the item is centred on; box of its captions or None; parts: [(box, kind)], kind
        'ink', 'txt', for a pointer's arrows and note 'pink', 'ptxt', and 'zone' for the empty place of an arrow's
        conditions, kept free of notes). A Group is centred on its main fragment (the
        largest by atom count), so a docked reagent or a caption does not pull the backbone off the reaction-arrow line.
        Pointer notes are left out of the full box: row() keeps them clear only of what stands level with them."""
        st = self.st
        if isinstance(it, RArrow) and not it.vertical:
            g = self._item_geo(it); P = np.array([q for pr in g["prims"] if pr[0] in ("line", "poly") for q in np.asarray(pr[1], float).reshape(-1, 2)])
            parts = [((*P.min(0), *P.max(0)), "ink")] + [(t["box"], "txt") for t in g["texts"]] + [(self._reagent_ext(m)[1], "ink") for m in it._mols]
            for sg in (1, -1):                                 # where it carries nothing, the place of its conditions: a
                if not any((t["box"][1] + t["box"][3]) / 2 * sg > it.y * sg for t in g["texts"]):   # note there would
                    lo, hi = sorted((it.y, it.y + sg * 1.2 * st.bond))                                # read as one
                    parts.append(((P[:, 0].min(), lo, P[:, 0].max(), hi), "zone"))
            return g["box"], g["box"], it.y, None, parts
        if isinstance(it, Text):
            g = self._item_geo(it); ext = [self.F.path(t, x, y, sz, i, bd).get_extents() for t, x, y, sz, i, bd in g["runs"]
                                            if t.strip() and not self.F.missing(t)]
            cy = (min(e.y0 for e in ext) + max(e.y1 for e in ext)) / 2 if ext else (g["box"][1] + g["box"][3]) / 2
            return g["box"], g["box"], cy, None, [(g["box"], "txt")]
        if not isinstance(it, Group): b = self._ibox(it); return b, b, (b[1] + b[3]) / 2, None, [(b, "ink")]
        sb, best = [], None
        for m in it.mols:
            self._prep_mol(m); sb.append(self._mol_box(m, curlies=False, prep=False)); seen = set()
            for a in m.atoms:
                if a.uid in seen: continue
                fr = m.frag(a); seen.update(x.uid for x in fr)
                if best is None or len(fr) > len(best[1]): best = (m, fr)
            if best is not None and best[0] is m: mb = self._mol_box(m, curlies=False, atoms=best[1], prep=False)
        bx = {id(m): b for m, b in zip(it.mols, sb)}
        sb += self._own_brackets(it.mols, bx)                # brackets belong to the structure they enclose
        S = np.array(sb); S = (S[:, 0].min(), S[:, 1].min(), S[:, 2].max(), S[:, 3].max()); fb = [S]; parts = [(S, "ink")]
        own = [c for c in self.curlies if self._owners(c) and self._owners(c) <= {id(m) for m in it.mols}]
        if own:
            self._solve(own, it.mols, quick=True)
            for c in own: P = np.vstack([c._pts, c._head]); fb.append((*P.min(0), *P.max(0))); parts.append((fb[-1], "ink"))
        caps = [cp for cp in self.captions if all(id(m) in bx for m in self._cmem(cp))]; C = None
        if caps:
            G = self._caps_geo(caps, lambda m: self._sbox(m, boxes=bx), lambda cp: [np.vstack([c._pts, c._head]) for c in own
                                                                                  if {id(m) for m in self._cmem(cp)} & self._owners(c)])
            cb = np.array([G[id(cp)]["box"] for cp in caps]); fb += list(cb); parts += [(tuple(b), "txt") for b in cb]
            C = (cb[:, 0].min(), cb[:, 1].min(), cb[:, 2].max(), cb[:, 3].max())
        for p in self.pointers:                              # a pointer and its note belong to their molecule's item
            if p.mols() and all(any(m is x for x in it.mols) for m in p.mols()):
                g = self._ptr_geo(p); P = np.array([q for pr in g["prims"] for q in np.asarray(pr[1], float).reshape(-1, 2)])
                parts += [((*P.min(0), *P.max(0)), "pink")] + [(t["box"], "ptxt") for t in g["texts"]]
        B = np.array(fb); return (B[:, 0].min(), B[:, 1].min(), B[:, 2].max(), B[:, 3].max()), S, (mb[1] + mb[3]) / 2, C, parts
    def _next_dx(self, rs, F, S, C, parts, dy, x0, gap):
        """Horizontal move (pt) that puts the next item of a row in place, rs being the row so far ({'prev', 'capR',
        'placed'}, updated): structures `gap` BU apart; full boxes 0.4 BU; captions side by side CAP_SEP; and a pointer's
        arrows and note clear of whatever stands level with them (0.4 BU; a note next to a caption or note CAP_SEP), so
        a note above a reaction arrow may reach over it, as in teaching figures."""
        B = self.st.bond
        if rs["prev"] is None: dx = x0 - F[0]
        else:
            dx = max(rs["prev"][0] + gap * B - S[0], rs["prev"][1] + 0.4 * B - F[0])
            if C is not None and rs["capR"] is not None: dx = max(dx, rs["capR"] + self.CAP_SEP * B - C[0])
            for q, qk in parts:
                for p_, pk in rs["placed"]:
                    if qk[0] != "p" and pk[0] != "p": continue
                    if "zone" in (qk, pk) and "ptxt" not in (qk, pk): continue
                    tt = qk.endswith("txt") and pk.endswith("txt"); vp = 0.0 if tt else 0.2 * B   # (text boxes hold their
                    if q[1] + dy < p_[3] + vp and p_[1] < q[3] + dy + vp:                          # own leading)
                        dx = max(dx, p_[2] + (self.CAP_SEP if tt else 0.4) * B - q[0])
        rs["prev"] = (S[2] + dx, F[2] + dx)
        if C is not None: rs["capR"] = C[2] + dx if rs["capR"] is None else max(rs["capR"], C[2] + dx)
        rs["placed"] += [((q[0] + dx, q[1] + dy, q[2] + dx, q[3] + dy), qk) for q, qk in parts]
        return dx
    def _main_box(self, it):
        """Box (pt) of an item's main structure: for a Group, its largest fragment (a merged reagent left out)."""
        if not isinstance(it, Group): return self._ibox(it)
        best = None
        for m in it.mols:
            self._prep_mol(m); seen = set()
            for a in m.atoms:
                if a.uid in seen: continue
                fr = m.frag(a); seen.update(x.uid for x in fr)
                if best is None or len(fr) > len(best[1]): best = (m, fr)
        return self._mol_box(best[0], curlies=False, atoms=best[1], prep=False)
    def _span_box(self, items):
        """Box (pt) round row items, with any caption written under several of their structures."""
        bs = [self._ibox(i) for i in items]
        ms = {id(m) for i in items if isinstance(i, Group) for m in i.mols}
        joint = [cp for cp in self.captions if isinstance(cp[0], tuple) and all(id(m) in ms for m in cp[0])]
        if joint:                                         # (below the captions of their own structures)
            sbx = {}
            def sbox(m):
                if id(m) not in sbx: sbx[id(m)] = self._sbox(m)
                return sbx[id(m)]
            G = self._caps_geo([cp for cp in self.captions if all(id(m) in ms for m in self._cmem(cp))], sbox)
            bs += [G[id(cp)]["box"] for cp in joint]
        bs = np.array(bs); return (bs[:, 0].min(), bs[:, 1].min(), bs[:, 2].max(), bs[:, 3].max())
    def _space_joint(self, items, idx):
        """A caption under several structures of a row is centred under them, so it can reach out towards a neighbour's
        caption on the same line: move the items from the right-hand caption's first structure on until the two stand
        CAP_SEP apart. idx: id(structure) -> item index."""
        rid = self._rowid; B = self.st.bond
        if not any(isinstance(cp[0], tuple) and cp[3] == rid for cp in self.captions): return
        for _ in range(2 * len(items)):
            caps = [cp for cp in self.captions if cp[3] == rid]; sbx = {}
            def sbox(m):
                if id(m) not in sbx: sbx[id(m)] = self._sbox(m)
                return sbx[id(m)]
            G = self._caps_geo(caps, sbox, align=True); worst = None
            for c1, c2 in itertools.combinations(caps, 2):
                k1, k2 = [sorted(idx[id(m)] for m in self._cmem(c)) for c in (c1, c2)]
                if set(k1) & set(k2): continue                 # one stack: placed one under the other
                if k1[-1] < k2[0]: lft, rgt, k = c1, c2, k2[0]
                elif k2[-1] < k1[0]: lft, rgt, k = c2, c1, k1[0]
                else: continue                                 # interleaved: nothing a shift can mend
                b1, b2 = G[id(lft)]["box"], G[id(rgt)]["box"]
                if b1[1] > b2[3] + 0.1 * B or b2[1] > b1[3] + 0.1 * B: continue   # on different lines
                need = b1[2] + self.CAP_SEP * B - b2[0]
                if need > 1e-6 and (worst is None or k < worst[0]): worst = (k, need)
            if worst is None: return
            for it in items[worst[0]:]: self._imove(it, worst[1], 0.0)
    def row(self, items, x=0.0, y=0.0, gap=None, under=None, align="center", row_gap=1.3, via=None):
        """Lay items left to right, each centred on the line y (BU): a structure on its main fragment, an arrow on its
        shaft, '+' on its glyph. Structures sit `gap` BU apart (measured between atoms, labels, marks and brackets);
        curly arrows and other overhangs only push the next item along when they would come closer than 0.4 BU to it,
        and captions side by side keep CAP_SEP (1 BU) of white, so neighbouring captions never read as one line.
        under=the items of an earlier row (what row() returned): this row goes beneath that one, row_gap BU below its
        lowest ink (captions included), centred under it (align "left" or "right" lines up those edges instead); x and
        y are then not used. Two rows of products after one arrow: row([..., arrow, p1, "+", p2]), then
        row([p3, "+", p4], under=[p1, "+", p2]). via=S.arrow(vertical=True, ...) turns a scheme's corner: the arrow
        runs down from the last item of the row above (align="right", the default with via; "left": the first) to the
        matching item of this one, centred under it; a scheme that snakes back reads right to left with kind="<-":
        row([product, S.arrow(kind="<-"), intermediate], under=row1, via=S.arrow(vertical=True))."""
        if under is not None:
            items = self.row(items, 0.0, 0.0, gap); B = self.st.bond
            up = [Group(i) if isinstance(i, Mol) else i for i in under if not isinstance(i, str)]
            ub = self._span_box(up); nb = self._span_box(items)
            if via is not None:
                al = "left" if align == "left" else "right"
                pick = lambda its: next(i for i in (its if al == "left" else its[::-1]) if not isinstance(i, (Text, str)))
                top, bot = pick(up), pick(items); tb, bb = self._ibox(top), self._ibox(bot)
                L0 = self._item_geo(via)["box"]; L0 = L0[3] - L0[1]          # the vertical arrow's own height
                tm, bm = self._main_box(top), self._main_box(bot)              # centred on the main structures
                dx = (tm[0] + tm[2]) / 2 - (bm[0] + bm[2]) / 2
                dy = min(ub[1] - row_gap * B - nb[3], tb[1] - (L0 + 0.6 * B) - bb[3])
                for it in items: self._imove(it, dx, dy)
                bb = self._ibox(bot); cx = (tm[0] + tm[2]) / 2
                via.length = (tb[1] - bb[3] - 0.6 * B) / B; via.x, via.y = cx, tb[1] - 0.3 * B
                for _ in range(8):                         # what it carries keeps clear of both rows' other items: a
                    g = self._item_geo(via); d = 0.0       # longer arrow lowers its middle by half what it lowers this row
                    sb = [t["box"] for t in g["texts"]] + [self._reagent_ext(m)[1] for m in via._mols]
                    for its_, sg in ((up, 1), (items, -1)):
                        for i in its_:
                            if i is top or i is bot or isinstance(i, str): continue
                            ob = self._ibox(i); pd_ = 0.3 * B
                            for q in sb:
                                if q[0] < ob[2] + pd_ and ob[0] < q[2] + pd_ and q[1] < ob[3] + pd_ and ob[1] < q[3] + pd_:
                                    d = max(d, q[3] + pd_ - ob[1] if sg > 0 else ob[3] + pd_ - q[1])
                    if d <= 1e-6: break
                    for it in items: self._imove(it, 0.0, -2 * d)
                    via.length += 2 * d / B
                return items
            dx = {"center": ((ub[0] + ub[2]) - (nb[0] + nb[2])) / 2, "left": ub[0] - nb[0], "right": ub[2] - nb[2]}[align]
            for it in items: self._imove(it, dx, ub[1] - row_gap * B - nb[3])
            return items
        st = self.st; gap = st.gap if gap is None else gap; B = st.bond
        items = [self._wrap(i) for i in items]; self._rowid += 1; rs = {"prev": None, "capR": None, "placed": []}
        for it in items:
            F, S, cy, C, parts = self._lboxes(it)
            dx = self._next_dx(rs, F, S, C, parts, y * B - cy, x * B, gap); self._imove(it, dx, y * B - cy)
        ids = {id(m): k for k, it in enumerate(items) if isinstance(it, Group) for m in it.mols}
        for cp in self.captions:                           # the row's captions share its lines (a joint one too)
            if all(id(m) in ids for m in self._cmem(cp)): cp[3] = self._rowid
        self._space_joint(items, ids)
        ref = next((i for i in items if isinstance(i, RArrow) and not i.vertical), None) or next((i for i in items if isinstance(i, Text)), None)
        if ref is not None:                                # remembered so _build can re-centre on the marks as drawn
            self.__dict__.setdefault("_rows", []).append(
                (ref, [[g, g.mols[0].atoms[0].xy[1] * B - self._line_y(ref)] for g in items if isinstance(g, Group) and g.mols and g.mols[0].atoms]))
        return items
    def flow(self, items, width_cm=16.0, x=0.0, y=0.0, gap=None, row_gap=1.3, squeeze=0.7):
        """Wrap items into rows no wider than width_cm, breaking before a reaction arrow when possible (a '+' stays with
        its neighbours). Before a lone trailing item would wrap onto a row of its own, the last row's gaps shrink (down
        to `squeeze` of normal; squeeze=None never) to keep it. A continuation row that has to begin with an arrow is
        indented so that arrow lines up under the first arrow of the row above, where it fits."""
        st = self.st; gap = st.gap if gap is None else gap; W = width_cm / 2.54 * 72; X0 = x * st.bond
        items = [self._wrap(i) for i in items]; B = st.bond
        LB = {id(i): (self._lboxes(i), self._ibox(i)) for i in items}
        def width(r, g):                                   # what row() would give, worked out without moving anything
            rs = {"prev": None, "capR": None, "placed": []}; right = -1e18
            for it in r:
                (F, S, cy, C, parts), ib = LB[id(it)]
                dx = self._next_dx(rs, F, S, C, parts, -cy, X0, g); right = max(right, ib[2] + dx)
            return right - X0
        rows, cur = [], []
        for k, it in enumerate(items):
            if not cur or width(cur + [it], gap) <= W: cur.append(it); continue
            if k == len(items) - 1 and squeeze and width(cur + [it], gap * squeeze) <= W: cur.append(it); continue
            kk = len(cur)
            if not isinstance(it, RArrow):
                for j in range(len(cur) - 1, 0, -1):
                    if isinstance(cur[j], RArrow): kk = j; break
            rows.append(cur[:kk]); cur = cur[kk:] + [it]
        rows.append(cur)
        gaps = []
        for r in rows:                                     # the widest gap (up to normal) at which the row fits
            g = gap
            if width(r, gap) > W and squeeze:
                lo, hi = gap * squeeze, gap
                for _ in range(20):
                    mid = (lo + hi) / 2
                    if width(r, mid) <= W: lo = mid
                    else: hi = mid
                g = lo if width(r, lo) <= W else gap
            gaps.append(g)
        yy = y * st.bond; prev = None
        for r, g in zip(rows, gaps):
            self.row(r, x, 0.0, g); bs = [self._ibox(i) for i in r]
            if prev is not None and isinstance(r[0], RArrow):  # indent a continuation that begins with an arrow
                arrows = [self._ibox(i)[0] for i in prev if isinstance(i, RArrow)]
                if arrows:
                    dx = arrows[0] - bs[0][0]
                    if dx > 0 and max(b[2] for b in bs) + dx - X0 <= W:
                        for i in r: self._imove(i, dx, 0)
                        bs = [self._ibox(i) for i in r]
            top = max(b[3] for b in bs); dy = yy - top
            for i in r: self._imove(i, 0, dy)
            yy = min(self._ibox(i)[1] for i in r) - row_gap * st.bond; prev = r
        return rows
    def place(self, item, at):
        """Centre an item (Mol, Group, arrow, text) at `at` (BU)."""
        it = self._wrap(item); b = self._ibox(it)
        self._imove(it, at[0] * self.st.bond - (b[0] + b[2]) / 2, at[1] * self.st.bond - (b[1] + b[3]) / 2); return it
    def cycle(self, stations, center=(0, 0), radius=4.0, cw=True, gap=0.3):
        """Catalytic cycle: stations [(item, angle_deg)], each centred on the circle; arcs join them in order.
        The radius grows automatically until stations are >= 1 BU apart and every arc shows enough free length
        (2 BU, plus 1.6 BU for each reagent/by-product branch later added with branch())."""
        stations = [((self._wrap(it) if isinstance(it, str) else it), a) for it, a in stations]   # a text station is made once
        cyc = {"stations": stations, "c": np.array(center, float), "r": radius, "cw": cw, "gap": gap, "branches": []}
        self._cycle_fit(cyc); self.extras.append(("cycle", cyc)); return cyc
    def _cycle_fit(self, cyc):
        st = self.st; radius = cyc["r"]; center = cyc["c"]; cw = cyc["cw"]; gap = cyc["gap"]
        need = {}
        for br in cyc["branches"]: need[br[0]] = need.get(br[0], 0) + 1
        for _ in range(80):
            its = []
            for it, a in cyc["stations"]:
                self.place(it, (center[0] + radius * math.cos(math.radians(a)), center[1] + radius * math.sin(math.radians(a))))
                its.append((self._wrap(it) if not isinstance(it, Mol) else Group(it), a))
            bx = [self._ibox(i) for i, _ in its]; ok = True
            for b1, b2 in itertools.combinations(bx, 2):
                dx = max(b1[0] - b2[2], b2[0] - b1[2]); dy = max(b1[1] - b2[3], b2[1] - b1[3])
                if max(dx, dy) < 1.0 * st.bond: ok = False
            if ok:
                C = np.array(center, float) * st.bond; R = radius * st.bond
                for k in range(len(its)):
                    a0, a1 = its[k][1], its[(k + 1) % len(its)][1]; sg = -1 if cw else 1
                    span = (a1 - a0) * sg % 360 or 360
                    free = [t for t in np.arange(0, span, 1.0) if all(box_dist(C + dirv(a0 + sg * t) * R, b) > gap * st.bond for b in bx)]
                    if len(free) * math.radians(1.0) * R < (2.0 + 1.6 * need.get(k, 0)) * st.bond: ok = False
            if ok: break
            radius += 0.25
        cyc["its"] = its; cyc["r"] = radius
    def branch(self, cyc, i, label=None, kind="in", at=None, reach=1.4, side=None, note=None):
        """Reagent ('in') or by-product ('out') arrow on arc i of a cycle, with a label at its free end: text, or a
        structure (Mol) drawn there. at: where it joins the arc (0 = start, 1 = end); None puts a lone branch at the
        middle and, when an arc has both, the reagents before the by-products (0.3 and 0.7), so they never share a
        junction or cross. Adding a branch re-fits the cycle radius so the arc has room for it.
        kind='step': the name of arc i's step ("oxidative addition") beside the arc, at its middle or the nearest place
        clear of its branches, outside the cycle or else inside (side= and at= to choose). kind='feed': a straight arrow from outside into station i (a precatalyst entering
        the cycle), the label (text or Mol) at its tail, note= beside the arrow ("activation"), reach= its length (BU).
        On a reaction arrow, branch(arrow, label, kind="in" | "out", at=, side=) draws a co-reagent curving into it or a
        by-product curving out of it (the iminium salt leaving a beta-hydride elimination), the label (text or Mol) at
        its end, on the side the arrow's own text leaves free unless side= says."""
        if isinstance(cyc, RArrow):                         # branch(arrow, label, kind=): a co-reagent ('in') curving
            arr, lab = cyc, (i if label is None else label)  # into a reaction arrow, or a by-product ('out') curving
            if kind not in ("in", "out"):                    # out of it, at= along it (0.5; 0.3 and 0.7 for one of each),
                raise ValueError(f"branch() on a reaction arrow: kind {kind!r}: use 'in' or 'out'")   # side= "above" |
            if arr.vertical: raise ValueError("branch() on a reaction arrow: horizontal arrows only")  # "below" (default:
            if isinstance(lab, Mol): self.add(lab); arr._bmols.append(lab)                              # the side its own
            arr.branches.append((lab, kind, at, side, 1.2 if reach == 1.4 else reach)); return arr   # text leaves free)
        if label is None: raise TypeError("branch(cycle, i, label): give the arc's index and the label")
        if isinstance(label, Mol): self.add(label)
        if kind == "step": cyc.setdefault("steps", []).append((i, label, at, side)); return
        if kind == "feed": cyc.setdefault("feeds", []).append((i, label, reach, note)); return
        if kind not in ("in", "out"): raise ValueError(f"branch(): kind {kind!r}: use 'in', 'out', 'step' or 'feed'")
        cyc["branches"].append((i, label, kind, at, reach, side or "out")); self._cycle_fit(cyc)

    # ---------- item geometry
    @staticmethod
    def _lines(s):
        """Text as a list of lines: a list/tuple of strings, or a string with \n between lines (structures skipped).
        <i> and <b> left open at a line break carry on into the next line."""
        if s is None or isinstance(s, Mol): return []
        if isinstance(s, (list, tuple)): ls = [l for x in s if not isinstance(x, Mol) for l in str(x).split("\n")]
        else: ls = str(s).split("\n")
        out, open_ = [], []
        for l in ls:
            l = "".join(open_) + l; open_ = []
            for tg in re.findall(r"<(/?)(i|b|c)(=[^<>]*)?>", l):
                if not tg[0]: open_.append(f"<{tg[1]}{tg[2]}>")
                else:
                    k = max((j for j, o in enumerate(open_) if o[1] == tg[1]), default=None)
                    if k is not None: open_.pop(k)
            out.append(l + "".join(f"</{o[1]}>" for o in reversed(open_)))
        return out
    @staticmethod
    def _parts(s):
        """What an arrow carries on one side, top to bottom: ('t', [lines]) for a block of text, ('m', Mol) for a
        structure."""
        out = []
        for x in (s if isinstance(s, (list, tuple)) else [s]):
            if isinstance(x, Mol): out.append(("m", x)); continue
            if x is None or x == "": continue
            ls = str(x).split("\n")
            if out and out[-1][0] == "t": out[-1][1].extend(ls)
            else: out.append(("t", ls))
        return out
    def _open_arrow(self, tail, tip):
        """Retrosynthesis (open) arrow from tail to tip: two parallel shafts and an open head with a sharp mitred point.
        Each shaft is cut along the centre line of the head's arm it runs into, so its end lies wholly under that arm's
        stroke (no square end showing inside the head)."""
        st = self.st; d = unit(tip - tail); n = perp(d); o = 0.11 * st.bond; hl = 0.45 * st.bond; sp = 2.2 * o
        k = sp / hl; h = st.lw / 2; out = []
        for sg in (1, -1):
            t1, t2 = (o + h) / k, (o - h) / k                  # where the shaft's outer and inner edges meet the arm
            out.append(("poly", [tail + n * sg * (o - h), tip - d * t2 + n * sg * (o - h),
                                 tip - d * t1 + n * sg * (o + h), tail + n * sg * (o + h)], st.ink))
        out.append(("line", [tip - d * hl + n * sp, tip, tip - d * hl - n * sp], st.lw, st.ink, None, "butt"))   # sharp mitred tip
        return out
    def _ink_shift(self, lines, side, size):
        """How far (pt) a block of text placed next to an arrow (or under a structure) must move away from it so that
        its nearest ink keeps the spacing of plain capitals: superscripts on the first line under an arrow (H₃O⁺, ²⁺)
        rise above the cap height, subscripts on the last line over it dip below the usual descender allowance."""
        ls = [_typo(l, chem=getattr(self.st, "chemtext", True)) for l in self._lines(lines)]
        if not ls: return 0.0
        if side == "below":
            up = self._line_ink(ls[0], size)[0] - self.F.cap * size
            return up if up > 0.05 * size else 0.0
        dn = self._line_ink(ls[-1], size)[1] - 0.22 * size
        return dn if dn > 0.05 * size else 0.0
    def _reagent_ext(self, m):
        """(structure box, box with its captions) of a reagent structure where it stands (pt): its captions ("17") hang
        under it and travel with it in an arrow's stack."""
        b = self._sbox(m); f = b
        caps = [cp for cp in self.captions if all(x is m for x in self._cmem(cp))]
        if caps:
            G = self._caps_geo(caps, lambda x: b, align=True); cb = np.array([G[id(cp)]["box"] for cp in caps])
            f = (min(b[0], cb[:, 0].min()), min(b[1], cb[:, 1].min()), max(b[2], cb[:, 2].max()), max(b[3], cb[:, 3].max()))
        return b, f
    def _reagent_box(self, m, cx, edge, side, ha="center"):
        """Move reagent structure m so its structure box is centred on cx (or starts at cx, ha='left'; ends there,
        'right') with its bottom (side 'above') or top (side 'below') at edge; return the box (pt)."""
        b, f = self._reagent_ext(m)
        dx = cx - {"center": (b[0] + b[2]) / 2, "left": f[0], "right": f[2]}[ha]; dy = edge - (f[1] if side == "above" else f[3])
        if abs(dx) > 1e-9 or abs(dy) > 1e-9: m.move(dx / self.st.bond, dy / self.st.bond)
        return (f[0] + dx, f[1] + dy, f[2] + dx, f[3] + dy)
    def _stack_h(self, parts, side):
        """Height (pt) a stack of text blocks and reagent structures takes, as _stack_parts lays it out."""
        st = self.st; h = 0.0
        for k, (kd, v) in enumerate(parts):
            if kd == "t":
                g = self._rich_geo(v, 0.0, 0.0, st.tsize, "left", "top"); h += g["box"][3] - g["box"][1] + self._ink_shift(v, side, st.tsize)
            else:
                f = self._reagent_ext(v)[1]; h += f[3] - f[1] + (0.22 * st.tsize if k == (len(parts) - 1 if side == "above" else 0) else 0.0)
        return h + 0.2 * st.bond * max(0, len(parts) - 1)
    def _stack_parts(self, parts, cx, edge, side, texts, boxes, ha="center"):
        """Stack an arrow's text blocks and reagent structures away from the shaft (side 'above': upwards from edge,
        the part nearest the arrow last in the list; 'below': downwards). Gaps: 0.2 BU between parts; a structure next
        to the shaft keeps the distance a line of text keeps from it."""
        st = self.st; seq = list(reversed(parts)) if side == "above" else list(parts); first = True
        for kd, v in seq:
            if kd == "t":
                sh = self._ink_shift(v, side, st.tsize) * (1 if side == "above" else -1)
                g = self._rich_geo(v, cx, edge + sh, st.tsize, ha, "bottom" if side == "above" else "top"); texts.append(g)
                edge = (g["box"][3] if side == "above" else g["box"][1]) + (1 if side == "above" else -1) * 0.2 * st.bond
            else:
                pad = 0.22 * st.tsize if first else 0.0
                b = self._reagent_box(v, cx, edge + (pad if side == "above" else -pad), side, ha); boxes.append(b)
                edge = (b[3] if side == "above" else b[1]) + (1 if side == "above" else -1) * 0.2 * st.bond
            first = False
        return edge
    def _line_ink(self, l, size):
        """(height above, depth below) the baseline (pt) of one line's ink, sub- and superscripts included."""
        memo = self.__dict__.setdefault("_inkc", {}); k = (l, size)
        if k not in memo:
            ext = [self.F.path(t, x, y, sz, i, bd).get_extents() for t, x, y, sz, i, bd in self._lay(parse_rich(l), 0.0, 0.0, size)
                   if t.strip() and not self.F.missing(t)]
            memo[k] = (max([e.y1 for e in ext] + [0.0]), max([-e.y0 for e in ext] + [0.0]))
        return memo[k]
    def _rich_geo(self, s, x, y, size, ha="center", va="center", color=None):
        """Rich text at (x, y). Several lines (a list, or \n) stack 1.2 em apart, further where a subscript would come
        within 0.2 em of a superscript or ascender on the next line (Et2O over H3O+); va places the whole block (top:
        first line's top at y, bottom: last line's bottom at y), ha aligns every line."""
        lines = [_typo(l, chem=getattr(self.st, "chemtext", True)) for l in self._lines(s)] or [""]; capH = self.F.cap * size; n = len(lines)
        out, boxes, offs = [], [], [0.0]
        for i in range(n - 1):
            gap = max(1.2 * size, self._line_ink(lines[i], size)[1] + self._line_ink(lines[i + 1], size)[0] + 0.2 * size)
            offs.append(offs[-1] + (1.2 * size if gap - 1.2 * size < 1e-6 else gap))
        yb0 = {"center": y - capH / 2 + offs[-1] / 2, "baseline": y, "top": y - capH,
               "bottom": y + 0.22 * size + offs[-1]}[va]
        for i, l in enumerate(lines):
            runs = parse_rich(l); w = self._runs_w(runs, size); yb = yb0 - offs[i]
            x0 = {"center": x - w / 2, "left": x, "right": x - w}[ha]
            out += self._lay(runs, x0, yb, size); boxes.append((x0, yb - 0.22 * size, x0 + w, yb + capH + 0.1 * size))
        B = np.array(boxes); return {"runs": out, "box": (B[:, 0].min(), B[:, 1].min(), B[:, 2].max(), B[:, 3].max()), "color": color}
    def _item_geo(self, it):
        st = self.st
        if isinstance(it, Text):
            return self._rich_geo(it.s, it.x, it.y, it.size or (st.size if it.s == "+" else st.tsize), it.ha, it.va, it.color)
        L = (it.length or st.rxn) * st.bond; x, y = it.x, it.y; prims, texts = [], []; rboxes = []
        tl = self._lines(it.above) + self._lines(it.below)
        if tl and not it.vertical:                    # the arrow is always longer than its text (0.4 BU spare each side)
            L = max(L, max(self._runs_w(parse_rich(_typo(l, chem=getattr(st, "chemtext", True))), st.tsize) for l in tl) + 0.8 * st.bond)
        if it._mols and not it.vertical:              # ... and than any reagent structure it carries
            L = max([L] + [(lambda b: b[2] - b[0])(self._mol_box(m, curlies=False)) + 0.8 * st.bond for m in it._mols])
        if getattr(it, "branches", None) and not it.vertical:   # room for its branches: one in and one out on a side
            L = max(L, 2.6 * st.bond)                          # keep their labels CAP_SEP apart
            for sd, bs in self._arrow_branch_sides(it).items():
                if {b[1] for b in bs} == {"in", "out"}:
                    w = lambda lab: (lambda b: b[2] - b[0])(self._mol_box(lab, curlies=False) if isinstance(lab, Mol) else self._rich_geo(lab, 0, 0, st.tsize)["box"])
                    wi = max(w(b[0]) for b in bs if b[1] == "in"); wo = max(w(b[0]) for b in bs if b[1] == "out"); Lr = max(b[4] for b in bs) * st.bond
                    L = max(L, 3.8 * st.bond, (wi + wo) / 2 + self.CAP_SEP * st.bond + st.head_len - 0.35 * Lr + 0.35 * st.bond)
        split = it.vertical and getattr(it, "split", False)
        one = it.vertical and not split and (it._mols or bool(self._parts(it.above)) != bool(self._parts(it.below)))
        if split or one:                              # ... and, vertical, than the stacks centred beside it: split, above=
            stacks = [self._parts(it.above), self._parts(it.below)] if split else [self._parts(it.above) + self._parts(it.below)]
            L = max([L] + [self._stack_h(pa, "below") + 0.8 * st.bond for pa in stacks if pa])   # left and below= right;
        elif tl and it.vertical:                      # else one stack to the right (structures, or one side alone); text
            hb = lambda t: (self._rich_geo(t, 0, 0, st.tsize, "left", "top")["box"][3] - self._rich_geo(t, 0, 0, st.tsize, "left", "top")["box"][1]) if t else 0.0
            L = max(L, 2 * (max(hb(it.above), hb(it.below)) + 0.3 * st.tsize) + 0.6 * st.bond)   # on both sides: either
                                                                                                    # side of its middle
        if it.vertical and (it.kind != "->" or getattr(it, "dashed", False)):   # the horizontal arrow turned to point down
            hz = RArrow(kind=it.kind, length=L / st.bond, dashed=getattr(it, "dashed", False)); hz._B = it._B
            def rot(P): P = np.asarray(P, float); return np.column_stack([x + P[..., 1], y - P[..., 0]]) if P.ndim == 2 else np.array([x + P[1], y - P[0]])
            for pr in self._item_geo(hz)["prims"]:
                prims.append((pr[0], [rot(q) for q in pr[1]] if pr[0] == "line" else rot(pr[1])) + tuple(pr[2:]))
        elif it.vertical:
            prims.append(("line", [np.array([x, y]), np.array([x, y - L + 0.5 * st.head_len])], st.lw, st.ink, None, "butt"))
            prims.append(self._head(np.array([x, y - L]), np.array([0, -1.0]), 1.0))
        if it.vertical:                                    # conditions to the right, clear of the widest part
            xs = [float(np.max(np.asarray(q, float).reshape(-1, 2)[:, 0])) for pr in prims if pr[0] in ("line", "poly")
                  for q in ([np.vstack(pr[1])] if pr[0] == "line" else [pr[1]])]
            tx = x + (0.35 * st.bond if it.kind == "->" and not getattr(it, "dashed", False) else max(0.35 * st.bond, (max(xs) - x if xs else 0.0) + 0.22 * st.bond))
            if split or one:                               # stacks centred on the shaft, clear of the head's barbs
                ms_ = [float(np.min(np.asarray(q, float).reshape(-1, 2)[:, 0])) for pr in prims if pr[0] in ("line", "poly")
                       for q in ([np.vstack(pr[1])] if pr[0] == "line" else [pr[1]])]
                lx = x - max(0.35 * st.bond, (x - min(ms_) if ms_ else 0.0) + 0.22 * st.bond)
                tx = x + max(0.35 * st.bond, (max(xs) - x if xs else 0.0) + 0.22 * st.bond)
                sides = [(self._parts(it.above), lx, "right"), (self._parts(it.below), tx, "left")] if split else \
                        [(self._parts(it.above) + self._parts(it.below), tx, "left")]
                for pa, xx, ha in sides:
                    if pa: self._stack_parts(pa, xx, y - L / 2 + self._stack_h(pa, "below") / 2, "below", texts, rboxes, ha)
            else:
                if it.above: texts.append(self._rich_geo(it.above, tx, y - L / 2 + 0.3 * st.tsize + self._ink_shift(it.above, "above", st.tsize), st.tsize, "left", "bottom"))
                if it.below: texts.append(self._rich_geo(it.below, tx, y - L / 2 - 0.3 * st.tsize - self._ink_shift(it.below, "below", st.tsize), st.tsize, "left", "top"))
        else:
            a, b = np.array([x, y]), np.array([x + L, y])
            if it.kind in ("->", "<->", "<-", "-/->", "<-/-"):
                lft, rgt = it.kind in ("<->", "<-", "<-/-"), it.kind not in ("<-", "<-/-")
                dash = (4.0 * st.lw, 2.5 * st.lw) if getattr(it, "dashed", False) else None
                prims.append(("line", [a + (0.8 * st.head_len if lft else 0, 0), b - (0.8 * st.head_len if rgt else 0, 0)], st.lw, st.ink, dash, "butt"))
                if rgt: prims.append(self._head(b, np.array([1.0, 0]), 1.0))
                if lft: prims.append(self._head(a, np.array([-1.0, 0]), 1.0))
                if it.kind in ("-/->", "<-/-"):        # the cross of a reaction that does not happen
                    c = (a + b) / 2 + ((0.4 if it.kind == "<-/-" else -0.4) * st.head_len, 0); r = 0.2 * st.bond
                    prims.append(("line", [c + (-r, -r), c + (r, r)], st.lw, st.ink, None, "butt"))
                    prims.append(("line", [c + (-r, r), c + (r, -r)], st.lw, st.ink, None, "butt"))
            elif it.kind in ("=>", "<="):
                prims += self._open_arrow(a, b) if it.kind == "=>" else self._open_arrow(b, a)
            elif it.kind == "->>":                     # several steps: two short arrows, stacked and staggered
                sy = 0.3 * st.bond; l = 0.58 * L
                for x0, yy in ((x, y + sy), (x + L - l, y - sy)):
                    p1 = np.array([x0 + l, yy])
                    prims.append(("line", [np.array([x0, yy]), p1 - (0.8 * st.head_len, 0)], st.lw, st.ink, None, "butt"))
                    prims.append(self._head(p1, np.array([1.0, 0]), 1.0))
            elif it.kind in ("<=>", "<=>>", "<<=>"):
                o = np.array([0, 0.09 * st.bond])        # each half drawn as one outline, shaft and barb together:
                c, h_ = (a + b) / 2, (b - a) * 0.25      # a separate shaft would meet the barb's straight side in
                fa, fb = (c - h_, c + h_) if it.kind == "<<=>" else (a, b)   # a visible step once the lines are
                ra, rb = (c - h_, c + h_) if it.kind == "<=>>" else (a, b)   # pixel-snapped; the disfavoured half of
                prims.append(self._half_arrow(fa + o, fb + o))               # an unequal one half as long, centred
                prims.append(self._half_arrow(rb - o, ra - o))
            off = {"<=>": 0.09, "<=>>": 0.09, "<<=>": 0.09, "=>": 0.11, "<=": 0.11, "-/->": 0.13, "<-/-": 0.13, "->>": 0.3}.get(it.kind, 0.0) * st.bond   # text keeps its distance from the
            yt = y + 0.32 * st.bond + off                                           # nearest line (an equilibrium arrow's
            yl = y - 0.32 * st.bond - off                                           # upper half sits 1.3 pt higher), and
            ink = [np.asarray(pr[1], float).reshape(-1, 2)[:, 1] for pr in prims if pr[0] in ("line", "poly")]
            if ink:                                                                 # clears the ink of any raised barb
                yt = max(yt, max(float(v.max()) for v in ink) + 0.125 * st.bond)   # (an ordinary head already clears
                yl = min(yl, min(float(v.min()) for v in ink) - 0.125 * st.bond)   # by that: only a raised barb moves it)
            if it.step is not None:
                r = 0.5 * st.tsize; cy = yt + r
                texts.append(self._rich_geo(str(it.step), x + L / 2, cy, st.tsize, "center", "center"))
                prims.append(("circle", np.array([x + L / 2, cy]), r, st.ink, False)); yt = cy + r + 0.15 * st.tsize
            if it._mols:
                self._stack_parts(self._parts(it.above), x + L / 2, yt, "above", texts, rboxes)
                self._stack_parts(self._parts(it.below), x + L / 2, yl, "below", texts, rboxes)
            else:
                if it.above: texts.append(self._rich_geo(it.above, x + L / 2, yt + self._ink_shift(it.above, "above", st.tsize), st.tsize, "center", "bottom"))
                if it.below: texts.append(self._rich_geo(it.below, x + L / 2, yl - self._ink_shift(it.below, "below", st.tsize), st.tsize, "center", "top"))
        pts = [p for pr in prims if pr[0] in ("line", "poly") for p in pr[1]]
        for pr in prims:
            if pr[0] == "circle": pts += [pr[1] - pr[2], pr[1] + pr[2]]
        for t in texts: pts += [np.array(t["box"][:2]), np.array(t["box"][2:])]
        for bx in rboxes: pts += [np.array(bx[:2]), np.array(bx[2:])]
        heads = []                                    # (shaft, head polygon) pairs, for check()'s arrowhead test
        for i, pr in enumerate(prims[:-1]):
            if pr[0] == "line" and prims[i + 1][0] == "poly" and len(prims[i + 1][1]) in (3, 4):
                heads.append((np.array(pr[1], float), np.array(prims[i + 1][1], float)))
        if getattr(it, "branches", None) and not it.vertical:   # co-reagents curving in, by-products curving out
            rgt = it.kind not in ("<-", "<-/-", "<="); Tv = np.array([1.0 if rgt else -1.0, 0.0])
            s0 = np.array([x, y]) if rgt else np.array([x + L, y])
            for sd, bs in self._arrow_branch_sides(it).items():
                both = {b[1] for b in bs} == {"in", "out"}
                for lab, kind, at, side, reach in bs:
                    grp = [b for b in bs if b[1] == kind]; k_ = grp.index((lab, kind, at, side, reach))
                    r = np.array([0.0, 1.0 if sd == "above" else -1.0]); Lr = reach * st.bond
                    lo_, hi_ = 0.9 * Lr + 0.15 * st.bond, L - st.head_len - 0.55 * Lr - 0.2 * st.bond   # an 'in' end within the
                    d0 = (at * L if at is not None else (lo_ if kind == "in" else hi_) if both else 0.5 * L) + 0.2 * st.bond * (k_ - (len(grp) - 1) / 2)
                    d0 = max(d0, lo_) if kind == "in" else min(d0, hi_)          # arrow; an 'out' leaves before the head
                    P0 = s0 + Tv * d0
                    if kind == "out":
                        end = P0 + r * Lr + Tv * 0.9 * Lr; bp = self._curve_arrow(bezier([P0, P0 + Tv * 0.55 * Lr, end - r * 0.55 * Lr, end], 48), log=False)
                    else:
                        end = P0 + r * Lr - Tv * 0.9 * Lr; bp = [("line", list(bezier([end, end - r * 0.55 * Lr, P0 - Tv * 0.55 * Lr, P0], 48)), st.lw, st.ink, None, "butt")]
                    prims += bp; pts += [q for pr in bp for q in np.asarray(pr[1], float).reshape(-1, 2)]
                    q = end + r * 0.25 * st.bond
                    if isinstance(lab, Mol): b_ = self._reagent_box(lab, q[0], q[1], sd); rboxes.append(b_); pts += [np.array(b_[:2]), np.array(b_[2:])]
                    else:
                        g_ = self._rich_geo(lab, q[0], q[1] + (self._ink_shift(lab, "above", st.tsize) if sd == "above" else -self._ink_shift(lab, "below", st.tsize)),
                                            st.tsize, "center", "bottom" if sd == "above" else "top"); texts.append(g_)
                        pts += [np.array(g_["box"][:2]), np.array(g_["box"][2:])]
        P = np.array(pts); return {"prims": prims, "texts": texts, "box": (*P.min(0), *P.max(0)), "heads": heads}
    @staticmethod
    def _half_head(tip, d, h, barb, k, a, s):
        """A half arrowhead (one barb, on side s of the direction d) whose straight side continues the shaft's outer edge
        on the other side into a sharp point: the shaft (half-width a) ends inside it, so neither a rounded shaft end
        nor a step shows at the tip."""
        n = perp(d); E, Q = tip - n * a * s, tip - d * k * h - n * a * s          # (five points: check() tells half heads
        return [E, tip - d * h + n * barb * s, tip - d * k * h + n * a * s, Q, (E + Q) / 2]   # from full ones by that)
    def _half_arrow(self, tail, tip):
        """One half of an equilibrium arrow as a single outline: a shaft of line width from tail to tip, its outer edge
        running straight into a sharp point, and one barb on the left of the direction of travel (outward)."""
        st = self.st; d = unit(tip - tail); n = perp(d); a = st.lw / 2; h, k = st.head_len, st.head_notch
        return ("poly", [tail - n * a, tip - n * a, tip - d * h + n * st.head_w * 1.6, tip - d * k * h + n * a, tail + n * a], st.ink)
    def _head(self, tip, d, scale=1.0, half=0):
        h, w, k = self.st.head_len * scale, self.st.head_w * scale, self.st.head_notch; n = perp(d)
        if half: return ("poly", self._half_head(tip, d, h, w * 1.6, k, self.st.lw / 2, half), self.st.ink)
        return ("poly", [tip, tip - d * h + n * w, tip - d * k * h, tip - d * h - n * w], self.st.ink)

    def _settle_marks(self, curlies, mols):
        """Pass 1 of placing curly arrows: they land where they fit best, ignoring the movable marks (charges, spare lone
        pairs) on their target atoms; each lone pair then turns to its arrow and the marks near the arrows are re-placed
        clear of them (pass 2, in _build, fits the arrows again around the final marks)."""
        st = self.st
        dst = set()
        for c in curlies:                                  # the target atom, or both atoms of a target bond (a cation's
            d = self._edst(c)                              # charge beside the bond its electrons come back to moves aside)
            if isinstance(d, Atom): dst.add(id(d))
            elif isinstance(d, (tuple, list)) and len(d) == 2 and all(isinstance(x, Atom) for x in d): dst.update(id(x) for x in d)
        self._solve(curlies, mols, skip=frozenset(dst), second=False)
        for c in curlies:
            if getattr(c, "_lpc", None): l, t, _ = c._lpc; l.angle = float(t); l._locked = True
        self._arrow_geo = [np.vstack([c._arc, self._head_samples(c._P, c._d, c.kind, c._barb)]) for c in curlies]
        self._arrow_src = [c.src if isinstance(c.src, Atom) else None for c in curlies]
        near = {}                                          # (by the atom itself: the steps of a mechanism are copies
        for c in curlies:                                  # whose atoms share uids)
            cd = self._edst(c)
            if isinstance(cd, Atom):
                cd._avoid.append(ang(c._P - self.P(cd))); near[id(cd)] = cd
                if c._class in ("hook", "shift"):           # an arrow curling round onto its atom: a charge or lone pair
                    A = self.P(cd)                          # in the loop reads as caught in it, so the stretch near the
                    cd._avoid += sorted({10.0 * round(ang(q - A) / 10.0) for q in c._arc[::3] if norm(q - A) < 1.2 * st.bond})   # atom is taken
            if getattr(c, "_lpc", None): near[id(c._lpc[0].atom)] = c._lpc[0].atom
        for m in mols:
            for a in m.atoms:
                if (a.charge or a.radical or a.partial or a.number or a.lps) and \
                   min(float(np.hypot(*(g - self.P(a)).T).min()) for g in self._arrow_geo) < 2.0 * st.bond:
                    near[id(a)] = a
        pre = {u: (a._chg, a._rad) for u, a in near.items()}
        mi = {id(m): i for i, m in enumerate(mols)}
        for u in sorted(near, key=lambda u: (near[u].uid, mi.get(id(near[u].mol), -1))):
            a = near[u]; self._marks(a); geo = self._ctx_geo(a); fix = {}
            for nm, k, mk in (("charge_at", 0, a._chg), ("radical_at", 1, a._rad)):
                if pre[u][k] is not None and mk is not None and getattr(a, nm) is None and self._clear(a, mk[0], mk[1], geo) < 1.2:
                    fix[nm] = pre[u][k][2]                # no room left beside the arrow: the charge keeps its place
            if fix:
                keep = self._arrow_geo; self._arrow_geo = []
                for nm, t in fix.items(): setattr(a, nm, t)
                self._marks(a)
                for nm in fix: setattr(a, nm, None)
                self._arrow_geo = keep

    def _line_y(self, ref):
        """y (pt) of the line a row is centred on: its reaction arrow's shaft, or the '+' glyph's centre."""
        return ref.y if isinstance(ref, RArrow) else self._lboxes(ref)[2]
    def _recentre(self):
        """Re-centre row structures on their arrow line once the marks are final. A lone pair turned towards its arrow, or a
        charge moved aside for one, can change a structure's box after the row was laid out; each structure that is still
        where row() put it (not moved by hand since) is moved back onto the line. True if anything moved."""
        moved = False; B = self.st.bond
        for ref, rec in getattr(self, "_rows", []):
            line = self._line_y(ref)
            for r in rec:
                g, off0 = r
                if not g.mols or not g.mols[0].atoms or abs(g.mols[0].atoms[0].xy[1] * B - line - off0) > 1e-6: continue
                best = None
                for m in g.mols:
                    seen = set()
                    for a in m.atoms:
                        if a.uid in seen: continue
                        fr = m.frag(a); seen.update(x.uid for x in fr)
                        if best is None or len(fr) > len(best[1]): best = (m, fr)
                if best is None: continue
                mb = self._mol_box(best[0], curlies=False, atoms=best[1], prep=False); dy = line - (mb[1] + mb[3]) / 2
                if abs(dy) > 0.01:
                    self._imove(g, 0, dy); r[1] = off0 + dy; moved = True
        return moved
    # ---------- render
    def _state_sig(self):
        """Everything a build depends on that the user can set (public attributes of the style, structures, items,
        arrows, captions, pointers and extras), as one hashable value: a build is reused while it is unchanged."""
        seen = {}
        def cv(v, d=0):
            if isinstance(v, (str, int, float, bool, type(None))): return v
            if isinstance(v, np.ndarray): return ("nd", v.shape, tuple(np.round(v.astype(float), 9).ravel()))
            if isinstance(v, (np.floating, np.integer)): return float(v)
            if isinstance(v, Atom): return ("A", v.uid, id(v))
            if isinstance(v, Mol): return ("M", id(v))
            if isinstance(v, (list, tuple)): return tuple(cv(x, d) for x in v)
            if isinstance(v, (set, frozenset)): return tuple(sorted(repr(cv(x, d)) for x in v))
            if isinstance(v, dict): return tuple(sorted((repr(k), cv(x, d)) for k, x in v.items()))
            if d > 3 or not hasattr(v, "__dict__"): return ("o", type(v).__name__, id(v))
            if id(v) in seen: return ("r", id(v))
            seen[id(v)] = True
            return (type(v).__name__,) + tuple((k, cv(x, d + 1)) for k, x in sorted(vars(v).items())
                                                if not k.startswith("_") and k != "mol")
        out = [cv(self.st)]
        for m in self.mols:
            out.append(("mol", id(m), tuple(cv(a) for a in m.atoms), tuple(cv(b) for b in m.bonds)))
        for name in ("items", "curlies", "captions", "pointers", "extras"):
            out.append((name, tuple(cv(x) for x in getattr(self, name, []))))
        return hash(repr(out))
    def _build(self, _again=False):
        if not _again:                                         # unchanged since the last build (save() then check()):
            sig = self._state_sig(); last = getattr(self, "_built", None)   # the same drawing, so reuse it
            if last is not None and last[0] == sig: return list(last[1]), list(last[2])
        try: r = self._build_(_again)
        finally: self._drawn = False
        self._built = (self._state_sig(), r[0], r[1]); return r
    def _build_(self, _again=False):
        st = self.st; prims, texts = [], []; self._arrow_log = []
        self._ctx = None; self._drawn = False
        for it in self.items:                                  # reagent structures sit on their arrows before anything
            if isinstance(it, RArrow) and (it._mols or getattr(it, "_bmols", None)): self._item_geo(it)   # is measured or drawn, and structures at the
        for ex in self.extras:                                 # ends of cycle branches at theirs
            if ex[0] == "cycle" and any(isinstance(b[1], Mol) for b in ex[1]["branches"] + ex[1].get("feeds", [])):
                log = list(self._arrow_log); self._cycle_prims(ex[1], [], place=True); self._arrow_log = log
        self._unify_sides()
        for m in self.mols: _linear_centres(m)
        sph = {id(e[1]): e[6] for e in self.extras if e[0] == "orbital" and e[2] == "s"}
        for m in self.mols:
            for a in m.atoms:
                if getattr(a, "_newman_r", None) is not None: a.clear = a._newman_r * self.st.bond
                if getattr(a, "_sclear", False): a.clear = None; a._sclear = False
                if id(a) in sph and a.clear is None and not (a.drawn and not a.hidden):   # an s orbital on a bare atom:
                    a.clear = 0.42 * self.st.bond * sph[id(a)]; a._sclear = True           # its bonds stop at the sphere
                a._g = self._label(a) if (a.drawn and not a.hidden) else None; a._hE = a._hEc = a._hEw = None
        for m in self.mols: self._untangle_labels(m)          # a terminal OH crowding another label reads another way
        for m in self.mols:
            for a in m.atoms:
                for l in a.lps: l._locked = False
        self._ctx = list(self.mols)
        for m in self.mols: self._prep_mol(m)
        if self.curlies:
            self._settle_marks(self.curlies, self.mols)
            if not _again and self._recentre(): self._arrow_geo = []; return self._build_(_again=True)
            self._solve(self.curlies, self.mols)                # pass 2: arrows fitted around the final marks
            self._arrow_geo = []
        self._drawn = True                                     # marks final: later boxes measure them as they are
        lobes_ = []
        for ex in self.extras:
            if ex[0] == "orbital":
                op = self._orbital_prims(ex); prims += op       # (a lobe of a labelled atom starts at its label's edge;
                lobes_ += [(id(ex[1]), np.asarray(q[1], float)) for q in op if q[0] == "orb"]   # other labels it reaches
                g_ = self._orbital_label(ex)                                                    # get a white halo)
                if g_: texts.append({"runs": g_["runs"], "color": (ex[7] if len(ex) > 7 else None) or st.ink, "box": g_["box"]})
            elif ex[0] == "wiggle": prims += self._wiggle_prims(ex)
            elif ex[0] == "dipole": prims += self._dipole_prims(ex)
        for m in self.mols:
            for a in getattr(m, "_newman", []):
                f, r = a; prims.append(("circle", self.P(f), r * st.bond, st.ink, False))
        self._trim = True                                      # the drawn lines: a single bond at a centred C=O stops at
        for m in self.mols:                                    # its near line, meeting it in a clean corner
            for b in m.bonds:
                if b.order == 1 and b.style in ("plain", "bold", "dash", "hash") and any(
                        self._centred_trim(X, b, unit(self.P(Y) - self.P(X))) for X, Y in ((b.a, b.b), (b.b, b.a))):
                    self._bond(b)
        self._trim = False
        for m in self.mols:
            for b in m.bonds: prims += b._prims
            prims += self._joints(m) + list(getattr(m, "_xprims", []))
            texts += list(getattr(m, "_xtexts", []))
            for a in m.atoms:
                col = a.color or (HETERO.get(re.sub(r"\d", "", label_units(a.label)[0]), st.ink) if st.hetero and a._g else st.ink)
                for g in a._g or []:
                    t_ = {"runs": g["runs"], "color": col, "atom": True}
                    if lobes_ and any(o_ != id(a) and _MPath(P_).intersects_bbox(_Bbox.from_extents(*g["box"]), filled=True)
                                      for o_, P_ in lobes_):
                        t_["halo"] = True                 # a label in front of an orbital lobe stays legible
                    texts.append(t_)
                mc = a.color or st.ink                       # charges, radicals, lone pairs and δ take the atom's own
                if a._chg:                                   # colour (set by Mol.color), not a hetero=True element colour
                    c, r, t = a._chg; sgn = "+" if a.charge > 0 else "-"; k = abs(a.charge)
                    if st.charge == "circle" and k == 1:
                        prims.append(("circle", c, r, mc, False, self.RING_W * st.lw)); s = 0.58 * r + self.MARK_W * st.lw / 2
                        prims.append(("line", [c - (s, 0), c + (s, 0)], self.MARK_W * st.lw, mc, None, "butt"))   # flat ends,
                        if sgn == "+": prims.append(("line", [c - (0, s), c + (0, s)], self.MARK_W * st.lw, mc, None, "butt"))  # as bonds
                    else:
                        txt = (str(k) if k > 1 else "") + ("+" if sgn == "+" else "−")
                        g = self._rich_geo(txt, c[0], c[1], 0.8 * st.size); texts.append({"runs": g["runs"], "color": mc, "box": g["box"]})
                if a._rad: prims.append(("dot", a._rad[0], a._rad[1], mc))
                if a._pc:
                    g = self._rich_geo("δ" + ("+" if a.partial == "+" else "−"), a._pc[0][0], a._pc[0][1], 0.85 * st.size)
                    texts.append({"runs": g["runs"], "color": mc})
                if getattr(a, "_amb", None):
                    t_, sz_ = self._textmark(a, "amb"); g = self._rich_geo(t_, a._amb[0][0], a._amb[0][1], sz_)
                    texts.append({"runs": g["runs"], "color": mc})
                if a._num:
                    g = self._rich_geo(str(a.number), a._num[0][0], a._num[0][1], 0.7 * st.size)
                    texts.append({"runs": g["runs"], "color": st.ink})
                if a.site:                                   # a site dot on the vertex, its label beside it
                    col = a.site[0] or "#E8261A"
                    prims.append(("dot", self.P(a), self._halo_r(a), col, "bottom") if a.drawn else   # a halo under it all
                                 ("dot", self._site_xy(a), self.SITE_R * st.lw, col, "top"))   # over the bond ends
                    if a._site:
                        g = self._site_geo(a)
                        texts.append({"runs": g["runs"], "color": a.site[2] or col, "box": g["box"]})
                for l, c, d1, d2 in a._lp:
                    prims += [("dot", d1, st.lp_r, mc), ("dot", d2, st.lp_r, mc)]
        for c in self.curlies:
            col = c.color or st.arrow
            dash = (3.0 * st.alw, 2.0 * st.alw) if getattr(c, "dashed", False) else None
            prims.append(("line", list(c._pts), st.alw, col, dash, "butt")); prims.append(("poly", list(c._head), col))
        for it in self.items:
            g = self._item_geo(it)
            if isinstance(it, Text):
                if getattr(it, "highlight", None): prims.append(self._highlight(g, it.highlight))
                texts.append({"runs": g["runs"], "color": g.get("color") or st.ink, "free": it.s}); continue
            prims += g["prims"]; texts += [dict(t, color=t.get("color") or st.ink) for t in g["texts"]]
            for sh, hd in g.get("heads", []): self._arrow_log.append({"shaft": sh, "poly": hd, "lw": st.lw, "curved": False})
        def cpts(m):
            ids = {id(x) for x in m} if isinstance(m, tuple) else {id(m)}
            return [np.vstack([c._pts, c._head]) for c in self.curlies if ids & self._owners(c)]
        sbx = {}
        def sbox(m):
            if id(m) not in sbx: sbx[id(m)] = self._sbox(m)
            return sbx[id(m)]
        for cp in self.captions:
            while len(cp) < 6: cp.append(None)
        G = self._caps_geo(self.captions, sbox, lambda cp: cpts(cp[0]), align=True)
        for cp in self.captions:
            g = G[id(cp)]; texts.append({"runs": g["runs"], "color": cp[5] or st.ink, "box": g["box"], "cap": cp})
            if len(cp) > 7 and cp[7]: prims.append(self._highlight(g, cp[7]))
        for p in self.pointers:
            g = self._ptr_geo(p); prims += g["prims"]
            texts += [dict(t, color=t.get("color") or st.ink, note=p) for t in g["texts"]]
        for ex in self.extras:
            if ex[0] == "ringtext":
                c = np.mean([self.P(a) for a in ex[1]], axis=0)
                g = self._rich_geo(ex[2], c[0], c[1], ex[4] or st.tsize, "center", "center", ex[3])
                texts.append(dict(g, color=ex[3] or st.ink))
            if ex[0] == "brackets": prims += self._brackets(ex, texts)
            elif ex[0] == "frame": prims += self._frame_prims(ex, texts)
            elif ex[0] == "cycle": prims += self._cycle_prims(ex[1], texts)
        self._ctx = None
        return prims, texts
