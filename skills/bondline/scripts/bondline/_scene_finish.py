"""bondline engine: Scene methods, part 3 of 3."""
# Part of the bondline engine (see __init__.py): the package's modules share one namespace, so names
# defined in any of them are visible here at call time, as in a single module.
from ._core import *  # noqa: F401,F403


class _SceneFinish:
    def _site_xy(self, a):
        """Centre (pt) of atom a's site dot: on the vertex, except at the bare end of an offset double bond (a terminal
        CH2=), where the vertex is the end of the main line only: there it sits on the middle of the pair's end and caps
        both lines rather than hanging off one."""
        P = self.P(a); nb = [x for x in a.mol.nbrs(a) if not getattr(x, "hidden", False)]
        if a.drawn or len(nb) != 1: return P
        b = a.mol.find(a, nb[0])
        if b is None or b.order != 2 or b.style not in ("plain", "bold"): return P
        s = self._dside(b)
        if s == 0: return P
        return P + perp(unit(self.P(b.b) - self.P(b.a))) * s * self.st.dbl / 2
    def _ptr_target(self, t, d):
        """(point, gap) for a pointer aimed at target t from direction d (unit, from the target towards the note): how far
        from the point its tip stops. The head, as drawn coming in along -d, keeps clear of the atom's label letters
        (1.5 pt: the tip stops at the letters it points at, not beyond the label's far end), or of the bonds round a
        bare atom (1.4 pt: a crowded vertex is pointed at, not run into), its site dot, or the bond's lines."""
        st = self.st; d = unit(np.asarray(d, float))
        if isinstance(t, Atom):
            T = self._site_xy(t) if t.site else self.P(t)
            lo = (self.SITE_R * st.lw + 2.2) if t.site else 0.2 * st.bond
            if t.drawn and t._g:
                hs = [h for g in t._g for h in g["hulls"]]
                E = np.vstack([np.hstack([h, np.roll(h, -1, axis=0)]) for h in hs]); s0 = np.cumsum([0] + [len(h) for h in hs[:-1]]).astype(int)
                return T, self._head_clear(T, d, lambda Q: _polys_dist(Q, E, s0).min(1), 1.5, 2.0)
            S = [np.r_[P[i], P[i + 1]] for b in t.mol.bonds if t in (b.a, b.b)
                 for pr in (getattr(b, "_prims", None) or [("line", [self.P(b.a), self.P(b.b)])])   # (not drawn yet: centre line)
                 if pr[0] in ("line", "poly") for P in [np.asarray(pr[1], float).reshape(-1, 2)] for i in range(len(P) - 1)]
            if not S: return T, lo
            S = np.array(S); return T, self._head_clear(T, d, lambda Q: _segs_dist(Q, S).min(1) - st.lw / 2, 1.4, lo, 0.6)
        if isinstance(t, Mid): return self._pt(t), 2.0
        a, b = t; T = (self.P(a) + self.P(b)) / 2; bd = a.mol.find(a, b)
        S = [np.r_[P[i], P[i + 1]] for pr in (bd._prims if bd is not None else []) if pr[0] in ("line", "poly")
             for P in [np.asarray(pr[1], float).reshape(-1, 2)] for i in range(len(P) - 1)]
        if not S: return T, 2.0 + st.lw / 2 + (st.dbl if bd is not None and bd.order >= 2 else 0.0)
        S = np.array(S); return T, self._head_clear(T, d, lambda Q: _segs_dist(Q, S).min(1) - st.lw / 2, 2.0, 0.0)
    PTR_HS = 0.8                  # a pointer's head is this size of a reaction arrow's
    def _head_clear(self, T, d, dist, need, lo, cap=1.5):
        """Least distance s >= lo (pt) along d from T at which a pointer's head, its tip at T + s d and pointing back along
        -d, keeps `need` pt from the ink that dist(points) measures; at most `cap` BU beyond lo (a pointer that comes in
        along a bond stops on it, where check() sees it, not past the next atom)."""
        st = self.st; tip0 = np.zeros(2); H = np.asarray(self._head(tip0, -d, self.PTR_HS)[1], float)
        rel = np.vstack([H[i] + (H[(i + 1) % len(H)] - H[i]) * f for i in range(len(H)) for f in np.linspace(0, 1, 5, endpoint=False)])
        S = lo + np.arange(0.0, cap * st.bond, 0.25); Q = (T + S[:, None, None] * d + rel[None]).reshape(-1, 2)
        ok = dist(Q).reshape(len(S), -1).min(1) >= need
        return float(S[int(np.argmax(ok))]) if ok.any() else float(S[-1])
    PTR_SIDES = {"below": (0, -1), "above": (0, 1), "left": (-1, 0), "right": (1, 0)}
    def _ptr_rays(self, p):
        """Segments (pt) along which pointer p's arrows will run, from its targets outward: known before the marks are
        placed (the arrows' lengths are not), so that marks keep out of their way."""
        B = self.st.bond; out = []
        u = np.array(self.PTR_SIDES[p.side], float) if isinstance(p.side, str) else dirv(float(p.side))
        far = 8.0 if isinstance(p.side, str) else p.length + 1.2
        for t in p.targets:
            T = self.P(t) if isinstance(t, Atom) else self._pt(t) if isinstance(t, Mid) else (self.P(t[0]) + self.P(t[1])) / 2
            out.append((T + u * 0.25 * B, T + u * far * B))
        return out
    def _ptr_ink(self, p):
        """What pointer p's note keeps clear of (pt): the line segments and the boxes of its molecules' ink (bonds,
        labels, charges and site labels, site dots, lone pairs)."""
        segs, boxes = [], []; st = self.st
        for m in p.mols():
            for b in m.bonds:
                prs = getattr(b, "_prims", None)               # (not drawn yet, when a row is laid out: its centre line)
                if prs is None: segs.append((self.P(b.a), self.P(b.b))); continue
                for pr in prs:
                    P = np.asarray(pr[1], float).reshape(-1, 2) if pr[0] in ("line", "poly") else None
                    if P is None: continue
                    if pr[0] == "line": segs += [(P[i], P[i + 1]) for i in range(len(P) - 1)]
                    else: boxes.append((*P.min(0), *P.max(0)))
            for a in m.atoms:
                boxes += [g["box"] for g in a._g or []]
                for k in ("_chg", "_rad", "_pc", "_num", "_site", "_amb"):
                    bk = self._mark_box(a, k)
                    if bk: boxes.append(bk)
                if a.site:
                    c, r = self._site_xy(a), self.SITE_R * st.lw; boxes.append((c[0] - r, c[1] - r, c[0] + r, c[1] + r))
                for _, c, d1, d2 in getattr(a, "_lp", []):
                    for d in (d1, d2): boxes.append((d[0] - st.lp_r, d[1] - st.lp_r, d[0] + st.lp_r, d[1] + st.lp_r))
            for pr in getattr(m, "_xprims", []):
                P = np.array(pr[1], float).reshape(-1, 2) if pr[0] in ("line", "poly") else np.array([pr[1] - pr[2], pr[1] + pr[2]])
                boxes.append((*P.min(0), *P.max(0)))
        return segs, boxes
    @staticmethod
    def _ink_extreme(segs, boxes, ax, lo, hi, sg):
        """Furthest reach (sg = +1: largest, -1: smallest coordinate on axis ax) of the ink lying within [lo, hi] on the
        other axis; None if none does."""
        o = 1 - ax; vals = []
        for p0, p1 in segs:
            a0, a1 = p0[o], p1[o]
            if max(a0, a1) < lo or min(a0, a1) > hi: continue
            for q in (p0, p1):
                if lo <= q[o] <= hi: vals.append(q[ax])
            if a1 != a0:
                for e in (lo, hi):
                    t = (e - a0) / (a1 - a0)
                    if 0 <= t <= 1: vals.append(p0[ax] + t * (p1[ax] - p0[ax]))
        for bx in boxes:
            if bx[2 + o] < lo or bx[o] > hi: continue
            vals.append(bx[2 + ax] if sg > 0 else bx[ax])
        return (max if sg > 0 else min)(vals) if vals else None
    def _note_clash(self, b1, b2):
        """Two notes too close: on one line closer than CAP_SEP (they would read as one), or corner to corner closer than
        0.5 BU across and 0.2 BU up or down."""
        B = self.st.bond; dx = max(b2[0] - b1[2], b1[0] - b2[2]); dy = max(b2[1] - b1[3], b1[1] - b2[3])
        return (dy < -0.1 * B and dx < self.CAP_SEP * B) or (dx < 0.5 * B and dy < 0.2 * B)
    def _ptr_geo(self, p, _memo=None):
        """Arrows and note of a pointer (pt): {'prims', 'texts', 'box'}. Its note keeps clear of the notes of the
        pointers declared before it: moved aside along the note line (its arrows then run from the note's nearest
        point), or, for a callout, further out."""
        st = self.st; B = st.bond; col = p.color or st.ink; prims, texts = [], []
        dash = (3.0 * st.lw, 2.0 * st.lw) if p.dashed else None; hs = self.PTR_HS
        memo = {} if _memo is None else _memo; avoid = []
        for q in self.pointers:
            if q is p: break
            if id(q) not in memo: memo[id(q)] = self._ptr_geo(q, memo)
            avoid += [t["box"] for t in memo[id(q)]["texts"]]
        def arrow(pts):
            P = np.asarray(pts, float); tip = P[-1]; h = st.head_len * hs; kk = st.head_notch
            far_ = np.flatnonzero(np.hypot(*(P - tip).T) >= kk * h)
            if len(P) > 2 and len(far_):                  # a curve: the head sits on the chord from where the shaft
                j = int(far_[-1]); a_, b_ = P[j], P[min(j + 1, len(P) - 1)]; lo_, hi_ = 0.0, 1.0   # meets it (as in
                for _ in range(30):                                                               # _curve_arrow): no
                    mid = (lo_ + hi_) / 2                                                         # kink, however short
                    if norm(a_ + (b_ - a_) * mid - tip) >= kk * h: lo_ = mid
                    else: hi_ = mid
                cut = a_ + (b_ - a_) * lo_; d = unit(tip - cut); shaft = list(P[:j + 1]) + [cut, tip - d * 0.8 * h]
            else: d = unit(tip - P[0]); shaft = [P[0], tip - d * 0.8 * h]
            prims.append(("line", shaft, st.lw, col, dash, "butt")); prims.append(("poly", self._head(tip, d, hs)[1], col))
        def tip_for(t, v):                                # (point, tip) on target t for an arrow coming in from v
            T_, g_ = self._ptr_target(t, v); return T_, T_ + v * g_
        def clash(bx): return any(self._note_clash(bx, b2) for b2 in avoid)
        segs, boxes = self._ptr_ink(p)
        if isinstance(p.side, str):
            s_ = np.array(self.PTR_SIDES[p.side], float); ax = 1 if s_[1] else 0; sg = s_[ax]; o = 1 - ax; eo = np.eye(2)[o]
            Ts = [self._ptr_target(t, s_)[0] for t in p.targets]; c_o = float(np.mean([T[o] for T in Ts])); plan = []
            for t, T in zip(p.targets, Ts):               # a bond lying along the pointer's direction is met from one
                v = s_                                    # side at 35 deg, not run along: the outer side of the
                if isinstance(t, tuple) and abs(float(unit(self.P(t[1]) - self.P(t[0])) @ s_)) > math.cos(math.radians(30)):
                    so = np.sign(T[o] - c_o) if abs(T[o] - c_o) > 1e-6 else self._ptr_room_side(p, T, eo)   # targets,
                    v = unit(s_ + eo * so * math.tan(math.radians(35)))                                       # or the
                plan.append((t, *tip_for(t, v), v))                                                           # freer one
            ends = [e for _, _, e, _ in plan]
            far = (min if sg < 0 else max)(e[ax] for e in ends) + sg * p.length * B
            tails_at = lambda f: [e + v * ((f - e[ax]) / v[ax]) for _, _, e, v in plan]
            g0 = self._rich_geo(p.text, 0.0, 0.0, st.tsize, "center", "center", col) if p.text else None
            for _ in range(2):                            # the tails stand clear of everything in the strip the
                tl = tails_at(far); c0 = float(np.mean([q[o] for q in tl]))     # arrows and the note take
                span = [min(q[o] for q in ends + tl), max(q[o] for q in ends + tl)]
                if g0: half = (g0["box"][2 + o] - g0["box"][o]) / 2; span = [min(span[0], c0 - half), max(span[1], c0 + half)]
                ext = self._ink_extreme(segs, boxes, ax, span[0] - 0.1 * B, span[1] + 0.1 * B, sg)
                if ext is not None: far = (min if sg < 0 else max)(far, ext + sg * 0.25 * B)
            tl = tails_at(far); c0 = float(np.mean([q[o] for q in tl])); moved = False
            if p.text:
                def note(c):
                    if ax == 1: return self._rich_geo(p.text, c, far + sg * 0.15 * B, st.tsize, "center", "top" if sg < 0 else "bottom", col)
                    return self._rich_geo(p.text, far + sg * 0.15 * B, c, st.tsize, "right" if sg < 0 else "left", "center", col)
                g = note(c0); hit = [b2 for b2 in avoid if self._note_clash(g["box"], b2)]
                if hit:                                   # aside from an earlier note on the same side: the way whose
                    tg_ = {id(x) for t in p.targets for x in ((t.a, t.b) if isinstance(t, Mid) else t if isinstance(t, tuple) else (t,))}
                    lab = [x["box"] for m in p.mols() for a in m.atoms if id(a) not in tg_ for x in (a._g or [])]
                    E_ = np.array([[*q0, *q1] for q0, q1 in segs], float).reshape(-1, 4)
                    def paths(bx):                        # arrows (from the note's nearest point) cross least ink
                        out = []
                        for t, T, e, v in plan:
                            q = np.array([min(max(T[0], bx[0]), bx[2]), min(max(T[1], bx[1]), bx[3])]); q = q + unit(T - q) * 1.5
                            out.append((q, tip_for(t, unit(q - T))[1]))
                        return out
                    Ts_ = [T for _, T, _, _ in plan]
                    def cost(bx):                         # bonds crossed or run along, labels and earlier notes entered
                        n = 0
                        for (q, e), T in zip(paths(bx), Ts_):
                            n += sum(1 for r in E_ if _seg_x(q, e, r[:2], r[2:])) + sum(1 for b2 in lab + avoid if seg_box_gap(q, e, b2) < 1.5)
                            Q = np.array([q + (e - q) * f for f in np.linspace(0, 1, 24)]); Q = Q[np.hypot(*(Q - T).T) > 0.45 * B]
                            if len(Q) and len(E_) and _segs_dist(Q, E_).min() < st.lw + 0.6: n += 1
                        return n
                    bc = (hit[0][o] + hit[0][2 + o]) / 2; pref = 1 if c_o >= bc else -1; best = None
                    for way in (pref, -pref):
                        c, gx = c0, g
                        for _ in range(len(avoid) + 1):
                            hx = [b2 for b2 in avoid if self._note_clash(gx["box"], b2)]
                            if not hx: break
                            b2 = max(hx, key=lambda b_: way * b_[2 + o] if way > 0 else -b_[o]); gb = gx["box"]
                            c += (b2[2 + o] + self.CAP_SEP * B - gb[o] + 0.01) if way > 0 else (b2[o] - self.CAP_SEP * B - gb[2 + o] - 0.01)
                            gx = note(c)
                        for k in range(3):                # at the nearest place clear of the notes, or a little further
                            gk = note(c + way * k * B)
                            if any(self._note_clash(gk["box"], b2) for b2 in avoid): continue
                            sc = (cost(gk["box"]), k)
                            if best is None or sc < best[0]: best = (sc, c + way * k * B, gk)
                    if best is not None: c0, g, moved = best[1], best[2], True
                texts.append(g)
            if moved:                                     # from the note's nearest point, straight to each target
                bx = g["box"]
                for t, T, e, v in plan:
                    q = np.array([min(max(T[0], bx[0]), bx[2]), min(max(T[1], bx[1]), bx[3])]); q = q + unit(T - q) * 1.5
                    arrow([q, tip_for(t, unit(q - T))[1]])
            else:
                for (t, T, e, v), q in zip(plan, tl): arrow([q, e])
        else:
            u = dirv(float(p.side)); T, gap = self._ptr_target(p.targets[0], u)
            ha = "left" if u[0] > 0.45 else "right" if u[0] < -0.45 else "center"
            va = "bottom" if u[1] > 0.45 else "top" if u[1] < -0.45 else "center"
            cv = bool(p.bend and p.text)                  # a curved callout: the note stands to one side, its near end
            if cv:                                        # short of the target, so the arrow can leave that end level
                sx = (1 if u[0] > 0 else -1) if abs(u[0]) > 0.05 else (1 if p.bend < 0 else -1)
                ha = "left" if sx > 0 else "right"
            E = np.array([[*q0, *q1] for q0, q1 in segs], float).reshape(-1, 4)
            labs = {tuple(g_["box"]) for m_ in p.mols() for a_ in m_.atoms for g_ in (a_._g or [])}
            def free(bx):                                                       # white between a note box and the ink;
                def lv(b2):                                                     # a label beside the note, on a line of
                    if tuple(b2) not in labs or max(b2[0] - bx[2], bx[0] - b2[2]) <= 0: return 0.0   # it, keeps twice
                    ov = min(bx[3], b2[3]) - max(bx[1], b2[1])                  # as much (OTMS  formation would read as
                    return self.TEXT_CLEAR * B / 14.4 if ov > 0.5 * (b2[3] - b2[1]) else 0.0   # one line)
                return min([box_gap(bx, b2) - lv(b2) for b2 in boxes] + [seg_box_gap(q[:2], q[2:], bx) - st.lw / 2 for q in E] + [99.0])
            g0 = self._rich_geo(p.text, 0.0, 0.0, st.tsize, ha, va, col) if p.text else None
            def at(k, da):                                                      # out along u, turned by da: the note
                ref = T + dirv(float(p.side) + da) * (gap + (p.length + 0.1 * k) * B)   # keeps TEXT_CLEAR from all
                if cv: ref[0] = max(ref[0], T[0] + 0.3 * B) if sx > 0 else min(ref[0], T[0] - 0.3 * B)   # ink (and is
                return ref                                                      # clear of earlier notes), at the
            ok = lambda ref: g0 is None or (free(tuple(np.add(g0["box"], [*ref, *ref]))) >= self.TEXT_CLEAR * B / 14.4
                                            and not clash(tuple(np.add(g0["box"], [*ref, *ref]))))
            cands = sorted((k + 2 * abs(da) / 15, abs(da), k, da) for k in range(25)   # nearest place, a turn of 15
                           for da in ((0,) if g0 is None else (0, -15, 15, -30, 30, -45, 45)))   # deg costing 0.2 BU out
            ref = next((r_ for _, _, k, da in cands for r_ in (at(k, da),) if ok(r_)), None)
            if ref is None:                                                     # no clear place: the freest one
                ref = max((at(k, da) for _, _, k, da in cands), key=lambda r_: free(tuple(np.add(g0["box"], [*r_, *r_]))) - 2.0 * clash(tuple(np.add(g0["box"], [*r_, *r_]))))
            g = self._rich_geo(p.text, ref[0], ref[1], st.tsize, ha, va, col) if p.text else None
            if p.text:
                texts.append(g); bx = g["box"]
                tail = np.array([min(max(T[0], bx[0]), bx[2]), min(max(T[1], bx[1]), bx[3])])
                tail = tail if norm(tail - T) > 0 else ref
                tail = tail + unit(T - tail) * 1.5
            else: tail = ref
            tip = tip_for(p.targets[0], unit(tail - T))[1]
            if p.bend and p.text:                         # a curved callout, as in teaching figures: it leaves the end of the
                bx = g["box"]; capH = self.F.cap * st.tsize   # note's line nearest the target, level, and sweeps round to it
                mids = sorted({round(r[2] + capH / 2, 3) for r in g["runs"] if r[0].strip()})
                yl = min(mids, key=lambda y: abs(y - T[1])) if mids else (bx[1] + bx[3]) / 2
                left = sx > 0                             # the end facing the target, heading towards it
                P0 = np.array([bx[0] - 1.5 if left else bx[2] + 1.5, yl]); t0 = np.array([-1.0 if left else 1.0, 0.0])
                D = norm(T - P0); t1 = unit(T - (P0 + t0 * 0.5 * D)); tp = tip_for(p.targets[0], -t1)[1]
                level = bezier([P0, P0 + t0 * 0.55 * D, tp - t1 * 0.45 * D, tp], 48)
                def room(P):                              # white between the curve and the ink it passes (its last
                    Q = [q for q in P if norm(q - T) > gap + 0.25 * B]          # stretch, into the target, aside)
                    if not Q: return 99.0
                    Q = np.array(Q); d = 99.0
                    for q0, q1 in segs:
                        ab = q1 - q0; t = np.clip(((Q - q0) @ ab) / max(ab @ ab, 1e-9), 0, 1)
                        d = min(d, float(np.hypot(*(q0 + t[:, None] * ab - Q).T).min()) - st.lw)
                    for b2 in boxes:
                        c = np.clip(Q, [b2[0], b2[1]], [b2[2], b2[3]]); d = min(d, float(np.hypot(*(Q - c).T).min()) - st.lw / 2)
                    return d
                plain = arc(tail, tip, p.bend)            # else from the note's nearest edge, an even arc, its tip set
                tq = tip_for(p.targets[0], unit(plain[-2] - plain[-1]))[1]        # for the way it comes in
                plain = arc(tail, tq, p.bend)
                best = level if room(level) >= 1.5 or room(level) >= room(plain) else plain
                cur = [np.vstack([c._pts, c._head]) for c in self.curlies
                       if getattr(c, "_pts", None) is not None and getattr(c, "_head", None) is not None]
                def curl_room(P):                         # white between the whole curve, tip included, and the curly
                    P = np.asarray(P, float); d = 99.0    # arrows (two heads at one atom must not run together)
                    for Cp in cur: d = min(d, float(_segs_dist(P, np.hstack([Cp[:-1], Cp[1:]])).min()) - st.lw)
                    return d
                tg = p.targets[0]; ats = [tg] if isinstance(tg, Atom) else []      # the target's own marks: its head
                mk = [bx for a_ in ats for k_ in ("_chg", "_rad", "_pc") for bx in [self._mark_box(a_, k_)] if bx]   # keeps
                mk += [(q[0] - st.lp_r, q[1] - st.lp_r, q[0] + st.lp_r, q[1] + st.lp_r)                        # off a
                       for a_ in ats for _, _c, d1, d2 in getattr(a_, "_lp", []) for q in (d1, d2)]           # charge too
                def tip_room(P):
                    if not mk: return 99.0
                    P = np.asarray(P, float); Q = P[np.hypot(*(P - P[-1]).T) <= 1.2 * st.head_len * hs]
                    return min(box_dist(q, bx) for q in Q for bx in mk) - st.head_w * hs / 2
                if (cur and curl_room(best) < 2.5) or tip_room(best) < 1.0:   # in from another side: the nearest turn
                    def level_at(r):                                          # that clears them
                        t1r = dirv(ang(t1) + r); tpr = tip_for(p.targets[0], -t1r)[1]
                        return bezier([P0, P0 + t0 * 0.55 * D, tpr - t1r * 0.45 * D, tpr], 48)
                    opts = [(0, best)] + [(abs(r), level_at(r)) for r in (10, -10, 20, -20, 35, -35, 50, -50, 75, -75)] + [(0, plain)]
                    def worst(P): return min((curl_room(P) - 2.5) if cur else 9.0, tip_room(P) - 1.0, 9.0)   # both kept,
                    best = max(opts, key=lambda o: (worst(o[1]) >= 0, room(o[1]) >= 1.5, -o[0] if worst(o[1]) >= 0   # else
                                                    else worst(o[1])))[1]                                   # the least bad
                arrow(best)
            elif p.bend: arrow(arc(tail, tip, p.bend))     # no note: an even arc
            else: arrow([tail, tip])
        pts = [q for pr in prims for q in (pr[1] if pr[0] in ("line", "poly") else [])]
        for t in texts: pts += [np.array(t["box"][:2]), np.array(t["box"][2:])]
        P = np.array(pts, float); return {"prims": prims, "texts": texts, "box": (*P.min(0), *P.max(0))}
    def _dipole_prims(self, ex):
        """The dipole arrow's shaft, head and crossbar (S.dipole): 0.8 of the bond long, 0.45 BU off it."""
        st = self.st; _, a, b, col, side = ex; A, B = self.P(a), self.P(b); u = unit(B - A); n = perp(u); L = norm(B - A)
        if side in ("left", "right"): sg = 1 if side == "left" else -1
        else:                                          # the side whose midpoint stands further from the other atoms
            pts = np.array([self.P(x) for x in a.mol.atoms if x is not a and x is not b and not x.hidden] or [[1e6, 1e6]])
            sg = max((1, -1), key=lambda k: float(np.hypot(*(pts - ((A + B) / 2 + n * k * 0.45 * st.bond)).T).min()))
        o = n * sg * 0.45 * st.bond; tail, tip = A + u * 0.1 * L + o, A + u * 0.9 * L + o; col = col or st.ink
        hd = self._head(tip, u, 0.5); hd = (hd[0], hd[1], col)
        c = tail + u * 0.16 * L; h = 0.11 * st.bond
        return [("line", [tail, tip - u * 0.5 * st.head_len * 0.5], st.lw, col, None, "butt"), hd,
                ("line", [c - n * h, c + n * h], st.lw, col, None, "butt")]
    def _wiggle_prims(self, ex):
        """The wiggle's line (S.wiggle): a sine wave (or a dashed line) square to the bond through the point at
        fraction at from a to b, 1.5 waves across a bond's middle, 2 at an end."""
        st = self.st; _, a, b, at, col, dashed, length = ex; A, B = self.P(a), self.P(b)
        u = unit(B - A); n = perp(u); end = at >= 0.99 or at <= 0.01
        Lh = 0.5 * (length if length is not None else (0.5 if end else 0.75)) * st.bond
        p = A + (B - A) * at; col = col or st.ink
        bd = a.mol.find(a, b); vis = getattr(bd, "_vis", None) if bd is not None else None
        if vis is not None and not end:                  # across the bond's drawn part (between its labels' edges)
            p0, p1 = (vis[0], vis[1]) if bd.a is a else (vis[1], vis[0]); p = np.asarray(p0) + (np.asarray(p1) - np.asarray(p0)) * at
        if dashed:
            return [("line", [p - n * Lh, p + n * Lh], st.lw, col, (3.0 * st.lw, 2.0 * st.lw), "butt")]
        t = np.linspace(-Lh, Lh, 80); waves = 2.0 if end else 1.5; amp = 0.075 * st.bond
        P = p + np.outer(t, n) + np.outer(amp * np.sin(2 * math.pi * waves * (t + Lh) / (2 * Lh)), u)
        return [("line", list(P), st.lw, col, None, "round")]
    def _ptr_room_side(self, p, T, eo):
        """+1 or -1 along eo: the side of point T with more room (its nearest atom further away)."""
        A = np.array([self.P(a) for m in p.mols() for a in m.atoms]); best = (-1.0, 1)
        for so in (1, -1):
            q = T + eo * so * 0.8 * self.st.bond; d = float(np.hypot(*(A - q).T).min()) if len(A) else 99.0
            if d > best[0] + 1e-6: best = (d, so)
        return best[1]
    def _orbital_label(self, ex):
        """Rich-text geometry of an orbital's name (S.orbital label=): just past the tip of its lobe along angle (the
        sphere's edge for s), on the side the lobe points to."""
        if len(ex) < 9 or not ex[8]: return None
        st = self.st; a, kind, angle, size = ex[1], ex[2], ex[3], ex[6]; A = self.P(a); u = dirv(angle)
        lab = a.drawn and not a.hidden and a._g
        off = self._lobe_tip(a, u, A, size) if (lab and kind != "s") else 0.0
        reach = 0.42 * st.bond * size if kind == "s" else off + self.LOBE_L * st.bond * size
        q = A + u * (reach + 0.12 * st.size)
        ha = "left" if u[0] > 0.3 else "right" if u[0] < -0.3 else "center"
        va = "bottom" if u[1] > 0.3 else "top" if u[1] < -0.3 else "center"
        col = ex[7] if len(ex) > 7 else None; T = A + u * reach; need = 0.15 * st.size
        g = self._rich_geo(ex[8], q[0], q[1], st.tsize, ha, va, col)
        for _ in range(3):                                 # its ink (a subscript's descent too) 1.5 pt clear of the tip
            ext = [self.F.path(t, x, y, sz, i, bd).get_extents() for t, x, y, sz, i, bd in g["runs"] if t.strip() and not self.F.missing(t)]
            if not ext: break
            x0, y0 = min(e.x0 for e in ext), min(e.y0 for e in ext); x1, y1 = max(e.x1 for e in ext), max(e.y1 for e in ext)
            gap = math.hypot(max(x0 - T[0], 0.0, T[0] - x1), max(y0 - T[1], 0.0, T[1] - y1))
            if gap >= need - 0.05: break
            q = q + u * (need - gap + 0.1); g = self._rich_geo(ex[8], q[0], q[1], st.tsize, ha, va, col)
        return g
    def _lobe_form(self, a, t):
        """(q, tip, cap) of the lobe along t on atom a: the full teardrop, slimmer where a bond of the atom leaves
        within 50 deg of its axis (a DMSO sulfur's p orbital between its methyls), so the two never run side by side."""
        d = min([abs(adiff(ang(self.P(n) - self.P(a)), t)) for n in a.mol.nbrs(a)] + [90.0])
        if d >= 50: return self.LOBE_Q, self.LOBE_TIP, self.LOBE_CAP
        tip = max(10.0, min(self.LOBE_TIP, d - 14)); return max(0.4, self.LOBE_Q - 0.009 * (50 - d)), tip, min(self.LOBE_CAP, tip - 4)
    def _lobe_half(self, L, q=None, soft=True, form=None):
        """Upper half of a lobe's outline (pt; x along its axis from the tip, y across), tip to far end: a cubic from the
        tip, leaving it LOBE_TIP deg off the axis, into a round far end q as wide as the lobe is long, met LOBE_CAP deg off
        the axis (soft=False: straight sides, tangent to the end; form: (q, tip, cap) from _lobe_form)."""
        if form is not None: q = form[0]
        q = self.LOBE_Q if q is None else q; R = q * L / 2; xc = L - R
        tip, cap = (form[1:] if form is not None else (self.LOBE_TIP, self.LOBE_CAP)) if soft else (math.degrees(math.asin(R / xc)),) * 2
        b = math.radians(cap); J = np.array([xc - R * math.sin(b), R * math.cos(b)])
        k = 0.36 * norm(J); P0, P3 = np.zeros(2), J
        P1, P2 = k * dirv(tip), J - k * np.array([math.cos(b), math.sin(b)])
        t = np.linspace(0.0, 1.0, 24, endpoint=False)[:, None]
        side = (1 - t) ** 3 * P0 + 3 * (1 - t) ** 2 * t * P1 + 3 * (1 - t) * t ** 2 * P2 + t ** 3 * P3
        ph = np.linspace(math.pi / 2 + b, 0.0, 36)
        return np.vstack([side, np.column_stack([xc + R * np.cos(ph), R * np.sin(ph)])])
    def _lobe_tip(self, a, u, A=None, size=1.0):
        """How far along u from atom a its lobe's tip sits (pt): at the nucleus of a bare vertex; on a labelled atom
        behind the letter, so the lobe shows from a bond's margin off the ink on, LOBE_OPEN line widths across there."""
        if not (a.drawn and not a.hidden and a._g): return 0.0
        A = self.P(a) if A is None else A; fm = self._lobe_form(a, ang(u)); H = self._lobe_half(self.LOBE_L * self.st.bond * size, form=fm)
        w = (self.LOBE_OPEN if fm[0] >= self.LOBE_Q else 1.5) * self.st.lw          # (a slim lobe beside a bond: just open)
        d = float(np.interp(w / 2, H[:24, 1], H[:24, 0]))
        return self._exit_hull(a, u, A, self.LOBE_GAP, whole=True) - d
    def _orbital_prims(self, ex):
        st = self.st; _, a, kind, angle, shade, ne, size = ex[:7]; col = (ex[7] if len(ex) > 7 else None) or st.ink
        A = self.P(a); out = []; L = self.LOBE_L * st.bond * size
        grey = "#B8B8B8" if col == st.ink else _mix(col, 0.38)
        lab = a.drawn and not a.hidden and a._g
        def lobe(t, k=1.0, q=None):                        # one lobe along t (a small back lobe has straight sides); a
            u = dirv(t); v = perp(u)                                                         # labelled atom's lobe
            H = self._lobe_half(k * L, q, k == 1.0, self._lobe_form(a, t) if k == 1.0 else None)
            o = self._lobe_tip(a, u, A, size) if lab else 0.0   # comes out from behind its letter, its tip hidden: cut there,
            if lab:                                        # the outline left open (no pinch of two strokes reads as a dot)
                c = self._exit_hull(a, u, A, self.LOBE_GAP, whole=True) - o
                y0 = float(np.interp(c, H[:, 0], H[:, 1])); H = np.vstack([[c, y0], H[H[:, 0] > c + 1e-6]])
            loc = np.vstack([H, (H * [1, -1])[::-1][1:]])
            return A + u * o + loc[:, :1] * u + loc[:, 1:] * v
        if kind == "s":                                    # a sphere, shaded for phase "up" and white for "down" (s + s in
            R = 0.42 * st.bond * size                      # and out of phase look different), its nucleus a dot on a bare
            P = A + R * np.array([dirv(t) for t in np.arange(0.0, 360.0, 5.0)])   # atom (not a bond's dash read as a minus)
            out.append(("orb", P, grey if shade == "up" else "white", col))
            if not lab: out.append(("dot", A, 1.2 * st.lp_r, st.ink, "top"))
        else:                                              # a hybrid lobe has its small back lobe, of the other phase
            back = [] if lab or any(abs(adiff(ang(self.P(n) - A), angle + 180)) < 40 for n in a.mol.nbrs(a)) else \
                   [(angle + 180, shade == "down", 0.36, 0.8)]   # (on a bond it reads as a stray mark, and off a label's
                                                                  # edge it floats free: only a bare vertex shows it)
            lobes = [(angle, shade == "up", 1.0, None)] + ([(angle + 180, shade == "down", 1.0, None)] if kind == "p" else back)
            for t, filled, k, q in lobes:
                P = lobe(t, k, q)
                out.append(("orb", P, grey if filled else "white", col) + (("open",) if lab else ()))
        re_ = 0.8 * st.lp_r; c0, v = self._lobe_pair(A, dirv(angle), self._lobe_tip(a, dirv(angle), A, size), size, re_,
                                                     self._lobe_form(a, angle))
        for k in range(ne):                                # (an orbital's electrons: a little smaller than a lone pair's, so
            out.append(("dot", c0 + (0 if ne == 1 else (2 * k - 1)) * v, re_, st.ink))   # two sit apart)
        return out
    def _lobe_pair(self, A, u, off, size, r=None, form=None):
        """Where a lobe's electron pair sits: its centre in the full part of the lobe (0.68 along from its tip, off along
        u) and the half-spacing vector: dots of radius r, as far apart as an atom's lone pair allows but 0.6 pt of white
        inside the outline."""
        st = self.st; L = self.LOBE_L * st.bond * size; f = 0.68; r = st.lp_r if r is None else r
        H = self._lobe_half(L, form=form); half = float(np.interp(f * L, H[:, 0], H[:, 1]))   # its half-width there
        h = max(r + 0.3, min(st.lp_sep / 2, half - st.lw / 2 - r - 0.6))
        return A + u * (off + f * L), perp(u) * h
    def _highlight(self, g, col):
        """A fill behind a block of text (pt): its ink with 1.5 pt to spare across and 1 pt up and down."""
        ext = [self.F.path(t, x, y, sz, i, bd).get_extents() for t, x, y, sz, i, bd in g["runs"] if t.strip() and not self.F.missing(t)]
        b = (min(e.x0 for e in ext), min(e.y0 for e in ext), max(e.x1 for e in ext), max(e.y1 for e in ext)) if ext else g["box"]
        x0, y0, x1, y1 = b[0] - 1.5, b[1] - 1.0, b[2] + 1.5, b[3] + 1.0
        return ("poly", [np.array([x0, y0]), np.array([x1, y0]), np.array([x1, y1]), np.array([x0, y1])], col)
    def _frame_prims(self, ex, texts):
        """A frame's rectangle (pt) round its items, or an empty box; its label in the top left corner."""
        st = self.st; B = st.bond; _, its, pad, label, size, at = ex
        if its:
            bs = np.array([self._ibox(Group(i) if isinstance(i, Mol) else i) for i in its if not isinstance(i, str)])
            x0, y0, x1, y1 = bs[:, 0].min() - pad * B, bs[:, 1].min() - pad * B, bs[:, 2].max() + pad * B, bs[:, 3].max() + pad * B
            if label: y1 += 1.2 * st.tsize                 # (room for its label above what it frames)
        else:
            cx, cy = at[0] * B, at[1] * B; w, h = size[0] * B, size[1] * B; x0, y0, x1, y1 = cx - w / 2, cy - h / 2, cx + w / 2, cy + h / 2
        C = [np.array(q) for q in ((x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0), (x1, y0))]   # round once more: every
        if label:                                                                                    # corner mitred
            g = self._rich_geo(f"<b>{label}</b>", x0 + 0.25 * B, y1 - 0.2 * B, st.tsize, "left", "top"); texts.append({"runs": g["runs"], "color": st.ink, "box": g["box"]})
        return [("line", C, st.lw, st.ink, None, "butt")]
    def _brackets(self, ex, texts):
        g = self._bracket_geo(ex); texts += g["texts"]; return g["prims"]
    def _bracket_geo(self, ex, boxes=None):
        """Brackets round ex's structures (pt): {'prims', 'texts', 'box'}; the box holds the charge (upper right) and the
        count (a subscript, lower right) set beside them. boxes: structure boxes already measured, by id."""
        st = self.st; p = ex[3] * st.bond; boxes = boxes or {}
        bs = np.array([boxes[id(m)] if id(m) in boxes else self._mol_box(m, curlies=False) for m in ex[1]])
        x0, y0, x1, y1 = bs[:, 0].min() - p, bs[:, 1].min() - p, bs[:, 2].max() + p, bs[:, 3].max() + p; t = 0.18 * st.bond
        out = [("line", [np.array([x0 + t, y0]), np.array([x0, y0]), np.array([x0, y1]), np.array([x0 + t, y1])], st.lw, st.ink, None, "butt"),
               ("line", [np.array([x1 - t, y0]), np.array([x1, y0]), np.array([x1, y1]), np.array([x1 - t, y1])], st.lw, st.ink, None, "butt")]
        texts = []; h = st.lw / 2; box = [x0 - h, y0 - h, x1 + h, y1 + h]
        if ex[2] == "‡":                                   # the double dagger drawn as strokes, a stem and two bars, at the
            H = 0.78 * 0.85 * st.size; W = 0.55 * H; w = self.DAGGER_W * st.lw   # weight of the lettering (a font's ‡ is ornamental)
            cx = x1 + 0.08 * st.bond + W / 2; top = y1 + H / 2
            out += [("line", [np.array([cx, top]), np.array([cx, top - H])], w, st.ink, None, "butt")]
            out += [("line", [np.array([cx - W / 2, top - f * H]), np.array([cx + W / 2, top - f * H])], w, st.ink, None, "butt") for f in (0.26, 0.64)]
            box = [box[0], box[1], max(box[2], cx + W / 2 + w / 2), max(box[3], top)]
        elif ex[2]: texts.append(self._rich_geo(ex[2], x1 + 0.08 * st.bond, y1, 0.85 * st.size, "left", "center"))
        if len(ex) > 4 and ex[4] is not None:              # a count: a subscript beside the lower end (dimer ₂, polymer ₙ)
            texts.append(self._rich_geo(str(ex[4]), x1 + 0.08 * st.bond, y0, 0.85 * st.size, "left", "center"))
        for g in texts:
            b = g["box"]; box = [min(box[0], b[0]), min(box[1], b[1]), max(box[2], b[2]), max(box[3], b[3])]
        return {"prims": out, "texts": [{"runs": g["runs"], "color": st.ink, "box": g["box"]} for g in texts], "box": tuple(box)}
    def _arrow_branch_sides(self, it):
        """A reaction arrow's branches by side: {'above': [...], 'below': [...]}; a branch without side= goes where the
        arrow's own text leaves room (above, unless only the space below is free)."""
        busy_a = bool(self._parts(it.above)) or it.step is not None; busy_b = bool(self._parts(it.below))
        out = {}
        for b in it.branches: out.setdefault(b[3] or ("below" if busy_a and not busy_b else "above"), []).append(b)
        return out
    def _curve_arrow(self, pts, scale=1.0, log=True):
        """Shaft and head of a curved black arrow drawn along the polyline pts (tail -> tip). The shaft is cut where its
        chord to the tip equals the head's notch depth, then runs straight along that chord a little way into the head
        (as straight arrows do), so nothing of it shows past the head; the head's axis is that chord, so the head sits
        symmetrically on the curve."""
        st = self.st; h = st.head_len * scale; k = st.head_notch; P = np.asarray(pts, float); tip = P[-1]
        dist = np.hypot(*(P - tip).T); far = np.flatnonzero(dist >= k * h)
        if not len(far):                                        # too short to curve: a straight arrow
            d = unit(tip - P[0]); shaft = [P[0], tip - d * 0.8 * h]
        else:
            j = int(far[-1]); a, b = P[j], P[min(j + 1, len(P) - 1)]; lo, hi = 0.0, 1.0
            for _ in range(40):                                 # chord |a + t(b - a) - tip| = k h, t in [0, 1]
                mid = (lo + hi) / 2
                if norm(a + (b - a) * mid - tip) >= k * h: lo = mid
                else: hi = mid
            cut = a + (b - a) * lo; d = unit(tip - cut)
            shaft = list(P[:j + 1]) + [cut, tip - d * 0.8 * h]
        if log: self._arrow_log.append({"shaft": np.array(shaft), "tip": tip, "d": d, "h": h, "w": st.head_w * scale, "k": k,
                                        "lw": st.lw, "curved": len(far) > 0, "half": 0})
        return [("line", shaft, st.lw, st.ink, None, "butt"), self._head(tip, d, scale)]
    def _cycle_prims(self, cy, texts, place=False):
        """Arcs, branches, step names and feeds of a cycle (pt). place=True (before a build) moves structures given as
        branch or feed labels to where they belong."""
        st = self.st; C = cy["c"] * st.bond; R = cy["r"] * st.bond; sgn = -1 if cy["cw"] else 1; out = []
        boxes = [self._ibox(it) for it, a in cy["its"]]; angs = [a for it, a in cy["its"]]; g = cy["gap"] * st.bond
        hulls = []                                    # a structure station's outline: the convex hull of its letters and
        for it, a in cy["its"]:                       # bond ends, so an arc runs on to within the gap of what is drawn
            ms = it.mols if isinstance(it, Group) else [it] if isinstance(it, Mol) else []   # (not of its box's empty
            P = [q for m in ms for at in m.atoms for x in (at._g or []) for q in ((x["box"][0], x["box"][1]),   # corners)
                 (x["box"][2], x["box"][1]), (x["box"][2], x["box"][3]), (x["box"][0], x["box"][3]))]          # but not
            P += [tuple(self.P(at)) for m in ms for at in m.atoms if not at.hidden and not at.drawn]          # into the
            hulls.append(_hull(np.array(P, float)) if len(P) >= 3 else None)                                 # gaps
        def outside(t, k):
            p = C + dirv(t) * R
            if hulls[k] is None or len(hulls[k]) < 3: return box_dist(p, boxes[k]) > g
            return _hull_dist(p, hulls[k]) > g
        arcs, heads = [], []
        for i in range(len(angs)):
            j = (i + 1) % len(angs); a0, a1 = angs[i], angs[j]
            span = (a1 - a0) * sgn % 360 or 360; t0 = 0.0
            while t0 < span and not outside(a0 + sgn * t0, i): t0 += 0.5
            t1 = span
            while t1 > t0 and not outside(a0 + sgn * t1, j): t1 -= 0.5
            ts = np.linspace(a0 + sgn * t0, a0 + sgn * t1, 80); pts = C + np.array([dirv(t) for t in ts]) * R
            ar = self._curve_arrow(pts); out += ar
            arcs.append(ts); heads.append(np.asarray(ar[1][1], float))
        def geom(i, kind, side):
            """Shape, label and clearance helpers for a branch on arc i."""
            ts = arcs[i]; Hh = _hull(heads[i])                # the arc's own head: a branch keeps 2.5 pt clear of it
            def shape(at, L):
                t = ts[int(min(max(at, 0.0), 1.0) * (len(ts) - 1))]; P0 = C + dirv(t) * R
                T = unit(perp(dirv(t)) * sgn); r = dirv(t) * (1 if side == "out" else -1)
                end = P0 + r * L + (T if kind == "out" else -T) * 0.9 * L
                if kind == "in": pts = bezier([end, end - r * 0.55 * L, P0 - T * 0.55 * L, P0], 60)
                else: pts = bezier([P0, P0 + T * 0.55 * L, end - r * 0.55 * L, end], 60)   # the mirror of an incoming branch
                return P0, r, end, pts
            def clear_of_head(pts): return min(_hull_dist(p, Hh) for p in pts) >= st.lw / 2 + 2.5
            return ts, shape, clear_of_head
        mb = {}
        def lab_geo(label, end, r):
            ha = "left" if r[0] > 0.3 else "right" if r[0] < -0.3 else "center"
            va = "bottom" if r[1] > 0.3 else "top" if r[1] < -0.3 else "center"
            q = end + r * 0.25 * st.bond
            if isinstance(label, Mol):                       # a structure: its box set where a text block would go
                if id(label) not in mb: mb[id(label)] = self._mol_box(label, curlies=False)
                b0 = mb[id(label)]; w, h = b0[2] - b0[0], b0[3] - b0[1]
                x0 = {"left": q[0], "right": q[0] - w, "center": q[0] - w / 2}[ha]; y0 = {"bottom": q[1], "top": q[1] - h, "center": q[1] - h / 2}[va]
                bx = (x0, y0, x0 + w, y0 + h); return {"runs": [], "box": bx, "mol": (label, x0 - b0[0], y0 - b0[1])}, bx
            if va == "top": q = q - (0, self._ink_shift(label, "below", st.tsize))      # a superscript leading the label
            elif va == "bottom": q = q + (0, self._ink_shift(label, "above", st.tsize))  # keeps the line's distance
            g2 = self._rich_geo(label, q[0], q[1], st.tsize, ha, va); return g2, g2["box"]
        def put(g2):                                          # a structure label moves to its place (before a build)
            if "mol" in g2 and place:
                m, dx, dy = g2["mol"]
                if abs(dx) > 1e-9 or abs(dy) > 1e-9: m.move(dx / st.bond, dy / st.bond); mb[id(m)] = g2["box"]
            elif g2["runs"]: texts.append({"runs": g2["runs"], "color": st.ink, "box": g2["box"]})
        def sep(gb, b): return max(gb[0] - b[2], b[0] - gb[2], gb[1] - b[3], b[1] - gb[3])
        brs = []                                              # where each branch joins: reagents in before by-products out
        for n_br, (i, label, kind, at, reach, side) in enumerate(cy["branches"]):
            auto = at is None
            if auto:
                same = [b for b in cy["branches"] if b[0] == i and b[3] is None and b[5] == side]
                ins = [b for b in same if b[2] == "in"]; outs = [b for b in same if b[2] == "out"]
                grp = ins if kind == "in" else outs; k = next(j for j, b in enumerate(grp) if b is cy["branches"][n_br])
                if not (ins and outs) and len(grp) == 1: at = 0.5
                elif kind == "in": at = 0.3 if len(grp) == 1 else 0.2 + 0.25 * k / (len(grp) - 1)
                else: at = 0.7 if len(grp) == 1 else 0.55 + 0.25 * k / (len(grp) - 1)
            ts = arcs[i]; Larc = abs(math.radians(ts[-1] - ts[0])) * R
            at = min(max(at, min(0.45, 0.35 * st.bond / Larc)), max(0.55, 1.0 - (1.6 * st.head_len + 0.35 * st.bond) / Larc))
            brs.append([i, label, kind, at, reach, side, auto])
        need = {}                                             # one reach for every branch on an arc: the longest any label needs
        for i, label, kind, at, reach, side, auto in brs:
            ts, shape, coh = geom(i, kind, side); L = reach * st.bond
            for _ in range(60):
                P0, r, end, pts = shape(at, L); g2, gb = lab_geo(label, end, r)
                if all(sep(gb, b) > 0.3 * st.bond for b in boxes): break
                L += 0.2 * st.bond
            need[(i, side)] = max(need.get((i, side), 0.0), L)
        for b in brs:                                         # a by-product leaving late starts earlier, clear of the arc's head,
            i, label, kind, at, reach, side, auto = b        # and an automatic reagent on the same arc mirrors it
            if kind != "out": continue
            ts, shape, coh = geom(i, kind, side); a2 = at
            while not coh(shape(a2, need[(i, side)])[3]) and a2 > 0.35: a2 = max(0.35, a2 - 0.02)
            b[3] = a2
            if auto:
                for b2 in brs:
                    if b2[0] == i and b2[5] == side and b2[2] == "in" and b2[6] and abs(b2[3] - (1 - at)) < 1e-9: b2[3] = 1 - a2
        placed_lbl, placed_paths, placed_junc = [], [], []
        for i, label, kind, at, reach, side, auto in brs:
            ts, shape, coh = geom(i, kind, side)
            def attempt(a_req):
                L = need[(i, side)]
                for _ in range(60):                   # lengthen the branch until its label is clear: arrow and label stay together
                    a2 = a_req
                    while True:                       # and start it earlier on the arc until it leaves the arc's head alone
                        P0, r, end, pts = shape(a2, L)
                        if coh(pts) or a2 <= 0.35: break
                        a2 = max(0.35, a2 - 0.02)
                    g2, gb = lab_geo(label, end, r)
                    if all(sep(gb, b) > 0.3 * st.bond for b in boxes) and all(sep(gb, b) > 1.2 * st.bond for b in placed_lbl): break
                    L += 0.2 * st.bond
                return pts, g2, gb, P0
            def crosses(pts, gb, P0):                 # 2: its line through another branch's label or its label on another's
                pad = st.lw / 2 + 1.5                 # line; 1: the same junction as another branch on this arc; 0: clear
                for b in placed_lbl:
                    if (((pts >= (b[0] - pad, b[1] - pad)) & (pts <= (b[2] + pad, b[3] + pad))).all(1)).any(): return 2
                for Pq in placed_paths:
                    if (((Pq >= (gb[0] - pad, gb[1] - pad)) & (Pq <= (gb[2] + pad, gb[3] + pad))).all(1)).any(): return 2
                return 1 if any(j == i and norm(P0 - Pj) < 0.8 * st.bond for j, Pj in placed_junc) else 0
            best = None                               # the requested place, or the nearest one that crosses nothing (failing
            for a_try in [at] + [at + sg * 0.05 * k for k in range(1, 12) for sg in (1, -1)]:   # that, shares a junction at worst)
                if not 0.15 <= a_try <= 0.9: continue
                res = attempt(a_try); sev = crosses(res[0], res[2], res[3])
                if best is None or sev < best[0]: best = (sev, res)
                if sev == 0: break
            pts, g2, gb, P0 = best[1]
            placed_paths.append(np.asarray(pts)); placed_junc.append((i, P0))
            if kind == "in": out.append(("line", list(pts), st.lw, st.ink, None, "butt"))
            else: out += self._curve_arrow(pts)
            placed_lbl.append(gb); put(g2)
        def busy(gb):                                         # what a label set here would run into
            pad = 0.3 * st.bond; n = sum(1 for b in boxes + placed_lbl if sep(gb, b) < pad)
            for Pq in placed_paths + [C + np.array([dirv(t) for t in ts_]) * R for ts_ in arcs]:
                if (((Pq >= (gb[0] - 2, gb[1] - 2)) & (Pq <= (gb[2] + 2, gb[3] + 2))).all(1)).any(): n += 1
            return n
        for i, label, at, side in cy.get("steps", []):      # step names beside their arcs, where nothing else stands
            ts = arcs[i]; best = None                         # (outside first, then inside, unless a side is given)
            for sd in ([side] if side else ["out", "in"]):
                for a_ in ([at] if at is not None else [0.5, 0.35, 0.65, 0.25, 0.75, 0.15, 0.85]):
                    t = ts[int(min(max(a_, 0.0), 1.0) * (len(ts) - 1))]; r = dirv(t) * (1 if sd == "out" else -1)
                    g2, gb = lab_geo(label, C + dirv(t) * R + r * st.lw / 2, r); n = busy(gb)
                    if best is None or n < best[0]: best = (n, g2, gb)
                    if n == 0: break
                if best[0] == 0: break
            put(best[1]); placed_lbl.append(best[2])
        for i, label, reach, note in cy.get("feeds", []):    # straight arrows into stations from outside
            b = boxes[i]; c = np.array([(b[0] + b[2]) / 2, (b[1] + b[3]) / 2]); u = unit(c - C) if norm(c - C) > 1e-6 else dirv(angs[i])
            ex_ = min((b[2] - b[0]) / 2 / abs(u[0]) if abs(u[0]) > 1e-9 else 1e9, (b[3] - b[1]) / 2 / abs(u[1]) if abs(u[1]) > 1e-9 else 1e9)
            tip = c + u * (ex_ + 0.2 * st.bond); tail = tip + u * (reach or 1.4) * st.bond
            out += self._curve_arrow([tail, tip]); put(lab_geo(label, tail, u)[0])
            gl = lab_geo(label, tail, u)[1]; placed_lbl.append(gl)
            if note:                                          # its note beside it, on the right (or above) unless that
                n0 = perp(u); n0 = n0 if (n0[0] > 1e-9 or (abs(n0[0]) <= 1e-9 and n0[1] > 0)) else -n0; best = None   # side
                for n in (n0, -n0):                                                                                 # is taken
                    q = (tail + tip) / 2 + n * 0.35 * st.bond
                    g3 = self._rich_geo(note, q[0], q[1], st.tsize, "left" if n[0] > 0.3 else "right" if n[0] < -0.3 else "center",
                                        "bottom" if n[1] > 0.3 else "top" if n[1] < -0.3 else "center")
                    k = busy(g3["box"])
                    if best is None or k < best[0]: best = (k, g3)
                    if k == 0: break
                texts.append({"runs": best[1]["runs"], "color": st.ink, "box": best[1]["box"]}); placed_lbl.append(best[1]["box"])
        return out
    def _draw(self, ax, prims, texts, grid=None):
        """grid=(x0, y0, px per pt across, px per pt up) for a bitmap: each run of text starts on a whole pixel, so
        the same letter is drawn with the same pixels wherever it stands (two L's of one complex, the H of every H⁺);
        a shift of half a pixel at most, so nothing else moves."""
        st = self.st
        if grid is not None:
            gx, gy, sx, sy = grid; snap = lambda x, y: (gx + round((x - gx) * sx) / sx, gy + round((y - gy) * sy) / sy)
            texts = [dict(t, runs=[(s_, *snap(x, y), sz, it, bd) for s_, x, y, sz, it, bd in t["runs"]]) for t in texts]
        for pr in prims:
            if pr[0] == "line":
                P = np.array(pr[1]); kw = {"dashes": pr[4]} if pr[4] else {}
                cap = pr[5] if len(pr) > 5 else "round"
                ax.plot(P[:, 0], P[:, 1], lw=pr[2], color=pr[3], solid_capstyle=cap, dash_capstyle=cap,   # (unsnapped: a
                        solid_joinstyle="miter" if cap == "butt" else "round", snap=False, **kw)   # level or upright line
                                                                    # snapped to the pixel grid leaves a step at its joins)
            elif pr[0] == "poly":
                if len(pr) > 3 and isinstance(pr[3], dict):   # a perspective ring's outline: a hairline of its own colour
                    ax.add_patch(Polygon(np.array(pr[1]), closed=True, fc=pr[2], ec=pr[2], lw=self.SEAM, joinstyle="miter"))
                else: ax.add_patch(Polygon(np.array(pr[1]), closed=True, fc=pr[2], ec="none", lw=0))
            elif pr[0] == "circle": ax.add_patch(Circle(pr[1], pr[2], fc="white" if not pr[4] else pr[3], ec=pr[3],
                                                        lw=pr[5] if len(pr) > 5 else 0.8 * st.lw))
            elif pr[0] == "dot": ax.add_patch(Circle(pr[1], pr[2], fc=pr[3], ec="none",
                                                     zorder={("top",): 2.5, ("bottom",): 0.5}.get(tuple(pr[4:]), 1)))
            elif pr[0] == "orb" and pr[4:] == ("open",):    # a lobe from behind a label: filled closed, its outline open,
                ax.add_patch(Polygon(np.array(pr[1]), closed=True, fc=pr[2], ec="none", lw=0))      # both ending square a
                ax.add_patch(Polygon(np.array(pr[1]), closed=False, fill=False, ec=pr[3], lw=st.lw, joinstyle="round", capstyle="butt"))   # hair off the letter
            elif pr[0] == "orb": ax.add_patch(Polygon(np.array(pr[1]), closed=True, fc=pr[2], ec=pr[3], lw=st.lw, joinstyle="round"))
        for t in texts:
            if t.get("halo"):
                for s, x, y, sz, it, bd in t["runs"]:
                    if s.strip(): ax.add_patch(PathPatch(self.F.path(s, x, y, sz, it, bd), fc="white", ec="white", lw=2.2 * st.lw, joinstyle="round"))
            for s, x, y, sz, it, bd in t["runs"]:
                if s.strip(): ax.add_patch(PathPatch(self.F.path(s, x, y, sz, it, bd), fc=getattr(s, "c", None) or t["color"], ec="none", lw=0,
                                                     snap=False))   # (Agg would round an H's or L's square stems to 3 or 4 px)
    def bbox(self):
        prims, texts = self._build(); pts = []
        for pr in prims:
            if pr[0] in ("line", "poly", "orb"): pts += [np.array(p) for p in pr[1]]
            else: pts += [pr[1] - pr[2], pr[1] + pr[2]]
        for t in texts:
            for s, x, y, sz, it, bd in t["runs"]:
                if s.strip():
                    e = self.F.path(s, x, y, sz, it, bd).get_extents(); pts += [np.array([e.x0, e.y0]), np.array([e.x1, e.y1])]
        P = np.array(pts); return (*P.min(0), *P.max(0)), prims, texts
    PAGES = {"A4": (595.276, 841.89), "letter": (612.0, 792.0), "A3": (841.89, 1190.551)}
    def save_pdf(self, path, page="A4", margin_cm=2.0, landscape=False, numbers=True):
        """Print-ready PDF: the figure at its true size (the ACS scale, 1 pt = 1 pt) centred on A4 pages ('letter' and
        'A3' too) inside margin_cm margins, vector throughout, text as outlines. A figure taller than a page is split
        between its rows, never through a structure, arrow or caption, over as many pages as it needs, numbered at
        the foot ('2 / 3'). A figure wider than the page is scaled down to fit, the same on every page (a note says so
        when its labels would print below 7.5 pt); landscape=True turns the pages instead. Lay long mechanisms out
        with flow() so their rows suit the page (its default 16 cm fits A4 portrait). For a tight PDF of just the
        figure (to place in a document), use save('figure.pdf'). Returns the number of pages."""
        from matplotlib.backends.backend_pdf import PdfPages
        (x0, y0, x1, y1), prims, texts = self.bbox()
        objs = []                                              # (y_lo, y_hi, kind, object) of everything drawn
        for pr in prims:
            if pr[0] in ("line", "poly", "orb"): P = np.asarray(pr[1], float); lo, hi = P[:, 1].min(), P[:, 1].max()
            else: lo, hi = pr[1][1] - pr[2], pr[1][1] + pr[2]
            w = pr[2] / 2 if pr[0] == "line" else (self.st.lw / 2 if pr[0] in ("circle", "orb") else 0.0)
            objs.append((lo - w, hi + w, "p", pr))
        for t in texts:
            ys = [e for s, x, y, sz, it, bd in t["runs"] if s.strip() for e in
                  [self.F.path(s, x, y, sz, it, bd).get_extents()] for e in (e.y0, e.y1)]
            if ys: objs.append((min(ys), max(ys), "t", t))
        blocks = []                                            # rows: stretches of height nothing crosses
        for lo, hi, _, _ in sorted(objs, key=lambda o: -o[1]):
            if blocks and hi > blocks[-1][0] - 1.0: blocks[-1][0] = min(blocks[-1][0], lo)
            else: blocks.append([lo, hi])
        pw0, ph0 = self.PAGES[page]; m = margin_cm / 2.54 * 72; W = x1 - x0
        pw, ph = (ph0, pw0) if landscape else (pw0, ph0)
        foot = 0.6 * m if numbers else 0.0
        s = min(1.0, (pw - 2 * m) / W, (ph - 2 * m) / max(hi - lo for lo, hi in blocks))
        if s * self.st.size < 7.5:
            print(f"bondline: {path} is scaled to {s:.0%} to fit the page, so its labels print at {s * self.st.size:.1f} pt: "
                  "wrap it with flow(width_cm=...) or use landscape=True for a larger print")
        pages, cur = [], None                                  # fill pages from the top, whole rows at a time
        for lo, hi in blocks:
            if cur and (cur[1] - lo) * s <= ph - 2 * m: cur[0] = lo
            else: cur = [lo, hi]; pages.append(cur)
        n = len(pages); xm = (x0 + x1) / 2
        with PdfPages(path, metadata={"CreationDate": None, "Producer": "bondline", "Creator": "bondline"}) as pdf:
            for i, (lo, hi) in enumerate(pages):
                fig = plt.figure(figsize=(pw / 72, ph / 72)); ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
                top = hi + m / s if n > 1 else (lo + hi) / 2 + ph / (2 * s)    # one page: centred; several: from the top
                ax.set_xlim(xm - pw / (2 * s), xm + pw / (2 * s)); ax.set_ylim(top - ph / s, top)
                keep = [o for o in objs if o[1] >= lo - 0.5 and o[0] <= hi + 0.5]
                self._draw(ax, [o[3] for o in keep if o[2] == "p"], [o[3] for o in keep if o[2] == "t"])
                if numbers and n > 1:
                    g = self._rich_geo(f"{i + 1} / {n}", xm, top - (ph - foot) / s, 8.0 / s)
                    self._draw(ax, [], [{"runs": g["runs"], "color": self.st.ink}])
                pdf.savefig(fig); plt.close(fig)
        return n
    DPI = 1200                    # default PNG resolution at the figure's true size: the line-art standard of ACS
    MAX_PX = 12000                # journals (the strictest major publisher); a huge figure is capped at this long edge
    def save(self, path, dpi=None, pad=4.0, debug=False, px=None):
        """Render at the fixed scale (1 data unit = 1 pt). Returns (width_pt, height_pt). Never rescale afterwards.
        A PNG comes out at DPI (1200) pixels per inch of its true size, sharp in print as well as on screen, with the
        long edge capped at MAX_PX; dpi= sets another resolution, px= an exact long edge in pixels. SVG and PDF are
        vector and ignore both."""
        (x0, y0, x1, y1), prims, texts = self.bbox()
        x0 -= pad; y0 -= pad; x1 += pad; y1 += pad; w, h = x1 - x0, y1 - y0
        if px is not None: dpi = px * 72.0 / max(w, h) * (1 + 1e-9)        # a hair over, so the edge is exactly px
        elif dpi is None: dpi = min(self.DPI, self.MAX_PX * 72.0 / max(w, h))
        fig = plt.figure(figsize=(w / 72, h / 72), dpi=dpi); ax = fig.add_axes([0, 0, 1, 1])
        ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.axis("off")
        ext = os.path.splitext(str(path))[1].lower()
        self._draw(ax, prims, texts, grid=(x0, y0, dpi / 72.0, dpi / 72.0) if ext not in (".svg", ".pdf") else None)
        if debug:
            for mi, m in enumerate(self.mols):
                for i, a in enumerate(m.atoms):
                    p = self.P(a); ax.text(p[0] + 2, p[1] + 2, f"{mi}.{a.name if a.name is not None else i}", fontsize=4.5, color="#1565C0")
                                                                # vector files carry no date and fixed ids: re-saves are identical
        kw = {".svg": {"metadata": {"Date": None}}, ".pdf": {"metadata": {"CreationDate": None}}}.get(ext, {})
        with matplotlib.rc_context({"svg.hashsalt": "bondline"}): fig.savefig(path, dpi=dpi, **kw)
        plt.close(fig); self._last = (path, dpi, x0, y1)
        return w, h

    # ---------- QA
    def check(self, quiet=False):
        """List drawing problems: bond lengths and angles, clashes and crossings, curly-arrow placement, label and lone-pair
        consistency, stereo depiction and valences. Returns the issues (notes for deliberate geometry are printed apart)."""
        st = self.st; prims_all, texts_all = self._build(); issues = []; notes = []
        react = set()
        for c in self.curlies:
            for x in (c.src, c.dst):
                for a in ([x.atom] if isinstance(x, LP) else [x] if isinstance(x, Atom) else [x.a, x.b] if isinstance(x, Mid)
                          else list(x) if isinstance(x, (tuple, list)) else []):
                    if isinstance(a, Atom): react.add(a.uid); react.update(n.uid for n in a.mol.nbrs(a))
        near_react = lambda a: a.uid in react or any(n.uid in react for n in a.mol.nbrs(a))
        for mi, m in enumerate(self.mols):
            for b in m.bonds:
                L = norm(b.b.xy - b.a.xy)
                if not b.free and b.style not in ("dash", "dot") and abs(L - 1) > 0.03:
                    seen, todo = {b.a.uid}, [b.a]                  # is it the only link between two fragments?
                    while todo:
                        x = todo.pop()
                        for y in m.nbrs(x):
                            if y.uid not in seen and not (x is b.a and y is b.b) and not (x is b.b and y is b.a):
                                seen.add(y.uid); todo.append(y)
                    hint = ": to connect two fragments use join(), which lays the new bond out at 1 BU" if b.b.uid not in seen else ""
                    if getattr(b, "_smiles", False) and hint:       # (from_smiles's own layout: join() does not apply)
                        hint = ": RDKit's layout of a crowded centre stretched it; pin the atoms round it with xy="
                    issues.append(f"mol{mi}: bond {b.a}-{b.b} length {L:.2f} BU{hint} (flag free=True if intended)")
        allb = [(m, b) for m in self.mols for b in m.bonds]
        vis = []
        for m, b in allb:
            if b.style != "plain" or b.order != 1 or b.free or (b.a.drawn == b.b.drawn): continue
            p, q = b._prims[0][1]; U, L = (b.b, b.a) if b.a.drawn else (b.a, b.b)
            vis.append((norm(q - p) + self._centred_trim(U, b, unit(self.P(L) - self.P(U))), b))   # (a C=O's near line: not a stub)
        for v, b in vis:
            if v < 0.40 * st.bond:
                L = b.a if b.a.drawn else b.b
                sup_only = b._limited and all(s for x in (L._g or []) for h, s in zip(x["hulls"], x["sup"])
                                              if h is not getattr(L, "_anchor_hull", None))
                if sup_only:                                   # R² with its bond leaving at 0-57 deg: starts past the ²
                    notes.append(f"bond {b.a}-{b.b} shows {v:.1f} pt: it leaves {L} past the superscript, as ChemDraw draws it; "
                                 "for a longer bond put the numbered group where its bond leaves sideways, up or down")
                    continue
                issues.append(f"bond {b.a}-{b.b} shows only {v:.1f} pt between labels (stub)" + (
                    "; the label's other characters stand in its way: write the label to read away from the bond" if b._limited else ""))
        from matplotlib.path import Path as _Pth          # nothing may run through a black arrowhead (cycle arcs, branches)
        for A in getattr(self, "_arrow_log", []):
            if "tip" not in A: continue
            tip, d, h, w, kk = A["tip"], A["d"], A["h"], A["w"], A["k"]; n = perp(d)
            H = _Pth(np.array([tip, tip - d * h + n * w, tip - d * kk * h, tip - d * h - n * w]))
            own = np.asarray(A["shaft"], float)
            for pr in prims_all:
                if pr[0] != "line": continue
                P = np.asarray(pr[1], float)
                if P.shape == own.shape and np.allclose(P, own): continue
                seg = np.diff(P, axis=0); Lg = np.hypot(*seg.T)
                Q = np.vstack([P[k] + seg[k] * t for k in range(len(seg)) for t in np.linspace(0, 1, max(2, int(Lg[k] / 0.3) + 1))] or [P])
                if H.contains_points(Q, radius=-(pr[2] / 2 + 0.6)).any() or H.contains_points(Q, radius=(pr[2] / 2 + 0.6)).any():
                    issues.append(f"a line runs into the arrowhead at ({tip[0]:.0f}, {tip[1]:.0f}) pt"); break
        bond_prims = {id(pr) for m, b in allb for pr in b._prims}      # arrows and cycle branches never run through text
        tlines = [pr for pr in prims_all if pr[0] == "line" and id(pr) not in bond_prims]
        if tlines:
            for t in texts_all:
                if t.get("atom") or not t["runs"]: continue
                hs = []
                for s_, x, y, sz, it, bd in t["runs"]:
                    for ch in s_:
                        H = self.F.hull(ch, it, bd)
                        if H is not None and len(H) >= 3: hs.append(H * sz + (x, y))
                        x += self.F.adv(ch, sz, it, bd)
                if not hs: continue
                allh = np.vstack(hs); tb0, tb1 = allh.min(0) - 2, allh.max(0) + 2
                E = np.vstack([np.hstack([h, np.roll(h, -1, axis=0)]) for h in hs]); s0 = np.cumsum([0] + [len(h) for h in hs[:-1]]).astype(int)
                for pr in tlines:
                    P = np.asarray(pr[1], float)
                    if (P.max(0) < tb0).any() or (P.min(0) > tb1).any(): continue
                    seg = np.diff(P, axis=0); Lg = np.hypot(*seg.T)
                    Q = np.vstack([P[k] + seg[k] * tt for k in range(len(seg)) for tt in np.linspace(0, 1, max(2, int(Lg[k] / 0.3) + 1))])
                    Q = Q[((Q >= tb0) & (Q <= tb1)).all(1)]
                    if len(Q) and _polys_dist(Q, E, s0).min() < pr[2] / 2 + 0.5:
                        whole = "".join(r[0] for r in t["runs"])
                        issues.append(f"a line runs through the text '{whole}'"); break
        for m, b in allb:                                   # white between a bond and the letters of its own labels
            if not b._prims: continue
            pieces = [(np.asarray(pr[1], float), pr[2]) for pr in b._prims if pr[0] == "line"]
            for L in (b.a, b.b):
                if not L._g or L.hidden or L.clear is not None: continue
                one = getattr(L, "_anchor1", False); ah = getattr(L, "_anchor_hull", None)
                hs = [h for x in L._g for h in x["hulls"] if not (one and h is ah)]     # the bonded letter: its circle
                if not hs: continue
                E = np.vstack([np.hstack([h, np.roll(h, -1, axis=0)]) for h in hs])
                s0 = np.cumsum([0] + [len(h) for h in hs[:-1]]).astype(int); worst = 99.0
                for P, w in pieces:
                    for i in range(len(P) - 1):
                        p0, q0 = P[i], P[i + 1]; Lp = norm(q0 - p0)
                        if Lp < 1e-6: continue
                        nv = perp((q0 - p0) / Lp) * w / 2; t = np.linspace(0, 1, max(2, int(Lp / 0.25) + 1))[:, None]
                        Q = np.vstack([p0 + nv + (q0 - p0) * t, p0 - nv + (q0 - p0) * t, p0 + np.linspace(-1, 1, 5)[:, None] * nv,
                                       q0 + np.linspace(-1, 1, 5)[:, None] * nv])          # the outline of the drawn line
                        worst = min(worst, float(_polys_dist(Q, E, s0).min()))
                if worst < 0.8 * st.bond / 14.4:
                    issues.append(f"bond {b.a}-{b.b} comes {worst:.1f} pt from the letters of {L} (keep 0.8 pt of white)")
        for m, b in allb:                                   # a dative arrow needs a shaft; an eta line must be readable
            if b.style == "dative" and b.a.drawn and b.b.drawn and not (b.a.hidden or b.b.hidden) and b._vis is not None:
                v = norm(np.asarray(b._vis[1]) - np.asarray(b._vis[0]))
                if v < self._dative_need() - 0.5:
                    issues.append(f"dative bond {b.a}-{b.b} shows {v:.1f} pt, too little for its arrowhead and a shaft; "
                                  f"lengthen it (length=1.6, free=True)")
            e = b.b if getattr(b.b, "_eta", None) else b.a if getattr(b.a, "_eta", None) else None
            if e is not None and not getattr(e, "_persp3d", False) and b._prims and b._prims[0][0] == "line":
                P = np.asarray(b._prims[0][1], float); v = norm(P[-1] - P[0])
                if v < 0.45 * st.bond:
                    M = b.a if e is b.b else b.b
                    issues.append(f"the eta line from {M} shows only {v:.1f} pt before it meets its pi system; put a ring's "
                                  f"centre about 1.75 BU from the metal (as sandwich() does)")
        for m, b in allb:                                   # label-to-label bonds (O-H, X-Pd, N=N) must stay visible
            if b.style != "plain" or not (b.a.drawn and b.b.drawn) or b.a.hidden or b.b.hidden or not b._prims: continue
            p, q = b._prims[0][1][0], b._prims[0][1][-1]
            if norm(np.array(q) - np.array(p)) < 0.18 * st.bond:
                issues.append(f"bond {b.a}-{b.b} between two labels shows only {norm(np.array(q) - np.array(p)):.1f} pt; "
                              "lengthen it (length=1.3, free=True) or write one condensed label")
        def own_pi(b1, b2):                                 # a metal's line to the centre of a pi system crosses that
            for x, y in ((b1, b2), (b2, b1)):               # system's own bonds by design (flat or in perspective), and a
                c = x.b if getattr(x.b, "_eta", None) else x.a if getattr(x.a, "_eta", None) else None   # floating
                if c is not None and y.a in c._eta and y.b in c._eta: return True                      # substituent its
                f = x.b if getattr(x.b, "_float", None) else x.a if getattr(x.a, "_float", None) else None   # ring edge
                if f is not None and {y.a.uid, y.b.uid} == {f._float[0].uid, f._float[1].uid}: return True
            return False
        for (m1, b1), (m2, b2) in itertools.combinations(allb, 2):
            if {b1.a.uid, b1.b.uid} & {b2.a.uid, b2.b.uid} or own_pi(b1, b2): continue
            if _seg_x(self.P(b1.a), self.P(b1.b), self.P(b2.a), self.P(b2.b)):
                if any(g[2] is b2 for g in getattr(b1, "_gaps", [])) or any(g[2] is b1 for g in getattr(b2, "_gaps", [])):
                    continue                                # drawn in perspective: the back bond is broken
                issues.append(f"bonds cross: {b1.a}-{b1.b} and {b2.a}-{b2.b}")
        stubs = [(norm(np.array(b._prims[0][1][-1]) - np.array(b._prims[0][1][0])), b) for m, b in allb
                 if b.a.drawn and b.b.drawn and not (b.a.hidden or b.b.hidden) and b.style == "plain" and not b.free
                 and b._prims and b._prims[0][0] == "line" and abs(norm(b.b.xy - b.a.xy) - 1) < 0.03]
        for mm in self.mols:                                # within one structure (S-O and Si-O in two molecules differ
            ss = [(s_, b) for s_, b in stubs if b.a.mol is mm]   # by their letters, not by drawing)
            if len(ss) < 2: continue
            med = float(np.median([s_ for s_, b in ss]))
            for s_, b in ss:
                if abs(s_ - med) > 0.4 and not b._limited:
                    issues.append(f"bond {b.a}-{b.b} between labels is {s_:.1f} pt, others {med:.1f} pt (should match)")
        for m in self.mols:
            for a in m.atoms:
                for l, c0, d1, d2 in a._lp:
                    if abs(norm(d1 - self.P(a)) - norm(d2 - self.P(a))) > 0.2: issues.append(f"lone pair on {a} is lopsided")
        segs_drawn = [(b, np.array(pr[1][i]), np.array(pr[1][i + 1])) for m, b in allb for pr in b._prims if pr[0] == "line"
                      for i in range(len(pr[1]) - 1) if b.style in ("plain", "bold", "dative")] + \
                     [(b, np.asarray(pr[3]["spine"][0], float), np.asarray(pr[3]["spine"][1], float)) for m, b in allb
                      for pr in b._prims if pr[0] == "poly" and len(pr) > 3 and isinstance(pr[3], dict)]  # perspective outlines
        for (b1, p1, q1), (b2, p2, q2) in itertools.combinations(segs_drawn, 2):
            if b1 is b2 or {b1.a.uid, b1.b.uid} & {b2.a.uid, b2.b.uid} or own_pi(b1, b2): continue
            if _seg_x(p1, q1, p2, q2): issues.append(f"bond {b1.a}-{b1.b} is drawn through bond {b2.a}-{b2.b}")
        boxes = [(g["box"], a) for m in self.mols for a in m.atoms for g in (a._g or [])]
        for (b1, a1), (b2, a2) in itertools.combinations(boxes, 2):
            if a1 is not a2 and b1[0] < b2[2] + 0.5 and b2[0] < b1[2] + 0.5 and b1[1] < b2[3] + 0.5 and b2[1] < b1[3] + 0.5:
                issues.append(f"labels overlap: {a1} / {a2}")
        def hull_gap(H1, H2):
            return min(min(_hull_dist(p, H2) for p in H1), min(_hull_dist(p, H1) for p in H2))
        glyphs = [(h, a) for m in self.mols for a in m.atoms for g in (a._g or []) for h in g["hulls"]]
        for (h1, a1), (h2, a2) in itertools.combinations(glyphs, 2):
            if a1 is a2: continue
            if norm(h1.mean(0) - h2.mean(0)) > 4 * st.size: continue
            if hull_gap(h1, h2) < 1.5: issues.append(f"labels of {a1} and {a2} are cramped (< 1.5 pt apart)")
        lbx = {}                                           # labels of two unbonded atoms side by side on one line, less
        for m in self.mols:                                # than a word space apart, read as one word (CH3CH3, OHS)
            for a in m.atoms:
                if a._g and not a.hidden:
                    B = np.array([g["box"] for g in a._g]); lbx[a] = (B[:, 0].min(), B[:, 1].min(), B[:, 2].max(), B[:, 3].max(), m)
        for (a1, b1), (a2, b2) in itertools.combinations(lbx.items(), 2):
            if b1[4] is not b2[4] or b1[4].find(a1, a2) is not None: continue
            if min(b1[3], b2[3]) - max(b1[1], b2[1]) < 0.6 * min(b1[3] - b1[1], b2[3] - b2[1]): continue
            gap = max(b2[0] - b1[2], b1[0] - b2[2])
            if 0 <= gap < 0.4 * st.size:
                issues.append(f"labels of {a1} and {a2} run together on one line ({gap:.1f} pt apart): they read as one; "
                              "turn one of the groups, or move the atoms apart")
        for m in self.mols:
            for a in m.atoms:
                for mk in (a._chg, a._rad):
                    if not mk: continue
                    for h, a2 in glyphs:
                        if a2 is not a and _hull_dist(mk[0], h) < mk[1] + 1.2:
                            issues.append(f"charge/radical on {a} is cramped against the label of {a2}")
        tb = []
        for t in texts_all:
            if t.get("atom"): continue
            pts = [(x, y - 0.25 * sz, x + self.F.adv(s, sz, it, bd), y + self.F.cap * sz) for s, x, y, sz, it, bd in t["runs"] if s.strip()]
            if pts: tb.append((min(p[0] for p in pts), min(p[1] for p in pts), max(p[2] for p in pts), max(p[3] for p in pts), "".join(r[0] for r in t["runs"]), t))
        texts_only = tb
        for m in self.mols:                                # a charge typed into a label is neither placed nor counted
            for a in m.atoms:
                if not a.hidden and a.drawn and re.search(r"[⊕⊖]|\^\{?\d*[+−-]", a.label or "") and \
                        (m.nbrs(a) or a.lps or a.uid in react):    # (a spectator ion's formula, Ru(bpy)3^{2+}, may carry it)
                    issues.append(f"charge written into the label '{a.label}' of {a}: set {a}.charge instead (it is placed "
                                  "clear of bonds and arrows and counted by validate())")
        for m in self.mols:                                # a one-bond group (Me, OMe, Ph, R) on a vertex with more bonds
            for a in m.atoms:                              # is a label on the wrong atom
                k = sum({"dash": 0.5, "dot": 0, "dative": 0}.get(b.style, 1) for b in m.bonds   # skeleton bonds: a TS's
                        if a in (b.a, b.b) and not (b.b if b.a is a else b.a).hidden)              # partial one half, an
                k = int(k) if k == int(k) else k                                                   # H-bond none
                on = (a.charge or 0) > 0 and _ucore(label_units(a.label)[0]) in ("N", "O", "S", "P", "Se")  # EtO+-H, R2NH2+
                if not a.hidden and a.label in ONE_BOND and k > 1 + on:
                    issues.append(f"{a} is labelled {a.label} but has {k} bonds: a one-bond group on the wrong atom?")
        for t in texts_all:                               # every glyph must exist in the house font (no tofu boxes)
            for r in t["runs"]:
                miss = self.F.missing(r[0])
                if miss:
                    whole = "".join(x[0] for x in t["runs"])
                    issues.append(f"text '{whole}' has characters the house font cannot draw: " + ", ".join(f"'{c}' (U+{ord(c):04X})" for c in dict.fromkeys(miss))); break
        def ov(b1, b2, pad=0.6): return b1[0] < b2[2] + pad and b2[0] < b1[2] + pad and b1[1] < b2[3] + pad and b2[1] < b1[3] + pad
        for t1, t2 in itertools.combinations(texts_only, 2):
            if ov(t1, t2): issues.append(f"texts overlap: '{t1[4]}' / '{t2[4]}'")
        def hulls_of(runs):                                # glyph outlines of laid-out runs (pt)
            out = []
            for s_, x, y, sz, it, bd in runs:
                for ch in s_:
                    H = self.F.hull(ch, it, bd)
                    if H is not None and len(H) >= 3: out.append(H * sz + (x, y))
                    x += self.F.adv(ch, sz, it, bd)
            return out
        def gap(A, B):                                      # white between two sets of convex outlines
            if not A or not B: return 99.0
            EA = np.vstack([np.hstack([h, np.roll(h, -1, axis=0)]) for h in A]); sA = np.cumsum([0] + [len(h) for h in A[:-1]]).astype(int)
            EB = np.vstack([np.hstack([h, np.roll(h, -1, axis=0)]) for h in B]); sB = np.cumsum([0] + [len(h) for h in B[:-1]]).astype(int)
            return float(min(_polys_dist(np.vstack(A), EB, sB).min(), _polys_dist(np.vstack(B), EA, sA).min()))
        for it in self.items:                              # a reaction arrow (placed by row() or by hand) runs through no
            if not isinstance(it, RArrow): continue        # atom label, and its outline (an equilibrium arrow's halves)
            g = self._item_geo(it); own = [t_["runs"] for t_ in g["texts"]]; QL, QP = [], []   # through no other text
            for pr in g["prims"]:
                if pr[0] not in ("line", "poly"): continue
                P = np.asarray(pr[1], float) if pr[0] == "line" else np.vstack([pr[1], pr[1][:1]]).astype(float)
                seg = np.diff(P, axis=0); Lg = np.hypot(*seg.T)
                (QL if pr[0] == "line" or len(pr[1]) != 5 else QP).extend(P[k] + seg[k] * tt for k in range(len(seg))   # (free text:
                                                                         for tt in np.linspace(0, 1, max(2, int(Lg[k] / 0.4) + 1)))  # half arrows only)
            for t in texts_all:
                if not t["runs"] or t["runs"] in own: continue
                Q = np.array(QL + QP if t.get("atom") else QP)
                if not len(Q): continue
                hs = hulls_of(t["runs"])
                if not hs: continue
                allh = np.vstack(hs)
                if (Q.max(0) < allh.min(0) - 2).any() or (Q.min(0) > allh.max(0) + 2).any(): continue
                E = np.vstack([np.hstack([h, np.roll(h, -1, axis=0)]) for h in hs]); s0 = np.cumsum([0] + [len(h) for h in hs[:-1]]).astype(int)
                if float(_polys_dist(Q, E, s0).min()) < st.lw / 2 + 0.5:
                    issues.append(f"a reaction arrow runs through {'the label' if t.get('atom') else 'the text'} '{''.join(r[0] for r in t['runs'])}'")
        for it in self.items:                              # a reaction arrow's text keeps clear of its shaft and head: the
            if not isinstance(it, RArrow): continue        # nearest ink at least 2.5 pt from the line's edge (capitals sit
            g = self._item_geo(it)                         # 4.6 pt off; a superscript once rode up to 1.3 pt)
            for kind, lim in (("line", 2.5), ("poly", 1.2)):   # shaft lines; heads (their barbs sit 2.75 pt off the shaft,
                segs = [np.array(pr[1], float) for pr in g["prims"] if pr[0] == kind and len(pr[1]) > 1]   # so text under
                if not segs: continue                                                                         # a head is
                E = np.vstack([np.hstack([P[:-1], P[1:]] if kind == "line" else [P, np.roll(P, -1, axis=0)]) for P in segs])  # closer)
                for t in g["texts"]:
                    H = hulls_of(t["runs"])
                    if not H: continue
                    dmin = float(_segs_dist(np.vstack(H), E).min()) - (st.lw / 2 if kind == "line" else 0.0)
                    if dmin < lim * st.bond / 14.4:
                        issues.append(f"arrow text '{''.join(r[0] for r in t['runs'])}' comes within {dmin:.1f} pt of its arrow")
        for t in texts_only:
            for bx, a in boxes:
                if not ov(t, bx): continue                  # boxes are only a first sieve: a minus or a dot is far smaller
                tg = hulls_of(t[5]["runs"]); lg = [h for g in (a._g or []) for h in g["hulls"]]
                if gap(tg, lg) < 0.6 * st.bond / 14.4: issues.append(f"text '{t[4]}' overlaps the label of {a}")
            for m, b in allb:
                for p0, p1, ex in self._bond_lines(b):
                    if any(box_dist(p, t[:4]) < 0.6 + ex for p in np.linspace(p0, p1, 14)):
                        issues.append(f"text '{t[4]}' runs into bond {b.a}-{b.b}"); break
        for m in self.mols:
            for a in m.atoms:
                for nm, t in [("charge", a._chg[2] if a._chg else None), ("radical", a._rad[2] if a._rad else None)] + \
                             [("lone pair", l.angle) for l in a.lps]:
                    if t is not None and self._in_ring(a, t): issues.append(f"{nm} on {a} is drawn inside a ring")
        for m in self.mols:                                # a bare vertex drawn on a bond it is not part of (a ring
            if getattr(m, "_newman", None): continue       # swung onto another): it reads as a junction of that bond
            V = [a for a in m.atoms if not a.drawn and not a.hidden and getattr(a, "_tmpl", None) != "chair"]
            if not V: continue
            P = np.array([self.P(a) for a in V]); hit = set()
            for b in m.bonds:                              # (its lines as drawn: a gap where it passes behind clears it)
                if b.a.hidden or b.b.hidden: continue
                for p0, p1, ex in self._bond_lines(b):
                    d_ = p1 - p0; L2 = float(d_ @ d_)
                    if L2 < 1e-9: continue
                    tt = np.clip(((P - p0) @ d_) / L2, 0, 1); D = np.hypot(*(P - (p0 + tt[:, None] * d_)).T) / st.bond
                    for k in np.nonzero(D < 0.15)[0]:
                        if V[k] in (b.a, b.b) or V[k].uid in hit: continue
                        hit.add(V[k].uid); issues.append(f"{V[k]} lies on bond {b.a}-{b.b} ({D[k]:.2f} BU from it): move one of them")
        for m, b in allb:
            for bx, a in boxes:
                if a in (b.a, b.b): continue
                for p0, p1, ex in self._bond_lines(b):
                    pts = np.linspace(p0, p1, 24); d = min(box_dist(p, bx) for p in pts)
                    if d >= st.margin + 1.0 + ex: continue      # near: the white between the line and the letters
                    hs = [h for g in (a._g or []) for h in g["hulls"]]   # themselves, at least a bond end's own margin
                    if d < 1.0 + ex or (hs and min(_hull_dist(p, h) for p in pts for h in hs) - st.lw / 2 - ex < st.margin):
                        issues.append(f"bond {b.a}-{b.b} runs into label {a}"); break   # (less reads as touching: OTBS|)
        for m, b in allb:                                  # a dashed or dotted bond ends where a solid one would at its
            if b.style not in ("dash", "dot"): continue    # own atoms' labels (no dot touching the H of N-H...Cl)
            for a, o in ((b.a, b.b), (b.b, b.a)):
                if not (a.drawn and not a.hidden and a._g): continue
                u = unit(self.P(o) - self.P(a))
                for pr in b._prims:
                    if pr[0] != "line" or len(pr[1]) != 2: continue
                    e = min(pr[1][0], pr[1][-1], key=lambda x: norm(np.asarray(x) - self.P(a)))
                    cap = pr[2] / 2 if pr[5:6] == ("round",) else 0.0
                    if float((np.asarray(e, float) - self.P(a)) @ u) - cap < self._exit(a, u) - 0.3:
                        issues.append(f"the {'dotted' if b.style == 'dot' else 'dashed'} bond {b.a}-{b.b} runs into label {a}: "
                                      f"it should end where a solid bond would"); break
        for m in self.mols:
            for a in m.atoms:
                marks = [(c, st.lp_r) for _, c, d1, d2 in a._lp] + [(mk[0], mk[1]) for mk in (a._chg, a._rad) if mk]
                for c, r in marks:
                    for _, b in allb:
                        for p0, p1, ex in self._bond_lines(b):
                            if min(norm(p - c) for p in np.linspace(p0, p1, 16)) < r + 0.6 + ex:
                                issues.append(f"mark on {a} touches bond {b.a}-{b.b}")
                    for bx, a2 in boxes:
                        if box_dist(c, bx) < r + 0.3 and a2 is not a: issues.append(f"mark on {a} touches label {a2}")
        for pt in self.pointers:                           # a pointer's arrows run clear of bonds and labels to their target
            g = self._ptr_geo(pt); Q = []
            for pr in g["prims"]:
                P = np.asarray(pr[1], float)
                if pr[0] == "line":
                    seg = np.diff(P, axis=0); Lg = np.hypot(*seg.T)
                    Q += [P[k] + seg[k] * t for k in range(len(seg)) for t in np.linspace(0, 1, max(2, int(Lg[k] / 0.4) + 1))]
                elif pr[0] == "poly": Q += list(P)
            if not Q: continue
            Q = np.array(Q); tname = ", ".join(str(t) for t in pt.targets)
            for m, b in allb:
                for p0, p1, ex in self._bond_lines(b):
                    if float(_segs_dist(Q, np.array([[*p0, *p1]])).min()) < st.lw + ex + 0.6:
                        issues.append(f"the pointer to {tname} runs into bond {b.a}-{b.b} (give it another side)"); break
            for m in self.mols:
                for a in m.atoms:
                    hs = [h for x in (a._g or []) for h in x["hulls"]]
                    if not hs or a in pt.targets: continue
                    E = np.vstack([np.hstack([h, np.roll(h, -1, axis=0)]) for h in hs]); s0 = np.cumsum([0] + [len(h) for h in hs[:-1]]).astype(int)
                    if float(_polys_dist(Q, E, s0).min()) < st.lw / 2 + 0.8:
                        issues.append(f"the pointer to {tname} runs through the label of {a} (give it another side)")
            heads = [np.asarray(pr[1], float)[0] for pr in g["prims"] if pr[0] == "poly"]
            for t, tip in zip(pt.targets, heads):          # each arrow ends at its target
                if isinstance(t, Atom) and t.drawn and t._g:
                    hs = [h for x in t._g for h in x["hulls"]]; E = np.vstack([np.hstack([h, np.roll(h, -1, axis=0)]) for h in hs])
                    miss = float(_polys_dist(tip[None], E, np.cumsum([0] + [len(h) for h in hs[:-1]]).astype(int)).min()) - 4.0
                elif isinstance(t, Atom): miss = norm(tip - self.P(t)) - 0.95 * st.bond
                elif isinstance(t, Mid): miss = norm(tip - self._pt(t)) - 4.0
                else: miss = float(_segs_dist(tip[None], np.array([[*self.P(t[0]), *self.P(t[1])]])).min()) - 6.0
                if miss > 0: issues.append(f"the pointer to {t} stops {miss:.1f} pt short of it")
            own = {id(m) for m in pt.mols()}
            for t in g["texts"]:                           # its note is not taken for another item's: not over a reaction
                nb = t["box"]; note = "".join(r[0] for r in t["runs"])   # arrow (where it reads as a condition), nor on
                for it in self.items:                                     # another structure
                    if not isinstance(it, RArrow) or it.vertical: continue
                    gi = self._item_geo(it); yy = it.y
                    P = np.array([q for pr in gi["prims"] if pr[0] in ("line", "poly") for q in np.asarray(pr[1], float).reshape(-1, 2)])
                    x0, x1 = P[:, 0].min(), P[:, 0].max(); up = nb[1] >= yy
                    between = any(min(nb[2], tb[2]) > max(nb[0], tb[0]) and (yy <= tb[1] and tb[3] <= nb[1] + 1.0 if up else nb[3] - 1.0 <= tb[1] and tb[3] <= yy)
                                  for tb in (t_["box"] for t_ in gi["texts"]))       # (the arrow's own text between them: no)
                    if not between and min(nb[2], x1) - max(nb[0], x0) > 0.3 * (nb[2] - nb[0]) and min(abs(nb[1] - yy), abs(nb[3] - yy)) < 1.2 * st.bond:
                        issues.append(f"the note '{note}' stands over a reaction arrow, where it reads as a condition: declare the "
                                      "pointer before the row() that places its structure, or give it another side")
                for m in self.mols:
                    if id(m) not in own and box_gap(nb, self._mol_box(m, curlies=False)) < 0.15 * st.bond:
                        issues.append(f"the note '{note}' runs into another structure ({m.atoms[0]}...): give it another side")
        clear = self.TEXT_CLEAR * st.bond / 14.4           # free text keeps a third of a bond from structures it is not
        def mol_ink(m):                                    # part of (a note under the F of a neighbour's aryl reads as
            L = [[*p0, *p1] for b in m.bonds for p0, p1, ex in self._bond_lines(b)]                                  # its)
            for e in getattr(self, "extras", []):          # (orbital lobes drawn on it are its ink too)
                if e[0] == "orbital" and e[1].mol is m:
                    for q in self._orbital_prims(e):
                        if q[0] == "orb": P = np.asarray(q[1], float); L += [[*P[i], *P[i + 1]] for i in range(len(P) - 1)]
            L = np.array(L, float).reshape(-1, 4)
            return L, [h for a in m.atoms if a.drawn and not a.hidden for g in (a._g or []) for h in g["hulls"]]
        inks = {id(m): mol_ink(m) for m in self.mols}; mbox = {id(m): self._mol_box(m, curlies=False, prep=False) for m in self.mols}
        def near(runs, own, what, txt):
            hs = hulls_of(runs)
            if not hs: return
            HP = np.vstack(hs); bx = (*HP.min(0), *HP.max(0))
            for m in self.mols:
                if id(m) in own or box_gap(bx, mbox[id(m)]) > clear: continue
                L, LH = inks[id(m)]
                d = min([float(_segs_dist(HP, L).min()) - st.lw / 2] if len(L) else [99.0])
                if LH: d = min(d, gap(hs, LH))
                if d < clear - 0.25:
                    issues.append(f"the {what} '{' '.join(str(txt).split())}' comes within {max(d, 0):.1f} pt of {m.atoms[0]}...'s "
                                  f"structure (keep {clear:.0f} pt): move it, or give it another side"); return
        for t in texts_all:                                # (as drawn) a pointer's note: clear of all ink, its own
            if t.get("note") is not None: near(t["runs"], set(), "note", t["note"].text)          # structure's too; a
            elif t.get("cap") is not None:                                                         # caption: of the
                cp = t["cap"]; near(t["runs"], {id(x) for x in (cp[0] if isinstance(cp[0], (tuple, list)) else (cp[0],))},   # ones
                                    "caption", cp[1])                                              # it is not under; free
            elif t.get("free") is not None and t["free"].strip() not in ("+", "⊕", "⊖"):            # text (not a + between
                near(t["runs"], set(), "text", t["free"])                                          # structures): of all
        issues += self._arrow_issues(); notes += self._arrow_notes
        issues += self._head_issues()
        issues += self._fragment_issues()
        metals = set("Li Na K Mg Ca B Al Ti V Cr Mn Fe Co Ni Cu Zn Ru Rh Pd Ag Ir Pt Au Hg Tc Re Os Mo W Sn M L".split())
        for m in self.mols:
            for a in m.atoms:
                el = re.sub(r"[^A-Za-z].*", "", label_units(a.label)[0]) if a.drawn else "C"
                nb = m.nbrs(a)
                if getattr(a, "_tmpl", None) == "cage" and len(nb) >= 2:       # a cage in perspective: still no near-overlap
                    ds = sorted(ang(n.xy - a.xy) % 360 for n in nb)
                    gp = min(((ds[(i + 1) % len(ds)] - ds[i]) % 360) for i in range(len(ds)))
                    if gp < 25: issues.append(f"bonds at {a} only {gp:.0f} deg apart")
                if getattr(a, "_tmpl", None) == "chair" and len(nb) >= 3:      # perspective template: still no near-overlap
                    ds = sorted(ang(n.xy - a.xy) % 360 for n in nb)
                    gp = min(((ds[(i + 1) % len(ds)] - ds[i]) % 360) for i in range(len(ds)))
                    if gp < 40: issues.append(f"bonds at {a} only {gp:.0f} deg apart (on a chair, use the axial and equatorial directions from chair())")
                if el in metals or len(nb) < 2 or getattr(a, "_tmpl", None): continue
                nbs = sorted(nb, key=lambda n: ang(n.xy - a.xy) % 360); ds = [ang(n.xy - a.xy) % 360 for n in nbs]
                gaps = [((ds[(i + 1) % len(ds)] - ds[i]) % 360) for i in range(len(ds))]
                tsb = [m.find(a, n).style in ("dash", "dot") for n in nbs]   # two bonds forming or breaking in a cyclic TS
                if len(nbs) >= 4 and all(m.ring_of(m.find(a, n)) for n in nbs): tsb = [True] * len(nbs)   # (four ring
                close = [g for i, g in enumerate(gaps) if g < (30 if tsb[i] and tsb[(i + 1) % len(ds)] else 55)]   # bonds
                # at one atom: the rings' polygons set the gaps, 36 deg in an angular triquinane)
                if close: issues.append(f"bonds at {a} only {min(close):.0f} deg apart")   # a 3- or 4-ring: 30 deg will do
                inring = any(m.ring_of(b) for b in m.bonds if a in (b.a, b.b))
                if m._lopsided(a):
                    gs = sorted(gaps); shape = "make a T" if gs[0] >= 80 and gs[2] <= 190 else "are spread unevenly"
                    issues.append(f"the three bonds at {a} {shape} ({', '.join(f'{g:.0f}' for g in gs)} deg): draw a three-bond "
                                  f"centre with its bonds 120 deg apart (release() redraws a centre it leaves like this)")
                if len(nb) == 2 and not inring and not any(e[0] == "orbital" and e[1] is a for e in self.extras):
                    # (an orbital picture's frame is in perspective: R2C: as 60 deg)
                    th = min(gaps); bs_ = [m.find(a, n) for n in nbs]; bo_ = sorted(b_.order for b_ in bs_)
                    if all(b_.style in ("plain", "bold") for b_ in bs_) and (bo_ == [2, 2] or bo_[-1] == 3) and th < 172:
                        issues.append(f"{a} has {'a triple bond' if bo_[-1] == 3 else 'two double bonds'}: an sp centre is "
                                      f"straight (180 deg), drawn at {th:.0f} deg; turn one side through to 180 (swing())")
                    elif not (110 <= th <= 130 or th >= 172):
                        (notes if (100 <= th <= 140 and near_react(a)) else issues).append(f"bond angle at {a} is {th:.0f} deg (use 120 or 180)")
        for m in self.mols:
            for a in m.atoms:
                el = re.sub(r"[^A-Za-z].*", "", label_units(a.label)[0]) if a.drawn else "C"
                if el in metals or getattr(a, "_tmpl", None): continue
                rb = [b for b in m.bonds if a in (b.a, b.b) and m.ring_of(b)]
                if len(rb) != 2: continue
                ex = [n for n in m.nbrs(a) if not any(n in (b.a, b.b) for b in rb)]
                r1, r2 = [ang((b.b if b.a is a else b.a).xy - a.xy) for b in rb]
                bis = ang(-(dirv(r1) + dirv(r2))) if norm(dirv(r1) + dirv(r2)) > 1e-6 else None
                if bis is None: continue
                cen = np.mean([x.xy for x in m.ring_of(rb[0])], axis=0)   # (away from the ring's centre: at a reflex
                if float(dirv(bis) @ (cen - a.xy)) > 0: bis = bis + 180     # vertex the bonds' sum points outwards)
                dev = abs(adiff(ang(ex[0].xy - a.xy), bis)) if len(ex) == 1 else 0
                gem_ok = (not a.drawn) and sum(b.order for b in m.bonds if a in (b.a, b.b)) + abs(a.charge) + a.radical < 4
                delib = getattr(a, "deliberate", False) or any(getattr(e, "deliberate", False) for e in ex)
                bucket = notes if ((near_react(a) and dev <= 35) or delib) else issues    # >35 deg only if marked deliberate
                if len(ex) == 1 and dev > 12 and not (gem_ok and abs(dev - 30) <= 12):
                    bucket.append(f"substituent {ex[0]} on ring atom {a} is {abs(adiff(ang(ex[0].xy - a.xy), bis)):.0f} deg off the exterior bisector")
                if len(ex) == 2:
                    d1, d2 = sorted(adiff(ang(e.xy - a.xy), bis) for e in ex)
                    if abs(d1 + d2) > 12: bucket.append(f"gem substituents on {a} not symmetric about the exterior bisector")
            for a in m.atoms:
                sb = [b for b in m.bonds if b.a is a and b.style in ("wedge", "hash")]
                if not sb or (a.drawn and re.sub(r"[^A-Za-z].*", "", label_units(a.label)[0]) in metals): continue   # (a
                plain = [b for b in m.bonds if a in (b.a, b.b) and b not in sb and b.style not in ("dash", "dot")]   # metal's
                # wedges show the perspective of its square plane or octahedron, not a stereocentre)
                if len(sb) + len(plain) == 4 and len(plain) == 3:
                    ds = sorted(ang((b.b if b.a is a else b.a).xy - a.xy) % 360 for b in plain)
                    gaps = [((ds[(i + 1) % 3] - ds[i]) % 360) for i in range(3)]
                    if max(gaps) >= 180: issues.append(f"stereo at {a} is ambiguous: its three plain bonds lie within 180 deg")
                if len(sb) >= 2 and len(plain) + len(sb) == 4:
                    ds = [ang(b.b.xy - a.xy) for b in sb]
                    if any(abs(abs(adiff(x, y)) - 180) < 15 for x, y in itertools.combinations(ds, 2)) and len({b.style for b in sb}) == 1:
                        issues.append(f"stereo at {a}: two {sb[0].style}s drawn opposite each other")
                for b in sb:
                    if getattr(b.a, "_tmpl", None) == "cage" and getattr(b.b, "_tmpl", None) == "cage": continue
                    if len(m.nbrs(b.b)) >= 3 and any(x.a is b.b and x.style in ("wedge", "hash") for x in m.bonds):
                        notes.append(f"{b.style} {b.a}-{b.b} joins two stereocentres (IUPAC prefers wedges on exocyclic bonds)")
            orb = {id(e[1]) for e in getattr(self, "extras", []) if e[0] == "orbital"}
            for b in m.bonds:
                if getattr(b.a, "_tmpl", None) == "persp": continue           # a perspective ring's near edges: depth, not stereo
                if id(b.a) in orb: continue                                    # an orbital picture's sigma frame in perspective
                if getattr(b.a, "_tmpl", None) == "cage" and getattr(b.b, "_tmpl", None) == "cage": continue   # a pinned
                                                                               # cage's front edge, drawn as a wedge: depth
                if b.style in ("wedge", "hash") and len(m.nbrs(b.a)) < 3 and not (b.a.drawn and re.sub(r"[^A-Za-z].*", "", label_units(b.a.label)[0]) in metals):
                    issues.append(f"{b.style} {b.a}-{b.b} starts at {b.a}, which is not a stereocentre (narrow end goes on the stereocentre)")
        for m in self.mols:
            for h in m.atoms:
                if h.label != "H" or len(m.nbrs(h)) != 1: continue
                x = m.nbrs(h)[0]; others = [n for n in m.nbrs(x) if n is not h]
                near = [y for y in m.atoms if y is not h and y is not x and y not in m.nbrs(h) and not y.hidden   # (its own structure:
                        and y not in others and norm(y.xy - h.xy) < 0.85              # a neighbour in the row is spacing)
                        and not (isinstance(y, Atom) and y.label == "H")]
                # (its parent's own neighbours are the bond-angle rule's business: a short fusion H in its free sector
                # stands 0.85 BU from them by geometry)
                if near: issues.append(f"H on {x} points into the structure (near {near[0]})")
        if self.curlies:                                   # P1-5: explicit X-H that nothing uses (standard depiction is OH, NH2)
            here = {id(y) for c in self.curlies for x in (c.src, c.dst) for y in self._xatoms(x)}
            anywhere = {y.uid for c in self.curlies for x in (c.src, c.dst) for y in self._xatoms(x)}
            for m in self.mols:
                stepm = any(id(a) in here for a in m.atoms)
                for X in m.atoms:
                    if not X.drawn or X.hidden or re.search(r"[_^<]", X.label): continue
                    if re.sub(r"\d+'*$", "", label_units(X.label)[0]) not in ("O", "N", "S", "P", "Se"): continue
                    hs = [h for h in m.nbrs(X) if h.label == "H" and len(m.nbrs(h)) == 1]
                    if not hs or any(id(h) in here for h in hs) or getattr(X, "deliberate", False) or any(getattr(h, "deliberate", False) for h in hs): continue
                    if stepm and any(h.uid in anywhere for h in hs): continue        # the H is used in another step
                    notes.append(f"explicit H on {X!r} with no arrow on it: fold it into the label (Mol.tidy()) unless it is shown on purpose")
        val = {"C": 4, "N": 3, "O": 2, "S": 2, "B": 3, "Al": 3, "F": 1, "Cl": 1, "Br": 1, "I": 1, "H": 1, "P": 3}
        VE = {"C": 4, "N": 5, "O": 6, "S": 6, "B": 3, "Al": 3, "F": 7, "Cl": 7, "Br": 7, "I": 7, "P": 5}   # valence electrons
        orb = {id(e[1]) for e in getattr(self, "extras", []) if e[0] == "orbital"}
        def subs(units, bo=1, onium=False):
            """Substituents written into a label after its anchor (H3, Me2, Et, Boc, tBu, F3, R2 ...), or None if the
            label holds anything else (O2 of CO2Me, a bracket), when its valence cannot be read from the label. A numbered
            group (R¹ of NR¹R²) is one substituent, a counted one (R₂ of NR₂) is several."""
            n = 0
            for k, x in enumerate(units):
                if k == 0: continue
                mm = re.fullmatch(r"(H|(?:tert|sec|iso|t|i|n|s|c)?(?:%s)|R|R'|X|F|Cl|Br|I)(\d*)" % "|".join(ABBR), x)
                if not mm: return None
                n += 1 if r_ident(units, k, bo, onium) else int(mm.group(2) or 1)
            return n
        for m in self.mols:
            for a in m.atoms:
                if a.hidden or not a.drawn or re.search(r"[_^<]", a.label): continue
                if not (m.nbrs(a) or a.lps or a.charge or a.radical): continue
                u = label_units(a.label)
                if len(u) < 2: continue
                u0 = re.sub(r"\d+'*$", "", u[0]); u0b = re.sub(r"^(?:tert|sec|iso|t|i|n|s|c)(?=[A-Z])", "", u0)
                het = [x for x in (re.sub(r"\d+'*$", "", y) for y in u[1:]) if x in HETMET]
                if a.label in ("OC", "ON", "SC") and any(re.sub(r"[^A-Za-z].*", "", n.label) in METALS for n in m.nbrs(a)):
                    issues.append(f"label '{a.label}' on {a} is anchored on its {a.label[0]}, but a carbonyl, nitrosyl or thiocarbonyl "
                                  f"binds the metal through {a.label[1]}: write '{a.label[::-1]}' (the engine flips it to read {a.label} on the left)")
                if (u0 == "H" or u0b in ABBR) and het and (m.nbrs(a) or a.lps or a.uid in react):   # (a spectator ion
                                                                    # with only its charge, HSO4-, reads as typed: no anchor in play)
                    issues.append(f"label '{a.label}' on {a} is anchored on '{u0}', so its bonds, charges and lone pairs attach "
                                  f"there and not to {het[0]}: write the {het[0]} first (the engine flips it to read correctly)")
        pz = {id(e[1]): 1 for e in getattr(self, "extras", []) if e[0] == "orbital" and e[2] == "p"}   # a p orbital drawn
        for m in self.mols:                                                     # on an atom is its pi bond (C=O as lobes)
            for a in m.atoms:
                if a.hidden or not a.label or getattr(a, "generic", False): continue
                u = label_units(a.label) if not re.search(r"[_^<]", a.label) else ["?"]
                el = re.sub(r"\d", "", u[0]); hs = 0
                if el not in val: continue
                bo = sum({"dative": 0, "dash": 0, "dot": 0}.get(b.style, 1.5 if b.style == "partial" else b.order)
                         for b in m.bonds if a in (b.a, b.b))   # a dashed bond is
                                                                       # a partial one (a transition state's), not a full bond
                if len(u) >= 2:
                    hs = subs(u, sum(b.order for b in m.bonds if a in (b.a, b.b)),
                              (a.charge or 0) > 0 or any(re.sub(r"[^A-Za-z].*", "", n.label) in METALS for n in m.nbrs(a)))
                    if hs is None: continue
                if any(re.sub(r"[^A-Za-z].*", "", n.label) in (METALS - {"B"}) | {"L"} for n in m.nbrs(a)): continue
                if re.fullmatch(r"[A-Z][a-z]?\d+", u[0]): continue   # N2Cl, N3: the anchor is several atoms (an ion pair)
                v = val[el] + (a.charge if el in ("N", "O", "S", "P", "F", "Cl", "Br", "I") else -a.charge if el in ("B", "Al")
                               else -abs(a.charge)) - (1 if a.radical else 0)   # R2Cl+ two bonds, BH4- and AlCl4- four
                if el == "C": v -= 2 * max(0, len(a.lps) + bool(getattr(a, "ambi", False)) - (1 if a.charge < 0 else 0))   # a
                                                                       # carbene's pair (:CH2, :CCl2, or its + -)
                plain = not any(b.style in ("dative", "dash", "dot", "partial") for b in m.bonds if a in (b.a, b.b)) and id(a) not in orb
                if a.drawn and a.lps and el in VE and plain:   # the drawn pairs fit the electrons left by bonds and charge
                    room = VE[el] - a.charge - (1 if a.radical else 0) - (bo + hs)
                    if 2 * len(a.lps) > room:
                        issues.append(f"{a} has {len(a.lps)} lone pair{'s' * (len(a.lps) > 1)} drawn, but its bonds and charge "
                                      f"leave room for {max(0, room) // 2}: one is on the wrong atom")
                if el == "C" and not a.drawn:
                    h0 = getattr(a, "_h0", None)             # its H count from the SMILES, while bonds and charge are as there
                    if bo > v:
                        issues.append(f"carbon {a} has {bo:g} bonds (max {v} with charge/radical)" if not a.lps else
                                      f"lone pair on {a}, which has {bo:g} bonds: no room for it (a carbanion needs its charge)")
                    elif h0 is not None and h0[1:] == (bo, a.charge) and (a.lps or a.radical) and bo + h0[0] > v and plain:
                        k = "lone pair" if a.lps else "radical"
                        issues.append(f"{k} on {a}: its bonds and its {h0[0]} H (from the SMILES) leave no room for it; it is on "
                                      f"the wrong atom, or write the {'carbanion or carbene' if a.lps else 'radical'} in the SMILES "
                                      f"({'[CH-], [C]' if a.lps else '[CH], [C]'})")
                elif bo + hs != v and not (el in "SP" and bo + hs > v) and bo + hs + pz.get(id(a), 0) != v and \
                        not (el in ("Cl", "Br", "I") and a.charge == 1 and bo + hs == 0) and \
                        all(abs(bo + hs + t * sum(1 for b in m.bonds if a in (b.a, b.b) and b.style in ("dash", "dot")) - v) > 1e-6
                            for t in (0.5, 1)):                # partial bonds: a transition state's halves (B in a four-centred
                                                               # TS, the H of O...H...A), or coordination either way (O...Li)
                    issues.append(f"valence: {a} uses {bo + hs:g}, expected {v} (charge {a.charge:+d})")
        for m in self.mols:                                # a charge or radical belongs right by its atom: one pushed
            for a in m.atoms:                              # out by crowding reads as floating free
                for nm, mk in (("charge", a._chg), ("radical", a._rad)):
                    if not mk or a.hidden: continue
                    c, r = np.asarray(mk[0], float), mk[1]
                    g, lim = self._mark_gap(a, c, r, mk[2], "chg" if nm == "charge" else "rad")   # as placement measures it
                    dots =[d for _, _, d1, d2 in getattr(a, "_lp", []) for d in (d1, d2)]   # set just outside its own lone
                    if dots and min(norm(d - c) for d in dots) - r - st.lp_r <= 1.6: g = 0.0    # pairs it is theirs (⊖:Ö–H),
                    if self._by_own_lobe(a, c, r): g = 0.0     # and beside its own orbital's lobe (C⊕ by its empty p orbital)
                    if g > lim:
                        issues.append(f"the {nm} on {a} sits {g:.1f} pt from it (at most {lim:.1f} pt): crowding pushed it out; give it "
                                      f"a clearer side ({nm}_at=deg) or more room round the atom")
        if not quiet:
            print("check: OK" if not issues else "check:\n  " + "\n  ".join(dict.fromkeys(issues)))
            if notes: print("notes (deliberate near reacting centres? confirm):\n  " + "\n  ".join(dict.fromkeys(notes)))
        self._notes = list(dict.fromkeys(notes))
        return list(dict.fromkeys(issues))
    def balance(self, before, after, quiet=False):
        """Mass balance between two steps (each a Mol or list of Mols): the atoms gained and lost as Hill formulas, and
        the net charge before and after. A loss or gain that no arrow text in the scene mentions is flagged, so a
        by-product dropped from the drawing (-H2O, -HBr) never vanishes silently. Returns {'lost', 'gained', 'charge'}."""
        fb, qb = formula(before); fa, qa = formula(after)
        lost, gained = fb - fa, fa - fb
        written = [f for it in self.items if isinstance(it, RArrow) for side in (it.above, it.below)
                   for t in (side if isinstance(side, (list, tuple)) else [side]) if t is not None
                   for f in ([formula(t)[0]] if isinstance(t, Mol) else formulas_in(t) if t else [])]
        said = lambda c: not c or any(+f == +c for f in written)
        lines = []
        if lost: lines.append(f"lost {hill(lost)}" + ("" if said(lost) else " (no arrow mentions it)"))
        if gained: lines.append(f"gained {hill(gained)}" + ("" if said(gained) else " (no arrow mentions it)"))
        if qb != qa: lines.append(f"charge {qb:+d} before, {qa:+d} after")
        if not quiet: print("balance: " + ("; ".join(lines) if lines else "OK"))
        return {"lost": hill(lost), "gained": hill(gained), "charge": (qb, qa)}
    def validate(self, before, after, quiet=False):
        """Replay the curly arrows drawn on `before` (Mol or list of Mols) and compare the predicted bonds,
        formal charges and radicals with `after`, atom by atom (atoms matched by uid via Mol.copy)."""
        B = before if isinstance(before, (list, tuple)) else [before]
        A = after if isinstance(after, (list, tuple)) else [after]
        batoms = {a.uid: a for m in B for a in m.atoms}; aatoms = {a.uid: a for m in A for a in m.atoms}
        def inb(x):
            if isinstance(x, LP): return x.atom.uid in batoms and batoms[x.atom.uid] is x.atom
            if isinstance(x, Atom): return x.uid in batoms and batoms[x.uid] is x
            if isinstance(x, Mid): return inb(x.a) and inb(x.b)
            if isinstance(x, (tuple, list)): return all(inb(y) for y in x)
            return False
        arrows = [c for c in self.curlies if inb(c.src) and inb(c.dst)]
        dN, dE, issues = {}, {}, []
        key = lambda x, y: frozenset((x.uid, y.uid))
        def addN(a, v): dN[a.uid] = dN.get(a.uid, 0) + v
        def addE(x, y, v): dE[key(x, y)] = dE.get(key(x, y), 0) + v
        hit = {}                                   # atoms each arrow attacks (its target atom, or the ends of a forming bond)
        for c in arrows: hit[id(c)] = {y.uid for y in self._hit(c)}
        for c in arrows:
            n = 2 if c.kind == "full" else 1; s, d = c.src, c.dst; srcatoms = []
            if isinstance(s, LP):
                if c.kind == "fish": issues.append(f"fishhook from a lone pair on {s.atom}: use a full arrow")
                addN(s.atom, -n); srcatoms = [s.atom]
            elif isinstance(s, Atom):
                if c.kind == "full" and s.radical and s.charge >= 0: issues.append(f"two-electron arrow from radical {s}: use fish()")
                addN(s, -n); srcatoms = [s]
            elif isinstance(s, (tuple, list)):
                if s[0].mol.find(*s) is None: issues.append(f"arrow tail on {s[0]}-{s[1]} but no bond is drawn there")
                addE(s[0], s[1], -n); srcatoms = list(s)
            else: issues.append("arrow tail must be a lone pair, atom or bond"); continue
            if isinstance(d, Atom):
                if isinstance(s, (tuple, list)) and d in s: addN(d, n)
                elif isinstance(s, (tuple, list)):
                    # a bond's electrons onto an atom bonded to one of its ends: a 1,2-shift, the group at the other end moving
                    # with them (C2-R onto C1+ makes R-C1); if another arrow takes that group (a base removing the H), they
                    # stay behind and make the double bond instead (E1). Otherwise the nearer atom bonds to d (a pi bond
                    # attacking an electrophile).
                    others = set().union(*[hit[id(x)] for x in arrows if x is not c])
                    e = self._shift_end(s, d); w = self._migrant(s, d, others) or (e[0] if e else None)
                    if w is None and d.mol is not s[0].mol:  # an atom in another item: distances mean nothing across
                        X = [z for z in s if z.uid in aatoms and d.uid in aatoms and aatoms[z.uid].mol is aatoms[d.uid].mol
                             and aatoms[z.uid].mol.find(aatoms[z.uid], aatoms[d.uid]) is not None]   # items, so the end the
                        w = X[0] if len(X) == 1 else None                                        # next step bonds to it
                    w = w or min(s, key=lambda z: norm(z.xy - d.xy)); addE(w, d, n)
                else: addE(srcatoms[0], d, n)
            elif isinstance(d, Mid): addE(d.a, d.b, n)
            elif isinstance(d, (tuple, list)): addE(d[0], d[1], n)
        bo = lambda m, x, y: (lambda b: 0 if b is None or b.style in ("dot",) else (0.5 if b.style == "dash" else b.order))(m.find(x, y))
        common = [u for u in batoms if u in aatoms]
        for u, v in itertools.combinations(common, 2):
            x, y = batoms[u], batoms[v]
            if x.mol is not y.mol and (x.mol.find(x, y) is None): before_o = 0
            else: before_o = bo(x.mol, x, y)
            e = dE.get(frozenset((u, v)), 0)
            if e % 2: issues.append(f"odd number of electrons moved into {x}-{y} (a fishhook is missing its partner)")
            pred = before_o + e / 2; X, Y = aatoms[u], aatoms[v]
            after_o = bo(X.mol, X, Y) if X.mol is Y.mol else 0
            if abs(pred - after_o) > 1e-6 and not (before_o % 1 or after_o % 1):
                issues.append(f"bond {x}-{y}: arrows give order {pred:g}, next step shows {after_o:g}")
        for u in common:
            x, X = batoms[u], aatoms[u]
            eb = sum(v for k, v in dE.items() if u in k)
            dfc = -dN.get(u, 0) - eb / 2
            if abs(x.charge + dfc - X.charge) > 1e-6:
                issues.append(f"charge on {x}: arrows give {x.charge + dfc:+g}, next step shows {X.charge:+d}")
            rad = bool(x.radical) ^ bool(dN.get(u, 0) % 2)
            if rad != bool(X.radical): issues.append(f"radical on {x}: arrows give {'a radical' if rad else 'none'}, next step shows {'a radical' if X.radical else 'none'}")
        if not quiet: print("validate: OK" if not issues else "validate:\n  " + "\n  ".join(dict.fromkeys(issues)))
        return list(dict.fromkeys(issues))
    @staticmethod
    def continuity(*steps, tol=0.25, quiet=False):
        """Atoms that persist between consecutive steps (same uid via Mol.copy) but move relative to the
        rest of the molecule. Each step is a Mol or a list of Mols. Returns a list of messages. Deliberate moves are
        not reported: atoms laid out by join(), swung by swing() or placed by dock()-style helpers (deprotonate,
        protonate, attack) in either of the two steps."""
        out = []
        def atoms(s): return {a.uid: a for m in (s if isinstance(s, (list, tuple)) else [s]) for a in m.atoms}
        def mols(s): return s if isinstance(s, (list, tuple)) else [s]
        for i in range(len(steps) - 1):
            A, B = atoms(steps[i]), atoms(steps[i + 1]); common = [u for u in A if u in B]
            if len(common) < 2: continue
            by = {}; moved = set().union(*[getattr(m, "_moved", set()) for m in mols(steps[i]) + mols(steps[i + 1])])
            for u in common:
                if u not in moved: by.setdefault((id(A[u].mol), id(B[u].mol)), []).append(u)
            for grp in by.values():
                d = np.array([B[u].xy - A[u].xy for u in grp]); shift = np.median(d, axis=0)
                for u, dv in zip(grp, d):
                    mv = norm(dv - shift)
                    if mv > tol: out.append(f"step {i + 1}->{i + 2}: {A[u]!r} moved {mv:.2f} BU")
        if not quiet: print("continuity: OK" if not out else "continuity:\n  " + "\n  ".join(out))
        return out
    def review(self, out="review.png", dpi=900):
        """Close-up sheet of every molecule and every curly arrow, for the mandatory visual check."""
        return self.zoom(out, dpi, mols=True)
    def zoom(self, out="zoom.png", dpi=900, pad=10.0, mols=False):
        """Contact sheet of close-ups of every curly arrow (and with mols=True every molecule) for visual QA."""
        from PIL import Image, ImageDraw
        tmp = out.replace(".png", "_full.png"); self.save(tmp, dpi=dpi); _, _, X0, Y1 = self._last
        im = Image.open(tmp); s = dpi / 72; tiles = []
        if mols:
            for mi, m in enumerate(self.mols):
                x0, y0, x1, y1 = self._mol_box(m, curlies=False, prep=False); x0 -= pad; y0 -= pad; x1 += pad; y1 += pad
                t = im.crop((int((x0 - X0) * s), int((Y1 - y1) * s), int((x1 - X0) * s), int((Y1 - y0) * s)))
                t.thumbnail((520, 520)); tiles.append((f"mol {mi}", t))
        for i, c in enumerate(self.curlies):
            P = np.vstack([c._pts, c._head, c._ctrl]); x0, y0 = P.min(0) - pad; x1, y1 = P.max(0) + pad
            t = im.crop((int((x0 - X0) * s), int((Y1 - y1) * s), int((x1 - X0) * s), int((Y1 - y0) * s)))
            t.thumbnail((520, 520)); tiles.append((f"#{i}", t))
        if not tiles: return tmp
        W = 540; cols = 3; rows = (len(tiles) + cols - 1) // cols
        sheet = Image.new("RGB", (W * cols, 560 * rows), "white"); d = ImageDraw.Draw(sheet)
        for k, (lab, t) in enumerate(tiles):
            x, y = (k % cols) * W + 10, (k // cols) * 560 + 30; sheet.paste(t, (x, y)); d.text((x, y - 22), lab, fill="blue")
        sheet.save(out); return out


def side_by_side(original, drawn, out="compare.png", crop=None, height=1000, gap=40):
    """The original (sketch, photo, slide) and the redrawn figure side by side at one height, as a PNG, to compare them
    element by element. crop=(x0, y0, x1, y1), fractions of the original, zooms into one part of it. Returns out."""
    from PIL import Image
    def load(p):
        im = Image.open(p)
        if im.mode in ("RGBA", "LA", "P"):
            im = im.convert("RGBA"); bg = Image.new("RGBA", im.size, "white"); bg.alpha_composite(im); im = bg
        return im.convert("RGB")
    a, b = load(original), load(drawn)
    if crop:
        w, h = a.size; a = a.crop((int(crop[0] * w), int(crop[1] * h), int(crop[2] * w), int(crop[3] * h)))
    a, b = (im.resize((max(1, round(im.size[0] * height / im.size[1])), height), Image.LANCZOS) for im in (a, b))
    c = Image.new("RGB", (a.size[0] + gap + b.size[0], height), "white"); c.paste(a, (0, 0)); c.paste(b, (a.size[0] + gap, 0))
    c.save(out); return out
