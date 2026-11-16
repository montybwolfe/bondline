# Curly arrows and mechanism operations (teaching-figure pass): SN2 attack straight anti to the leaving group; a
# bond's electrons onto the neighbouring atom read as a 1,2-shift (the group moves with them; validate() agrees) unless
# a base takes that group (E1 onto the cation); such arrows start on the side facing their target and curl over the
# outside; a hook onto a crowded vertex leans away and curls back; a charge stays out of a hook's loop in every copied
# step; a caption leaves check() measuring the marks as drawn.
from bondline import *
import math

def deg(m, a, b): return round(math.degrees(math.atan2(*(m[b].xy - m[a].xy)[::-1])) % 360)

S = Scene(Style())
s = Mol(); c = s.atom("C", (0, 0), name="c", show=True); s.attach(c, 0, "Br", name="br"); s.attach(c, 100, "Me")
s.stereo(c, s.attach(c, 200, "Et"), "hash"); s.stereo(c, s.attach(c, 260, "H"), "wedge")
ho = hydroxide(); ho.key.name = "o"; s.attack(ho, s["c"]); S.curly(s["o"].lp(), s["c"])
# [1,2]-hydride shift: isobutyl cation to tert-butyl cation
h1 = Mol(); c1 = h1.atom("C", (0, 0), name="c1"); c2 = h1.attach(c1, -30, name="c2"); h1.attach(c2, 90, "H", name="h")
h1.attach(c2, -30); h1.attach(c2, -90); h1.attach(c1, 150); c1.charge = 1
a1 = S.curly((c2, h1["h"]), c1)
h2 = h1.copy(); h2.unbond(h2["c2"], h2["h"]); h2.join(h2["c1"], h2["h"]); h2["c1"].charge = 0; h2["c2"].charge = 1
S.validate(h1, h2)
# [1,2]-methyl shift onto a secondary cation, drawn into the cation as teaching figures do; the same atoms with a base taking the
# H make an E1 step: the C-H electrons stay and make the C=C
m = Mol(); p = m.atom("C", (0, 0), name="p"); p.charge = 1; m.attach(p, 150); q = m.attach(p, 30, name="q")
m.attach(q, 90, name="me"); m.attach(q, -30); m.attach(q, -90)
a2 = S.curly((m["q"], m["me"]), m["p"])
n = m.copy(); n.unbond(n["q"], n["me"]); n.join(n["p"], n["me"]); n["p"].charge = 0; n["q"].charge = 1
S.validate(m, n)
e = Mol(); p2 = e.atom("C", (0, 0), name="p"); p2.charge = 1; e.attach(p2, 150); q2 = e.attach(p2, 30, name="q")
hq = e.attach(q2, 90, "H", name="h"); e.attach(q2, -30); w = water(explicit=False); w.key.name = "w"; e.deprotonate(w, hq)
S.curly(e["w"].lp(), e["h"]); a3 = S.curly((e["q"], e["h"]), e["p"])
f = e.copy(); f.unbond(f["q"], f["h"]); f.join(f["w"], f["h"]); f["w"].charge = 1; f.find(f["p"], f["q"]).order = 2; f["p"].charge = 0
S.validate(e, f)
S.row([s, h1, m, e]); print(S.save("sn2_shifts.png", dpi=300)); S.check()
print("SN2: HO- at", deg(s, "c", "o"), "deg (Br at 0, Et at 200: the straight run in clears Et)")
print("hydride shift:", a1._class, "| tail on the C1 side:", bool(a1._T[0] < S.P(c2)[0]),
      "| curls over the top:", bool(a1._arc[:, 1].max() > S.P(h1["h"])[1]), "| methyl shift:", a2._class, "| E1:", a3._class)

S = Scene(Style())
# chloroform: the C-H electrons stay on C, and every quadrant but one is taken
k = Mol(); kc = k.atom("C", (0, 0), name="c")
for t, lab in ((90, "Cl"), (180, "Cl"), (270, "Cl")): k.attach(kc, t, lab)
kh = k.attach(kc, 0, "H", name="h"); a4 = S.curly((kc, kh), kc)
# a short O-H between two labels: the hook curls over and the charge keeps out of its loop, in both copies of the step;
# a caption does not move the charge check() measures
t1 = Mol(); cc = t1.atom("C", (0, 0), name="c"); o = t1.attach(cc, 90, "O", 2, name="o"); t1.attach(cc, 210); t1.attach(cc, -30, "OMe")
t1.attach(o, 30, "H", name="h", deliberate=True); o.charge = 1
t2 = t1.copy()
for t in (t1, t2): S.curly((t["h"], t["o"]), t["o"])
S.caption(t2, "<i>the same step, captioned</i>")
S.row([k, t1, S.arrow(kind="<->"), t2]); print(S.save("hooks.png", dpi=300)); S.check()
inloop = lambda t, c_: math.hypot(*(t["o"]._chg[0] - c_._arc.mean(0))) < math.hypot(*(S.P(t["o"]) - c_._arc.mean(0)))
print("hook on CHCl3:", a4._class, "| lands between H and the upper Cl:", 0 < math.degrees(math.atan2(*(a4._P - S.P(kc))[::-1])) < 90,
      "| O-H hooks: charges outside the loops:", [not inloop(t, c_) for t, c_ in zip((t1, t2), S.curlies[1:])],
      "| same place in both copies:", t1["o"]._chg[2] == t2["o"]._chg[2])

S = Scene(Style())
# cyclic arrays: benzene's Kekule arrows are rotated copies of one short arc; a 4-centre array (a beta-hydride
# elimination drawn U-shaped) has room inside for one arrow, the other runs outside, each turning with the flow
bz = Mol(); bz.ring(6, (0, 0), rot=90, kekule=True)
for b in [b for b in bz.bonds if b.order == 2]:
    nb = [x for x in bz.bonds if x is not b and x.order == 1 and b.b in (x.a, x.b)][0]; S.curly((b.a, b.b), (nb.a, nb.b))
be = from_smiles("[H]CC[*]", explicit_h=True, labels={3: "PdL_{3}"}, xy={1: (0, 0), 2: (1, 0), 0: (-0.5, 0.87), 3: (1.5, 0.87)})
S.curly((be[1], be[0]), be[3]); S.curly((be[2], be[3]), (be[1], be[2]))
S.row([bz, be], gap=1.5); print(S.save("arrays.png", dpi=300)); S.check()
sw = [round(c._sweep, 1) for c in S.curlies[:3]]; fr = [S._frame(c) for c in S.curlies[:3]]
loc = [tuple(round(float(v), 1) for v in ((c._T - o) @ ex, (c._T - o) @ ey, (c._P - o) @ ex, (c._P - o) @ ey)) for c, (o, ex, ey) in zip(S.curlies[:3], fr)]
print("benzene arrows alike:", len(set(sw)) == 1 and len(set(loc)) == 1, "| sweep", sw[0], "| beta-H arrows turn",
      ["ccw" if c._ccw else "cw" for c in S.curlies[3:]])

# released groups read as free species; a lone diatomic shows its carbon; a pinned charge turns with its molecule
t = Mol(); c = t.atom("C", (0, 0)); t.attach(c, 90, "O", 2); t.attach(c, 180); om = t.attach(c, 0, "OMe"); w = t.attach(c, 270, "OH2")
meo = t.release(om, c); h2o = t.release(w, c); om.charge = -1
k = from_smiles("CC(C)([O-])C#N"); cn = k.release(k[4], k[1])
h3 = Mol(); o = h3.atom("O", (0, 0)); o.charge = 1; o.charge_at = 90; hh = h3.attach(o, 150, "H"); h3.attach(o, 30, "H"); h3.attach(o, 270, "H"); h3.key = hh
al = from_smiles("CC=C"); al.protonate(h3, al[2], angle=30)
S = Scene(Style())
print("free species read:", [S._ldir(x) for x in (om, w)], "| cyanide C shown:", any(a.show for a in cn.atoms),
      "| hydronium charge pinned at", round(o.charge_at % 360), "deg (its O-H bonds turned by", round((math.degrees(math.atan2(*(hh.xy - o.xy)[::-1])) - 150) % 360), "deg)")

# layout of the mechanism operations: a fused ring goes on the free side; a departing cation and a tidied cross are
# trigonal; a linear centre taking a third bond opens to 120 deg (one branch flipped if the two would meet)
spread = lambda m, x: sorted(round((math.degrees(math.atan2(*(n.xy - x.xy)[::-1])) - math.degrees(math.atan2(*(m.nbrs(x)[0].xy - x.xy)[::-1]))) % 360) for n in m.nbrs(x))
r = Mol(); rr = r.ring(6, (0, 0), rot=90, labels={3: "O"}); tri = r.ring_fuse(rr[4], rr[5], 3)
e = from_smiles("CC(=O)OC(C)(C)C", abbrev=False); tb = e.release(e[4], e[3])
t = Mol(); c = t.atom("C", (0, 0)); t.attach(c, 90, "O", 2, name="o"); t.attach(c, 210, "Ph"); t.attach(c, 330, "H")
am = Mol(); am.attach(am.atom("R", (0, 0)), 0, "NH2", name="n"); t.attack(am, c); t.join(c, t["n"]); t["n"].charge = 1
t.find(c, t["o"]).order = 1; t.tidy(carbon=True)
d = from_smiles("C1CCC(CC1)N=C=NC1CCCCC1", align=(6, 8, 0)); ac = from_smiles("CC(=O)[O-]"); o = [a for a in ac.atoms if a.charge][0]
d.merge(ac); d.attack(o, d[7]); d.join(d[7], o)
print("fused ring outside:", math.hypot(*tri[0].xy) > 1.0, "| tBu+ angles", spread(tb, tb.atoms[0]),
      "| tidied cross", spread(t, c), "| carbodiimide C after join", spread(d, d[7]))

# an arrow from a bare H is refused by check(); one with its target as its source is refused at once
S = Scene(Style())
hb = Mol(); c1 = hb.atom("C", (0, 0)); c2 = hb.attach(c1, 0, order=2); hh = hb.atom("H", (0.5, 1.6)); S.curly(hh, (c1, c2))
S.row([hb]); S.save("from_h.png", dpi=300); print("arrow from a bare H:", [p for p in S.check(quiet=True) if "no electrons" in p] != [])
try: S.curly(c1, c1)
except ValueError as x: print("arrow onto its own source:", str(x)[:60])
