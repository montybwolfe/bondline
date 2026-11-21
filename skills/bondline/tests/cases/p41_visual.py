# Visual pass (fixer "visual"): orbital lobes as soft textbook teardrops, a lobe on a labelled atom coming out from
# behind its letter with its outline open (no pinch beside the letter), electrons apart in the full part of a lobe, a
# carbocation's charge beside its empty p orbital; bond lines unsnapped, so an upright bond meets its ring joins clean.
import numpy as np
from bondline import *
from bondline import _segs_dist

print("== a p orbital: two soft teardrops about 0.6 as wide as long, meeting at a bare carbon's nucleus")
S = Scene(Style())
m = Mol(); c = m.atom("C", (0, 0)); o = m.atom("O", (1, 0)); m.bond(c, o)
m.attach(c, 150, style="hash"); m.attach(c, 210, style="bold")
S.orbital(c, "p", 90, shade="up"); S.orbital(o, "p", 90, shade="up", size=1.4)
k = Mol(); x = k.atom("C", (0, 0)); x.charge = 1; k.attach(x, 0, "R"); k.attach(x, 150, "R", style="hash")
k.attach(x, 210, "R", style="wedge"); S.orbital(x, "p", 90, shade="down")
S.row([m, k], gap=2.0); print(S.save("pi_cation.png", dpi=300)); S.check()
for ex in S.extras:
    if ex[0] != "orbital": continue
    a = ex[1]; A = S.P(a); u = dirv(ex[3]); P = [np.asarray(p[1], float) for p in S._orbital_prims(ex) if p[0] == "orb"]
    up = P[0] - A; along = up @ u; across = up @ perp(u)
    print(a, "lobe", round(float(along.max() - along.min()), 1), "pt long,", round(float(across.max() - across.min()), 1),
          "pt wide (", round(float((across.max() - across.min()) / along.max()), 2), "of its reach ); starts",
          round(float(along.min()), 1), "pt out;", "open at the label" if S._orbital_prims(ex)[0][4:5] == ("open",) else "closed at the nucleus")
cut = [p for ex in S.extras if ex[0] == "orbital" and ex[1] is o for p in S._orbital_prims(ex) if p[0] == "orb"][0]
P = np.asarray(cut[1], float); print("O lobe's open end:", round(float(norm(P[0] - P[-1])), 2), "pt across (line width", S.st.lw, ")")
print("C+ charge at", round(x._chg[2]), "deg,", round(float(norm(x._chg[0] - S.P(x)) / S.st.bond), 2), "BU out; beside its lobe:",
      S._by_own_lobe(x, x._chg[0], x._chg[1]), "; nearer its lobe than R:",
      bool(S._lobe_dist(x, x._chg[0]) < norm(x._chg[0] - S.P(k.atoms[1]))))

print("== a lobe with a bond of its atom within 50 deg of its axis is slimmer, so the two never run side by side")
S = Scene(Style())
d = from_smiles("CS(=O)C", align=(1, 2, 0)); S.orbital(d[1], "p", 90, shade="up", size=1.6)
S.row([d]); print(S.save("slim.png", dpi=300)); S.check()
for ex in S.extras:
    a = ex[1]; P = np.vstack([np.asarray(p[1], float) for p in S._orbital_prims(ex) if p[0] == "orb"])
    L = np.array([[*p0, *p1] for b in d.bonds if b.order == 1 for p0, p1, e in S._bond_lines(b)])
    print("S lobes (size", ex[6], ") form", tuple(round(float(v), 2) for v in S._lobe_form(a, 90)), "; white to the methyl bonds",
          round(float(_segs_dist(P, L).min()) - S.st.lw, 2), "pt")

print("== electrons: two dots apart in the full part of a hybrid lobe, on a labelled atom and on a bare carbanion")
S = Scene(Style())
am = Mol(); n = am.atom("N", (0, 0)); am.attach(n, 200, "Me", style="wedge"); am.attach(n, -20, "H", style="hash"); am.attach(n, 260, "H")
S.orbital(n, "lobe", 90, electrons=2)
cb = Mol(); h = cb.atom("C", (0, 0)); cb.attach(h, 210); cb.attach(h, 330); h.charge = -1; S.orbital(h, "lobe", 90, electrons=2)
S.row([am, cb], gap=2.0); print(S.save("pairs.png", dpi=300)); S.check()
for ex in S.extras:
    pr = S._orbital_prims(ex); d = [np.asarray(p[1], float) for p in pr if p[0] == "dot"]
    lob = np.asarray([p for p in pr if p[0] == "orb"][0][1], float)
    print(ex[1], "lobes", sum(1 for p in pr if p[0] == "orb"), "dots", round(float(norm(d[0] - d[1])), 2), "pt apart,",
          round(float(min(np.hypot(*(lob - q).T).min() for q in d)) - S.st.lw / 2 - 0.8 * S.st.lp_r, 2), "pt inside the outline")

print("== an upright bond meets its ring joins without a step (lines drawn unsnapped)")
S = Scene(Style()); r = Mol(); r.ring(6, (0, 0), rot=90, kekule=True); S.row([r]); print(S.save("ring.png", dpi=300)); S.check()

print("== dotted bonds (H-bonds): a solid bond's white at a label, dots evenly spaced and centred, alike when alike")
S = Scene(Style())
m = Mol(); c = m.atom("C", (0, 0)); m.attach(c, 90, "O", 2); n1 = m.attach(c, 210, "N"); n2 = m.attach(c, 330, "N")
m.attach(n1, 150); m.attach(n2, 30); h1 = m.attach(n1, 270, "H"); h2 = m.attach(n2, 270, "H")
y = m.atom("O", (0, -2.3)); y.charge = -1; m.attach(y, 270, "Me")
m.bond(h1, y, style="dot", free=True); m.bond(h2, y, style="dot", free=True)
S.row([m]); print(S.save("hbonds.png", dpi=300)); S.check()
for b in [b for b in m.bonds if b.style == "dot"]:
    P = np.asarray(b._prims[0][1], float); u = unit(P[1] - P[0]); on, off = b._prims[0][4]
    n = int(round((norm(P[1] - P[0]) + off) / (on + off)))
    ends = [round(float((P[0] - S.P(b.a)) @ u) - S.st.lw / 2 - S._exit(b.a, u), 2), round(float((S.P(b.b) - P[1]) @ u) - S.st.lw / 2 - S._exit(b.b, -u), 2)]
    print(b.a, "...", b.b, n, "dots", round(on + off, 2), "pt apart; white past a solid bond's end at each label:", ends)

print("== two H-bonds to one Cl-: as many dots on each though one meets the C and the other the l, all clear of the letters")
from bondline import _hull_dist
S = Scene(Style())
m = Mol(); c = m.atom("C", (0, 0)); m.attach(c, 90, "S", 2); n1 = m.attach(c, 210, "N"); n2 = m.attach(c, 330, "N")
m.attach(n1, 150); m.attach(n2, 30); h1 = m.attach(n1, 270, "H"); h2 = m.attach(n2, 270, "H")
x = m.atom("Cl", (0, -2.45)); x.charge = -1
m.bond(h1, x, style="dot", free=True); m.bond(h2, x, style="dot", free=True)
S.row([m]); print(S.save("hbonds_cl.png", dpi=300)); S.check()
for b in [b for b in m.bonds if b.style == "dot"]:
    P = np.asarray(b._prims[0][1], float); u = unit(P[1] - P[0]); on, off = b._prims[0][4]
    n = int(round((norm(P[1] - P[0]) + off) / (on + off))); D = [P[0] + u * (i * (on + off) + on / 2) for i in range(n)]
    wh = [round(min(_hull_dist(d, h) for d in D for g in a._g for h in g["hulls"]) - S.st.lw / 2, 2) for a in (b.a, b.b)]
    print(b.a, "...", b.b, n, "dots", round(on + off, 2), "pt apart; white from the dots to the letters:", wh)
print("== free text keeps a third of a bond from structures it is not part of: a callout's note clear of all ink; check()")
S = Scene(Style())
k = from_smiles("O=C(OC)c1ccc(F)cc1", align=(1, 0, 90))
S.pointer([k[2]], "ester O", side=-30, bend=-40); S.row([k]); print(S.save("note.png", dpi=300)); S.check()
nb = S._ptr_geo(S.pointers[0])["texts"][0]["box"]
L = [(p0, p1) for b in k.bonds for p0, p1, ex in S._bond_lines(b)]
w = min([seg_box_gap(p0, p1, nb) - S.st.lw / 2 for p0, p1 in L] + [box_gap(nb, g["box"]) for x in k.atoms for g in (x._g or [])])
print("the note's white to the ink:", round(w, 1), "pt (at least", S.TEXT_CLEAR, ")")
S = Scene(Style()); k = from_smiles("CC(=O)C"); S.row([k]); S.text("a free note", at=(k[2].xy[0] + 0.5, k[2].xy[1]), ha="left"); S.check()
