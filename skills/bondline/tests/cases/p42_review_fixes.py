# 2.0.1 review fixes: two wedges in a row fill as one shape; single bonds meet a centred C=O at its near line (no
# line running into the gap between its two lines); a lobe on a labelled atom shows from just off the letter; cycle
# arcs end near a station's outline, not its box; a callout's arrow comes in clear of curly arrows at its target;
# a terminal OH crowding another label reads another way.
import numpy as np
from bondline import *
from bondline import box_gap, _hull, _hull_dist, _segs_dist

print("== two wedges in a row: the first's wide end runs on into the second's edges")
S = Scene(Style())
m = Mol(); a = m.atom("C", (0, 0)); b = m.atom("C", (0.866, 0.5)); c = m.atom("C", (1.732, 0))
w1 = m.bond(a, b, style="wedge"); w2 = m.bond(b, c, style="wedge"); m.attach(a, 210); m.attach(a, 270); m.attach(b, 90); m.attach(c, 330)
S.row([m]); print(S.save("wedge_row.png", dpi=1200)); S.check()
P1 = np.asarray(w1._prims[0][1], float)
print("first wedge's wide end joined to the second (outline points):", len(P1), "(a square end has 4)")

print("== a centred C=O: the single bonds stop at its near lines, meeting them in a corner")
S = Scene(Style())
acid = from_smiles("CC(=O)O"); ket = from_smiles("CC(C)=O")
S.row([acid, ket], gap=1.5); print(S.save("centred_co.png", dpi=1200)); S.check()
for mol in (acid, ket):
    C = [a for a in mol.atoms if any(mol.find(a, n).order == 2 for n in mol.nbrs(a))][0]
    for n in mol.nbrs(C):
        bd = mol.find(C, n)
        if bd.order != 1: continue
        P = np.asarray(bd._prims[0][1], float); e = min(P, key=lambda q: np.linalg.norm(q - S.P(C)))
        dl = [np.asarray(pr[1], float) for pr in mol.find(C, [k for k in mol.nbrs(C) if mol.find(C, k).order == 2][0])._prims]
        meets = min(np.linalg.norm(e - L[i]) for L in dl for i in (0, -1))
        print(C, "-", n, "stops", round(float(np.linalg.norm(e - S.P(C))), 2), "pt short of the carbon; a C=O line ends there:", meets < 0.05)

print("== lobes on a labelled atom show from just off the letter")
S = Scene(Style())
m = Mol(); c = m.atom("C", (0, 0)); n = m.atom("N", (1, 0)); m.bond(c, n, order=2); m.attach(c, 150); m.attach(c, 210); m.attach(n, 300, "H")
S.orbital(c, "p", 90, shade="up"); S.orbital(n, "p", 90, shade="up")
S.row([m]); print(S.save("lobe_letter.png", dpi=1200)); S.check()
for ex in S.extras:
    if ex[1] is not n: continue
    for pr in S._orbital_prims(ex):
        if pr[0] != "orb": continue
        P = np.asarray(pr[1], float); H = [hh for g in n._g for hh in g["hulls"]]
        print("lobe outline's end to the N's ink:", round(min(_hull_dist(q, hh) for q in (P[0], P[-1]) for hh in H), 2), "pt")

print("== cycle arcs end the gap from a station's outline (its letters and bond ends), not from its box")
S = Scene(Style())
def pd(e, w):
    q = Mol(); p = q.atom("Pd", (0, 0)); q.ligate(p, "sqpl", N="L", S="L", E=e, W=w); return q
pd0 = Mol(); p = pd0.atom("Pd^0", (0, 0)); pd0.ligate(p, "lin", E="L", W="L")
cy = S.cycle([(pd0, 90), (pd("Br", "Ar"), -30), (pd("Ar'", "Ar"), 210)], radius=4)
print(S.save("cycle_ends.png", dpi=300)); S.check()
boxes = [x["box"] for q in (pd0,) + tuple(it.mols[0] if hasattr(it, "mols") else it for it, _ in cy["its"][1:]) for at in q.atoms for x in (at._g or [])]
ends = [np.asarray(pr[1], float)[k] for pr in S._cycle_prims(cy, []) if pr[0] == "line" and len(pr[1]) > 10 for k in (0, -1)]
print("arc ends to the nearest letter:", sorted(round(min(box_gap((*e, *e), bx) for bx in boxes), 1) for e in ends)[:3], "pt (gap",
      round(cy["gap"] * S.st.bond, 1), ")")

print("== a callout's arrow comes in clear of a curly arrow ending at the same atom")
S = Scene(Style())
k = from_smiles("C1=[O+]CCc2ccccc12", align=(0, 1, 30)); c1, o1 = k[0], k[1]
nu = from_smiles("*C(*)=C1CCCC1", labels={0: "OTMS", 2: "OMe"}, align=(1, 3, 150, 2)); cb, ca, ot = nu[3], nu[1], nu[0]
k.merge(nu); k.dock(cb, c1, -90, 1.8)
S.curly(ot.lp(), (ot, ca)); S.curly((ca, cb), c1); S.curly((c1, o1), o1)
S.pointer([o1], "oxocarbenium\nformation", side=-20, bend=-40)
S.row([k]); print(S.save("callout_curly.png", dpi=300)); S.check()
g = S._ptr_geo(S.pointers[0]); P = np.vstack([np.asarray(pr[1], float) for pr in g["prims"]])
C = np.vstack([np.vstack([c._pts, c._head]) for c in S.curlies])
print("callout arrow to the curly arrows:", round(float(_segs_dist(P, np.hstack([C[:-1], C[1:]])).min()) - S.st.lw, 1) >= 2.0)

print("== a terminal OH crowding a C=O label reads another way")
S = Scene(Style())
m = Mol(); c1 = m.atom("C", (0, 0)); c2 = m.atom("C", (0, -1)); m.bond(c1, c2); m.attach(c1, 120); m.attach(c1, 210)
o1 = m.atom("OH", (1.0, 0.0)); m.bond(c1, o1, style="hash")
oe = m.atom("O", (0.87, -1.5)); m.bond(c2, oe, style="wedge"); m.attach(c2, 210)
ca = m.atom("C", (1.73, -1.0)); m.bond(oe, ca); oc = m.atom("O", (1.9, -0.02)); m.bond(ca, oc, order=2); m.attach(ca, -30)
S.row([m]); print(S.save("oh_untangled.png", dpi=600)); S.check()
print("OH reads", o1.ldir, "; white to the C=O's O:", round(min(box_gap(g1["box"], g2["box"]) for g1 in o1._g for g2 in oc._g), 1), "pt")
