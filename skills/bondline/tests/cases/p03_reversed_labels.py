# P0-3: labels typed the way they read with the bond on the right (HO, H2N, H3C, MeO, EtO, AcO, BrMg, TsO, F3C, Me2N)
# are stored anchor-first, so bonds, lone pairs, charges and radicals attach to the right letter; a lone species keeps
# the order it was typed in. check() names any label that still anchors on H or a group.
from bondline import *
st = Style()
print("== substituents typed reversed, bonds from either side")
S = Scene(st); mols = []
for lab in ("HO", "H2N", "H3C", "MeO", "EtO", "AcO", "TsO", "F3C", "Me2N", "BocHN", "TBSO"):
    m = Mol(); c = m.atom("C", (0, 0)); c2 = m.attach(c, 30); x = m.attach(c2, -30, lab); y = m.attach(c, 150, lab)
    print(f"{lab:6s} -> {x.label}"); mols.append(m)
S.row(mols[:6]); S.save("substituents1.png")
S2 = Scene(st); S2.row(mols[6:]); S2.save("substituents2.png"); S.check(); S2.check()
print("== lone species: charges, lone pairs and radicals on the heteroatom or carbon")
S = Scene(st)
oh = Mol(); o = oh.atom("HO", (0, 0)); o.charge = -1; o.add_lps(3)
eo = Mol(); e = eo.atom("EtO", (0, 0)); e.charge = -1; e.add_lps(2)
me = Mol(); r = me.atom("H3C", (0, 0)); r.radical = True
gr = Mol(); g = gr.atom("BrMg", (0, 0)); c = gr.attach(g, 0, "Me", length=1.3, free=True)
w = Mol(); h2o = w.atom("H2O", (0, 0))
S.row([oh, eo, me, gr, w]); print(S.save("species.png")); S.check()
for a in (o, e, r, g, h2o): print(repr(a), "anchor", label_units(a.label)[0])
print("== negative: a label that still anchors on H; valence read from substituents in labels")
S = Scene(st)
m = Mol(); c = m.atom("C", (0, 0)); m.attach(c, 0, "HOCH2"); n = m.attach(c, 120, "OEt"); m.attach(n, 90, "Me")
k = Mol(); c = k.atom("C", (0, 0)); nn = k.attach(c, 0, "NMe2"); k.attach(nn, 60, "Me")
S.row([m, k]); S.check()
