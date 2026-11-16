# Teaching-figure pass: the generic base B:⁻ (stock base(), read as a generic symbol, not boron); lone ions
# keep spare lone pairs on the axes and a one-unit label's charge after its last letter (:Base⊖, Nu⊖, Na⊕); a fourth
# substituent on a chain carbon makes a gem pair (an explicit alpha-H next to the carbonyl), not a 60 deg squeeze; a dotted
# O...Li coordination does not upset valence; arrow text that rides into its arrow is reported.
from bondline import *
st = Style()

print("== enolisation of an isopropyl ketone by B:⁻ (alpha-H placed by attach(), base by deprotonate())")
S = Scene(st)
k = Mol(); ca = k.atom("C", (0, 0), name="ca"); cc = k.attach(ca, 30, name="cc"); k.attach(cc, 90, "O", 2, name="o")
k.attach(cc, -30, "R"); k.attach(ca, 150); k.attach(ca, None, "H", name="h"); k.attach(ca, None, name="me")
print("H at", round(ang(k["h"].xy - k["ca"].xy)) % 360, "deg, Me at", round(ang(k["me"].xy - k["ca"].xy)) % 360, "deg")
b = base(); k.merge(b); k.deprotonate(b.key, k["h"])
S.curly(b.key.lps[0], k["h"]); S.curly((k["ca"], k["h"]), (k["ca"], k["cc"])); S.curly((k["cc"], k["o"]), k["o"])
e = k.copy(); e.unbond(e["ca"], e["h"]); e.join(e[b.key], e["h"]); e[b.key].charge = 0
e.find(e["ca"], e["cc"]).order = 2; e.find(e["cc"], e["o"]).order = 1; e["o"].charge = -1
bh = e.split(e[b.key])
S.validate(k, [e, bh]); S.balance(k, [e, bh])
S.row([k, S.arrow(above=base(), kind="<=>"), e, "+", bh])
print(S.save("enolisation.png", dpi=600)); S.check()

print("== lone ions: spare pairs on the axes, a one-unit label's charge after it, a radical by its atom")
S = Scene(st)
def ion(l, q=0, rad=False, lp=0):
    m = Mol(); a = m.atom(l, (0, 0)); a.charge = q; a.radical = rad
    for _ in range(lp): a.lp()
    return m
S.row([base(), base("Base"), base("Nu"), ion("Na", 1), ion("Br", -1, lp=4), ion("Br", rad=True), alkoxide("Et"), ion("NMe4", 1)], gap=1.0)
print(S.save("ions.png", dpi=600)); S.check()

print("== a gem pair on a chain vertex, either order; OH then Me on a zigzag")
S = Scene(st)
g1 = Mol(); c = g1.atom("C", (0, 0)); g1.attach(c, 150); g1.attach(c, 30); g1.attach(c, None, "OH"); g1.attach(c, None, "Me")
g2 = Mol(); c = g2.atom("C", (0, 0)); x = g2.attach(c, 150); g2.attach(x, 210); g2.attach(c, 30, "CO2Et"); g2.attach(c, None); g2.attach(c, None)
S.row([g1, g2], gap=1.5); print(S.save("gem.png", dpi=600)); S.check()

print("== lithium enolate formation through a cyclic arrangement: O...Li and N...H dotted, no valence complaint")
S = Scene(st)
t = Mol(); r = t.ring(6, rot=90, labels={0: "Li", 1: "O", 3: "C", 4: "H", 5: "N"})
for i, a in enumerate(r): a.name = f"r{i}"
t.find(t["r0"], t["r1"]).style = "dot"; t.find(t["r4"], t["r5"]).style = "dot"
t.find(t["r1"], t["r2"]).order = 2; t.attach(t["r2"], None, "R"); t.attach(t["r5"], None, "iPr"); t.attach(t["r5"], None, "iPr")
S.curly(t["r5"].lp(), t["r4"]); S.curly((t["r3"], t["r4"]), (t["r2"], t["r3"])); S.curly((t["r1"], t["r2"]), t["r1"])
S.row([t]); print(S.save("cyclic.png", dpi=600)); S.check()

print("== negative: superscripts under an arrow, spaced from the cap height (the old way), are reported")
S = Scene(st)
S._ink_shift = lambda *a, **k: 0.0
S.row([from_smiles("CC=O"), S.arrow(below="[Ru(bpy)3]2+, NH4+, PF6-"), from_smiles("CCO")]); S.check()

print("== reversed condensed labels anchored on carbon or nitrogen (typed with the bond on the right)")
S = Scene(st); ms = []
for l in ("R2CH", "R2C", "EtO2C", "HO2C", "O2N", "PhCH2", "HOCH2", "NC"):
    m = Mol(); c = m.atom("C", (0, 0), show=True); o = m.attach(c, 90, "O", 1 if l == "R2C" else 2)
    o.charge = -1 if l == "R2C" else 0; m.attach(c, 210, l, 2 if l == "R2C" else 1); m.attach(c, -30, "Me")
    print(f"{l:6} -> {m.atoms[2].label}"); ms.append(m)
S.row(ms, gap=1.0); print(S.save("reversed.png", dpi=600)); S.check()
