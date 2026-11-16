# Detail gallery: every kind of element the engine draws (marks on labels and bare carbons, every bond style and join,
# labels, reaction and curly arrows of every kind, brackets), for the zoomed element-by-element review that every
# release gets (SKILL.md, review checklist). Each is drawn the way the instructions say to draw it.
from bondline import *
st = Style()

print("== marks: charges, radicals, lone pairs, partial charges, numbers")
S = Scene(st); ms = []
m = Mol(); n = m.atom("N", (0, 0)); [m.attach(n, t, "Me") for t in (30, 150, 270)]; m.attach(n, 90, "H"); n.charge = 1; ms.append(m)   # ammonium
m = Mol(); o = m.atom("O", (0, 0)); m.attach(o, 210); o.charge = -1; o.lp(); o.lp(); o.lp(); ms.append(m)                         # alkoxide with 3 LPs
m = Mol(); c = m.atom("C", (0, 0)); [m.attach(c, t) for t in (90, 210, 330)]; c.charge = 1; ms.append(m)                         # bare carbocation
m = Mol(); c = m.atom("C", (0, 0)); [m.attach(c, t) for t in (150, 30)]; c.charge = -1; c.lp(); ms.append(m)                     # bare carbanion
m = Mol(); c = m.atom("C", (0, 0)); [m.attach(c, t) for t in (90, 210, 330)]; c.radical = True; ms.append(m)                     # bare radical
m = Mol(); b = m.atom("Br", (0, 0)); b.radical = True; ms.append(m)
m = Mol(); c = m.atom("C", (0, 0)); o = m.attach(c, 30, "O", 2); m.attach(c, 150); m.attach(c, 270); c.partial = "+"; o.partial = "-"; ms.append(m)
m = Mol(); r = m.ring(6, kekule=True)
for k, a in enumerate(r): a.number = k + 1
ms.append(m)
m = Mol(); x = m.atom("NH3", (0, 0)); x.charge = 1; m.attach(x, 330); ms.append(m)                                                # labelled charge
m = Mol(); x = m.atom("O", (0, 0)); m.attach(x, 210); m.attach(x, 330); m.attach(x, 90, "H"); x.charge = 1; x.lp(); ms.append(m)   # oxonium
S.flow(ms, width_cm=17); print(S.save("marks.png")); S.check()

print("== bonds")
S = Scene(st); ms = []
ms.append(from_smiles("CC=CC")); ms.append(from_smiles("CC#CC"))
m = Mol(); a = m.atom("C", (0, 0)); m.attach(a, 30, "O", 2); m.attach(a, 150); m.attach(a, 270, "OH"); ms.append(m)
m = Mol(); a = m.atom("C", (0, 0)); m.attach(a, 90, style="wedge"); m.attach(a, 210, style="hash"); m.attach(a, 330, "Br"); ms.append(m)
m = Mol(); a = m.atom("C", (0, 0)); b = m.attach(a, 30); m.attach(b, -30, style="wavy"); m.attach(a, 150, style="bold"); ms.append(m)
m = Mol(); a = m.atom("C", (0, 0)); b = m.attach(a, 30, style="dash"); c = m.attach(b, -30, style="partial", order=2); ms.append(m)
m = Mol(); n = m.atom("N", (0, 0)); [m.attach(n, t, "Me") for t in (90, 210, 270)]; bb = m.attach(n, 0, "B", style="dative", length=1.6, free=True)
[m.attach(bb, t, "F") for t in (90, 330, 270)]; ms.append(m)
m = Mol(); r = m.ring(6); m.deloc(r); ms.append(m)
m = from_smiles("c1ccc2ccccc2c1"); ms.append(m)
m = from_smiles("O=C1CCCC1"); ms.append(m)
S.flow(ms, width_cm=17); print(S.save("bonds.png")); S.check()

print("== labels")
S = Scene(st); ms = []
for smi in ("OCC(=O)OC", "NC1=CC=CC=C1", "CN(C)C(F)(F)F", "c1ccncc1", "c1cc[nH]c1", "CS(=O)(=O)Oc1ccccc1", "[NH3+]CC([O-])=O", "C[Si](C)(C)OC"):
    ms.append(from_smiles(smi, abbrev=False))          # every atom drawn: this gallery is about the labels themselves
S.flow(ms, width_cm=17); print(S.save("labels.png")); S.check()

print("== reaction arrows and text")
S = Scene(st)
a = from_smiles("CC=O"); b = from_smiles("CC(O)C#N"); c = from_smiles("CC(O)C(=O)O")
S.flow([a, S.arrow(above=["HCN, KCN (cat.)", "H_2O, 0 °C"], below="1 h, 92%"), b, S.arrow(kind="<=>", above="H_3O^+", below="Δ"), c,
        S.arrow(kind="<->"), from_smiles("CC([O-])C(=O)O"), S.arrow(kind="=>", above="retro")], width_cm=17)
print(S.save("arrows.png")); S.check()

print("== curly arrows of every class")
S = Scene(st)
m = Mol(); c = m.atom("C", (0, 0), name="C"); o = m.attach(c, 90, "O", 2, name="O"); m.attach(c, 210); m.attach(c, 330)
nu = m.attack(hydroxide(), c)
S.curly(nu.lp(), c); S.curly((c, o), o)
m2 = Mol(); c1 = m2.atom("C", (0, 0)); c2 = m2.attach(c1, 30, order=2); c3 = m2.attach(c2, -30); x = m2.attach(c3, 30, "Br")
S.curly((c1, c2), (c2, c3)); S.curly((c3, x), x); 
m3 = Mol(); a1 = m3.atom("O", (0, 0)); a1.charge = -1; a2 = m3.attach(a1, 30); a3 = m3.attach(a2, -30, order=2)
S.curly(a1.lp(), (a1, a2)); S.curly((a2, a3), a3)
S.row([m, m2, m3]); print(S.save("curly.png")); S.check()

print("== brackets, TS, resonance")
S = Scene(st)
t = Mol(); c = t.atom("C", (0, 0)); t.attach(c, 90, "H"); t.attach(c, 240, "H", style="wedge"); t.attach(c, 300, "H", style="hash")
nu = t.attach(c, 180, "HO", style="dash", length=1.4, free=True); lg = t.attach(c, 0, "Br", style="dash", length=1.4, free=True)
nu.partial = "-"; lg.partial = "-"; S.brackets([t], charge="‡")
S.row([t]); print(S.save("ts.png")); S.check()
