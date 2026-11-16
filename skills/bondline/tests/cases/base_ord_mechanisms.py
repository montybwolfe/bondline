# Ordinary mechanisms: SN1 ionisation, Grignard addition, radical chain step, resonance, TS brackets
from bondline import *
st = Style()
print("== SN1 ionisation")
S = Scene(st)
m = Mol(); c = m.atom("C", (0, 0), name="c"); [m.attach(c, t, "Me") for t in (90, 180, 270)]
br = m.attach(c, 0, "Br", name="br")
S.curly((c, br), br)
p = m.copy(); p.unbond(p["c"], p["br"]); p["c"].charge = 1; p["br"].charge = -1
S.validate(m, p)
S.row([m, S.arrow(), p]); print(S.save("sn1.png")); S.check()
print("== Grignard addition")
S = Scene(st)
m = Mol(); c = m.atom("C", (0, 0), name="c"); m.attach(c, 150, "Me"); m.attach(c, 30, "Me"); O = m.attach(c, -90, "O", 2, name="O")
C = m.atom("C", (0, 0), name="Cg"); Mg = m.attach(C, 90, "MgBr", name="Mg"); m.dock(C, c, 90)
S.curly((C, Mg), c); S.curly((c, O), O)
p = m.copy(); p.unbond(p["Cg"], p["Mg"]); p.bond(p["Cg"], p["c"]); p.find(p["c"], p["O"]).order = 1; p["O"].charge = -1; p["Mg"].charge = 1
S.validate(m, p)
S.row([m, S.arrow(above="Et_2O"), p]); print(S.save("grignard.png")); S.check()
print("== radical abstraction")
S = Scene(st)
m = Mol(); C = m.atom("C", (0, 0), name="C"); [m.attach(C, t, "Me") for t in (90, 210)]; H = m.attach(C, -30, "H", name="H")
Br = m.atom("Br", (0, 0), name="Br"); Br.radical = True; m.dock(Br, H, -30)
S.fish(Br, Mid(H, Br)); S.fish((C, H), Mid(H, Br)); S.fish((C, H), C)
p = m.copy(); p.unbond(p["C"], p["H"]); p.bond(p["H"], p["Br"]); p["Br"].radical = False; p["C"].radical = True
S.validate(m, p)
S.row([m, S.arrow(), p]); print(S.save("radical.png")); S.check()
print("== allyl resonance")
S = Scene(st)
a = Mol(); c1 = a.atom("C", (0, 0), name="c1"); c2 = a.attach(c1, 30, name="c2", order=2); c3 = a.attach(c2, -30, name="c3"); c3.charge = 1
S.curly((c1, c2), (c2, c3))
b = a.copy(); b.find(b["c1"], b["c2"]).order = 1; b.find(b["c2"], b["c3"]).order = 2; b["c3"].charge = 0; b["c1"].charge = 1
S.validate(a, b)
S.brackets([a, b])
S.row([a, S.arrow(kind="<->"), b]); print(S.save("allyl.png")); S.check()
print("== SN2 transition state")
S = Scene(st)
t = Mol(); c = t.atom("C", (0, 0)); t.attach(c, 90, "H"); t.attach(c, 210, "H"); t.attach(c, 330, "H")
nu = t.attach(c, 180, "HO", style="dash", length=1.4, free=True); lg = t.attach(c, 0, "Br", style="dash", length=1.4, free=True)
nu.partial = "-"; lg.partial = "-"
S.brackets([t], charge="‡")
S.row([t]); print(S.save("ts.png")); S.check()
