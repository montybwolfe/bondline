# Radical chain (bromination of methane): one row per step, fishhooks (S.fish) for single electrons. A radical that
# abstracts an atom is docked in line with the bond it attacks; Mid() is the middle of the bond being made. Br-Br is
# drawn 1.3 BU long (free=True) so the bond between the two labels shows; the step names stand in a column at left.
from bondline import *
S = Scene(Style())
def bromine():
    m = Mol(); a = m.atom("Br", (0, 0), name="Br1"); m.attach(a, 0, "Br", name="Br2", length=1.3, free=True); return m
# initiation
i0 = bromine(); S.fish((i0["Br1"], i0["Br2"]), i0["Br1"]); S.fish((i0["Br1"], i0["Br2"]), i0["Br2"])
i1 = i0.copy(); i1.unbond(i1["Br1"], i1["Br2"]); b2 = i1.split(i1["Br2"]); i1["Br1"].radical = b2["Br2"].radical = True
S.validate(i0, [i1, b2])
S.row([i0, S.arrow(above="hν"), i1, "+", b2], y=0); S.text("initiation", at=(-1.2, 0), ha="right")
# propagation 1: Br* takes H from methane
p0 = Mol(); c = p0.atom("CH3", (0, 0), name="C"); h = p0.attach(c, 0, "H", name="H")
br = p0.atom("Br", (0, 0), name="Br", radical=True); p0.dock(br, h, 0)          # in line with the C-H bond
S.fish(br, Mid(h, br)); S.fish((c, h), Mid(h, br)); S.fish((c, h), c)
p1 = p0.copy(); p1.unbond(p1["C"], p1["H"]); p1.join(p1["Br"], p1["H"]); p1["Br"].radical = False; p1["C"].radical = True
hbr = p1.split(p1["Br"])
S.validate(p0, [p1, hbr])
S.row([p0, S.arrow(), p1, "+", hbr], y=-3.2); S.text("propagation", at=(-1.2, -3.2), ha="right")
# propagation 2: the methyl radical takes Br from Br2, giving back Br*
q0 = Mol(); c = q0.atom("CH3", (0, 0), name="C", radical=True); q0.merge(bromine()); q0.dock(q0["Br1"], c, 0)
S.fish(c, Mid(c, q0["Br1"])); S.fish((q0["Br1"], q0["Br2"]), Mid(c, q0["Br1"])); S.fish((q0["Br1"], q0["Br2"]), q0["Br2"])
q1 = q0.copy(); q1.unbond(q1["Br1"], q1["Br2"]); q1.join(q1["C"], q1["Br1"]); q1["C"].radical = False
q1["Br2"].radical = True; brr = q1.split(q1["Br2"])
S.validate(q0, [q1, brr])
S.row([q0, S.arrow(), q1, "+", brr], y=-6.4)
print(S.save("radical_chain.png")); S.check()
