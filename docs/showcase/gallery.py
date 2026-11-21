# Orbitals, conformations and an organometallic: a π bond as two p orbitals (larger on O), a carbocation's empty p
# orbital, a singlet carbene, the nitrogen lone pair of methylamine in its sp3 lobe, a chair with an equatorial
# methyl, anti butane in a Newman projection, and (η⁶-benzene)Cr(CO)₃.
from bondline import *
S = Scene(Style())
co = Mol(); c = co.atom("C", (0, 0)); o = co.atom("O", (1, 0)); co.bond(c, o)
co.attach(c, 150, style="hash"); co.attach(c, 210, style="bold")
S.orbital(c, "p", 90, shade="up"); S.orbital(o, "p", 90, shade="up", size=1.6); S.caption(co, "π (C=O)")
cat = Mol(); k = cat.atom("C", (0, 0)); k.charge = 1; cat.attach(k, 0, "Me"); cat.attach(k, 150, "Me", style="hash")
cat.attach(k, 210, "Me", style="wedge"); S.orbital(k, "p", 90, shade="down"); S.caption(cat, "<i>tert</i>-butyl cation")
cb = Mol(); q = cb.atom("C", (0, 0)); cb.attach(q, 150, "Cl", style="hash"); cb.attach(q, 210, "Cl", style="wedge")
S.orbital(q, "p", 90, shade="down"); S.orbital(q, "lobe", 0, electrons=2); S.caption(cb, "singlet CCl_2")
am = Mol(); n = am.atom("N", (0, 0)); am.attach(n, 210, "CH3", style="wedge"); am.attach(n, 330, "H", style="hash"); am.attach(n, 270, "H")
S.orbital(n, "lobe", 90, electrons=2, size=1.4, label="n_{N(sp^3)}"); S.caption(am, "methylamine")
S.row([co, cat, cb, am], gap=2.5)
eq = Mol(); at, ax, eqd = eq.chair(); eq.attach(at[0], eqd[0], "Me"); eq.attach(at[0], ax[0], "H")
S.caption(eq, "methylcyclohexane")
anti = Mol(); anti.newman(["Me", "H", "H"], ["Me", "H", "H"], back_rot=270); S.caption(anti, "butane, anti")
pst = Mol(); cr = pst.atom("Cr", (0, 0)); ar = pst.persp_ring(6, (0, 1.75), toward=cr)
pst.deloc(ar); pst.eta(cr, ar); pst.ligate(cr, "tet", left="CO", fr="CO", br="CO"); S.caption(pst, "(η^6-C_6H_6)Cr(CO)_3")
S.row([eq, anti, pst], gap=3.0, y=-5.5)
print(S.save("gallery.png", dpi=260)); S.check()
