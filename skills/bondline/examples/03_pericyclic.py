# Pericyclic reaction (Diels-Alder). Draw the product ring first and make the reactants as a copy with the two new
# sigma bonds removed, so no atom moves between the two structures. The three arrows form one cyclic array: each starts
# on a pi bond and ends on the bond it makes (Mid for a new sigma bond), all turning the same way round the ring.
from bondline import *
S = Scene(Style())
P = Mol(); ring = P.ring(6, rot=30)                     # c0-c3 from the diene, c4 and c5 from the dienophile
for k, a in enumerate(ring): a.name = f"c{k}"
P.find(P["c1"], P["c2"]).order = 2; P.attach(P["c4"], None, "CO2Me")
R = P.copy(); R.unbond(R["c3"], R["c4"]); R.unbond(R["c5"], R["c0"])
for i, j, o in ((0, 1, 2), (1, 2, 1), (2, 3, 2), (4, 5, 2)): R.find(R[f"c{i}"], R[f"c{j}"]).order = o
S.curly((R["c0"], R["c1"]), (R["c1"], R["c2"]))
S.curly((R["c2"], R["c3"]), Mid(R["c3"], R["c4"]))
S.curly((R["c4"], R["c5"]), Mid(R["c5"], R["c0"]))
S.validate(R, P)
S.row([R, S.arrow(above="Δ"), P]); print(S.save("diels_alder.png")); S.check()
