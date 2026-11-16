# Organometallics: sandwich() draws metallocenes (stagger, n=6 arenes, bend= for bent metallocenes), ligate() adds
# ligands (the "bent" sites up/mid/down sit in a bent metallocene's opening), persp_ring(toward=) + eta() builds a
# half-sandwich by hand, and deloc() draws the delocalised ring.
from bondline import *
S = Scene(Style())
fc = Mol(); fc.sandwich("Fe")                                                     # ferrocene
cp2 = Mol(); ti, _ = cp2.sandwich("Ti", bend=130); cp2.ligate(ti, "bent", up="Cl", down="Cl")   # titanocene dichloride
pst = Mol(); cr = pst.atom("Cr", (0, 0)); ar = pst.persp_ring(6, (0, 1.75), toward=cr)
pst.deloc(ar); pst.eta(cr, ar); pst.ligate(cr, "tet", left="CO", fr="CO", br="CO")  # (benzene)Cr(CO)3
for m, name in ((fc, "ferrocene"), (cp2, "Cp_2TiCl_2"), (pst, "(η^6-C_6H_6)Cr(CO)_3")): S.caption(m, name)
S.row([fc, cp2, pst]); print(S.save("organometallics.png")); S.check()
