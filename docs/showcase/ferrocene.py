# Ferrocene and acetylferrocene: freshly cracked cyclopentadiene is deprotonated by KOH to potassium
# cyclopentadienide, which iron(II) chloride turns into ferrocene; Friedel-Crafts acylation then gives acetylferrocene.
from bondline import *
S = Scene(Style())
cpd = from_smiles("C1C=CC=C1")
cpk = Mol(); r = cpk.ring(5); cpk.deloc(r, charge=-1)
fc = Mol(); fc.sandwich("Fe")
acfc = Mol(); M, (top, bot) = acfc.sandwich("Fe")
co = acfc.attach(top[4]); acfc.attach(co, label="O", order=2); acfc.attach(co)
S.caption(cpd, "cyclopentadiene"); S.caption(cpk, "K^+ C_5H_5^−")
S.caption(fc, "ferrocene"); S.caption(acfc, "acetylferrocene")
S.flow([cpd, S.arrow(above="KOH", below="DME"), cpk, S.arrow(above="FeCl_2·4H_2O", below="DMSO"), fc,
        S.arrow(above="Ac_2O, H_3PO_4", below="Δ"), acfc], width_cm=16)
print(S.save("ferrocene.png", dpi=260)); S.check()
