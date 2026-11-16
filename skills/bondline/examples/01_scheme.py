# Reaction scheme: structures from SMILES, conditions above and below the arrow (numbered lines, symbols rather than
# words), compound names as captions, and a layout that wraps at a reaction arrow when the row gets too wide.
from bondline import *
S = Scene(Style())
# align=(i, j, deg): the bond from atom i to atom j (SMILES order) points along deg; flip=True mirrors left-right.
# Here every C-O points up with the ring on the left, so the three structures read alike.
ald = from_smiles("O=Cc1ccccc1", align=(1, 0, 90), flip=True)
alc = from_smiles("CC(O)c1ccccc1", align=(1, 2, 90), flip=True)
ket = from_smiles("CC(=O)c1ccccc1", align=(1, 2, 90), flip=True)
for m, name in ((ald, "benzaldehyde"), (alc, "1-phenylethanol"), (ket, "acetophenone")): S.caption(m, name)
S.flow([ald, S.arrow(above=["1. MeMgBr, Et_2O, 0 °C", "2. H_3O^+"]), alc,
        S.arrow(above="PCC", below="CH_2Cl_2, rt"), ket], width_cm=16)
print(S.save("scheme.png")); S.check()
