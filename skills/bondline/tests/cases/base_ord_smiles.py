# from_smiles layouts
from bondline import *
st = Style()
S = Scene(st)
asp = from_smiles("CC(=O)Oc1ccccc1C(=O)O"); caf = from_smiles("Cn1cnc2c1c(=O)n(C)c(=O)n2C")
men = from_smiles("CC(C)[C@@H]1CC[C@@H](C)C[C@H]1O"); stereo_report(men)
S.caption(asp, "aspirin"); S.caption(caf, "caffeine"); S.caption(men, "menthol")
S.row([asp, caf, men]); print(S.save("smiles.png")); S.check()
