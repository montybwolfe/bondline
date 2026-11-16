# Molfiles: to_molfile() writes a structure for ChemDraw (File > Open), keeping the drawn coordinates, wedges and
# abbreviations; from_molfile() reads one back with ChemDraw's layout kept exactly (a .mol saved from ChemDraw reads the
# same way). Here (R)-carvone goes out and comes back.
from bondline import *
S = Scene(Style())
m = from_smiles("CC1=CC[C@H](CC1=O)C(C)=C", align=(1, 7, 30))
to_molfile(m, "carvone.mol")
back = from_molfile("carvone.mol")
S.caption(m, "drawn"); S.caption(back, "read back from carvone.mol")
S.row([m, back]); print(S.save("molfile.png")); S.check(); stereo_report(back)
