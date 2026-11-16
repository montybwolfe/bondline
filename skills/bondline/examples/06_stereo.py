# Stereochemistry: wedges and hashes come from the SMILES (or attach(..., style="wedge"/"hash")); stereo_report() reads
# the CIP labels back from the drawing, so the captions say what is drawn. An enantiomer is the other SMILES (@ for @@);
# flip=True only mirrors the layout, the configuration stays. Descriptors are italic: (<i>R</i>).
from bondline import *
S = Scene(Style())
r = from_smiles("C[C@@H](O)CC", align=(1, 0, 210))                   # butan-2-ol
s = from_smiles("C[C@H](O)CC", align=(1, 0, 210), flip=True)         # its enantiomer, drawn as the mirror image
e = from_smiles("C/C=C/C(=O)O")                                       # but-2-enoic acid
for m in (r, s, e):
    rep = stereo_report(m, quiet=True)
    S.caption(m, " ".join(f"(<i>{x}</i>)" for x in list(rep["atoms"].values()) + list(rep["bonds"].values())))
S.row([r, s, e]); print(S.save("stereo.png")); S.check()
