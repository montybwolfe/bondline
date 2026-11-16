# Takasago synthesis of (-)-menthol: Li-catalysed hydroamination of myrcene, Rh-BINAP asymmetric 1,3-hydrogen shift
# of the allylic amine to the (R)-citronellal (E)-enamine, hydrolysis, ZnBr2 carbonyl-ene cyclisation, hydrogenation.
from bondline import *
from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

SMI = {  # isomeric SMILES (PubChem; the enamine is citronellal's with CH2CHO -> (E)-CH=CH-NEt2)
    "myrcene": "C=CC(=C)CCC=C(C)C",
    "amine": "CC(=CCC/C(=C/CN(CC)CC)/C)C",               # geraniol with OH -> NEt2
    "enamine": "C[C@H](CCC=C(C)C)/C=C/N(CC)CC",
    "citronellal": "C[C@H](CCC=C(C)C)CC=O",
    "isopulegol": "C[C@@H]1CC[C@H]([C@@H](C1)O)C(=C)C",
    "menthol": "C[C@@H]1CC[C@H]([C@@H](C1)O)C(C)C",
}
S = Scene(Style())
myr = from_smiles(SMI["myrcene"])
ami = from_smiles(SMI["amine"], like=myr)
ena = from_smiles(SMI["enamine"], like=ami)
cit = from_smiles(SMI["citronellal"], like=ena)
iso = from_smiles(SMI["isopulegol"])
men = from_smiles(SMI["menthol"], like=iso)
mols = dict(myrcene=myr, amine=ami, enamine=ena, citronellal=cit, isopulegol=iso, menthol=men)
for k, m in mols.items():
    rm = Chem.MolFromSmiles(SMI[k]); rdCIPLabeler.AssignCIPLabels(rm)
    ref = {a.GetIdx(): a.GetProp("_CIPCode") for a in rm.GetAtoms() if a.HasProp("_CIPCode")}
    refb = sorted(b.GetProp("_CIPCode") for b in rm.GetBonds() if b.HasProp("_CIPCode"))
    rep = stereo_report(m, quiet=True); got = {int(a.name): v for a, v in rep["atoms"].items()}
    print(k, "stereo", "OK" if (got, sorted(rep["bonds"].values())) == (ref, refb) else f"DIFF {got} {ref} {rep['bonds']} {refb}")
S.caption(myr, "myrcene"); S.caption(ami, "<i>N</i>,<i>N</i>-diethylgeranylamine")
S.caption(ena, "(<i>R</i>)-citronellal (<i>E</i>)-enamine"); S.caption(cit, "(<i>R</i>)-citronellal")
S.caption(iso, "(−)-isopulegol"); S.caption(men, "(−)-menthol")
S.flow([myr, S.arrow(above="Et_2NH", below="<i>n</i>BuLi (cat.)"), ami,
        S.arrow(above=["[Rh((<i>S</i>)-BINAP)(cod)]ClO_4", "(cat.)"], below="THF, Δ"), ena,
        S.arrow(above="H_3O^+"), cit,
        S.arrow(above="ZnBr_2"), iso,
        S.arrow(above="H_2", below="Raney Ni"), men], width_cm=16)
print(S.save("synthesis.png", dpi=260)); S.check()
