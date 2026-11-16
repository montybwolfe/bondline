# Gallery of famous natural products with full stereochemistry, drawn straight from PubChem isomeric SMILES.
# stereo_report() reads every stereocentre back from the drawing and is compared with RDKit's CIP labels.
from bondline import *
from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

NP = [  # (caption, PubChem isomeric SMILES, from_smiles options)
    ("morphine", "CN1CC[C@]23[C@@H]4[C@H]1CC5=C2C(=C(C=C5)O)O[C@H]3[C@H](C=C4)O", {}),
    ("strychnine", "C1CN2CC3=CCO[C@H]4CC(=O)N5[C@H]6[C@H]4[C@H]3C[C@H]2[C@@]61C7=CC=CC=C75", {}),
    ("reserpine", "CO[C@H]1[C@@H](C[C@@H]2CN3CCC4=C([C@H]3C[C@@H]2[C@@H]1C(=O)OC)NC5=C4C=CC(=C5)OC)OC(=O)C6=CC(=C(C(=C6)OC)OC)OC", {}),
    ("quinine", "COC1=CC2=C(C=CN=C2C=C1)[C@H]([C@@H]3C[C@@H]4CCN3C[C@@H]4C=C)O", {}),
    ("penicillin G", "CC1([C@@H](N2[C@H](S1)[C@@H](C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C", {}),
    ("cholesterol", "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C", dict(align=(19, 20, 90, 24))),  # C5->C10 up, ring A on the left: level, so ring-fusion H are drawn
]

def cip(smi):
    rm = Chem.MolFromSmiles(smi); rdCIPLabeler.AssignCIPLabels(rm)
    atoms = {a.GetIdx(): a.GetProp("_CIPCode") for a in rm.GetAtoms() if a.HasProp("_CIPCode")}
    bonds = sorted(b.GetProp("_CIPCode") for b in rm.GetBonds() if b.HasProp("_CIPCode"))
    return atoms, bonds

S = Scene(Style())
mols = []
for name, smi, kw in NP:
    m = from_smiles(smi, **kw)
    if name == "cholesterol":                      # C17: wedge the side chain (beta), not the ring bond C16-C17
        m.find(m[8], m[9]).style = "plain"; m.stereo(m[8], m[1], "wedge")
    rep = stereo_report(m, quiet=True)
    got = {int(a.name): v for a, v in rep["atoms"].items()}
    ref_atoms, ref_bonds = cip(smi)
    print(f"{name}: {len(got)} stereocentres {'OK' if got == ref_atoms else f'DIFF {got} vs {ref_atoms}'};"
          f" double bonds drawn {ref_bonds}, read back {sorted(rep['bonds'].values())}")
    S.caption(m, name); mols.append(m)
r1 = S.row(mols[:3])
S.row(mols[3:], under=r1)
print(S.save("structures.png", dpi=260)); S.check()
