# Stage 7: from_molfile / to_molfile. Round trip keeps layout, wedges and CIP; aliases and contracted superatoms (SUP)
# come back as labels; export draws abbreviations out as real atoms around fixed drawn atoms.
from bondline import *
from rdkit import Chem
from rdkit.Chem import rdDepictor, rdMolDescriptors

print("== round trip")
for smi in ["C[C@H](N)C(=O)O", "O=C(OC)c1ccc([C@@H](O)CC)cc1", "C[C@@H]1CC[C@H](O)CC1"]:
    m = from_smiles(smi); back = from_molfile(to_molfile(m))
    A = np.array([a.xy for a in m.atoms]); Bk = np.array([a.xy for a in back.atoms])
    st0 = sorted((b.a.name, b.b.name, b.style) for b in m.bonds if b.style != "plain")
    st1 = sorted((b.a.name, b.b.name, b.style) for b in back.bonds if b.style != "plain")
    cip = lambda mm: sorted((a.name, v) for a, v in stereo_report(mm, quiet=True)["atoms"].items())
    print(f"  {smi}: atoms {len(m.atoms)}->{len(back.atoms)}, layout {float(np.abs((A - A.mean(0)) - (Bk - Bk.mean(0))).max()) < 0.005},"
          f" wedges kept {st0 == st1}, CIP kept {cip(m) == cip(back)} {cip(back)}")

print("== abbreviations out and back")
q = Mol(); c = q.atom("C", (0, 0)); q.attach(c, 30, "OMe"); q.attach(c, 150, "Ph"); q.attach(c, 270, "CO2Et")
r = Chem.MolFromMolBlock(to_molfile(q, "abbr.mol")); print("  drawn out:", rdMolDescriptors.CalcMolFormula(r), "= formula():", hill(formula(q)[0]))
print("  as aliases and back:", [a.label for a in from_molfile(to_molfile(q, expand=False)).atoms])

print("== contracted superatoms (SUP) from a drawing program")
mol = Chem.MolFromSmiles("COc1ccc(C(=O)OC(C)(C)C)cc1"); rdDepictor.Compute2DCoords(mol)
for atoms, label in [((0, 1), "OMe"), ((9, 10, 11, 12), "tBu")]:
    sg = Chem.CreateMolSubstanceGroup(mol, "SUP"); sg.SetProp("LABEL", label)
    for a in atoms: sg.AddAtomWithIdx(a)
for v3 in (False, True):
    m = from_molfile(Chem.MolToMolBlock(mol, forceV3000=v3))
    print("  V3000" if v3 else "  V2000", [a.label for a in m.atoms if a.drawn], len(m.atoms), "atoms; drawn out:",
          len(from_molfile(Chem.MolToMolBlock(mol, forceV3000=v3), contract=False).atoms))
S = Scene(); S.add(m); S.check(); S.save("sup.png")
