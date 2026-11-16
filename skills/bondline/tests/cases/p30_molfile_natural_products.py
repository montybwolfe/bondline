# Teaching-figure pass (natural products): a molfile round trip keeps every stereocentre (RDKit could lay a
# bridged system such as morphine out as its mirror image, which written with the drawn wedges was the enantiomer) and
# writes drawn-out abbreviations as contracted SUP groups, so the copy read back has the same labels (silvestrol).
from bondline import *
NP = {"morphine": "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5",
      "silvestrol": "CO[C@H]1[C@@H](O[C@H](CO1)[C@@H](CO)O)OC2=CC3=C(C(=C2)OC)[C@@]4([C@@H]([C@@H]([C@H]([C@@]4(O3)"
                    "C5=CC=C(C=C5)OC)C6=CC=CC=C6)C(=O)OC)O)O"}
for name, smi in NP.items():
    m = from_smiles(smi, abbrev=True); t = to_molfile(m); m2 = from_molfile(t)
    a, b = stereo_report(m, quiet=True)["atoms"], stereo_report(m2, quiet=True)["atoms"]
    same = {x.name: v for x, v in a.items()} == {x.name: v for x, v in b.items()}
    labs = sorted(x.label for x in m.atoms if x.drawn) == sorted(x.label for x in m2.atoms if x.drawn)
    sup = sum(1 for l in t.split("\n") if l.startswith("M  SMT"))
    print(f"{name}: atoms {len(m.atoms)} -> {len(m2.atoms)}, R/S identical {same}, labels identical {labs}, SUP groups {sup}")
    print("  ", "".join(v for _, v in sorted(((x.name, v) for x, v in a.items()))))
