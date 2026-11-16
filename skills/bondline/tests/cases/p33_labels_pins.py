# Natural-products pass: C6H5 is one label (never reversed to H5C6) and a contraction choice of its own; check()
# sees a wedge stroked through a label; with pins, a cage's ring bonds lose their wedges, wedges="exo" keeps the
# substituents' own, and a bridgehead's substituent is drawn plain (the cage fixes it: the methyl of camphor).
from bondline import *
print(label_units("C6H5"), label_units("C6F5"))
m = from_smiles("CC(=O)Nc1ccccc1", abbrev=["C6H5"], rotate=180)
print("labels:", [a.label for a in m.atoms if a.drawn])
w = Mol(); c = w.atom("C", (0, 0)); w.attach(c, 90, "CH3"); w.attach(c, 210); w.attach(c, -30)
w.attach(c, 72, style="wedge", length=1.25)
S = Scene(Style()); S.row([m, w]); print(S.save("labels.png", dpi=300)); S.check()
smi = "C[C@@]12CC[C@@H](C1(C)C)C[C@H]2O"
pins = {1: (0, 0), 9: (0.5, -0.5), 8: (1.5, -0.5), 4: (2.0, 0), 3: (1.32, 0.22), 2: (0.68, 0.22), 5: (1.0, 1.1)}
for wd in (None, "exo"):
    p = from_smiles(smi, xy=pins, wedges=wd)
    print(f"wedges={wd!r}:", sorted((int(b.a.name), int(b.b.name), b.style) for b in p.bonds if b.style != "plain"))
b = p.find(p[9], p[8]); b.style = "wedge"                # a front edge drawn as a wedge: depth, not a stereocentre
S = Scene(Style()); S.row([p]); print(S.save("camphor_pins.png", dpi=300)); S.check()
q = from_smiles("C[C@@]12CCCC[C@H]1CCCC2", xy={1: (0, 0.5), 2: (-0.866, 1.0), 3: (-1.732, 0.5), 4: (-1.732, -0.5),
                5: (-0.866, -1.0), 6: (0, -0.5), 7: (0.866, -1.0), 8: (1.732, -0.5), 9: (1.732, 0.5), 10: (0.866, 1.0)},
                wedges="exo")                            # a flat ring fusion: its methyl keeps the wedge
print("flat ring fusion keeps its wedge:", sorted((int(b.a.name), int(b.b.name), b.style) for b in q.bonds if b.style != "plain"))
