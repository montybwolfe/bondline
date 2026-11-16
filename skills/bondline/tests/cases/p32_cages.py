# Bridged and caged ring systems (F37): RDKit's ring-system template where one exists (cubane, adamantane), its cage
# keeping its own perspective (bond lengths and angles not held to the flat rules, substituents at 1 BU); a perspective
# drawing pinned atom by atom (xy=), plain inside the cage, with the bond that passes behind broken where it crosses
# (behind=True, else from a 3D conformer fitted to the drawing); a crossing next to an atom is never broken.
from bondline import *
S = Scene(Style())
cub = from_smiles("COC(=O)C12C3C4C1C1C2C3C41C(=O)OC", abbrev=True)
ada = from_smiles("OC12CC3CC(CC(C3)C1)C2")
pins = {6: (0, 0), 1: (0.78, 0.2), 0: (1.65, 0.69), 2: (0.46, -0.58), 3: (-0.5, -0.3), 4: (-1.48, -0.59), 5: (-0.94, 0.29),
        7: (0, -0.96), 10: (-0.54, -1.85), 11: (0.43, -1.54), 8: (-0.94, -1.22), 9: (0.89, -1.41)}
bcn = from_smiles("O=C1C2CCCC1C(C)(O)CC2", xy=pins)
bcn.find(bcn[6], bcn[7]).behind = True; bcn.find(bcn[7], bcn[9]).behind = True
S.row([cub, ada, bcn], gap=1.5); print(S.save("cages.png", dpi=300)); S.check()
for name, m in (("cubane", cub), ("adamantanol", ada), ("pinned bicyclononanone", bcn)):
    ext = [round(norm(b.b.xy - b.a.xy), 2) for b in m.bonds if (getattr(b.a, "_tmpl", None) == "cage") != (getattr(b.b, "_tmpl", None) == "cage")]
    print(f"{name}: cage atoms {sum(getattr(a, '_tmpl', None) == 'cage' for a in m.atoms)}, broken bonds "
          f"{sum(bool(getattr(b, '_gaps', None)) for b in m.bonds)}, wedges {sum(b.style in ('wedge', 'hash') for b in m.bonds)}, "
          f"bonds out of the cage {sorted(set(ext))}")
print("pins kept:", all(norm(bcn[i].xy - bcn[6].xy - np.array(v)) < 1e-6 for i, v in pins.items()))
try: from_smiles("CCO", xy={5: (0, 0)})
except ValueError as e: print("raised:", e)
