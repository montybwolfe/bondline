# Curated skeleton layouts (TEMPLATES) and unambiguous wedges.
# 1. a penam (penicillin G) is laid out as textbooks draw it: the beta-lactam square upright (edges level and
#    vertical), the thiazolidine to its right with S on top; align= still turns it.
# 2. strychnine takes the strychnos template: every bond 1 BU, nothing crossed, stereo as the SMILES has it (RDKit's own
#    template stretched two bonds to 2.07 BU).
# 3. a stereocentre whose wedge left its three other bonds within 180 deg (a T) has the wedge moved to a bond that
#    leaves them spread: artemisinin's ketal carbon, strychnine's spiro carbon.
from bondline import *
import numpy as np
from rdkit import Chem
def cip(smi):
    return sorted(Chem.FindMolChiralCenters(Chem.MolFromSmiles(smi), useLegacyImplementation=False))
def got(m):
    return sorted((a.name, v) for a, v in stereo_report(m, quiet=True)["atoms"].items())
def edges(m, ids):
    sq = [m[i] for i in ids]
    return [round(float(np.degrees(np.arctan2(*(b.xy - a.xy)[::-1]))) % 180) for a, b in zip(sq, sq[1:] + sq[:1])]
def lengths(m):
    L = [norm(b.b.xy - b.a.xy) for b in m.bonds]
    return f"{min(L):.2f}-{max(L):.2f} BU"

pen = "CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O"
a = from_smiles(pen, abbrev=["CO2H", "Ph"]); b = from_smiles(pen, abbrev=["CO2H", "Ph"], align=(4, 19, 0))
print("penam square edges:", edges(a, (4, 5, 16, 18)), "S above N:", bool(a[3].xy[1] > a[18].xy[1]),
      "aligned:", edges(b, (4, 5, 16, 18)), "stereo", got(a) == cip(pen))
S = Scene(Style()); S.row([a, b], gap=2.0); print(S.save("penam.png", dpi=300)); S.check()

st = "C1CN2CC3=CCO[C@H]4CC(=O)N5[C@H]6[C@H]4[C@H]3C[C@H]2[C@@]61C7=CC=CC=C75"
art = "C[C@@H]1CC[C@H]2[C@H](C(=O)O[C@H]3[C@@]24[C@H]1CC[C@](O3)(OO4)C)C"
s, t = from_smiles(st), from_smiles(art)
print("strychnine bonds", lengths(s), "stereo", got(s) == cip(st), got(s))
print("artemisinin stereo", got(t) == cip(art), got(t))
S = Scene(Style()); S.row([s, t], gap=2.0); print(S.save("strychnos_cage.png", dpi=300)); S.check()

# 4. morphinans: A/B/C flat on the hexagon grid, the 4,5-ether and the N-bridge in perspective (free bonds), the bridge
#    crossing ring C with the back bond broken (natural morphine's bridge in front, dextromethorphan's behind); wedges on
#    exocyclic bonds or bonds to non-stereocentres (C13 on its bridge bond, C14 on C14-C8)
mo = "CN1CC[C@]23[C@@H]4[C@H]1CC5=C2C(=C(C=C5)O)O[C@H]3[C@H](C=C4)O"
dx = "CN1CC[C@@]23CCCC[C@@H]2[C@@H]1CC4=C3C=C(C=C4)OC"
p, q = from_smiles(mo), from_smiles(dx, abbrev=["OMe"])
S = Scene(Style()); S.row([p, q], gap=2.0); print(S.save("morphinans.png", dpi=300)); S.check()
for nm, m, smi in (("morphine", p, mo), ("dextromethorphan", q, dx)):
    print(nm, "stereo", got(m) == cip(smi), "marks", sorted((b.a.name, b.b.name, b.style) for b in m.bonds if b.style != "plain"),
          "broken", sorted((b.a.name, b.b.name) for b in m.bonds if getattr(b, "_gaps", None)))
