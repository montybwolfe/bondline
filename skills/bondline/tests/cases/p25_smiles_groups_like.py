# Teaching-figure pass: structures from SMILES write silyl ethers, nitro groups, O-Li and the like as
# labels (Me3SiO, O2N, OLi; abbrev=True adds OMe, CO2Et ...), and like= lays each product out as its starting material is
# drawn, so a sequence reads in one orientation (gingerol from its protected ketone; retro-aldol to zingerone).
from bondline import *
st = Style()

print("== gingerol by a kinetic enolate aldol, then retro-aldol to zingerone")
S = Scene(st)
sm = from_smiles("COc1cc(CCC(C)=O)ccc1O[Si](C)(C)C", abbrev=True)
en = from_smiles("COc1cc(CCC(=C)O[Li])ccc1O[Si](C)(C)C", like=sm, abbrev=True)
hx = from_smiles("CCCCCC=O")
gi = from_smiles("CCCCCC(O)CC(=O)CCc1ccc(O)c(OC)c1", like=sm, abbrev=True)
zi = from_smiles("CC(=O)CCc1ccc(O)c(OC)c1", like=sm, abbrev=True)
for m in (sm, en, gi, zi): print(sorted(a.label for a in m.atoms if a.label not in ("C",)))
S.caption(en, "<i>kinetic enolate</i>"); S.caption(gi, "gingerol, 57%"); S.caption(zi, "zingerone")
S.flow([sm, S.arrow(above="<i>LiHMDS</i>", below=["THF", "−78 °C"]), en, S.arrow(above=hx, below=["add electrophile", "then H3O+"]),
        gi, S.arrow(above="Δ", below="<i>retro</i>-aldol"), zi], width_cm=17)
print(S.save("gingerol.png", dpi=300)); S.check()

print("== automatic labels: nitro, silyl, stannyl, boronic acid, sulfonates, CF3; none with abbrev=False")
S = Scene(st); ms = []
for smi in ("O=[N+]([O-])c1ccc(C=O)cc1", "C=CC[Sn](CCCC)(CCCC)CCCC", "OB(O)c1ccc(C=O)cc1", "FC(F)(F)S(=O)(=O)Oc1ccccc1",
            "CC(C)[Si](C(C)C)(C(C)C)OCC=C", "Cc1ccc(cc1)S(=O)(=O)N1CC=CC1"):
    m = from_smiles(smi); ms.append(m); print(smi, "->", sorted(a.label for a in m.atoms if a.label not in ("C",)))
print("abbrev=False:", sorted(a.label for a in from_smiles("O=[N+]([O-])c1ccccc1", abbrev=False).atoms if a.label != "C"))
S.flow(ms, width_cm=17); print(S.save("groups.png", dpi=300)); S.check()
