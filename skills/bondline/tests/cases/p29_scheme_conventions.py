# Teaching-figure pass: conventions a scheme needs. A reagent multiple (2×Me3SiCl); from_smiles
# align=(i, j, deg, k) putting atom k on the bond's left (a level chain, C=O up); a dashed crossed arrow pointing left (a
# ring that does not form); like= re-laying a pinned substituent when its ring carbon becomes a gem pair; ring-fusion
# stereocentres drawn with their H straight up or down (the steroid convention) when the drawing is level; letters in
# rings; a coloured caption.
from bondline import *
RED = "#E8261A"

print("== level chain, C=O up; a reagent multiple")
S = Scene(Style())
zing = from_smiles("COc1cc(CCC(C)=O)ccc1O", align=(5, 8, 0, 9), abbrev=True)
print("C=O points", round(ang(zing[9].xy - zing[7].xy)) % 360, "deg; chain end-to-end", round(ang(zing[8].xy - zing[5].xy)) % 360, "deg")
S.row([zing, S.arrow(above="2xLDA", below="2×Me3SiCl")]); print(S.save("multiples.png", dpi=300)); S.check()
try: from_smiles("C=CC(=O)OC", align=(0, 1, 30, 3), abbrev=True)
except ValueError as e: print("raised:", e)

print("== the ring that does not form (dashed crossed arrow, pointing left) and the one that does")
S = Scene(Style())
dk = from_smiles("CC(=O)CCCCC(C)=O", align=(0, 8, 0, 2))
six = from_smiles("CC(=O)C1=C(C)CCCC1", align=(1, 3, 270)); six.find(six[3], six[4]).color = RED
eight = from_smiles("O=C1CCCCCC(C)=C1", align=(1, 0, 90)); eight.find(eight[7], eight[9]).color = RED
S.caption(six, "<i>85% yield</i>"); S.caption(eight, "<i>NOT FORMED</i>")
S.row([eight, S.arrow(kind="<-/-", above="H+", dashed=True), dk, S.arrow(above="H2SO4", below="H+"), six])
print(S.save("not_formed.png", dpi=300)); S.check()

print("== Michael addition at a substituted ring carbon: the product's gem pair is symmetric (like=)")
S = Scene(Style())
mc = from_smiles("O=C1CCCC1C", align=(1, 0, 90)); pr = from_smiles("O=C1CCCC1(C)CCC(=O)OC", like=mc, abbrev=["OMe"])
c5 = pr[5]; ts = sorted(round(ang(n.xy - c5.xy)) % 360 for n in pr.nbrs(c5) if n not in (pr[1], pr[4]))
print("exocyclic bonds at the quaternary carbon:", ts)
S.row([mc, S.arrow(above="t-BuOK"), pr]); print(S.save("gem_like.png", dpi=300)); S.check()

print("== steroids: ring-fusion H up or down when the drawing is level; letters in the rings; coloured caption")
S = Scene(Style())
testo = from_smiles("C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@]34C)[C@@H]1CC[C@@H]2O", labels={0: "Me", 15: "Me"})
hs = [a for a in testo.atoms if a.label == "H"]
print("fusion H atoms:", len(hs), "directions", sorted(round(ang(h.xy - testo.nbrs(h)[0].xy)) % 360 for h in hs),
      "styles", sorted(testo.find(h, testo.nbrs(h)[0]).style for h in hs))
print(stereo_report(testo)["atoms"].values())
sk = from_smiles("C[C@]12CCC3C(CCC4CCCC[C@]34C)C1CC[C@@H]2[*]", labels={0: "Me", 14: "Me", 19: "R"})
for letter, ring in (("A", (12, 11, 10, 9, 8, 13)), ("B", (8, 7, 6, 5, 4, 13)), ("C", (5, 4, 3, 2, 1, 15)), ("D", (1, 15, 16, 17, 18))):
    S.ring_text([sk[i] for i in ring], letter)
S.caption(testo, "<b>Testosterone</b>"); S.caption(sk, "Steroid skeleton", color=RED)
S.row([sk, testo], gap=1.5); print(S.save("steroids.png", dpi=300)); S.check()
