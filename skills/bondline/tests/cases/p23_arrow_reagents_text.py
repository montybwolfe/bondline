# Teaching-figure pass: reagent structures carried on reaction arrows; left, retrosynthetic-left and
# crossed arrows; chemical typography in running text (formulas, S_N2, pKa, italic prefixes, circled charges).
from bondline import *
st = Style()

print("== LDA drawn over the arrow, conditions under it; an aldehyde electrophile over the next arrow")
S = Scene(st)
k = Mol(); c = k.atom("C", (0, 0)); k.attach(c, 90, "O", 2); k.attach(c, 210, "Ph"); ca = k.attach(c, -30); k.attach(ca, 30)
lda = Mol(); n = lda.atom("N", (0, 0)); lda.attach(n, 270, "Li")
for t in (150, 30):
    ch = lda.attach(n, t); lda.attach(ch, 90); lda.attach(ch, t - 60 if t == 30 else t + 60)
e = Mol(); c = e.atom("C", (0, 0)); e.attach(c, 90, "OLi"); e.attach(c, 210, "Ph"); ca = e.attach(c, -30, order=2); e.attach(ca, 30)
al = Mol(); c = al.atom("C", (0, 0)); al.attach(c, 90, "O", 2); al.attach(c, 210, "H"); x = al.attach(c, -30); al.attach(x, 30)
p = from_smiles("CCC(O)C(C)C(=O)c1ccccc1", align=(6, 7, 90))
S.row([k, S.arrow(above=[lda, "<i>LDA</i>"], below=["THF", "−78 °C"]), e,
       S.arrow(above=al, below=["add electrophile", "then H3O+"]), p])
print(S.save("reagents.png", dpi=600)); S.check()
S.balance([k], [p])

print("== arrow kinds: <- (with B:⊖ in text), -/-> (does not happen), <= (retrosynthesis, pointing left)")
S = Scene(st)
a = from_smiles("CC(=O)C"); b = from_smiles("CC(=O)C"); c = from_smiles("CC(O)C")
S.row([a, S.arrow(above="B:⊖", kind="<-"), b, S.arrow(above="NaBH4", kind="-/->"), c, S.arrow(kind="<="), from_smiles("CC=O")])
print(S.save("kinds.png", dpi=600)); S.check()
try: S.arrow(kind="-->"); print("no error (should not happen)")
except ValueError as err: print("raised:", err)

print("== chemical typography in conditions and captions")
S = Scene(st)
m = from_smiles("c1ccccc1"); S.caption(m, "(<i>S</i>)-proline, (R)-aldol, pKa ≈ 20")
S.row([m, S.arrow(above=["1. Pd(PPh3)4, CuI, Et3N", "2. tBuOK, CH2Cl2, −H2O"], below=["[Ru(bpy)3]2+, NH4+, PF6-", "SN2, SNAr, E1cB, E2"]),
       from_smiles("c1ccccc1"), S.arrow(above=["O-silylation, p-TsOH, m-CPBA", "Fe2+, Pd0, R1R2NH, kcal mol-1"], below="H2, Pd/C; \\C2")])
print(S.save("typography.png", dpi=600)); S.check()

print("== vertical arrows of every kind (a scheme that turns a corner: equilibria, retrosynthesis), text clear of the pair")
S = Scene(Style())
for i, k in enumerate(("->", "<=>", "<->", "=>", "-/->", "<-")):
    S.arrow(kind=k, vertical=True, at=(i * 2.2, 0), above="H+")
print(S.save("vertical_kinds.png", dpi=300)); S.check()
