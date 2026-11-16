# Teaching-figure pass: teaching marks as teaching figures draw them. Site dots on vertices (red Nu, blue E) drawn
# over the bond ends, their labels in the widest gap the bonds and any pointer leave; dashed pointers with one note for
# several targets, their tails beyond the molecule's ink in the note's strip; a curved callout ("C–C bond formed") at a
# coloured bond; pointers sent through their molecule are reported. Four products after one arrow: a second row
# under the first (row(..., under=)), one caption centred under each pair (caption([m1, m2], text)).
from bondline import *
RED, BLUE, GREEN = "#E8261A", "#1F4FD1", "#2E7D32"

print("== crossed aldol: two enolisable aldehydes, two sites for enolate formation")
S = Scene(Style())
def ald(n):
    m = Mol(); c = m.atom("C", (0, 0), name="c"); m.attach(c, 90, "O", 2); m.attach(c, -30, "H")
    a = m.attach(c, 210, name="a")
    if n == 3: m.attach(a, 150)
    return m
pa, ac = ald(3), ald(2)
for m in (pa, ac):
    S.dot(m["a"], RED, "Nu"); S.dot(m["c"], BLUE, "E")
S.pointer([pa["a"], ac["a"]], "<i>2 possible sites for\nenolate formation</i>", side="below", color=GREEN, dashed=True)
x1 = from_smiles("CC(O)C(C)C=O", align=(3, 5, 30))
x1.find(x1[1], x1[3]).color = RED
S.pointer((x1[1], x1[3]), "C–C bond formed", side=40, bend=-60, color=RED, length=1.4)
S.row([pa, "+", ac, S.arrow(above="NaOH"), x1])
for m in (pa, ac):
    for a in (m["a"], m["c"]): print(a.site[1], "label at", round(a._site[2]) % 360, "deg")
print(S.save("crossed.png", dpi=300)); S.check()

print("== one enolisable partner: notes from above and below, a bond aligned at 180 - deg, then flipped")
S = Scene(Style())
ap = from_smiles("CC(=O)c1ccccc1", align=(3, 1, 150), flip=True)
nb = from_smiles("[H]C(=O)c1ccc(cc1)[N+](=O)[O-]", align=(1, 3, -30), explicit_h=True)
print("acetophenone: O at", round(ang(ap[2].xy - ap[1].xy)) % 360, "deg, CH3 at", round(ang(ap[0].xy - ap[1].xy)) % 360, "deg")
S.dot(ap[0], RED, "Nu"); S.dot(nb[1], BLUE, "E")
S.pointer(ap[0], "<i>1 site for enolate\nformation</i>", side="below", color=GREEN, dashed=True)
S.pointer(ap[3], "<i>no α protons</i>", side="above", color=GREEN, dashed=True)
S.pointer(nb[3], "<i>no α protons</i>", side="above", color=GREEN, dashed=True)
ch = from_smiles("O=C(/C=C/c1ccc(cc1)[N+](=O)[O-])c1ccccc1", align=(13, 1, 150), flip=True)
ch.find(ch[2], ch[3]).color = RED
S.row([ap, "+", nb, S.arrow(above="NaOH", below="H2O/ EtOH"), ch])
S.caption(ch, "<i>99% yield\nCrossed aldol adduct</i>")
print(S.save("one_site.png", dpi=300)); S.check()

print("== four products after one arrow: a second row under the first, one caption under each pair")
S = Scene(Style())
pa, ac = ald(3), ald(2)
c1 = from_smiles("CC(O)C(C)C=O", align=(3, 5, 30)); c2 = from_smiles("CCC(O)CC=O", like=c1)
s1 = from_smiles("CCC(O)C(C)C=O", like=c1); s2 = from_smiles("CC(O)CC=O", like=c1)
for m, (i, j) in ((c1, (1, 3)), (c2, (2, 4)), (s1, (2, 4)), (s2, (1, 3))): m.find(m[i], m[j]).color = RED
S.row([pa, "+", ac, S.arrow(above="NaOH"), c1, "+", c2])
S.caption([c1, c2], "<i>Crossed aldol adducts</i>")
S.row([s1, "+", s2], under=[c1, "+", c2])
S.caption([s1, s2], "<i>Self condensation adducts</i>")
top, low = [S._span_box([Group(a), Group(b)]) for a, b in ((c1, c2), (s1, s2))]
print("rows centred:", round(abs((top[0] + top[2]) - (low[0] + low[2])) / 2, 1), "pt apart; row gap",
      round((top[1] - low[3]) / S.st.bond, 2), "BU")
print(S.save("two_rows.png", dpi=300)); S.check()

print("== negative: pointers sent through their molecule are reported")
S = Scene(Style())
ap = from_smiles("CC(=O)c1ccccc1", align=(3, 1, 30))
nb = from_smiles("[H]C(=O)c1ccc(cc1)[N+](=O)[O-]", align=(1, 3, -60), explicit_h=True)
S.dot(ap[0], RED, "Nu")
S.pointer(ap[0], "<i>1 site</i>", side="below", color=GREEN, dashed=True)
S.pointer(nb[3], "<i>no α protons</i>", side="above", color=GREEN, dashed=True)
S.row([ap, "+", nb]); S.check()
