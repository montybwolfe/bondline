# Scheme layout: brackets with a count and a caption under them, spaced in a row as part of their structure; captions
# side by side kept apart; several captions on one structure stacked; a caption under a group below the members' own;
# a reagent structure's caption carried with it above the arrow.
from bondline import *

def white(b1, b2):                                         # horizontal white between two boxes (pt)
    return max(b2[0] - b1[2], b1[0] - b2[2])

S = Scene(Style())
a = from_smiles("CCC#C", align=(1, 3, 0)); cu = from_smiles("CCC#C[Cu]", align=(1, 3, 0)); d = from_smiles("CCC#CC#CCC", align=(1, 3, 0))
S.caption(a, "alkyne"); S.caption(cu, "copper acetylide"); S.caption(d, "diyne"); S.brackets([cu], count=2)
ar1, ar2 = S.arrow(above="CuCl"), S.arrow(above="O2")
S.row([a, ar1, cu, ar2, d]); print(S.save("bracket_count.png", dpi=300)); S.check()
prims, texts = S._build()
bg = S._bracket_geo([e for e in S.extras if e[0] == "brackets"][0])
cap = next(t for t in texts if "".join(r[0] for r in t["runs"]) == "copper acetylide")
cnt = next(t for t in bg["texts"])
print("caption under the brackets:", cap["box"][3] < bg["box"][1] + 0.5,
      "| count is a subscript at the lower right:", round(cnt["box"][1] - bg["box"][1], 1) <= 0.0,
      "| arrows clear the brackets (pt):", round(float(white(S._item_geo(ar1)["box"], bg["box"])), 1), round(float(white(bg["box"], S._item_geo(ar2)["box"])), 1))

S = Scene(Style())
ms = [from_smiles(s) for s in ("c1ccccc1", "c1ccncc1", "c1ccc2ncccc2c1")]
for m, t in zip(ms, ("all positions equal", "three kinds of position", "seven kinds of position")): S.caption(m, t)
S.row(ms, gap=2.5); print(S.save("caption_gap.png", dpi=300)); S.check()
prims, texts = S._build(); cb = sorted([t["box"] for t in texts if not t.get("atom")], key=lambda b: b[0])
print("white between neighbouring captions (BU):", [round(float(white(cb[i], cb[i + 1])) / 14.4, 2) for i in range(len(cb) - 1)])

S = Scene(Style())
p, q = from_smiles("CC(=O)c1ccc(O)cc1"), from_smiles("N#Cc1ccc(O)cc1")
S.caption(p, "<b>4</b>"); S.caption(p, "65%"); S.caption(q, "<b>5</b>"); S.caption(q, "62%"); S.caption([p, q], "from water")
r = S.row([p, q], gap=1.5)
x, y = from_smiles("CCO"), from_smiles("CCN")
S.caption(x, "one line"); S.caption([x, y], "a caption under both")
S.row([x, y], under=r); print(S.save("caption_stack.png", dpi=300)); S.check()
prims, texts = S._build(); tb = {"".join(r_[0] for r_ in t["runs"]): t["box"] for t in texts if not t.get("atom")}
print("stacked in order:", tb["4"][1] > tb["65%"][3] - 0.5 and tb["65%"][1] > tb["from water"][3],
      "| first captions on one line:", abs(tb["4"][3] - tb["5"][3]) < 0.01, "| joint below the own caption:", tb["one line"][1] > tb["a caption under both"][3])

S = Scene(Style())
s, t = from_smiles("Brc1ccc(C)cc1"), from_smiles("Cc1ccc(cc1)C1CCCO1")
bor = from_smiles("CC1(C)OB(C2CCCO2)OC1(C)C", abbrev=False); S.caption(bor, "<b>3</b>")
ar = S.arrow(above=bor, below="[Pd]"); S.row([s, ar, t]); print(S.save("reagent_caption.png", dpi=300)); S.check()
prims, texts = S._build(); c3 = next(t_ for t_ in texts if "".join(r_[0] for r_ in t_["runs"]) == "3")
print("reagent caption sits above the shaft (pt):", round(float(c3["box"][1] - ar.y), 1), "| structure above its caption:",
      S._mol_box(bor, curlies=False)[1] > c3["box"][3])

# A corner arrow (via=) carrying a reagent structure grows to fit it, the stack centred beside the shaft; split=True
# sets above= to the left and below= to the right.
for split in (False, True):
    S = Scene(Style())
    u, v = from_smiles("OCc1ccccc1"), from_smiles("O=Cc1ccccc1")
    w, z = from_smiles("O=C(O)c1ccccc1"), from_smiles("COC(=O)c1ccccc1")
    rg = from_smiles("BrCC(=O)c1ccccc1"); S.caption(rg, "<b>6</b>")
    r1 = S.row([u, S.arrow(above="MnO2"), v])
    vi = S.arrow(vertical=True, above=rg, below=["K2CO3", "DMF"], split=split)
    S.row([z, S.arrow(kind="<-", above="MeI"), w], under=r1, via=vi)
    print(S.save(f"via_reagent_{'split' if split else 'right'}.png", dpi=300)); S.check()
    g = S._item_geo(vi); top, tip = vi.y, vi.y - vi.length * 14.4
    parts = [t["box"] for t in g["texts"]] + [S._reagent_ext(rg)[1]]
    print("stack within the arrow:", all(b[3] < top and b[1] > tip for b in parts),
          "| sides (structure, text):", "left" if S._reagent_ext(rg)[1][2] < vi.x else "right",
          "right" if all(t["box"][0] > vi.x for t in g["texts"] if "".join(r[0] for r in t["runs"]) != "6") else "left")

# Pointers: the tip stops at the letters it points at (a callout to a wide label at an angle reaches it), short of
# the bonds round a crowded vertex, and on the side of a bond lying along the pointer; notes of later pointers keep
# clear of earlier ones; a note above an arrow's conditions may reach over the arrow; a curved callout has no kink.
S = Scene(Style())
ac = from_smiles("OC(=O)C1CCCCC1", abbrev=["CO2H"], align=(3, 1, 0))
sp = from_smiles("C1CCC2(CC1)CCCC2")
S.pointer(ac[1], "acid", side=20, bend=-60); S.pointer(sp[3], "spiro centre", side=90)
S.row([ac, sp], gap=3.0); print(S.save("pointer_tips.png", dpi=300)); S.check()
lab = [h for g in ac[1]._g for h in g["hulls"]]; tips = [np.asarray(pr[1])[0] for p in S.pointers for pr in S._ptr_geo(p)["prims"] if pr[0] == "poly"]
print("callout tip to its label (pt):", round(min(float(np.hypot(*(h - tips[0]).T).min()) for h in lab), 1),
      "| spiro tip to its vertex (BU):", round(float(norm(tips[1] - S.P(sp[3]))) / 14.4, 2))
S = Scene(Style())
d = from_smiles("CC(C)C[C@]([H])(O)CC(C)=O", explicit_h=True, align=(0, 8, 0))
S.row([d]); S.pointer(d[6], "hydroxy up", side="below"); S.pointer(d[5], "hydrogen down", side="below")
print(S.save("pointer_notes.png", dpi=300)); S.check()
nb_ = [t["box"] for p in S.pointers for t in S._ptr_geo(p)["texts"]]
print("notes apart (BU):", round(float(max(nb_[1][0] - nb_[0][2], nb_[0][0] - nb_[1][2])) / 14.4, 2))
S = Scene(Style())
k0 = from_smiles("O=C1CCCCC1C", align=(1, 0, 90)); k1 = from_smiles("O=C1CCCCC1(C)CC=C", like=k0)
S.pointer((k1[6], k1[8]), "bond formed", side=115, bend=-60, color="#E8261A")
ar = S.arrow(above="LDA", below="allyl bromide"); S.row([k0, ar, k1]); print(S.save("pointer_overhang.png", dpi=300)); S.check()
print("note reaches over the arrow:", S._ptr_geo(S.pointers[0])["texts"][0]["box"][0] < S._item_geo(ar)["box"][2],
      "| structure gap after the arrow (BU):", round(float(S._mol_box(k1, curlies=False)[0] - S._item_geo(ar)["box"][2]) / 14.4, 2))
S = Scene(Style())
im = from_smiles("[*]/C=[N+]1CCCC1", labels={0: "R"}, align=(2, 1, 270))
S.pointer(im[2], "no empty p orbital", side=-15, bend=60); S.row([im]); S.save("pointer_short_callout.png", dpi=300); S.check()
g = S._ptr_geo(S.pointers[0]); sh = np.asarray([pr for pr in g["prims"] if pr[0] == "line"][0][1]); hd = np.asarray([pr for pr in g["prims"] if pr[0] == "poly"][0][1])
print("shaft meets head straight (deg):", round(abs(adiff(ang(sh[-1] - sh[-2]), ang(hd[0] - hd[2]))), 1))
S = Scene(Style())
m1, m2 = from_smiles("CCCC", align=(0, 3, 0)), from_smiles("CCCC", align=(0, 3, 0))
S.row([m1, S.arrow(), m2], gap=2.0); S.pointer(m1[3], "late note", side="right"); S.save("pointer_late.png", dpi=300)
print("a note declared after the row, over the arrow:", S.check(quiet=True))

# Cycles: a structure at the end of a branch (drawn, not dropped), step names beside their arcs, a feed from outside
# into a station with its note.
S = Scene(Style())
st0 = Mol(); st0.atom("L_{n}Pd^{0}", (0, 0)); st1 = Mol(); st1.atom("<b>B</b>", (0, 0)); st2 = Mol(); st2.atom("<b>D</b>", (0, 0))
c = S.cycle([(st0, 90), (st1, -30), (st2, 210)], radius=4)
prod = from_smiles("CC(=O)c1ccc(C)cc1"); S.branch(c, 2, prod, kind="out")
S.branch(c, 0, "<i>oxidative addition</i>", kind="step"); S.branch(c, 0, "Pd(OAc)2", kind="feed", note="reduction")
print(S.save("cycle_features.png", dpi=300)); S.check()
pb = S._mol_box(prod, curlies=False); cb = (c["c"] * 14.4, c["r"] * 14.4)
print("product drawn outside the cycle:", prod in S.mols and min(float(norm(np.array([x, y]) - cb[0])) for x in pb[::2] for y in pb[1::2]) > cb[1],
      "| texts:", sorted("".join(r[0] for r in t["runs"]) for t in S._build()[1] if not t.get("atom")))

# Negative: an equilibrium arrow placed by hand through a label and a text is reported; text clears a crossed arrow.
S = Scene(Style())
h = hydronium(explicit=False); S.row([h])
S.place(S.arrow(vertical=True, kind="<=>", length=1.6), (h.key.xy[0] - 0.3, h.key.xy[1] + 0.6)); S.text("HCl", (h.key.xy[0] - 0.3, h.key.xy[1] + 1.4))
S.save("placed_arrow.png", dpi=300); print("placed arrow:", S.check(quiet=True))
S = Scene(Style()); S.row([from_smiles("Brc1ccccc1"), S.arrow(kind="-/->", above="SN2"), S.text("no reaction")])
S.save("crossed_text.png", dpi=300); print("crossed arrow text:", S.check(quiet=True))

# Branches on reaction arrows: a by-product (a structure) curving out of one arrow; a co-reagent curving into and a
# by-product curving out of another, below its conditions' free side, their labels kept apart.
S = Scene(Style())
a, b, c = from_smiles("CCO"), from_smiles("CC=O"), from_smiles("CC(=O)O")
ar1 = S.arrow(below="<i>oxidation</i>"); S.branch(ar1, from_smiles("O=C=O"), kind="out")
ar2 = S.arrow(above="[Cu]"); S.branch(ar2, "O2", kind="in"); S.branch(ar2, "H2O", kind="out")
S.row([a, ar1, b, ar2, c]); print(S.save("arrow_branches.png", dpi=300)); S.check()
g2 = S._item_geo(ar2); tb = sorted([t["box"] for t in g2["texts"] if "".join(r[0] for r in t["runs"]) != "[Cu]"], key=lambda q: q[0])
print("branch labels below the arrow:", all(q[3] < ar2.y for q in tb), "| apart (BU):", round(float(tb[1][0] - tb[0][2]) / 14.4, 2),
      "| CO2 drawn above arrow 1:", S._mol_box(S.mols[-1] if False else ar1._bmols[0], curlies=False)[1] > ar1.y)

# Headings over structures (caption side="above"), a highlight behind a caption and a text, a frame round an equation
# with its label, and an empty labelled answer box.
S = Scene(Style())
m, n = from_smiles("Brc1ccccc1", align=(1, 0, 30)), from_smiles("OB(O)c1ccccc1", align=(1, 3, 150))
S.caption(m, "X = Cl, Br, I"); S.caption(m, "aryl halide", side="above")
S.caption(n, "boronic acid", side="above"); S.caption(n, "by-product: B(OH)3", highlight="#FFFC79")
r = S.row([m, "+", n], gap=1.5); S.frame(r, label="A"); S.frame(size=(6, 3), at=(4, -6), label="B")
S.text("activation", (12, -6), highlight="#FFFC79")
print(S.save("frames_headings.png", dpi=300)); S.check()
prims, texts = S._build(); tb = {"".join(q[0] for q in t["runs"]): t["box"] for t in texts if not t.get("atom") and "box" in t}
print("headings above, on one line:", tb["aryl halide"][1] > S._mol_box(m, curlies=False)[3], abs(tb["aryl halide"][1] - tb["boronic acid"][1]) < 0.01,
      "| fills:", sum(1 for p in prims if p[0] == "poly" and p[2] == "#FFFC79"), "| frames:", sum(1 for p in prims if p[0] == "line" and len(p[1]) == 6))

# Unequal equilibria: the disfavoured half is half as long, centred.
S = Scene(Style())
a1, a2 = S.arrow(kind="<=>>", above="fast"), S.arrow(kind="<<=>", below="slow")
S.row([from_smiles("CCO"), a1, from_smiles("CCOC"), a2, from_smiles("CCOCC")]); print(S.save("unequal_eq.png", dpi=300)); S.check()
half = lambda ar, k: (lambda P: round(float(P[:, 0].max() - P[:, 0].min()) / 14.4, 2))(np.asarray(S._item_geo(ar)["prims"][k][1], float))
print("halves (BU), forward then reverse:", half(a1, 0), half(a1, 1), "|", half(a2, 0), half(a2, 1))
