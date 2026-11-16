# Teaching-figure pass: orbital lobes on a labelled atom start at the letter's edge (whole,
# the label clear), a caption clears the lobes, a p orbital counts as the pi bond in check(), a sigma frame drawn with
# wedges round an orbital atom is perspective; dashed curly arrows; tidy() folds CH2-H into CH3; the stock H-Br keeps
# a readable bond wherever protonate() turns it; halonium and borate valences; generic HY/HNu/HX labels bond through
# Y; a colon before a word is punctuation; nested sub/superscripts; text stations in a cycle; the several-steps arrow;
# arrow text clear of an equilibrium arrow's raised barb.
from bondline import *

S = Scene(Style())
m = Mol(); c = m.atom("C", (0, 0)); o = m.atom("O", (1, 0)); m.bond(c, o)
m.attach(c, 150, style="hash"); m.attach(c, 210, style="bold")
S.orbital(c, "p", 90, shade="up"); S.orbital(o, "p", 90, shade="up", size=1.6)
S.caption(m, "π (C=O)")
n = Mol(); c2 = n.atom("C", (0, 0)); n2 = n.atom("N", (1, 0)); n.bond(c2, n2); n.attach(n2, 0, "H", style="hash")
n.attach(c2, 150, style="hash"); n.attach(c2, 210, style="bold")
S.orbital(c2, "p", 90, shade="up", color="#00897B"); S.orbital(n2, "p", 90, shade="up", size=1.25, color="#00897B")
S.caption(n, "π (C=N)")
S.row([m, n]); print(S.save("orbitals.png", dpi=300)); S.check()
low = lambda x: min(np.asarray(pr[1])[:, 1].min() for ex in S.extras if ex[0] == "orbital" and ex[1].mol is x
                    for pr in S._orbital_prims(ex) if pr[0] == "orb")
print("structure boxes hold their lobes (captions go below them):", all(S._mol_box(x)[1] <= low(x) + 1e-6 for x in (m, n)))

S = Scene(Style())
k = from_smiles("CC(C)=O", align=(1, 3, 90))
nu = base("Nu", charge=0); k.attack(nu, k[1])
S.curly(nu.key.lps[0], k[1], dashed=True); S.curly((k[1], k[3]), k[3], dashed=True)
al = Mol(); a0 = al.atom("C", (0, 0)); al.attach(a0, 210); al.attach(a0, -30); x = al.attach(a0, 90, "CH2", 2)
hb = acid("Br"); al.protonate(hb, x); S.curly((a0, x), hb.key); S.curly((hb.key, hb.a), hb.a)
S.row([k, al]); print(S.save("dashed_hbr.png", dpi=300)); S.check()
t = Mol(); c = t.atom("C", (0, 0)); t.attach(c, 210); t.attach(c, -30)
y = t.attach(c, 90, "CH2"); t.attach(y, 30, "H"); c.charge = 1; t.tidy(carbon=True)
print("tidy:", [a.label for a in t.atoms])

S = Scene(Style())
br = Mol(); c1 = br.atom("C", (0, 0)); c2 = br.attach(c1, 0); b = br.atom("Br", (0, 0)); br.dock(b, c1, 60, 1.0)
br.bond(c1, b); br.bond(c2, b); b.charge = 1; br.attach(c1, 210); br.attach(c2, -30)
bh = Mol(); B = bh.atom("B", (0, 0)); B.charge = -1
for t_ in (45, 135, 225, 315): bh.attach(B, t_, "H")
hy = []
for lab in ("HY", "HNu", "HX"):
    g = Mol(); q = g.atom(lab, (0, 0)); g.attach(q, 0, "R2", length=1.3, free=True); hy.append(g)
dz = from_smiles("O=[N+]([O-])c1ccc(cc1)C", labels={9: "N2Cl"}, align=(3, 6, 0))
S.row([br, bh] + hy + [dz]); print(S.save("valence_labels.png", dpi=300)); S.check()
print("labels:", [g.atoms[0].label for g in hy])

S = Scene(Style())
S.text("M = SnBu3: Stille coupling\nY = NH, O: Buchwald-Hartwig\nB: + H-A", ha="left")
S.text("<i>K</i>_{eq} = 10^{(p<i>K</i>_a(H_3O^+) − p<i>K</i>_a(HA))}", (0, -3.2))
print(S.save("text.png", dpi=300)); S.check()

S = Scene(Style())
c = S.cycle([("L_nPd^0", 90), ("<b>B</b>", -30), ("<b>D</b>", 210)], radius=3)
print("cycle radius:", c["r"], "texts:", sum(isinstance(i, Text) for i in S.items))
print(S.save("cycle_text.png", dpi=300)); S.check()

S = Scene(Style())
S.row(["HA + H2O", S.arrow(above="<i>K</i>eq", kind="<=>"), "H3O⊕ + A⊖"])
print(S.save("eq_text.png", dpi=300)); S.check()
S = Scene(Style())
a = from_smiles("O=C1CCCCC1"); b = from_smiles("OC1CCCCC1")
S.row([a, S.arrow(kind="->>", above="several", below="steps"), b]); print(S.save("multistep.png", dpi=300)); S.check()
