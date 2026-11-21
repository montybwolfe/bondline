# 1.1.0 fixes: every O-H between two labels shows the same length, arrows or not; a centred double bond (C=O or
# C=CH2) stops its neighbours at its near line even at a wide angle, and goes offset beside the wide end of a hash;
# the + and - of text carry the charge marks' weight at any size (an operator keeps the font's); text starts on
# whole pixels in a bitmap, so the same letter has the same pixels anywhere; an arrow no arc can draw clear of the
# ink gets a free curve.
import numpy as np
from bondline import *

def vis(m):                                                   # drawn length (pt) of each label-to-label bond
    return sorted(round(float(np.linalg.norm(b._prims[0][1][-1] - b._prims[0][1][0])), 2)
                  for b in m.bonds if b.a.drawn and b.b.drawn)

print("== water: the O-H an arrow leaves shows as much line as the other")
S = Scene(Style())
m = Mol(); h = m.atom("H", (0, 0)); o = m.attach(h, 40, "O"); m.attach(o, None, "H")
S.curly((h, o), o); S.row([m]); print(S.save("water.png", dpi=1200)); S.check()
print("O-H lines (pt):", vis(m))

print("== centred double bonds: the neighbours stop at the near line, never in the gap")
S = Scene(Style())
k = Mol(); c = k.atom("C", (0, 0)); k.attach(c, 90); k.attach(c, 210); k.attach(c, 300, "O", order=2)
ib = from_smiles("CC(=C)C"); iso = from_smiles("C[C@@H]1CC[C@H]([C@@H](C1)O)C(=C)C")
S.row([k, ib, iso], gap=1.5); print(S.save("centred.png", dpi=1200)); S.check()
for mol in (k, ib):
    C = [a for a in mol.atoms if len(mol.nbrs(a)) == 3][0]
    db = [b for b in mol.bonds if b.order == 2][0]
    ends = [np.asarray(pr[1][i], float) for pr in db._prims for i in (0, -1)]
    for n in mol.nbrs(C):
        bd = mol.find(C, n)
        if bd.order != 1: continue
        P = np.asarray(bd._prims[0][1], float); e = min(P, key=lambda q: np.linalg.norm(q - S.P(C)))
        print(C, "-", n, "meets a line of the double bond:", min(np.linalg.norm(e - q) for q in ends) < 0.05)
db = [b for b in iso.bonds if b.order == 2][0]
print("isopulegol's C=CH2 beside the hash is offset:", S._dside(db) != 0)

print("== signs in text: never lighter than the circled charges' bars")
S = Scene(Style())
for s, sz in (("+", S.st.size), ("+", 0.7 * S.st.size), ("−", 0.7 * S.st.size)):
    e = S.F.path(s, 0, 0, sz).get_extents(); bar = min(e.width, e.height) if s == "−" else None
    P = S.F._sign(s, 0, 0, sz); ys = P.vertices[:4, 1]
    print(s, round(sz, 1), "pt: bar", round(float(ys.max() - ys.min()), 3), "pt (marks:", round(S.MARK_W * S.st.lw, 3), ")")
P = S.F._sign("−", 0, 0, 0.7 * S.st.size, op=True); ys = P.vertices[:4, 1]
print("an operator (a − b in an exponent) keeps the font's weight:", round(float(ys.max() - ys.min()), 3), "pt")
S = Scene(Style())
a = Mol(); a.atom("H^+", (0, 0)); b = Mol(); b.atom("Na^+", (0, 0))
S.row([a, S.arrow(below="H^+ transfer"), b]); print(S.save("signs.png", dpi=1200)); S.check()

print("== text on whole pixels: one letter, the same pixels anywhere")
S = Scene(Style())
p = Mol(); pd = p.atom("Pd", (0, 0)); p.ligate(pd, "sqpl", N="L", S="L", E="Br", W="Ar")
S.row([p]); S.save("ls.png", dpi=260)
from PIL import Image
im = np.asarray(Image.open("ls.png").convert("L")); ys = np.flatnonzero((im < 128).any(1))
rows = np.split(ys, np.flatnonzero(np.diff(ys) > 1) + 1)
def glyph(r):
    sub = im[r[0]:r[-1] + 1]; xs = np.flatnonzero((sub < 250).any(0)); return sub[:, xs[0]:xs[-1] + 1]
top, bot = glyph(rows[0]), glyph(rows[-1])
print("top and bottom L identical:", top.shape == bot.shape and bool((top == bot).all()))

print("== free curves: the fallback for an arrow no arc draws clear of the ink")
S = Scene(Style())
m = Mol(); c = m.atom("C", (0, 0)); o = m.attach(c, 90, "O", order=2)
for t in (210, 330): m.attach(c, t)
nu = m.atom("HO", (0.0, -1.6)); nu.charge = -1
S.curly(nu.lp(90), c); S.row([m]); print(S.save("free.png", dpi=1200)); S.check()
F = S._free_cands(S.curlies[0])
print("free curves offered:", len(F) > 0)
print("all bend one way, 25-270 deg, no tighter than the arcs:",
      all(25 <= k["sw"] <= 270 and k["R"] >= 0.15 * S.st.bond for k in F))
print("placed clear:", S.curlies[0]._ok)
