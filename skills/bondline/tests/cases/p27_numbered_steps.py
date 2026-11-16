# Workflow C (build from the user's step list): a numbered reference drawing, atoms found by their numbers in later
# steps (copies keep them, a merged reagent keeps its own), and a step built exactly from a list ("lp on 7 -> 6;
# bond 1-6 -> 1-2; bond 2-3 -> 3"). Workflow B's side_by_side() puts an original beside the redraw at one height.
from bondline import *
from PIL import Image

print("== numbered reference drawings")
S = Scene(Style())
k = Mol(); ca = k.atom("C", (0, 0)); cc = k.attach(ca, 30); k.attach(cc, 90, "O", 2); k.attach(cc, -30, "OEt")
k.attach(ca, 150); k.attach(ca, None, "H")
b = base()
S.row([k.numbered(), b.numbered(7)], gap=1.5)
print(S.save("numbered.png", dpi=300)); S.check()
n2 = from_smiles("CCOc1ccc(cc1)[N+](=O)[O-]", align=(2, 3, 0))
S = Scene(Style()); S.row([n2.numbered()]); print(S.save("numbered_smiles.png", dpi=300)); S.check()

print("== the user's steps: lp on 7 -> 6; bond 1-6 -> 1-2; bond 2-3 -> 3")
S = Scene(Style())
m = k.copy(); m.merge(b); m.deprotonate(b.key, m.num(6))
S.curly(m.num(7).lps[0], m.num(6)); S.curly((m.num(1), m.num(6)), (m.num(1), m.num(2))); S.curly((m.num(2), m.num(3)), m.num(3))
e = m.copy(); e.unbond(e.num(1), e.num(6)); e.join(e.num(7), e.num(6)); e.num(7).charge = 0
e.find(e.num(1), e.num(2)).order = 2; e.find(e.num(2), e.num(3)).order = 1; e.num(3).charge = -1
bh = e.split(e.num(7))
print("atoms 1-7 found in the second step:", all(e.num(i) for i in range(1, 6)), bh.num(7).label, bh.num(6).label)
S.validate(m, [e, bh]); S.balance(m, [e, bh])
S.row([m, S.arrow(kind="<=>"), e, "+", bh]); print(S.save("steps.png", dpi=300)); S.check()
try: e.num(42)
except KeyError as x: print("num(42):", x)

print("== side_by_side: an original beside the redraw, one height, optional crop")
out = side_by_side("numbered.png", "steps.png", "compare.png", height=300)
print(Image.open(out).size[1], "px high")
out = side_by_side("numbered.png", "steps.png", "compare_crop.png", crop=(0.0, 0.0, 0.5, 1.0), height=300)
print(Image.open(out).size[1], "px high (cropped)")
