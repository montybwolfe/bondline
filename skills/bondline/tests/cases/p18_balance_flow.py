# Stage 7: Scene.balance (atoms gained and lost between steps, matched against the arrow text), formula(), and flow():
# accurate row widths, a lone trailing item kept by shrinking gaps, a continuation that starts with an arrow indented
# under the first arrow of the row above.
from bondline import *

print("== formulas")
for name, m in [("ethanol", from_smiles("CCO")), ("acetate", from_smiles("CC(=O)[O-]")), ("benzene", from_smiles("c1ccccc1"))]:
    c, q = formula(m); print(f"  {name}: {hill(c)} ({q:+d})")
f = Mol(); f.sandwich("Fe"); print("  ferrocene:", hill(formula(f)[0]))
cp = Mol(); cp.deloc(cp.ring(5, rot=90), charge=-1); c, q = formula(cp); print(f"  Cp-: {hill(c)} ({q:+d})")
for lab in ["Et3N", "iPr2NEt", "PPh3", "LiAlH4", "tBuOK", "NHBoc", "Pd(PPh3)4", "HBr", "NH3", "LDA"]:
    print(f"  {lab}: {hill(label_formula(lab))}")

print("== balance")
S = Scene(); a = from_smiles("CC(O)C"); b = from_smiles("C=CC"); S.add(a); S.arrow("H2SO4", "−H_2O"); S.add(b)
print(" ", S.balance(a, b))
S = Scene(); S.add(a); S.arrow("H2SO4"); S.add(b); S.balance(a, b)
S = Scene(); h = from_smiles("CC(Br)(C)C"); e = from_smiles("C=C(C)C"); S.add(h); S.arrow("EtOH, Δ"); S.add(e)
S.balance(h, e); S.balance(h, [e, halide("Br"), hydronium(explicit=False)])

print("== flow: a lone trailing item stays on the row by closing the gaps a little")
def seq(S):
    its = []
    for k in range(3): its += [from_smiles("CC(=O)C"), S.arrow("step")]
    return its + [from_smiles("CCO")]
S = Scene(); its = S.row(seq(S)); bx = [S._ibox(i) for i in its]
natural = (max(b[2] for b in bx) - min(b[0] for b in bx)) / 72 * 2.54          # cm, one row at normal gaps
S = Scene(); rows = S.flow(seq(S), width_cm=0.97 * natural, squeeze=None); print("  without squeeze:", [len(r) for r in rows])
S = Scene(); rows = S.flow(seq(S), width_cm=0.97 * natural); print("  with squeeze:", [len(r) for r in rows])
S.save("squeeze.png")
print("== flow: a continuation that must start with an arrow lines up under the first arrow above")
S = Scene(); its = []
for k in range(4):
    its += [from_smiles("c1ccccc1C(=O)C"), S.arrow(step=k + 1)]
its += [from_smiles("CCO"), "+", halide("Br")]
rows = S.flow(its, width_cm=12)
arr = [[round(S._ibox(i)[0], 2) for i in r if isinstance(i, RArrow)] for r in rows]
print("  rows:", [len(r) for r in rows], "first arrow x per row:", [a_[0] if a_ else None for a_ in arr])
S.check(); S.save("indent.png")
