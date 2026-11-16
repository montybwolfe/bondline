# P0-4: reaction-arrow text is measured (the arrow is always longer than its text), stacked lines, symbols in the
# house font, vector arrows in text, and no tofu: a glyph the font lacks is reported by check() and raised by save().
from bondline import *
st = Style()
def ring():
    m = Mol(); m.ring(6, kekule=True); return m
print("== long and stacked conditions")
S = Scene(st)
a1 = S.arrow(above="1. LiAlH_4, THF, 0 °C → rt, 12 h", below="2. H_3O^+ workup\nthen NaHCO_3")
a2 = S.arrow(above=["1. LDA, THF, -78 °C", "2. MeI"], below="Δ, hν", step=1)
S.row([ring(), a1, ring(), a2, ring()]); print(S.save("conditions.png")); S.check()
for k, a in enumerate((a1, a2)):
    g = S._item_geo(a); L = g["prims"][0][1][-1][0] - g["prims"][0][1][0][0]
    print(f"arrow {k}: text fits inside the shaft:", all(t["box"][0] >= a.x and t["box"][2] <= a.x + L + 0.8 * st.head_len for t in g["texts"]))
print("== symbols, equilibrium glyph, vertical arrow with stacked text")
S = Scene(st)
x = ring(); y = ring(); z = ring()
S.place(x, (0, 0)); S.arrow(vertical=True, above="hν\nCH_2Cl_2", below="−78 °C\n2 h", at=(0, -1.8)); S.place(y, (0, -8.5))
S.text("<i>K</i> ⇌ ‡ δ+ δ− 5 ℃ H₃O⁺", at=(7.5, -4))
print(S.save("symbols.png")); S.check()
print("== negative: a glyph the house font lacks")
S = Scene(st); m = Mol(); n = m.atom("N⊕", (0, 0)); m.attach(n, 0, "Me"); S.row([m, S.arrow(above="x ⇄ y")]); S.check()
try: S.save("tofu.png"); print("saved (should not happen)")
except ValueError as e: print("save raised:", str(e)[:60])
