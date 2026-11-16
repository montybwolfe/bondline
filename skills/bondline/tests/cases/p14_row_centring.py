# P1-4: rows centre each structure on its main fragment (captions, curly arrows and docked reagents do not shift
# the backbone off the reaction-arrow line) and space structures by a constant gap.
from bondline import *
st = Style()
def ethanal():
    m = Mol(); c = m.atom("C", (0, 0)); c2 = m.attach(c, 30); m.attach(c2, -30, "O", 2); return m
print("== captions")
S = Scene(st); a, b, p = ethanal(), ethanal(), ethanal()
S.caption(a, "ethanal"); S.caption(p, "a much longer caption")
arw = S.arrow(above="NaOH (aq)")
S.row([a, "+", b, arw, p]); print(S.save("captions.png"))
for m in (a, b, p):
    bx = S._mol_box(m, curlies=False); print("structure centre - arrow line:", round((bx[1] + bx[3]) / 2 - arw.y, 3), "pt")
S.check()
print("== docked reagent below, curly arrows")
S = Scene(st)
m = Mol(); c1 = m.atom("C", (0, 0), name="c1"); c2 = m.attach(c1, 30, name="c2"); O = m.attach(c2, -30, "O", 2, name="O")
H = m.attach(c1, -90, "H", name="H"); B = m.atom("O", (0, 0), name="B"); B.charge = -1; m.attach(B, 210, "H"); m.dock(B, H, -90)
S.curly(B.lp(), H); S.curly((H, c1), (c1, c2)); S.curly((c2, O), O)
arw = S.arrow(); q = ethanal(); S.row([m, arw, q]); print(S.save("docked.png"))
fr = m.frag(c1); bx = S._mol_box(m, curlies=False, atoms=fr)
print("main fragment centre - arrow line:", round((bx[1] + bx[3]) / 2 - arw.y, 3), "pt")
sm = S._mol_box(m, curlies=False); g = S._item_geo(arw)["box"]; sq = S._mol_box(q, curlies=False)
print("gaps (BU):", round((g[0] - sm[2]) / st.bond, 2), round((sq[0] - g[2]) / st.bond, 2))
S.check()
