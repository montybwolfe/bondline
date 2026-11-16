# Aldol figure 3 (crossed aldol), script as first written
from bondline import *
st = Style(); S = Scene(st)

# product: Ph-CH(O)-CH2-C(=O)-CH3
P = Mol()
C1b = P.atom("C", (0, 0), name="C1b"); O1b = P.attach(C1b, 90, "O", name="O1b")
P.ring_out(C1b, 210, 6, kekule=True)
Ca = P.attach(C1b, -30, "C", name="Ca"); C1 = P.attach(Ca, 30, "C", name="C1")
O1 = P.attach(C1, 90, "O", 2, name="O1"); P.attach(C1, -30, "C", name="Me")

s3 = P.copy(); s3["O1b"].charge = -1
Hw = s3.atom("H", s3["O1b"].xy + dirv(90) * 2.2, name="Hw")
Ow = s3.attach(Hw, 90, "O", name="Ow"); s3.attach(Ow, None, "H")

s2 = P.copy(); s2.unbond(s2["Ca"], s2["C1b"])
s2.bond(s2["C1b"], s2["O1b"], 2); s2.bond(s2["Ca"], s2["C1"], 2); s2.bond(s2["C1"], s2["O1"], 1)
s2["O1"].charge = -1
s2.dock(s2["C1b"], s2["Ca"], 150)

s1 = s2.copy(); s1.remove_frag(s1["C1b"])
s1.bond(s1["Ca"], s1["C1"], 1); s1.bond(s1["C1"], s1["O1"], 2); s1["O1"].charge = 0
Ha = s1.attach(s1["Ca"], -90, "H", name="Ha")
Ob = s1.atom("O", Ha.xy + dirv(-90) * 2.2, name="Ob"); Ob.charge = -1
s1.attach(Ob, 210, "H")

s4 = s3.copy(); s4.transfer(s4["Hw"], s4["O1b"], 90); s4["O1b"].charge = 0
s4["Ow"].charge = -1; s4.remove_frag(s4["Ow"])
hx = Mol(); ox = hx.atom("O", (0, 0)); ox.charge = -1; hx.attach(ox, 210, "H")

S.curly(s1["Ob"].lp(), s1["Ha"]); S.curly((s1["Ca"], s1["Ha"]), (s1["Ca"], s1["C1"])); S.curly((s1["C1"], s1["O1"]), s1["O1"])
S.curly(s2["O1"].lp(), (s2["O1"], s2["C1"])); S.curly((s2["Ca"], s2["C1"]), s2["C1b"]); S.curly((s2["C1b"], s2["O1b"]), s2["O1b"])
S.curly(s3["O1b"].lp(), s3["Hw"]); S.curly((s3["Hw"], s3["Ow"]), s3["Ow"])

for a, b in ((s1, s2), (s2, s3), (s3, s4)): S.validate(a, b)
S.flow([s1, S.arrow(below="−H_2O", step=1), s2, S.arrow(step=2), s3, S.arrow(step=3), s4, "+", hx], width_cm=17)
print(S.save("crossed_aldol_mechanism.png")); S.check()
