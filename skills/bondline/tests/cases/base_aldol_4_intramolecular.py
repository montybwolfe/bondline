# Aldol figure 4 (intramolecular aldol condensation), final script as first written
from bondline import *
st = Style(); S = Scene(st)

# ring product skeleton (3-hydroxycyclohexanone / alkoxide): vertices at 30..330 deg
R = Mol(); ang_of = {"C1": 30, "C2": 90, "C3": 150, "C4": 210, "C5": 270, "C6": 330}
at = {k: R.atom("C", dirv(t), name=k) for k, t in ang_of.items()}
for a, b in (("C1", "C2"), ("C2", "C3"), ("C3", "C4"), ("C4", "C5"), ("C5", "C6"), ("C6", "C1")): R.bond(at[a], at[b])
R.attach(at["C2"], None, "O", 2, name="O2"); R.attach(at["C6"], None, "O", name="O6")

TH = 20
def open_chain(m):
    m.unbond(m["C1"], m["C6"]); m.bond(m["C6"], m["O6"], 2)
    m.swing(m["C5"], m["C4"], -10); m.swing(m["C6"], m["C5"], -10); m.swing(m["C1"], m["C2"], TH); return m

# step 3: alkoxide + water
s3 = R.copy(); s3["O6"].charge = -1
Hw = s3.atom("H", s3["O6"].xy + dirv(330) * 2.2, name="Hw"); Ow = s3.attach(Hw, 330, "O", name="Ow"); s3.attach(Ow, None, "H")
# step 2: open-chain enolate
s2 = open_chain(R.copy()); s2.bond(s2["C1"], s2["C2"], 2); s2.bond(s2["C2"], s2["O2"], 1); s2["O2"].charge = -1
# step 1: 5-oxohexanal + hydroxide
s1 = open_chain(R.copy())
H1 = s1.attach(s1["C1"], None, "H", name="H1"); Ob = s1.atom("O", H1.xy + unit(H1.xy - s1["C1"].xy) * 2.2, name="Ob")
Ob.charge = -1; s1.attach(Ob, None, "H")
# step 4: aldol (3-hydroxycyclohexanone) + hydroxide at the alpha-H between C=O and C-OH
s4 = s3.copy(); s4.transfer(s4["Hw"], s4["O6"], 330); s4["O6"].charge = 0; s4["Ow"].charge = -1
H1r = s4.attach(s4["C1"], None, "H", name="H1r"); s4.dock(s4["Ow"], H1r, 30)
# step 5: ring enolate (water gone)
s5 = s4.copy(); s5.transfer(s5["H1r"], s5["Ow"]); s5["Ow"].charge = 0; s5.remove_frag(s5["Ow"])
s5.bond(s5["C1"], s5["C2"], 2); s5.bond(s5["C2"], s5["O2"], 1); s5["O2"].charge = -1
# step 6: cyclohex-2-enone + hydroxide
s6 = s5.copy(); s6.bond(s6["C1"], s6["C2"], 1); s6.bond(s6["C2"], s6["O2"], 2); s6["O2"].charge = 0
s6.bond(s6["C1"], s6["C6"], 2)
s6.unbond(s6["C6"], s6["O6"]); s6.remove_frag(s6["O6"])
hx = Mol(); ox = hx.atom("O", (0, 0)); ox.charge = -1; hx.attach(ox, 210, "H")

S.curly(s1["Ob"].lp(), s1["H1"]); S.curly((s1["C1"], s1["H1"]), (s1["C1"], s1["C2"])); S.curly((s1["C2"], s1["O2"]), s1["O2"])
S.curly(s2["O2"].lp(), (s2["O2"], s2["C2"])); S.curly((s2["C1"], s2["C2"]), s2["C6"]); S.curly((s2["C6"], s2["O6"]), s2["O6"])
S.curly(s3["O6"].lp(), s3["Hw"]); S.curly((s3["Hw"], s3["Ow"]), s3["Ow"])
S.curly(s4["Ow"].lp(), s4["H1r"]); S.curly((s4["C1"], s4["H1r"]), (s4["C1"], s4["C2"])); S.curly((s4["C2"], s4["O2"]), s4["O2"])
S.curly(s5["O2"].lp(), (s5["O2"], s5["C2"])); S.curly((s5["C1"], s5["C2"]), (s5["C1"], s5["C6"])); S.curly((s5["C6"], s5["O6"]), s5["O6"])

for a, b in ((s1, s2), (s2, s3), (s3, s4), (s4, s5), (s5, s6)): S.validate(a, b)
Scene.continuity(s1, s2, s3, s4, s5, s6)
S.flow([s1, S.arrow(below="−H_2O", step=1), s2, S.arrow(step=2), s3, S.arrow(step=3), s4,
        S.arrow(below="−H_2O", step=4), s5, S.arrow(step=5), s6, "+", hx], width_cm=17)
print(S.save("intramolecular_aldol_mechanism.png")); S.check()
