# Workflow B (redraw a sketch): a handwritten Blanc chloromethylation of 4-nitrophenetole, rebuilt step by step:
# EAS on protonated formaldehyde (arenium cation beside OEt; the leaving ring H drawn next to it, the CH2OH on the
# roomier side), re-aromatisation by chloride, protonation of the benzylic alcohol by HCl (lone pair above the O, clear
# of its H), SN2 by chloride (near anti: OEt blocks the exact line), water leaving; a snake layout (row under row via a
# vertical arrow, the second row read right to left). tidy(carbon=True), unfold(), unbond() returning a ring substituent
# to the bisector, join() choosing the roomier orientation.
from bondline import *
S = Scene(Style())
# 1. phenetole
ph = from_smiles("CCOc1ccccc1", align=(3, 2, 90), abbrev=True)
# 2. 4-nitrophenetole, built so its atoms carry through every step (ring C1 top ... C4 bottom, NO2 on C4)
s2 = Mol(); r = s2.ring(6, rot=90)
for nm, a in zip(("c1", "c6", "c5", "c4", "c3", "c2"), r): a.name = nm
for x, y in (("c1", "c2"), ("c3", "c4"), ("c5", "c6")): s2.find(s2[x], s2[y]).order = 2
s2.attach(s2["c1"], 90, "OEt", name="oet"); s2.attach(s2["c4"], 270, "NO2", name="no2")
f = Mol(); fc = f.atom("C", (0, 0), name="fc"); fo = f.attach(fc, 90, "O", 2, name="fo"); fo.charge = 1
f.attach(fo, 30, "H", name="foh").deliberate = True; f.attach(fc, 210, "H", name="fh1"); f.attach(fc, 330, "H", name="fh2")
s2.merge(f); s2.attack(s2["c2"], s2["fc"])
S.curly((s2["c1"], s2["c2"]), s2["fc"]); S.curly((s2["fc"], s2["fo"]), s2["fo"])
# 3. arenium ion (cation at C1, next to OEt); chloride takes the ring proton
s3 = s2.copy()
s3.find(s3["c1"], s3["c2"]).order = 1; s3["c1"].charge = 1
s3.find(s3["fc"], s3["fo"]).order = 1; s3["fo"].charge = 0
s3["foh"].deliberate = False; s3.tidy(carbon=True)                                   # H2C=O+H drawn with its H atoms becomes a CH2-OH chain
s3.join(s3["c2"], s3["fc"]); hr = s3.attach(s3["c2"], None, "H", name="hr")
cl = halide("Cl"); s3.merge(cl); s3.deprotonate(cl, s3["hr"])
S.curly((s3["c2"], s3["hr"]), (s3["c1"], s3["c2"])); S.curly(cl.key.lps[0] if cl.key.lps else cl.key, s3["hr"])
S.validate(s2, s3)
# 4. the benzylic alcohol, its OH taking the proton of the HCl just formed
s4 = s3.copy()
s4.find(s4["c1"], s4["c2"]).order = 2; s4["c1"].charge = 0
clx = [a for a in s4.atoms if a.label == "Cl"][0]; clx.charge = 0
s4.unbond(s4["c2"], s4["hr"]); s4.join(clx, s4["hr"])
s4.protonate(s4["hr"], s4["fo"])
S.validate(s3, s4)
S.curly(s4["fo"].lp(), s4["hr"]); S.curly((s4["hr"], clx), clx)
# 5. the oxonium ion; chloride attacks the CH2 (SN2), water leaves
s5 = s4.copy(); ox = s5["fo"]; s5.unfold(ox)
cl5 = [a for a in s5.atoms if a.label == "Cl"][0]
s5.unbond(s5["hr"], cl5); s5.join(ox, s5["hr"]); ox.charge = 1; cl5.charge = -1
s5.attack(cl5, s5["fc"])
S.validate(s4, s5)
S.curly(cl5.lps[0] if cl5.lps else cl5, s5["fc"]); S.curly((s5["fc"], ox), ox)
# 6. 2-(chloromethyl)-1-ethoxy-4-nitrobenzene (+ water)
s6 = s5.copy(); c6 = [a for a in s6.atoms if a.label == "Cl"][0]
s6.unbond(s6["fc"], s6["fo"]); s6.join(s6["fc"], c6); c6.charge = 0; s6["fo"].charge = 0
w = s6.split(s6["fo"])
S.validate(s5, [s6, w])
r1 = S.row([ph, S.arrow(above=["HNO3", "H2SO4"], below="(standard nitration)"), s2, S.arrow(), s3, S.arrow(), s4])
S.row([s6, S.arrow(kind="<-", below="−H2O"), s5], under=r1, via=S.arrow(vertical=True))
S.balance(s5, s6)
print(S.save("chloromethylation.png", dpi=300)); S.check()
ag = lambda m, a, b: round(ang(m[b].xy - m[a].xy)) % 360
print("arenium: CH2 at", ag(s3, "c2", "fc"), "deg, its H at", ag(s3, "c2", "hr"), "deg (exterior bisector 30; H beside C1)")
print("alcohol: CH2 at", ag(s4, "c2", "fc"), "deg (back on the bisector once the H has gone)")
print("SN2: Cl- comes in at", round(ang(cl5.xy - s5["fc"].xy)) % 360, "deg; exactly anti would be", round(ang(s5["fc"].xy - ox.xy)) % 360, "(OEt in the way)")
lp = s4["fo"].lps[0]; print("alcohol O lone pair at", round(lp.angle) % 360, "deg (above the O, clear of its H)")
