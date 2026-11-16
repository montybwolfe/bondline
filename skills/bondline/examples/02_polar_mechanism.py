# Polar mechanism, written forwards with one Mol per step: ester hydrolysis by hydroxide (saponification).
# attack() and deprotonate() put each reagent at the standard distance and angle; join() and release() make and break
# bonds and lay the result out; validate() replays every step's arrows; the products are tidied (OH labels).
from bondline import *
S = Scene(Style())
s1 = Mol(); c = s1.atom("C", (0, 0), name="C"); s1.attach(c, 210, "Me"); s1.attach(c, 90, "O", 2, name="Oc")
Or = s1.attach(c, 330, "O", name="Or"); s1.attach(Or, 30, "Me")                     # methyl acetate
oh = hydroxide(explicit=True); oh.key.name = "Ob"; oh.atoms[1].name = "H"           # its H is used in step 3
s1.attack(oh, c); S.curly(s1["Ob"].lp(), c); S.curly((c, s1["Oc"]), s1["Oc"])
s2 = s1.copy(); s2.join(s2["C"], s2["Ob"]); s2["Ob"].charge = 0; s2.find(s2["C"], s2["Oc"]).order = 1; s2["Oc"].charge = -1
S.curly(s2["Oc"].lp(), (s2["Oc"], s2["C"])); S.curly((s2["C"], s2["Or"]), s2["Or"])  # tetrahedral intermediate collapses
s3 = s2.copy(); meo = s3.release(s3["Or"], s3["C"]); s3.find(s3["C"], s3["Oc"]).order = 2; s3["Oc"].charge = 0
meo["Or"].charge = -1; s3.deprotonate(meo["Or"], s3["H"])                          # methoxide takes the acid's proton
S.curly(s3["Or"].lp(), s3["H"]); S.curly((s3["H"], s3["Ob"]), s3["Ob"])
s4 = s3.copy(); s4.unbond(s4["H"], s4["Ob"]); s4.join(s4["Or"], s4["H"]); s4["Ob"].charge = -1; s4["Or"].charge = 0
meoh = s4.split(s4["Or"]); s4.tidy(); meoh.tidy()
S.validate(s1, s2); S.validate(s2, [s3, meo]); S.validate(s3, [s4, meoh]); Scene.continuity(s1, s2, s3, [s4, meoh])
S.flow([s1, S.arrow(), s2, S.arrow(), s3, S.arrow(), s4, "+", meoh], width_cm=16)
print(S.save("saponification.png")); S.check()
