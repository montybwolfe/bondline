# The examples in SKILL.md sections 5 and 8 of the pre-revision skill
from bondline import *
print("== proline")
st = Style(); S = Scene(st)
pro = from_smiles("OC(=O)[C@@H]1CCCN1", align=(3, 1, 0)); stereo_report(pro)
S.caption(pro, "(<i>S</i>)-proline"); S.row([pro]); print(S.save("proline.png")); S.check()

print("== lactone protonation")
st = Style(); S = Scene(st)
s1 = Mol(); r = s1.ring(6, rot=90, labels={5: "O"}); Oc = s1.attach(r[0], 90, "O", 2, name="Oc")
Ht = s1.atom("H", Oc.xy + dirv(150) * 2.2, name="Ht")
Oa = s1.attach(Ht, 150, "O", name="Oa"); Oa.charge = 1; s1.attach(Oa, None, "H"); s1.attach(Oa, None, "Me")
S.curly(Oc.lp(), Ht); S.curly((Ht, Oa), Oa)
s2 = s1.copy(); s2.transfer(s2["Ht"], s2["Oc"], 150); s2["Oc"].charge = 1; s2.remove_frag(s2["Oa"])
S.validate(s1, s2); S.row([s1, S.arrow(above="MeOH"), s2]); print(S.save("step.png")); S.check()

print("== Pd cycle from one template")
st = Style(); S = Scene(st)
def pd_complex(X):
    m = Mol(); P = m.atom("Pd", (0, 0)); m.ligate(P, "sqpl", N="L", E=X)
    m.ring_out(P, 180, 6, kekule=True, length=1.3); return m
pd0 = Mol(); P0 = pd0.atom("Pd^0", (0, 0)); pd0.ligate(P0, "lin", E="L", W="L")
oxadd = pd_complex("Br"); amido = pd_complex("NHPh")
cyc = S.cycle([(pd0, 90), (oxadd, -30), (amido, 210)], radius=4.2)
S.branch(cyc, 0, "<b>A</b>  Ar–Br", kind="in")
S.branch(cyc, 1, "<i>t</i>BuOH + NaBr", kind="out", at=0.8)
print(S.save("cycle.png")); S.check()
