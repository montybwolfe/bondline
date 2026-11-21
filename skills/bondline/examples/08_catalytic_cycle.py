# Catalytic cycle (Suzuki coupling of 4-bromoanisole with phenylboronic acid): stations placed round a circle at the
# angles given, arcs drawn in order, reagents in and products out as branches (a branch can be a structure), and each
# step named beside its arc. ligate() puts the PPh3 ligands on standard square-planar sites; join() adds the phenyl
# cis to the aryl group, ready for reductive elimination.
from bondline import *
S = Scene(Style())
def anisyl_pd():                                      # 4-MeOC6H4-Pd, its C-Pd bond pointing right
    smi = "COc1ccc(cc1)[Pd]"; m = from_smiles(smi); pd = m.atoms[-1]; ip = m.nbrs(pd)[0]
    m = from_smiles(smi, align=(m.atoms.index(ip), m.atoms.index(pd), 0)); return m, m.atoms[-1]
oa, p1 = anisyl_pd(); oa.ligate(p1, "sqpl", N="PPh_3", S="PPh_3", E="Br")
tm, p2 = anisyl_pd(); tm.ligate(p2, "sqpl", E="PPh_3", S="PPh_3")
ph = from_smiles("c1ccccc1"); ipso = ph.atoms[0]; tm.merge(ph); tm.join(p2, ipso, angle=90)
pd0 = Mol(); p = pd0.atom("Pd^0", (0, 0)); pd0.ligate(p, "lin", E="PPh_3", W="PPh_3")
c = S.cycle([(pd0, 90), (oa, -30), (tm, 210)], radius=5.5)
S.branch(c, 0, from_smiles("COc1ccc(Br)cc1"), kind="in"); S.branch(c, 1, "PhB(OH)_2, K_2CO_3", kind="in")
S.branch(c, 2, from_smiles("COc1ccc(cc1)c1ccccc1"), kind="out")
S.branch(c, 0, "oxidative\naddition", kind="step"); S.branch(c, 1, "transmetallation", kind="step")
S.branch(c, 2, "reductive\nelimination", kind="step")
print(S.save("suzuki_cycle.png")); S.check()
