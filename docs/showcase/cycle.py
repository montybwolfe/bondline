# The skill's examples/08_catalytic_cycle.py with every PPh3 drawn out: ring_out() hangs three phenyls off each P,
# fanned 80 deg apart, and the two cis phosphines of the transmetallation product sit further out, their fans turned
# away from each other. Suzuki coupling of 4-bromoanisole with phenylboronic acid; the stations go round a circle at
# the angles given, reagents in and products out as branches (a branch can be a structure), each step named by its arc.
from bondline import *
S = Scene(Style())
def ph3(m, sites, turn={}):                       # three drawn phenyls on each P, fanned 80 deg apart (turn: per site)
    for k, P in sites.items():
        if P.label != "P": continue
        t = ang(P.xy - m.nbrs(P)[0].xy) + turn.get(k, 0)
        for d in (-80, 0, 80): m.ring_out(P, t + d, kekule=True)
def anisyl_pd():                                      # 4-MeOC6H4-Pd, its C-Pd bond pointing right
    smi = "COc1ccc(cc1)[Pd]"; m = from_smiles(smi); pd = m.atoms[-1]; ip = m.nbrs(pd)[0]
    m = from_smiles(smi, align=(m.atoms.index(ip), m.atoms.index(pd), 0)); return m, m.atoms[-1]
oa, p1 = anisyl_pd(); ph3(oa, oa.ligate(p1, "sqpl", length=1.8, N="P", S="P", E="Br"))
tm, p2 = anisyl_pd(); ph3(tm, tm.ligate(p2, "sqpl", length=2.5, E="P", S="P"), {"E": 20, "S": -20})
ph = from_smiles("c1ccccc1"); ipso = ph.atoms[0]; tm.merge(ph); tm.join(p2, ipso, angle=90)
pd0 = Mol(); p = pd0.atom("Pd^0", (0, 0)); ph3(pd0, pd0.ligate(p, "lin", length=1.5, E="P", W="P"))
c = S.cycle([(pd0, 90), (oa, -30), (tm, 210)], radius=8)
S.branch(c, 0, from_smiles("COc1ccc(Br)cc1"), kind="in"); S.branch(c, 1, "PhB(OH)_2, K_2CO_3", kind="in")
S.branch(c, 2, from_smiles("COc1ccc(cc1)c1ccccc1"), kind="out")
S.branch(c, 0, "oxidative\naddition", kind="step"); S.branch(c, 1, "transmetallation", kind="step", side="in")
S.branch(c, 2, "reductive\nelimination", kind="step")
print(S.save("cycle.png", dpi=260)); S.check()
