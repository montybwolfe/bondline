# Catalytic cycle (Suzuki coupling): stations placed round a circle at the angles given, arcs drawn in order, reagents
# in and products out as branches on the arcs. ligate() puts ligands on standard sites at even label spacing.
from bondline import *
S = Scene(Style())
def pd(e, w):
    m = Mol(); p = m.atom("Pd", (0, 0)); m.ligate(p, "sqpl", N="L", S="L", E=e, W=w); return m
pd0 = Mol(); p = pd0.atom("Pd^0", (0, 0)); pd0.ligate(p, "lin", E="L", W="L")
c = S.cycle([(pd0, 90), (pd("Br", "Ar"), -30), (pd("Ar'", "Ar"), 210)], radius=4)
S.branch(c, 0, "Ar–Br", kind="in"); S.branch(c, 1, "Ar'B(OH)_2, base", kind="in"); S.branch(c, 2, "Ar–Ar'", kind="out")
print(S.save("suzuki_cycle.png")); S.check()
