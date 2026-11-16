# Catalytic cycle with branches
from bondline import *
st = Style(); S = Scene(st)
def pd(X, Y):
    m = Mol(); P = m.atom("Pd", (0, 0)); m.ligate(P, "sqpl", N="L", S="L", E=X, W=Y); return m
p0 = Mol(); P = p0.atom("Pd^0", (0, 0)); p0.ligate(P, "lin", E="L", W="L")
c = S.cycle([(p0, 90), (pd("Br", "Ar"), -30), (pd("Ar'", "Ar"), 210)], radius=4)
S.branch(c, 0, "Ar–Br", kind="in"); S.branch(c, 1, "Ar'B(OH)_2 + base", kind="in"); S.branch(c, 2, "Ar–Ar'", kind="out")
print(S.save("suzuki.png")); S.check()
