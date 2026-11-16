# P0-5: heads on curved black arrows (cycle arcs, in/out branches) sit on the curve with nothing of the shaft past them.
from bondline import *
import numpy as np
st = Style()
def pd(X, Y):
    m = Mol(); P = m.atom("Pd", (0, 0)); m.ligate(P, "sqpl", N="L", S="L", E=X, W=Y); return m
def cycle(S):
    p0 = Mol(); P = p0.atom("Pd^0", (0, 0)); p0.ligate(P, "lin", E="L", W="L")
    c = S.cycle([(p0, 90), (pd("Br", "Ar"), -30), (pd("Ar'", "Ar"), 210)], radius=4)
    S.branch(c, 0, "Ar–Br", kind="in"); S.branch(c, 1, "Ar'B(OH)_2 + base", kind="in"); S.branch(c, 2, "Ar–Ar'", kind="out")
    S.branch(c, 1, "X^−", kind="out", at=0.3)
print("== cycle with in and out branches, 1200 dpi")
S = Scene(st); cycle(S); print(S.save("cycle1200.png", dpi=1200)); S.check()
print("== negative: the old construction (shaft running up to the tip, head on a two-sample chord) is caught")
def old_curve(self, pts, scale=1.0):
    st = self.st; P = np.asarray(pts, float); tip = P[-1]; d = unit(P[-1] - P[-3])
    shaft = list(P[:-1]) + [tip - d * 0.55 * st.head_len]
    self._arrow_log.append({"shaft": np.array(shaft), "tip": tip, "d": d, "h": st.head_len, "w": st.head_w, "k": st.head_notch,
                            "lw": st.lw, "curved": True, "half": 0})
    return [("line", shaft, st.lw, st.ink, None), self._head(tip, d)]
keep = Scene._curve_arrow; Scene._curve_arrow = old_curve
S = Scene(st); cycle(S); S.check()
Scene._curve_arrow = keep
