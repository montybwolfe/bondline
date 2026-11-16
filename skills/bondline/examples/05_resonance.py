# Resonance: contributors in one bracket joined by double-headed arrows (kind="<->"), each built from the last by
# moving electrons exactly as its curly arrows say (validate() replays them). Enolate of acetone.
from bondline import *
S = Scene(Style())
a = Mol(); c1 = a.atom("C", (0, 0), name="ca"); c2 = a.attach(c1, 30, name="c"); a.attach(c2, 90, "O", 2, name="o")
a.attach(c2, -30); c1.charge = -1
S.curly(c1.lp(), (c1, c2)); S.curly((c2, a["o"]), a["o"])
b = a.copy(); b.find(b["ca"], b["c"]).order = 2; b.find(b["c"], b["o"]).order = 1; b["ca"].charge = 0; b["o"].charge = -1
S.validate(a, b); S.brackets([a, b])
S.row([a, S.arrow(kind="<->"), b]); print(S.save("enolate_resonance.png")); S.check()
