# Style presets and bond styles
from bondline import *
print("== acs preset")
st = Style(preset="acs"); S = Scene(st)
m = from_smiles("CC(=O)N[C@@H](Cc1ccccc1)C(=O)OC", align=(1, 2, 30)); stereo_report(m)
S.row([m]); print(S.save("acs.png")); S.check()
print("== bond styles")
st = Style(); S = Scene(st)
m = Mol(); a = m.atom("C", (0, 0)); styles = ["plain", "bold", "dash", "dot", "wavy", "dative", "partial"]
prev = a
for i, s in enumerate(styles):
    nx = m.atom("C", (i * 1.6 + 1.6, 0)); 
x = Mol(); a0 = x.atom("C", (0, 0))
row = []
for s in styles:
    q = Mol(); u = q.atom("C", (0, 0)); v = q.attach(u, 30); w = q.attach(v, -30, style=s); row.append(q)
S.row(row); print(S.save("bondstyles.png")); S.check()
