# Chair and Newman templates
from bondline import *
st = Style()
print("== chair")
S = Scene(st)
m = Mol(); at, ax, eq = m.chair()
m.attach(at[0], eq[0], "Me"); m.attach(at[0], ax[0], "H"); m.attach(at[3], eq[3], "tBu")
S.row([m]); print(S.save("chair.png")); S.check()
print("== newman")
S = Scene(st)
n = Mol(); n.newman(["Me", "H", "H"], ["Me", "H", "H"], back_rot=90)
n2 = Mol(); n2.newman(["Me", "H", "H"], ["H", "Me", "H"])
S.row([n, S.arrow(kind="<=>"), n2]); print(S.save("newman.png")); S.check()
