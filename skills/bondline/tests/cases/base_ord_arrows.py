# Reaction-arrow kinds, numbered steps, vertical arrows, text items
from bondline import *
st = Style()
print("== arrow kinds")
S = Scene(st)
def ket():
    m = Mol(); c = m.atom("C", (0, 0)); m.attach(c, 150, "Me"); m.attach(c, 30, "Me"); m.attach(c, -90, "O", 2); return m
a, b, c, d = ket(), ket(), ket(), ket()
S.row([a, S.arrow(kind="<=>", above="H_3O^+"), b, S.arrow(kind="=>"), c, S.arrow(above="NaBH_4", below="MeOH", step=1), d])
print(S.save("kinds.png")); S.check()
print("== vertical and text")
S = Scene(st)
e = ket(); f = ket()
S.place(e, (0, 0)); S.arrow(vertical=True, above="LDA", below="THF", at=(0, -1.6)); S.place(f, (0, -6))
S.text("<i>trans</i> only", at=(4.5, -3))
print(S.save("vertical.png")); S.check()
