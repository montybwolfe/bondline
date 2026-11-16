# Conformations: chair() returns the ring atoms with their axial and equatorial directions, so substituents go on the
# right bonds; chair_flip() is the other chair, where axial and equatorial swap on the same face. newman() takes the
# three front and three back substituents and the rotation of each set.
from bondline import *
S = Scene(Style())
eq = Mol(); at, ax, eqd = eq.chair(); eq.attach(at[0], eqd[0], "Me"); eq.attach(at[0], ax[0], "H")
axm = Mol(); at2, ax2, eq2 = axm.chair_flip(); axm.attach(at2[0], ax2[0], "Me"); axm.attach(at2[0], eq2[0], "H")
S.caption(eq, "equatorial Me"); S.caption(axm, "axial Me")
S.row([eq, S.arrow(kind="<=>"), axm], y=0)
anti = Mol(); anti.newman(["Me", "H", "H"], ["Me", "H", "H"], back_rot=270)       # butane, anti
gau = Mol(); gau.newman(["Me", "H", "H"], ["H", "Me", "H"], back_rot=270)         # gauche
S.caption(anti, "anti"); S.caption(gau, "gauche")
S.row([anti, S.arrow(kind="<=>"), gau], y=-5)
print(S.save("conformations.png")); S.check()
