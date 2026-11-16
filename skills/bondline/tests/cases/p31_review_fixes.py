# Fixes from a review of the teaching-figure renders: a site dot on the bare end of an offset double bond
# (a terminal CH2=) caps both lines instead of hanging off the main one; curved callouts leave the end of their note
# level and sweep round to the target in an even curve (no hook); a site label of two lines ("E" over "more δ+") is
# placed as its own box, aligned on the side facing its dot; a one-bond label (Me) on an atom with more bonds is refused.
from bondline import *
BLUE, GREEN, MAG = "#1F4FD1", "#2E7D32", "#D81B9A"

S = Scene(Style())
en = Mol(); c1 = en.atom("C", (0, 0)); en.attach(c1, 90, "O", 2); en.attach(c1, -30, "R'")
c2 = en.attach(c1, 210); c3 = en.attach(c2, 150, order=2)
S.dot(c3, BLUE, "E"); S.dot(c1, BLUE, "E\nmore δ+")
S.pointer(c3, "1,4-addition\n<i>thermodynamic</i>", side=120, bend=-70, color=GREEN, length=1.2)
S.pointer(c1, "1,2-addition\n<i>kinetic</i>", side=-110, bend=70, color=MAG, length=1.2)
S.row([en]); print(S.save("callouts.png", dpi=300)); S.check()
b = en.find(c2, c3); off = S._site_xy(c3) - S.P(c3)
print("terminal dot: off the vertex by", round(norm(off) / (S.st.dbl / 2), 3), "x half the line spacing, at",
      round(abs(adiff(ang(off), ang(S.P(c3) - S.P(c2))))), "deg to the bond")
print("ring/chain vertex dot stays on the vertex:", norm(S._site_xy(c1) - S.P(c1)) < 1e-9)
for p in S.pointers:
    pr = [q for q in S._ptr_geo(p)["prims"] if q[0] == "line"][0][1]; P = np.asarray(pr, float)
    t = [ang(P[i + 1] - P[i]) for i in range(len(P) - 1)]; turn = [abs(adiff(t[i + 1], t[i])) for i in range(len(t) - 1)]
    print(f"callout: leaves level {abs(adiff(t[0], 0)) < 1 or abs(adiff(t[0], 180)) < 1}, sharpest turn per step "
          f"{max(turn):.1f} deg (even: < 8)")

try: from_smiles("C[C@]12CCC3C(CCC4CCCC[C@]34C)C1CC[C@@H]2[*]", labels={0: "Me", 15: "Me"})
except ValueError as e: print("raised:", e)
m = Mol(); c = m.atom("C", (0, 0)); me = m.attach(c, 0, "Me"); m.attach(me, 60, "C")
S = Scene(Style()); S.row([m]); S.check()
