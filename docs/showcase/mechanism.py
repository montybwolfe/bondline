# Robinson annulation to the Wieland-Miescher ketone, base-mediated (hydroxide/water), one Mol per step:
# Michael addition of the 2-methylcyclohexane-1,3-dione enolate to methyl vinyl ketone, proton transfer, enolate
# formation at the methyl ketone, intramolecular aldol addition, proton transfer, and E1cB dehydration.
# The bicyclic skeleton is built once; the open-chain steps are copies with the new C-C bond broken and the chain swung.
from bondline import *
S = Scene(Style())

G = Mol()
P = {"q": (0, 0.5), "k2": (0, -0.5), "k1": (-0.866, 1), "r6": (-1.732, 0.5), "r5": (-1.732, -0.5), "r4": (-0.866, -1),
     "b": (0.866, 1), "a": (1.732, 0.5), "c": (1.732, -0.5), "m": (0.866, -1),
     "o1": (-0.866, 2), "me": (0, 1.5), "o3": (2.598, -1.0), "o2": (0, -1.5)}
A = {k: G.atom("O" if k[0] == "o" else "C", xy, name=k) for k, xy in P.items()}
for u, v in (("q", "k1"), ("k1", "r6"), ("r6", "r5"), ("r5", "r4"), ("r4", "k2"), ("k2", "q"), ("q", "b"), ("b", "a"),
             ("a", "c"), ("c", "m"), ("m", "k2"), ("k1", "o1"), ("q", "me"), ("c", "o3"), ("k2", "o2")):
    G.bond(A[u], A[v])

def order(m, *pairs):
    for u, v, o in pairs: m.find(m[u], m[v]).order = o

def opened(m):                     # break the aldol bond and swing the side chain clear of the ring ketone
    m.unbond(m["m"], m["k2"]); m.swing(m["a"], m["b"], 10); m.swing(m["c"], m["a"], 10); m.swing(m["m"], m["c"], 10)
    m["o2"].xy = m["k2"].xy + dirv(300); m["o2"].deliberate = True   # ring C=O tilted 30 deg to make room for the chain
    return m

def water(m, H_at, ang, name):     # H2O with its H 2 BU from the acceptor
    h = m.atom("H", H_at, name="H" + name); o = m.attach(h, ang, "O", name="O" + name); m.attach(o, None, "H"); return h, o

# s2: Michael adduct (enolate) + water at the alpha carbon
s2 = opened(G.copy()); order(s2, ("k1", "o1", 2), ("a", "c", 2), ("k2", "o2", 2)); s2["o3"].charge = -1
water(s2, s2["a"].xy + dirv(40) * 2.0, 40, "w")
S.curly(s2["o3"].lp(), (s2["o3"], s2["c"])); S.curly((s2["a"], s2["c"]), s2["Hw"]); S.curly((s2["Hw"], s2["Ow"]), s2["Ow"])

# s1: dione enolate + methyl vinyl ketone (the adduct with the new bond broken, MVK moved off)
s1 = s2.copy(); s1.remove_frag(s1["Ow"]); s1.unbond(s1["q"], s1["b"])
order(s1, ("k1", "o1", 1), ("k1", "q", 2), ("b", "a", 2), ("a", "c", 1), ("c", "o3", 2))
s1["o1"].charge = -1; s1["o3"].charge = 0; s1["me"].xy = s1["q"].xy + dirv(30)   # sp2 enolate carbon: Me on the bisector
s1.dock(s1["b"], s1["q"], 330, 2.0)
S.curly(s1["o1"].lp(), (s1["o1"], s1["k1"])); S.curly((s1["k1"], s1["q"]), s1["b"])
S.curly((s1["b"], s1["a"]), (s1["a"], s1["c"])); S.curly((s1["c"], s1["o3"]), s1["o3"])

# s3: triketone + hydroxide at the methyl ketone
s3 = s2.copy(); s3.transfer(s3["Hw"], s3["a"], 40); s3.remove_frag(s3["Ow"])
order(s3, ("a", "c", 1), ("c", "o3", 2)); s3["o3"].charge = 0
hm1 = s3.attach(s3["m"], 300, "H", name="hm1")
ob = s3.atom("O", hm1.xy + dirv(300) * 2.0, name="Ob"); ob.charge = -1; s3.attach(ob, None, "H")
S.curly(ob.lp(), hm1); S.curly((s3["m"], hm1), (s3["m"], s3["c"])); S.curly((s3["c"], s3["o3"]), s3["o3"])

# s4: enolate attacks the ring ketone (intramolecular aldol)
s4 = s3.copy(); s4.transfer(s4["hm1"], s4["Ob"]); s4.remove_frag(s4["Ob"])
order(s4, ("c", "m", 2), ("c", "o3", 1)); s4["o3"].charge = -1; s4.tidy(carbon=True)
S.curly(s4["o3"].lp(), (s4["o3"], s4["c"])); S.curly((s4["c"], s4["m"]), s4["k2"]); S.curly((s4["k2"], s4["o2"]), s4["o2"])

# s5: bicyclic alkoxide + water
s5 = G.copy(); order(s5, ("k1", "o1", 2), ("c", "o3", 2)); s5["o2"].charge = -1
water(s5, s5["o2"].xy + dirv(210) * 2.0, 210, "x")
S.curly(s5["o2"].lp(), s5["Hx"]); S.curly((s5["Hx"], s5["Ox"]), s5["Ox"])

# s6: aldol + hydroxide at the CH2 between C=O and C-OH
s6 = s5.copy(); s6.transfer(s6["Hx"], s6["o2"], 210); s6["o2"].charge = 0; s6["Ox"].charge = -1
hm2 = s6.attach(s6["m"], None, "H", name="hm2"); s6.dock(s6["Ox"], hm2, 300)
S.curly(s6["Ox"].lp(), hm2); S.curly((s6["m"], hm2), (s6["m"], s6["c"])); S.curly((s6["c"], s6["o3"]), s6["o3"])

# s7: enolate expels hydroxide (E1cB)
s7 = s6.copy(); s7.transfer(s7["hm2"], s7["Ox"]); s7.remove_frag(s7["Ox"])
order(s7, ("c", "m", 2), ("c", "o3", 1)); s7["o3"].charge = -1
S.curly(s7["o3"].lp(), (s7["o3"], s7["c"])); S.curly((s7["c"], s7["m"]), (s7["m"], s7["k2"])); S.curly((s7["k2"], s7["o2"]), s7["o2"])

# s8: Wieland-Miescher ketone + hydroxide
s8 = s7.copy(); order(s8, ("c", "o3", 2), ("c", "m", 1), ("m", "k2", 2)); s8["o3"].charge = 0
s8.unbond(s8["k2"], s8["o2"]); s8["o2"].charge = -1; ho = s8.split(s8["o2"]); s8.tidy(); ho.tidy()

steps = [s1, s2, s3, s4, s5, s6, s7, [s8, ho]]
for x, y in zip(steps, steps[1:]): S.validate(x, y)
Scene.continuity(*steps)
S.caption(s8, "Wieland–Miescher ketone")
S.flow([s1, S.arrow(below="Michael"), s2, S.arrow(below="H^+ transfer"), s3, S.arrow(below="enolisation"), s4,
        S.arrow(below="aldol"), s5, S.arrow(below="H^+ transfer"), s6, S.arrow(below="enolisation"), s7,
        S.arrow(below="E1cB"), s8, "+", ho], width_cm=16, row_gap=1.3)
print(S.save("mechanism.png", dpi=260)); S.check()
