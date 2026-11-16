# P1-3: standard docking distances, deprotonate / protonate / attack choose the spot by rule, stock reagents are always
# drawn the same way; P1-5: tidy() gives the standard depiction and check() notes explicit X-H that nothing uses.
# Consistency: the same reaction with different substituents gives the same placement and the same arrows.
from bondline import *
st = Style()
def arrows(S):
    for i, c in enumerate(S.curlies):
        L = float(np.hypot(*np.diff(c._arc, axis=0).T).sum()) / st.bond
        print(f"  #{i} {c._class:5s} sweep {c._sweep:3.0f}  shaft {L:.2f} BU")
def ald(n=""):
    m = Mol(); ca = m.atom("C", (0, 0), name="Ca" + n); c1 = m.attach(ca, 30, name="C1" + n); m.attach(c1, -30, "O", 2, name="O1" + n); return m
# (the aldol written forwards with these helpers is p20a_rebuilt_aldol, with its acceptance metrics)
d = lambda a, b: norm(a.xy - b.xy)
print("== consistency: SN2 by the same rule on three substrates gives the same placement and the same arrows")
for X, R in (("Br", "H"), ("I", "Me"), ("OTs", "Et")):
    S = Scene(st); m = Mol(); c = m.atom("C", (0, 0), name="c"); lg = m.attach(c, 0, X, name="lg"); m.attach(c, 120, R); m.attach(c, 240, "H")
    nu = m.attack(halide("I") if X != "I" else hydroxide(), c); S.curly(nu.lp(), c); S.curly((c, lg), lg)
    S.row([m]); S._build(); iss = S.check(quiet=True)
    sw = [round(x._sweep / 5) * 5 for x in S.curlies]
    print(f"  {X:3s}/{R:2s}: nucleophile at {ang(nu.xy - c.xy) % 360:.0f} deg, {d(nu, c):.2f} BU; sweeps {sw}; issues {len(iss)}")
print("== consistency: deprotonation by HO-, EtO- and Et3N places each base the same way")
for base in (hydroxide(), alkoxide("Et"), amine("Et3N")):
    S = Scene(st); m = Mol(); c = m.atom("C", (0, 0)); o = m.attach(c, 90, "O", 2); m.attach(c, 210, "Me"); a = m.attach(c, -30)
    h = m.attach(a, -90, "H"); m.attach(a, 30, "Me"); B = m.deprotonate(base, h)
    S.curly(B.lp(), h); S.curly((h, a), (a, c)); S.curly((c, o), o); S.row([m]); S._build(); iss = S.check(quiet=True)
    print(f"  {B.label:5s}: {ang(B.xy - h.xy) % 360:.0f} deg from H, {d(B, h):.2f} BU; base arrow sweep {S.curlies[0]._sweep:.0f}; issues {len(iss)}")
print("== acid catalysis: hydronium protonates a ketone where the new O-H will point")
S = Scene(st); m = Mol(); c = m.atom("C", (0, 0)); o = m.attach(c, 90, "O", 2, name="O"); m.attach(c, 210, "Me"); m.attach(c, 330, "Me")
H = m.protonate(hydronium(), o); S.curly(o.lp(), H); S.curly((H, m.nbrs(H)[0]), m.nbrs(H)[0]); S.row([m]); S.save("acid.png"); S.check(); arrows(S)
print("== stock reagents")
for f in (hydroxide(), hydroxide(explicit=True), water(), water(explicit=False), hydronium(), hydronium(explicit=False), alkoxide("Me"), halide("Cl"), amine("iPr2NEt"), acid("Cl")):
    print("  ", " ".join(f"{a.label}{'+' if a.charge > 0 else '-' if a.charge < 0 else ''}" for a in f.atoms), "key", f.key.label)
print("== dock() keeps its signature; its default distance follows the rule")
m = Mol(); c = m.atom("C", (0, 0)); m.attach(c, 90, "O", 2); n = m.atom("C", (5, 5)); m.attach(n, 0, "MgBr"); m.dock(n, c, 210)
o2 = m.atom("O", (9, 9)); o2.charge = -1; m.dock(o2, c, 330); print(f"  C to C {d(n, c):.2f} BU, O- to C {d(o2, c):.2f} BU")
print("== negative: explicit O-H that nothing uses is noted; tidy() clears it")
S = Scene(st); m = Mol(); c = m.atom("C", (0, 0)); o = m.attach(c, 90, "O", name="O"); m.attach(o, 30, "H"); br = m.attach(c, -30, "Br")
m.attach(c, 210, "Me"); S.curly((c, br), br); S.row([m]); S.check()
m.tidy(); S.check()
