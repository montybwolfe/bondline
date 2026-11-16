# P1-1 short arrows keep a real shaft, P1-2 every arrow is curved, P1-7 a cyclic array turns one way; negative tests
from bondline import *
st = Style()
def metrics(S):
    for i, c in enumerate(S.curlies):
        L = float(np.hypot(*np.diff(c._arc, axis=0).T).sum()) / st.bond
        print(f"  #{i} {c._class:5s} sweep {c._sweep:3.0f}  shaft {L:.2f} BU  {'ccw' if c._ccw else 'cw '}  ok={c._ok}")
print("== E2: base, C-H into C-C (adj), C-Br onto Br (hook)")
S = Scene(st); m = Mol(); c1 = m.atom("C", (0, 0)); c2 = m.attach(c1, 30); c3 = m.attach(c2, -30); Br = m.attach(c2, 90, "Br")
H = m.attach(c3, -90, "H"); B = m.atom("O", (0, 0)); B.charge = -1; m.attach(B, -90, "Et"); m.dock(B, H, -90)
S.curly(B.lp(), H); S.curly((H, c3), (c3, c2)); S.curly((c2, Br), Br); S.row([m]); S.save("e2.png"); S.check(); metrics(S)
print("== enolate attacking ethanal: lone pair into C-O (lpadj)")
S = Scene(st); m = Mol(); Ca = m.atom("C", (0, 0)); Cc = m.attach(Ca, 30); O = m.attach(Cc, -30, "O"); m.find(Ca, Cc).order = 2; O.charge = -1
E1 = m.atom("C", (0, 0)); m.attach(E1, 210); EO = m.attach(E1, 90, "O", 2); m.dock(E1, Ca, 180, 1.8)
S.curly(O.lp(), (O, Cc)); S.curly((Ca, Cc), E1); S.curly((E1, EO), EO); S.row([m]); S.save("enolate.png"); S.check(); metrics(S)
print("== radical abstraction: hook onto an unlabelled carbon stays compact and never wraps round a bond")
S = Scene(st); m = Mol(); C = m.atom("C", (0, 0), name="C"); [m.attach(C, t, "Me") for t in (90, 210)]; H = m.attach(C, -30, "H", name="H")
Br = m.atom("Br", (0, 0), name="Br"); Br.radical = True; m.dock(Br, H, -30)
S.fish(Br, Mid(H, Br)); S.fish((C, H), Mid(H, Br)); S.fish((C, H), C); S.row([m]); S.save("radical.png"); S.check(); metrics(S)
print("== Diels-Alder: one cyclic array, all inside the ring, all turning the same way")
S = Scene(st); m = Mol(); r = m.ring(6, rot=30)
m.find(r[0], r[1]).order = 2; m.find(r[2], r[3]).order = 2; m.find(r[4], r[5]).order = 2; m.unbond(r[3], r[4]); m.unbond(r[5], r[0])
S.curly((r[0], r[1]), (r[1], r[2])); S.curly((r[2], r[3]), Mid(r[3], r[4])); S.curly((r[4], r[5]), Mid(r[5], r[0]))
p = m.copy()
for i, o in enumerate([1, 2, 1, 1, 1, 1]): p.bond(p[r[i]], p[r[(i + 1) % 6]], o)
S.validate(m, p); S.row([m, S.arrow(above="Δ"), p]); S.save("da.png"); S.check(); metrics(S)
print("== negative tests: the checks fire on tampered arrows")
S._build(); c = S.curlies[0]
keep = (c._arc, c._sweep); c._arc = np.vstack([c._arc[0], c._arc[-1]]); c._sweep = 8.0
print("  straight:", [x for x in S._arrow_issues() if "straight" in x])
c._arc, c._sweep = keep; S.curlies[1]._ccw = not S.curlies[1]._ccw
print("  mixed loop:", [x for x in S._arrow_issues() if "cyclic" in x]); S.curlies[1]._ccw = not S.curlies[1]._ccw
S = Scene(st); m = Mol(); c0 = m.atom("C", (0, 0)); X = m.attach(c0, 0, "Br"); m.attach(c0, 120); m.attach(c0, 240)
S.curly((c0, X), X); S.row([m]); S._build(); c = S.curlies[0]; P = S.P(X); u = unit(S.P(c0) - P)
c._arc = np.array([P + perp(u) * 9 + u * 2, P + u * 4.5, P - perp(u) * 9 + u * 2])     # slips between Br and its bond
print("  crossing:", [x for x in S._arrow_issues() if "crosses a bond" in x])
