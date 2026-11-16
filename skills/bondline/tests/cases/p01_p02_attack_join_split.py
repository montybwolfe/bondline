# P0-1: an arrow to Mid(a, b) across a wide gap points at the electrophile atom; a concerted pair still meets nose to
# nose; a head that ends in empty space is reported. P0-2: join() lays the joined fragment out on standard angles and
# split()/release() send a departing group off as its own item; check() reports leaving groups left in place and strays.
from bondline import *
st = Style()
def enolate_plus_aldehyde(dist):
    m = Mol(); Ca = m.atom("C", (0, 0), name="Ca"); Cc = m.attach(Ca, 30, name="Cc"); O = m.attach(Cc, -30, "O", name="O")
    m.find(Ca, Cc).order = 2; O.charge = -1
    E = m.atom("C", (0, 0), name="E"); m.attach(E, 210, name="Em"); EO = m.attach(E, 90, "O", 2, name="EO"); m.dock(E, Ca, 180, dist)
    return m
print("== Mid across 2.2 BU: aimed at the carbonyl carbon")
S = Scene(st); m = enolate_plus_aldehyde(2.2)
S.curly(m["O"].lp(), (m["O"], m["Cc"])); S.curly((m["Ca"], m["Cc"]), Mid(m["Ca"], m["E"])); S.curly((m["E"], m["EO"]), m["EO"])
p = m.copy(); p.find(p["Ca"], p["Cc"]).order = 1; p.find(p["Cc"], p["O"]).order = 2; p["O"].charge = 0
p.join(p["Ca"], p["E"]); p.find(p["E"], p["EO"]).order = 1; p["EO"].charge = -1
S.validate(m, p); Scene.continuity(m, p)
S.row([m, S.arrow(), p]); print(S.save("attack.png")); S.check()
print("== concerted pair (radical recombination) still meets nose to nose")
S = Scene(st); r = Mol(); a = r.atom("C", (0, 0), name="a"); r.attach(a, 120, "Me"); r.attach(a, 240, "Me"); a.radical = True
b = r.atom("C", (2.0, 0), name="b"); r.attach(b, 60, "Me"); r.attach(b, -60, "Me"); b.radical = True
S.fish(a, Mid(a, b)); S.fish(b, Mid(a, b))
q = r.copy(); q["a"].radical = False; q["b"].radical = False; q.join(q["a"], q["b"])
S.validate(r, q); S.row([r, S.arrow(), q]); print(S.save("recombination.png")); S.check()
print("== SN2 with join and release: product on standard angles, bromide as its own item")
S = Scene(st)
m = Mol(); C = m.atom("C", (0, 0), name="C"); Br = m.attach(C, -30, "Br", name="Br"); m.attach(C, 90, "Me")
O = m.atom("HO", (0, 0), name="O"); O.charge = -1; m.dock(O, C, 210, 2.0)
S.curly(O.lp(), C); S.curly((C, Br), Br)
p = m.copy(); p["O"].charge = 0; br = p.release(p["Br"], p["C"]); br["Br"].charge = -1; p.join(p["C"], p["O"])
S.validate(m, [p, br]); Scene.continuity(m, [p, br])
S.row([m, S.arrow(), p, "+", br]); print(S.save("sn2.png")); S.check()
print("== SN1 with release: the carbocation left as a T is redrawn as a Y (one bond kept, the others turned least)")
S = Scene(st)
m = Mol(); C = m.atom("C", (0, 0), name="C"); [m.attach(C, t, "Me") for t in (90, 180, 270)]; Br = m.attach(C, 0, "Br", name="Br")
S.curly((C, Br), Br)
p = m.copy(); br = p.release(p["Br"], p["C"]); br["Br"].charge = -1; p["C"].charge = 1
th = sorted(round(ang(n.xy - p["C"].xy) % 360) for n in p.nbrs(p["C"])); print("cation bonds at", th); assert th == [60, 180, 300], th
S.validate(m, [p, br]); Scene.continuity(m, [p, br])
S.row([m, S.arrow(), p, "+", br]); print(S.save("sn1.png")); S.check()
print("== negative: a leaving group left where it was, a stray by-product, and a head pointing at nothing")
S = Scene(st)
m = Mol(); C = m.atom("C", (0, 0)); [m.attach(C, t, "Me") for t in (90, 180, 270)]; Br = m.attach(C, 0, "Br"); m.unbond(C, Br)
Br.charge = -1; C.charge = 1
w = Mol(); c = w.atom("C", (0, 0)); w.attach(c, 30); o = w.atom("O", (5, 0)); w.attach(o, 0, "H"); w.attach(o, 120, "H")
z = Mol(); zc = z.atom("C", (0, 0)); zo = z.attach(zc, 90, "O", 2); S.curly(zo.lp(), (3.0, 3.0))
S.row([m, w, z]); S.check()
