# Charges, radicals, lone pairs, partial charges, numbering
from bondline import *
st = Style()
print("== ions")
S = Scene(st)
cat = Mol(); c = cat.atom("C", (0, 0)); cat.attach(c, 90, "Me"); cat.attach(c, 210, "Me"); cat.attach(c, -30, "Me"); c.charge = 1
an = Mol(); o = an.atom("O", (0, 0)); an.attach(o, 210); o.charge = -1; o.add_lps(3)
am = Mol(); n = am.atom("N", (0, 0)); [am.attach(n, t, "Me") for t in (90, 210, 330)]; am.attach(n, 0, "H"); n.charge = 1
nit = Mol(); c = nit.atom("C", (0, 0)); nn = nit.attach(c, 30, "N"); nn.charge = 1; o1 = nit.attach(nn, 90, "O", 2); o2 = nit.attach(nn, -30, "O"); o2.charge = -1
S.row([cat, an, am, nit]); print(S.save("ions.png")); S.check()
print("== radicals and marks")
S = Scene(st)
rd = Mol(); c = rd.atom("C", (0, 0)); rd.attach(c, 210, "Me"); rd.attach(c, -30, "Me"); c.radical = True
pc = Mol(); c = pc.atom("C", (0, 0)); o = pc.attach(c, 30, "O", 2); pc.attach(c, 150, "Me"); pc.attach(c, -90, "H"); c.partial = "+"; o.partial = "-"
nb = Mol(); r = nb.ring(6, rot=90)
for i, a in enumerate(r): a.number = i + 1
S.orbital(rd.atoms[0], "p", 90, electrons=1)
S.row([rd, pc, nb]); print(S.save("marks.png")); S.check()
print("== plain charges and hetero colours")
st2 = Style(hetero=True, charge="plain"); S = Scene(st2)
ac = Mol(); c = ac.atom("C", (0, 0)); o1 = ac.attach(c, 90, "O", 2); o2 = ac.attach(c, -30, "O"); o2.charge = -1; ac.attach(c, 210, "Me")
su = Mol(); s = su.atom("S", (0, 0)); su.attach(s, 90, "O", 2); su.attach(s, -90, "O", 2); su.attach(s, 0, "Cl"); su.attach(s, 180, "Ph")
S.row([ac, su]); print(S.save("plain.png")); S.check()
