import sys, os, traceback, hashlib, time
from bondline import *
st = Style()
def run(name, f):
    print("=====", name)
    try:
        S = Scene(st); f(S); print(S.save(name + ".png", dpi=300)); S.check()
    except Exception: traceback.print_exc(limit=3)

def sn2(S):
    m = Mol(); C = m.atom("C",(0,0),name="C"); Br = m.attach(C,0,"Br",name="Br"); m.attach(C,120,"H"); m.attach(C,240,"H"); m.attach(C,180,"H",name="Hx")
    m.remove(m["Hx"]); m.attach(C,90,"H") ; 
    O = m.atom("HO",(0,0),name="O"); O.charge=-1; m.dock(O,C,180)
    S.curly(O.lp(),C); S.curly((C,Br),Br)
    p = m.copy(); p.bond(p["O"],p["C"]); p["O"].charge=0; p.unbond(p["C"],p["Br"]); p["Br"].charge=-1
    S.validate(m,p); Scene.continuity(m,p)
    S.row([m,S.arrow(above="S_N2"),p])
def e2(S):
    m = Mol(); c1 = m.atom("C",(0,0),name="c1"); c2 = m.attach(c1,30,name="c2"); c3=m.attach(c2,-30,name="c3"); Br=m.attach(c2,90,"Br",name="Br")
    H = m.attach(c3,-90,"H",name="H"); B = m.atom("EtO",(0,0),name="B"); B.charge=-1; m.dock(B,H,-90)
    S.curly(B.lp(),H); S.curly((H,c3),(c3,c2)); S.curly((c2,Br),Br)
    p = m.copy(); p.transfer(p["H"],p["B"],90); p["B"].charge=0; p.find(p["c2"],p["c3"]).order=2; p.unbond(p["c2"],p["Br"]); p["Br"].charge=-1
    S.validate(m,p); S.row([m,S.arrow(),p])
def eas(S):
    m = Mol(); r = m.ring(6,rot=30,kekule=True)
    E = m.atom("Br",(0,0),name="E"); E.charge=1; m.dock(E,r[0],30)
    S.curly((r[0],r[1]),E)
    a = m.copy(); a.find(a[r[0]],a[r[1]]).order=1; a.bond(a[r[0]],a["E"]); a["E"].charge=0; a[r[1]].charge=1
    a["E"].xy = a[r[0]].xy + dirv(60); H = a.attach(a[r[0]],0,"H",name="H")
    S.validate(m,a)
    b = a.copy(); b.find(b[r[2]],b[r[3]]).order=1; b.find(b[r[1]],b[r[2]]).order=2; b[r[1]].charge=0; b[r[3]].charge=1
    S.curly((a[r[2]],a[r[3]]),(a[r[1]],a[r[2]]))
    S.validate(a,b)
    S.brackets([a,b])
    S.row([m,S.arrow(),a,S.arrow(kind="<->"),b])
def carbonyl_add(S):
    m = Mol(); c = m.atom("C",(0,0),name="c"); m.attach(c,150,"Ph"); m.attach(c,30,"Me"); O=m.attach(c,-90,"O",2,name="O")
    Nu = m.atom("H_3C",(0,0),name="Nu"); Mg = m.attach(Nu,180,"MgBr",name="Mg"); m.dock(Nu,c,90)
    S.curly((Nu,Mg),c); S.curly((c,O),O)
    S.row([m])
def da(S):
    m = Mol(); r = m.ring(6, rot=30)
    # diene r0..r3 ; dienophile r4,r5; break r3-r4 and r5-r0
    m.find(r[0],r[1]).order=2; m.find(r[2],r[3]).order=2; m.find(r[4],r[5]).order=2
    m.unbond(r[3],r[4]); m.unbond(r[5],r[0])
    S.curly((r[0],r[1]),(r[1],r[2])); S.curly((r[2],r[3]),Mid(r[3],r[4])); S.curly((r[4],r[5]),Mid(r[5],r[0]))
    p = m.copy(); 
    for i,o in enumerate([1,2,1,1,1,1]): p.bond(p[r[i]],p[r[(i+1)%6]],o)
    S.validate(m,p); S.row([m,S.arrow(above="Δ"),p])
def radical(S):
    m = Mol(); C = m.atom("H_3C",(0,0),name="C"); H = m.attach(C,0,"H",name="H"); Cl = m.atom("Cl",(0,0),name="Cl"); Cl.radical=True; m.dock(Cl,H,0)
    S.fish(Cl,Mid(H,Cl)); S.fish((C,H),Mid(H,Cl)); S.fish((C,H),C)
    p = m.copy(); p.unbond(p["C"],p["H"]); p.bond(p["H"],p["Cl"]); p["Cl"].radical=False; p["C"].radical=True
    S.validate(m,p); S.row([m,S.arrow(),p])
def smiles(S):
    a = from_smiles("CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O"); b = from_smiles("C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C")
    S.caption(a,"ibuprofen"); S.caption(b,"cholesterol"); S.row([a,b])
def chairs(S):
    m = Mol(); at, ax, eq = m.chair(); m.attach(at[0], ax[0], "OH"); m.attach(at[3], eq[3], "tBu"); m.attach(at[1], ax[1], "H"); m.attach(at[1], eq[1], "Me")
    n = Mol(); n.newman(["H","H","Me"],["Me","H","H"])
    S.row([m, S.arrow(kind="<=>"), n])
def ester(S):
    m = Mol(); c = m.atom("C",(0,0),name="c"); m.attach(c,210,"R"); O = m.attach(c,90,"O",2,name="O"); m.attach(c,-30,"OH",name="OH")
    H = m.atom("H",(0,0),name="H"); A = m.attach(H,90,"OH2",name="A"); A.charge=1; m.dock(H,O,90)
    S.curly(O.lp(),H); S.curly((H,A),A)
    s2 = m.copy(); s2.transfer(s2["H"],s2["O"],90); s2["O"].charge=1; s2.remove_frag(s2["A"])
    Me = s2.atom("O",(0,0),name="Om"); s2.attach(Me,0,"Me"); s2.attach(Me,90,"H"); s2.dock(Me,s2["c"],-110)
    S.curly(s2["Om"].lp(), s2["c"]); S.curly((s2["c"],s2["O"]), s2["O"])
    S.validate(m,s2)
    S.row([m,S.arrow(kind="<=>"),s2])
def cyc(S):
    def pd(X, Y):
        m = Mol(); P = m.atom("Pd",(0,0)); m.ligate(P,"sqpl",N="L",S="L",E=X,W=Y); return m
    p0 = Mol(); P = p0.atom("Pd^0",(0,0)); p0.ligate(P,"lin",E="L",W="L")
    c = S.cycle([(p0,90),(pd("Br","Ar"),-30),(pd("Ar'","Ar"),210)], radius=4)
    S.branch(c,0,"Ar–Br",kind="in"); S.branch(c,1,"Ar'B(OH)_2 + base",kind="in"); S.branch(c,2,"Ar–Ar'",kind="out")
for n,f in [("sn2",sn2),("e2",e2),("eas",eas),("add",carbonyl_add),("da",da),("radical",radical),("smiles",smiles),("chairs",chairs),("ester",ester),("cyc",cyc)]: run(n,f)
