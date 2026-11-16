import sys, os, traceback, hashlib, time
from bondline import *
import time
# --- scheme
st = Style(); S = Scene(st)
def ethanal():
    m = Mol(); c1 = m.atom("C",(0,0)); c2 = m.attach(c1, 30); m.attach(c2, -30, "O", 2); return m
a, b = ethanal(), ethanal()
p = Mol(); c = p.atom("C",(0,0)); c2=p.attach(c,30); p.attach(c2,90,"OH"); c3=p.attach(c2,-30); c4=p.attach(c3,30); p.attach(c4,-30,"O",2)
S.caption(a,"ethanal"); S.caption(b,"ethanal"); S.caption(p,"3-hydroxybutanal")
S.row([a,"+",b,S.arrow(above="NaOH (aq)", below="5 °C"),p]); print(S.save("aldol.png", dpi=300)); S.check()

# --- mechanism, built as the skill instructs
t=time.time()
S = Scene(st)
s1 = Mol(); Ca = s1.atom("C",(0,0),name="Ca"); Cc = s1.attach(Ca,30,name="Cc"); O = s1.attach(Cc,-30,"O",2,name="O")
H = s1.attach(Ca,-90,"H",name="H")
Ob = s1.atom("O",(0,0),name="Ob"); Ob.charge=-1; s1.attach(Ob,180,"H",name="Hb"); s1.dock(Ob,H,-90)
S.curly(Ob.lp(),H); S.curly((H,Ca),(Ca,Cc)); S.curly((Cc,O),O)
# step 2: enolate + second ethanal
s2 = s1.copy(); s2.remove_frag(s2["Ob"]); s2.remove(s2["H"]); s2.find(s2["Ca"],s2["Cc"]).order=2; s2.find(s2["Cc"],s2["O"]).order=1; s2["O"].charge=-1
e = Mol(); E1 = e.atom("C",(0,0),name="Ec"); E2 = e.attach(E1,150,name="Em"); EO = e.attach(E1,90,"O",2,name="EO")
for x in e.atoms: s2.atoms.append(x); x.mol=s2
s2.bonds += e.bonds
s2.dock(s2["Ec"], s2["Ca"], 180)
S.curly(s2["O"].lp(), (s2["O"],s2["Cc"])); S.curly((s2["Ca"],s2["Cc"]), Mid(s2["Ca"],s2["Ec"])); S.curly((s2["Ec"],s2["EO"]), s2["EO"])
s3 = s2.copy(); s3.find(s3["Ca"],s3["Cc"]).order=1; s3.find(s3["Cc"],s3["O"]).order=2; s3["O"].charge=0
s3.bond(s3["Ca"],s3["Ec"]); s3.find(s3["Ec"],s3["EO"]).order=1; s3["EO"].charge=-1
S.validate(s1,s2); S.validate(s2,s3)
Scene.continuity(s1,s2,s3)
S.flow([s1,S.arrow(below="−H_2O",step=1),s2,S.arrow(step=2),s3], width_cm=18)
print(S.save("mech_naive.png", dpi=300)); S.check(); print("timing:", "ok" if time.time() - t < 120 else "slow")
