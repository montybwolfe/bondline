import sys, os, traceback, hashlib, time
from bondline import *
st = Style()
# 1 arrow text overflow / multiline
S = Scene(st); a = Mol(); a.ring(6,kekule=True); b = a.copy()
ar = S.arrow(above="1. LiAlH_4, THF, 0 °C → rt, 12 h", below="2. H_3O^+ workup\nthen NaHCO_3")
S.row([a, ar, b]); print(S.save("arrowtext.png", dpi=300)); S.check()
g = S._item_geo(ar); print("arrow len pt", st.rxn*st.bond, "text widths", [t["box"][2]-t["box"][0] for t in g["texts"]])
# 2 caption shifts centreline
S = Scene(st); m = Mol(); c = m.atom("C",(0,0)); c2 = m.attach(c,30); m.attach(c2,-30,"O",2); n = m.copy()
S.caption(m, "ethanal"); arw = S.arrow(); S.row([m, arw, n])
bm = S._mol_box(m, curlies=False); bn = S._mol_box(n, curlies=False)
print("captioned mol centre y", (bm[1]+bm[3])/2, "arrow y", arw.y, "uncaptioned mol centre y", (bn[1]+bn[3])/2)
# 3 arrow metrics for E2 + aldol
def metrics(S):
    S._build()
    for i,c in enumerate(S.curlies):
        A = c._arc; L = float(np.hypot(*np.diff(A,axis=0).T).sum())
        print(f"  #{i} class={c._class} sweep={c._sweep:.0f} shaft={L/st.bond:.2f} BU head={st.c_head_len/st.bond:.2f} BU ok={c._ok}")
S = Scene(st); m = Mol(); c1 = m.atom("C",(0,0)); c2 = m.attach(c1,30); c3=m.attach(c2,-30); Br=m.attach(c2,90,"Br")
H = m.attach(c3,-90,"H"); B = m.atom("O",(0,0)); B.charge=-1; m.attach(B,-90,"Et"); m.dock(B,H,-90)
S.curly(B.lp(),H); S.curly((H,c3),(c3,c2)); S.curly((c2,Br),Br); S.row([m]); print("E2"); metrics(S); S.save("e2b.png",dpi=400)
# enolate attack: Mid vs atom, and dock distances
for tag, dist, useMid in (("mid22",2.2,True),("atom22",2.2,False),("atom16",1.6,False)):
    S = Scene(st); m = Mol(); Ca = m.atom("C",(0,0)); Cc = m.attach(Ca,30); O = m.attach(Cc,-30,"O"); m.find(Ca,Cc).order=2; O.charge=-1
    E1 = m.atom("C",(0,0)); m.attach(E1,150); EO = m.attach(E1,90,"O",2); m.dock(E1,Ca,180,dist)
    S.curly(O.lp(),(O,Cc)); S.curly((Ca,Cc), Mid(Ca,E1) if useMid else E1); S.curly((E1,EO),EO); S.row([m])
    print(tag, S.save(tag+".png",dpi=400)); iss = S.check(quiet=True); print("  check issues:", iss); metrics(S)
# 4 determinism
h = [hashlib.md5(open(f,"rb").read()).hexdigest() for f in ("atom22.png",)]; S.save("atom16b.png",dpi=400)
print("deterministic:", hashlib.md5(open("atom16.png","rb").read()).hexdigest()==hashlib.md5(open("atom16b.png","rb").read()).hexdigest())
print("font used:", S.F.family)
