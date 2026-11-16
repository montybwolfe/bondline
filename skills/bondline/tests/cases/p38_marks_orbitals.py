# Marks on atoms: lone pairs, charges, site dots, locants; orbital pictures; small notations (fixer "marks").
import numpy as np
from bondline import *

print("== two or more lone pairs on one atom stand apart as pairs (120 deg for two, 90 for more), on the open side")
S = Scene(Style())
k = Mol(); c = k.atom("C", (0, 0)); k.attach(c, 210); o = k.attach(c, 90, "O", 2); k.attach(c, 330)
w = water(); w.key.lp(); src = w.key.lp(); k.attack(w, c)                # water adds to a ketone C=O
S.curly(src, c)
et = Mol(); x = et.atom("O", (0, 0)); et.attach(x, 150); et.attach(x, 30); x.add_lps(2)       # ether O in a zigzag
hx = Mol(); y = hx.atom("O", (0, 0)); y.charge = -1; hx.attach(y, 0, "H"); y.add_lps(3)       # ⊖:Ö–H
al = Mol(); z = al.atom("O", (0, 0)); al.attach(z, 210); z.charge = -1; z.add_lps(3)          # alkoxide
br = Mol(); b = br.atom("Br", (0, 0)); b.charge = -1; b.add_lps(4)                             # bromide
S.row([k, et, hx, al, br], gap=1.6)
print(S.save("lone_pairs.png", dpi=300)); S.check()
for nm, a in (("water", w.key), ("ether", x), ("hydroxide", y), ("alkoxide", z), ("bromide", b)):
    print(nm, "pairs", sorted(round(l.angle) for l in a.lps), "charge", a._chg[2] if a._chg else None)

print("== check(): a lone pair or radical with no room is a slip; carbenes and cyclic transition states are not")
S = Scene(Style())
p = from_smiles("C1CCCN1C"); p[3].lp()                                  # lone pair on a ring CH2 (meant for the N)
r = from_smiles("Cc1ccccc1"); r[3].radical = True                       # radical on an arene CH
n = Mol(); nn = n.atom("N", (0, 0)); [n.attach(nn, t) for t in (90, 210, 330)]; nn.add_lps(2)   # amine N with two pairs
S.row([p, r, n], gap=1.5); S.check()
S = Scene(Style())
k = Mol(); c = k.atom("C", (0, 0), show=True); k.attach(c, 150, "H"); k.attach(c, 30, "H"); c.lp()      # :CH2
d = from_smiles("Cl[C]Cl"); d[1].radical = False; d[1].lp()                                         # singlet :CCl2
ts = from_smiles("C1=CCCCC1", align=(0, 1, 90), flip=True)                # Simmons-Smith butterfly transition state
ta, tb = (ts[0], ts[1]) if ts[0].xy[1] > ts[1].xy[1] else (ts[1], ts[0])
tc = ts.atom("C", (0, 0), show=True); ts.dock(tc, ta, -20, 1.5)
ts.bond(ta, tc, style="dash"); ts.bond(tb, tc, style="dash"); ts.find(ta, tb).style = "partial"
ti = ts.attach(tc, 40, "I", style="dash", length=1.2, free=True); tz = ts.attach(tc, -40, "Zn", style="dash", length=1.2, free=True)
ts.bond(tz, ti, style="dash"); ts.attach(tz, -60, "I"); ts.attach(tc, 100, "H"); ts.attach(tc, 260, "H")
S.brackets([ts], charge="‡"); S.row([k, d, ts], gap=1.5)
print(S.save("carbenes_ts.png", dpi=300)); S.check()
print("formula :CH2", hill(formula(k)[0]), " :CCl2", hill(formula(d)[0]))

print("== site labels: their text box beside the dot, clear of a triple bond; a dot stays with its step")
S = Scene(Style())
dn = from_smiles("N#CCCCC(=O)C", align=(2, 5, 0))
S.dot(dn[1], "#1F4FD1", "E"); S.dot(dn[4], "#E8261A", "2^{nd} Nu"); S.dot(dn[5], "#1F4FD1", "1^{st} E")
nx = dn.copy()
S.row([dn, S.arrow(), nx]); print(S.save("site_labels.png", dpi=300)); S.check()
for a in (dn[1], dn[4], dn[5]):
    bx = S._mark_box(a, "_site"); O = S._site_xy(a)
    print(a.site[1], "at", a._site[2], "deg, box", round(box_dist(O, bx), 1), "pt from its dot")
print("dots after copy():", sum(1 for a in nx.atoms if a.site))

print("== charges: the default place passes check(); a charge by its own p orbital belongs to it; a perspective frame")
S = Scene(Style())
nn = Mol(); n = nn.atom("N", (0, 0)); n.charge = 1; nn.attach(n, 90, "O", 2); nn.attach(n, 270, "O", 2)   # O=N⊕=O
cc = Mol(); c = cc.atom("C", (0, 0)); c.charge = 1; cc.attach(c, 0, "R"); cc.attach(c, 150, "R", style="hash")
cc.attach(c, 210, "R", style="wedge"); S.orbital(c, "p", 90, shade="down")                       # R3C⊕ with its empty p
cb = Mol(); k = cb.atom("C", (0, 0)); cb.attach(k, 150, "R", style="hash"); cb.attach(k, 210, "R", style="wedge")
S.orbital(k, "p", 90, shade="down"); S.orbital(k, "lobe", 0, electrons=2)                        # singlet carbene
S.row([nn, cc, cb], gap=2.0); print(S.save("charges.png", dpi=300)); S.check()
for a in (n, c):
    g, lim = S._mark_gap(a, a._chg[0], a._chg[1], a._chg[2], "chg"); print(a, "charge at", a._chg[2], "deg,", round(g, 1), "pt (allowed", round(lim, 1), ")")

print("== a ring carbanion's pair goes outside the ring; its charge beside the pair, not by a charged neighbour")
S = Scene(Style())
yl = from_smiles("[C-]1=[N+](C)C=CS1"); yl[0].lp()                                    # azolium ylide
dt = from_smiles("C[C-]1SCCCS1", align=(1, 6, 0, 3)); dt[1].lp(); S.dot(dt[1], "#E8261A", "Nu")   # dithiane anion
S.row([yl, dt], gap=2.0); print(S.save("ring_anions.png", dpi=300)); S.check()
for a in (yl[0], dt[1]): print(a, "pair", round(a.lps[0].angle) % 360, "charge", round(a._chg[2]) % 360, "in ring:", S._in_ring(a, a.lps[0].angle))

print("== locants: in a ring's corner where the outside is full; clear of their own atom's label")
S = Scene(Style())
tq = from_smiles("CC12CCCC1(C)CCC2", labels={0: "Me", 6: "Me"})          # a fused bicycle, both fusion carbons numbered
for k, a in enumerate(x for x in tq.atoms): a.number = k + 1
S.row([tq]); print(S.save("locants.png", dpi=300)); S.check()
print("locants inside a ring:", sorted(a.number for a in tq.atoms if a._num and S._in_ring(a, a._num[2])))

print("== a lone pair clear of a wide lone label (Nu), and in the empty sp2 lobe drawn for it")
S = Scene(Style())
e = Mol(); c = e.atom("C", (0, 0)); e.attach(c, 90, "R"); e.attach(c, 150, "H"); e.attach(c, 270, "H"); i = e.attach(c, 30, "I")
nu = base("Nu", charge=0); nu.dock(nu.key, c, 210, 2.4); e.merge(nu); S.curly(nu.key.lps[0], c)
cb = Mol(); k = cb.atom("C", (0, 0)); cb.attach(k, 150, "Br", style="hash"); cb.attach(k, 210, "Br", style="wedge")
S.orbital(k, "p", 90, shade="down"); S.orbital(k, "lobe", 0); lp = k.lp()
S.row([e, cb], gap=2.0); print(S.save("lp_label_lobe.png", dpi=300)); S.check()
H = np.vstack([h for g in nu.key._g for h in g["hulls"]]); d1, d2 = nu.key._lp[0][2:]
print("Nu pair clear of its letters:", round(min(float(np.hypot(*(H - d).T).min()) for d in (d1, d2)) - S.st.lp_r, 1), "pt;",
      "carbene pair at", round(lp.angle), "deg (its lobe at 0)")

print("== orbitals: s spheres show their phase (bonds stop at them), a hybrid has its back lobe where no bond runs, an orbital can carry its name")
S = Scene(Style())
ss = Mol(); s1 = ss.atom("C", (0, 0)); s2 = ss.atom("C", (1.6, 0)); ss.bond(s1, s2, style="dash")
S.orbital(s1, "s", 0, shade="up"); S.orbital(s2, "s", 0, shade="down"); S.caption(ss, "σ* (s + s)")
am = Mol(); n = am.atom("N", (0, 0)); am.attach(n, 200, "CH3", style="wedge"); am.attach(n, -20, "H", style="hash"); am.attach(n, 260, "H")
S.orbital(n, "lobe", 90, electrons=2, color="#00807F", label="n_{N(sp^3)}")   # (back lobe: on an N-H bond, so left out)
hb = Mol(); h0 = hb.atom("C", (0, 0)); hb.attach(h0, 210); hb.attach(h0, 330); h0.charge = -1   # a carbanion's sp3
S.orbital(h0, "lobe", 90, electrons=2)                                                       # pair: its back lobe shows
S.row([ss, am, hb], gap=2.0); print(S.save("orbital_shapes.png", dpi=300)); S.check()
print("lobes drawn:", [sum(1 for p in S._orbital_prims(e) if p[0] == "orb") for e in S.extras if e[0] == "orbital"],
      "fills:", [p[2] for e in S.extras if e[0] == "orbital" and e[2] == "s" for p in S._orbital_prims(e) if p[0] == "orb"])

print("== a floating substituent across a ring edge; the wiggly line of a disconnection and of an attachment point")
S = Scene(Style())
ax = from_smiles("Brc1ccccc1", labels={0: "X"}, align=(1, 0, 30))
ed = min((b for b in ax.bonds if b.order == 2), key=lambda b: b.a.xy[0] + b.b.xy[0]); r1 = ax.floating(ed.a, ed.b, "R^{1}")
es = from_smiles("CCOC(C)=O", align=(2, 3, 0)); S.wiggle((es[2], es[3]), color="#00918F")
dg = from_smiles("Cc1ccccn1", align=(1, 0, 30)); S.wiggle((dg[1], dg[0]), at=1)
S.row([ax, es, dg], gap=1.5); print(S.save("floating_wiggle.png", dpi=300)); S.check()
gaps = [pr for pr in ed._prims if pr[0] == "line"]
print("formula", hill(formula(ax)[0]), "| crossed edge drawn in", len(gaps), "lines (inner line broken)")

print("== partial charges: kept clear by their text box; charge_at steers one")
S = Scene(Style())
ee = from_smiles("C=C[*]", labels={2: "EDG"}, align=(0, 2, 0)); ee[1].partial = "+"; ee[1].charge_at = 270
S.pointer(ee[1], "<b>E</b>", side="above", length=1.0)
ts = Mol(); c = ts.atom("C", (0, 0), show=True); ts.attach(c, 90, "H"); ts.attach(c, 240, "H"); ts.attach(c, 300, "H")
br = ts.attach(c, 180, "Br", style="dash", length=1.5); io = ts.attach(c, 0, "I", style="dash", length=1.5)
br.partial = "-"; io.partial = "-"; c.partial = "+"
S.row([ee, ts], gap=2.0); print(S.save("partials.png", dpi=300)); S.check()
print("steered δ+ at", ee[1]._pc[2], "deg; the TS carbon's δ+ at", c._pc[2], "deg")

print("== a halo behind an H label (S.dot on a labelled atom); crossed dipole arrows beside polar bonds")
S = Scene(Style())
hm = Mol(); c = hm.atom("C", (0, 0)); h = hm.attach(c, 90, "H"); hm.attach(c, 210); hm.attach(c, -30); S.dot(h, color="#D9D9D9")
fb = from_smiles("Fc1ccccc1", align=(1, 0, 90)); S.dipole((fb[1], fb[0]))
S.row([hm, fb], gap=2.0); print(S.save("halo_dipole.png", dpi=300)); S.check()

print("== the + - notation of a carbene (a.ambi): one mark, net charge 0, the carbene's pair counted")
S = Scene(Style())
cb2 = Mol(); c = cb2.atom("C", (0, 0)); cb2.attach(c, 150, "R"); cb2.attach(c, 30, "R"); c.ambi = True
S.row([cb2]); print(S.save("carbene_pm.png", dpi=300)); S.check()
print("+ - at", c._amb[2], "deg; formula", hill(formula(cb2)[0]), "charge", formula(cb2)[1])
