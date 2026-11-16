# Labels as chemists write them: NC bonded through N on the right of a bond (an isocyanide) but a cyano group written
# leftwards; acronyms (FG, PG, EWG) and lowercase words (Nuc, Bpin, Phth) kept whole; Ln with a subscript n; ligand-first
# metal labels ((PPh3)4Pd^0, L2(Ph)Pd^{II}, LnNi^{II}, PdII) anchored on the metal, printed as typed, with automatic
# subscripts beside explicit markup; EtOH and iPr2NEt anchored on their O and N; a whole-label <b>[M]</b> keeps its
# anchor; a lone formula ion (PF6-) is never written backwards; NiPr2 / PtBu3 read as N(iPr)2 / P(tBu)3.
from bondline import *
from bondline import _anchor, label_units

def printed(S, a):
    """The label as drawn, left to right: _x a subscript, ^x a superscript, `x an italic prefix."""
    out = []
    for g in a._g:
        for t, x, y, sz, it, bd in g["runs"]:
            mark = "" if sz >= S.st.size - 1e-6 else ("`" if it else "^" if y > g["y"] + 0.01 else "_")
            out.append((round(x, 2), mark + t))
    return "".join(t for x, t in sorted(out))

for l in ("NC", "FG", "EWG", "Nuc", "Ln", "Bpin", "NPhth", "(PPh3)4Pd^0", "L2(Ph)Pd^{II}", "L_nNi^{II}", "PdII(Ph)L2",
          "LnPdX", "EtOH", "iPr2NEt", "NiPr2", "PtBu3", "SnBu3", "PF6", "(CH2)11CH3"):
    st, tl, pv = _anchor(l)
    print(f"{l:14} stored {st:16} {'typed leftwards' if tl else 'pivot %d' % pv if pv else ''}  units {label_units(st)}")

S = Scene(Style())
iso = Mol(); q = iso.atom("C", (0, 0)); nc = iso.attach(q, 0, "NC"); iso.attach(q, 120); iso.attach(q, 240); iso.attach(q, 180)
cy = Mol(); c = cy.atom("C", (0, 0)); cn = cy.attach(c, 150, "NC"); cy.attach(c, 30, "CO2Et")
pd = Mol(); p = pd.atom("L2(Ph)Pd^{II}", (0, 0)); pd.attach(p, 0, "OR", length=1.6, free=True); pd.attach(p, 90, "Ar", length=1.5, free=True)
ni = Mol(); n = ni.atom("L_nNi^{II}", (0, 0)); ni.ligate(n, "sqpl", S="Ar")
bm = from_smiles("Cc1ccccc1", labels={0: "<b>[M]</b>"}, align=(1, 0, 30))
fg = Mol(); f0 = fg.atom("C", (0, 0)); f = fg.attach(f0, 150, "FG"); nu = fg.attach(f0, 30, "Nuc")
et = Mol(); o = et.atom("EtOH", (0, 0)); o.lp(90)
hb = Mol(); h = hb.atom("iPr2NEt", (0, 0)); h.lp(90)
cat = from_smiles("C[P+](C)(C)C"); pf = cat.atom("PF6", (0, 0)); pf.charge = -1; cat.dock(pf, cat[1], 150, 2.8)
S.row([iso, cy, pd, ni, bm, fg, et, hb, cat], gap=1.2)
print(S.save("labels.png", dpi=300)); S.check()
print("NC on the right:", nc.label, "| NC on the left:", cn.label)
for a in (nc, cn, p, n, bm[0], f, nu, o, h, pf):                 # where the atom (the bonds' end) sits in its label
    A = S.P(a); X = np.vstack([h_ for g in a._g for h_ in g["hulls"]])[:, 0]
    print(f"{printed(S, a):18} ink {A[0] - X.min():5.1f} pt left of the atom, {X.max() - A[0]:5.1f} pt right")

print("== text: formulas beside markup, tags, glued charges and radicals, unicode signs, pKa', hyphens, adducts, colour")
from bondline import _typo, parse_rich, _hull_dist
for t in ("H3O⊕", "H3O^+", "N3- I3- Fe3+", "NMe3⁺ N3⁻", "R2CH•", "<i>−CO2</i>", "<b>Ru(bpy)3.Cl2</b>", "pKa' (K2CO3)",
          "PdCl2[P(o-tolyl)3]2", "NaOt-Bu", "LnPd(0)", "MgSO4.7H2O"):
    print(f"  {t:22} {_typo(t, chem=True)}")
print("  colours:", [(r[0], getattr(r[0], "c", None)) for r in parse_rich("<c=#E8261A>Ligand</c>, Base")])
S = Scene(Style()); a = from_smiles("Brc1ccccc1"); b = from_smiles("c1ccc(cc1)-c1ccccc1")
S.caption(b, "<c=#1F4FD1><b>biphenyl</b>\nC12H10</c>\nmp 69 °C")
S.row([a, S.arrow(above=["<c=#E8261A>Transition metal</c>", "<c=#2E7D32>Ligand</c>, Base"], below="NaOt-Bu"), b])
print(S.save("colour_text.png", dpi=300)); S.check()

print("== check(): partial bonds of a transition state, H-bonds, oniums, X as a heteroatom, attach(length=), formula ions")
S = Scene(Style())
ts = Mol(); c = ts.atom("C", (0, 0), show=True); ts.attach(c, 90, "H"); ts.attach(c, 240, "H"); ts.attach(c, 300, "H")
ts.attach(c, 180, "Br", style="dash"); ts.attach(c, 0, "I", style="dash")                 # SN2 TS: C 3 + 2 x 1/2
fc = Mol(); pd = fc.atom("Pd", (0, 0)); o = fc.attach(pd, 0, "OR", length=1.6, style="dash")   # four-centred TS: B 2 + 2 x 1/2
b = fc.attach(o, 90, "B", length=1.5, style="dash"); fc.attach(b, 30, "OH"); fc.attach(b, -30, "OH")
ar = fc.attach(b, 180, "Ar", length=1.6, style="dash"); fc.bond(pd, ar, style="dash")
ox = Mol(); q = ox.atom("C", (0, 0)); ox.attach(q, 90, "OH"); ox.attach(q, 210); e = ox.attach(q, -30, "OEt"); ox.attach(e, 270, "H"); e.charge = 1
am = Mol(); q = am.atom("C", (0, 0)); am.attach(q, 90, "O").charge = -1; n = am.attach(q, 210, "NH2"); n.charge = 1; am.attach(n, 150, "R")
fu = from_smiles("c1ccoc1", labels={3: "X"}, align=(3, 0, 90))
pa = from_smiles("C1CCC2c3ccccc3CCC2C1")
for i, j in ((11, 12), (12, 13)): bd = pa.find(pa[i], pa[j]); bd.order = 2; bd.style = "partial"
ln = Mol(); q = ln.atom("C", (0, 0)); ln.attach(q, 30, "R1", length=1.4); ln.attach(q, 150, "C6H13")
ion = Mol(); ion.atom("Ru(bpy)_3^{2+}", (0, 0))
S.brackets([ts], charge="‡"); S.brackets([fc], charge="‡")
S.row([ts, fc, ox, am, fu, pa, ln, ion], gap=1.2); print(S.save("check_rules.png", dpi=300)); S.check()
print("  C6H13 on the left of its bond:", printed(S, ln.atoms[2]))
S = Scene(Style()); bad = Mol(); q = bad.atom("C", (0, 0)); bad.attach(q, 0, "OMe"); bad.attach(bad.atoms[1], 90)
nh = Mol(); q = nh.atom("C", (0, 0)); n = nh.attach(q, 0, "NH2"); nh.attach(n, 60, "R")             # negatives: still refused
S.row([bad, nh], gap=1.5); S.save("check_neg.png", dpi=300); S.check()

print("== bonds between labels show at least 6 pt (a terminal label slides out); a level bond starts beyond a superscript")
def vis(b): P = np.asarray(b._prims[0][1], float); return float(np.hypot(*(P[-1] - P[0])))
S = Scene(Style())
hcl = Mol(); h = hcl.atom("H", (0, 0)); cl = hcl.attach(h, 0, "Cl")
pdh = Mol(); p = pdh.atom("Pd^{II}", (0, 0)); pdh.ligate(p, "sqplx", NW="PPh3", SW="PPh3", SE="Cl"); hh = pdh.atom("H", (2, 2)); pdh.join(p, hh)
ni = Mol(); n = ni.atom("Ni^{III}", (0, 0)); ar = ni.attach(n, 10, "Ar", length=1.3); ni.attach(n, 90, "R", length=1.3)
st = Mol(); c = st.atom("C", (0, 0), show=True); st.attach(c, 120, "R3"); st.attach(c, 0, "Br"); r2 = st.attach(c, 240, "R2"); r1 = st.attach(c, 300, "R1")
st.stereo(c, r2, "hash"); st.stereo(c, r1, "wedge")
lg = Mol(); q = lg.atom("Pd", (0, 0)); cs = lg.ligate(q, "sqpl", N="Cl", S="Cl"); lg.ligate(q, "sqpl", E="P", W="P", length=3.4)
lb = Mol(); q2 = lb.atom("L_nNi^{II}", (0, 0)); bb = lb.ligate(q2, "sqpl", E="Br", S="Ar"); q2.number = 8
S.row([hcl, pdh, ni, st, lg, lb], gap=1.2); print(S.save("label_bonds.png", dpi=300)); S.check()
print(f"  H-Cl {vis(hcl.find(h, cl)):.1f} pt; Pd-H {vis(pdh.find(p, hh)):.1f}; Ni(III)-Ar {vis(ni.find(n, ar)):.1f}; "
      f"C-R2 hash {norm(r2.xy - c.xy):.2f} BU, C-R1 wedge {norm(r1.xy - c.xy):.2f} BU; Pd-Cl {norm(cs['N'].xy - q.xy):.2f} BU")
sup = [h_ for g in q2._g for h_, s_ in zip(g["hulls"], g["sup"]) if s_][0]
br = lb.find(q2, bb["E"]); P0 = np.asarray(br._prims[0][1][0], float)
print(f"  Ni(II)-Br starts {P0[0] - sup[:, 0].max():.1f} pt beyond the superscript; locant clear of it: "
      f"{_hull_dist(q2._num[0], sup) - q2._num[1] > 0.5}")
