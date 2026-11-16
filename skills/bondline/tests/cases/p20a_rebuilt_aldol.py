# Acceptance : the four aldol figures and repro_1 rebuilt with the intent helpers (figures 1 and 2).
# No hand-placed coordinates: reagents are placed by deprotonate / protonate / attack at the standard distances, bonds
# are made and broken with join / release, products are tidied to the standard depiction (OH labels) and by-products
# stand behind a "+". Printed: check(), validate(), continuity(), every arrow's shaft and sweep, the reagent distances
# against the rule, and how far each structure's main fragment sits from its row's arrow line.
from bondline import *
st = Style()

def report(S, rows, pairs=()):
    sh = []
    for c in S.curlies:
        L = float(np.hypot(*np.diff(c._arc, axis=0).T).sum()) / st.bond; sh.append((L, c._sweep))
    if sh: print(f"  arrows: {len(sh)}, shortest shaft {min(L for L, _ in sh):.2f} BU, smallest sweep {min(w for _, w in sh):.0f} deg,"
                 f" all >= 0.8 BU and >= 25 deg: {all(L >= 0.8 and w >= 25 for L, w in sh)}")
    off = []
    for r in rows:
        line = [it.y for it in r if isinstance(it, RArrow) and not it.vertical]
        if not line: line = [S._lboxes(it)[2] for it in r if isinstance(it, Text)]
        off += [abs(S._lboxes(it)[2] - line[0]) for it in r if isinstance(it, Group)]
    print(f"  structures off their row's arrow line: at most {max(off):.2f} pt")
    for what, a, b in pairs:
        print(f"  {what}: {norm(a.xy - b.xy):.2f} BU (rule {Mol._dock_dist(a, b):.1f})")

def ethanal(n=""):
    m = Mol(); me = m.atom("C", (0, 0), name="Ca" + n); c = m.attach(me, 30, name="C1" + n)
    m.attach(c, -30, "O", 2, name="O1" + n); return m

def aldol_steps(enolisable, electrophile, S):
    """Base-catalysed aldol addition, written forwards: s1 deprotonation, s2 C-C bond formation, s3 protonation,
    s4 product + hydroxide."""
    s1 = enolisable; Ha = s1.attach(s1["Ca"], None, "H", name="Ha")
    B = s1.deprotonate(hydroxide(), Ha); B.name = "B"
    S.curly(B.lp(), Ha); S.curly((s1["Ca"], Ha), (s1["Ca"], s1["C1"])); S.curly((s1["C1"], s1["O1"]), s1["O1"])
    s2 = s1.copy(); s2.release(s2["Ha"], s2["Ca"]); s2.remove_frag(s2["B"])
    s2.find(s2["Ca"], s2["C1"]).order = 2; s2.find(s2["C1"], s2["O1"]).order = 1; s2["O1"].charge = -1
    s2.merge(electrophile); s2.attack(s2["Ca"], s2["C1b"])
    S.curly(s2["O1"].lp(), (s2["O1"], s2["C1"])); S.curly((s2["Ca"], s2["C1"]), s2["C1b"]); S.curly((s2["C1b"], s2["O1b"]), s2["O1b"])
    s3 = s2.copy(); s3.find(s3["Ca"], s3["C1"]).order = 1; s3.find(s3["C1"], s3["O1"]).order = 2; s3["O1"].charge = 0
    s3.join(s3["Ca"], s3["C1b"]); s3.find(s3["C1b"], s3["O1b"]).order = 1; s3["O1b"].charge = -1
    w = water(); Hw = s3.protonate(w.h, s3["O1b"]); Hw.name = "Hw"; w.key.name = "Ow"
    S.curly(s3["O1b"].lp(), Hw); S.curly((Hw, s3["Ow"]), s3["Ow"])
    s4 = s3.copy(); oh = s4.release(s4["Ow"], s4["Hw"]); s4.join(s4["O1b"], s4["Hw"]); s4["O1b"].charge = 0
    oh["Ow"].charge = -1; s4.tidy(); oh.tidy()
    S.validate(s1, s2); S.validate(s2, s3); S.validate(s3, [s4, oh]); Scene.continuity(s1, s2, s3, [s4, oh])
    pairs = [("base to H", s1["B"], s1["Ha"]), ("enolate C to carbonyl C", s2["Ca"], s2["C1b"]), ("O- to water H", s3["O1b"], s3["Hw"])]
    print("  product O label:", s4["O1b"].label, "| by-product:", " ".join(a.label for a in oh.atoms), oh["Ow"].charge)
    return s1, s2, s3, s4, oh, pairs

print("== figure 1: aldol scheme")
S = Scene(st); a1, a2 = ethanal(), ethanal("x")
p = Mol(); c4 = p.atom("C", (0, 0)); c3 = p.attach(c4, 30); p.attach(c3, 90, "OH"); c2 = p.attach(c3, -30)
c1 = p.attach(c2, 30); p.attach(c1, -30, "O", 2)
S.caption(a1, "ethanal"); S.caption(a2, "ethanal"); S.caption(p, "3-hydroxybutanal")
rows = [S.row([a1, "+", a2, S.arrow(above="NaOH (aq)", below="5 °C"), p])]
print(S.save("aldol_scheme.png")); S.check(); report(S, rows)

print("== figure 2: aldol mechanism (also repro_1 rebuilt)")
S = Scene(st)
s1, s2, s3, s4, oh, pairs = aldol_steps(ethanal(), ethanal("b"), S)
rows = S.flow([s1, S.arrow(below="−H_2O", step=1), s2, S.arrow(step=2), s3, S.arrow(step=3), s4, "+", oh], width_cm=16)
print(S.save("aldol_mechanism.png")); S.check(); report(S, rows, pairs)

