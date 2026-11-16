# Stage 7: R-group numbers and superscripts, chair templates, colour, charges, output (PNG resolution, paged PDF)
import hashlib
from bondline import *
from PIL import Image

def shown(a, S):
    """The label as drawn: superscripts and subscripts written with Unicode, from the laid-out runs."""
    sup = dict(zip("0123456789+−n", "⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻ⁿ")); sub = dict(zip("0123456789", "₀₁₂₃₄₅₆₇₈₉"))
    runs = sorted((r for g in a._g for r in g["runs"]), key=lambda r: r[1]); out = ""
    for t, x, y, sz, it, bd in runs:
        if sz < S.st.size and not it:
            up = y > min(r[2] for r in runs) + 0.01
            out += "".join((sup if up else sub).get(c, c) for c in t)
        else: out += t
    return out

print("== numbered groups are superscripts, counts are subscripts")
S = Scene()
k1 = Mol(); c = k1.atom("C", (0, 0)); b = k1.attach(c, 30); k1.attach(b, 90, "OR1"); k1.attach(b, -30, "CO2R2")
k1.attach(c, 150, "Pd0"); k1.attach(c, 270, "NR1R2")
k2 = Mol(); c = k2.atom("C", (0, 0)); b = k2.attach(c, 30); k2.attach(c, 150, "NR1R2"); k2.attach(b, -30, "NR2")
k2.attach(c, 270, "CO2R3"); k2.attach(b, 90, "SiR3")
k3 = Mol(); pd = k3.atom("Pd", (0, 0)); k3.attach(pd, 0, "PAr3", length=1.3, free=True); k3.attach(pd, 180, "Cl", length=1.3, free=True)
S.row([k1, k2, k3]); S.check(); S.save("labels.png")
print("  " + "  ".join(shown(a, S) for m in (k1, k2, k3) for a in m.atoms if a.drawn))

print("== a bond leaving R3 between 0 and 90 deg: up to about 57 deg it starts past the 3, steeper the 3 moves aside")
S = Scene(); ms = []
for t in (0, 15, 30, 45, 54, 60, 75, 90):
    m = Mol(); c = m.atom("C", (0, 0)); r = m.attach(c, t + 180, "R3"); m.attach(c, t + 60); ms.append(m)
S.row(ms); iss = S.check(quiet=True); S.save("rgroups.png")
for t, m in zip((0, 15, 30, 45, 54, 60, 75, 90), ms):
    r = m.atoms[1]; R = [g for g in r._g][0]["runs"]; xs = [x for t_, x, y, sz, it, bd in R]
    b = m.bonds[0]; p, q = b._prims[0][1]
    print(f"  {t:3d} deg: superscript {xs[1] - xs[0]:.2f} pt after the R, bond shows {norm(q - p):.1f} pt")
print("  issues:", len(iss))

print("== chair and ring-flipped chair: axial and equatorial swap on the same face")
S = Scene(); m = Mol(); at, ax, eq = m.chair((0, 0)); m.attach(at[0], eq[0], "Me"); m.attach(at[3], eq[3], "tBu")
m2 = Mol(); at2, ax2, eq2 = m2.chair_flip((0, 0)); m2.attach(at2[0], ax2[0], "Me"); m2.attach(at2[3], ax2[3], "tBu")
S.row([m, S.arrow(kind="<=>"), m2]); S.check(); S.save("chairs.png")
print("  axial  ", [round(x) for x in ax], " flipped", [round(x) for x in ax2])
print("  equat. ", [round(x) for x in eq], " flipped", [round(x) for x in eq2])

print("== colour: a whole molecule, a fragment (split joint), an atom; plain charges in ink")
S = Scene(Style(hetero=True, charge="plain"))
en = Mol(); c1 = en.atom("C", (0, 0)); c2 = en.attach(c1, 30, "C", 2); o = en.attach(c2, 90, "O"); o.charge = -1; o.add_lps(3)
en.attach(c1, 210); en.attach(c2, -30, "Me"); en.color("#1F4FD1")
ch = Mol(); a0 = ch.atom("C", (0, 0)); a1 = ch.attach(a0, 30); a2 = ch.attach(a1, -30); a3 = ch.attach(a2, 30); a4 = ch.attach(a3, -30, "OH")
ch.colour("#D01C1C", [a2, a3, a4])
ac = Mol(); c = ac.atom("C", (0, 0)); ac.attach(c, 90, "O", 2); o2 = ac.attach(c, -30, "O"); o2.charge = -1; ac.attach(c, 210, "Me")
S.row([en, ch, ac]); S.check(); S.save("colour.png")
pr, tx = S._build()
print("  ink colours used:", sorted({p[3] if p[0] != "poly" else p[2] for p in pr}), sorted({t["color"] for t in tx}))

print("== output: PNG at 1200 ppi of its true size; print-ready PDF pages, split between rows")
S = Scene(); m = Mol(); r = m.ring(6); m.attach(r[0], 90, "OH"); S.row([m]); w, h = S.save("ppi.png")
im = Image.open("ppi.png"); print(f"  {w:.1f} x {h:.1f} pt -> {im.size[0]} x {im.size[1]} px at {round(im.info['dpi'][0])} ppi")
S.save("ppi_px.png", px=1000); print("  px=1000 ->", Image.open("ppi_px.png").size)
S = Scene(); items = []
for i in range(60):
    if i: items.append(S.arrow())
    mm = Mol(); rr = mm.ring(6); mm.attach(rr[0], 90, ["OH", "Cl", "OMe"][i % 3]); items.append(mm)
S.flow(items)
n1 = S.save_pdf("long.pdf"); n2 = S.save_pdf("long2.pdf")
same = hashlib.sha256(open("long.pdf", "rb").read()).digest() == hashlib.sha256(open("long2.pdf", "rb").read()).digest()
print(f"  {n1} pages; second save identical: {same}")

print("== a reagent and a by-product on one arc: the reagent joins before the by-product leaves, never crossing")
S = Scene()
def pdx(X):
    m = Mol(); P = m.atom("Pd", (0, 0)); m.ligate(P, "sqpl", N="L", S="L", E=X, W="Ar"); return m
p0 = Mol(); P = p0.atom("Pd^0", (0, 0)); p0.ligate(P, "lin", E="L", W="L")
c = S.cycle([(p0, 90), (pdx("Br"), -30), (pdx("Ar'"), 210)], radius=4)
S.branch(c, 0, "Ar–Br", kind="in"); S.branch(c, 1, "Ar'B(OH)_2 + base", kind="in"); S.branch(c, 1, "Br^−", kind="out")
S.branch(c, 2, "Ar–Ar'", kind="out")
S.save("cycle_inout.png"); S.check()
