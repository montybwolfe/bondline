# Depiction pass: structures from SMILES and molfiles.
# 1. stereo_report reads labels as the atoms they stand for (OMs, OSiMe3, OTBS, CH2OH, CO2Bn...: a label it could not
#    read ranked below H and flipped R/S); a generic R counts as an alkyl and X as a halogen, '?' where it depends on
#    which; formulas of from_smiles's own labels (TIPS, TBDPS, Bpin, ONf).
from bondline import *
from rdkit import Chem
from rdkit.Chem import rdCIPLabeler
def rd(smi):
    m = Chem.MolFromSmiles(smi); rdCIPLabeler.AssignCIPLabels(m)
    return "".join(a.GetProp("_CIPCode") for a in m.GetAtoms() if a.HasProp("_CIPCode"))
for smi in ("CS(=O)(=O)OC[C@H](C)N", "C[Si](C)(C)OC[C@H](C)N", "C[C@H](O[Si](C)(C)C(C)(C)C)CO", "OC[C@@H](C)C(=O)OCc1ccccc1"):
    for ab in ("auto", True, False):
        m = from_smiles(smi, abbrev=ab)
        got = "".join(stereo_report(m, quiet=True)["atoms"].values())
        print(f"{smi} abbrev={ab}: {got} (RDKit {rd(smi)})", [a.label for a in m.atoms if a.drawn and len(a.label) > 2])
m = Mol(); c = m.atom("C", (0, 0)); m.attach(c, 90, "OTr"); m.attach(c, 210, "CH2OTBS"); x = m.attach(c, 330); m.stereo(c, x, "wedge")
stereo_report(m)
for smi, labs in (("C/C=C/*", None), ("O[C@@H](C)*", None), ("O[C@@H](C(C)C)*", None), ("O[C@@H]([1*])[2*]", {2: "R", 3: "R'"}),
                  ("*[C@@H](C)C(=O)O", {0: "X"})):
    stereo_report(from_smiles(smi, labels=labs))
print({k: hill(label_formula(k)) for k in ("OTIPS", "OTBDPS", "Bpin", "ONf", "OMOM", "OTHP")})
a = from_smiles("CS(=O)(=O)OC[C@H](C)N"); b = from_smiles("C[C@H](O[Si](C)(C)C(C)(C)C)CO")
S = Scene(Style()); S.row([a, b], gap=2.0); print(S.save("stereo_labels.png", dpi=300)); S.check()

# 2. an sp carbon between two double bonds in line (CO2, an isocyanate, a ketene, an allene) is drawn as C, so N=C=O
#    no longer reads N=O; a group in line with a double bond counts on neither side of it (the allene's second lines
#    go to the methyls, the ketene's C=C is centred)
S = Scene(Style())
ms = [from_smiles("CN=C=O", align=(1, 3, 0)), from_smiles("CC(C)=C=O", align=(3, 4, 0)), from_smiles("CC=C=CC", align=(1, 3, 0))]
S.row(ms, gap=2.0); print(S.save("cumulenes.png", dpi=300)); S.check()
print("C drawn:", [[a.name for a in m.atoms if a.drawn and a.label == "C"] for m in ms],
      "sides:", [[S._dside(b) for b in m.bonds if b.order == 2] for m in ms])

# 3. slips get a message in chemical terms: labels= on an atom the SMILES lacks, align= on an atom inside an
#    abbreviation, pins RDKit cannot lay the rest out around (and RDKit's own banner stays quiet); labels= win over
#    abbrev (the user's Me's keep the tBu drawn, so the O-ethyl contracts instead)
for kw in (dict(labels={9: "Me"}), dict(align=(2, 4, 0), abbrev=True)):
    try: from_smiles("CC(C)(C)OCC", **kw)
    except ValueError as e: print("ValueError:", e)
try: from_smiles("Cc1ccc(cc1)S(Cl)(=O)=O", xy={5: (0, 0), 8: (0, 1), 9: (0.87, 0.5), 10: (0.87, -0.5), 7: (-1, 0)})
except ValueError as e: print("ValueError:", str(e)[:120])
print([a.label for a in from_smiles("CC(C)(C)OCC", labels={0: "Me", 2: "Me", 3: "Me"}, abbrev=True, align=(1, 4, 30)).atoms])

# 4. where bonds cross, the deeper one is broken only if 3 pt of it keeps showing either side of the gap; with the
#    crossing 2 pt from an atom the other bond is broken instead (an atom never looks unbonded); with no depth known a
#    dashed partial bond passes under a solid one, its pieces holding whole dashes
import numpy as np
def bow(t):
    X = np.array([t, 0.0]); d = np.array([0.5, 0.866])
    m = from_smiles("CCCCC", xy={0: (0, 0), 1: (1, 0), 2: (1.2, 0.98), 3: tuple(X + 0.5 * d), 4: tuple(X - 0.5 * d)})
    for a, z in zip(m.atoms, (-1, -1, 0, 1, 1)): a.z = z
    return m
ms = [bow(0.14), bow(0.5), bow(0.5)]
for a in ms[2].atoms: a.z = None                            # no depth known: the dashed bond is the one broken
ms[2].find(ms[2][3], ms[2][4]).style = "dash"
S = Scene(Style()); S.row(ms, gap=1.5); print(S.save("crossings.png", dpi=300)); S.check()
for m in ms:
    print("broken:", [(b.a.name, b.b.name, [round(float(g[3]), 1) for g in b._gaps]) for b in m.bonds if b._gaps],
          "dashed pieces:", [len(b._prims) for b in m.bonds if b.style == "dash"])

# 5. like= keeps every E/Z of the SMILES: the methyl that takes a (Z)-vinyl bromide's Br stands where the Br stood
vb = from_smiles("CCCC/C=C\\Br")
print([sorted(stereo_report(from_smiles(s, like=vb), quiet=True)["bonds"].values()) for s in ("CCCC/C=C\\C", "CCCC/C=C/C", "CCCC/C=C\\C=C\\C(=O)OC")])

# 6. a ring atom whose two exocyclic groups sit as a gem pair (explicit H on a small ring) takes a wedge and a hash
#    where one mark would leave its three plain bonds within 180 deg (check() called from_smiles's own drawing ambiguous)
for smi in ("BrC1(Br)[C@@]([H])(C)[C@]1([H])C", "[H][C@@]1(CCC#N)CC[C@]1([H])CCC#N"):
    m = from_smiles(smi, explicit_h=True); S = Scene(Style()); S.row([m]); S.check()
    print(sorted(b.style for b in m.bonds if b.style != "plain"), "".join(stereo_report(m, quiet=True)["atoms"].values()), rd(smi))

# 7. ring-fusion H's (the steroid convention) on any fusion with room across its near-vertical bond (4/7, 5/7, not only
#    6/6), in the middle of the free arc; a like= product takes them as its starting material does; check() does not
#    call a fusion H "into the structure" for standing 0.85 BU from its own carbon's ring neighbours
fa = from_smiles("O=C1CC[C@@H]2CC[C@@H]2CC1", align=(4, 7, 90))
sub = from_smiles("CC(=O)[C@@H]1Cc2ccccc2[C@H]2OCOC[C@H]21", rotate=180)
prd = from_smiles("CC(=O)[C@@H]1[C@H](C)c2ccccc2[C@H]2OCOC[C@H]21", like=sub)
S = Scene(Style()); S.row([fa, sub, S.arrow(), prd]); print(S.save("fusion_h.png", dpi=300)); S.check()
print([sum(1 for a in m.atoms if a.label == "H") for m in (fa, sub, prd)])

# 8. a label clash is relieved without breaking 120 deg: an NH between two bonds has its H stacked (not running into
#    the next label), a group on a chain centre goes to the zigzag's other side, only a four-bond centre swings one
def angles(m, k):
    c = m[k]; return sorted(round(ang(n.xy - c.xy) % 360) for n in m.nbrs(c))
ms = [from_smiles("O=C(N*)c1ccccc1", labels={3: "R"}, align=(4, 1, -30)), from_smiles("CNC(=O)OCC", labels={0: "R", 6: "R'"}),
      from_smiles("CCC(CO)(CO)C=O", align=(7, 1, 180, 3)), from_smiles("CCOC=C(C(=O)OCC)C(=O)c1ccccc1", abbrev=True, align=(10, 11, 90))]
S = Scene(Style()); S.row(ms, gap=1.5); print(S.save("angles.png", dpi=300)); S.check()
print(angles(ms[0], 2), angles(ms[1], 1), angles(ms[2], 3), angles(ms[2], 5), angles(ms[3], 3))

# 9. abbreviations: protecting-group ethers with abbrev=True (OMOM before its OMe); alkyl and acyl names when named
#    (tBu, Ac on carbon, Et never the end of a propyl); a ring N-methyl as NMe, the N keeping its ring bonds
ms = []
for smi, ab in (("CC(C)(C)c1ccc(C(C)=O)cc1", ["tBu", "Ac"]), ("COc1ccccc1OCOC", True), ("CN1CCC(CC1)C(=O)OC", ["NMe", "CO2Me"]),
                ("CCCC1CCCCC1", ["Et"]), ("CCc1ccccc1", ["Et"])):
    m = from_smiles(smi, abbrev=ab); ms.append(m); print(smi, ab, sorted(a.label for a in m.atoms if a.drawn))
S = Scene(Style()); S.row(ms[:3], gap=1.5); print(S.save("abbrev.png", dpi=300)); S.check()
print("formula NMe ring:", hill(formula(ms[2])[0]))

# 10. check(): labels of unbonded atoms on one line under a word space apart read as one (reported); a metal's wedges
#     show perspective, not a stereocentre; four ring bonds at one atom (an angular triquinane) may meet at 36 deg
for dx in (1.15, 1.6):
    m = from_smiles("OCCCCl", xy={0: (0, 1), 1: (0, 0), 2: (0.6, -0.5), 3: (dx, 0), 4: (dx, 1)})
    S = Scene(Style()); S.row([m]); S.check()
pd = Mol(); M_ = pd.atom("Pd", (0, 0))
for t, lab, sty in ((90, "PPh3", "plain"), (210, "Cl", "plain"), (330, "Me", "hash"), (270, "Ph", "bold")):
    pd.find(M_, pd.attach(M_, t, lab, style=sty, length=1.3)).free = True
tq = from_smiles("C1CC2CCC3CCCC23C1")
S = Scene(Style()); S.row([pd, tq], gap=2.0); print(S.save("checks.png", dpi=300)); S.check()

# 11. check(): a bond passing a label keeps at least a bond end's own margin (1.6 pt) of white from its letters ("OTBS|"
#     read as touching); a bare vertex drawn on a bond it is not part of is reported (a ring swung onto another)
k5 = from_smiles("CCC(=O)C[C@@H](O)CC(C)C", align=(1, 9, 0, 3))
lb = from_smiles("CC=C(OS(=O)(=O)C(F)(F)F)C[C@@H](O[Si](C)(C)C(C)(C)C)CC(C)C", like=k5)
vb = from_smiles("CC1CCCCC1"); vb[0].xy = (vb[3].xy + vb[4].xy) / 2; vb.find(vb[0], vb[1]).free = True
for m in (lb, vb):
    S = Scene(Style()); S.row([m]); S.check()

# 12. check(): an sp centre (two double bonds, or a triple bond) off a small ring must be straight: O=N+=O left at
#     120 deg after the leaving group goes is reported; a straight nitrile is fine
no2 = Mol(); n = no2.atom("N", (0, 0)); n.charge = 1; no2.attach(n, 90, "O", 2); no2.attach(n, 210, "O", 2)
nit = from_smiles("CC#N")
for m in (no2, nit):
    S = Scene(Style()); S.row([m]); S.check()

# 13. like= maps the shared part the right way round: the embedding whose bond orders agree (the alkyl bromide end onto
#     the alkyl end, not the vinyl one), and a labelled group (CF3) counted as the atoms it stands for, so the ring is
#     not mirrored and the chain not unpinned
sm = from_smiles("BrCCC/C=C/Br"); pr = from_smiles("[*]CCC/C=C/Br", labels={0: "Nu"}, like=sm)
side = lambda m, a, b: "same" if (m[a].xy[0] - m[b].xy[0]) > 0 else "reversed"
tg = from_smiles("OCCCc1cccc(c1)C(F)(F)F", align=(0, 4, 0)); fr = from_smiles("NCCCc1cccc(c1)C(F)(F)F", like=tg)
cf = lambda m: [a for a in m.atoms if a.label == "CF3"][0]
print("vinyl end:", side(sm, 6, 0), side(pr, 6, 0), "| CF3 offset:", [tuple(round(float(v), 2) for v in cf(m).xy - m[4].xy) for m in (tg, fr)])

# 14. like= never makes a drawing worse than its own layout: where fitting the shared atoms of a cage bends the rest out
#     of shape, the molecule's own layout is turned onto them instead
sm = from_smiles("O=CCCC12CC(CC(=O)C1=C)C(=C)C2"); pr = from_smiles("O=C1CCC23CC(CC(=O)C2C1)C(=C)C3", like=sm)
L = [norm(b.a.xy - b.b.xy) for b in pr.bonds]; print("like= bridged product bonds: %.2f-%.2f BU" % (min(L), max(L)))

# 15. pins: an unpinned group RDKit lays onto a ring bond or onto its neighbour's group is turned into the free sector
#     round its pinned atom with the most room (a lone group anywhere 45 deg clear of the bonds; a gem pair split)
pp = from_smiles("OC1(CC#CC)CCC=C1", xy={1: (0, 0), 2: (0.3, -0.95), 3: (1.25, -0.6), 4: (2.1, -0.2), 5: (2.95, 0.2),
                                        6: (-0.95, 0.3), 7: (-0.6, 1.25), 8: (0.4, 1.3), 9: (0.6, 0.35)})
S = Scene(Style()); S.row([pp]); print(S.save("pins_free.png", dpi=300)); S.check()
c1 = pp[1]; print("OH at", round(ang(pp[0].xy - c1.xy) % 360), "deg; bonds at", sorted(round(ang(n.xy - c1.xy) % 360) for n in pp.nbrs(c1) if n is not pp[0]))

# 16. stereo_report keeps an alkene's E/Z when the neighbouring ring stereocentre's wedge sits on the bond to it (RDKit
#     keeps one direction a bond, so the alkene's geometry is read off the drawing)
for smi in ("CCC/C=C/[C@H]1CC[C@@H](O)[C@@H]1C", "CCC/C=C\\[C@H]1CC[C@@H](O)[C@@H]1C"):
    print(smi, sorted(stereo_report(from_smiles(smi), quiet=True)["bonds"].values()), rd(smi))

# 17. a salt's counter-ion stands beside the opposite charge, 2 BU or more out along its free direction (not 1 BU off,
#     where check() reads it as a group left behind)
sa = [from_smiles("[Cl-].C[N+](C)(C)Cc1ccccc1"), from_smiles("CC(=O)[O-].[Na+]")]
S = Scene(Style()); S.row(sa, gap=1.5); print(S.save("salts.png", dpi=300)); S.check()
for m in sa:
    q = [a for a in m.atoms if a.charge]; print([(a.label, round(float(norm(q[0].xy - q[1].xy)), 2)) for a in q[:1]])

# 18. a like= product keeps its drawing's wedge bond at each shared stereocentre (the ether wedged on C-O, like the
#     alcohol), the mark chosen to keep the configuration
for m in (k5, lb):
    c = next(iter(stereo_report(m, quiet=True)["atoms"]))
    print([(n.label, m.find(c, n).style) for n in m.nbrs(c)], stereo_report(m, quiet=True)["atoms"][c])

# 19. check(): at a reflex ring vertex the exterior bisector points away from the ring's centre (the C=O of a ring opened
#     with its old geometry kept is not "180 deg off"); a bond from_smiles's own layout stretched gets advice that fits
ro = from_smiles("CCC1=C2C(CC1=O)CC1(O)CCCC21")
ro.find(ro[9], ro[10]).order = 2; ro[10].label = "O"; ro.unbond(ro[9], ro[14]); ro[14].charge = -1
st_ = from_smiles("CCCC"); st_[3].xy = st_[2].xy + (st_[3].xy - st_[2].xy) * 1.3
for m in (ro, st_):
    S = Scene(Style()); S.row([m]); S.check()
