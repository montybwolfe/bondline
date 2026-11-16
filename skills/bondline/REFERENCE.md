# bondline API reference

Read the part you need while building; SKILL.md says how to work. Units: positions in bond units (BU, 1 BU = 14.4 pt
bond), angles in degrees counter-clockwise from east. Every script starts `from bondline import *` and draws into one
`S = Scene(Style())`. The ten scripts in `examples/` are complete starting points.

## Structures

- `m = Mol()`; `a = m.atom(label="C", xy=(0, 0), name=None, show=False)`; `m.bond(a, b, order=1, style="plain", free=False)`;
  styles `plain`, `wedge`, `hash`, `bold`, `dash` (partial), `dot` (coordination, forming bonds in a TS), `dative`,
  `wavy`, `partial`; `free=True` for a bond drawn at a deliberate length (attach(length=) sets it). check() counts a
  dashed or dotted bond as a partial one (a transition state's half bond, an H-bond). A bond between two labels shows
  at least 6 pt of line: a label with no other bond (O-H, H-Cl, Pd-H, C-R¹) is moved out along it at drawing time
  (not one a curly arrow uses); a level bond leaves a label past its superscript (Pd^II-OR). `m.find(a, b)` gives the bond;
  `m[name]` the atom (merge() keeps names unique: a merged atom whose name is taken gets a prime, 8 → "8'");
  `m.nbrs(a)`; `m.unbond(a, b)`; `m.remove(a)`.
- `m.attach(a, angle=None, label="C", order=1, style="plain", length=1.0, name=None, toward=None)`: a new atom on the
  standard outward direction (zigzag, or the exterior bisector on a ring); a fourth substituent on a chain carbon
  makes a gem pair. `toward=partner` keeps a direction free for an incoming reagent. `m.chain(a, angles, labels,
  orders)` for runs of atoms.
- Rings: `m.ring(n, center, rot=90, labels={i: "O"}, kekule=False)`; `m.ring_out(a, angle, n)` (hanging off a);
  `m.ring_at(a, angle, n)` (a spiro vertex); `m.ring_fuse(a, b, n)` (sharing bond a-b, on the side with more room; `side=1 | -1` for left or right of a->b).
  A floating substituent (any ring position, generic schemes): `m.floating(a, b, "R1")`, its bond across the ring edge a-b.
- `from_smiles(smi, align=(i, j, deg), rotate=0, flip=False, labels=None, explicit_h=False, abbrev="auto",
  like=None)`. Atom names are RDKit indices (`m[3]`). `align` turns the line from atom i to atom j (a bond, or two
  ends of a chain for a level zigzag) to deg; `align=(i, j, deg, k)` also puts atom k on its left (above it when it
  points right), e.g. a ketone's O up. `flip` then mirrors left-right. Ring-fusion stereocentres get their H drawn
  straight up or down on a short wedge or hash (the steroid convention) when the drawing is level. `abbrev="auto"` writes NO2, silyl groups, OTf/OTs/OMs, Ts, CF3, Bpin, OLi, MgBr as labels; `True` also esters,
  ethers (OMOM, OTHP, OPMB too), NHBoc, NMe2, CN...; a list picks names, alkyl and acyl ones among them (tBu, iPr,
  Et, Bn, Ac, Bz, Piv, PMB, Cy) and `"NMe"` for a ring N-methyl (the N keeps its ring bonds); `False` draws every atom; an atom named in `labels=` is never
  contracted. `like=other` lays the part shared with a
  drawn Mol out exactly as there, so products match starting materials. `explicit_h=True` keeps [H] atoms.
- Curated skeletons (`TEMPLATES`): a molecule containing a penam (penicillins: the beta-lactam square upright, S on
  top), a morphinan (morphine, codeine, dextromethorphan: A/B/C flat, the N-bridge and 4,5-ether in perspective, the
  bridge crossing ring C with the back bond broken) or the strychnos ring system (strychnine, brucine: flat, 1 BU bonds)
  is laid out on the textbook drawing of it, the rest round it; `xy=` and `like=` take precedence, `align=`, `rotate=`,
  `flip=` still turn it. On a template its wedges sit on bonds to non-stereocentres where they can; anywhere, a
  stereocentre whose one wedge would leave its other bonds within 180 deg gets it on a bond that leaves them spread.
- Bridged and caged polycycles: `from_smiles` uses RDKit's ring-system template where one exists (norbornane,
  camphor, adamantane, cubane) and keeps the cage's perspective; wherever two bonds cross, the
  one behind is broken (depth read from a 3D conformer fitted to the drawing; `bond.behind = True` decides it). For any
  other perspective drawing (a slide's twin chair, taxol, a sketch) pin the atoms where the picture has them:
  `from_smiles(smi, xy={i: (x, y)})` in BU (y up, 1 BU = one bond of the picture), the rest laid out around them;
  bonds inside the pinned part are drawn plain (the perspective carries the depth); `wedges="exo"` keeps the wedges
  of pinned substituents (the OH or OAc on a pinned cage, as a textbook taxol shows them), `wedges=True` keeps all.
  Pins also fix any conformation `align=` cannot (two rotamers at once, a chelating amide): pin the few atoms that
  set it; an unpinned group that lands on a bond is turned into the free sector beside it.
- `from_molfile(path_or_text, labels=None, align=..., contract=True)` keeps a ChemDraw export's own layout (median
  bond scaled to 1 BU, wedges, charges, abbreviations); `to_molfile(m, path=None)` writes the drawing back out.
- Labels: write them as read (OMe, CO2Me, NHBoc, tBu, R1, Ar); anchoring, flipping (MeO/OMe), subscripts and italic
  prefixes are automatic. Rich text anywhere: `<i>`, `<b>`, `<c=#E8261A>..</c>` (a colour: one word of the
  conditions, the first line of a caption), `_{}`, `^{}`. Abbreviations are single atoms. A metal
  written after its ligands (L2Pd, LnNi^{II}, (PPh3)4Pd^0, L2(Ph)PdII) is anchored on the metal; EtOH, iPr2NEt and
  LnPdX on the atom in the middle; NC on the right of a bond is an isocyanide (write CN for a cyano group there);
  `\Nuc` prints a label exactly as typed; `a.ldir = "left"` (or "right", "up", "down") sets the side a label's
  other letters take when the automatic one collides (CO2H read as HO2C beside an N-H). A locant can be coloured:
  `a.number = "<c=#00897B>3</c>"`.
- Stereo: `m.stereo(centre, other, "wedge" | "hash" | "wavy" | "bold")`, narrow end on the centre;
  `stereo_report(m)` reads R/S and E/Z back, labels counted as the atoms they stand for (OTBS, CH2OH, CO2Bn); R counts
  as an alkyl group and X as a halogen, and a descriptor that depends on which is `?`. `m.flip()` mirrors but keeps
  configurations (it swaps wedges).
- Charges and electrons: `a.charge = -1`; `a.radical = True`; `a.lp()` adds a lone pair (dots only where you add
  them; `a.add_lps(n)`); `a.partial = "+"` (δ+); `a.number = 3` (a locant); `a.charge_at = deg` pins a charge (or,
  with no charge, a δ); `a.ambi = True` the + − of a carbene (empty p, filled sp2; net charge 0).
- Colour: `m.color("#1F4FD1", atoms)` colours a fragment (labels, bonds, marks) through every later step;
  `bond.color = "#E8261A"` one bond (a bond formed, highlighted); `a.color` one atom. `Style(hetero=True)` colours
  heteroatom labels.
- Shapes: `m.chair()`, `m.chair_flip()`, `m.newman(front, back)`; metals `m.ligate(M, "sqpl" | "sqplx" | "oct" |
  "tbp" | "tet" | "lin" | "bent", N="PPh3", ...)`, `m.bridge(a, b, k)` (chelates), `m.eta(M, atoms)`,
  `m.sandwich("Fe", n=5, bend=)`, `m.persp_ring(...)`, `m.deloc(atoms, charge)` (Cp⁻, allyl as one hybrid).
- `m.numbered(start=1)`: a copy with every atom numbered, as a reference drawing for a step list; `m.num(k)` is atom
  k in that molecule and in every step copied from it. `m.copy()`, `m.move(dx, dy)`, `m.rotate(deg)`, `m.center()`.

## Mechanisms

- One Mol per step: `nxt = prev.copy()` then edit; atoms keep uid and position. `Scene.continuity(*steps)` lists any
  that wander.
- Reagents: `hydroxide()`, `water()`, `hydronium()`, `alkoxide("Et")`, `halide("Br")`, `amine("Et3N")`,
  `acid("Cl")` (key atom: its H), `base("B" | "Base" | "Nu", charge=-1)` (a generic B:⁻, `.key`); `m.merge(other)`
  takes any Mol in.
- Placing partners (the added species moves): `m.deprotonate(base, H)` (in line with X-H, 2.0 BU beyond H);
  `m.protonate(H, acceptor)` (H: the acid's H atom, or a reagent Mol: water() gives its `.h`, a hand-built acid its one
  O-H or N-H proton; anything else is refused); `m.attack(nu, E)` (SN2: straight anti to a leaving group; Bürgi–Dunitz
  to C=O; always where its arrow can run straight in, clear of other groups); `m.dock(atom, target, angle, dist)` for
  anything else (an atom of another Mol is merged in, its molecule moving with it).
- Bond changes: `m.join(a, b)` (a bond between fragments on standard angles); `m.release(atom, partner)` (break, and
  return the leaving fragment as its own Mol); `bond`/`unbond`, `split(atom)`, `transfer(H, new_parent, angle)`,
  `remove_frag(atom)`, `swing(atom, pivot, deg)` (open a ring for a ring-forming step: build the product, copy, break,
  swing). Finish with `tidy()` (folds unused explicit H into OH, NH labels, leaving any atom a curly arrow of this step
  uses; `carbon=True` also makes spare C-H
  atoms implicit, e.g. a protonated formaldehyde once it is a CH2 in a chain). `unfold(atom)` is the reverse for one
  label (OH → O-H, NH2 → N with two H) when a step needs those bonds. Breaking the bond to one of a gem pair on a ring
  returns the other substituent to the exterior bisector.
- Curly arrows: `S.curly(src, dst)`, `S.fish(src, dst)`; src `atom.lp()`, `(a, b)` (a bond's electrons) or an atom
  (anion, radical); dst an atom (intermolecular attack: the electrophile), `(a, b)`, or `Mid(a, b)` (a bond forming in
  a concerted step), or a point `(x, y)` in BU for a free end (an arrow drawn wrong on purpose; check() reports it). A bond's electrons onto an atom bonded to one end of it, `S.curly((C2, R), C1)`, is a 1,2-shift:
  R moves to C1 with them (validate() replays it so; `Mid(R, C1)` says the same); if another arrow takes R (a base
  removing that H), they stay and make the C1=C2 bond (E1 drawn onto the cation). Such arrows, and a bond's electrons
  onto its own atom, curl over the outside, leaning away from where they turn back to. Placement is automatic; only if
  `check()` complains, steer with `side="left" | "right"`,
  `bend=deg`, or move the partner. `color=` for a coloured arrow; `dashed=True` for a step that does not happen
  (or a minor pathway).
- Cyclic transition states: build the ring with `m.ring()`, forming and breaking bonds `style="dot"` or `"dash"`;
  `S.brackets([ts], charge="‡")`. The arrows of a cyclic array keep to one side of the ring and turn with the flow
  round it (a 4-centre array has room inside for one: the other runs outside); arrows that are rotated copies of each
  other (benzene's Kekulé arrows) are drawn alike. `count=2` (or `"n"`) sets a subscript at the lower right: a dimer [ ]₂, a polymer;
  the structure's captions go under the brackets, and rows space neighbours from the brackets.
- Bookkeeping: `S.validate(step_k, step_k1)` replays the arrows and compares bonds, charges and radicals (fix the
  chemistry until `validate: OK`); `S.balance(before, after)` flags atoms gained or lost that no arrow text mentions.

## Arrows, text and layout

- `S.arrow(above=None, below=None, kind="->", length=None, step=None, dashed=False, vertical=False)`: kinds `->`,
  `<-`, `<=>` (equilibrium; `<=>>`, `<<=>` unequal), `<->` (resonance), `=>`, `<=` (retrosynthesis), `-/->` and
  `<-/-` (no reaction), `->>` (several steps);
  `dashed=True` for a step that does not happen or is only proposed; any kind can be vertical. above/below take a string ("\n" between lines), a
  list of lines, a Mol (a reagent drawn as a structure) or a list mixing them: `above=[catalyst, "(5 mol%)"]`; the arrow
  grows to fit, vertical ones too. `step=n` puts a circled step number on the arrow. `S.branch(arrow, label,
  kind="out")` curves a by-product (text or Mol) out of a horizontal arrow, `kind="in"` a co-reagent into it, on the
  side its text leaves free (`side=`, `at=` to choose). A vertical arrow sets what it
  carries to its right (a structure, or one side alone, centred on it); `split=True` puts above= on its left and
  below= on its right (a reagent structure on one side, the base and solvent on the other).
- Text typography is automatic, inside <i>, <b> and beside markup too: H2O → H₂O, CH2Cl2, Pd(PPh3)4, LnPd(0), SN2 →
  S_N2, pKa′, N3- (azide), Et3N.HCl, ⊕ ⊖ as circled charges (raised when glued to a formula: H3O⊕), a colon against
  a symbol (B:, :NEt3) as lone-pair dots; numbers, temperatures, equiv. and percentages stay as written. `\word`
  keeps a word literal. Symbols, not words: Δ, hν, −78 °C, ‡, δ+.
- `S.caption(mol, text, color=None)` under a structure; `S.caption([m1, m2], text)` one caption under several.
  Declare before layout. Several captions on one structure stack in the order written (a bold number, then the yield);
  a caption under several goes below their own; neighbouring captions in a row keep 1 BU apart. `side="above"` sets a
  heading over the structure; `highlight="#FFFC79"` a fill behind a caption (also `S.text(..., highlight=)`).
  `S.frame(items, label="A")` draws a box round laid-out items (an equation on an exam page); `S.frame(size=(w, h),
  at=(x, y), label="B")` an empty answer box. `S.text(s, at)` free text; `S.plus()`; `S.ring_text(ring_atoms, "A")` a letter in a ring.
  A leading multiple in conditions ("2xLDA") prints as 2×.
- `S.row([...], gap=None)`: left to right, structures centred on the arrow line. `S.row(items, under=earlier_items,
  align="center" | "left" | "right", row_gap=0.9)`: a second row of products beneath the first. A scheme that turns a
  corner and snakes back: `r1 = S.row([...])`, then `S.row([product, S.arrow(kind="<-"), intermediate], under=r1,
  via=S.arrow(vertical=True))` (the vertical arrow runs down from r1's last item). `S.flow([...], width_cm=16)` wraps
  at arrows. `S.place(item, at)`.
- Cycles: `cyc = S.cycle([(species, deg), ...], radius=4)`, `S.branch(cyc, i, "label", kind="in" | "out", at=)`; the
  label may be a structure (a Mol drawn at the branch's end). `kind="step"` names arc i's step beside it ("oxidative
  addition"); `kind="feed"` is a straight arrow from outside into station i (`note="activation"` beside it).
- `S.orbital(atom, "p" | "s" | "lobe", angle, shade=, electrons=, size=, color=, label=)`: `size=1.6` for the larger
  coefficient of a polarised pi bond, `color=` for the outline (a shaded lobe takes a light tint of it). `shade="up"`
  shades the lobe along angle (an s sphere: shaded), `"down"` the other (an s sphere: white), so phases show; a
  `"lobe"` (hybrid) has its small back lobe; `label="n_{N(sp3)}"` names the orbital past its tip. Lobes are soft
  teardrops meeting at a bare atom's nucleus; on a labelled atom they come out from behind the letter, open there. An
  orbital's electrons sit as two dots in the full part of its lobe. A lone pair `a.lp()` on an atom with an empty lobe
  sits in it; a charge may stand off beside its own atom's lobe (C⊕ by its empty p orbital). A pair of p orbitals is the pi bond,
  so draw the sigma bond alone (check() counts it); a sigma frame drawn with wedges round an orbital atom is read as
  perspective, not stereochemistry.

## Teaching marks (lecture and textbook figures)

- `S.dot(atom, color="#E8261A", label="Nu")`: a site dot on a bare vertex, drawn over the bond ends, its label in the
  widest free gap (red nucleophilic, blue `#1F4FD1` electrophilic). It belongs to its step: `copy()` leaves it behind.
  On a labelled atom (an H marking a C-H position) it is a halo, a disc behind the label: `color="#D9D9D9"`.
- `S.dipole((a, b))`: the crossed dipole arrow beside bond a-b, pointing to its negative end b (`side=` to choose).
- `S.pointer(targets, text, side="below" | "above" | "left" | "right", color=, dashed=True, length=1.2)`: straight
  arrows from one note to one or more atoms or bonds ("2 possible sites for enolate formation"); the note stands
  beyond the molecule. `side=deg, bend=-60` instead draws a curved callout to one target ("C–C bond formed").
  `check()` reports a pointer that runs into a bond or a label: give it another side. Declare pointers, like captions,
  before the row() that places their structure: the row keeps their notes clear (a note above an arrow's conditions
  may reach over the arrow), and the structure's captions go below them. A later pointer's note moves aside from an
  earlier one's; a bond lying along the pointer is met from its side.
- `S.wiggle((a, b), color=)`: the wiggly line across a bond, the disconnection of a retrosynthesis; `at=1` crosses b's
  free end (an attachment point of a fragment); `dashed=True` a straight dashed line (a bond broken).
- Highlights: `bond.color`, `m.color(c, atoms)`. Teaching colours: red #E8261A, blue #1F4FD1,
  green #2E7D32.

## Output and checks

- `S.save("fig.png")` → (width_pt, height_pt), 1200 ppi of true size; `S.save("fig.svg")` only when asked;
  `S.save_pdf("fig.pdf")` print-ready A4. Insert at natural size (width = w_pt / 72 in).
- `S.check()`, `S.validate()`, `S.balance()`, `Scene.continuity()`, `stereo_report()`; `S.review("review.png")`
  and `S.zoom()` write close-ups for the visual audit.
- `side_by_side(original, "fig.png", "compare.png", crop=(x0, y0, x1, y1))`: a sketch, photo or slide beside the
  redraw at one height (crop in fractions of the original) for the element-by-element comparison.
- `Style(preset="acs")` strict ACS 0.6 pt lines; `Style(charge="plain")`; `Style(chemtext=False)` turns text
  typography off.
