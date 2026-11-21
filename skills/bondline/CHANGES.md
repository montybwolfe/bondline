# bondline changes

## 1.1.0

### Fixes

- Bonds between two labels (O–H, H–Cl) show the same length everywhere. An O–H that a curly arrow starts or ends on
  used to keep its short default stub, so the two O–H bonds of one water could differ.
- At a centred double bond (a ketone's C=O, an isopropenyl's C=CH₂), the neighbouring bonds stop at its near line,
  even when one of them is almost in line with it. Previously a C=CH₂ didn't trim them at all and a wide-angled C=O
  didn't either, so a stub of the neighbour poked into the gap between the two lines.
- A double bond beside the wide end of a hash is drawn offset, as beside a wedge, so its line no longer runs into the
  last hash bar.
- The + and − of text (H⁺, Na⁺, superscript charges) carry the same bar weight as the signs in circled charges. The
  font's own superscript plus was about half the weight of a letter's stem. A sign followed by a space, the minus of
  "a − b" in an exponent, is an operator and keeps the font's weight.
- Text in PNGs is drawn at its true weight everywhere. The rasteriser used to round the square stems of H, L, I, E
  and T to a whole number of pixels, so the same letter came out 3 or 4 pixels wide depending on where it stood (two
  L ligands, a bonded H beside the H of H₃PO₄). Each run of text also starts on a whole pixel, so a letter is drawn
  with identical pixels wherever it stands.
- An orbital lobe drawn on a labelled atom (the p orbital of a C=O's oxygen) is shaded right up to the letter. A thin
  white line used to separate the shading from the letter. An orbital's name (`label=`) stands 5 pt clear of its
  lobe instead of 1.5 pt.
- Rows of a scheme (`flow()`, `row(under=)`) sit 1.3 BU apart instead of 0.9, so captions no longer crowd the
  conditions of the row below.

### Curly arrows

- When no single or joined arc can reach its target without touching or crossing a bond or label, the arrow is drawn
  as a free curve instead: a smooth one-way bend that leaves its source over a wider fan of directions and swings
  round whatever is in the way. Arrows that fit as arcs are unchanged.

### Examples

- The Suzuki cycle (example 08 and the README) draws its real aryl groups (4-bromoanisole and phenylboronic acid to
  4-methoxybiphenyl) on Pd(PPh₃)₂, with the phenyl cis to the aryl for reductive elimination, and names each step.
- New README figures: ferrocene and acetylferrocene, and a gallery of orbitals, a chair, a Newman projection and
  (η⁶-benzene)Cr(CO)₃.
- The skill says plainly that molfiles go in with `from_molfile()` and out with `to_molfile()`.

### Style decisions

- Signs in text are drawn at the charge marks' weight (0.79 pt by default) rather than the font's, so a superscript
  plus is a little heavier than in ChemDraw.

## 1.0.0

First release.

### Style decisions

House choices that differ from ChemDraw's defaults or common practice; each could become an option.

- Circled charges are as tall as a capital letter, with a lighter ring, because ChemDraw's larger ones crowd letters.
- Greek letters (Δ, δ, ν) come from a Times-style serif, as journals set them.
- A superscript slides up to 2 pt aside for a steep bond, rather than the bond starting beyond it.
- Lopsided and narrow letters (L, F, P, I, T) meet bonds at their ink, not their box.
- An aldehyde's H is left implicit unless you draw it, to keep schemes uncluttered.
- Numbered groups are superscripts (R¹) and counts subscripts (R₂), so the two never look alike.
- Ring-fusion stereocentres show a wedged H (the steroid convention) rather than a wedged ring bond.
- Perspective cages use plain bonds with the back bond broken, keeping wedges for stereochemistry.
- Decorative wedges on non-stereocentres (a gem-dimethyl) are drawn plain.
- Electrons in orbitals are dots, not half-arrows, and orbitals are grey unless given `color=`.
- Arrowheads are ChemDraw's default size, a little larger against the bonds than in most textbooks; it shows most on
  the double-headed resonance arrow.
- H⁺ and H₃O⁺ in conditions take a plain superscript plus; type ⊕ for the circled one.
- The two fishhooks of a homolysis sit on opposite sides of the bond, rather than mirrored on one side.
- Spectator ions are separate ions with circled charges, the charge after the last letter.
- Conditions are centred on their arrow, and captions under their structure.
- `abbrev="auto"` contracts only what schemes always write as labels, larger groups first, unless `labels=` says
  otherwise.
- Outside `like=`, the bond RDKit picks to wedge is kept (move it with `m.stereo()`).

### Known limitations

- Octahedral tris-chelates (Ru(bpy)₃²⁺, fac-Ir(ppy)₃) and other chelate complexes come out tangled from SMILES. Pin
  every atom with `xy=`, or import a molfile.
- Taxanes and tropanes have no curated template yet and come out crowded, and in 14-hydroxymorphinans (oxycodone,
  naloxone) the 14-OH crowds the bridge. Pin the skeleton with `xy=`, or import a molfile.
- Macrocycles come out as regular polygons, and the alkynes round a dehydroannulene bend; pin the ring with `xy=`.
- Very crowded centres (an aryl-laden quaternary stereocentre, a bulky biaryl phosphine, a carbon shared by three
  rings) can give uneven bonds or tight angles, which `check()` reports. Pinning the atoms round the centre fixes it.
- `like=` doesn't fold an open chain onto its product's ring, so a ring-closing metathesis substrate stays extended
  unless you pin the chain.
- There is no fan-out from one starting material to several products, no square layout closed by a second vertical
  arrow, and no elbow arrow under a row. Use one row per product, or place items with `S.place()`.
- Two catalytic cycles that share a step, or a reversible pair of arcs between two stations, can't be joined; draw
  the two cycles side by side with `center=`.
- There are no orbitals along a bond (σ ellipses, σ* lobes) and no wavy interaction lines; p orbitals are the nearest
  stand-in.
- Reaction arrows can't be coloured (their text can), highlights go behind captions and free text only, there are no
  braces, hollow block arrows or angle arcs, and ≡ between structures is text-sized.
- Polymer brackets take their n, but the backbone bonds stop at the brackets instead of crossing them.
- In very crowded places a curly arrow can still go the long way round or touch a bond, and the placing helpers can
  struggle with bulky partners or a hypervalent centre. `check()` says so; steer with `side=` or `bend=`, or place the
  partner with `dock()`.
- After an SN2, `join()` lays the new bond in the widest gap rather than along the line of attack; turn it with
  `swing()`.
- Only a curved callout's note moves itself to keep 5 pt from structures; notes on a named side, captions and free
  text are checked and reported, not moved.
- A ring bond between two labels (N–Pd in a palladacycle) can stay short; neighbouring labels that collide (CO₂H
  beside an N–H) need their side set with `a.ldir`; an H-bond acceptor between two donors is placed by hand.
