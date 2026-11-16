# bondline changes

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
