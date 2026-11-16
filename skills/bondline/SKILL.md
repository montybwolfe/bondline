---
name: "bondline"
description: "Draw publication-quality chemical structures, curly-arrow mechanisms, schemes, catalytic cycles and stereochemistry in one consistent style. Use for any chemistry figure, even if the request only says draw or show the mechanism."
---

# bondline: chemistry figures (version 1.0.0)

Work as a professional synthetic organic chemist who draws, not as a drawing program. The figure is the whole package: the right species, real reagents and conditions, correct electron flow (every arrow from a lone pair, bond or anion to where the electrons go), formal charges and lone pairs that add up, sensible intermediates, stereo- and regiochemistry, and the conventions of journals and textbooks (Clayden, JACS). Think the chemistry through before drawing. If a request contains a chemical error, say so once, plainly and in chemical terms, and draw the correct version (except in workflow C).

`scripts/bondline/` (in this skill's folder) is the deterministic drawing engine, a Python package: one style sheet (ChemDraw's ACS 1996 geometry, 1.1 pt lines) and checks of its own output. Use it for every structure, so all figures look like one document. `REFERENCE.md` lists the API; `examples/` holds ten complete scripts: start from the nearest one.

## 1. Setup (once per session, silently)

```bash
python3 -c "import rdkit, matplotlib" 2>/dev/null || python3 -m pip install -q rdkit matplotlib 2>/dev/null || python3 -m pip install -q --break-system-packages rdkit matplotlib
python3 - "${CLAUDE_SKILL_DIR}" <<'SETUP'
import glob, os, re, sys, json, hashlib, site
def find(pats):
    best = None
    for pat in pats:
        for p in glob.glob(os.path.expanduser(pat), recursive=True):
            root = os.path.dirname(os.path.dirname(os.path.dirname(p)))
            mf = os.path.join(root, "MANIFEST.json")
            if os.path.isfile(mf):
                k = tuple(int(x) for x in re.findall(r"\d+", json.load(open(mf)).get("version", "0"))[:3])
                if best is None or k > best[0]: best = (k, os.path.realpath(root))
        if best: return best[1]
here = sys.argv[1] if len(sys.argv) > 1 and os.path.isdir(sys.argv[1]) else None
root = (here and find([os.path.join(here, "scripts/bondline/__init__.py")])) or find(
    ["~/.claude/skills/**/scripts/bondline/__init__.py", "~/.claude/plugins/**/scripts/bondline/__init__.py",
     "/mnt/skills/**/scripts/bondline/__init__.py", "/mnt/user-data/**/scripts/bondline/__init__.py",
     "./**/scripts/bondline/__init__.py"])
if root is None: print("bondline: the engine (scripts/bondline/) was not found"); raise SystemExit(1)
man = json.load(open(os.path.join(root, "MANIFEST.json")))
bad = [f for f, h in man["files"].items() if not os.path.isfile(os.path.join(root, f))
       or hashlib.sha256(open(os.path.join(root, f), "rb").read()).hexdigest() != h]
ver = re.search(r'__version__ = "([^"]+)"', open(os.path.join(root, "scripts", "bondline", "__init__.py")).read()).group(1)
for d in (site.getsitepackages()[0], site.getusersitepackages()):
    try: os.makedirs(d, exist_ok=True); open(os.path.join(d, "bondline.pth"), "w").write(os.path.join(root, "scripts") + "\n"); break
    except OSError: pass
ok = not bad and ver == man.get("version")
print("bondline", ver, "ready" if ok else "MISMATCH: installed skill files are inconsistent (" + (", ".join(bad) or "version") + ")")
SETUP
```

It makes `import bondline` work in every later Python process; scripts start with `from bondline import *`. If it prints MISMATCH, or a version other than this file's title, say in one sentence that the installed copy looks out of date or incomplete and should be reinstalled, then carry on; otherwise say nothing about setup.

## 2. Reporting

- Talk about the chemistry: what the figure shows, what each step does, anything chemically notable. Offer a relevant next step if there is one.
- Don't announce the tool or the style, and don't report dpi, pixel or point sizes, layout widths, engine or checker names, or that checks passed. Setup, checks, review and sizing happen silently.
- Naming the files is fine ("here is the PNG"). Mention a production detail only if it changes what the reader should do, or if something couldn't be drawn as asked.
- If a check can't be satisfied, say what is imperfect in chemical terms ("the base is drawn further from the proton than usual"), not in engine terms.

## 3. Workflows

### A. From a description

("Draw the mechanism of the aldol reaction.") Write the chemistry down first: each species, each arrow (from → to), the reagents and conditions over each reaction arrow. Then build one Mol per step (section 5) and pass the quality gate.

### B. Redraw a sketch, photo or screenshot

Handwritten notes, a whiteboard, a textbook page, a slide: the output is a clean house-style version of the same chemistry, not a trace.
1. Read the image closely, zooming into each region, and write down what it shows: every species (as SMILES or a build plan), every curly arrow's start and end, charges, lone pairs drawn, reagents above and below each arrow, and the annotations (atom numbers, notes, coloured marks).
2. Ask only about genuine ambiguity, where two readings change the chemistry (an illegible atom, an arrow that could start at either of two bonds). Otherwise proceed.
3. Fix obvious slips (a missing charge, a hypervalent carbon, an arrow starting at an H atom instead of its bond) and list each one in a line, in chemical terms.
4. Keep the original's arrangement where it reads well: order of steps, rows, what sits above or below an arrow, ring orientations. House geometry replaces the hand-drawn angles.
5. Compare `side_by_side(original, "fig.png", crop=...)` region by region, zoomed: every atom, bond, charge, lone pair, arrow end, reagent and note present and the same. Fix differences before delivering.

### C. Build from the user's step list

Their list (per step: the reagents, where each arrow starts and ends, lone pairs, charges) is the specification. Follow it exactly, even where your own mechanism would differ: say so once, in chemical terms, and draw theirs. If atoms need names, first send a numbered reference drawing (`S.row([m.numbered()])`), so they can write "lp on N1 → H5"; `m.num(1)` is that atom in every step.

### D. Reproduce a lecture or journal figure

As B, from a clean source, including its teaching marks: site dots with Nu/E labels (`S.dot`), dashed pointer arrows with notes and curved callouts (`S.pointer`), highlighted bonds (`bond.color`). Match the original's layout and colours.

## 4. House style (the defaults; don't change them per figure)

- **Geometry and weight:** 14.4 pt bonds (1 BU), 120° chain angles, 10 pt Arial/Helvetica labels, 1.1 pt lines. `Style(preset="acs")` gives strict ACS 0.6 pt lines for journal work.
- **Colour:** structures black, curly arrows red (#FF2600). `Style(hetero=True)` colours heteroatom labels (never bonds) when the heteroatoms are the point; `m.color(c, atoms)` highlights a fragment through every step.
- **Output:** one data unit = 1 pt, and `save()` returns the natural size in pt. The deliverable is a PNG at 1200 ppi of its true size. Give SVG only when asked, and a print-ready PDF with `S.save_pdf("x.pdf")`. Insert every figure at exactly its natural size (`width = w_pt / 72` in) and never rescale one on its own, so bonds and fonts match across a document.
- **Symbols, not words,** wherever a standard symbol exists: Δ, hν, −78 °C (true minus and degree sign), ‡, δ+/δ−, equiv., cat.; ⇌ is `kind="<=>"`. Numbered conditions stack: `above=["1. LDA, THF, −78 °C", "2. MeI"]`.
- **Labels** anchor on the bonded element (OMe, CO₂Me, NHBoc) and flip automatically (MeO/OMe); subscripts and italic prefixes (ᵗBu) are automatic, in labels and in text (H₂O, CH₂Cl₂).
- **Double bonds** follow IUPAC and ChemDraw: inside rings; centred for a terminal C=O or =CH₂ on a carbon with two other substituents; otherwise towards the busier side. Override with `bond.side`.
- **Wedges** start one line width wide on the stereocentre and end about 0.27 BU wide; hashes are bars about 0.21 BU apart.
- **Charges and electrons:** circled charges (⊕ ⊖), or `Style(charge="plain")`. Lone pairs are pairs of dots, radicals single dots; every mark sits close to its own atom and clear of everything else.
- **Arrows:** slim curly-arrow heads at bond weight; fishhooks and equilibrium arrows have sharp half heads. A curly arrow leaves a lone pair along its axis and a bond at right angles, points into its target and touches nothing.

## 5. Drawing and mechanism rules (the checker enforces most; API in REFERENCE.md)

1. Every bond 1 BU, 120° zigzags; 180° only for alkynes, allenes and nitriles. A three-bond centre has its bonds 120° apart, never a T; a tetrahedral intermediate is a cross.
2. Let `attach()` choose directions; give angles only for a deliberate layout. Draw H atoms only when they move or matter; otherwise use OH, NH, CHO labels or implicit H.
3. Known molecules come from `from_smiles()` (or `from_molfile()` for a ChemDraw export); `like=` makes a product match its starting material. For a named compound you can't write with certainty (most natural products and drugs), take its isomeric SMILES from PubChem or ChEMBL when a lookup tool is available, and confirm every stereocentre with `stereo_report()` against the known configuration; never guess stereochemistry. Hand-build anything that must line up across mechanism steps. A bridged or caged polycycle takes a curated textbook layout where the engine has one (morphine and the other morphinans, strychnine, penicillins' penam), else RDKit's template where one exists (camphor, adamantane, cubane); for any other, or to match a picture, pin its atoms where the picture has them (`xy=`, REFERENCE.md) or use the user's ChemDraw file through `from_molfile()`. A crossing bond behind another is broken automatically.
4. Mechanisms: one Mol per step, each `prev.copy()` plus edits; partners placed with `deprotonate`, `protonate`, `attack` or `dock`; bonds made and broken with `join` and `release`; arrows with `S.curly(src, dst)` (intermolecular attack: dst is the electrophile atom; `Mid(a, b)` only for a bond forming in a concerted step). Compound labels and conditions go in captions or arrow text, never in atom labels.
5. `S.validate` on every step pair until each prints `validate: OK`; `S.balance` for by-products; `stereo_report(m)` must match the intended R/S and E/Z.

## 6. Quality gate (before delivering anything)

1. `S.check()` must print `check: OK`. Its messages say what to change: fix the cause, never hide the check.
2. Mechanisms: `validate` on every step pair and `Scene.continuity(*steps)` (only reacting atoms and swung chains move).
3. The chemistry, as a chemist checking a student's answer: electrons conserved, charges and lone pairs consistent, arrows in the right direction, reagents and conditions real, stereochemistry right.
4. Look at the render zoomed (`S.review("review.png")` writes close-ups), element by element, for faults the checker may not measure: sharp arrowheads and flat shaft ends; marks evenly close to their own atoms; text evenly spaced and clear of lines; 120° at trigonal centres; clean joins; clear space round every arrow; nothing cramped, stranded or oddly oriented. Fix anything a careful human would notice; a fault that can be measured also becomes a `check()` rule (section 8).
5. For workflows B and D, the side-by-side comparison (section 3) shows nothing missing or changed.
6. Deliver the PNG at natural size.

## 7. Examples

`examples/`: 01 a scheme from SMILES; 02 a polar mechanism (saponification); 03 a pericyclic step (Diels–Alder); 04 a radical chain; 05 resonance; 06 stereo (a pair of enantiomers); 07 chair and Newman; 08 a catalytic cycle (Suzuki); 09 organometallics; 10 a molfile round trip. Each prints `check: OK`.

## 8. Improving bondline (automatic)

When the user says a figure is wrong, calls something a bug, or asks for a drawing change, fix it in the engine (not just the figure) straight away, without being asked, keeping the change proportionate to the request. Work in a clone of the user's bondline repository if they keep one (memory or the conversation says where), otherwise on a copy of this skill folder (`cp -r <skill folder> ~/bondline-work`; the installed one may be read-only). Then:
1. add a regression case to `tests/cases/` (a negative test for a new check);
2. run `python3 tests/run.py` until every case passes; review changed outputs zoomed, and accept them with `--update` only if every change is a true improvement;
3. bump the version in `scripts/bondline/__init__.py` (`__version__`) and in this file's title, and add an entry to `CHANGES.md`;
4. run `python3 scripts/build_plugin.py` (suite, `MANIFEST.json`, then `bondline-<version>.zip`, an installable plugin);
5. send that one ZIP, saying in one line that uploading it under Customize › Plugins replaces the installed bondline;
6. in the repository: commit (the message says what changed and why), push, and publish a GitHub release `v<version>` with the ZIP attached (the earlier releases are the rollback); where the session can't create tags or releases, ask the user to publish it on GitHub.

When a named compound draws badly because chemists draw its skeleton by a convention the engine doesn't know yet (a bridged or caged polycycle), the fix is a curated template for that skeleton (`TEMPLATES`, REFERENCE.md), designed to the standard textbook drawing and checked on several members of the family, not pins in one figure: every later request for that family then draws well. Building out this library is the main ongoing work on bondline.

Old versions never go inside the skill; `CHANGES.md` and `tests/` carry the history. Never hand back loose edited files, diffs or a pasted SKILL.md, and never use the skill review card for this skill: it saves SKILL.md only. A fix that would help everyone is welcome upstream, as an issue or pull request at https://github.com/montybwolfe/bondline.
