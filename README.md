# bondline

<img src="https://raw.githubusercontent.com/montybwolfe/bondline/main/docs/images/icon.png" alt="" width="120" align="right">

A natural-language chemistry engine for Claude.

bondline lets Claude draw clean, consistent chemistry figures from a plain description: structures, curly-arrow mechanisms, reaction schemes and catalytic cycles. Every figure is checked before it's shown, including valences, charges, stereochemistry and the electron flow in mechanisms.

![Morphine, strychnine, reserpine, quinine, penicillin G and cholesterol](https://raw.githubusercontent.com/montybwolfe/bondline/main/docs/images/structures.png)

## Features

- Structures from a name, SMILES or molfile, including natural products and bridged ring systems
- Mechanisms with curly arrows, checked step by step
- Reaction schemes and catalytic cycles with reagents and conditions
- Orbitals, chairs, Newman projections, transition states and callouts
- Redraws figures from photos of notes or whiteboards
- Molfile import/export for ChemDraw
- PNG at 1200 ppi by default, SVG or PDF on request
- Deterministic output

## Usage

Once it's installed, ask Claude for the figure you want, e.g.

- "Draw morphine"
- "Show the mechanism of the Robinson annulation"
- "Draw the catalytic cycle for a Suzuki coupling"

## Examples

Takasago synthesis of (−)-menthol:

![Takasago synthesis of (−)-menthol](https://raw.githubusercontent.com/montybwolfe/bondline/main/docs/images/synthesis.png)

Robinson annulation to the Wieland–Miescher ketone:

![Robinson annulation mechanism](https://raw.githubusercontent.com/montybwolfe/bondline/main/docs/images/mechanism.png)

Suzuki coupling:

<img src="https://raw.githubusercontent.com/montybwolfe/bondline/main/docs/images/cycle.png" alt="Suzuki coupling catalytic cycle" width="420">

The scripts for these are in [docs/showcase](https://github.com/montybwolfe/bondline/tree/main/docs/showcase).

## Installation

### claude.ai and the desktop app

Download `bondline-<version>.zip` from [Releases](https://github.com/montybwolfe/bondline/releases) and upload it under Customize › Plugins › Add › Upload plugin. Uploading a newer zip replaces the installed version.

### Claude Code

```
/plugin install bondline --marketplace montybwolfe/bondline
```

On older versions of Claude Code:

```
/plugin marketplace add montybwolfe/bondline
/plugin install bondline@bondline
```

You can also clone the repo and run `claude --plugin-dir path/to/bondline`.

## Requirements

Python 3 with matplotlib and RDKit. Missing packages are installed with pip on first use, and a `bondline.pth` file is added to site-packages so `import bondline` works. Rendering runs offline. If Claude needs to confirm a structure, it may look the compound up on PubChem or ChEMBL.

## Tests

```
python3 skills/bondline/tests/run.py
```

This redraws about 200 reference figures and compares them with stored output. Release notes and known issues are in [CHANGES.md](https://github.com/montybwolfe/bondline/blob/main/skills/bondline/CHANGES.md).

## Licence

MIT
