"""Build the bondline plugin ZIP: the one deliverable of any change to the skill.

    python3 scripts/build_plugin.py [--out DIR] [--skip-tests]

1. checks that SKILL.md's title line and scripts/bondline/__init__.py (__version__) carry the same version;
2. runs the regression suite (tests/run.py) and stops if any case fails;
3. writes MANIFEST.json (version + SHA-256 of every file of the skill, which the setup step verifies);
4. writes DIR/bondline-<version>.zip, an installable plugin: .claude-plugin/plugin.json, README.md, LICENSE and
   skills/bondline/ (sorted entries, fixed timestamps, so the same files always give the same ZIP). In the plugin's
   repository these come from the repository root, and the manifest's version is brought up to date; in an installed
   copy, where only the skill folder exists, they are generated;
5. re-opens the ZIP and checks every file against the manifest.
DIR defaults to the folder that holds the repository (or the skill folder), so the ZIP never lands inside it.
"""
import os, re, sys, json, hashlib, zipfile, subprocess, argparse

SKILL = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
REPO = os.path.normpath(os.path.join(SKILL, "..", ".."))
IN_REPO = (os.path.basename(os.path.dirname(SKILL)) == "skills"
           and os.path.isfile(os.path.join(REPO, ".claude-plugin", "plugin.json")))
HOME = "https://github.com/montybwolfe/bondline"

def skill_files():
    out = []
    for d, ds, fs in os.walk(SKILL):
        ds[:] = sorted(x for x in ds if x != "__pycache__" and not x.startswith("."))
        for f in sorted(fs):
            rel = os.path.relpath(os.path.join(d, f), SKILL).replace(os.sep, "/")
            if not f.startswith(".") and not f.endswith((".pyc", ".zip")) and rel != "MANIFEST.json": out.append(rel)
    return out

def read(p): return open(p, "rb").read()

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out"); ap.add_argument("--skip-tests", action="store_true")
    a = ap.parse_args()
    v = re.search(r'__version__ = "([^"]+)"', open(os.path.join(SKILL, "scripts", "bondline", "__init__.py")).read()).group(1)
    skill_md = open(os.path.join(SKILL, "SKILL.md"), encoding="utf-8").read()
    m = re.search(r"^# .*\(version ([^)]+)\)", skill_md, re.M)
    if not m or m.group(1) != v:
        sys.exit(f"version mismatch: __init__.py says {v}, SKILL.md title says {m.group(1) if m else None}")
    if not a.skip_tests and subprocess.run([sys.executable, os.path.join(SKILL, "tests", "run.py")]).returncode:
        sys.exit("tests failed: fix them (or review and --update) before packaging")
    fl = skill_files()
    man = {"name": "bondline", "version": v, "files": {f: hashlib.sha256(read(os.path.join(SKILL, f))).hexdigest() for f in fl}}
    open(os.path.join(SKILL, "MANIFEST.json"), "w").write(json.dumps(man, indent=1, sort_keys=True) + "\n")
    if IN_REPO:                                            # the repository's own files, version brought up to date
        pj = os.path.join(REPO, ".claude-plugin", "plugin.json")
        plug = json.load(open(pj, encoding="utf-8")); plug["version"] = v
        open(pj, "w", encoding="utf-8").write(json.dumps(plug, indent=2, ensure_ascii=False) + "\n")
        readme, licence = read(os.path.join(REPO, "README.md")), read(os.path.join(REPO, "LICENSE"))
    else:                                                  # an installed copy: only the skill folder is here
        desc = re.search(r'^description: "(.*)"$', skill_md, re.M).group(1)
        plug = {"name": "bondline", "version": v, "description": desc, "author": {"name": "Monty Wolfe"},
                "homepage": HOME, "license": "MIT"}
        readme = (f"# bondline {v}\n\nA natural-language chemistry engine for Claude, rebuilt from an installed copy. "
                  f"Source, documentation and releases: {HOME}\n").encode()
        licence = read(os.path.join(SKILL, "LICENSE.txt"))
    entries = [(".claude-plugin/plugin.json", (json.dumps(plug, indent=2, ensure_ascii=False) + "\n").encode()),
               ("LICENSE", licence), ("README.md", readme)]
    entries += [("skills/bondline/" + f, read(os.path.join(SKILL, f))) for f in fl + ["MANIFEST.json"]]
    out = a.out or os.path.dirname(REPO if IN_REPO else SKILL); os.makedirs(out, exist_ok=True)
    zp = os.path.join(out, f"bondline-{v}.zip")
    with zipfile.ZipFile(zp, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as z:
        for name, data in entries:
            zi = zipfile.ZipInfo(name, date_time=(2026, 1, 1, 0, 0, 0))
            zi.external_attr = 0o644 << 16; zi.compress_type = zipfile.ZIP_DEFLATED
            z.writestr(zi, data)
    with zipfile.ZipFile(zp) as z:
        bad = [f for f, h in man["files"].items() if hashlib.sha256(z.read("skills/bondline/" + f)).hexdigest() != h]
        pj = json.loads(z.read(".claude-plugin/plugin.json"))
        if bad or pj.get("name") != "bondline" or pj.get("version") != v: sys.exit(f"ZIP check failed: {bad or pj}")
    print(f"{zp}: bondline {v}, {len(entries)} files, {os.path.getsize(zp) // 1024} kB, verified")

if __name__ == "__main__": main()
