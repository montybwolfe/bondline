"""Regression suite for bondline.

    python3 tests/run.py            run every case, compare with tests/expected, print PASS/FAIL and the total time
    python3 tests/run.py -k cycle   only cases whose name contains "cycle"
    python3 tests/run.py --twice    render every case twice and require byte-identical PNGs (determinism)
    python3 tests/run.py --keep D   keep the rendered files in D for inspection
    python3 tests/run.py --update   accept the current output as expected (only after reviewing every change)

On Linux the cases run in processes forked from one that has already imported the engine and matplotlib (each case
still gets a fresh process; this saves the half-second import per case); elsewhere, or with --spawn, each case starts a
new interpreter. Both run with PYTHONHASHSEED=0.

A case is tests/cases/<name>.py: an ordinary figure script (`from bondline import *`) that saves PNGs into the
working directory and prints what check(), validate(), continuity() and stereo_report() print. The expected output
is tests/expected/<name>.txt (exact text) plus a SHA-256 of every PNG's decoded pixels (tests/expected/hashes.json).
Pixels depend on the rendering stack (matplotlib, FreeType, font file), so hashes are compared only where
tests/expected/envfp.json matches this machine; hashing pixels rather than file bytes keeps the PNG encoder (Pillow and
its zlib) out of it. The printed check/validate output is always compared.
"""
import sys, os, json, time, hashlib, subprocess, tempfile, shutil, argparse, difflib, platform
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
ENGINE = os.path.normpath(os.path.join(HERE, "..", "scripts"))
CASES, EXP = os.path.join(HERE, "cases"), os.path.join(HERE, "expected")
BOOT = ("import sys, os, runpy; sys.path.insert(0, sys.argv[1]); os.chdir(sys.argv[3]); "
        "import matplotlib; matplotlib.use('Agg'); sys.argv = [sys.argv[2]]; runpy.run_path(sys.argv[0], run_name='__main__')")

def envfp():
    import matplotlib, numpy
    from matplotlib import ft2font
    from matplotlib.font_manager import findfont, FontProperties
    f = findfont(FontProperties(family="Liberation Sans"), fallback_to_default=True)
    return {"python": platform.python_version(), "matplotlib": matplotlib.__version__, "numpy": numpy.__version__,
            "freetype": ft2font.__freetype_version__, "font": os.path.basename(f),
            "font_sha": hashlib.sha256(open(f, "rb").read()).hexdigest()[:16]}

def sha(p):
    """SHA-256 of a PNG's decoded pixels (mode, size, pixel bytes): the image, whichever zlib compressed the file."""
    from PIL import Image
    with Image.open(p) as im:
        h = hashlib.sha256(f"{im.mode} {im.size[0]}x{im.size[1]}\n".encode()); h.update(im.tobytes())
    return h.hexdigest()

def render(case, out):
    os.makedirs(out, exist_ok=True)
    env = dict(os.environ, PYTHONHASHSEED="0")
    t = time.time()
    p = subprocess.run([sys.executable, "-c", BOOT, ENGINE, os.path.join(CASES, case + ".py"), out],
                       capture_output=True, text=True, env=env)
    txt = p.stdout + (("\n[exit %d]\n" % p.returncode + p.stderr[-3000:]) if p.returncode else "")
    open(os.path.join(out, "stdout.txt"), "w").write(txt)
    files = {f: sha(os.path.join(out, f)) for f in sorted(os.listdir(out)) if f.endswith(".png")}
    return txt, files, time.time() - t

def render_forked(cases, j):
    """Run each (case, out) in a child forked from this process, which has imported the engine already: as isolated as a
    new interpreter, without the import cost. At most j children at a time. Returns {out: (txt, files, seconds)}."""
    import runpy, traceback
    sys.path.insert(0, ENGINE)
    import matplotlib; matplotlib.use("Agg")
    import bondline  # noqa: F401  (imported once here, inherited by every child)
    todo = list(cases); running = {}; res = {}
    def start(case, out):
        os.makedirs(out, exist_ok=True); fo, fe = os.path.join(out, "stdout.txt"), os.path.join(out, "stderr.txt")
        sys.stdout.flush(); sys.stderr.flush(); pid = os.fork()
        if pid == 0:                                       # child: run the case with its output in files
            code = 0
            try:
                o = os.open(fo, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644); e = os.open(fe, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)
                os.dup2(o, 1); os.dup2(e, 2); os.chdir(out); sys.argv = [os.path.join(CASES, case + ".py")]
                runpy.run_path(sys.argv[0], run_name="__main__")
            except SystemExit as x: code = x.code if isinstance(x.code, int) else (0 if x.code is None else 1)
            except BaseException: traceback.print_exc(); code = 1
            finally:
                try: sys.stdout.flush(); sys.stderr.flush()
                finally: os._exit(code)
        running[pid] = (case, out, time.time())
    while todo or running:
        while todo and len(running) < j: start(*todo.pop(0))
        pid, st = os.waitpid(-1, 0); case, out, t = running.pop(pid); code = os.waitstatus_to_exitcode(st)
        txt = open(os.path.join(out, "stdout.txt")).read(); err = open(os.path.join(out, "stderr.txt")).read()
        os.remove(os.path.join(out, "stderr.txt"))
        if code: txt += "\n[exit %d]\n" % code + err[-3000:]
        open(os.path.join(out, "stdout.txt"), "w").write(txt)
        res[out] = (txt, {f: sha(os.path.join(out, f)) for f in sorted(os.listdir(out)) if f.endswith(".png")}, time.time() - t)
    return res

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("-k", default=""); ap.add_argument("-j", type=int, default=max(2, os.cpu_count() or 2))
    ap.add_argument("--twice", action="store_true"); ap.add_argument("--keep"); ap.add_argument("--update", action="store_true")
    ap.add_argument("--spawn", action="store_true", help="start a new interpreter for every case instead of forking")
    ap.add_argument("-v", action="store_true", help="show diffs of printed output")
    a = ap.parse_args()
    fork = hasattr(os, "fork") and sys.platform.startswith("linux") and not a.spawn
    if fork and os.environ.get("PYTHONHASHSEED") != "0":   # the forked cases inherit this process's hash seed
        os.execve(sys.executable, [sys.executable] + sys.argv, dict(os.environ, PYTHONHASHSEED="0"))
    names = sorted(f[:-3] for f in os.listdir(CASES) if f.endswith(".py") and a.k in f)
    hp = os.path.join(EXP, "hashes.json"); hashes = json.load(open(hp)) if os.path.exists(hp) else {}
    fp = envfp(); fpp = os.path.join(EXP, "envfp.json")
    same_stack = os.path.exists(fpp) and json.load(open(fpp)) == fp
    root = a.keep or tempfile.mkdtemp(prefix="cdtests_"); t0 = time.time()
    if fork:
        jobs = [(n, os.path.join(root, n)) for n in names] + ([(n, os.path.join(root, n + "__2")) for n in names] if a.twice else [])
        R = render_forked(jobs, a.j)
        results = [(n, R[os.path.join(root, n)], R.get(os.path.join(root, n + "__2"))) for n in names]
    else:
        def job(n):
            r1 = render(n, os.path.join(root, n))
            r2 = render(n, os.path.join(root, n + "__2")) if a.twice else None
            return n, r1, r2
        with ThreadPoolExecutor(a.j) as ex: results = list(ex.map(job, names))
    fails = 0
    for n, (txt, files, dt), r2 in results:
        probs = []
        ep = os.path.join(EXP, n + ".txt")
        if a.update:
            open(ep, "w").write(txt); hashes[n] = files
        else:
            exp = open(ep).read() if os.path.exists(ep) else None
            if exp is None: probs.append("no expected output (new case: review it, then --update)")
            elif exp != txt:
                probs.append("printed output differs")
                if a.v: probs += ["   " + l for l in difflib.unified_diff(exp.splitlines(), txt.splitlines(), "expected", "now", lineterm="", n=0)]
            if same_stack and n in hashes and hashes[n] != files:
                ch = sorted(f for f in set(files) | set(hashes[n]) if files.get(f) != hashes[n].get(f))
                probs.append("PNG differs: " + ", ".join(ch))
        if r2 is not None and (r2[1] != files or r2[0] != txt): probs.append("not deterministic (second render differs)")
        if "[exit" in txt: probs.append("script failed")
        fails += bool(probs)
        print(f"{'FAIL' if probs else 'PASS'}  {n:38s} {dt:5.1f}s")
        for p in probs: print("      " + p)
    if a.update:
        json.dump(hashes, open(hp, "w"), indent=1, sort_keys=True); json.dump(fp, open(fpp, "w"), indent=1)
    if not a.keep: shutil.rmtree(root, ignore_errors=True)
    note = "" if same_stack or a.update else "  (PNG hashes not compared: different rendering stack)"
    print(f"{len(names) - fails}/{len(names)} passed in {time.time() - t0:.1f} s{note}")
    return 1 if fails else 0

if __name__ == "__main__": sys.exit(main())
