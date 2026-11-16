# Every script in examples/ runs as written (each prints its check, validate and stereo results; its figure is compared
# like any other). The examples are the models a figure is written from, so they must stay correct.
import os, runpy, glob
here = os.path.dirname(os.path.abspath(__file__))
for f in sorted(glob.glob(os.path.join(here, "..", "..", "examples", "*.py"))):
    print("==", os.path.basename(f)); runpy.run_path(f, run_name="__main__")
