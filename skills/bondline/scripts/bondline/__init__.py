"""bondline: consistent, coordinate-based structure and mechanism drawings.

Molecule coordinates are in bond lengths (BU). Everything renders at one fixed physical
scale set by Style.bond (pt per bond), so every figure drawn with the same Style matches.
"""
__version__ = "1.0.0"
from ._core import *  # noqa: F401,F403
from ._scene_build import _SceneBuild
from ._scene_arrows import _SceneArrows
from ._scene_finish import _SceneFinish


class Scene(_SceneBuild, _SceneArrows, _SceneFinish):
    """A figure: structures, arrows, marks and text, laid out and drawn."""


def _share():                                     # one namespace across the package, as in a single module
    from . import _core, _scene_build, _scene_arrows, _scene_finish
    mods = (_core, _scene_build, _scene_arrows, _scene_finish)
    names = {k: v for k, v in globals().items() if not k.startswith("__")}
    for m in mods:
        for k, v in names.items(): vars(m).setdefault(k, v)
_share()
