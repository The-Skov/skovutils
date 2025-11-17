"""Top-level package for skovutils utilities."""

from importlib import metadata as _metadata

__author__ = "Mads Skov"

try:
    __version__ = _metadata.version("skovutils")
except _metadata.PackageNotFoundError:  # pragma: no cover - during editable installs
    __version__ = "0.0.0"

# Re-export frequently used helpers at the package root for backwards compatibility.
from .tree import *  # noqa: F401,F403
from .rm_correlation import *  # noqa: F401,F403
from .filtering import *  # noqa: F401,F403

__all__ = [
    "filtering",
    "rm_correlation",
    "tree",
]
