try:
    from importlib.metadata import version as get_version
except ImportError:
    # for python < 3.8
    try:
        from importlib_metadata import version as get_version
    except ImportError:
        def get_version(package):
            return "unknown"

from ._reconstruct import reconstruct
from ._solve import solve, solve_m
from ._ut_constants import (
    constit_index_dict,
    cycles_per_hour,
    hours_per_cycle,
    ut_constants,
)


try:
    __version__ = get_version("UTide")
except Exception:
    __version__ = "unknown"

__all__ = [
    "solve",
    "solve_m",
    "reconstruct",
    "ut_constants",
    "constit_index_dict",
    "hours_per_cycle",
    "cycles_per_hour",
]
