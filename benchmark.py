import time
import utide as ut
import xarray as xr
from pandas import date_range
from utide.harmonics import FUV, linearized_freqs
from utide._time_conversion import _normalize_time

# 配置参数
DATA_PATH = "../harmonic_test/h572a.nc"
START_TIME = "2018-02-01 15:00:00"
PERIODS = 8000
FREQ = "h"
LATITUDE = 46.208
CONSTITUENTS = ["M2", "S2", "K1", "O1", "P1", "Q1", "N2", "K2"]

def load_ssh_data(path, periods):
    ds = xr.open_dataset(path)
    sea_level = ds["sea_level"].values
    ssh = sea_level[0, -periods:] / 1000.0
    return ssh


def benchmark_solve(t, ssh):
    # start = time.time()
    for _ in range(10):
        out = ut.solve(
            t,
            ssh,
            lat=LATITUDE,
            method="ols",
            conf_int="none",
            constit=CONSTITUENTS,
            order_constit="frequency",
            trend=False,
            nodal=True,
            verbose=False,
        )
    # elapsed = time.time() - start
    # print(f"ut.solve: total time is {elapsed:.3f} seconds")
    return out


def benchmark_solve_m(t, ssh, out1):
    # Use the frequency-ordered constituent names from the first solve call
    # This ensures the cache parameters match the internal constituent selection
    fuv_cache = FUV(
        _normalize_time(t, None),
        out1["aux"]["reftime"],
        out1["aux"]["lind"],
        30.0,
        [1, 0, False, False],
    )
    freqs = linearized_freqs(out1["aux"]["reftime"])

    # Use the frequency-ordered constituent names to match the cache
    ordered_constituents = list(out1["name"])
    print(ordered_constituents)

    # start = time.time()
    for _ in range(10):
        out = ut.solve_m(
            t,
            ssh,
            lat=LATITUDE,
            method="ols",
            conf_int="none",
            constit=ordered_constituents,  # Use ordered constituents to match cache
            order_constit="frequency",
            fuv_cache=fuv_cache,
            freqs=freqs,
            trend=False,
            nodal=True,
            verbose=False,
        )
    # elapsed = time.time() - start
    # print(f"ut.solve_m: total time is {elapsed:.3f} seconds")
    return out


if __name__ == "__main__":
    t = date_range(start=START_TIME, periods=PERIODS, freq=FREQ)
    ssh = load_ssh_data(DATA_PATH, PERIODS)
    out1 = benchmark_solve(t, ssh)
    out2 = benchmark_solve_m(t, ssh, out1)
    print(out1["A"])
    print(out1["name"])
    print(out2["A"])
    print(out2["name"])
