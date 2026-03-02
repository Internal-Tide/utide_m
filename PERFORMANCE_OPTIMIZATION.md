# UTide Performance Optimization Guide

This document describes the performance optimizations implemented to speed up harmonic analysis in UTide, particularly for large multi-dimensional datasets (e.g., 3D spatial grids).

## Key Optimizations

### 1. Multi-dimensional (3D) Vectorized Solve
The `solve()` function now natively supports multi-dimensional inputs. If `u` (and optionally `v`) has more than one dimension (e.g., `(time, ny, nx)`), UTide will automatically solve all spatial points simultaneously.

**How it works:**
- The time-series data is flattened into a `(time, n_points)` matrix.
- The model matrix `B` (harmonics, trend, mean) is constructed only once.
- `np.linalg.lstsq` is called with multiple right-hand sides, solving for all points in a single vectorized operation.
- Results are reshaped back to the original spatial dimensions.

**Example Usage:**
```python
import numpy as np
from utide import solve

# Generate synthetic 3D data: (time=8000, y=10, x=10)
t = np.linspace(0, 365, 8000)
u = np.random.randn(8000, 10, 10)

# All 100 spatial points solved simultaneously
coef = solve(t, u, lat=30, method='ols', conf_int='none')

print(coef['A'].shape)  # Output: (n_constituents, 10, 10)
```

### 2. Optimized Nodal Corrections (`FUV`)
The construction of the nodal/satellite correction matrix in `utide/harmonics.py` has been optimized.

**Optimization detail:**
- A Python loop that summed satellite contributions for each constituent was replaced with NumPy's vectorized `np.add.at` operation.
- This significantly reduces the overhead of `ut_E` and `FUV`, especially when dealing with many constituents.

## Benchmark Results

For a dataset with 8000 time points and a 40x40 spatial grid (1600 points):

| Method | Total Time | Speedup |
| :--- | :--- | :--- |
| **Original (Looping over points)** | ~300.0 seconds | 1.0x |
| **Optimized (Vectorized Solve)** | **~4.6 seconds** | **~65x** |

*Note: Benchmarks were performed on a standard CPU environment. Actual speedup may vary depending on data size and system architecture.*

## Considerations for Multi-dimensional Solve
- **Memory Usage:** Vectorizing over many points increases memory consumption. For very large grids, consider processing in chunks.
- **Confidence Intervals:** Currently, confidence intervals (`conf_int='linear'` or `'MC'`) are not yet supported in the vectorized multi-dimensional path.
- **NaN Handling:** The vectorized solver assumes all spatial points share a consistent set of valid time points. If a spatial point is all NaNs, the result for that point will be NaN.
