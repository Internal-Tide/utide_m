# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

UTide is a Python re-implementation of the Matlab UTide package for tidal analysis and prediction. It calculates astronomical harmonics separately from the solve harmonic analysis, making `solve_m` more efficient for multiple time series with the same time vector.

## Core Architecture

The package is organized around several key modules:

- `_solve.py`: Main tidal harmonic analysis with `solve()` and `solve_m()` functions
- `_reconstruct.py`: Signal reconstruction from solved coefficients
- `astronomy.py`: Astronomical calculations and variable computations
- `harmonics.py`: Harmonic analysis functions including `ut_E()` and `ut_E_m()`
- `constituent_selection.py`: Constituent selection algorithms
- `confidence.py`: Confidence interval calculations
- `ellipse_params.py`: Conversion to ellipse parameters
- `utilities.py`: Utility classes including the `Bunch` data structure
- `_time_conversion.py`: Time normalization utilities

## Common Development Commands

### Package Management
```bash
# Install in development mode
pip install -e .

# Install dependencies
pip install -r requirements.txt
pip install -r requirements-dev.txt
```

### Testing and Quality
```bash
# Run tests (if test suite exists)
pytest

# Run linting
flake8 utide/
black --check utide/
isort --check-only utide/

# Format code
black utide/
isort utide/
```

### Build and Distribution
```bash
# Build package
python -m build

# Check package
check-manifest

# Run benchmark (using benchmark.py)
python benchmark.py
```

## Key Functions and Usage

The main public API consists of:
- `solve()`: Primary harmonic analysis function
- `solve_m()`: Optimized for multiple time series with same time vector
- `reconstruct()`: Reconstruct tidal signal from coefficients

All functions expect time series data and return `Bunch` objects (dict-like with attribute access) containing coefficients and metadata.

## Data Files

The package includes essential data files in `utide/data/`:
- `ut_constants.npz`: Tidal constituent constants
- `FUV0.npz`: Astronomical coefficients

These are automatically included via package data configuration in `setup.cfg`.

## Configuration

Default options for analysis are defined in `_solve.py:16-31` and can be customized via the `opts` parameter in solve functions. The package supports various analysis methods including OLS and robust fitting.