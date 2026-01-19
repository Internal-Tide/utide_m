UTide
=====
|travis| |license| |conda| |downloads| |anaconda_cloud| |appveyor|

.. |travis| image:: https://travis-ci.org/wesleybowman/UTide.svg?branch=master
   :target: https://travis-ci.org/wesleybowman/UTide

.. |license| image:: https://anaconda.org/conda-forge/utide/badges/license.svg
   :target: https://choosealicense.com/licenses/mit/

.. |conda| image:: https://anaconda.org/conda-forge/utide/badges/installer/conda.svg
   :target: https://anaconda.org/conda-forge/utide

.. |downloads| image:: https://anaconda.org/conda-forge/utide/badges/downloads.svg
   :target: https://anaconda.org/conda-forge/utide

.. |anaconda_cloud| image:: https://anaconda.org/conda-forge/utide/badges/version.svg
   :target: https://anaconda.org/conda-forge/utide

.. |appveyor| image:: https://ci.appveyor.com/api/projects/status/4o163ma4ehhr3q48/branch/master?svg=true
   :target: https://ci.appveyor.com/project/wesleybowman/utide/branch/master


Python re-implementation of the Matlab package UTide.

Calculates astronomy seperately from the solve harmonic analysis, so the function solve_m
is more efficient for multiple time series with the same time vector.
Still in heavy development--everything is subject to change!

Note: the user interface differs from the Matlab version, so
consult the Python function docstrings to see how to specify
parameters. Some functionality from the Matlab version is
not yet available. For more information see:

::

    Codiga, D.L., 2011. Unified Tidal Analysis and Prediction Using the
    UTide Matlab Functions. Technical Report 2011-01. Graduate School
    of Oceanography, University of Rhode Island, Narragansett, RI.
    59pp. ftp://www.po.gso.uri.edu/pub/downloads/codiga/pubs/
    2011Codiga-UTide-Report.pdf

    UTide v1p0 9/2011 d.codiga@gso.uri.edu
     http://www.po.gso.uri.edu/~codiga/utide/utide.htm

Installation
============

.. code:: shell

    pip install utide

If you are using conda,

.. code:: shell

    conda install utide -c conda-forge


The public functions can be imported using

.. code:: python

    from utide import solve, reconstruct

A sample call would be

.. code:: python

    from utide import solve

    coef = solve(time, time_series_u, time_series_v,
                 lat=30,
                 nodal=False,
                 trend=False,
                 method='ols',
                 conf_int='linear',
                 Rayleigh_min=0.95,)

if use model data you can doing this to speed up the calculation  

.. code:: python

    from utide import solve_m

def benchmark_solve(t, ssh, repeat=10):
    out = None
    for _ in range(repeat):
        out = ut.solve(
            t, ssh, lat=LATITUDE, method="ols", conf_int="none",
            constit=CONSTITUENTS, order_constit="frequency",
            trend=False, nodal=True, verbose=False
        )
    return out

def benchmark_solve_m(t, ssh, out1, repeat=10):
    F,U,V = FUV(_normalize_time(t, None), out1["aux"]["reftime"], out1["aux"]["lind"], 30.0, [1, 0, False, False])
    fuv_cache = (F, U, V, out1["aux"]["lind"])
    freqs = linearized_freqs(out1["aux"]["reftime"])
    for _ in range(repeat):
        out = ut.solve_m(
            t, ssh, lat=LATITUDE, method="ols", conf_int="none",
            constit=CONSTITUENTS, order_constit="frequency",
            fuv_cache=fuv_cache, freqs=freqs,
            trend=False, nodal=True, verbose=False
        )
    return out  
    
For more examples see the
`notebooks <https://nbviewer.jupyter.org/github/wesleybowman/UTide/tree/master/notebooks/>`__
folder.
