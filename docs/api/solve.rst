solve
=====

.. autofunction:: utide.solve

.. note::

   `utide.reconstruct` 需要通过 `solve` 计算得到的置信区间。

参数详解
--------

* **t**: 时间数组。支持天数（相对于 `epoch`）、`numpy.datetime64` 或 `pandas.datetime`。
* **u**: 观测序列（如海面高度、流速分量）。
* **v**: 可选。如果是矢量流速，则为正交分量。
* **lat**: 纬度（度），必选参数。
* **constit**: 调和常数选择，默认为 'auto'（根据时间跨度自动选择）。
* **conf_int**: 置信区间计算方法：'linear', 'MC', 或 'none'。
* **method**: 求解方法：'ols' (普通最小二乘) 或 'robust' (鲁棒回归)。
* **trend**: 是否包含线性趋势，默认 True。
