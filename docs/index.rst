UTide 官方文档 (中文)
=====================

.. toctree::
   :maxdepth: 2
   :caption: 目录

   principles
   usage
   api/index

简介
----

UTide 是潮汐调和分析和预测包 UTide 的 Python 实现。它可以分析一维（标量，如海平面高度）或二维（矢量，如流速分量）的时间序列。

主要特点
--------

* 包含天文修正常量（Nodal and Satellite Corrections）。
* 自动根据观测序列的时间长度进行调和常数的选择（瑞利判据）。
* 提供普通最小二乘（OLS）和鲁棒回归（Robust Fit）两种求解算法。
* 支持多维数据（如 3D 空间数据）的向量化解算。

索引与搜索
----------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
