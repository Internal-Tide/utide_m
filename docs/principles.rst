原理说明
========

UTide 是一个用于潮汐调和分析和预测的工具，它是 Matlab 包 UTide 的 Python 实现。

潮汐调和分析
------------

潮汐调和分析的基本思想是，潮汐信号可以表示为一系列具有特定频率的正弦和余弦分量的叠加：

.. math::

   h(t) = H_0 + \sum_{j=1}^M [A_j \cos(\omega_j t + \phi_j)]

其中：
- :math:`H_0` 是平均海平面（mean sea level）。
- :math:`A_j` 是第 :math:`j` 个分量的振幅。
- :math:`\omega_j` 是第 :math:`j` 个分量的频率（通常由天文计算得出）。
- :math:`\phi_j` 是第 :math:`j` 个分量的相位。

瑞利判据 (Rayleigh Criterion)
-----------------------------

在选择调和常数时，UTide 使用瑞利判据来决定哪些分量可以从给定的时间序列中可靠地提取出来。瑞利判据要求两个分量的频率之差至少应为观测时间长度的倒数的一定倍数（默认为 1.0）：

.. math::

   |\omega_1 - \omega_2| \geq \frac{1}{T}

其中 :math:`T` 是观测序列的时间长度。

天文修正 (Nodal and Satellite Corrections)
------------------------------------------

UTide 包含节点修正（nodal correction）和卫星修正（satellite correction），以考虑长周期（如 18.61 年的月亮交点周期）对潮汐常数的影响。

解算方法
--------

UTide 提供两种主要的解算方法：
- **OLS (Ordinary Least Squares)**：普通最小二乘法，适用于数据分布良好的情况。
- **Robust Fit**：鲁棒回归方法，能够有效降低异常值对结果的影响。
