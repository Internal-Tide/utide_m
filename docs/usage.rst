使用指南与示例
==============

本节将向您展示如何使用 UTide 进行基本的潮汐分析和重构。

安装
----

.. code-block:: shell

    pip install utide

基本用法
--------

首先，让我们看看如何从一个时间序列中求解调和常数：

.. code-block:: python

    from utide import solve
    import numpy as np

    # 准备示例数据（时间、海平面高度、纬度）
    time = np.linspace(0, 100, 1000) # 时间，单位：天
    ssh = np.sin(2 * np.pi * time / 0.5) + np.random.normal(0, 0.1, 1000) # 简化的潮汐模拟
    lat = 30.0 # 纬度

    # 求解调和常数
    coef = solve(time, ssh, lat=lat, method='ols', conf_int='linear')

    # 输出结果
    print(coef.name) # 提取的分量名称
    print(coef.A)    # 提取的振幅

重构潮汐信号
------------

得到调和常数后，我们可以重构潮汐信号：

.. code-block:: python

    from utide import reconstruct

    # 使用求得的 coef 重构信号
    tide = reconstruct(time, coef)

    # 访问重构后的高度
    print(tide.h)

进阶用法：处理多个时间序列
--------------------------

如果您有多个时间序列（例如，在同一位置的多个年份，或者相同时间点的多个位置数据），可以使用 `solve_m` 进行更高效的计算。

.. code-block:: python

    from utide import solve_m
    # 更多关于 solve_m 的示例见 API 参考。
