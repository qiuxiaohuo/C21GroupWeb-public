# 满足 Hill--Mandel 条件的边界类型

均质化方法要求在给定宏观应变 $\boldsymbol{\varepsilon}_{0}$ 或宏观应力 $\boldsymbol{\sigma}_{0}$ 时，在单胞 $\Theta$ 或代表体积元 $V$（以下均用 $V$ 表示）求解如下无体力项的弹性力学方程：

$$
\begin{equation}
\boldsymbol{\sigma} \cdot \nabla = \boldsymbol{0}, \quad
\boldsymbol{\sigma} = \mathbb{L} : \boldsymbol{\varepsilon}, \quad
\langle \boldsymbol{\varepsilon} \rangle = \boldsymbol{\varepsilon}_{0} \quad or \quad
\langle \boldsymbol{\sigma} \rangle = \boldsymbol{\sigma}_{0}.
\label{eq:bvp}
\end{equation}
$$

式中，$\langle \cdot \rangle$ 表示在区域 $V$ 内的体积积分。上式还未组成一个良定义的 BVP 问题，为在区域 $V$ 内定解该问题，需要施加指定类型的边界条件。这些边界条件应对任意给定的宏观应变或应力，满足 Hill-Mandel 条件 $\langle \boldsymbol{\sigma}:\boldsymbol{\varepsilon} \rangle = \langle \boldsymbol{\sigma} \rangle : \langle \boldsymbol{\varepsilon} \rangle$。能够在有限尺寸区域内自动满足 Hill-Mandel 条件的边界条件一般是如下三种，或是这三种边界条件的混合：

* **线性位移边界条件**。边界点处的位移指定为

$$
\boldsymbol{u}|_{\partial V} = \boldsymbol{\varepsilon}_{0} \cdot \boldsymbol{y}, \quad \boldsymbol{y} \in \partial V.
$$

* **均匀应力边界条件**。设边界点处的外法向量为 $\boldsymbol{n}$，该点处的力边界条件指定为

$$
\boldsymbol{t}|_{\partial V} = \boldsymbol{\sigma}_{0} \cdot \boldsymbol{n}.
$$

* **周期性边界（约束）条件**。施加周期性边界条件要求单胞 $\Theta$ 能在空间中密铺，因此一般设置单胞形状为长方体。此时，单胞在边界处有两个相对的点，分别记为 $\boldsymbol{y}^{-}$ 和 $\boldsymbol{y}^{+}$，在周期性边界意义下实为同一个点 $\boldsymbol{x}$。周期性边界条件要求，在相对点处位移场相等：

$$
\boldsymbol{u}^{\prime}(\boldsymbol{y}^{-})
= \boldsymbol{u}^{\prime}(\boldsymbol{y}^{+}), \quad \boldsymbol{y} \in \partial \Theta.
$$

