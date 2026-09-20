# Hashin--Shtrikman 变分原理的构造方法

本文档关注如何构造 Hashin--Shtrikman (HS) 变分原理中的极化应力场，并得到各向同性非均质材料等效模量上下界的估计。

## HS 变分原理的陈述

HS 变分原理如下所示：

$$
\begin{equation}\label{eq:hs_vari}
\begin{aligned}
\mathcal{H}(\boldsymbol{p},\tilde{\boldsymbol{\varepsilon}}\mid\mathbb{L}_{0})
&= \int_{\Omega} \left( \frac{1}{2}\boldsymbol{p}:\tilde{\boldsymbol{\varepsilon}}
+ \boldsymbol{p}:\boldsymbol{\varepsilon}_{0} 
- \frac{1}{2} \boldsymbol{p}:(\mathbb{L} - \mathbb{L}_{0})^{-1}:\boldsymbol{p} \right) \mathrm{d}\Omega\\
&+ \frac{1}{2} \boldsymbol{\varepsilon}_{0} : \mathbb{L}_{0} : \boldsymbol{\varepsilon}_{0}.
\end{aligned}
\end{equation}
$$

泛函 $\eqref{eq:hs_vari}$ 含有两个依赖变量：极化应力 $\boldsymbol{p}$ 和应变场扰动 $\tilde{\boldsymbol{\varepsilon}}$，以及一个参数：参考模量 $\mathbb{L}_{0}$。依赖变量之间并不独立，而是通过如下辅助方程相互联系：

$$
\begin{equation}
\left( \mathbb{L}_{0} : \nabla_{s}\tilde{\boldsymbol{u}} + \boldsymbol{p} \right)\cdot\nabla = \boldsymbol{0}, \quad
\tilde{\boldsymbol{u}}|_{\partial\Omega} = \boldsymbol{0}.
\label{eq:subsd}
\end{equation}
$$

极化应力的可行域是**可积函数空间**；应变场扰动的可行域空间是**所有相容且在区域 $\Omega$ 内体积平均等于零的应变场**；宏观应变 $\boldsymbol{\varepsilon}_{0}$ 是给定值。泛函 $\eqref{eq:hs_vari}$ 在如下位置取**驻值**：

$$
\begin{equation}
\boldsymbol{p} = \Delta \mathbb{L}:(\boldsymbol{\varepsilon}_{0} + \tilde{\boldsymbol{\varepsilon}}).
\label{eq:stat}
\end{equation}
$$

当参考模量 $\mathbb{L}_{0}^{-}$ 使得 $\Delta\mathbb{L} \triangleq \mathbb{L} - \mathbb{L}_{0}$ **正定**时，泛函 $\eqref{eq:hs_vari}$ 在驻点处取**最大值**；而当参考模量 $\mathbb{L}_{0}^{+}$ 使 $\Delta\mathbb{L}$ 负定时，泛函在驻点处取**最小值**。因此有如下不等式成立：

$$
\mathcal{H}(\hat{\boldsymbol{p}}^{-},\tilde{\boldsymbol{\varepsilon}}\mid\mathbb{L}_{0}^{-})
\leq \frac{1}{2} \boldsymbol{\varepsilon}_{0} : \mathbb{L}^{c} : \boldsymbol{\varepsilon}_{0}
\leq \mathcal{H}(\hat{\boldsymbol{p}}^{+},\tilde{\boldsymbol{\varepsilon}}\mid\mathbb{L}_{0}^{+}),
$$

式中，$\hat{\boldsymbol{p}}$ 是任意满足可积条件的极化应力场。接下来将**构造分片常值分布的极化应力场**，并显式地写出泛函取值。

## 极化应力场的构造

考虑材料相均为**各向同性材料**，因此将任意二阶张量 $\boldsymbol{\varepsilon}$ 分解成体积部分和剪切部分：

$$
\boldsymbol{\varepsilon} = \varepsilon_{h}\boldsymbol{I} + \boldsymbol{\varepsilon}_{s}, \quad 
\varepsilon_{h} \triangleq\frac{1}{3}\varepsilon_{kk},
$$

那么泛函 $\eqref{eq:hs_vari}$ 就可以简化为

$$
\begin{equation}
\mathcal{H}(\boldsymbol{p},\tilde{\boldsymbol{\varepsilon}}\mid K_{0},G_{0})
= U_{0} + \tilde{U} + U_{\boldsymbol{p}},
\label{eq:hs_iso}
\end{equation}
$$

其中，由应力场扰动项贡献的能量 $\tilde{U}$，以及关于极化应力的二次型 $U_{p}$ 分别为

$$
\begin{aligned}
\tilde{U} &= \frac{1}{2}\int_{\Omega} \left( 
3p_{h}\tilde{\varepsilon}_{h} + \boldsymbol{p}_{s}:\tilde{\boldsymbol{\varepsilon}}_{s}
\right)\mathrm{d}\Omega, \\
U_{\boldsymbol{p}} &= \frac{1}{2}\int_{\Omega} \left(6p_{h}\varepsilon_{h}^{0} + 2\boldsymbol{p}_{s}:\boldsymbol{\varepsilon}^{0}_{s}
- \frac{p_{h}^{2}}{K-K_{0}} - \frac{\boldsymbol{p}_{s}:\boldsymbol{p}_{s}}{2(G-G_{0})} 
\right)\mathrm{d}\Omega.
\end{aligned}
$$

构造分片常值分布的极化应力场 $\hat{\boldsymbol{p}}$：

$$
\hat{\boldsymbol{p}}(\boldsymbol{x}) = \sum_{\alpha=1}^{N} \boldsymbol{p}^{(\alpha)}\chi^{(\alpha)}(\boldsymbol{x}),
$$

并代入到泛函 $U_{\boldsymbol{p}}$ 中，得到

$$
\begin{equation}
U_{\boldsymbol{p}}(p_{h}^{(\alpha)},\boldsymbol{p}_{s}^{(\alpha)})
= 3\varepsilon_{h}^{0} \left\langle p_{h} \right\rangle
- \frac{1}{2} \Big\langle \frac{p_{h}^{2}}{K-K_{0}} \Big\rangle
+ \boldsymbol{\varepsilon}_{s}^{0}: \left\langle \boldsymbol{p}_{s} \right\rangle
- \frac{1}{2} \Big\langle \frac{\boldsymbol{p}_{s}:\boldsymbol{p}_{s}}{2(G-G_{0})} \Big\rangle.
\label{eq:up_disc}
\end{equation}
$$

比较麻烦的是 $\tilde{U}$ 项，可通过 Fourier 方法获得这一相的表达式为

$$
\begin{equation}
2 \tilde{U}(p_{h}^{(\alpha)},\boldsymbol{p}_{s}^{(\alpha)})
= \alpha_{0}
\left\{ \left\langle p_{h}^{2} \right\rangle - \left\langle p_{h} \right\rangle^{2} \right\}
+ \beta_{0} \left\{ \left\langle \boldsymbol{p}_{s}:\boldsymbol{p}_{s} \right\rangle 
- \left\langle \boldsymbol{p}_{s} \right\rangle:\left\langle \boldsymbol{p}_{s} \right\rangle \right\},
\label{eq:uf_disc}
\end{equation}
$$

式中，

$$
\begin{equation}
\alpha_{0}\triangleq -\frac{3}{3 K^0 + 4\mu^0}, \quad
\beta_{0} \triangleq -\frac{3 K^0 + 6\mu^0}{5\mu^0 (3 K^0 + 4\mu^0)}.
\label{eq:hs_factor}
\end{equation}
$$

结合式 $\eqref{eq:up_disc}$，$\eqref{eq:uf_disc}$ 和 $\eqref{eq:hs_factor}$，代入各向同性材料的 HS 变分式 $\eqref{eq:hs_iso}$ 中，就得到关于有限维依赖变量 $(p_{h}^{(\alpha)},\boldsymbol{p}_{s}^{(\alpha)})$ 的**二次型表达式**：

$$
\begin{equation}
\begin{aligned}
&\hat{\mathcal{H}} (p_{h}^{(\alpha)},\boldsymbol{p}_{s}^{(\alpha)}; K_{0},G_{0}) - U_{0}\\
=&\, 3\varepsilon_{h}^{0} \left\langle p_{h} \right\rangle
- \frac{1}{2} \Big\langle \big( \frac{1}{K-K_{0}}-\alpha_{0} \big) p_{h}^{2} \Big\rangle
- \frac{1}{2}\alpha_{0} \left\langle p_{h} \right\rangle^{2} \\
+&\, \boldsymbol{\varepsilon}_{s}^{0}: \left\langle \boldsymbol{p}_{s} \right\rangle
- \frac{1}{2} \Big\langle \big( \frac{1}{2(G-G_{0})} - \beta_{0} \big) \boldsymbol{p}_{s}:\boldsymbol{p}_{s} \Big\rangle
- \frac{1}{2}\beta_{0}\left\langle \boldsymbol{p}_{s} \right\rangle:\left\langle \boldsymbol{p}_{s} \right\rangle.
\end{aligned}
\label{eq:hs_nb}
\end{equation}
$$

上式关于体积和剪切部分是**解耦**的，因此可以分别考虑。为使得式 $\eqref{eq:hs_nb}$ 取极值，对变量 $p_{h}^{(\alpha)}$ 求偏导数之后等于零，就得到 $N$ 个关于 $p_{h}^{(\alpha)}$ 的方程：

$$
\begin{aligned}
\frac{\partial \hat{\mathcal{H}}}{\partial p_{h}^{(\alpha)}}=0 
\quad\Rightarrow\quad
\big( \frac{1}{K^{(\alpha)}-K_{0}}-\alpha_{0} \big) p_{h}^{(\alpha)}
+ \alpha_{0} \sum_{\beta=1}^{N} c^{(\beta)} p_{h}^{(\beta)} = 3\varepsilon_{h}^{0}, \\
\alpha=1,2,\ldots,N.
\end{aligned}
$$

这个方程系统是对角矩阵加秩一形式，如果记列向量

$$
(\mathbf{m})_{i} = \frac{1}{K^{(i)}-K_{0}}-\alpha_{0}, \quad
(\mathbf{c})_{i} = c^{(i)}, \quad
\mathbf{a} = \alpha_0 \mathbf{1},
$$

那么方程系数的逆就等于

$$
(\mathrm{diag}\ \mathbf{m} + \mathbf{a}\mathbf{c}^{\top})^{-1}
= \mathrm{diag}\ \mathbf{m}^{-1} 
- \frac{\alpha_{0}\mathbf{m}^{-1} \mathbf{c}^{\top} \mathrm{diag}\ \mathbf{m}^{-1}}
{1+\alpha_{0}\mathbf{c}^{\top} \mathbf{m}^{-1}},
$$

再乘恒一向量 $\mathbf{1}$，就得到

$$
\mathbf{p}_{h} = 3\varepsilon_{h}^{0} (\mathrm{diag}\ \mathbf{m} + \mathbf{a}\mathbf{c}^{\top})^{-1} \mathbf{1}
= \frac{3\varepsilon_{h}^{0}}{1+\alpha_{0}\mathbf{c}^{\top} \mathbf{m}^{-1}}\mathbf{m}^{-1}
$$

因此，如果记

$$
K_{hs}^{(\alpha)} = (\mathbf{m})_{\alpha}^{-1}
= \frac{1}{K^{(\alpha)}-K_{0}}-\alpha_{0}, \quad
\bar{K}_{hs} =\mathbf{c}^{\top} \mathbf{m}^{-1}
= \sum_{\alpha=1}^{N} c^{(\alpha)} K_{hs}^{(\alpha)},
$$

那么方程的解等于

$$
p_{h}^{(\alpha)} = \frac{3K_{hs}^{(\alpha)}}{ 1+\alpha_{0} \bar{K}_{hs} } \varepsilon_{h}^{0}
$$

类似的，关于剪切部分的极化应力满足方程

$$
\begin{aligned}
\frac{\partial \hat{\mathcal{H}}}{\partial \boldsymbol{p}_{s}^{(\alpha)}}=\boldsymbol{0}
\quad\Rightarrow\quad
\big( \frac{1}{2(G^{(\alpha)}-G_{0})}-\beta_{0} \big) \boldsymbol{p}_{s}^{(\alpha)}
+ \beta_{0} \sum_{\beta=1}^{N} c^{(\beta)} \boldsymbol{p}_{s}^{(\beta)} = \boldsymbol{\varepsilon}_{s}^{0}, \\
\alpha=1,2,\ldots,N.
\end{aligned}
$$

如果记

$$
2G_{hs}^{(\alpha)} = \frac{1}{2(G^{(\alpha)}-G_{0})}-\beta_{0}, \quad
\bar{G}_{hs} = \sum_{\alpha=1}^{N} c^{(\alpha)} G_{hs}^{(\alpha)},
$$

那么方程的解为

$$
\boldsymbol{p}_{s}^{(\alpha)} = \frac{2G_{hs}^{(\alpha)}}{ 1+2\beta_{0} \bar{G}_{hs} } \boldsymbol{\varepsilon}_{s}^{0}.
$$

最后，将方程的解代入到  中，就得到

$$
\begin{aligned}
2\hat{\mathcal{H}} (p_{h}^{(\alpha)},\boldsymbol{p}_{s}^{(\alpha)}; K_{0},G_{0})
&= 9(\varepsilon_{h}^{0})^{2}\Big( K_{0} 
+ \sum_{\alpha=1}^{N}\frac{c^{(\alpha)} K_{hs}^{(\alpha)}}{ 1+\alpha_{0} \bar{K}_{hs} }
\Big) \\
&+ 2\boldsymbol{\varepsilon}_{s}^{0}:\boldsymbol{\varepsilon}_{s}^{0}
\Big( G_{0}
+ \sum_{\alpha=1}^{N}\frac{c^{(\alpha)}G_{hs}^{(\alpha)}}{ 1+2\beta_{0} \bar{G}_{hs} }
\Big)
\end{aligned}
$$

## 更新日志

### 2026/07/24

1. 创建了文档 `Hashin--Shtrikman 变分原理的构造方法.md`
