# 线上组装降阶方程

在线上阶段，ROH 求解的非线性方程组为：

$$
\Delta \boldsymbol{\varepsilon}^{(\alpha)}
= \mathbb{E}^{(\alpha)}:\Delta \boldsymbol{\varepsilon}^{(c)}
+ \sum_{\beta=1}^{N}\mathbb{P}^{(\alpha\beta)} : \Delta\boldsymbol{\mu}^{(\beta)}, \quad
\alpha=1,2,\ldots,N
$$

方程的未知量是分块的应变增量 $\{\Delta\boldsymbol{\varepsilon}^{(\alpha)}\}$，四阶张量 $\{ \mathbb{E}^{(\alpha)}, \mathbb{P}^{(\alpha\beta)} \}$ 是线下阶段计算弹性方程，或使用细观力学方法估计得到的影响张量，由非均质材料各组分的弹性常数和几何构型决定。本征应变增量 $\Delta\boldsymbol{\mu}^{(\beta)}$ 由分块的非弹性段本构确定的演化方程，以及状态变量确定，记成增量形式有

$$
\Delta\boldsymbol{\mu}^{(\beta)} = \Delta\boldsymbol{\varepsilon}^{(\beta)} 
- \mathbb{M}^{(\beta)}: \frac{\partial\Delta\boldsymbol{\sigma}^{(\beta)} }{\partial\Delta\boldsymbol{\varepsilon}^{(\beta)} } : \Delta\boldsymbol{\varepsilon}^{(\beta)}
$$

因此，本征应变增量 $\Delta\boldsymbol{\mu}^{(\beta)}$ 是关于相同材料相上应变增量 $\Delta\boldsymbol{\varepsilon}^{(\beta)}$ 的**非线性函数**。宏观应变增量 $\Delta\boldsymbol{\varepsilon}^{c}$ 与分块的应变增量之间满足如下体积平均关系：

$$
\Delta\boldsymbol{\varepsilon}^{c} = \sum_{\gamma=1}^{N} c^{(\gamma)} \Delta\boldsymbol{\varepsilon}^{(\gamma)}
$$

在方程组中代入上述关系，移项后得到

$$
\sum_{\gamma=1}^{N} \big( \delta^{(\alpha\gamma)}\mathbb{I}-c^{(\gamma)} \mathbb{E}^{(\alpha)} \big):\Delta\boldsymbol{\varepsilon}^{(\gamma)} 
= \sum_{\beta=1}^{N}\mathbb{P}^{(\alpha\beta)} : \Delta\boldsymbol{\mu}^{(\beta)}, \quad
\alpha=1,2,\ldots,N
$$

将非线性方程组写成矩阵形式为：

$$
\boldsymbol{\mathsf{E}}\Delta \boldsymbol{\mathsf{x}} = \boldsymbol{\mathsf{P}} \Delta \boldsymbol{\mathsf{u}}
$$

式中，维度为 $6N\times6N$ 的矩阵 $\boldsymbol{\mathsf{E}}$ 和 $\boldsymbol{\mathsf{P}}$ 分别由影响张量 $\{ \mathbb{E}^{(\alpha)} \}$ 和 $\{ \mathbb{P}^{(\alpha\beta)} \}$ 组装得到：

$$
\boldsymbol{\mathsf{E}}
\triangleq 
\boldsymbol{\mathsf{I}}-
\begin{pmatrix}
c^{(1)}\{ \mathbb{E}^{(1)} \} & c^{(2)}\{ \mathbb{E}^{(1)} \} & \cdots & c^{(N)}\{ \mathbb{E}^{(1)} \} \\
c^{(1)}\{ \mathbb{E}^{(2)} \} & c^{(2)}\{ \mathbb{E}^{(2)} \} & \cdots & c^{(N)}\{ \mathbb{E}^{(2)} \} \\
\vdots & \vdots & \ddots & \vdots \\
c^{(1)}\{ \mathbb{E}^{(N)} \} & c^{(2)}\{ \mathbb{E}^{(N)} \} & \cdots & c^{(N)}\{ \mathbb{E}^{(N)} \}
\end{pmatrix}, \quad
\boldsymbol{\mathsf{P}}
\triangleq \begin{pmatrix}
\{ \mathbb{P}^{(11)} \} & \{ \mathbb{P}^{(12)} \} & \cdots & \{ \mathbb{P}^{(1N)} \} \\
\{ \mathbb{P}^{(21)} \} & \{ \mathbb{P}^{(22)} \} & \cdots & \{ \mathbb{P}^{(2N)} \} \\
\vdots & \vdots & \ddots & \vdots \\
\{ \mathbb{P}^{(N1)} \} & \{ \mathbb{P}^{(N2)} \} & \cdots & \{ \mathbb{P}^{(NN)} \}
\end{pmatrix}
$$

式中，$\boldsymbol{\mathsf{I}}$ 是维度为 $6N\times6N$ 的单位矩阵。维度为 $6N\times 1$ 的未知量 $\Delta \boldsymbol{\mathsf{x}}$ 和 $\Delta \boldsymbol{\mathsf{u}}$ 分别由应变和本征应变增量组装得到：

$$
\Delta \boldsymbol{\mathsf{x}} = \begin{pmatrix}
\{ \Delta\boldsymbol{\varepsilon}^{(1)} \} \\ \{ \Delta\boldsymbol{\varepsilon}^{(2)} \} \\ \vdots \\ \{ \Delta\boldsymbol{\varepsilon}^{(N)} \}
\end{pmatrix}, \quad
\Delta \boldsymbol{\mathsf{u}} = \begin{pmatrix}
\{ \Delta\boldsymbol{\mu}^{(1)} \} \\ \{ \Delta\boldsymbol{\mu}^{(2)} \} \\ \vdots \\ \{ \Delta\boldsymbol{\mu}^{(N)} \}
\end{pmatrix}
$$

矩阵 $\boldsymbol{\mathsf{E}}$ 可以更进一步地表示为

$$
\boldsymbol{\mathsf{E}} = \boldsymbol{\mathsf{I}} - \boldsymbol{\mathsf{e}} \boldsymbol{\mathsf{c}}^{\top}, \quad
\underbrace{\boldsymbol{\mathsf{e}}}_{6N\times 6} \triangleq \begin{pmatrix}
\{ \mathbb{E}^{(1)} \} \\ \{ \mathbb{E}^{(2)} \} \\ \vdots \\ \{ \mathbb{E}^{(N)} \}
\end{pmatrix}, \quad
\underbrace{\boldsymbol{\mathsf{c}}}_{6N\times6} \triangleq \begin{pmatrix}
c^{(1)}\{ \mathbb{I} \} \\ c^{(2)}\{ \mathbb{I} \} \\ \vdots \\ c^{(N)}\{ \mathbb{I} \}
\end{pmatrix}
$$

由影响张量的恒等式

$$
\sum_{\alpha=1}^{N} c^{(\alpha)} \mathbb{E}^{(\alpha)} = \mathbb{I}
$$

可以得到

$$
\boldsymbol{\mathsf{c}}^{\top} \boldsymbol{\mathsf{e}} = \{ \mathbb{I} \} 
\Rightarrow \boldsymbol{\mathsf{E}} \boldsymbol{\mathsf{e}} 
= \underbrace{\boldsymbol{\mathsf{0}}}_{6N\times6}
$$

因此 $\boldsymbol{\mathsf{e}}$ 的每一个列向量属于矩阵 $\boldsymbol{\mathsf{E}}$ 的零空间。事实上，列向量组 $\boldsymbol{\mathsf{e}}$ 张成的空间**恰好等于**矩阵 $\boldsymbol{\mathsf{E}}$ 的零空间，

$$
\mathcal{N}(\boldsymbol{\mathsf{E}}) = \mathrm{span} \ \boldsymbol{\mathsf{e}}
$$

Fish et al. (2013) 给出的修正方法，如果设置 6 个变形模式全部为基体主导，修正后的 $\boldsymbol{\mathsf{P}}$ 实际上等于矩阵 $\boldsymbol{\mathsf{E}}$：

$$
\begin{equation}
\boldsymbol{\mathsf{P}} = \boldsymbol{\mathsf{E}}
\label{eq:infl_mod}
\end{equation}
$$

因此非线性方程组可以写成

$$
\boldsymbol{\mathsf{E}}(\Delta \boldsymbol{\mathsf{x}} - \Delta \boldsymbol{\mathsf{u}})
= \underbrace{\boldsymbol{0}}_{6N\times1}
\Rightarrow \Delta \boldsymbol{\mathsf{x}} - \Delta \boldsymbol{\mathsf{u}} \in \mathcal{N}(\boldsymbol{\mathsf{E}})\\
\Rightarrow \Delta \boldsymbol{\mathsf{x}} - \Delta \boldsymbol{\mathsf{u}} \in \mathrm{span} \ \boldsymbol{\mathsf{e}}
$$

所以

$$
\begin{equation}
\Delta \boldsymbol{\mathsf{x}} = \Delta \boldsymbol{\mathsf{u}} 
+ \underbrace{\boldsymbol{\mathsf{e}}}_{6N\times6}\ \underbrace{\boldsymbol{\lambda}}_{6\times1}
\label{eq:strain_fld}
\end{equation}
$$

式中的向量 $\boldsymbol{\lambda} \in \mathbb{R}^{6}$，其物理含义可在对上式左乘 $\boldsymbol{\mathsf{c}}^{\top}$ 后得到：

$$
\boldsymbol{\mathsf{c}}^{\top}\Delta \boldsymbol{\mathsf{x}} 
= \boldsymbol{\mathsf{c}}^{\top}\Delta \boldsymbol{\mathsf{u}} 
+ \boldsymbol{\mathsf{c}}^{\top}\boldsymbol{\mathsf{e}} \boldsymbol{\lambda}
\quad\Rightarrow \quad \boldsymbol{\lambda} = \Delta\boldsymbol{\varepsilon}^{c} - \boldsymbol{\mathsf{c}}^{\top}\Delta \boldsymbol{\mathsf{u}}
= \mathbb{M}^{c} : \Delta \boldsymbol{\sigma}^{c}
$$

考虑某一相（或某一个分块）完全破坏，也即本征应变增量等于应变增量：

$$
\Delta\boldsymbol{\varepsilon}^{(\alpha)} = \Delta \boldsymbol{\mu}^{(\alpha)}
$$

这可能对应为：

1. 损伤本构材料完全损伤，损伤因子 $d=1$;
2. 理想塑性材料进入塑性流动状态

而这一相（分块）的应变集中张量 $\mathbb{E}^{(\alpha)}$ 非奇异，也即

$$
\mathbb{E}^{(\alpha)} : \boldsymbol{\lambda} = \boldsymbol{0}
\iff \boldsymbol{\lambda} = \boldsymbol{0}
$$

那么根据式 $\eqref{eq:strain_fld}$，就得到

$$
\Delta \boldsymbol{\mathsf{x}} = \Delta \boldsymbol{\mathsf{u}} 
\quad\Rightarrow\quad \Delta \boldsymbol{\sigma}^{c} = \boldsymbol{0}
$$

这意味着：只要**某一**相（分块）完全破坏，那么无论这一相（分块）的体积分数有多小，都将强制所有其它材料相（分块）的**应变增量等于本征应变增量**，这有两种可能：

1. 其它材料相（分块）也恰好完全破坏，因此也有不等于零的应变增量
   

$$
   \Delta\boldsymbol{\varepsilon}^{(\alpha)} = \Delta\boldsymbol{\mu}^{(\alpha)}\neq\boldsymbol{0}
   $$

2. 由于本构的限制（比如说硬化模量不等于零的塑性本构，或还未完全损伤的分块），当本征应变增量不等于零时，应变增量不可能等于本征应变增量，此时分块处于“凝滞”状态，既不加载也不卸载：
   

$$
   \Delta\boldsymbol{\varepsilon}^{(\alpha)} = \Delta\boldsymbol{\mu}^{(\alpha)}=\boldsymbol{0}
   $$

无论其它分块出现什么样的状态，此时宏观的应力增量总是等于零的。
