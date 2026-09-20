考虑非均质的 RVE，包含任意多个夹杂相。RVE 体力项 $\boldsymbol{f}=\boldsymbol{0}$，边界条件为周期性边界条件。首先构造两组容许空间中的应力和应变场，分别为 $(\boldsymbol{\sigma}^{(1)}, \boldsymbol{\varepsilon}^{(1)})$ 和 $(\boldsymbol{\sigma}^{(2)}, \boldsymbol{\varepsilon}^{(2)})$。第一组考虑温度均匀变化，由此在每一相上产生热应变：

$$
\langle \boldsymbol{\sigma}^{(1)} \rangle=\boldsymbol{0}, \quad
\langle \boldsymbol{\varepsilon}^{(1)} \rangle = \boldsymbol{m}^{c} \Delta\theta, \quad
\boldsymbol{\varepsilon}^{(1)} = \mathbb{M}:\boldsymbol{\sigma}^{(1)} + \boldsymbol{m} \Delta \theta.
$$

式中，$\boldsymbol{m}$ 是材料的热膨胀二阶张量，在空间中可以非均匀分布；$\boldsymbol{m}^{c}$ 是均匀化后得到的宏观热膨胀张量。第二组考虑在 RVE 上施加宏观应力 $\boldsymbol{\sigma}_{0}$，无温度变化，由此得到的单胞内应变与应力场的响应：

$$
\langle \boldsymbol{\sigma}^{(2)} \rangle=\boldsymbol{\sigma}_{0}, \quad
\langle \boldsymbol{\varepsilon} \rangle^{(2)} = \mathbb{M}^{c} : \boldsymbol{\sigma}_{0}, \quad
\boldsymbol{\sigma}^{(2)} = \mathbb{F}:\boldsymbol{\sigma}_{0}, \quad
\boldsymbol{\varepsilon}^{(2)} = \mathbb{M}:\mathbb{F}:\boldsymbol{\sigma}_{0}.
$$

式中，四阶张量 $\mathbb{F}$ 是**应力影响函数**，通过应力影响函数和宏观应力，可以还原细观尺度下的应力场。

在构造完满足平衡方程的应力场和满足相容条件的应变场之后，接下来将应用虚功原理，得到宏观热膨胀张量 $\boldsymbol{m}^{c}$ 的表达式。注意到应变场 $\boldsymbol{\varepsilon}^{(1)}$ 和 $\boldsymbol{\varepsilon}^{(2)}$ 体积平均后不等于零，所以不能用做虚位移表示的虚功原理中应变场的扰动，因此使用由**虚应力表示的虚功原理**：

$$
\begin{equation}
\int_{\Omega} \delta\boldsymbol{\sigma} : \boldsymbol{\varepsilon} \,\mathrm{d}V
- |\Omega| \langle \delta\boldsymbol{\sigma} \rangle: \langle \boldsymbol{\varepsilon} \rangle = 0,
\label{eq:vs}
\end{equation}
$$

其中，虚应力要求满足平衡方程 $\delta\boldsymbol{\sigma}\cdot\nabla=\boldsymbol{0}$，以及周期性边界条件；应变场只需要满足相容条件即可。首先将 $\boldsymbol{\sigma}^{(1)}$ 作为虚应力，$\boldsymbol{\varepsilon}^{(2)}$ 作为试验函数，代入虚功原理 $\eqref{eq:vs}$ 中得到：

$$
\int_{\Omega} \boldsymbol{\sigma}^{(1)} : \mathbb{M} : \boldsymbol{\sigma}^{(2)} \,\mathrm{d}V = 0.
$$

然后再将 $\boldsymbol{\sigma}^{(2)}$ 和 $\boldsymbol{\varepsilon}^{(1)}$ 代入得到

$$
|\Omega|\boldsymbol{\sigma}_{0} : \boldsymbol{m}^{c}\Delta\theta
= \int_{\Omega} (\mathbb{F}:\boldsymbol{\sigma}_{0}) : \mathbb{M}:\boldsymbol{\sigma}^{(1)} \,\mathrm{d}V
+ \Delta\theta\int_{\Omega} (\mathbb{F}:\boldsymbol{\sigma}_{0}) : \boldsymbol{m} \,\mathrm{d}V.
$$

因为 $\boldsymbol{\sigma}^{(1)}$ 的体积平均等于零，所以上式等号右侧第一项积分等于零（这一结论同样来源于 $\eqref{eq:vs}$）。再根据 $\boldsymbol{\sigma}_{0}$ 的任意性，就得到

$$
\boldsymbol{m}^c = \frac{1}{|\Omega|} \int_{\Omega} \boldsymbol{m}: \mathbb{F}\,\mathrm{d}V.
$$

如果令 $\boldsymbol{\mu} = \boldsymbol{m} \Delta \theta$，$\boldsymbol{\mu}$ 是本征应变，那么上式就可以写作

$$
\begin{equation}
\boxed{\boldsymbol{\mu}^c = \frac{1}{|\Omega|} \int_{\Omega} \boldsymbol{\mu}: \mathbb{F}\,\mathrm{d}V.}
\label{eq:levin}
\end{equation}
$$

进一步的，如果本征应变场分片常值分布，那么上式又可以写成求和的形式：

$$
\begin{equation}
\boxed{\boldsymbol{\mu}^c = \sum_{\alpha=1}^{N} c^{(\alpha)} \boldsymbol{\boldsymbol{\mu}}^{(\alpha)} : \mathbb{F}^{(\alpha)}.}
\label{eq:levin_sum}
\end{equation}
$$

式中，$c^{(\alpha)}$ 是体积分数，$\boldsymbol{\mu}^{(\alpha)}$ 是区域 $\Omega^{(\alpha)}$ 内常值分布的本征应变，$\mathbb{F}^{(\alpha)}$ 是在区域 $\Omega^{(\alpha)}$ 内对应力影响函数的体积平均。以下是一些评论：

* 宏观尺度的应力和应变等于细观尺度中应力和应变场的**体积平均**，因此符号 $\boldsymbol{\sigma}^{c}$ 和 $\langle \boldsymbol{\sigma} \rangle$，$\boldsymbol{\varepsilon}^{c}$ 和 $\langle \boldsymbol{\varepsilon} \rangle$ 是可以不加区别地混用的；然而，**宏观本征应变并不是本征应变场的体积平均**，而是需要应力影响函数进行加权，如公式 $\eqref{eq:levin}$ 和 $\eqref{eq:levin_sum}$ 所示。
* 公式的推导过程只用到虚功原理，因此对于任意非弹性变形，在用本征应变进行表示之后，Levin 总是成立。