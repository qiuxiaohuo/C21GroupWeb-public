# 细观力学


统计平均（statistically homogeneous）是指：对描述非均质材料内各组分 $\alpha$ 分布的 $n$ 点概率函数 $\mathcal{S}_{n}^{(\alpha)}(\boldsymbol{x}_{1}, \boldsymbol{x}_{2}, \ldots,\boldsymbol{x}_{n})$，满足**平移不变性**，也即对空间中的任意向量 $\boldsymbol{y}$，概率函数满足

$$
\mathcal{S}_{n}^{(\alpha)}(\boldsymbol{x}_{1}+\boldsymbol{y}, \boldsymbol{x}_{2}+\boldsymbol{y}, \ldots,\boldsymbol{x}_{n}+\boldsymbol{y})
= \mathcal{S}_{n}^{(\alpha)}(\boldsymbol{x}_{1}, \boldsymbol{x}_{2}, \ldots,\boldsymbol{x}_{n}).
$$

一个直接的推论是，满足统计平均的非均质材料，其一维概率函数 $\mathcal{S}_{1}^{(\alpha)}$ 是常数。一维概率函数又可以解释为**材料组分的体积分数**，统计上平均的非均质材料的体积分数是常数：

$$
c^{(\alpha)} = \mathcal{S}_{1}^{(\alpha)}.
$$

这意味着，当在使用不变的体积分数描述非均质材料时，就已经默认该材料在统计上是平均的了。

遍历性假设（ergodicity assumption）是指：可以用无穷大非均质材料 $V_{\infty}$ 内的体积平均，代替有限区域 $V$ 内采样得到的**系综平均**（ensemble average）。以下使用符号

$$
\begin{equation}
\langle \cdot \rangle_{\infty} \triangleq \lim_{V\to \infty} \frac{1}{|V|} \int_{V} \cdot \ \mathrm{d}V
\label{eq:inf_vol_integral}
\end{equation}
$$

表示在无穷大区域 $V_{\infty}$ 内对某一物理量的体积平均。以下想要说明，在满足统计平均条件的非均质材料内，应用遍历性假设，可以自然导出 Hill-Mandel 条件的结论：

$$
\begin{equation}
\left\langle \boldsymbol{\sigma}: \boldsymbol{\varepsilon} \right\rangle_{\infty}
= \left\langle \boldsymbol{\sigma}\right\rangle_{\infty}: \left\langle \boldsymbol{\varepsilon} \right\rangle_{\infty}.
\label{eq:hill_condition}
\end{equation}
$$

为证明上述结论，首先考虑在一个有限区域内 $V$ 的体积平均，也即

$$
\left\langle \boldsymbol{\sigma}:\boldsymbol{\varepsilon} \right\rangle
= \frac{1}{|V|} \int_{V} \boldsymbol{\sigma} : \boldsymbol{\varepsilon} \ \mathrm{d}V.
\tag{$a$}
$$

如果将应力与应变场分解成平均常值和扰动场：

$$
\boldsymbol{\varepsilon} = \langle \boldsymbol{\varepsilon} \rangle + \boldsymbol{\varepsilon}^{\prime}, \quad
\boldsymbol{\sigma} = \langle \boldsymbol{\sigma} \rangle + \boldsymbol{\sigma}^{\prime}, \quad
\left\langle \boldsymbol{\varepsilon}^{\prime} \right\rangle
=\left\langle \boldsymbol{\sigma}^{\prime} \right\rangle
= \boldsymbol{0},
$$

那么式 $(a)$ 可以重新写成

$$
\begin{aligned}
\left\langle \boldsymbol{\sigma}:\boldsymbol{\varepsilon} \right\rangle
&= \left\langle \left( \langle \boldsymbol{\sigma} \rangle + \boldsymbol{\sigma}^{\prime} \right)
:\left( \langle \boldsymbol{\varepsilon} \rangle + \boldsymbol{\varepsilon}^{\prime} \right) \right\rangle \\
&= \langle \boldsymbol{\sigma} \rangle : \langle \boldsymbol{\varepsilon} \rangle
+ \langle \boldsymbol{\sigma}^{\prime} : \boldsymbol{\varepsilon}^{\prime} \rangle
+ \underbrace{\langle \boldsymbol{\sigma}^{\prime} \rangle : \langle \boldsymbol{\varepsilon} \rangle}_{=0}
+ \underbrace{\langle \boldsymbol{\sigma} \rangle : \langle \boldsymbol{\varepsilon}^{\prime} \rangle}_{=0},
\end{aligned}
$$

也即

$$
\begin{equation}
\left\langle \boldsymbol{\sigma}:\boldsymbol{\varepsilon} \right\rangle
= \langle \boldsymbol{\sigma} \rangle : \langle \boldsymbol{\varepsilon} \rangle
+ \langle \boldsymbol{\sigma}^{\prime} : \boldsymbol{\varepsilon}^{\prime} \rangle.
\label{eq:energy_product}
\end{equation}
$$

如果将位移场也相应地分解为线性位移和扰动项：

$$
\boldsymbol{u} = \left\langle \boldsymbol{\varepsilon} \right\rangle \cdot \boldsymbol{x} + \boldsymbol{u}^{\prime}, \quad
\boldsymbol{\varepsilon}^{\prime} = \nabla_{s} \boldsymbol{u}^{\prime}.
$$

代入式 $\eqref{eq:energy_product}$ 等式右边的最后一项得到

$$
\begin{equation}
\langle \boldsymbol{\sigma}^{\prime}:\boldsymbol{\varepsilon}^{\prime} \rangle
= \frac{1}{|V|} 
\left( \int_{\partial V} \boldsymbol{\sigma}^{\prime} : (\boldsymbol{u}^{\prime}\otimes\boldsymbol{n}) \ \mathrm{d}S
- \int_{V} (\boldsymbol{\sigma}^{\prime}\cdot\nabla) \cdot \boldsymbol{u}^{\prime} \ \mathrm{d}V
\right)
= \frac{1}{|V|} 
\int_{\partial V} \boldsymbol{\sigma}^{\prime} : (\boldsymbol{u}^{\prime}\otimes\boldsymbol{n}) \ \mathrm{d}S.
\label{eq:energy_fluc}
\end{equation}
$$

上式中分布积分后的体积积分等于零是因为区域内满足平衡方程 $\boldsymbol{\sigma}^{\prime} \cdot \nabla = \boldsymbol{0}$。注意到在推导式 $\eqref{eq:energy_product}$ 和式 $\eqref{eq:energy_fluc}$ 时，并没有引入额外的假设，也即对任何区域 $V$ 都是成立的。而**遍历性假设要求，关于扰动项的边界积分必须等于零**：

$$
\begin{equation}
\langle \boldsymbol{\sigma}^{\prime}:\boldsymbol{\varepsilon}^{\prime} \rangle_{\infty}
= \lim_{V\to \infty}\frac{1}{|V|} 
\int_{\partial V} \boldsymbol{\sigma}^{\prime} : (\boldsymbol{u}^{\prime}\otimes\boldsymbol{n}) \ \mathrm{d}S
= 0,
\label{eq:ergodicity_assump}
\end{equation}
$$

将上式代入到式 $\eqref{eq:energy_product}$ 中，就得到式 $\eqref{eq:hill_condition}$ 的结论。

> **Remark**
>
> 1. 除式 $\eqref{eq:ergodicity_assump}$ 外，遍历性假设要求扰动项的边界积分等于零还有如下例子：
>
>    * 边界处位移场扰动项等于零，

$$
\langle \boldsymbol{\varepsilon}^{\prime} \rangle_{\infty} = \boldsymbol{0}
\Rightarrow\lim_{V\to\infty} \frac{1}{|V|} \int_{\partial V} \boldsymbol{u}^{\prime} \otimes \boldsymbol{n} \ \mathrm{d}S = \boldsymbol{0}.
$$

>   * 边界处应力扰动项等于零，

$$
\langle \boldsymbol{\sigma}^{\prime} \rangle_{\infty} = \boldsymbol{0}
\Rightarrow\lim_{V\to\infty} \frac{1}{|V|} \int_{\partial V} \boldsymbol{x} \otimes ( \boldsymbol{\sigma}^{\prime} \cdot \boldsymbol{n} ) \ \mathrm{d}S = \boldsymbol{0}.
$$

>然而，文档作者确实没有找到推导出式 $\eqref{eq:ergodicity_assump}$ 的方法。
>
>2. Hill 在他的经典论文中给出代表体积元的说明如下，式 $\eqref{eq:ergodicity_assump}$ 可以看作是对标红文字的数学表述：“通过取足够大体积的样本，在边界处的（外力或位移扰动项）对整体平均可以忽略不计”。
>
>   <img src="../../../../../../../../../../../../assets/image-20251215193714676.png" alt="image-20251215193714676" style="zoom:50%;" />
