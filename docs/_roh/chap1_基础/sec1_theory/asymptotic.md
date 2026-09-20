# 渐进均质化

## 双尺度渐进展开

设位移场 $\boldsymbol{u}^{\zeta}$ 和本征应变场 $\boldsymbol{\mu}^{\zeta}$ 可展开成如下形式：

$$
\begin{equation}\label{eq:asymptotic_expansion}
\begin{aligned}
\boldsymbol{u}^{\zeta}(\boldsymbol{x}) 
&= \boldsymbol{u}^{(0)}(\boldsymbol{x},\, \boldsymbol{x}/\zeta)
+ \zeta \boldsymbol{u}^{(1)}(\boldsymbol{x},\, \boldsymbol{x}/\zeta)
+ \zeta^2 \boldsymbol{u}^{(2)}(\boldsymbol{x},\, \boldsymbol{x}/\zeta)
+ \cdots, \\
\boldsymbol{\mu}^{\zeta}(\boldsymbol{x}) 
&= \boldsymbol{\mu}^{(0)}(\boldsymbol{x},\, \boldsymbol{x}/\zeta)
+ \zeta \boldsymbol{\mu}^{(1)}(\boldsymbol{x},\, \boldsymbol{x}/\zeta)
+ \zeta^2 \boldsymbol{\mu}^{(2)}(\boldsymbol{x},\, \boldsymbol{x}/\zeta)
+ \cdots,
\end{aligned}
\end{equation}
$$

式中，$\boldsymbol{u}^{(k)}$ 和 $\boldsymbol{\mu}^{(k)}$ 是定义在 $\Omega\times\Theta$ 的双尺度函数，$k=0,1,\ldots$ 是渐进展开式的阶数。对双尺度函数 $f$，当固定宏观尺度的坐标 $\boldsymbol{x}_{\mathtt{fix}}$ 时，$\tilde{f}: \boldsymbol{y} \mapsto f(\boldsymbol{x}_{\mathtt{fix}},\,  \boldsymbol{y})$ 是以单胞 $\Theta$ 为周期的函数，与定义式 $\eqref{eq:fine_scale_stiffness}$ 给出的弹性刚度张量 $\mathbb{L}$ 类似。将坐标变换式 $\eqref{eq:coarse_fine_scale_trans}$ 代入到双尺度函数中，也即 $f^{\zeta}(\boldsymbol{x})\triangleq f(\boldsymbol{x},\, \boldsymbol{x}/\zeta)$，再应用链式法则，函数 $f^{\zeta}$ 的偏导数等于

$$
\begin{equation}\notag
\frac{\partial f^{\zeta}}{\partial x_{i}}
= \frac{\partial f}{\partial x_{i}} + \frac{1}{\zeta} \frac{\partial f}{\partial y_{i}}.
\end{equation}
$$

记关于双尺度函数宏观和细观坐标的偏导数算子分别为 $\nabla_{\mathtt{x}}$ 和 $\nabla_{\mathtt{y}}$，那么上式可写成如下算子的形式，

$$
\begin{equation}\label{eq:two_scale_diff}
\nabla f^{\zeta}
= \nabla_{\mathtt{x}} f + \zeta^{-1} \nabla_{\mathtt{y}} f.
\end{equation}
$$

将位移和本征应变的渐进展开式代入到几何方程 $\eqref{eq:geom}$ 和本构关系 $\eqref{eq:cstt_eigenstrain}$ 中，得到应变和应力场的表达式为

$$
\begin{equation}\label{eq:asymp_strain_stress}
\begin{aligned}
\boldsymbol{\varepsilon}^{\zeta}(\boldsymbol{x})
&= \zeta^{-1} \underbrace{\nabla_{\mathtt{y}}^{\mathtt{s}} \boldsymbol{u}^{(0)}}_{\triangleq \boldsymbol{\varepsilon}^{(-1)}}
+ \zeta^{0}\underbrace{\big( \nabla_{\mathtt{x}}^{\mathtt{s}} \boldsymbol{u}^{(0)} + \nabla_{\mathtt{y}}^{\mathtt{s}} \boldsymbol{u}^{(1)} \big)}_{\triangleq\boldsymbol{\varepsilon}^{(0)}}
+ \zeta\underbrace{\big( \nabla_{\mathtt{x}}^{\mathtt{s}} \boldsymbol{u}^{(1)} + \nabla_{\mathtt{y}}^{\mathtt{s}} \boldsymbol{u}^{(2)} \big)}_{\triangleq\boldsymbol{\varepsilon}^{(1)}} 
+ \cdots,\\
\boldsymbol{\sigma}^{\zeta}(\boldsymbol{x})
&= \zeta^{-1} \underbrace{\mathbb{L} : \nabla_{\mathtt{y}}^{\mathtt{s}} \boldsymbol{u}^{(0)}}_{\triangleq\boldsymbol{\sigma}^{(-1)}}
+ \zeta^{0}\underbrace{\mathbb{L} : \big( \nabla_{\mathtt{x}}^{\mathtt{s}} \boldsymbol{u}^{(0)} + \nabla_{\mathtt{y}}^{\mathtt{s}} \boldsymbol{u}^{(1)} - \boldsymbol{\mu}^{(0)} \big)}_{\triangleq\boldsymbol{\sigma}^{(0)}} \\
&+ \zeta \underbrace{\mathbb{L} : \big( \nabla_{\mathtt{x}}^{\mathtt{s}} \boldsymbol{u}^{(1)} + \nabla_{\mathtt{y}}^{\mathtt{s}} \boldsymbol{u}^{(2)} - \boldsymbol{\mu}^{(1)} \big)}_{\triangleq\boldsymbol{\sigma}^{(1)}} + \cdots,
\end{aligned}
\end{equation}
$$

继续把渐进展开式 $\eqref{eq:asymptotic_expansion}$ 代入到平衡方程 $\eqref{eq:govern}$ 当中，然后匹配关于 $\zeta$ 不同阶次的方程项，由此给出一系列定解在单胞区域 $\Theta$ 内的平衡方程。这一过程繁琐且容易出错，因此本论文定义如下微分算子，将渐进展开的过程抽象化：

$$
\begin{equation}\label{eq:def_op}
\begin{aligned}
\mathcal{A}^{\zeta} (\bullet) &\triangleq - \nabla\cdot\big[ \mathbb{L}^{\zeta} : \nabla(\bullet) \big]
= - \nabla\cdot\big[ \mathbb{L}^{\zeta} : \nabla^{\mathtt{s}}(\bullet) \big], \\
\mathcal{B}^{\zeta} (\bullet) &\triangleq - \nabla\cdot\big[ \mathbb{L}^{\zeta} : (\bullet) \big],
\end{aligned}
\end{equation}
$$

式中，关于算子 $\mathcal{A}^{\zeta}$ 定义式中第二个等号应用到四阶弹性刚度张量的次对称性，也即 $L_{ijkl}^{\zeta}=L_{ijlk}^{\zeta}=L_{jikl}^{\zeta}$。将式 $\eqref{eq:cstt_eigenstrain}$ 代入到平衡方程 $\eqref{eq:govern}$ 中，并应用式 $\eqref{eq:def_op}$ 定义的算子，得到

$$
\begin{equation}\label{eq:govern_op}
\mathcal{A}^{\zeta}\boldsymbol{u}^{\zeta} = \mathcal{B}^{\zeta}\boldsymbol{\mu}^{\zeta} + \boldsymbol{f},
\end{equation}
$$

当上述算子作用在双尺度函数上时，将微分法则 $\eqref{eq:two_scale_diff}$ 代入上式中，得到

$$
\begin{equation}\label{eq:op_two_scale}
\mathcal{A}^{\zeta} 
= \zeta^{-2} \mathcal{A}_{\mathtt{yy}}
+ \zeta^{-1} \big( \mathcal{A}_{\mathtt{xy}} + \mathcal{A}_{\mathtt{yx}} \big)
+ \mathcal{A}_{\mathtt{xx}}, \quad
\mathcal{B}^{\zeta} 
= \zeta^{-1} \mathcal{B}_{\mathtt{y}}
+ \mathcal{B}_{\mathtt{x}},
\end{equation}
$$

式中，

$$
\begin{equation}\notag
\begin{cases}
\mathcal{A}_{\mathtt{yy}} (\bullet) 
\triangleq - \nabla_{\mathtt{y}} \cdot 
\big[ \mathbb{L} : \nabla_{\mathtt{y}} (\bullet) \big], \\
\mathcal{A}_{\mathtt{xy}} (\bullet) 
\triangleq - \nabla_{\mathtt{x}} \cdot 
\big[ \mathbb{L} : \nabla_{\mathtt{y}} (\bullet) \big], \\
\mathcal{A}_{\mathtt{yx}} (\bullet) 
\triangleq - \nabla_{\mathtt{y}} \cdot 
\big[ \mathbb{L} : \nabla_{\mathtt{x}} (\bullet) \big], \\
\mathcal{A}_{\mathtt{xx}} (\bullet) 
\triangleq - \nabla_{\mathtt{x}} \cdot 
\big[ \mathbb{L} : \nabla_{\mathtt{x}} (\bullet) \big],
\end{cases} \quad\quad
\begin{cases}
\mathcal{B}_{\mathtt{y}} (\bullet)  
\triangleq - \nabla_{\mathtt{y}} \cdot \big[ \mathbb{L} : (\bullet) \big], \\
\mathcal{B}_{\mathtt{x}} (\bullet)  
\triangleq - \nabla_{\mathtt{x}} \cdot \big[ \mathbb{L} : (\bullet) \big].
\end{cases}
\end{equation}
$$

记 $\mathcal{A}_{\mathtt{xy}}^{\mathtt{s}}\triangleq \mathcal{A}_{\mathtt{xy}} + \mathcal{A}_{\mathtt{yx}}$，将渐进展开式 $\eqref{eq:asymptotic_expansion}$ 代入到算子形式的控制方程 $\eqref{eq:govern_op}$ 中，应用式 $\eqref{eq:op_two_scale}$，得到

$$
\begin{equation}\notag
\begin{aligned}
\big[ \zeta^{-2} \mathcal{A}_{\mathtt{yy}}
+ \zeta^{-1} \mathcal{A}_{\mathtt{xy}}^{\mathtt{s}}
+ \mathcal{A}_{\mathtt{xx}} \big] 
&\big( \boldsymbol{u}^{(0)} + \zeta\boldsymbol{u}^{(1)} + \zeta^{2}\boldsymbol{u}^{(2)} + \cdots \big) \\
= \big[ \zeta^{-1} \mathcal{B}_{\mathtt{y}}
+ \mathcal{B}_{\mathtt{x}} \big]
&\big( \boldsymbol{\mu}^{(0)} + \zeta\boldsymbol{\mu}^{(1)} + \zeta^{2}\boldsymbol{\mu}^{(2)} + \cdots \big)  + \boldsymbol{f},
\end{aligned}
\end{equation}
$$

展开后，整理成类似于乘法表的形式为：

![image-20260801133047917](../../../assets/image-20260801133047917.png)

可以看到，在上式给出的矩阵中，从右上到左下对角线方程项的阶次是一样的。按照关于 $\zeta$ 的阶次 $q$ 整理方程，得到每一阶方程如下：

$$
\begin{equation}\label{eq:asymptotic_eqs}
\begin{aligned}
\mathcal{O}(\zeta^{q}):\quad\mathcal{A}_{\mathtt{yy}}\boldsymbol{u}^{(q+2)}
&= \boldsymbol{f}^{(q+1)},\\
\boldsymbol{f}^{(q+1)} 
&\triangleq - \mathcal{A}_{\mathtt{xy}}^{\mathtt{s}}\boldsymbol{u}^{(q+1)}
- \mathcal{A}_{\mathtt{xx}}\boldsymbol{u}^{(q)} \\
&+ \mathcal{B}_{\mathtt{y}}\boldsymbol{\mu}^{(q+1)} 
+ \mathcal{B}_{\mathtt{x}}\boldsymbol{\mu}^{(q)}
+ \delta_{0q} \boldsymbol{f},
\end{aligned}
\end{equation}
$$

式中，整数 $q\geq -2$，$\boldsymbol{u}^{(k)}$ 和 $\boldsymbol{\mu}^{(k)}$ 在 $k<0$ 时置为零。$\delta_{0q}$ 只在 $q=0$ 时等于 1，其它情况等于 0。

## 渐进方程组的适定性条件

式 $\eqref{eq:asymptotic_eqs}$ 中待求解的未知量是满足周期性边界条件的位移场 $\boldsymbol{u}^{(q+2)}$，由低阶的信息组成方程的右端项 $\boldsymbol{f}^{(q+1)}$（此处假设本征应变场是已知项），其中就包含低阶位移场 $\boldsymbol{u}^{(q)}$ 和 $\boldsymbol{u}^{(q+1)}$，因此应按照阶次由低到高的顺序依次求解。在实际求解这一组方程之前，本论文将引入方程解存在且唯一的条件。这不仅是为了理论推导的严谨，更是之后求解渐进方程组 $\eqref{eq:asymptotic_eqs}$ 的重要工具。

为刻画方程解存在且唯一条件，需要引入恰当的函数空间。定义 $L^{2}(\Theta)$ 为单胞域内平方可积函数组成的空间，$H^{1}(\Theta)$ 为一阶弱导数平方可积函数组成的空间。之后默认给出的函数空间定义在单胞区域 $\Theta$ 内，因此在函数空间记号中略去定义域。$H^{1}$ 空间中满足周期性条件的函数组成的子空间为 $H_{\mathtt{per}}^{1}$，也即运动学允许的位移场，定义为

$$
\begin{equation}\label{eq:kinematic_disp_space}
H_{\mathtt{per}}^{1} \triangleq \big\{ \boldsymbol{u}\in H^{1} \mid
\langle \boldsymbol{u} \rangle=\boldsymbol{0},\, \boldsymbol{u}|_{\partial\Theta^{-}} = \boldsymbol{u}|_{\partial\Theta^{+}}  \big\}
\end{equation}
$$

式中，$\boldsymbol{u}|_{\partial\Theta^{-}} = \boldsymbol{u}|_{\partial\Theta^{+}}$ 表示位移场 $\boldsymbol{u}$ 在相对边界点处取值相同，条件 $\langle \boldsymbol{u} \rangle=\boldsymbol{0}$ 排除了刚体平移。记 $(H_{\mathtt{per}}^{1})'$ 为 $H_{\mathtt{per}}^{1}$ 的对偶空间，由定义在 $H_{\mathtt{per}}^{1}$ 的所有连续线性泛函组成。方程 $\eqref{eq:asymptotic_eqs}$ 对应的弱形式可表述为：寻找函数 $\boldsymbol{u}^{(q+2)}\in H_{\mathtt{per}}^{1}$，使得对任意的测试函数 $\delta\boldsymbol{v}\in H_{\mathtt{per}}^{1}$，有如下等式恒成立：

$$
\begin{equation}\label{eq:asymptotic_eq_weak}
\int_{\Theta} \nabla_{\mathtt{y}} \delta\boldsymbol{v}
: \big( \mathbb{L}: \nabla_{\mathtt{y}}\boldsymbol{u}^{(q+2)} \big) \,\mathrm{d}\boldsymbol{y}
= \int_{\Theta} \delta\boldsymbol{v} \cdot \boldsymbol{f}^{(q+1)} \,\mathrm{d}\boldsymbol{y}.
\end{equation}
$$

现在不加证明地给出方程解存在且唯一的条件：设函数 $\boldsymbol{f}^{(q+1)}\in (H_{\mathtt{per}}^{1})'$， 那么弱形式 $\eqref{eq:asymptotic_eq_weak}$ 解存在且唯一的充分必要条件是：对任意定义在单胞域 $\Theta$ 的常值函数 $\boldsymbol{c}$，总有

$$
\begin{equation}\label{eq:theorem_exist_unique}
\int_{\Theta} \boldsymbol{c} \cdot \boldsymbol{f}^{(q+1)} \,\mathrm{d} \boldsymbol{y} = 0.
\end{equation}
$$

## 宏细观尺度分离

接下来本论文将逐阶次分析渐进展开方程组 $\eqref{eq:asymptotic_eqs}$ 的适定性。渐进展开方程组中 $\mathcal{O}(\zeta^{-2})$ 阶方程为

$$
\begin{equation}\label{eq:o-2}
\mathcal{A}_{\mathtt{yy}} \boldsymbol{u}^{(0)} = \boldsymbol{0},
\end{equation}
$$

方程右端项显然满足适定性条件 $\eqref{eq:theorem_exist_unique}$，并且常值函数是该方程的解，因此记方程的零阶位移解 $\boldsymbol{u}^{(0)}(\boldsymbol{x},\, \boldsymbol{y}) = \boldsymbol{u}^{(0)}(\boldsymbol{x})$。$\mathcal{O}(\zeta^{-2})$ 阶方程组的解总结如下：

$$
\begin{equation}\label{eq:o-2_sln}
\boldsymbol{u}^{(0)} = \boldsymbol{u}^{(0)}(\boldsymbol{x}), \quad
\boldsymbol{\varepsilon}^{(-1)} = \boldsymbol{0}, \quad
\boldsymbol{\sigma}^{(-1)} = \boldsymbol{0}.
\end{equation}
$$

渐进展开方程组中 $\mathcal{O}(\zeta^{-1})$ 阶方程为

$$
\begin{equation}\label{eq:o-1}
\begin{aligned}
\mathcal{A}_{\mathtt{yy}} \boldsymbol{u}^{(1)} 
&= -\mathcal{A}_{\mathtt{xy}}^{\mathtt{s}} \boldsymbol{u}^{(0)}
+ \mathcal{B}_{\mathtt{y}} \boldsymbol{\mu}^{(0)} \\
&= \nabla_{\mathtt{y}}\cdot \big[ \mathbb{L} : \big( \nabla_{\mathtt{x}}^{\mathtt{s}} \boldsymbol{u}^{(0)}
- \boldsymbol{\mu}^{(0)}\big) \big],
\end{aligned}
\end{equation}
$$

需要注意的是，式中的散度运算是定义在分布意义下的，因为方程右端项的极化应力场 $\boldsymbol{\lambda}\triangleq\mathbb{L} : \big( \nabla_{\mathtt{x}}^{\mathtt{s}}  \boldsymbol{u}^{(0)} - \boldsymbol{\mu}^{(0)}\big)\in L^{2}$ 在不同材料相的界面处一般不连续，求偏导之后得到集中于材料界面上的 Dirac 型分布。分布意义下的散度运算满足，对任意的满足周期性条件的无穷次可微测试函数 $\delta\boldsymbol{v}^{\infty}\in C_{\mathtt{per}}^{\infty}$，有

$$
\begin{equation}\notag
\int_{\Theta} \big( \nabla_{\mathtt{y}}\cdot\boldsymbol{\lambda} \big) \cdot \delta\boldsymbol{v}^{\infty}
\,\mathrm{d} \boldsymbol{y}
= -\int_{\Theta}\boldsymbol{\lambda} : \nabla_{\mathtt{y}}\delta\boldsymbol{v}^{\infty} \,\mathrm{d} \boldsymbol{y}
\end{equation}
$$

恒成立。选取测试函数 $\delta\boldsymbol{v}^{\infty}$ 为常值函数，可看到方程 $\eqref{eq:o-1}$ 右端项满足适定性条件 $\eqref{eq:theorem_exist_unique}$。该方程的解可通过线弹性方程的 Green 函数表示，将在下一节给出。

渐进展开方程组中 $\mathcal{O}(\zeta^{0})$ 阶方程为

$$
\begin{equation}\label{eq:o-0}
\begin{aligned}
\mathcal{A}_{\mathtt{yy}} \boldsymbol{u}^{(2)} 
&= -\mathcal{A}_{\mathtt{xy}}^{\mathtt{s}} \boldsymbol{u}^{(1)}
-\mathcal{A}_{\mathtt{xx}} \boldsymbol{u}^{(0)}
+ \mathcal{B}_{\mathtt{y}} \boldsymbol{\mu}^{(1)}
+ \mathcal{B}_{\mathtt{x}} \boldsymbol{\mu}^{(0)}
+ \boldsymbol{f}\\
&= \nabla_{\mathtt{y}} \cdot \big[ \mathbb{L} : \big( \nabla_{\mathtt{x}} \boldsymbol{u}^{(1)} 
- \boldsymbol{\mu}^{(1)} \big)\big]
+ \nabla_{\mathtt{x}} \cdot \big[ \mathbb{L} : \big( \nabla_{\mathtt{x}}^{\mathtt{s}}  \boldsymbol{u}^{(0)} + \nabla_{\mathtt{y}}^{\mathtt{s}}  \boldsymbol{u}^{(1)}
- \boldsymbol{\mu}^{(0)} \big)\big]
+ \boldsymbol{f}.
\end{aligned}
\end{equation}
$$

将上式右端项代入到适定性条件中，交换单胞域内的积分与对宏观坐标的微分运算顺序，就得到

$$
\begin{equation}\label{eq:temp1}
\nabla_{\mathtt{x}} \cdot \big\langle \mathbb{L} : \big( \nabla_{\mathtt{x}}^{\mathtt{s}}  \boldsymbol{u}^{(0)} + \nabla_{\mathtt{y}}^{\mathtt{s}}  \boldsymbol{u}^{(1)}
- \boldsymbol{\mu}^{(0)} \big) \big\rangle
+ \big\langle \boldsymbol{f} \big\rangle
= \boldsymbol{0},
\end{equation}
$$

式中的 $\langle \bullet \rangle$ 是公式 $\eqref{eq:def_volume_avg}$ 定义的体积平均运算。注意到在陈述控制方程 $\eqref{eq:govern}$ 时已经假设外力项不包含快速振荡项，因此 $\langle \boldsymbol{f} \rangle=\boldsymbol{f}$。此时定义宏观应力

$$
\begin{equation}\label{eq:def_macro_stress}
\boldsymbol{\sigma}^{\mathtt{c}} \triangleq \big\langle \mathbb{L} : \big( \nabla_{\mathtt{x}}^{\mathtt{s}}  \boldsymbol{u}^{(0)} + \nabla_{\mathtt{y}}^{\mathtt{s}}  \boldsymbol{u}^{(1)}
- \boldsymbol{\mu}^{(0)} \big) \big\rangle,
\end{equation}
$$

再代入到方程 $\eqref{eq:temp1}$ 中，形式上就得到了宏观尺度下的平衡方程：

$$
\begin{equation}\label{eq:macro_equilibrium}
\nabla_{\mathtt{x}} \cdot \boldsymbol{\sigma}^{\mathtt{c}} + \boldsymbol{f} = \boldsymbol{0}.
\end{equation}
$$



## 草稿

**误差估计**

数学上可以严格证明，当单胞特征长度 $l$ 与宏观问题特征长度 $L$ 的比值 $\zeta = l/L \to 0$ 时，考虑微结构分布的真实解 $\boldsymbol{u}^{\zeta}$ 与经式 $\eqref{eq:modulus_homo}$ 给出的等效模量均匀化后的方程解 $\boldsymbol{u}_{0}$ 满足如下误差估计：

$$
\| \boldsymbol{u}^{\zeta} - \boldsymbol{u}_{0} \|_{L^{2}} = \mathcal{O}(\sqrt{\zeta}).
$$




> 1. 这里应该给出一些均质化的一维算例，帮助读者理解推导过程，以及下面的两个 Remark。
> 2. 这里添加两个 Remarks，一个是说当 $\zeta \ll 1$ 但仍然是有限值的时候，给出一个误差估计；同时讨论如果存在边界层，裂纹等等，对渐进表达式的影响。这里可以参考 Yuan Zheng Partition of Unity，以及其它加入的分析高阶渐进展开和边界层影响的论文。在宏观尺度有没有必要加密网格呢？
> 3. 这里还应该指出，当 $\zeta \to 0^{+}$，也即尺度分离假设成立，此时给出一阶计算均质化的框架，式 $\eqref{eq:def_macro_stress}$ 给出细观尺度更新宏观应力的公式。Hill-Mandel 条件和边界条件在下一节给出。