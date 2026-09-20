# 均匀场变换分析

从上一节中可以看到，宏观尺度的应力 $\boldsymbol{\sigma}^{\mathtt{c}}$ 是通过求解细观尺度 $\mathcal{O}(\zeta^{-1})$ 阶的单胞控制方程 $\eqref{eq:o-1}$ 之后，再在单胞区域内作体积平均 $\eqref{eq:def_macro_stress}$ 得到的。本节关注如何求解该单胞方程。

## 单胞方程的弱形式

这一小节将重新描述渐进均质化得到的单胞方程，并给出对应的弱形式，由此看到和一阶计算均质化能量一致性条件的关联。为了记法上的清晰，将渐进均质化得到的不同阶场变量重写为

$$
\begin{equation}\notag
\left\{\begin{aligned}
&\boldsymbol{u}^{(0)}(\boldsymbol{x})\to \boldsymbol{u}^{\mathtt{c}}(\boldsymbol{x}), \quad
\boldsymbol{u}^{(1)}(\boldsymbol{x},\, \boldsymbol{y})\to \tilde{\boldsymbol{u}}^{\mathtt{f}}(\boldsymbol{x},\, \boldsymbol{y}), \\
&\nabla_{\mathtt{x}}^{\mathtt{s}}\boldsymbol{u}^{(0)}(\boldsymbol{x})\to \boldsymbol{\varepsilon}^{\mathtt{c}}(\boldsymbol{x}), \quad
\nabla_{\mathtt{y}}^{\mathtt{s}}\boldsymbol{u}^{(1)}(\boldsymbol{x},\, \boldsymbol{y})\to \tilde{\boldsymbol{\varepsilon}}^{\mathtt{f}}(\boldsymbol{x},\, \boldsymbol{y}), \\
&\boldsymbol{\mu}^{(0)}(\boldsymbol{x},\, \boldsymbol{y})\to \boldsymbol{\mu}^{\mathtt{f}}(\boldsymbol{x},\, \boldsymbol{y}),
\end{aligned}\right.
\end{equation}
$$

式中，上标 $\square^{\mathtt{f}}$ 表示细观尺度的场变量，与细观尺度的坐标 $\boldsymbol{y}$ 相关；$\tilde{\square}$ 描述场变量的扰动项，其体积积分 $\langle \tilde{\square} \rangle$ 等于零，之后将进一步说明。由于单胞方程的定解区域在单胞 $\Theta$ 内，宏观坐标 $\boldsymbol{x}$ 可视作参数，本节将暂时略去场变量对 $\boldsymbol{x}$ 的依赖关系。定义细观尺度的应变和应力场分别为

$$
\begin{equation}\label{eq:fine_scale_cstt}
\boldsymbol{\varepsilon}^{\mathtt{f}}(\boldsymbol{y}) 
\triangleq \boldsymbol{\varepsilon}^{\mathtt{c}} + \tilde{\boldsymbol{\varepsilon}}^{\mathtt{f}}(\boldsymbol{y}),\quad
\boldsymbol{\sigma}^{\mathtt{f}}(\boldsymbol{y})
\triangleq \mathbb{L}(\boldsymbol{y}):\big( \boldsymbol{\varepsilon}^{\mathtt{c}}
+ \tilde{\boldsymbol{\varepsilon}}^{\mathtt{f}}(\boldsymbol{y}) - \boldsymbol{\mu}^{\mathtt{f}}(\boldsymbol{y}) \big),
\end{equation}
$$

于是单胞方程可重新表示为无体力项的平衡方程，

$$
\begin{equation}\label{eq:fine_scale_equili}
\nabla_{\mathtt{y}} \cdot \boldsymbol{\sigma}^{\mathtt{f}} = \boldsymbol{0}, \quad \boldsymbol{y} \in \Theta.
\end{equation}
$$

限定位移 $\tilde{\boldsymbol{u}}^{\mathtt{f}}$ 和应力 $\boldsymbol{\sigma}^{\mathtt{f}}$ 在边界处分别满足周期和反周期边界条件：

$$
\begin{equation}\label{eq:fine_scale_bc}
\tilde{\boldsymbol{u}}^{\mathtt{f}}|_{\partial \Theta^{-}} = \tilde{\boldsymbol{u}}^{\mathtt{f}}|_{\partial \Theta^{+}},\quad
\boldsymbol{\sigma}^{\mathtt{f}}\cdot\boldsymbol{n}|_{\partial \Theta^{-}} 
= -\boldsymbol{\sigma}^{\mathtt{f}}\cdot\boldsymbol{n}|_{\partial \Theta^{+}},
\end{equation}
$$

此时式 $\eqref{eq:fine_scale_equili}$、$\eqref{eq:fine_scale_cstt}$ 和 $\eqref{eq:fine_scale_bc}$ 共同组成细观尺度的单胞方程。接下本论文将给出与单胞方程等价的弱形式，读者将看到，这一弱形式恰好是 Hill--Mandel 宏细观能量一致性原理的数学表达。选取运动学允许的测试函数 $\delta\tilde{\boldsymbol{v}}$，与控制方程 $\eqref{eq:fine_scale_equili}$ 相乘后在单胞域内积分，并应用分部积分公式，就得到

$$
\begin{equation}\notag
\int_{\Theta} \big( \nabla_{\mathtt{y}}\cdot \boldsymbol{\sigma}^{\mathtt{f}} \big) \cdot \delta\tilde{\boldsymbol{v}}
\,\mathrm{d}\boldsymbol{y}
= - \int_{\Theta} \boldsymbol{\sigma}^{\mathtt{f}} : \big( \nabla_{\mathtt{y}}\delta\tilde{\boldsymbol{v}} \big)
\,\mathrm{d}\boldsymbol{y}
+ \int_{\partial\Theta} \delta \tilde{\boldsymbol{v}} 
\cdot \big( \boldsymbol{\sigma}^{\mathtt{f}} \cdot \boldsymbol{n} \big) \,\mathrm{d} S,
\end{equation}
$$

由边界条件 $\eqref{eq:fine_scale_bc}$ 可知，上式右侧的边界积分项等于零。弱形式因而可陈述为：

$$
\begin{equation}\label{eq:fine_scale_weak_form}
\int_{\Theta} \boldsymbol{\sigma}^{\mathtt{f}} : \big( \nabla_{\mathtt{y}}\delta\tilde{\boldsymbol{v}} \big)
\,\mathrm{d}\boldsymbol{y} = 0, \quad
\forall \delta \tilde{\boldsymbol{v}} \in H_{\mathtt{per}}^{1}.
\end{equation}
$$

在周期性边界条件设置下，扰动应变场的体积平均等于零，也即

$$
\begin{equation}\label{eq:zero_mean_fluc}
\langle \tilde{\boldsymbol{\varepsilon}}^{\mathtt{f}} \rangle
= \frac{1}{|\Theta|} \int_{\Theta} \nabla_{\mathtt{y}}^{\mathtt{s}} \tilde{\boldsymbol{u}}^{\mathtt{f}}
\,\mathrm{d}\boldsymbol{y}
= \frac{1}{2|\Theta|} \int_{\partial\Theta} \boldsymbol{n}\otimes\tilde{\boldsymbol{u}}^{\mathtt{f}}
+ \tilde{\boldsymbol{u}}^{\mathtt{f}}\otimes\boldsymbol{n} \,\mathrm{d}S = \boldsymbol{0},
\end{equation}
$$

式中，$\boldsymbol{n}$ 是单胞边界的外法向，第二个等号由 Gauss 定理得到，第三个等号由周期性边界条件 $\eqref{eq:fine_scale_bc}$ 得到。因此，对式 $\eqref{eq:fine_scale_cstt}$ 给出的细观尺度应变和应力场作体积平均，就得到

$$
\begin{equation}\label{eq:macro_strain_stress}
\langle \boldsymbol{\varepsilon}^{\mathtt{f}} \rangle = \boldsymbol{\varepsilon}^{\mathtt{c}}
+ \langle \tilde{\boldsymbol{\varepsilon}}^{\mathtt{f}} \rangle
= \boldsymbol{\varepsilon}^{\mathtt{c}}, \quad
\langle \boldsymbol{\sigma}^{\mathtt{f}} \rangle = \boldsymbol{\sigma}^{\mathtt{c}}.
\end{equation}
$$

由式 $\eqref{eq:macro_strain_stress}$ 可以看到，在本论文处理的单胞方程中，宏观应力和应变总是等于细观尺度应力和应变场的体积平均，以后本论文将不再区分宏观应力/应变和平均应力/应变。若给定宏观应变增量 $\Delta \boldsymbol{\varepsilon}^{\mathtt{c}}$，在细观尺度下，容许应变场 $\delta\boldsymbol{\varepsilon}^{\mathtt{f}}\triangleq\delta\boldsymbol{\varepsilon}^{\mathtt{c}}+\nabla_{\mathtt{y}}\delta\tilde{\boldsymbol{v}}$ 所做虚功等于

$$
\begin{equation}\label{eq:hill_mandel_cond}
\begin{aligned}
\langle \boldsymbol{\sigma}^{\mathtt{f}} : \delta\boldsymbol{\varepsilon}^{\mathtt{f}} \rangle
&= \frac{1}{|\Theta|} \Big( \int_{\Theta} \boldsymbol{\sigma}^{\mathtt{f}} \,\mathrm{d}\boldsymbol{y} \Big)
: \delta\boldsymbol{\varepsilon}^{\mathtt{c}}  
+ \frac{1}{|\Theta|} \underbrace{\int_{\Theta} \boldsymbol{\sigma}^{\mathtt{f}} 
: \nabla_{\mathtt{y}}\delta\tilde{\boldsymbol{v}} \,\mathrm{d}\boldsymbol{y}}_{=0,\, \text{ref. }\eqref{eq:fine_scale_weak_form}} \\
&= \boldsymbol{\sigma}^{\mathtt{c}} : \delta\boldsymbol{\varepsilon}^{\mathtt{c}},
\end{aligned}
\end{equation}
$$

这就给出了 Hill--Mandel 条件，即宏观尺度的虚功密度等于细观尺度虚功的体积平均。

## 单胞方程的 Green 函数

虽然上一节给出的单胞方程可直接通过有限元方法计算，但是当如下情景共同出现时，直接求解单胞方程的计算量将是难以承受的：

1. 单胞几何形貌复杂，例如通过 CT 扫描得到的图片、拓扑优化得到的复杂结构、生成的含多个夹杂的 RVE 等。这将使用大量的单元表征复杂夹杂形状的边界。
2. 单胞组分具有材料非线性。这意味着在求解单胞问题时，将使用多个增量步、在每个增量步中多次迭代、在每次迭代的每一个积分点上多次求解非线性方程以更新应力。

此外，观察宏观应力 $\boldsymbol{\sigma}^{\mathtt{c}}$ 的表达式 $\eqref{eq:def_macro_stress}$，其中细观尺度的贡献来自位移场梯度 $\nabla_{\mathtt{y}}^{\mathtt{s}}  \boldsymbol{u}^{(1)}$ 和本征应变 $\boldsymbol{\mu}^{(0)}$，并不直接需要单胞内的位移场分布。这似乎意味着一种直接计算单胞内应变和应力场的数值方法，可能比以节点位移自由度为未知量的有限元方法更加高效。Moulinec 和 Suquet (1994) 提出的 FFT 方法印证了这一观点。在这一小节中将通过弹性力学的 Green 函数表示单胞方程的应变场，读者将会在下一节中看到，方程解的 Green 函数表示为单胞问题的降阶计算提供基础。

> 补充一下最近的 FE2 的文献，以及 FFT 的论述

以下推导过程将使用张量的指标记法，并默认张量的分量在三维空间中笛卡尔坐标系下给出。Green 函数 $G_{k}^{m}$ 是如下方程的解：

$$
\begin{equation}\label{eq:green_func}
L_{ijkl}G_{k,lj}^{m} + \delta_{im} \delta(\boldsymbol{y}-\boldsymbol{y}') = 0,
\end{equation}
$$

式中，$\delta(\boldsymbol{y}-\boldsymbol{y}')$ 是三维空间中的 Dirac 分布，$\delta_{im}$ 是 Kronecker Delta 记号，四阶张量 $L_{ijkl}$ 是弹性刚度张量，具有主对称性和次对称性，也即 $L_{ijkl}=L_{klij}$ 和 $L_{ijkl}=L_{jikl}=L_{ijlk}$。Green 函数 $G_{k}^{m}(\boldsymbol{y},\, \boldsymbol{y}')$ 的物理含义是：作用在点 $\boldsymbol{y}'$，方向为 $m$ 的单位集中力在弹性介质内产生的在点 $\boldsymbol{y}$，方向为 $k$ 的位移分量。需要指出的是，式 $\eqref{eq:green_func}$ 给出的 Green 函数还需要匹配边界条件才能完全定解。特别的，在无穷大各向同性均质弹性区域内，此处定义的 Green 函数又称为 Kelvin 基本解。根据弹性力学的互易原理，Green 函数具有如下互易性：

$$
\begin{equation}\label{eq:recip_green}
G_{k}^{m}(\boldsymbol{y},\,  \boldsymbol{y}') = G_{m}^{k}(\boldsymbol{y}',\,  \boldsymbol{y}).
\end{equation}
$$

即使本论文求解的单胞问题 $\eqref{eq:fine_scale_equili}$ 定解在有限尺寸的非均质单胞 $\Theta$ 内，其 Green 函数仍满足上述性质。不过，只有在无穷大均质弹性区域内（可以是各向异性），Green 函数才具有如下平移不变性

$$
\begin{equation}\label{eq:trans_invar}
G_{k}^{m}(\boldsymbol{y},\,  \boldsymbol{y}') 
= G_{k}^{m}(\boldsymbol{y}+\boldsymbol{a},\,  
\boldsymbol{y}'+\boldsymbol{a}) , \quad \forall \boldsymbol{a}\in \mathbb{R}^{3}
\end{equation}
$$

和中心反演对称性

$$
\begin{equation}\label{eq:central_symm}
G_{k}^{m}(\boldsymbol{y}, \boldsymbol{y}') = G_{m}^{k}(\boldsymbol{y}',\boldsymbol{y}),
\end{equation}
$$

此时才可以将 Green 函数记作 $G_{k}^{m}(\boldsymbol{y},\, \boldsymbol{y}')=\hat{G}_{k}^{m}(\boldsymbol{y} - \boldsymbol{y}')$。现在将细观应力场的表达式 $\eqref{eq:fine_scale_cstt}$ 代入到平衡方程 $\eqref{eq:fine_scale_equili}$ 中，整理之后得到

$$
\begin{equation}\notag
\nabla_{\mathtt{y}}\cdot\big( \mathbb{L}:\tilde{\boldsymbol{\varepsilon}}^{\mathtt{f}} \big)
= - \nabla_{\mathtt{y}}\cdot \boldsymbol{\lambda},\quad
\boldsymbol{\lambda} \triangleq \mathbb{L}:\big( \boldsymbol{\varepsilon}^{\mathtt{c}}  - \boldsymbol{\mu}^{\mathtt{f}} \big)
\end{equation}
$$

式中，$\boldsymbol{\lambda}$ 称为本征应力。应用 Green 函数，单胞内的位移可以表示为

$$
\begin{equation}\label{eq:green_disp}
\begin{aligned}
\tilde{u}_{k}^{\mathtt{f}}(\boldsymbol{y}) 
&= \int_{\Theta}
G_{k}^{m} ( \boldsymbol{y},\,  \boldsymbol{y}' ) \frac{\partial \lambda_{mn}(\boldsymbol{y}')}{\partial y_{n}^\prime} \, \mathrm{d} \boldsymbol{y}' \\
&+ \int_{\partial\Theta} G_{k}^{m} ( \boldsymbol{y},\,  \boldsymbol{y}' ) t_{m}(\boldsymbol{y}') \,\mathrm{d}S'
- \int_{\partial\Theta} \tilde{u}_{m}^{\mathtt{f}}(\boldsymbol{y}' ) T_{k}^{m}(\boldsymbol{y},\, \boldsymbol{y}') \,\mathrm{d}S',
\end{aligned}
\end{equation}
$$

式中，$t_{m}$ 和 $T_{k}^{m}$ 分别是弹性位移解 $\tilde{\boldsymbol{u}}^{\mathtt{f}}$ 和 Green 函数贡献的边界处的力：

$$
\begin{equation}\notag
t_{m} = L_{mnij}\tilde{u}_{i,j}^{\mathtt{f}} n_{n}, \quad
T_{k}^{m} = L_{klij}G_{i,j}^{m} n_{l}
\end{equation}
$$

在周期性边界条件设置下，$T_{k}^{m}|_{\partial\Theta^{-}} = -T_{k}^{m}|_{\partial\Theta^{+}}$，因此关于 $T_{k}^{m}$ 的面积分等于零。对式 $\eqref{eq:green_disp}$ 应用分部积分，得到

$$
\begin{equation}\notag
\begin{aligned}
\tilde{u}_{k}^{\mathtt{f}}(\boldsymbol{y}) 
= &-\int_{\Theta}
G_{k}^{m,n} ( \boldsymbol{y},\, \boldsymbol{y}' ) \lambda_{mn}(\boldsymbol{y}') \, \mathrm{d} \boldsymbol{y}'\\
&+ \int_{\partial\Theta} G_{k}^{m} ( \boldsymbol{y},\, \boldsymbol{y}' )
\big[ t_{m}(\boldsymbol{y}') + \lambda_{mn}(\boldsymbol{y}')n_{n}(\boldsymbol{y}') \big]
\,\mathrm{d}S',\\
&\quad \quad G_{k}^{m,n} ( \boldsymbol{y},\, \boldsymbol{y}' )
\triangleq \frac{1}{2} \Big( \frac{\partial G_{k}^{m} ( \boldsymbol{y},\, \boldsymbol{y}' ) }{\partial y_{n}'} 
+ \frac{\partial G_{k}^{n} ( \boldsymbol{y},\, \boldsymbol{y}' ) }{\partial y_{m}'} \Big),
\end{aligned}
\end{equation}
$$

式中出现的 $G_{k}^{m,n}$ 的定义用到了本征应力的对称性，$\lambda_{mn}=\lambda_{nm}$。注意到上式中边界积分项中 $t_{m}+\lambda_{mn}n_{n}=\sigma_{mn}^{\mathtt{f}}n_{n}$，其中细观应力场满足反周期性条件 $\eqref{eq:fine_scale_bc}$，因此边界项积分后等于零。继续对 $\tilde{\boldsymbol{u}}^{\mathtt{f}}$ 关于坐标 $\boldsymbol{y}$ 求偏导数，并代入本征应力场的表达式，得到单胞内应变扰动场的 Green 函数表示：

$$
\begin{equation}\label{eq:stran_green_repre}
\tilde{\boldsymbol{\varepsilon}}^{\mathtt{f}}(\boldsymbol{y}) 
= -\int_{ \Omega }
\boldsymbol{\Gamma} ( \boldsymbol{y},\, \boldsymbol{y}' ) 
: \mathbb{L}(\boldsymbol{y}')
: \big[ \boldsymbol{\varepsilon}^{\mathtt{c}}  
- \boldsymbol{\mu}^{\mathtt{f}}(\boldsymbol{y}') \big]
\, \mathrm{d} \boldsymbol{y}',
\end{equation}
$$

其中，Green 算子 $\boldsymbol{\Gamma}$（取值为四阶张量）定义为

$$
\begin{equation}\label{eq:def_green_op}
\begin{aligned}
\boldsymbol{\Gamma}&=\Gamma_{kl}^{mn} ( \boldsymbol{y},\, \boldsymbol{y}' ) \\
&\triangleq \frac{1}{4}\Big( 
\frac{\partial^{2}   G_{k}^{m} ( \boldsymbol{y},\, \boldsymbol{y}' ) }{\partial y_{n}'\partial y_{l}} 
+ \frac{\partial^{2} G_{k}^{n} ( \boldsymbol{y},\, \boldsymbol{y}' ) }{\partial y_{m}'\partial y_{l}} 
+ \frac{\partial^{2} G_{l}^{m} ( \boldsymbol{y},\, \boldsymbol{y}' ) }{\partial y_{n}'\partial y_{k}} 
+ \frac{\partial^{2} G_{l}^{n} ( \boldsymbol{y},\, \boldsymbol{y}' ) }{\partial y_{m}'\partial y_{k}} 
\Big).
\end{aligned}
\end{equation}
$$

通过 Green 算子 $\boldsymbol{\Gamma}$，式 $\eqref{eq:stran_green_repre}$ 给出了单胞方程应变扰动场。有如下对 Green 算子 $\boldsymbol{\Gamma}$ 的评注：

1. Green 算子 $\boldsymbol{\Gamma}$ 继承了 Green 函数的互易性，也即 $\Gamma_{kl}^{mn}(\boldsymbol{y},\, \boldsymbol{y}')=\Gamma_{mn}^{kl}(\boldsymbol{y}',\, \boldsymbol{y})$。
2. Green 算子 $\boldsymbol{\Gamma}$ 具有次对称性，也即 $\Gamma_{kl}^{mn}=\Gamma_{lk}^{mn}=\Gamma_{kl}^{nm}$，但因为 Green 函数可能不具备平移不变性和中心反演对称性，因此一般不具有主对称性，也即 $\Gamma_{kl}^{mn} \neq \Gamma_{mn}^{kl}$。
3. Green 算子 $\boldsymbol{\Gamma}$ 与（1）定解区域 $\Theta$、（2）弹性刚度张量分布函数 $\mathbb{L}(\boldsymbol{y})$、（3）边界条件 $\mathcal{B}$ 三者密切相关。在第 $\ref{}$ 节讨论不同降阶方法的区别时将详细说明。

## 场变换算子

式 $\eqref{eq:stran_green_repre}$ 事实上给出了本征应变场到应变扰动场的映射，因此定义如下场变换算子 $\mathcal{D}$，将本征应变场 $\boldsymbol{\mu}$ 变换为单胞域内的影响函数 $\mathcal{D}*\boldsymbol{\mu}$：

$$
\begin{equation}\label{eq:def_trans_op}
\mathcal{D}*\boldsymbol{\mu}(\boldsymbol{y}) 
\triangleq \int_{\Theta} \boldsymbol{\Gamma}(\boldsymbol{y},\, \boldsymbol{y}') 
: \mathbb{L}(\boldsymbol{y}') : \boldsymbol{\mu}(\boldsymbol{y}')\,\mathrm{d}\boldsymbol{y},
\end{equation}
$$

那么，单胞内的应变场 $\eqref{eq:stran_green_repre}$ 可以用场变换算子重新表示为

$$
\begin{equation}\label{eq:fine_scale_strain_fluc_transop}
\boldsymbol{\varepsilon}^{\mathtt{f}}
= \boldsymbol{\varepsilon}^{\mathtt{c}}
+ \tilde{\boldsymbol{\varepsilon}}^{\mathtt{f}}
= \boldsymbol{\varepsilon}^{\mathtt{c}}-\mathcal{D}*\boldsymbol{\varepsilon}^{\mathtt{c}}
+ \mathcal{D}*\boldsymbol{\mu}^{\mathtt{f}}.
\end{equation}
$$

通过定义式 $\eqref{eq:def_trans_op}$ 以及 Green 算子 $\boldsymbol{\Gamma}$ 的性质，可得到场变换算子具有如下性质：

1. 算子 $\mathcal{D}$ 是线性算子。

2. 对任意的本征应变场 $\boldsymbol{\mu}$，影响函数 $\mathcal{D}*\boldsymbol{\mu}$ 是一个相容的应变扰动场，也即

$$
\begin{equation}\label{eq:transop_zero_mean}
\langle \mathcal{D}*\boldsymbol{\mu} \rangle = \boldsymbol{0}.
\end{equation}
$$

3. 对任意相容的应变扰动场 $\tilde{\boldsymbol{\varepsilon}}$，算子 $\mathcal{D}$ 将该应变场映射至自身，
   
$$
\begin{equation}\label{eq:transop_fixed_point}
\mathcal{D}*\tilde{\boldsymbol{\varepsilon}} = \tilde{\boldsymbol{\varepsilon}}.
\end{equation}
$$

其中第三条性质的证明可参考 Fish et al. (2014)。通过场变换根据定义 $\eqref{eq:def_trans_op}$，以及算子 $\mathcal{D}$ 的线性，式 $\eqref{eq:fine_scale_strain_fluc_transop}$ 由宏观应变 $\boldsymbol{\varepsilon}^{\mathtt{c}}$ 贡献的影响函数可分解为

$$
\begin{equation}\label{eq:def_elastic_infl_func}
\begin{aligned}
\boldsymbol{\varepsilon}^{\mathtt{c}}-\mathcal{D}*\boldsymbol{\varepsilon}^{\mathtt{c}}
&= \boldsymbol{\varepsilon}^{\mathtt{c}}-\int_{\Theta} \boldsymbol{\Gamma}(\boldsymbol{y},\, \boldsymbol{y}') 
: \mathbb{L}(\boldsymbol{y}') : \boldsymbol{\varepsilon}^{\mathtt{c}}\,\mathrm{d}\boldsymbol{y} \\
&= \mathbb{I}(\boldsymbol{y}) : \boldsymbol{\varepsilon}^{\mathtt{c}}-\Big( \int_{\Theta} \boldsymbol{\Gamma}(\boldsymbol{y},\, \boldsymbol{y}') 
: \mathbb{L}(\boldsymbol{y}') : \mathbb{I}(\boldsymbol{y}') \,\mathrm{d}\boldsymbol{y} \Big)
: \boldsymbol{\varepsilon}^{\mathtt{c}} \\
&\triangleq \mathbb{E}(\boldsymbol{y}) : \boldsymbol{\varepsilon}^{\mathtt{c}},
\end{aligned}
\end{equation}
$$

式中，$\mathbb{I}(\boldsymbol{y})$ 是在单胞区域 $\Theta$ 内取值为四阶单位张量的常值函数；四阶单位张量定义为 $(\mathbb{I})_{ij}^{kl}\triangleq\frac{1}{2}(\delta_{ik}\delta_{jl}+\delta_{il}\delta_{jk})$，具有主对称性和次对称性；四阶张量场 $\mathbb{E}(\boldsymbol{y})$ 称为弹性应变影响函数。当不考虑材料非线性时，也即本征应变恒等于零，根据式 $\eqref{eq:fine_scale_cstt}$、$\eqref{eq:fine_scale_strain_fluc_transop}$、$\eqref{eq:def_elastic_infl_func}$，细观尺度的应变和应力场分别等于

$$
\begin{equation}\label{eq:fins_scale_strain_stress_infl}
\boldsymbol{\varepsilon}^{\mathtt{f}}(\boldsymbol{y})
= \mathbb{E}(\boldsymbol{y}) : \boldsymbol{\varepsilon}^{\mathtt{c}},\quad
\boldsymbol{\sigma}^{\mathtt{f}}(\boldsymbol{y})
= \mathbb{L}(\boldsymbol{y}) : \mathbb{E}(\boldsymbol{y}) : \boldsymbol{\varepsilon}^{\mathtt{c}},
\end{equation}
$$

在体积平均之后，宏观应力与宏观应变之间的关系为

$$
\begin{equation}\label{eq:def_macro_elastic_stiffness_tensor}
\boldsymbol{\sigma}^{\mathtt{c}}
= \langle \mathbb{L}:\mathbb{E} \rangle : \boldsymbol{\varepsilon}^{\mathtt{c}}
\triangleq \mathbb{L}^{\mathtt{c}} : \boldsymbol{\varepsilon}^{\mathtt{c}},
\end{equation}
$$

式中，四阶张量 $\mathbb{L}^{\mathtt{c}}$ 称为宏观弹性刚度张量，具有主对称性、次对称性和正定性（证明见 Fish (2014)），通过弹性应变影响函数和材料相的弹性刚度张量获得：

$$
\begin{equation}\label{eq:elastic_infl_tensor}
\begin{aligned}
\mathbb{L}^{\mathtt{c}} 
&= \frac{1}{|\Theta|} \sum_{\alpha=1}^{N_{\mathtt{phs}}} 
\int_{\Theta^{(\alpha)}} \mathbb{L}^{(\alpha)} : \mathbb{E}(\boldsymbol{y}) \,\mathrm{d}\boldsymbol{y}\\
&= \sum_{\alpha=1}^{N_{\mathtt{phs}}} c^{(\alpha)} \mathbb{L}^{(\alpha)} : \mathbb{E}^{(\alpha)},\quad
\mathbb{E}^{(\alpha)} \triangleq \langle \mathbb{E} \rangle^{(\alpha)},
\end{aligned}
\end{equation}
$$

式中，运算符 $\langle \bullet \rangle^{(\alpha)}$ 表示在材料相 $\alpha$ 内的体积平均；四阶张量 $\mathbb{E}^{(\alpha)}$ 称为材料相 $\alpha$ 的弹性应变影响张量，由影响函数 $\mathbb{E}$ 在材料相 $\alpha$ 内体积平均后得到。

## 分片常值分布的本征应变场

在单胞内应变场表达式 $\eqref{eq:fine_scale_strain_fluc_transop}$ 中，由本征应变场贡献的部分为 $\mathcal{D}*\boldsymbol{\mu}^{\mathtt{f}}$。然而表征非弹性变形的本征应变场一般在单胞域内并不均匀，无法直接给出类似于式 $\eqref{eq:def_elastic_infl_func}$ “细观影响函数-宏观物理量”的分离变量形式。类似有限元的思路，可以定义一组有限数量的形函数 $\boldsymbol{\mu}^{(k)}(\boldsymbol{y})$，$k=1,2,\ldots,N_{\mathtt{sp}}$，以及与形函数匹配的“自由度” $\xi^{k}$，并假设本征应变场可写成

$$
\begin{equation}\label{eq:ansatz_eigenstrain}
\boldsymbol{\mu}^{\mathtt{f}}(\boldsymbol{y})
= \sum_{k=1}^{N_{\mathtt{sp}}} \xi^{(k)} \boldsymbol{\mu}^{(k)}(\boldsymbol{y}),
\end{equation}
$$

在场变换算子作用下，可得本征应变场贡献的单胞内应变扰动场为

$$
\mathcal{D}*\boldsymbol{\mu}^{\mathtt{f}}
= \sum_{k=1}^{N_{\mathtt{sp}}} \xi^{(k)} \big( \mathcal{D}*\boldsymbol{\mu}^{(k)} \big),
$$

上式应用到场变换算子的线性性质，并在形式上得到了类似于式 $\eqref{eq:def_elastic_infl_func}$ 的宏细观变量分离的形式。然而，紧接而来的问题是：（1）如何选取形函数 $\boldsymbol{\mu}^{(k)}$？以及（2）如何计算系数 $\xi^{(k)}$？Dvorak 等人提出的均匀场变换分析（TFA），针对上述问题作出如下假设：

1. 假设本征应变场分片常值分布。若设本征应变场在同一材料相上均匀分布（One Phase One Partition，OPOP），那么有

$$
\begin{equation}\label{eq:eigenstrain_opop}
\boldsymbol{\mu}^{\mathtt{f}}(\boldsymbol{y}) 
= \sum_{\alpha=1}^{N_{\mathtt{phs}}} \mathbb{I}^{(\alpha)}(\boldsymbol{y}) : \boldsymbol{\mu}^{(\alpha)},\quad
\mathbb{I}^{(\alpha)}(\boldsymbol{y})\triangleq\chi^{(\alpha)}(\boldsymbol{y})\mathbb{I},
\end{equation}
$$

​	式中，$\boldsymbol{\mu}^{(\alpha)}$ 是分块 $\alpha$ 的本征应变。将假设的本征应变式 $\eqref{eq:eigenstrain_opop}$ 代入到单胞内应变场表达式 $\eqref{eq:fine_scale_strain_fluc_transop}$ 中，得到

$$
\begin{equation}\label{eq:def_eigenstrain_infl_func}
  \boldsymbol{\varepsilon}^{\mathtt{f}}(\boldsymbol{y})
  = \mathbb{E}(\boldsymbol{y}) : \boldsymbol{\varepsilon}^{\mathtt{c}}
  + \sum_{\alpha=1}^{N_{\mathtt{phs}}} \mathbb{P}^{(\alpha)}(\boldsymbol{y})
  : \boldsymbol{\mu}^{(\alpha)},\quad
  \mathbb{P}^{(\alpha)}(\boldsymbol{y})
  \triangleq \mathcal{D}*\mathbb{I}^{(\alpha)}(\boldsymbol{y}),
  \end{equation}
$$

​	式中，$\mathbb{P}^{(\alpha)}$ 是分块 $\alpha$ 的本征应变影响函数。

2. 假设分块 $\alpha$ 的本征应变由同一分块的**平均应变**驱动。对细观应变场 $\eqref{eq:def_eigenstrain_infl_func}$ 在每一个分块 $\alpha$ 作体积平均，得到

$$
\begin{equation}\label{eq:reduced_eqs}
\begin{aligned}
\boldsymbol{\varepsilon}^{(\beta)}
&= \mathbb{E}^{(\beta)}:\boldsymbol{\varepsilon}^{\mathtt{c}}
+ \sum_{\alpha=1}^{N_{\mathtt{phs}}} \mathbb{P}^{(\beta\alpha)} : \boldsymbol{\mu}^{(\alpha)}, \quad
\beta=1,2,\ldots,N_{\mathtt{phs}},\\
\boldsymbol{\varepsilon}^{(\beta)}
&\triangleq\langle \boldsymbol{\varepsilon}^{\mathtt{f}} \rangle^{(\beta)}, \quad
\mathbb{E}^{(\beta)}\triangleq\langle \mathbb{E} \rangle^{(\beta)}, \quad
\mathbb{P}^{(\beta\alpha)} \triangleq \langle \mathbb{P}^{(\alpha)} \rangle^{(\beta)},
\end{aligned}
\end{equation}
$$

​	式中，$\boldsymbol{\varepsilon}^{(\beta)}$ 是分块 $\beta$ 的平均应变；四阶张量 $\mathbb{E}^{(\beta)}$ 是分块 $\beta$ 的弹性应变影响张量；四阶张量 $\mathbb{P}^{(\beta\alpha)}$ 是分块 $\alpha$ 对分块 $\beta$ 的本征应变影响张量。本征应变 $\boldsymbol{\mu}^{(\alpha)}$ 依赖于分块应变 $\boldsymbol{\varepsilon}^{(\alpha)}$ 的加载历史，因此式 $\eqref{eq:reduced_eqs}$ 得到的是关于分块应变（增量）的非线性方程组。之后将在 ROH 线上阶段这一节中详细说明。

> 对第二点假设，给出 Suquet 文献的引用。

TFA 允许进一步将各个材料相区域 $\Theta^{(\alpha)}$ 进一步分成若干个分块 $\Theta^{(\alpha_{k})}$，$k=1,2,\ldots,N_{\mathtt{ptn}}^{(\alpha)}$，不同分块之间满足

$$
\Theta^{(\alpha_{k})} \cap \Theta^{(\alpha_{j})} = \varnothing, \quad \forall \alpha_{k} \neq \alpha_{j}
\quad \text{and} \quad 
\bigcup_{k=1}^{N_{\mathrm{ptn}}^{(\alpha)}} \Theta^{(\alpha_{k})} = \Theta^{(\alpha)},
$$

式中，$\alpha_{k}$ 是材料相 $\alpha$ 的分块序号，$N_{\mathtt{ptn}}^{(\alpha)}$ 是材料相 $\alpha$ 的分块数量。由此，假设在分块上分片常值分布的本征应变场表示为（One Phase Multiple Partitions，OPMP）

$$
\begin{equation}\label{eq:C0124__eigstran_tfa_opmp}
\boldsymbol{\mu}^{\mathtt{f}}(\boldsymbol{y}) = \sum_{\alpha=1}^{N_{\mathtt{phs}}}\sum_{k=1}^{N_{\mathtt{pth}}^{(\alpha)}}
\mathbb{I}^{(\alpha_{k})}(\boldsymbol{y}) : \boldsymbol{\mu}^{(\alpha_{k})},\quad
\mathbb{I}^{(\alpha_{k})}(\boldsymbol{y})\triangleq\chi^{(\alpha_{k})}(\boldsymbol{y})\mathbb{I},
\end{equation}
$$

式中，$\chi^{(\alpha_{k})}$ 是分块 $\alpha_{k}$ 的特征函数，$\boldsymbol{\mu}^{(\alpha_{k})}$ 是分块 $\alpha_{k}$ 的本征应变。若单胞使用有限元网格离散，每个实体单元所占区域用 $\Theta^{(e)}$ 表示，其中 $e$ 是单元序号，那么 TFA 最精细的分块划分方案，是将每个实体单元视作独立的分块（One Element One Partition，OEOP），此时本征应变场表示为

$$
\begin{equation}\label{eq:C0124__eigstran_tfa_oeop}
\boldsymbol{\mu}^{\mathtt{f}}(\boldsymbol{y}) = \sum_{e=1}^{N_{\mathtt{elem}}} \mathbb{I}^{(e)}(\boldsymbol{y}) : \boldsymbol{\mu}^{(e)},\quad
\mathbb{I}^{(e)}(\boldsymbol{y})\triangleq\chi^{(e)}(\boldsymbol{y})\mathbb{I},
\end{equation}
$$

式中，$\chi^{(e)}$ 是单元 $e$ 的特征函数，$N_{\mathtt{elem}}$ 单胞网格离散后的单元总数，$\boldsymbol{\mu}^{(e)}$ 是单元的本征应变。本论文将在第三章证明，当 ROH 使用 OEOP 分块方案时，其数值计算结果和使用常应变单元的有限元方法是相同的。
