# Eshelby 张量的表示方法

这一部分将给出弹性介质内的 Green 函数的定义，并用于表示 Eshelby 张量。需注意，此处只给出**形式上**的 Green 函数，弹性介质可以是**任意的各向异性**，夹杂也可以是**任意形状**。

## Green 函数的定义与性质

考虑三维空间 $\mathbb{R}^{3}$ 内的弹性力学方程

$$
\boldsymbol{\sigma} \cdot \nabla + \boldsymbol{f} = \boldsymbol{0}, \quad 
\boldsymbol{\sigma} = \mathbb{L}: \boldsymbol{\varepsilon}, \quad
\boldsymbol{\varepsilon} = \frac{1}{2} ( \nabla \boldsymbol{u} + \boldsymbol{u} \nabla ),
$$

式中，$\mathbb{L}$ 是弹性刚度张量，$\boldsymbol{f}$ 是体力项。上式又可以写成分量形式为

$$
\begin{equation}
L_{ijkl}\ u_{k,lj} + f_{i} = 0.
\label{eq:elastic}
\end{equation}
$$

方程的解可通过 Green 函数[^1] $G_{k}^{m}$ 表示成与体力项 $f_{i}$ 的卷积：

$$
\begin{equation}
u_{k}(\boldsymbol{x}) = \int_{\mathbb{R}^{3}}
G_{k}^{m} ( \boldsymbol{x} - \boldsymbol{x}' ) f_{m}(\boldsymbol{x}') \ \mathrm{d} \boldsymbol{x}',
\label{eq:sln_green}
\end{equation}
$$

式中使用了 Einstein 求和约定。Green 函数是如下方程的解：

$$
\begin{equation}
L_{ijkl}G_{k,lj}^{m} + \delta_{im} \delta(\boldsymbol{x}-\boldsymbol{x}') = 0,
\label{eq:green}
\end{equation}
$$

式中，$\delta(\boldsymbol{x}-\boldsymbol{x}')$ 是 Dirac 函数。因此，Green 函数 $G_{k}^{m}(\boldsymbol{x}-\boldsymbol{x}')$ 的物理含义是：作用在点 $\boldsymbol{x}'$，方向为 $m$ 的单位集中力在弹性介质内产生的在点 $\boldsymbol{x}$，方向为 $k$ 的位移分量。Green 函数有着如下性质：

1. 平移不变性：

$$
\begin{equation}
G_{k}^{m}(\boldsymbol{x}, \boldsymbol{x}') = G_{k}^{m}(\boldsymbol{x}-\boldsymbol{x}') = G_{k}^{m}(\boldsymbol{x}',\boldsymbol{x}).
\label{eq:trans_invar}
\end{equation}
$$


2. 互易原理：

$$
\begin{equation}
G_{k}^{m}(\boldsymbol{x}, \boldsymbol{x}') = G_{m}^{k}(\boldsymbol{x}', \boldsymbol{x}) = G_{m}^{k}(\boldsymbol{x}, \boldsymbol{x}').
\label{eq:recip}
\end{equation}
$$

## Eshelby 张量

考虑如下夹杂问题：

$$
\boldsymbol{\sigma} \cdot \nabla= \boldsymbol{0}, \quad 
\boldsymbol{\sigma} = \mathbb{L}: (\boldsymbol{\varepsilon} - \chi_{\Omega} \boldsymbol{\mu})
= \mathbb{L}: \boldsymbol{\varepsilon} + \chi_{\Omega}\boldsymbol{\lambda},
$$

式中，$\boldsymbol{\mu}$，$\boldsymbol{\lambda}$ 分别是夹杂 $\Omega$ 内的本征应变和本征应力场，示性函数 $\chi_{\Omega}$ 表示本征应变和本征应力场在夹杂外取值为零。夹杂的本征应力和本征应变可以视作是方程 $\eqref{eq:elastic}$ 中的体力项：

$$
\boldsymbol{f} = - (\mathbb{L}:\chi_{\Omega}\boldsymbol{\mu}) \cdot\nabla = \chi_{\Omega}\boldsymbol{\lambda} \cdot\nabla,
$$

或者写成指标形式：

$$
f_{m} = -\chi_{\Omega} \frac{\partial }{\partial x_{n}}(L_{mnij} \mu_{ij}) 
= \chi_{\Omega}\frac{\partial \lambda_{mn}}{\partial x_{n}}
$$

代入式 $\eqref{eq:sln_green}$ 中，由于示性函数 $\chi_{\Omega}$ 限制积分项只在夹杂内取值，所以积分区域从 $\mathbb{R}^{3}$ 转换为 $\Omega$：

$$
\begin{aligned}
u_{k}(\boldsymbol{x}) 
&= -\int_{\Omega}
G_{k}^{m} ( \boldsymbol{x} - \boldsymbol{x}' ) 
\frac{\partial }{\partial x_{n}^\prime}\big[ L_{mnij} \mu_{ij}(\boldsymbol{x}') \big] \ \mathrm{d} \boldsymbol{x}' \\
&= \int_{\Omega}
G_{k}^{m} ( \boldsymbol{x} - \boldsymbol{x}' ) \frac{\partial }{\partial x_{n}^\prime}\lambda_{mn}(\boldsymbol{x}') \ \mathrm{d} \boldsymbol{x}',
\end{aligned}
$$

对上式应用分部积分公式（严谨的分析应考虑偏导数在夹杂边界的间断，这里略过），并根据平移不变性，将 Green 函数关于 $\boldsymbol{x}'$ 的偏导数转换为 $\partial G_{k}^{m}/ \partial x_{n}' = -\partial G_{k}^{m}/ \partial x_{n}$：

$$
\begin{aligned}
u_{k}(\boldsymbol{x}) 
&= -\int_{ \Omega }
G_{k,n}^{m} ( \boldsymbol{x} - \boldsymbol{x}' ) \big[ L_{mnij} \mu_{ij}(\boldsymbol{x}') \big] \ \mathrm{d} \boldsymbol{x}' \\
&= \int_{ \Omega }
G_{k,n}^{m} ( \boldsymbol{x} - \boldsymbol{x}' ) \lambda_{mn}(\boldsymbol{x}') \ \mathrm{d} \boldsymbol{x}'.
\end{aligned}
$$

继续对位移场求偏导数得到应变场，根据二阶偏导数和指标 $m,n$ 的可交换性，有 $G_{k,nl}^{m}=G_{m,nl}^{k} = G_{k,ln}^{m}$，可以定义如下**具有次对称性**的四阶张量算子 $\Gamma_{mnkl}$：

$$
\Gamma_{mnkl}
= -\frac{1}{4} \big( G_{k,nl}^{m} + G_{l,nk}^{m} + G_{k,ml}^{n} + G_{l,mk}^{n} \big),
$$

注意到根据互易原理 $\eqref{eq:recip}$，四阶张量 $\Gamma_{mnij}$ 同时还**具有主对称性**，$\Gamma_{mnij} = \Gamma_{ijmn}$。最终，应变场可以表示为

$$
\begin{aligned}
\varepsilon_{kl}(\boldsymbol{x}) 
&= \int_{ \Omega }
\Gamma_{klmn} ( \boldsymbol{x} - \boldsymbol{x}' ) L_{mnij} \mu_{ij}(\boldsymbol{x}') \ \mathrm{d} \boldsymbol{x}' \\
&= -\int_{ \Omega }
\Gamma_{klmn} ( \boldsymbol{x} - \boldsymbol{x}' ) \lambda_{mn}(\boldsymbol{x}') \ \mathrm{d} \boldsymbol{x}'.
\end{aligned}
$$

若夹杂内的本征应变或本征应力是**均匀的**，那么就可以将 $\mu_{ij}$ 或 $\lambda_{ij}$ 提出到积分之外，得到用四阶张量 $S_{ijkl}$ 或 $P_{ijkl}$ [^2]表示的应变场：

$$
\begin{gathered}
\varepsilon_{ij}(\boldsymbol{x}) = S_{ijkl}(\boldsymbol{x}) \mu_{kl}
= -P_{ijkl}(\boldsymbol{x}) \lambda_{kl}, \newline
S_{ijkl}(\boldsymbol{x})
= \int_{\Omega} \Gamma_{ijmn} (\boldsymbol{x} - \boldsymbol{x}') L_{mnkl} \ \mathrm{d} \boldsymbol{x}'
= P_{ijmn}(\boldsymbol{x}) L_{mnkl}.
\end{gathered}
$$

特别的，当夹杂形状为椭球，并且点 $\boldsymbol{x}$ 在椭球内部时，此时 $S_{ijkl}$ 和 $P_{ijkl}$ 为**常数**，也即 **Eshelby 张量**和 **Hill 极化张量**。这一结论对任意各向异性的弹性介质均成立。

## 更新日志

### 2026/07/24

1. 创建了文档 `Eshelby 张量的表示方法.md`

[^1]: 数学上严谨的叫法应该是基本解（fundamental solution），Green 函数一般是需要考虑边界条件的。
[^2]: Eshelby 张量 $S_{ijkl}$ 具有次对称性，但一般不具有主对称性；而 Hill 张量 $P_{ijkl}$ 同时具有主对称性和次对称性。

