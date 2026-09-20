# Papkovich--Neuber 通解求解椭球夹杂问题

在空间中某一点 $\boldsymbol{\xi}=(\xi,\eta,\zeta)$ 中有一质量为 $M$ 的质点，对空间中另外一点处 $\boldsymbol{x} = (x,y,z)$ 的单位质量的引力可以表示为

$$
\boldsymbol{F} = -\frac{GM}{\| \boldsymbol{r} \|^{2}}\hat{\boldsymbol{r}}, \quad 
\boldsymbol{r} \triangleq \boldsymbol{x} - \boldsymbol{\xi}, \quad
\hat{\boldsymbol{r}} = \boldsymbol{r}/ \| \boldsymbol{r} \|.
$$

如果定义引力势函数 $\varphi$

$$
\varphi(\boldsymbol{x}) = -GM \frac{1}{\|\boldsymbol{r}\|},
$$

可以验证，势函数 $\varphi$ 的梯度等于负引力场：

$$
\nabla \varphi = -\boldsymbol{F}.
$$

上式可通过如下对各向同性函数的微分运算公式得到：

$$
\nabla \|\boldsymbol{r}\|^{-n} = -n \|\boldsymbol{r}\|^{-(n+1)}\hat{\boldsymbol{r}}, \quad
\nabla \cdot ( \nabla \|\boldsymbol{r}\|^{-n} )
= \nabla^{2} \|\boldsymbol{r}\|^{-n}
= n(n-1) \| \boldsymbol{r} \|^{-(n+2)}.
$$

值得注意的是，当公式中 $n=1$ 时，$\nabla^{2} \|\boldsymbol{r}\|^{-1} = 0$，因此 $\|\boldsymbol{r}\|^{-1}$ 是在除点 $\boldsymbol{\xi}$ 之外的调和函数。Poisson 方程的基本解：

$$
\nabla^{2} G = \delta(\boldsymbol{x} - \boldsymbol{y}),
$$

在各维度中由如下公式给出：

$$
G(\boldsymbol{x} - \boldsymbol{y}) = \begin{cases}
\dfrac{1}{2}|x - y|, & d=1, \\
\dfrac{1}{2\pi} \ln\|\boldsymbol{x} - \boldsymbol{y}\|, & d=2, \\
\dfrac{1}{(2-d) \omega_{d-1}} \| \boldsymbol{x}-\boldsymbol{y} \|^{2-d}, & d\geq 3.
\end{cases}
$$

式中，$\omega_{d-1}$ 是在 $d-1$ 维空间中单位球的表面积公式。当 $d=3$ 时，方程基本解又称为 Newton 位势 $N$：

$$
N(\boldsymbol{x}-\boldsymbol{y}) = -\frac{1}{4\pi^{2}} \| \boldsymbol{x}-\boldsymbol{y} \|^{-1}.
$$

因此，对任意三维空间中的函数 $f$，Poisson 方程 $\nabla^{2} u = f$ 的解就可以通过 Newton 位势表示为

$$
u(\boldsymbol{x}) = \int_{V} N(\boldsymbol{x} - \boldsymbol{y}) f(\boldsymbol{y}) \ \mathrm{d} \boldsymbol{y}
= \frac{1}{4\pi} \int_{V} \frac{f(\boldsymbol{y})}{\| \boldsymbol{x}-\boldsymbol{y} \|} \ \mathrm{d} \boldsymbol{y},
$$

可以将以 Newton 位势为 kernal 的积分作为 Poisson 方程的逆算子，定义为 $\mathcal{F}$：

$$
\mathcal{F} (f) = u, \quad where \quad \nabla^{2} u = f.
$$


Navier 方程，也即各向同性弹性力学方程可以表示为

$$
\nabla^{2} \boldsymbol{u} + \frac{1}{1-2\nu} \nabla ( \nabla \cdot \boldsymbol{u} ) + \frac{1}{\mu} \boldsymbol{f} = \boldsymbol{0}, \quad
or \quad \mathcal{L} \boldsymbol{u} + \frac{1}{\mu} \boldsymbol{f} = \boldsymbol{0}.
$$

现考虑无体力项的弹性力学方程 $\mathcal{L} \boldsymbol{u}  = \boldsymbol{0}$，可以验证，如果有两个调和函数 $\boldsymbol{P}$ 和 $P_{0}$，那么可以直接构造出满足该方程的解 $\hat{\boldsymbol{u}}$：

$$
\begin{cases}
\hat{\boldsymbol{u}} = \boldsymbol{P} - \dfrac{1}{4(1-\nu)}\nabla( P_{0} + \boldsymbol{r}\cdot\boldsymbol{P} ), \\
\nabla^{2} \boldsymbol{P} = \boldsymbol{0}, \quad \nabla^{2}P_{0} = 0.
\end{cases}
$$

其中 $\boldsymbol{r}$ 是位置矢量。以上通解形式又称 Papkovich--Neuber（P-N）通解。以下验证函数 $\hat{\boldsymbol{u}}$ 确实是无体力弹性力学方程的解，首先将解代入到方程中得到

$$
\begin{aligned}
\nabla^{2} \hat{\boldsymbol{u}}
&=  - \dfrac{1}{4(1-\nu)} \nabla\nabla^{2}(\boldsymbol{r}\cdot\boldsymbol{P}), \\
\nabla ( \nabla \cdot \hat{\boldsymbol{u}} )
&= \nabla ( \nabla \cdot \boldsymbol{P} )
- \dfrac{1}{4(1-\nu)} \nabla\nabla^{2}(\boldsymbol{r}\cdot\boldsymbol{P}),
\end{aligned}
$$

再将恒等式

$$
\nabla^{2}( \boldsymbol{r}\cdot\boldsymbol{P} )
= \boldsymbol{r}\cdot\nabla^{2}\boldsymbol{P} + 2\nabla\cdot\boldsymbol{P}
$$

代入得到

$$
\nabla^{2} \hat{\boldsymbol{u}}
= - \dfrac{1}{2(1-\nu)} \nabla (\nabla\cdot\boldsymbol{P}),\quad
\nabla ( \nabla \cdot \hat{\boldsymbol{u}} )
= \dfrac{1-2\nu}{2(1-\nu)} \nabla (\nabla\cdot\boldsymbol{P}),
$$

这就可以看到，构造的解 $\hat{\boldsymbol{u}}$ 满足

$$
\nabla^{2} \hat{\boldsymbol{u}} + \frac{1}{1-2\nu} \nabla ( \nabla \cdot \hat{\boldsymbol{u}} )
= \boldsymbol{0} = \mathcal{L} \hat{\boldsymbol{u}}.
$$

另一方面，解的完备性要求，是否对任意满足 $\mathcal{L} \boldsymbol{u} = \boldsymbol{0}$ 方程的解，都可以按照式用两个调和函数进行表示？可以验证，如果 $\boldsymbol{u}$ 是无体力弹性力学方程的解，那么可以构造如下两个调和函数：

$$
\hat{\boldsymbol{P}} = \boldsymbol{u} + \frac{1}{1-2\nu} \nabla ( \nabla \cdot \mathcal{F}\boldsymbol{u} ), \quad
\hat{P}_{0}  = \frac{4(1-\nu)}{1-2\nu} \nabla \cdot \mathcal{F}\boldsymbol{u} - \boldsymbol{r}\cdot\hat{\boldsymbol{P}},
$$

可以验证，这两种构造确实是调和函数，因为

$$
\nabla^{2} \hat{\boldsymbol{P}}
= \nabla^{2} \boldsymbol{u} 
+ \frac{1}{1-2\nu} \nabla( \nabla \cdot \boldsymbol{u} ) = \mathcal{L} \boldsymbol{u} = \boldsymbol{0}, \\
\nabla^{2} \hat{P}_{0}
= \frac{4(1-\nu)}{1-2\nu} \nabla \cdot \boldsymbol{u} - 2\nabla\cdot\boldsymbol{P}
= \frac{4(1-\nu)}{1-2\nu} \nabla \cdot \boldsymbol{u} - 2\nabla\cdot\boldsymbol{u}
- \frac{2}{1-2\nu} \nabla \cdot \boldsymbol{u} = 0.
$$

并且满足通解表达式。

另外一种通解的表达方式，可以只借助一个矢量调和函数 $\boldsymbol{F}$：

$$
\boldsymbol{u} = \boldsymbol{F} - \frac{1}{2(1-\nu)} \nabla \big( \mathcal{F}(\nabla\cdot\boldsymbol{F}) \big)
$$

以下通过 P-N 通解，给出含体力项的 Kelvin 特解形式。此时为了定解，函数 $P_{0}$ 和 $\boldsymbol{P}$ 不再是调和函数，但满足 Poisson 方程，因此可以通过 Newton 位势给出。还是假设特解 $\boldsymbol{u}$ 的形式同式，代入到含体力项的 Navier 方程中得到：

$$
\nabla^{2}\boldsymbol{P} 
- \frac{1}{2(1-2\nu)} \nabla ( \nabla^{2}P_{0} + \boldsymbol{r}\cdot\nabla^{2}\boldsymbol{P} )
+ \frac{1}{\mu} \boldsymbol{f} = \boldsymbol{0}.
$$

可以看到，如果函数 $\boldsymbol{P}$ 和 $P_{0}$ 分别满足如下 Poisson 方程：

$$
\nabla^{2}\boldsymbol{P} = - \frac{1}{\mu} \boldsymbol{f}, \quad
\nabla^{2} P_{0} = \frac{1}{\mu} \boldsymbol{r} \cdot \boldsymbol{f},
$$

那么含体力项的 Navier 方程自动满足。可以通过 Newton 位势得到函数 $\boldsymbol{P}$ 和 $P_{0}$ 的表示：

$$
\boldsymbol{P} = -\frac{1}{\mu} \mathcal{F}(\boldsymbol{f}), \quad
P_{0}  =  \frac{1}{\mu} \mathcal{F}(\boldsymbol{r}\cdot\boldsymbol{f}), \quad
$$

代回位移场中得到

$$
\boldsymbol{u}= -\frac{1}{\mu} \mathcal{F}(\boldsymbol{f}) 
+ \dfrac{1}{4(1-\nu)\mu}\nabla \Big( \boldsymbol{r}\cdot \mathcal{F}(\boldsymbol{f}) - \mathcal{F}(\boldsymbol{r}\cdot\boldsymbol{f}) \Big)
$$

将 Newton 位势的具体表达式代入得到

$$
\boldsymbol{u} = \frac{3-4\nu}{16\pi(1-\nu)\mu}
\int_{V} \frac{\boldsymbol{f}(\boldsymbol{y})}{\| \boldsymbol{x}-\boldsymbol{y} \|} \ \mathrm{d} \boldsymbol{y}
+ \dfrac{1}{16\pi(1-\nu)\mu}
\int_{V}\frac{\boldsymbol{x}-\boldsymbol{y}}{\| \boldsymbol{x}-\boldsymbol{y} \|^{3}} \boldsymbol{f}(\boldsymbol{y})\cdot(\boldsymbol{x}-\boldsymbol{y}) \ \mathrm{d} \boldsymbol{y}
$$

为了使记号更加紧凑，定义

$$
\boldsymbol{\rho} = \boldsymbol{x} - \boldsymbol{y}, \quad \rho = \| \boldsymbol{x} - \boldsymbol{y} \|,
$$

代入方程中，并隐藏函数的自变量，得到

$$
\boldsymbol{u}(\boldsymbol{x}) = \frac{3-4\nu}{16\pi(1-\nu)\mu}
\int_{V} \frac{\boldsymbol{f}}{\rho} \ \mathrm{d} \boldsymbol{y}
+ \dfrac{1}{16\pi(1-\nu)\mu}
\int_{V} (\boldsymbol{f}\cdot\boldsymbol{\rho}) \frac{\boldsymbol{\rho}}{\rho^{3}} \ \mathrm{d} \boldsymbol{y}
$$

特别的，如果选取体力项为在空间点 $\boldsymbol{x}_{0}$ 处的单位集中力 $\boldsymbol{f}_{\delta}^{(m)}$，定义为

$$
(\boldsymbol{f}_{\delta}^{(m)})_{i}(\boldsymbol{x}) \triangleq \delta_{mi} \delta(\boldsymbol{x} - \boldsymbol{x}_{0}),
$$

代入后得到

$$
u_{i}^{(m)}(\boldsymbol{x}; \boldsymbol{x}_{0}) = \frac{3-4\nu}{16\pi(1-\nu)\mu} 
\frac{\delta_{mi}}{\rho}
+ \dfrac{1}{16\pi(1-\nu)\mu}
\frac{\rho_{m}\rho_{i}}{\rho^{3}}, \quad
\boldsymbol{\rho} = \boldsymbol{x} - \boldsymbol{x}_{0}, \quad \rho = \| \boldsymbol{x} - \boldsymbol{x}_{0} \|
$$

方程解的坐标具有对称性：

$$
u_{i}^{(m)}(\boldsymbol{x}; \boldsymbol{y}) = u_{i}^{(m)}(\boldsymbol{y}; \boldsymbol{x})
$$

Betti 互易公式

$$
\int_{\Omega} \mathcal{L}\boldsymbol{u}^{(2)} \cdot \boldsymbol{u}^{(1)}
-\mathcal{L}\boldsymbol{u}^{(1)} \cdot \boldsymbol{u}^{(2)} \ \mathrm{d}\Omega 
= \int_{\partial\Omega} \boldsymbol{t}^{(2)} \cdot \boldsymbol{u}^{(1)}
- \boldsymbol{t}^{(1)} \cdot \boldsymbol{u}^{(2)}\ \mathrm{d}\Gamma
$$

如果考虑基本解，那么边界项的积分可以忽略不计（$\boldsymbol{u}^{(m)}$ 以 $\mathcal{O}(1/\rho)$ 的速率衰减）。那么根据 Betti 互易原理，有

$$
u_{i}^{(j)}(\boldsymbol{x}_{2};\boldsymbol{x}_{1} ) = u_{j}^{(i)}(\boldsymbol{x}_{1};\boldsymbol{x}_{2} )
= u_{j}^{(i)}(\boldsymbol{x}_{2};\boldsymbol{x}_{1} )
$$

Somigliana 公式：

$$
\begin{aligned}
u_{i}(\boldsymbol{x})
&= \int_{V} \boldsymbol{u}^{(i)}( \boldsymbol{x}, \boldsymbol{y} ) \cdot \boldsymbol{f}(\boldsymbol{y}) \ \mathrm{d} V_{\boldsymbol{y}}\\
&+ \int_{\partial V} \boldsymbol{u}^{(i)}( \boldsymbol{x}, \boldsymbol{y} ) \cdot \boldsymbol{t}(\boldsymbol{y})
- \boldsymbol{t}^{(i)}( \boldsymbol{x}, \boldsymbol{y} ) \cdot \boldsymbol{u}(\boldsymbol{y})\ \mathrm{d}S_{\boldsymbol{y}}
\end{aligned}
$$


以下在推导 Eshelby 夹杂问题时，做出如下符号上的约定：

1. 方程中的张量使用**指标形式**。

2. 各向同性弹性力学方程的基本解用符号 $u_{i}^{(m)}(\boldsymbol{x},\boldsymbol{y})$ 表示，根据互易定理，基本解有如下性质：

$$
   u_{i}^{(m)}(\boldsymbol{x},\boldsymbol{y}) = u_{m}^{(i)}(\boldsymbol{y},\boldsymbol{x})
   = u_{m}^{(i)}(\boldsymbol{x},\boldsymbol{y}).
$$

   此外，基本解有如下量级上的估计：

$$
   u_{i}^{(m)} \to \mathcal{O}( \|\boldsymbol{x}-\boldsymbol{y}\|^{-1} ),\ \boldsymbol{x}\to\boldsymbol{y}, \quad
   u_{i}^{(m)} \to \mathcal{O}( \|\boldsymbol{x}\|^{-1} ),\ \boldsymbol{x}\to\infty.
$$

   因此，在靠近奇点处，基本解有一阶奇异性；在远场处，如果三维空间中的面积分项量级在 $\mathcal{o}( \|\boldsymbol{x}\|^{-2} )$，那么认为该积分结果等于零。

3. 函数的自变量分为两类，一类是 $\xi_{i}$，专用于积分项中的空间坐标，在积分运算后，位移场是关于空间坐标 $x_{i}$ 的函数。函数 $f$ 关于坐标 $\xi_{i}$ 的偏导数简记作 $f_{,i}$，关于坐标 $x_{i}$ 的偏导数写作 ${\partial f}/{\partial x_{i}}$。基本解对不同类型的空间坐标偏导数有如下关系：

$$
\frac{\partial f}{\partial x_{i}} = - \frac{\partial f}{\partial \xi_{i}}.
$$

考虑无穷大无体力项弹性区域 $\Omega_{\infty}$，在区域 $\Omega$ 内施加了一个均匀本征应变场 $\boldsymbol{\mu}$，控制方程为

$$
L_{ijkl} u_{k,lj} - L_{ijkl}\mu_{kl,j} = 0
$$

将 $- L_{ijkl}\mu_{kl,j}$ 视作体力项，代入 Somigliana 公式中得到（无穷远边界处的积分被忽略不计）

$$
u_{m}(\boldsymbol{x})  = - \int_{\Omega_{\infty}} u_{i}^{(m)}(\boldsymbol{x},\boldsymbol{\xi}) L_{ijkl}\mu_{kl,j} \ \mathrm{d}\Omega
$$

再应用分部积分，同样略去无穷远边界处的积分，得到

$$
u_{m}(\boldsymbol{x}) = L_{ijkl}\mu_{kl} \int_{\Omega} u_{i,j}^{(m)}(\boldsymbol{x},\boldsymbol{\xi})\ \mathrm{d}\Omega
$$

如果应用 Gauss 公式，上式转化为在夹杂边界处的面积分：

$$
u_{m}(\boldsymbol{x}) 
= \int_{\partial\Omega} 
u_{i}^{(m)}(\boldsymbol{x},\boldsymbol{\xi}) t_{i}^{\mu} \ \mathrm{d}\Gamma, \quad
t_{i}^{\mu} \triangleq L_{ijkl}\mu_{kl} n_{j}.
$$

接下来考虑各向同性对简化，有本构关系

$$
L_{ijkl}\mu_{kl} = \lambda \mu_{kk} \delta_{ij} + 2G \mu_{ij}, \quad
\lambda = \frac{2\nu}{1-2\nu} G
$$

以及基本解

$$
u_{i}^{(m)}(\boldsymbol{x}, \boldsymbol{\xi}) = (3-4\nu)\alpha \frac{\delta_{mi}}{\rho}
+ \alpha\frac{\rho_{m}\rho_{i}}{\rho^{3}}, \quad
\alpha = \dfrac{1}{16\pi(1-\nu)\mu}, \quad
\quad \rho = \| \boldsymbol{x} - \boldsymbol{\xi} \|
$$

基本解的偏导数有

$$
u_{i,j}^{(m)} = -\frac{\partial u_{i}^{(m)}\!\!\!\!\!}{\partial x_{j}}
= (3-4\nu)\alpha \frac{\rho_{j}\delta_{mi}}{\rho^{3}}
- \alpha\frac{\rho_{i}\delta_{mj} + \rho_{m}\delta_{ij}}{\rho^{3}}
+ 3\alpha\frac{\rho_{m}\rho_{i}\rho_{j}}{\rho^{5}}
$$

因此有

$$
u_{m}(\boldsymbol{x}) 
= 2\alpha G \mu_{ij}
\int_{\Omega} (1-2\nu) \Big(
\frac{\rho_{i}\delta_{mj}}{\rho^{3}}
+ \frac{\rho_{j}\delta_{mi}}{\rho^{3}}
- \frac{\rho_{m}\delta_{ij}}{\rho^{3}}
\Big)
+ 3\frac{\rho_{m}\rho_{i}\rho_{j}}{\rho^{5}}
\ \mathrm{d}\Omega
$$

引入三阶张量 $g_{mij}(\boldsymbol{n})$

$$
g_{mij} = (1-2\nu)(n_{i}\delta_{mj} + n_{j}\delta_{mi} - n_{m}\delta_{ij}) + 3n_{m}n_{i}n_{j}, \quad
n_{i} = \rho_{i}/\rho
$$

注意到 $g_{mij} = g_{mji}$，代回得到

$$
u_{m}(\boldsymbol{x}) 
=  \frac{1}{8\pi(1-\nu)} \mu_{ij}
\int_{\Omega}\frac{1}{\rho^{2}} g_{mij}(\boldsymbol{n})
\ \mathrm{d}\Omega
$$

在获取 $u_{m}$ 的解析表达式时，一个核心问题是怎么对形如 $n_{i}/\rho^{2}$ 或 $n_{i}n_{j}n_{k}/\rho^2$ 的被积函数在 $\Omega$ 内进行体积积分。如果区域 $\Omega$ 是椭球形状，并且点 $\boldsymbol{x}$ 在椭球内，那么可以借用椭圆球引力势的求解方法。首先将体积分微元根据当前问题进行简化，表示为沿方向 $\boldsymbol{n}$，长度为 $\rho$，立体角大小为 $\mathrm{d} \omega$ 的圆锥体积：

$$
\mathrm{d}\Omega = \rho^{2} \ \mathrm{d}\rho \mathrm{d}\omega
$$

变量 $\rho$ 的取值范围为 $[0,t(\boldsymbol{n})]$，$t(\boldsymbol{n})$ 为从点 $\boldsymbol{x}$ 出发，沿方向 $\boldsymbol{n}$ 与椭球面 $\partial \Omega$ 相交点之间的距离。在给出上述体积微元之后，原体积积分就可以表示为：

$$
u_{m}(\boldsymbol{x}) 
=  \frac{1}{8\pi(1-\nu)} \mu_{ij}
\int_{S} \int_{0}^{t(\boldsymbol{n})}\frac{1}{\rho^{2}} g_{mij}(\boldsymbol{n})
\rho^{2} \ \mathrm{d}\rho \mathrm{d}\omega
$$

 可以看到，在立体角体积微元表示下，积分项大大简化。首先是 $\rho^{-2}$ 项与体积微元中出现的 $\rho^{2}$ 项抵消，其次是最内层的积分固定方向 $\boldsymbol{n}$ 对 $\rho$ 作积分，因此 $g_{mij}$ 可以直接提到内层积分号之外，最终得到

$$
u_{m}(\boldsymbol{x}) 
=  \frac{1}{8\pi(1-\nu)} \mu_{ij}
\int_{S} t(\boldsymbol{n})g_{mij}(\boldsymbol{n})\ \mathrm{d}\omega
$$

而点 $\boldsymbol{x} - t\boldsymbol{n}$ 在椭球面 $\partial \Omega$ 上，代入到椭球面方程中，得到

$$
\frac{(x_{1} - tn_{1})^{2}}{a_{1}^{2}}
+ \frac{(x_{2} - tn_{2})^{2}}{a_{2}^{2}}
+ \frac{(x_{3} - tn_{3})^{2}}{a_{3}^{2}} = 1
$$

展开之后，得到关于 $t$ 的二次方程：

$$
\underbrace{\frac{x_{1}^{2}}{a_{1}^{2}}
+ \frac{x_{2}^{2}}{a_{2}^{2}}
+ \frac{x_{3}^{2}}{a_{3}^{2}}-1}_{\triangleq e}
- 2\Big(\underbrace{\frac{x_{1}n_{1}}{a_{1}^{2}}
+ \frac{x_{2}n_{2}}{a_{2}^{2}}
+ \frac{x_{3}n_{3}}{a_{3}^{2}}}_{\triangleq f}\Big)t
+ \Big(\underbrace{\frac{n_{1}^{2}}{a_{1}^{2}}
+ \frac{n_{2}^{2}}{a_{2}^{2}}
+ \frac{n_{3}^{2}}{a_{3}^{2}}}_{\triangleq g}
\Big)t^{2}= 0
$$

椭球内过一点的直线与椭球面有两个交点，因此 $t$ 的取值一正一负，$u_{m}$ 可以进一步化简为

$$
u_{m}(\boldsymbol{x}) 
=  \frac{1}{16\pi(1-\nu)} \mu_{ij}
\int_{S} \big(t_{1}(\boldsymbol{n})+t_{2}(\boldsymbol{n})\big)g_{mij}(\boldsymbol{n})\ \mathrm{d}\omega
$$

式中，$t_{1}$ 和 $t_{2}$ 对应同一方向 $\boldsymbol{n}$ 二次方程的两个根，因此有 $t_{1}+t_{2}=2f/g$。重新定义变量 $r$，为坐标系原点沿 $\boldsymbol{n}$ 方向到椭球面的距离，因此 $gr^2=1$，所以积分简化为

$$
u_{m}(\boldsymbol{x}) 
=  \frac{\mu_{ij} x_{k}}{16\pi(1-\nu)} 
\int_{S} 2\lambda_{k} g_{mij}r^{2} \ \mathrm{d}\omega, \quad
\lambda_{k} = n_{k}/a_{k}^{2}
$$

可以看到，椭球内的位移场仅是空间坐标 $x_{k}$ 的一次函数，因此，再次求偏导数之后，应变场为常值函数：

$$
\varepsilon_{mn} = S_{mnij}\mu_{ij}, \quad 
S_{mnij} = \frac{1}{16\pi(1-\nu)} 
\int_{S} (\lambda_{m} g_{nij} + \lambda_{n} g_{mij})r^{2} \ \mathrm{d}\omega
$$

式中，$\mathbb{S}$ 是 Eshelby 张量，是具有次对称性的四阶张量，$S_{mnij}=S_{nmij}=S_{mnji}$。Eshelby 张量的分量是形如 $n_{i}n_{j}\delta_{mn}r^2$ 或 $n_{i}n_{j}n_{m}n_{n}r^2$  作为积分项的积分结果，因此只有积分项指标出现偶数次才能得到非零值。因此，在 Eshelby 张量分量中，所有可能出现的积分为（以下指标均不求和）：

$$
\begin{aligned}
I &= \int_{S} r^{2} \ \mathrm{d}\omega
= 2\pi a_1 a_2 a_3 \int_{0}^{\infty} \frac{\mathrm{d}u}{\Delta}, \\
I_{i} &= \int_{S} \frac{n_{i}^{2} r^{2}}{a_{i}^{2}} \ \mathrm{d}\omega
= 2\pi a_{1} a_{2} a_{3} \int_{0}^{\infty} \frac{1}{(a_{i}^{2}+u)} \frac{\mathrm{d}u}{\Delta}, \\
I_{ii} &= \int_{S} \frac{n_{i}^{4} r^{2}}{a_{i}^{4}} \ \mathrm{d}\omega
= 2\pi a_{1} a_{2} a_{3} \int_{0}^{\infty} \frac{1}{(a_{i}^{2}+u)^2}\frac{\mathrm{d}u}{\Delta}, \\
I_{ij} &= 3\int_{S} \frac{n_{i}^{2}n_{j}^{2} r^{2}}{a_{i}^{2}a_{j}^{2}} \ \mathrm{d}\omega
= 2\pi a_{1} a_{2} a_{3} \int_{0}^{\infty} \frac{1}{(a_{i}^{2}+u)(a_{j}^{2}+u)}\frac{\mathrm{d}u}{\Delta}
\end{aligned}
$$

$$
\Delta^2 \triangleq ( a_{1}^{2} + u )( a_{2}^{2} + u )( a_{3}^{2} + u )
$$

这些积分值有如下代数上的关系：

$$
\sum_{i=1}^{3} a_{i}^{2}I_{i} = I, \quad
\sum_{i=1}^{3} I_{i} = 4\pi, \quad
\sum_{j=1}^{3} a_{j}^{2} I_{ij} = 3I_{i} -2 I_{ii}
$$


**立体角**

考虑空间中的一点 $\boldsymbol{x}$，以及面 $S$，面 $S$ 相对于点 $\boldsymbol{x}$ 的**立体角**定义为这个面投影到以 $\boldsymbol{x}$ 为球心的单位球上的面积：

$$
\Omega = \int_{S} \frac{\hat{\boldsymbol{r}}\cdot \ \mathrm{d}\boldsymbol{A}_{\boldsymbol{\xi}}}{ r^{2}}, \quad
\boldsymbol{r} = \boldsymbol{\xi} - \boldsymbol{x}, \quad r = \|\boldsymbol{\xi} - \boldsymbol{x}\|, \quad \hat{\boldsymbol{r}} = \boldsymbol{r}/r.
$$

因此，如果 $S$ 是一个闭曲面，并且点 $\boldsymbol{x}$ 在闭曲面内，那么曲面 $S$ 相对于点 $\boldsymbol{x}$ 的立体角就等于 $4\pi$。

考虑如下积分

$$
I \triangleq \frac{1}{2\pi} \int_{S} r^{2} \ \mathrm{d} \omega,
$$

其中 $S$ 是椭球面，$\mathrm{d} \omega$ 是立体角微元，$r$ 是椭球面 $S$ 到面内某一点 $\boldsymbol{x}$ 之间距离大小。点 $\boldsymbol{x}$ 是固定且任意的，只要保证在椭球 $S$ **面内**，然而，积分值 $I$ 却与 $\boldsymbol{x}$ 的位置无关。考虑在参数坐标系下描述椭球面方程：

$$
\frac{1}{r^{2}} = \sin^{2}\varphi\Big(
\frac{\cos^{2}\vartheta}{a_{1}^{2}}
+ \frac{\sin^{2}\vartheta}{a_{2}^{2}}\Big)
+ \frac{\cos^{2}\varphi}{a_{3}^{2}}
$$

以及立体角微元在球坐标系下表述为 $\mathrm{d} \omega = \sin\varphi\ \mathrm{d} \varphi\mathrm{d} \vartheta$，再利用该积分在空间中的对称性，可以只选用 $1/8$ 象限内的积分区域，因此有

$$
2\pi I = 8 \int_{0}^{\frac{\pi}{2}}\sin\varphi \ \mathrm{d} \varphi \int_{0}^{\frac{\pi}{2}}
\frac{\mathrm{d} \vartheta}
{\sin^{2}\varphi\Big(
{\cos^{2}\vartheta}/{a_{1}^{2}}
+ {\sin^{2}\vartheta}/{a_{2}^{2}}\Big)
+ \frac{1}{a_{3}^{2}}\cos^{2}\varphi},
$$

被积分项关于 $\vartheta$ 作为积分变量时是偶次项的三角函数有理式，利用常规的换元方法 $t = \tan\vartheta$，可以得到

$$
2\pi I = 8 \int_{0}^{\frac{\pi}{2}}\sin\varphi \ \mathrm{d} \varphi 
\int_{0}^{\infty}
\frac{\mathrm{d} t}
{\sin^{2}\varphi/a_{1}^{2} + \cos^{2}\varphi/a_{3}^{2}
+ t^{2}(\sin^{2}\varphi/a_{2}^{2}+ \cos^{2}\varphi/a_{3}^{2})},
$$

因此，应用积分公式

$$
\int_{0}^{\infty} \frac{\mathrm{d} t}{a_{1}^{2} + a_{2}^{2} t^{2}} = \frac{\pi}{2a_{1}a_{2}},
$$

得到

$$
I = 2a_1 a_2 a_3^2 \int_{0}^{\frac{\pi}{2}}
\frac{\sec^{2}\varphi\sin\varphi \ \mathrm{d} \varphi }
{ (a_{1}^{2} + a_{3}^{2}\tan^{2}\varphi)^{1/2}
(a_{2}^{2} + a_{3}^{2}\tan^{2}\varphi)^{1/2}},
$$

再进行换元

$$
u = a_{3}^{2} \tan^{2}\varphi, \quad
\mathrm{d} u = 2 a_{3}^{2} \sin\varphi\sec^{3}\varphi \ \mathrm{d}\varphi.
$$

就得到

$$
I = a_1 a_2 a_3 \int_{0}^{\infty} \frac{\mathrm{d}u}{\Delta}, \quad
\Delta^2 \triangleq ( a_{1}^{2} + u )( a_{2}^{2} + u )( a_{3}^{2} + u )
$$

进一步的，如果定义 $n_{i}$ 为矢径 $\boldsymbol{r}$ 方向在坐标轴 $i$ 的分量，那么有如下积分公式成立：

$$
\int_{S}r^{2}n_{i}^{2} \ \mathrm{d} \omega 
= 2\pi a_{i}^{2} A_{i}, \quad
A_{i}\triangleq a_1 a_2 a_3 \int_{0}^{\infty} \frac{\mathrm{d}u}{(a_{i}^{2}+u)\Delta}
$$

 也可以直接看到 $A_{i}$ 与积分值 $I$ 之间的关系为

$$
\sum_{i=1}^{3}  a_{i}^{2}A_{i} = I
$$

另外可以直接验证

$$
\frac{\partial I}{\partial a_{i}}
= \frac{I}{a_{i}} - a_{i}A_{i}
$$

and

$$
\sum_{i=1}^{3} A_{i} = 2.
$$

and

$$
\int_{S} r^{4}l_{i}^{2} \ \mathrm{d} \omega = \pi a_{i}^{3} \frac{\partial I}{\partial a_{i}}.
$$

## 更新日志

### 2026/07/24

1. 邱俊淞创建了文档 `Papkovich--Neuber 通解求解椭球夹杂问题.md`

