# Gurson 塑性模型

这篇文档旨在通过 McClintock 解和 Gurson 模型的推导，理解如何将细观力学方法应用在损伤问题中。

## McClintock 解

McClintock 最先给出理想塑性设置下空心圆柱或圆球体积扩张的解析解[^1]。这里主要沿用李少凡等人的推导过程[^2]。考虑如下具有轴向对称性的圆柱模型，其中空心圆柱的半径为 $a$，其它区域的材料设置为理想 J2 塑性，应力的偏张量部分 $\boldsymbol{s}$ 和应变率 $\dot{\boldsymbol{\varepsilon}}$ 之间的关系为

$$
\begin{equation}\label{eq:cstt}
\dot{\boldsymbol{\varepsilon}} = \frac{3}{2} \frac{ \dot{\varepsilon}_{\mathrm{eq}} }{ \sigma_{y} } \boldsymbol{s},
\end{equation}
$$

式中，$\dot{\varepsilon}_{\mathrm{eq}} \triangleq \sqrt{\frac{2}{3} \dot{\varepsilon}_{ij} \dot{\varepsilon}_{ij}}$ 是等效塑性应变率，$\sigma_{y}$ 是屈服应力。上式给出了应力应变之间的非线性本构关系。

<img src="https://raw.githubusercontent.com/qiuxiaohuo/picgo-bed/main/figs/20260702103705550_McClintock.png" alt="image-20260702103705504" style="zoom:50%;" />

考虑问题具有反对称性，并且是在面 1-2 的广义平面应变问题，因此速度场 $\dot{\boldsymbol{u}}$ 在柱坐标系下的分量可以设为如下函数形式：

$$
(\dot{\boldsymbol{u}})_{r} = v(r), \quad
(\dot{\boldsymbol{u}})_{z} = e_{z} z, \quad 
(\dot{\boldsymbol{u}})_{\theta} = 0,
$$

式中，$e_{z}$ 是常值。在给定速度场的函数形式之后，对应的应变率分量可以表示为

$$
\begin{equation}\label{eq:rate_e}
\dot{\varepsilon}_{rr} = \frac{\mathrm{d} v}{\mathrm{d} r}, \quad
\dot{\varepsilon}_{\theta\theta} = \frac{1}{r} v, \quad
\dot{\varepsilon}_{zz} = e_{z}, \quad
\dot{\varepsilon}_{r\theta} = \dot{\varepsilon}_{rz} = \dot{\varepsilon}_{\theta z} = 0.
\end{equation}
$$

由于材料设置为理想塑性，因此应变率的体积部分应恒等于零，也即

$$
\dot{\varepsilon}_{rr} + \dot{\varepsilon}_{\theta\theta} + \dot{\varepsilon}_{zz} \equiv 0,
$$

这就得到关于径向函数 $v$ 的微分方程：

$$
\frac{\mathrm{d} (rv)}{\mathrm{d} r} + re_{z} = 0,
$$

对上式在区间 $(a,\rho)$ 内作积分得到

$$
\rho v(\rho) - av_{a} = \int_{a}^{\rho} \,\mathrm{d}(rv) = -\int_{a}^{\rho} e_{z}r\,\mathrm{d}r
= -\frac{e_{z}}{2}( \rho^{2} - a^{2} ),
$$

式中，$v_{a}$ 是空心圆柱边界的扩展速度。速度场的表达式由此给出：

$$
v(r)
= - \frac{e_{z}}{2 }r + \Big( \frac{e_{z}}{2}a^{2} + av_{a} \Big) \frac{1}{r} ,
$$

将速度场表达式代入到应变率分量 $\eqref{eq:rate_e}$ 得到

$$
\begin{equation}\label{eq:rate_e_explicit}
\begin{aligned}
\dot{\varepsilon}_{rr} 
&= - \frac{e_{z}}{2 } - \Big( \frac{e_{z}}{2}a^{2} + av_{a} \Big) \frac{1}{r^{2}}, \\
\dot{\varepsilon}_{\theta\theta} 
&= - \frac{e_{z}}{2 } + \Big( \frac{e_{z}}{2}a^{2} + av_{a} \Big) \frac{1}{r^{2}}, \\
\dot{\varepsilon}_{zz} 
&= e_{z},
\end{aligned}
\end{equation}
$$

以及等效应变率

$$
\dot{\varepsilon}_{\mathrm{eq}}
= \sqrt{\frac{2}{3} \big( 
\dot{\varepsilon}_{rr}^{2} + \dot{\varepsilon}_{\theta\theta}^{2} + \dot{\varepsilon}_{zz}^{2}
\big)}
= \mathrm{sign}\,(e_{z}) e_{z}\sqrt{1 + \frac{1}{3r^{4}} \big( a^2 + 2a\lambda \big)^{2}},\quad \lambda \triangleq v_{a}/ e_{z}.
$$

将应变率和等效应变率代入到本构关系 $\eqref{eq:cstt}$ 中得到

$$
\begin{equation}\label{eq:cstt_expl}
\frac{1}{\sigma_{y}}\begin{pmatrix}
\sigma_{rr} - \sigma \\
\sigma_{\theta\theta} - \sigma \\
\sigma_{zz} - \sigma
\end{pmatrix}
= \frac{\mathrm{sign}\,(e_{z})}{\sqrt{9 + 3 \big[ ( a^{2} + 2a\lambda )/r^{2} \big]^{2} }}
\begin{pmatrix}
- 1 - \Big( a^{2} + 2a\lambda \Big)/r^{2} \\
- 1 + \Big( a^{2} + 2a\lambda \Big)/r^{2} \\
2
\end{pmatrix}.
\end{equation}
$$

上述给出的三个方程并不是完全独立的，因此还需要额外的平衡方程进行补充，

$$
\begin{equation}\label{eq:equili}
\frac{\mathrm{d} \sigma_{rr}}{\mathrm{d} r}
+ \frac{\sigma_{rr} - \sigma_{\theta\theta}}{r} = 0,
\end{equation}
$$

其中，径向应力减环向应力分量可通过式 $\eqref{eq:cstt_expl}$ 给出

$$
\sigma_{rr} - \sigma_{\theta\theta}
= -\frac{2\mathrm{sign}\,(e_{z})\sigma_{y}( a^{2} + 2a\lambda )/r^{2}}
{\sqrt{9 + 3 \big[ ( a^{2} + 2a\lambda )/r^{2} \big]^{2} }}.
$$

将上式代入到平衡方程 $\eqref{eq:equili}$ 中得到

$$
\frac{\sigma_{rr}^{\infty}}{\sigma_{y}}
= \mathrm{sign}\,(e_{z})\int_{a}^{\infty} \frac{( a^{2} + 2a\lambda )/r^{2}}
{\sqrt{9 + 3 \big[ ( a^{2} + 2a\lambda )/r^{2} \big]^{2} }}\frac{2}{r}\,\mathrm{d} r,
$$

式中，$\sigma_{rr}^{\infty}$ 是径向应力分量在 $r\to\infty$ 的取值。记 $x= ( a^{2} + 2a\lambda )/r^{2}$，那么 $2\,\mathrm{d}r/r=-\,\mathrm{d}x/x$，以及 $r=a$ 对应 $x=1+2\lambda/a\triangleq x_{a}$，$r\to\infty$ 对应 $x\to0$，所以上述积分转化为

$$
\begin{equation}\label{eq:stress_radial}
\frac{\sigma_{rr}^{\infty}}{\sigma_{y}}
= \mathrm{sign}\,(e_{z})\int_{0}^{x_{a}} \frac{1}{\sqrt{9 + 3 x^{2} }}\,\mathrm{d} x
= \mathrm{sign}\,(e_{z})\frac{1}{\sqrt{3}} \mathrm{arsinh}\,\big( x_{a}/\sqrt{3} \big),
\end{equation}
$$

因而扩展速度 $v_{a}$ 可表示为关于横向远场应力 $\sigma_{rr}^{\infty}$ 的函数：

$$
\begin{equation}\label{eq:mcclintock}
v_{a}/a = \frac{\sqrt{3}|e_{z}|}{2}\sinh\frac{\sqrt{3}\sigma_{rr}^{\infty}}{\sigma_{y}} - \frac{e_{z}}{2}
\end{equation}
$$

式 $\eqref{eq:mcclintock}$ 给出 McClintock 解。

## Gurson 模型

McClintock 解 $\eqref{eq:mcclintock}$ 给出了无穷大区域内应变和偏应力的场分布，接下来将通过对细观局域场在**有限区域** $V$ 内平均化的方法，得到在宏观应力空间中的屈服面方程[^3]。设宏观应力为 $\boldsymbol{\Sigma}$，不变量 $\Sigma_{m}=\Sigma_{ii}$ 和 $\Sigma_{\mathrm{eq}}=\sqrt{3J_{2}(\boldsymbol{\Sigma})}$，空相的体积分数为 $c$。McClintock 求解模型问题时取边界 $r$ 无穷大，这里考虑空相在**有限体积的弹性介质** $r=b>a$ 内，因此在上一节中的积分将限制在有限的区域内。设宏观应力与远场应力之间的关系为

$$
\begin{equation}\label{eq:macro_stress}
\Sigma_{11} = \Sigma_{22} = \sigma_{rr}|_{\partial V}, \quad
\Sigma_{33} = \sigma_{zz}|_{\partial V}.
\end{equation}
$$

首先考察宏观应力分量 $\Sigma_{11}$，将积分式 $\eqref{eq:stress_radial}$ 修改在有限的区间 $r=[a,b]$ 内，得到

$$
\begin{equation}\label{eq:sigma11_macro}
\Sigma_{11}/\sigma_{y}
= \frac{\mathrm{sign}\,(e_{z})}{\sqrt{3}} \Big( 
\mathrm{arsinh}\,\big( x_{a}/\sqrt{3} \big)
- \mathrm{arsinh}\, \big( x_{b}/\sqrt{3} \big)\Big),
\end{equation}
$$

式中，$x_{b}=c x_{a}$。绘制归一化后的宏观应力分量 $\Sigma_{11}/\sigma_{y}$ 关于半径的图像如图所示。

<img src="https://raw.githubusercontent.com/qiuxiaohuo/picgo-bed/main/figs/20260703_gurson_sigma11.png" alt="gurson" style="zoom: 67%;" />

接下来考察轴向的宏观应力分量 $\Sigma_{33}$。对式 $\eqref{eq:cstt_expl}$ 的第三个方程在整个区域 $V$ 内求体积平均，并利用宏观应力等于细观应力场的体积平均，以及对宏观应力分量的假设 $\eqref{eq:macro_stress}$，得到

$$
\frac{1}{V}\int_{V} \frac{\sigma_{zz} - \sigma}{\sigma_{y}} \,\mathrm{d}V
= \frac{2}{3\sigma_{y}} ( \Sigma_{33} - \Sigma_{11} )
= \frac{2\times2\pi h}{\pi b^{2} h}\int_{a}^{b}
\frac{\mathrm{sign}\,(e_{z})}{\sqrt{9 + 3 \big[ ( a^{2} + 2a\lambda )/r^{2} \big]^{2} }}
r\,\mathrm{d}r,
$$

整理之后，并应用换元公式，得到

$$
\big( \Sigma_{33} - \Sigma_{11} \big) / \sigma_{y}
= \frac{x_{b}}{ \mathrm{sign}\,(e_{z}) }   
\int_{x_{b}}^{x_{a}} \frac{3}{x^{2}\sqrt{9 + 3 x^{2} }} \,\mathrm{d} x\\
= \mathrm{sign}\,(e_{z})
\Big( \sqrt{1+x_{b}^{2}/3} - c\sqrt{1+x_{a}^{2}/3} \Big).
$$

在式 $\eqref{eq:macro_stress}$ 约定的宏观应力状态下，上式的绝对值等于宏观等效应力，

$$
\big( \Sigma_{\mathrm{eq}}  / \sigma_{y} \big)^2
= \Big( \sqrt{1+x_{b}^{2}/3} - c\sqrt{1+x_{a}^{2}/3} \Big)^{2} \\
= 1 + c^{2} + \frac{2}{3}fx_{b}x_{b} - 2f\sqrt{1+x_{b}^{2}/3}\sqrt{1+x_{a}^{2}/3}
$$

注意到式 $\eqref{eq:sigma11_macro}$ 可以转化为

$$
\Sigma_{11}/\sigma_{y}
=  \cosh (\sqrt{3}\Sigma_{11}/\sigma_{y})
= \sqrt{1 + x_{a}^{2}/3 }\sqrt{1 + x_{b}^{2}/3 } - \frac{1}{3}x_{a}x_{b},
$$

代入宏观等效应力的表达式当中，就得到 Gurson 屈服面模型：

$$
\begin{equation}\label{eq:gurson}
\big( \Sigma_{\mathrm{eq}}  / \sigma_{y} \big)^2
+ 2f \cosh (\sqrt{3}\Sigma_{11}/\sigma_{y}) - (1 + c^{2}) = 0.
\end{equation}
$$

Gurson 屈服面绘制在 $\Sigma_{m}$-$\Sigma_{\mathrm{eq}}$ 坐标平面内的图像如图所示。

<img src="https://raw.githubusercontent.com/qiuxiaohuo/picgo-bed/main/figs/20260703gurson_yld.png" alt="gurson_yld" style="zoom:67%;" />

## Gurson-Tvergaard-Needleman 模型

Gurson-Tvergaard-Needleman (GTN) 模型[^4]在式 $\eqref{eq:gurson}$ 中引入额外的参数 $q_{1}$ 和修正的体积分数 $c^{*}$，以更好地和实验数据贴合，

$$
\begin{equation}\label{eq:gtn}
\Phi=\big( \Sigma_{\mathrm{eq}}  / \sigma_{y} \big)^2
+ 2q_{1}c^{*} \cosh (\frac{3\Sigma_{m}}{2\sigma_{y}}) - \big( 1 + (q_{1}c^{*})^{2} \big) \equiv 0,
\end{equation}
$$

式中，引入的无量纲参数 $q_{1}$ 可以更好地贴合周期性空洞的结果；$c^{*}$ 是关于体积分数 $c$ 的函数，在 $c$ 小于阈值 $c_{C}$ 时与 Gurson 模型保持一致，而在大于阈值 $c_{C}$ 时使用断裂体积分数 $c_{F}$ 修正：

$$
c^{*} = \begin{cases}
c, & c \leq c_{C}, \\
c_{C} + \dfrac{1/q_{1} - c_{C}}{c_{F} - c_{C}} (c - c_{C}) , & c>c_{C}.
\end{cases}
$$

当 $c=c_{F}$，此时 $c^{*}=1/q_{1}$，代入屈服面函数 $\eqref{eq:gtn}$ 后表示完全损伤状态。宏观流动法则为

$$
\dot{\boldsymbol{E}}^{p} = \dot{\lambda} \frac{\partial \Phi}{\partial \boldsymbol{\Sigma}}
= \dot{\lambda} \Big( \frac{1}{3} \frac{\partial \phi}{\partial \Sigma_{m}} \boldsymbol{I}
+ \frac{\partial \phi}{\partial \Sigma_{\mathrm{eq}}} \boldsymbol{N} \Big)
$$

式中，$\dot{\boldsymbol{E}}^{p}$ 时宏观塑性应变率，$\dot{\lambda}$ 是塑性乘子，二阶张量 $\boldsymbol{I}$ 和 $\boldsymbol{N}$ 分别是单位张量和单位偏张量。可以看到，宏观塑性流动法则不再是 $J_{2}$ 类型。

## 更新日志

### 2026/07/24

1. 邱俊淞创建了文档 `Gurson 塑性模型.md`

[^1]: McClintock, F.A., 1968. A Criterion for Ductile Fracture by the Growth of Holes. Journal of Applied Mechanics 35, 363–371. https://doi.org/10.1115/1.3601204
[^2]: Li, S., Wang, G., 2018. Introduction to Micromechanics and Nanomechanics, 2nd ed. WORLD SCIENTIFIC. https://doi.org/10.1142/8644
[^3]: Gurson, A.L., 1977. Continuum Theory of Ductile Rupture by Void Nucleation and Growth: Part I—Yield Criteria and Flow Rules for Porous Ductile Media. Journal of Engineering Materials and Technology 99, 2–15. https://doi.org/10.1115/1.3443401
[^4]: Tvergaard, V., Needleman, A., 1984. Analysis of the cup-cone fracture in a round tensile bar. Acta Metallurgica 32, 157–169. https://doi.org/10.1016/0001-6160(84)90213-X
