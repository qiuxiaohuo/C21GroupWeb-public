Hashin-Shtrikman 变分原理中给出的泛函表达式为

$$
U^p (p_{ij}, \widetilde\varepsilon_{ij}) 
= \underbrace{\frac{1}{2}\int \bar\sigma_{ij} \bar\varepsilon_{ij}\mathrm{d} y}_{\triangleq \bar U}
+ \underbrace{\frac{1}{2} \int p_{ij}(y) \widetilde\varepsilon_{ij}(y) \mathrm{d} y}_{\triangleq \widetilde U}
- \frac{1}{2} \int \left( p_{ij}(y)~\delta L_{ijkl}^{-1}(y)~p_{kl}(y) - 2p_{ij}(y) \bar\varepsilon_{ij}\right)  \mathrm{d} y
$$

对二阶张量进行 dev 分解

$$
p_{ij} = p \delta_{ij} + f_{ij},\quad
\bar\varepsilon_{ij} = \bar\varepsilon \delta_{ij} + \bar e_{ij}
$$

代入到泛函表达式当中

$$
U^p = \bar U + \widetilde U - \frac{1}{2} \int
\left(\frac{p^2}{K-K^0} + \frac{f_{ij} f_{ij}}{2(G-G^0)}
-6 p \bar\varepsilon-2 f_{i j} \bar e_{i j}\right) \mathrm{d} y
$$

希望得到使用极化应力表示的 $\widetilde U$ 项表达式, 但是 $\widetilde U$ 项考量的是极化应力与应变扰动量之间的相互作用, 所以包含未知量 $\widetilde \varepsilon_{ij}$​. 幸运的是, 考察的 $\widetilde U$ 是能量项, 而根据 Fourier 分析中的 Parseval 公式, 在相空间中的能量内积等于物理空间中的能量内积: 

$$
\widetilde U = \frac{1}{2} \int p_{ij}(y) \widetilde \varepsilon_{ij}(y) ~\mathrm{d} y
\equiv \frac{1}{2} \int \widehat p_{ij}(\xi) \widehat{\widetilde{u}}_{i,j}^{*}(\xi) ~\mathrm{d} \xi
\tag{0}
$$

并且, 在 Fourier 变换下, 可以找到应变与极化应力的代数表达式. 式 (0) 中, $\widehat{p}_{ij}$, $\widehat{\widetilde{u}}_{i,j}$ 表示场变量 $p_{ij}$, $\widetilde{u}_{i,j}$ 的 Fourier 变换, $\square^{*}$​ 表示共轭复数, 三维的 Fourier 变换定义为

$$
\mathscr{F}[f(y)](\xi) \triangleq \left( \frac{1}{2 \pi} \right)^{3/2} 
\int_{\mathbb{R}^3} f(y) e^{-i \xi \cdot y} \mathrm{d} y
$$

现在希望找到相空间中 $\widehat{p}_{ij}$, $\widehat{\widetilde{u}}_{i,j}$ 之间的关系, 不过首先给出单胞上所满足的微分方程的 statement. 单胞上的控制方程为

$$
\left( L_{ijkl}^{0} \widetilde \varepsilon_{kl} + \widetilde p_{ij} \right)_{,j} = 0
\tag{1}
$$

以及小变形的几何方程

$$
\widetilde \varepsilon_{kl} = \frac{1}{2} (\widetilde{u}_{k,l} + \widetilde{u}_{l,k})
\tag{2}
$$

和齐次位移边界条件:

$$
\widetilde{u}_k(\partial Y) = 0
\tag{3}
$$

方程 (1) 中的 $\widetilde{p}_{ij}$ 表征极化应力的扰动, 定义为

$$
\widetilde{p}_{ij} \triangleq p_{ij} - \bar{p}_{ij}
\tag{4}
$$

接下来, 对控制方程 (1) 几何方程 (2) 作 Fourier 变换得到

$$
i\xi_j \left( L_{ijkl}^{0} \widehat{ \widetilde{\varepsilon} }_{kl}   + \widehat{ \widetilde{p} }_{ij} \right) = 0
$$

和

$$
\widehat{\widetilde{\varepsilon} }_{kl} = \frac{1}{2} i(\xi_l \widehat{\widetilde{u}}_{k} + \xi_k \widehat{\widetilde{u}}_{l})
$$

将各向同性材料的弹性模量

$$
L_{ijkl}^{0} = \lambda^0 \delta_{ij} \delta_{kl} + \mu^{0} \delta_{ik} \delta_{jl} + \mu^{0} \delta_{il} \delta_{jk}
$$

和应变的 Fourier 表达式代入到控制方程的 Fourier 变换式中, 得到关于 $\widehat{\widetilde{u}}_j$ 的代数方程:

$$
(\lambda^0 + \mu^0) \xi_i \xi_k \widehat{\widetilde{u}}_k + \mu^0 \xi^2 \widehat{\widetilde{u}}_i
= i \xi_j \widehat{\widetilde{p}}_{ij}
$$

式中, $\xi^2 = \xi_{kk}$. 上述方程左端的系数矩阵为

$$
A \triangleq \mu^0 \xi^{T} \xi ~\mathbf{I} + (\lambda^0 + \mu^0)\xi \xi^{T}
= \mu^0 \xi^2 \mathbf{I}
+ (\lambda^0 + \mu^0)
\begin{pmatrix}
\xi_1 \xi_1 & \xi_1 \xi_2 & \xi_1 \xi_3 \\
\xi_2 \xi_1 & \xi_2 \xi_2 & \xi_2 \xi_3 \\
\xi_3 \xi_1 & \xi_3 \xi_2 & \xi_3 \xi_3 
\end{pmatrix}
$$

式中, $\mathbf{I}$ 为 3 阶单位矩阵. 对于形式为 $I + \alpha vv^T$ 的矩阵, 其逆矩阵的形式为 $I + \beta vv^T$, 矩阵 $A$ 的逆求解得到

$$
A^{-1} = \frac{1}{\mu^0} \frac{1}{\xi^2} \mathbf{I}
- \frac{\lambda^0 + \mu^0}{\mu^0 (\lambda^0 + 2\mu^0)} \frac{1}{(\xi^2)^2} \xi \xi^{T}
$$

所以, $\widehat{\widetilde{u}}_{i}$ 可以被显式地用 $\widehat{\widetilde{p}}_{ij}$ 表示为

$$
\widehat{\widetilde{u}}_{i} = i \frac{1}{\mu^0} \frac{\xi_k}{\xi^2} \widehat{\widetilde{p}}_{ik}
- i\frac{\lambda^0 + \mu^0}{\mu^0 (\lambda^0 + 2\mu^0)} \frac{\xi_i \xi_k \xi_l}{(\xi^2)^2} \widehat{\widetilde{p}}_{kl}
$$

根据 Fourier 对导函数的变换规则, 只需要对上式两边乘 $i \xi_j$, 就可以得到关于 $\widehat{\widetilde{u}}_{i,j}$ 的 Fourier 变换

$$
\widehat{\widetilde{u}}_{i,j} = i \xi_j \widehat{\widetilde{u}}_{i} 
= - \frac{1}{\mu^0} \frac{\xi_j\xi_k}{\xi^2} \widehat{\widetilde{p}}_{ik}
+ \frac{\lambda^0 + \mu^0}{\mu^0 (\lambda^0 + 2\mu^0)} \frac{\xi_i \xi_j \xi_k \xi_l}{(\xi^2)^2} \widehat{\widetilde{p}}_{kl}
$$


代入到式 (0) 中得到

$$
\widetilde U 
= \frac{1}{2} \int \widehat{\widetilde{p}}_{ij} \widehat{\widetilde{u}}_{i,j}^{*} ~\mathrm{d} \xi
= \frac{1}{2} \int \left( 
- \frac{1}{\mu^0} \frac{\xi_j\xi_k}{\xi^2} \widehat{\widetilde{p}}_{ij} \widehat{\widetilde{p}}_{ik}^{*}
+ \frac{\lambda^0 + \mu^0}{\mu^0 (\lambda^0 + 2\mu^0)} \frac{\xi_i \xi_j \xi_k \xi_l}{(\xi^2)^2} \widehat{\widetilde{p}}_{ij} \widehat{\widetilde{p}}_{kl}^{*}
\right)~\mathrm{d} \xi
\tag{5}
$$

如果假设极化应力是各向同性函数, 也就是说, 极化应力场 $p_{ij}(y)$ 可以通过函数 $P_{ij}(r): \mathbb{R} \mapsto \mathbb{S}^2$ 表示:

$$
\widetilde p_{ij}(y) \equiv \widetilde P_{ij}(|y|)
$$

如果函数在物理空间中是各向同性函数, 那么在相空间中也是各向同性函数

$$
\widehat{\widetilde{p}}_{ij}(\xi) \equiv \widehat{\widetilde{P}}_{ij}(|\xi|)
$$

代入到式 (5) 中, 得到用各向同性函数 $P_{ij}(r)$ 表示的积分

$$
2 \widetilde U 
= \int \left( 
- \frac{1}{\mu^0} \frac{\xi_j\xi_k}{\xi^2} \widehat{\widetilde{P}}_{ij}(|\xi|) \widehat{\widetilde{P}}_{ik}^{*}(|\xi|)
+ \frac{\lambda^0 + \mu^0}{\mu^0 (\lambda^0 + 2\mu^0)} \frac{\xi_i \xi_j \xi_k \xi_l}{(\xi^2)^2} \widehat{\widetilde{P}}_{ij}(|\xi|) \widehat{\widetilde{P}}_{kl}^{*}(|\xi|)
\right)~\mathrm{d} \xi
\tag{6}
$$

如果将式 (6) 转换到球坐标系下进行积分

$$
\xi_1 = r \sin\theta \cos\phi, \quad
\xi_2 = r \sin\theta \sin\phi, \quad
\xi_3 = r \cos\theta , \quad
\mathrm{d} \xi_1 ~\mathrm{d} \xi_2 ~\mathrm{d} \xi_3
= r^2 \sin \theta ~\mathrm{d} r ~\mathrm{d} \theta ~\mathrm{d} \phi
$$

就可以在积分项中消去关于 $\xi_i$ 的有理分式, 得到只关于径向微分 $\mathrm{d} r$ 的积分式

$$
2 \widetilde U 
=  - \frac{1}{\mu^0} \frac{4\pi}{3} \int_{0}^{\infty}
\widehat{\widetilde{P}}_{ij}(r) \widehat{\widetilde{P}}_{ij}^{*}(r) r^2 \mathrm{d} r
+ \frac{\lambda^0 + \mu^0}{\mu^0 (\lambda^0 + 2\mu^0)} \frac{4\pi}{15} 
\int_{0}^{\infty} \left( \widehat{\widetilde{P}}_{mm}(r) \widehat{\widetilde{P}}_{kk}^{*}(r)
+ 2\widehat{\widetilde{P}}_{ij}(r) \widehat{\widetilde{P}}_{ij}^{*}(r) \right) r^2 \mathrm{d} r
\tag{7}
$$

式 (7) 已经和最终的结果很接近了! 只要能够确定式中的径向积分, 这可以通过 Parseval 等式与物理空间中的能量积分关联为

$$
\int_{0}^{\infty} \widehat{\widetilde{P}}_{ij}(r) \widehat{\widetilde{P}}_{ij}^{*}(r) r^2 \mathrm{d} r
= \frac{1}{4\pi} \int_{Y} \widetilde{p}_{ij}(y) \widetilde{p}_{ij}(y) ~\mathrm{d}y, \quad
\int_{0}^{\infty} \widehat{\widetilde{P}}_{mm}(r) \widehat{\widetilde{P}}_{kk}^{*}(r) r^2 \mathrm{d} r
= \frac{9}{4\pi} \int_{Y} \widetilde{p}(y) \widetilde{p}(y) ~\mathrm{d}y
\tag{a}
$$

将式 (4) 代入上式左端项当中, 就得到使用 $p_{ij}$ 表示的积分表达式:

$$
\int_{Y} \widetilde{p}_{ij}(y) \widetilde{p}_{ij}(y) ~\mathrm{d}y
= \int_{Y}  p_{ij}(y) p_{ij}(y)~\mathrm{d}y
- 2 \bar{p}_{ij} \int_{Y} p_{ij}(y)  ~\mathrm{d}y
+ \bar{p}_{ij}\bar{p}_{ij}
= \int_{Y}  p_{ij}(y) p_{ij}(y)~\mathrm{d}y
- \bar{p}_{ij}\bar{p}_{ij}
$$

以及

$$
\int_{Y} \widetilde{p}(y) \widetilde{p}(y) ~\mathrm{d}y
= \int_{Y}  p(y) p(y)~\mathrm{d}y - \bar{p}^2
$$

将分块常值的极化应力表达式

$$
p_{ij} := \sum_{\alpha}^{N} p_{ij}^{\alpha} \chi^{\alpha}
$$

以及极化应力的 dev 分解代入式中得到

$$
\begin{gathered}
\int_{Y} \widetilde{p}_{ij}(y) \widetilde{p}_{ij}(y) ~\mathrm{d}y
= \underbrace{\left( \sum_{\alpha} c^{\alpha} f_{ij}^{\alpha} f_{ij}^{\alpha} - \bar{f}_{ij}\bar{f}_{ij} \right)}_{\triangleq \clubsuit}
+ 3\underbrace{\left(\sum_{\alpha} c^{\alpha} (p^{\alpha})^2 - \bar{p}^2 \right)}_{\triangleq \spadesuit}
= \clubsuit + 3\spadesuit \\
\int_{Y} \widetilde{p}(y) \widetilde{p}(y) ~\mathrm{d}y = \left(\sum_{\alpha} c^{\alpha} (p^{\alpha})^2 - \bar{p}^2 \right) = \spadesuit
\end{gathered}\tag{b}
$$

最终, 代入到式 (7) 中, 得到

$$
\begin{aligned}
2 \widetilde U 
&\stackrel{(a)}{=} - \frac{1}{3\mu^0}
\int_{Y} \widetilde{p}_{ij}(y) \widetilde{p}_{ij}(y) ~\mathrm{d}y
+ \frac{\lambda^0 + \mu^0}{\mu^0 (\lambda^0 + 2\mu^0)} \frac{2}{15} 
\int_{Y} \widetilde{p}_{ij}(y) \widetilde{p}_{ij}(y) ~\mathrm{d}y
+ \frac{\lambda^0 + \mu^0}{\mu^0 (\lambda^0 + 2\mu^0)} \frac{3}{5}
\int_{Y} \widetilde{p}(y) \widetilde{p}(y) ~\mathrm{d}y \\
&\stackrel{(b)}{=} - \frac{1}{3\mu^0} ( \clubsuit + 3\spadesuit )
+ \frac{\lambda^0 + \mu^0}{\mu^0 (\lambda^0 + 2\mu^0)} \frac{2}{15} ( \clubsuit + 3\spadesuit )
+ \frac{\lambda^0 + \mu^0}{\mu^0 (\lambda^0 + 2\mu^0)} \frac{3}{5} \spadesuit \\
&= \left(-\frac{1}{3\mu^0} + \frac{\lambda^0 + \mu^0}{\mu^0 (\lambda^0 + 2\mu^0)} \frac{2}{15} \right)\clubsuit
+ \left(-\frac{1}{\mu^0} + \frac{\lambda^0 + \mu^0}{\mu^0 (\lambda^0 + 2\mu^0)}\right)\spadesuit \\
&= \frac{-3\lambda^0 -8\mu^0}{15\mu^0 (\lambda^0 + 2\mu^0)} \clubsuit
+ \frac{-1}{\lambda^0 + 2\mu^0}\spadesuit
\end{aligned}
$$

如果将

$$
\lambda^0 = K^0 - \frac{2}{3}\mu^0
$$

代入, 就可以得到

$$
2 \widetilde U 
= \frac{-3 K^0 - 6\mu^0}{5\mu^0 (3 K^0 + 4\mu^0)} \left( \sum_{\alpha} c^{\alpha} f_{ij}^{\alpha} f_{ij}^{\alpha} - \bar{f}_{ij}\bar{f}_{ij} \right)
+ \frac{-3}{3 K^0 + 4\mu^0} \left(\sum_{\alpha} c^{\alpha} (p^{\alpha})^2 - \bar{p}^2 \right)
\tag{$\bigstar$}
$$

## 更新日志

### 2026/07/24

1. 创建了文档 `Papkovich--Neuber 通解求解椭球夹杂问题.md`

