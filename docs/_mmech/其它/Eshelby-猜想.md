# Eshelby 猜想

这篇文档推导 Liu (2007)[^1] 给出的 Eshelby conjectures 的解。对于空间 $\mathbb{R}^{d}$ 中的偏微分方程：

$$
\nabla\cdot( \mathbb{L} : \nabla\boldsymbol{v} + \boldsymbol{P}\chi_{\Omega} ) = \boldsymbol{0},
$$

区域 $\Omega$ 是空间 $\mathbb{R}^{d}$ 的夹杂，$\chi_{\Omega}$ 是关于夹杂区域 $\Omega$ 的特征函数，$\boldsymbol{P}$ 是极化应力。方程的解可通过 Fourier 变换得到得到

$$
(\nabla \boldsymbol{v})_{pi}
= \frac{-1}{(2\pi)^{n}} \int_{\mathbb{R}^{d}} N_{pq}(\boldsymbol{\xi}) P_{qj}\xi_{j}\xi_{i}
\int_{\Omega} e^{i\boldsymbol{\xi}\cdot(\boldsymbol{x}-\boldsymbol{x}')} \,\mathrm{d}\boldsymbol{x}'\,\mathrm{d}\boldsymbol{\xi},
$$

式中，$N_{pq}$ 是 $L_{piqj}\xi_{i}\xi_{j}$ 的逆。（见 Mura, 1987）如果弹性刚度张量 $\mathbb{L}$ 表示为如下各向同性张量的组合，

$$
L_{piqj} = \lambda \delta_{pi}\delta_{qj} 
+ \mu_{1} \delta_{ij}\delta_{pq} + \mu_{2} \delta_{pj}\delta_{qi},
$$

那么 $N_{pq}$ 等于

$$
N_{pq} = \frac{1}{\mu_{1}|\xi|^2}\delta_{pq}
- \frac{\mu_{2}+\lambda}{\mu_{1}(\lambda+\mu_{1}+\mu_{2})} \frac{\xi_{p}\xi_{q}}{|\xi|^{4}},
$$

式中，$|\xi| = \sqrt{\xi_{i}\xi_{i}}$。而对于在 $\mathbb{R}^{d}$ 空间中的 Poisson 方程，

$$
\Delta u = -\chi_{\Omega},
$$

方程解的二阶导数可通过 Newtonian 位势和 Fourier 变换表示为

$$
(\nabla \nabla u)_{ij}
= \frac{-1}{(2\pi)^{n}} \int_{\mathbb{R}^{d}} \frac{\xi_{i}\xi_{j}}{|\xi|^{2}}
\int_{\Omega} e^{i\boldsymbol{\xi}\cdot(\boldsymbol{x}-\boldsymbol{x}')} \,\mathrm{d}\boldsymbol{x}'\,\mathrm{d}\boldsymbol{\xi},
$$

如果 $P_{qj} = \delta_{qj}$，那么有

$$
N_{pq}(\boldsymbol{\xi}) \xi_{q}\xi_{i}
= \frac{1}{\lambda+\mu_{1}+\mu_{2}}  \frac{\xi_{p}\xi_{i}}{|\xi|^{2}},
$$

此时原问题和 Newtonian 位势问题的联系为

$$
\nabla \boldsymbol{v} = \frac{1}{\lambda+\mu_{1}+\mu_{2}} \nabla\nabla u.
$$

而如果材料参数满足 $\lambda+\mu_{1}=0$，此时 $N_{pq}$ 中只保留 $\delta_{pq}$ 项，同样可以和 Newtonian 位势关联：

$$
\nabla \boldsymbol{v} = \boldsymbol{P} \nabla\nabla u / \mu_{1}
$$

构造变分原理

$$
G_{r}(u_{r}) = \inf_{v\in K_{r}} \int_{B_{r}} \frac{1}{2}|\nabla v|^2 \,\mathrm{d}\boldsymbol{x}
$$

## 更新日志

### 2026/07/24

1. 创建了文档 `Eshelby 猜想.md`

[^1]: Liu, L.P., 2008. Solutions to the Eshelby conjectures. Proc. R. Soc. A. 464, 573–594. https://doi.org/10.1098/rspa.2007.0219
