# Fourier 方法求解椭球夹杂问题

对于各向同性材料, 方程左端的系数矩阵为

$$
\hat K \triangleq \mu \xi^{T} \xi ~\mathbf{I} + (\lambda + \mu)\xi \xi^{T}
= \mu |\xi|^2 \mathbf{I}
+ (\lambda + \mu)
\begin{pmatrix}
\xi_1 \xi_1 & \xi_1 \xi_2 & \xi_1 \xi_3 \\
\xi_2 \xi_1 & \xi_2 \xi_2 & \xi_2 \xi_3 \\
\xi_3 \xi_1 & \xi_3 \xi_2 & \xi_3 \xi_3 
\end{pmatrix}
$$

式中, $\mathbf{I}$ 为 3 阶单位矩阵, $\xi = (\xi_1, \xi_2, \xi_3)^{\top}$, $|\xi|^2 = \xi_1^2 + \xi_2^2 + \xi_3^2$. 对于形式为 $I + \alpha vv^T$ 的矩阵, 其逆矩阵的形式依旧为 $I + \beta vv^T$, 因此矩阵 $\hat{K}$ 的逆为

$$
\hat{K}^{-1}  = \frac{1}{\mu} \frac{1}{|\xi|^2} \mathbf{I}
- \frac{\lambda + \mu}{\mu (\lambda + 2\mu)} \frac{1}{|\xi|^4} \xi \xi^{\top}, \quad
\det(\hat{K}) = (\lambda + 2\mu)\mu^2 |\xi|^6
$$

各项同性材料的 Green 函数为

$$
G(x) = \frac{1}{(2\pi)^3} \frac{1}{\mu} \int_{\mathbb{R}^3} \frac{1}{|\xi|^2} \mathbf{I} ~\mathrm{d} \xi
- \frac{1}{(2\pi)^3} \frac{\lambda + \mu}{\mu (\lambda + 2\mu)} 
\int_{\mathbb{R}^3} \frac{1}{|\xi|^4} \xi \xi^{\top} \mathrm{d} \xi
$$

将积分分成两部分进行计算

$$
I_0(x) = \int_{\mathbb{R}^3} \frac{1}{|\xi|^2} e^{i\xi \cdot x} ~\mathrm{d} \xi, \quad
I_{ij}(x) = - \frac{\partial^2}{\partial x_i \partial x_j} 
\int_{\mathbb{R}^3} \frac{1}{|\xi|^4} e^{i\xi \cdot x} \mathrm{d} \xi
$$

一波操作得到

$$
I_0(x) = 2\pi^2 \frac{1}{|x|}, \quad
I_{ij} = \pi^2 \frac{\partial^2 |x|}{\partial x_i \partial x_j}
$$

代入方程中得到

$$
G_{ij}(x) = \frac{1}{4\pi \mu} \frac{1}{|x|} \delta_{ij}
- \frac{\lambda + \mu}{8\pi \mu (\lambda + 2\mu)} 
\frac{\partial^2}{\partial x_i \partial x_j}  |x|
$$

对于径向对称幂函数 $1/|\xi|^a$, Fourier 逆变换依旧得到径向对称幂函数 $C_a/|x|^{d-a}$, 这时需要确定常数 $C_a$, 当 $a < d$ 时, 有

$$
C_a=(2 \pi)^{\frac{n}{2}} \frac{2^{\frac{n-a}{2}} \Gamma\left(\frac{n-a}{2}\right)}{2^{\frac{a}{2}} \Gamma\left(\frac{a}{2}\right)} .
$$

[Radial functions and the Fourier transform](https://math.arizona.edu/~faris/methodsweb/hankel.pdf)
对于 $a > d$ 的情况, 积分是不收敛的, 但注意到 $I_{ij}$ 还是关于 $x_i$ 的二阶导数, 这就抵消了 2 幂次, 积分是收敛的. 尝试计算 $I_{11}(x)$, 并令 $x=(0,0,1)$, 最终得到

$$
C_{4} = 4\pi \int_{0}^{\infty} \frac{1}{r^3} \sin r - \frac{1}{r^2} \cos r ~\mathrm{d}r
$$

有意思的是, 如果直接通过计算积分 $\int_{\mathbb{R}^3} \frac{1}{|\xi|^4} e^{i\xi \cdot x} \mathrm{d} \xi$ 求解常数 $C_4$, 最终会得到上式积分项中关于 $\sin r$ 部分

$$
C_4 = 4\pi \int_{0}^{\infty} \frac{1}{r^3} \sin r ~\mathrm{d}r
$$

而这个积分是不收敛的。

## 更新日志

### 2026/07/24

1. 创建了文档 `Fourier 方法求解椭球夹杂问题.md`
