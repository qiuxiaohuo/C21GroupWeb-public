# 非均匀场变换分析

这篇文档希望理解 NTFA 是如何拓展 TFA 的工作，如何通过广义标准材料获得假设本征应变模式的系数的演化方程。

## 单胞问题

### 自由能和耗散势

设复合材料 $\Theta$ 由 $N$ 相组分构成，组分材料可以使用 Helmholtz 自由能 $\omega$ 和耗散势函数 $\varphi$ 描述本构关系（包括状态和演化方程），在单胞区域内，自由能和耗散势是关于细观尺度坐标 $\boldsymbol{y}$ 的标量值场变量

$$
\begin{equation}\label{eq:micro_potential}
\omega(\boldsymbol{y},\, \boldsymbol{\varepsilon},\, \boldsymbol{\mu})
=\sum_{\alpha=1}^{N} \chi^{(\alpha)}(\boldsymbol{y}) \omega^{(\alpha)}( \boldsymbol{\varepsilon},\, \boldsymbol{\mu}),\quad
\varphi(\boldsymbol{y},\, \dot{\boldsymbol{\mu}})
=\sum_{\alpha=1}^{N} \chi^{(\alpha)}(\boldsymbol{y}) \varphi^{(\alpha)}( \dot{\boldsymbol{\mu}} ),
\end{equation}
$$

式中，$\omega^{(\alpha)}$ 和 $\varphi^{(\alpha)}$ 分别是材料相 $\alpha$ 单位体积的 Helmholtz 自由能和耗散势函数，$\boldsymbol{\varepsilon}$ 和 $\boldsymbol{\mu}$ 分别是材料在点 $\boldsymbol{y}$ 处的应变和状态变量，记号 $\dot{\square}$ 表示物理量 $\square$ 的物质导数，$\chi^{(\alpha)}$ 是材料相 $\alpha$ 的特征函数。以上预设的函数形式可参考以凸分析为基础的广义标准材料模型。单胞域内的应力平衡方程和几何方程分别为

$$
\begin{equation}\label{eq:uc_eqs}
\left\{\begin{aligned}
&\nabla \cdot\boldsymbol{\sigma} 
= \boldsymbol{0},\\
&\boldsymbol{\varepsilon} 
= \boldsymbol{\varepsilon}^{\mathtt{c}} + \tilde{\boldsymbol{\varepsilon}},\quad
\tilde{\boldsymbol{\varepsilon}} \in \mathscr{E}_{0},
\end{aligned}\right.
\end{equation}
$$

式中，$\boldsymbol{\varepsilon}^{\mathtt{c}}$ 是施加的宏观应变，$\tilde{\boldsymbol{\varepsilon}}$ 是运动学允许的应变扰动场，$\langle \tilde{\boldsymbol{\varepsilon}} \rangle=\boldsymbol{0}$，其中 $\langle \bullet \rangle$ 是单胞区域内的体积平均算子。同样还需要补充的是局部的状态方程和演化方程，由自由能和耗散势函数的偏导数分别给出：

$$
\begin{equation}\label{eq:stat_dissip_strong_local}
\begin{cases}
\boldsymbol{\sigma}
= \dfrac{\partial \omega}{\partial \boldsymbol{\varepsilon}}, & \text{state law},\\
\dfrac{\partial \omega}{\partial \boldsymbol{\mu}}
+ \dfrac{\partial \varphi}{\partial \dot{\boldsymbol{\mu}}} = 0, &\text{evolution law}.
\end{cases}
\end{equation}
$$

宏观的应力响应为细观尺度应力场的体积平均：

$$
\begin{equation}\label{eq:def_macro_stress}
\boldsymbol{\sigma}^{\mathtt{c}} \triangleq \langle \boldsymbol{\sigma} \rangle.
\end{equation}
$$

### 增量变分原理

现在考虑将分析步时间离散为为 $t_{0}=0,\, t_{1},\, t_{2},\, \ldots,t_{n},\ldots$，并已知在 $t_{n}$ 以及之前所有时刻的场变量取值，那么定义 $t_{n}$ 时刻的增量泛函

$$
\begin{equation}\label{eq:micro_inc_funcl}
\mathcal{F}_{\Delta}(\boldsymbol{y},\, \boldsymbol{\varepsilon},\, \boldsymbol{\mu})
\triangleq \sum_{\alpha=1}^{N} \chi^{(\alpha)}(\boldsymbol{y})
\Big( \omega^{(\alpha)}(\boldsymbol{\varepsilon},\, \boldsymbol{\mu})
+ \Delta t \varphi^{(\alpha)}\big( \tfrac{\boldsymbol{\mu}-\boldsymbol{\mu}_{n}}{\Delta t} \big)  \Big),
\end{equation}
$$

那么在点 $\boldsymbol{y}$ 处，和式 $\eqref{eq:stat_dissip_strong_local}$ 等价的增量变分原理为：

$$
\begin{equation}\label{eq:stat_dissip_vari_local}
\omega_{\Delta} (\boldsymbol{y},\, \boldsymbol{\varepsilon})
\triangleq \inf_{\boldsymbol{\mu}\in \mathscr{U}} \mathcal{F}_{\Delta}(\boldsymbol{y},\, \boldsymbol{\varepsilon},\, \boldsymbol{\mu}),\quad
\boldsymbol{\sigma} = \frac{\partial \omega_{\Delta}}{\partial \boldsymbol{\varepsilon}},
\end{equation}
$$

式中，有限维空间 $\mathscr{U}$ 表示 $\boldsymbol{\mu}$ 所属的向量空间，其维度记作 $K$。泛函取极值的欧拉--拉格朗日方程则对应为演化方程 $\eqref{eq:stat_dissip_strong_local}$。Suquet 等人通过单胞域内的变分问题，定义宏观的应变能

$$
\begin{equation}\label{eq:macro_helmholtz}
W_{\Delta}(\boldsymbol{\varepsilon}^{\mathtt{c}}) 
\triangleq 
\inf_{\tilde{\boldsymbol{\varepsilon}}\in \mathscr{E}_{0}} \big\langle 
\omega_{\Delta}(\boldsymbol{\varepsilon}^{\mathtt{c}}
+\tilde{\boldsymbol{\varepsilon}}) \big\rangle
=\inf_{\tilde{\boldsymbol{\varepsilon}}\in \mathscr{E}_{0}}\inf_{\boldsymbol{\mu}\in L^{2}}  \big\langle \mathcal{F}_{\Delta}(\boldsymbol{\varepsilon}^{\mathtt{c}}
+\tilde{\boldsymbol{\varepsilon}},\, \boldsymbol{\mu}) \big\rangle,
\end{equation}
$$

以及宏观应力等于

$$
\boldsymbol{\sigma}^{\mathtt{c}}
= \frac{\partial W_{\Delta}}{\partial \boldsymbol{\varepsilon}^{\mathtt{c}}}.
$$

其中，$L^{2}$ 是在单胞内平方可积的本征应变场空间。

## 非均匀变形模式

NTFA 假设单胞域内的状态变量场假设由如下 $M$ 个非均匀变形模式 $\boldsymbol{\mu}^{(k)}$，以及随时间演化的系数 $\mathsf{q}^{(k)}$ 组成：

$$
\begin{equation}\label{eq:ansatz_sv_field}
\begin{aligned}
\text{Ansatz}: \boldsymbol{\mu}(\boldsymbol{y},\, t) 
&= \sum_{k=1}^{M}
\mathsf{q}^{(k)}(t)\boldsymbol{\mu}^{(k)}(\boldsymbol{y})
\triangleq \hat{\boldsymbol{\mu}}(\boldsymbol{y},\, \boldsymbol{\mathsf{q}}(t)),\\
\boldsymbol{\mathsf{q}}(t) &\triangleq \begin{pmatrix}
\mathsf{q}^{(1)}(t)&\mathsf{q}^{(2)}(t)&\cdots\mathsf{q}^{(M)}(t)
\end{pmatrix}^{\top},
\end{aligned}
\end{equation}
$$

并假设模式 $\boldsymbol{\mu}^{(k)}$ 之间线性无关。因此非均匀变形模式在有限维空间 $\text{span}\, \{ \boldsymbol{\mu}^{(k)} \}$ 中，因此式 $\eqref{eq:macro_helmholtz}$ 定义的宏观应变能 $W_{\Delta}$ 有如下不等式约束关系：

$$
\begin{equation}\label{eq:ineq}
\begin{aligned}
W_{\Delta}(\boldsymbol{\varepsilon}^{\mathtt{c}})
&=\inf_{\tilde{\boldsymbol{\varepsilon}}\in \mathscr{E}_{0}} \inf_{\boldsymbol{\mu}\in L^{2}} \big\langle \mathcal{F}_{\Delta}(\boldsymbol{\varepsilon}^{\mathtt{c}}
+\tilde{\boldsymbol{\varepsilon}},\, \boldsymbol{\mu}) \big\rangle\\
&\leq \inf_{\tilde{\boldsymbol{\varepsilon}}\in \mathscr{E}_{0}} \inf_{\boldsymbol{\mathsf{q}}\in\mathbb{R}^{M}}
\big\langle \mathcal{F}_{\Delta}\big(\boldsymbol{\varepsilon}^{\mathtt{c}}
+\tilde{\boldsymbol{\varepsilon}},\, \hat{\boldsymbol{\mu}}(\boldsymbol{\mathsf{q}})\big) \big\rangle\\
&= \inf_{\boldsymbol{\mathsf{q}}\in \mathbb{R}^{M}} 
\Big( \Pi(\boldsymbol{\varepsilon}^{\mathtt{c}},\, \boldsymbol{\mathsf{q}})
+ \Delta t \Phi\big( (\boldsymbol{\mathsf{q}}-\boldsymbol{\mathsf{q}}_{n})/\Delta t \big) \Big)
\triangleq W_{\Delta}^{+}(\boldsymbol{\varepsilon}^{\mathtt{c}}),
\end{aligned}
\end{equation}
$$

其中，

$$
\begin{equation}\label{eq:macro_potential}
\begin{aligned}
\Pi(\boldsymbol{\varepsilon}^{\mathtt{c}},\, \boldsymbol{\mathsf{q}})
&\triangleq \inf_{\tilde{\boldsymbol{\varepsilon}}\in \mathscr{E}_{0}} 
\big\langle 
\omega\big(\boldsymbol{\varepsilon}^{\mathtt{c}} +\tilde{\boldsymbol{\varepsilon}},\, \hat{\boldsymbol{\mu}}(\boldsymbol{\mathsf{q}})\big)
\big\rangle,\\
\Phi\big( (\boldsymbol{\mathsf{q}}-\boldsymbol{\mathsf{q}}_{n})/\Delta t \big)
&\triangleq \big\langle \varphi\big[ (\hat{\boldsymbol{\mu}}(\boldsymbol{\mathsf{q}})-\hat{\boldsymbol{\mu}}(\boldsymbol{\mathsf{q}}_{n}))/\Delta t \big] \big\rangle.
\end{aligned}
\end{equation}
$$

观察 $W_{\Delta}^{+}$ 定义式 $\eqref{eq:ineq}$，可以看到，其保留了细观尺度增量变分原理相似的数学结构 $\eqref{eq:stat_dissip_vari_local}$。同样的，式 $\eqref{eq:macro_potential}$ 也和细观尺度的势函数具有相似的结构（函数 $\Phi$ 应按照 $\Delta t\to 0^{+}$ 的极限情况理解），状态变量从无穷维的状态变量场 $\boldsymbol{\mu}$ 转变为有限维的状态向量 $\boldsymbol{\mathsf{q}}$。同时，式 $\eqref{eq:macro_potential}$ 的变分问题 $\Pi$ 是固定状态变量 $\boldsymbol{\mathsf{q}}$ 后的弹性问题。NTFA 使用 $W_{\Delta}^{+}$ 代替 $W_{\Delta}$，宏观应力和式 $\eqref{eq:ineq}$ 提供的变分原理的极值条件（欧拉--拉格朗日方程）分别对应为

$$
\begin{equation}\label{eq:stat_dissip_strong_macro}
\begin{cases}
\boldsymbol{\sigma}^{\mathtt{c}} 
\simeq \dfrac{\partial W_{\Delta}^{+}}{\partial \boldsymbol{\varepsilon}^{\mathtt{c}}}
= \dfrac{\partial \Pi}{\partial \boldsymbol{\varepsilon}^{\mathtt{c}}}(\boldsymbol{\varepsilon}^{\mathtt{c}},\, \boldsymbol{\mathsf{q}}),\\
\dfrac{\partial \Pi}{\partial \boldsymbol{\mathsf{q}}}(\boldsymbol{\varepsilon}^{\mathtt{c}},\, \boldsymbol{\mathsf{q}})
+ \dfrac{\partial \Phi}{\partial \dot{\boldsymbol{\mathsf{q}}}}\big( \dot{\boldsymbol{\mathsf{q}}} \big)
= 0,
\end{cases}
\end{equation}
$$

将上式与式 $\eqref{eq:stat_dissip_strong_local}$ 的结构进行对比，可以看到两者之间结构的类似。状态变量 $(\boldsymbol{\varepsilon}^{\mathtt{c}},\, \boldsymbol{\mathsf{q}})$，宏观 Helmholtz 自由能 $\Pi$ 和耗散势 $\Phi$，共同组成了非线性复合材料的宏观本构。

### 系数的演化方程

接下来进一步探索 NTFA 中系数 $\boldsymbol{\mathsf{q}}$ 的演化方程。式 $\eqref{eq:macro_potential}$ 可视作含参数 $\boldsymbol{\mathsf{q}}$ 的方程，分别对 $\mathsf{q}^{(k)}$ 求偏导数得到

$$
\frac{\partial \Pi}{\partial \mathsf{q}^{(k)}}(\boldsymbol{\varepsilon}^{\mathtt{c}},\, \boldsymbol{\mathsf{q}})
= \Big\langle 
\frac{\partial \omega}{\partial \boldsymbol{\mu}}:
\boldsymbol{\mu}^{(k)}
\Big\rangle,\quad
\frac{\partial \Phi}{\partial \dot{\mathsf{q}}^{(k)}}( \dot{\boldsymbol{\mathsf{q}}} )
= \Big\langle 
\frac{\partial \varphi}{\partial \dot{\boldsymbol{\mu}}}:
\boldsymbol{\mu}^{(k)}
\Big\rangle,
$$

结合式 $\eqref{eq:stat_dissip_strong_macro}$，实际上得到了如下关于 $\boldsymbol{\mathsf{q}}$ 的常微分方程组：

$$
\begin{equation}\label{eq:reduced_odes}
\begin{aligned}
\frac{\partial \Pi}{\partial \mathsf{q}^{(k)}}(\boldsymbol{\varepsilon}^{\mathtt{c}},\, \boldsymbol{\mathsf{q}})
+ \frac{\partial \Phi}{\partial \dot{\mathsf{q}}^{(k)}}( \dot{\boldsymbol{\mathsf{q}}} )
&= \Big\langle 
\Big( \frac{\partial \omega}{\partial \boldsymbol{\mu}}(\boldsymbol{\varepsilon},\, \boldsymbol{\mu})
+ \frac{\partial \varphi}{\partial \dot{\boldsymbol{\mu}}}(\dot{\boldsymbol{\mu}}) \Big)
: \boldsymbol{\mu}^{(k)} \Big\rangle
\\&= 0,\quad k=1,2,\ldots,M.
\end{aligned}
\end{equation}
$$

### 案例 1：Maxwell 粘弹性本构

以下将以粘弹性材料为例，说明如何通过预先计算给出 $\eqref{eq:reduced_odes}$ 的常微分方程组。设单胞组分满足 Maxwell 粘弹性本构关系：

$$
\begin{equation}\label{eq:consti_visco}
\boldsymbol{\sigma}(\boldsymbol{y}) = \mathbb{L}(\boldsymbol{y})
:\big( \boldsymbol{\varepsilon}(\boldsymbol{y})-\boldsymbol{\varepsilon}_{\mathtt{v}}(\boldsymbol{y}) \big),\quad
\dot{\boldsymbol{\varepsilon}}_{\mathtt{v}}(\boldsymbol{y}) = \mathbb{M}_{\mathtt{v}}(\boldsymbol{y}) : \boldsymbol{\sigma}(\boldsymbol{y}),
\end{equation}
$$

式中，$\boldsymbol{\varepsilon}_{\mathtt{v}}$ 表示粘性应变，$\mathbb{L}$ 是四阶弹性刚度张量，$\mathbb{M}_{\mathtt{v}}$ 也是四阶张量，关联粘性应变率 $\dot{\boldsymbol{\varepsilon}}_{\mathtt{v}}$ 和应力。$\mathbb{L}$ 和 $\mathbb{M}_{\mathtt{v}}$ 是在各个材料相上取常值的分布函数，用材料相 $\alpha$ 的特征函数 $\chi^{(\alpha)}$ 可以表示为

$$
\mathbb{L}(\boldsymbol{y})
= \sum_{\alpha=1}^{N}\chi^{(\alpha)}(\boldsymbol{y})\mathbb{L}^{(\alpha)}, \quad
\mathbb{M}_{\mathtt{v}}(\boldsymbol{y})
= \sum_{\alpha=1}^{N}\chi^{(\alpha)}(\boldsymbol{y})\mathbb{M}_{\mathtt{v}}^{(\alpha)}.
$$

与本构关系 $\eqref{eq:consti_visco}$ 对应的材料相 $\alpha$ 的 Helmholtz 自由能和耗散势函数写成

$$
\begin{equation}\label{eq:visco_potential}
\begin{aligned}
\omega^{(\alpha)}(\boldsymbol{\varepsilon},\, \boldsymbol{\varepsilon}_{\mathtt{v}})
&:= \frac{1}{2} ( \boldsymbol{\varepsilon} - \boldsymbol{\varepsilon}_{\mathtt{v}} )
:\mathbb{L}^{(\alpha)}:( \boldsymbol{\varepsilon} - \boldsymbol{\varepsilon}_{\mathtt{v}} ),\\
\varphi^{(\alpha)}(\dot{\boldsymbol{\varepsilon}}_{\mathtt{v}})
&:= \frac{1}{2} \dot{\boldsymbol{\varepsilon}}_{\mathtt{v}} : \mathbb{L}_{\mathtt{v}}
: \dot{\boldsymbol{\varepsilon}}_{\mathtt{v}},
\end{aligned}
\end{equation}
$$

式中，四阶张量 $\mathbb{L}_{\mathtt{v}}\triangleq(\mathbb{M}_{\mathtt{v}})^{-1}$。在选择式 $\eqref{eq:ansatz_sv_field}$ 作为粘性应变场 $\boldsymbol{\varepsilon}_{\mathtt{v}}$ 的降阶表示后，式 $\eqref{eq:macro_potential}$ 中变分问题

$$
\begin{equation}\label{eq:variational_thermo_elas}
\Pi(\boldsymbol{\varepsilon}^{\mathtt{c}},\, \boldsymbol{\mathsf{q}})
= \inf_{\tilde{\boldsymbol{\varepsilon}}\in \mathscr{E}_{0}} 
\big\langle 
\omega\big(\boldsymbol{\varepsilon}^{\mathtt{c}} +\tilde{\boldsymbol{\varepsilon}},\, \hat{\boldsymbol{\mu}}(\boldsymbol{\mathsf{q}})\big)
\big\rangle
\end{equation}
$$

对应的热弹性方程的解可通过线性的场变换算子 $\mathcal{D}$ 表示为

$$
\begin{equation}\label{eq:fine_scale_strain}
\begin{aligned}
\boldsymbol{\varepsilon}^{\star}(\boldsymbol{y})
&= \boldsymbol{\varepsilon}^{\mathtt{c}}
+ \tilde{\boldsymbol{\varepsilon}}^{\star}(\boldsymbol{y})\\
&= \mathbb{E}(\boldsymbol{y}) : \boldsymbol{\varepsilon}^{\mathtt{c}}
+ \mathcal{D}*\Big( \sum_{k=1}^{M} \mathsf{q}^{(k)}\boldsymbol{\mu}^{(k)} \Big)\\
&= \mathbb{E}(\boldsymbol{y}) : \boldsymbol{\varepsilon}^{\mathtt{c}}
+ \sum_{k=1}^{M} \mathsf{q}^{(k)} \boldsymbol{P}^{(k)}(\boldsymbol{y}),  \quad
\boldsymbol{P}^{(k)}(\boldsymbol{y})\triangleq \mathcal{D}*\boldsymbol{\mu}^{(k)},
\end{aligned}
\end{equation}
$$

式中，$\boldsymbol{P}^{(k)}$ 是变形模式 $\boldsymbol{\mu}^{(k)}$ 在场变换算子作用后得到的影响函数（应变场）。将上述应变场表达式代入到式 $\eqref{eq:consti_visco}$ 中，即得到应力场表达式：

$$
\begin{equation}\label{eq:stress_field}
\begin{aligned}
\boldsymbol{\sigma}(\boldsymbol{y})
=\mathbb{L}(\boldsymbol{y}):\mathbb{E}(\boldsymbol{y}):\boldsymbol{\varepsilon}^{\mathtt{c}}
+ \sum_{k=1}^{M} \mathsf{q}^{(k)} \boldsymbol{A}^{(k)}(\boldsymbol{y}), \\
\boldsymbol{A}^{(k)}(\boldsymbol{y}) \triangleq 
\mathbb{L}(\boldsymbol{y}) : \big( \boldsymbol{P}^{(k)}(\boldsymbol{y})
- \boldsymbol{\mu}^{(k)}(\boldsymbol{y})\big).
\end{aligned}
\end{equation}
$$

结合粘弹性材料具体的势函数表达式 $\eqref{eq:visco_potential}$，将应变场 $\eqref{eq:fine_scale_strain}$ 和本征应变场 $\eqref{eq:ansatz_sv_field}$ 代入到宏观势函数 $\Pi$ 的表达式 $\eqref{eq:macro_potential}$ 中，得到

$$
\begin{equation}\label{eq:reduced_macro_helmholtz}
\begin{aligned}
\Pi(\boldsymbol{\varepsilon}^{\mathtt{c}},\, \boldsymbol{\mathsf{q}})
&= \big\langle 
\omega\big(\boldsymbol{\varepsilon}^{\mathtt{c}} +\tilde{\boldsymbol{\varepsilon}}^{\star},\, \sum_{k=1}^{M} \mathsf{q}^{(k)}\boldsymbol{\mu}^{(k)}\big)
\big\rangle\\
&= \frac{1}{2} \Big\langle \sum_{\alpha=1}^{N}\chi^{(\alpha)}
\Big( \boldsymbol{\varepsilon}^{\mathtt{c}} +\tilde{\boldsymbol{\varepsilon}}^{\star} - \sum_{k=1}^{M} \mathsf{q}^{(k)}\boldsymbol{\mu}^{(k)} \Big)
:\mathbb{L}^{(\alpha)}:\Big( \boldsymbol{\varepsilon}^{\mathtt{c}} +\tilde{\boldsymbol{\varepsilon}}^{\star} - \sum_{k=1}^{M} \mathsf{q}^{(k)}\boldsymbol{\mu}^{(k)} \Big) 
\Big\rangle,
\end{aligned}
\end{equation}
$$

若使用 $\{\bullet\}_{\mathtt{vgt}}$ 表示使用 Voigt 记法将张量映射为列向量或矩阵，并给出如下矩阵记号：

$$
\begin{equation}\label{eq:def_matrix}
\begin{aligned}
\underbrace{\mathbb{L}^{\mathtt{c}}}_{6\times 6} \triangleq \{\langle \mathbb{E}:\mathbb{L}:\mathbb{E} \rangle\}_{\mathtt{vgt}},\quad
\underbrace{\boldsymbol{\mathsf{A}}}_{6\times M}(:,\, k) \triangleq
\{ \langle \boldsymbol{\mu}^{(k)} : \mathbb{L} : \mathbb{E} \rangle \}_{\mathtt{vgt}},\\
\underbrace{\boldsymbol{\mathsf{D}}}_{M\times M}(k,\, l) 
\triangleq \langle \boldsymbol{\mu}^{(k)} : \mathbb{L} : \boldsymbol{P}^{(l)} \rangle, \quad
\underbrace{\boldsymbol{\mathsf{L}}}_{M\times M}(k,\, l) 
\triangleq \langle \boldsymbol{\mu}^{(k)} : \mathbb{L} : \boldsymbol{\mu}^{(l)} \rangle.
\end{aligned}
\end{equation}
$$

为简便起见还是使用张量记号表示四阶张量 $\mathbb{L}^{\mathtt{c}}$，以及之后用到的所有的二阶对称张量，但四阶张量均为 $6\times 6$ 的矩阵，二阶对称张量均为 $6 \times 1$ 的列向量。由此，宏观势函数 $\eqref{eq:reduced_macro_helmholtz}$ 可表示为

$$
\Pi(\boldsymbol{\varepsilon}^{\mathtt{c}},\, \boldsymbol{\mathsf{q}})
= \frac{1}{2}\boldsymbol{\varepsilon}^{\mathtt{c},\top}\mathbb{L}^{\mathtt{c}}\boldsymbol{\varepsilon}^{\mathtt{c}}
- \boldsymbol{\varepsilon}^{\mathtt{c},\top} \boldsymbol{\mathsf{A}}\boldsymbol{\mathsf{q}}
+ \frac{1}{2} \boldsymbol{\mathsf{q}}^{\top} ( \boldsymbol{\mathsf{L}} - \boldsymbol{\mathsf{D}} )\boldsymbol{\mathsf{q}}.
$$

类似的，式 $\eqref{eq:macro_potential}$ 定义的宏观耗散势在代入具体的粘弹性本构 $\eqref{eq:visco_potential}$ 和假设的本征应变场 $\eqref{eq:ansatz_sv_field}$ 之后，得到

$$
\begin{equation}\label{eq:macro_dissip_potential}
\Phi(\dot{\boldsymbol{\mathsf{q}}})
= \frac{1}{2} \dot{\boldsymbol{\mathsf{q}}}^{\top} \boldsymbol{\mathsf{L}}_{\mathtt{v}}\dot{\boldsymbol{\mathsf{q}}}, \quad
\underbrace{\boldsymbol{\mathsf{L}}_{\mathtt{v}}}_{M\times M}(k,\, l)
\triangleq \langle \boldsymbol{\mu}^{(k)} : \mathbb{L}_{\mathtt{v}} : \boldsymbol{\mu}^{(l)} \rangle.
\end{equation}
$$

将 $\boldsymbol{\mathsf{f}}$ 记作与 $\boldsymbol{\mathsf{q}}$ 功共轭的热力学力，那么有

$$
\begin{equation}\label{eq:stat_dissip_reduced_macro}
\begin{cases}
\boldsymbol{\mathsf{f}} = - \dfrac{\partial \Pi}{\partial\boldsymbol{\mathsf{q}} }
= \boldsymbol{\mathsf{A}} \boldsymbol{\varepsilon}^{\mathtt{c}}+(\boldsymbol{\mathsf{D}}-\boldsymbol{\mathsf{L}})\boldsymbol{\mathsf{q}},&\text{state law},\\
\boldsymbol{\mathsf{f}} = \dfrac{\partial \Phi}{\partial\dot{\boldsymbol{\mathsf{q}}} }
= \boldsymbol{\mathsf{L}}_{\mathtt{v}}\dot{\boldsymbol{\mathsf{q}}}, &\text{evolution law},
\end{cases}
\end{equation}
$$

结合上式，得到关于 $\boldsymbol{\mathsf{q}}$ 的常微分方程组为

$$
\begin{equation}\label{eq:ode_macro_stat_vars}
\boldsymbol{\mathsf{L}}_{\mathtt{v}}\dot{\boldsymbol{\mathsf{q}}}
= \boldsymbol{\mathsf{A}} \boldsymbol{\varepsilon}^{\mathtt{c}}+(\boldsymbol{\mathsf{D}}-\boldsymbol{\mathsf{L}})\boldsymbol{\mathsf{q}},
\end{equation}
$$

或者等价地写成关于 $\boldsymbol{\tau}$ 的常微分方程组

$$
\dot{\boldsymbol{\mathsf{f}}}=\boldsymbol{\mathsf{A}} \dot{\boldsymbol{\varepsilon}}^{\mathtt{c}}+(\boldsymbol{\mathsf{D}}-\boldsymbol{\mathsf{L}})(\boldsymbol{\mathsf{L}}_{\mathtt{v}})^{-1}\boldsymbol{\mathsf{f}}.
$$

## 宏观耗散势的近似

### 对偶势函数

在上面给出的粘弹性材料的例子中，式 $\eqref{eq:visco_potential}$ 的耗散势 $\Phi$ 是简单的二次型，因而最后能解析地给出宏观状态变量的演化方程。但在一般的非弹性问题中，宏观耗散势可能无法写成二次型，即使用二次型作近似，类似于式 $\eqref{eq:macro_dissip_potential}$ 中 $\boldsymbol{\mathsf{L}}_{\mathtt{v}}$ 也可能与 $\dot{\boldsymbol{\mathsf{q}}}$ 相关，在求偏导数之后得到的常微分方程组为 $\boldsymbol{g}(\dot{\boldsymbol{\mathsf{q}}})=\boldsymbol{f}(\boldsymbol{\mathsf{q}},\, t)$，其中函数 $\boldsymbol{g}$ 是非线性的。因此接下来的问题是使用耗散势的对偶函数 $\hat{\Phi}^{*}$ 代替 $\Phi$，由此对功共轭的热力学力 $\boldsymbol{\mathsf{f}}$ 求偏导之后，可直接得到 $\dot{\boldsymbol{\mathsf{q}}}$ 的方程组：

$$
\begin{equation}\label{eq:ode_q}
\dot{\boldsymbol{\mathsf{q}}} = \frac{\partial \hat{\Phi}^{*}}{\partial \boldsymbol{\mathsf{f}}}.
\end{equation}
$$

为得到上述函数形式，首先用 Legendre 变换表示耗散势 $\varphi$：

$$
\varphi(\dot{\boldsymbol{\mu}})
= \sup_{\boldsymbol{\tau}\in L^{2}} \big( \dot{\boldsymbol{\mu}} : \boldsymbol{\tau} 
- \varphi^{*}(\boldsymbol{\tau})\big).
$$

再将上式代入到增量变分原理中，求体积平均后得到

$$
\langle \mathcal{F}_{\Delta}(\boldsymbol{\varepsilon},\, \boldsymbol{\mu}) \rangle
= \sup_{\boldsymbol{\tau}\in L^{2}}\Big( \langle \omega(\boldsymbol{\varepsilon},\, \boldsymbol{\mu}) \rangle
+ \langle (\boldsymbol{\mu}-\boldsymbol{\mu}_{n}) : \boldsymbol{\tau} \rangle
-\Delta t \langle \varphi^{*}(\boldsymbol{\tau}) \rangle
\Big),
$$

将假设的本征应变变形模式代入，再求极值，得到（注意其中 $\tilde{\boldsymbol{\varepsilon}}$ 和 $\boldsymbol{\tau}$ 是完全解耦的，所以可以交换 $\inf$ 和 $\sup$ 的顺序）

$$
\begin{equation}\label{eq:temp_a}
\begin{aligned}
&\inf_{\tilde{\boldsymbol{\varepsilon}}\in \mathscr{E}_{0}}
\inf_{\boldsymbol{\mu}\in L^{2}} \langle \mathcal{F}_{\Delta}(\boldsymbol{\varepsilon},\, \boldsymbol{\mu}) \rangle \\
\leq& \inf_{\boldsymbol{\mathsf{q}}\in \mathbb{R}^{M}}
\sup_{\boldsymbol{\tau}\in L^{2}}\Big( \Pi(\boldsymbol{\varepsilon}^{\mathtt{c}},\, \boldsymbol{\mathsf{q}})
+ \sum_{k=1}^{M} (\mathsf{q}^{(k)}-\mathsf{q}_{n}^{(k)}) \langle \boldsymbol{\mu}^{(k)} : \boldsymbol{\tau} \rangle
-\Delta t \langle \varphi^{*}(\boldsymbol{\tau}) \rangle
\Big).
\end{aligned}
\end{equation}
$$

当场变量 $\boldsymbol{\tau}$ 取任意值时，上式不等号的右端项与式 $\eqref{eq:ineq}$ 中不等号的右端项是完全等价的。Fritzen 等人假设场变量 $\boldsymbol{\tau}$ 取值为求解热弹性问题 $\eqref{eq:variational_thermo_elas}$ 得到的自平衡应力场：

$$
\begin{equation}\label{eq:assume_tau}
\boldsymbol{\tau}(\boldsymbol{\varepsilon}^{\mathtt{c}},\, \boldsymbol{\mathsf{p}}) 
:= -\frac{\partial \omega}{\partial \boldsymbol{\mu} }(\boldsymbol{\varepsilon}^{\mathtt{c}}+\tilde{\boldsymbol{\varepsilon}}^{\star},\, \sum_{k=1}^{M}\mathsf{p}^{(k)}\boldsymbol{\mu}^{(k)}),
\end{equation}
$$

注意，这一在有限维空间中的选择意味着，在将上式代入式 $\eqref{eq:temp_a}$ 右端项之后，不再有不等关系。代入之后得到

$$
\begin{equation}\notag
\begin{aligned}
\text{R. H. S. of Eq. }\eqref{eq:temp_a} = \\
\inf_{\boldsymbol{\mathsf{q}}\in \mathbb{R}^{M}}
\sup_{\boldsymbol{\mathsf{p}}\in \mathbb{R}^{M}}\Big( \Pi(\boldsymbol{\varepsilon}^{\mathtt{c}},\, \boldsymbol{\mathsf{q}})
+ \sum_{k=1}^{M} (\mathsf{q}^{(k)}-\mathsf{q}_{n}^{(k)}) \langle \boldsymbol{\mu}^{(k)} : \boldsymbol{\tau}(\boldsymbol{\varepsilon}^{\mathtt{c}},\, \boldsymbol{\mathsf{p}}) \rangle
-\Delta t \big\langle \varphi^{*}\big(\boldsymbol{\tau}(\boldsymbol{\varepsilon}^{\mathtt{c}},\, \boldsymbol{\mathsf{p}})\big) \big\rangle
\Big).
\end{aligned}
\end{equation}
$$

上式关于 $\mathsf{q}^{(k)}$ 的极值条件为（极值点记作 $\check{\boldsymbol{\mathsf{q}}}$）

$$
\begin{equation}\label{eq:cond_min_xi}
\frac{\partial \Pi}{\partial \mathsf{q}^{(k)}}\Big|_{\check{\mathsf{q}}^{(k)}} + \langle \boldsymbol{\mu}^{(k)} : \boldsymbol{\tau} \rangle
= 0,\quad k=1,2,\ldots,M,
\end{equation}
$$

以及关于 $\mathsf{p}^{(l)}$ 的极值条件为

$$
\begin{equation}\label{eq:cond_min_p}
\sum_{k=1}^{M} (\mathsf{q}^{(k)}-\mathsf{q}_{n}^{(k)}) 
\big\langle \boldsymbol{\mu}^{(k)} 
: \frac{\partial \boldsymbol{\tau}}{\partial \mathsf{p}^{(l)}} \big\rangle
=\Delta t \big\langle \frac{\partial \varphi^{*}}{\partial \mathsf{p}^{(l)}} \big\rangle, \quad
l=1,2,\ldots,M.
\end{equation}
$$

接下来将说明当场变量 $\boldsymbol{\tau}$ 按照式 $\eqref{eq:assume_tau}$ 由 $\boldsymbol{\mathsf{p}}$ 给出时，并在进一步假设泛函 $\Pi$ 关于 $\boldsymbol{\mathsf{q}}$ 是**严格凸**的情况下，$\check{\boldsymbol{\mathsf{q}}} = \boldsymbol{\mathsf{p}}$。一方面，根据式 $\eqref{eq:stat_dissip_strong_macro}$，以及假设的 $\eqref{eq:assume_tau}$，可得到广义热力学力

$$
\mathsf{f}^{(k)}(\boldsymbol{\varepsilon}^{\mathtt{c}},\, \boldsymbol{\mathsf{p}})
= -\frac{\partial \Pi}{\partial \mathsf{q}^{(k)}}
= \langle \boldsymbol{\mu}^{(k)} : \boldsymbol{\tau} \rangle,
$$

另一方面，根据极值条件 $\eqref{eq:cond_min_xi}$，

$$
\begin{equation}\label{eq:thermo_force}
\mathsf{f}^{(k)}(\boldsymbol{\varepsilon}^{\mathtt{c}},\, \check{\boldsymbol{\mathsf{q}}})
= -\frac{\partial \Pi}{\partial \mathsf{q}^{(k)}}
= \langle \boldsymbol{\mu}^{(k)} : \boldsymbol{\tau} \rangle,
\end{equation}
$$

应用泛函 $\Pi$ 关于 $\boldsymbol{\mathsf{q}}$ 是严格凸的假设，即可得到要证明的结论。对广义热力学力 $\boldsymbol{\mathsf{f}}$ 进一步关于 $\mathsf{q}^{(l)}$ 求偏导数，结合极值条件 $\eqref{eq:cond_min_p}$，得到

$$
\begin{equation}\label{eq:temp_b}
\frac{\partial \Phi^{*}}{\partial \mathsf{q}^{(l)}}(\boldsymbol{\varepsilon}^{\mathtt{c}},\, \boldsymbol{\mathsf{q}})
= \sum_{k=1}^{M} \dot{\mathsf{q}}^{(k)} \frac{\partial \mathsf{f}^{(k)}}{\partial \mathsf{q}^{(l)}},\quad
\Phi^{*} \triangleq \langle \varphi^{*}(\boldsymbol{\tau}(\boldsymbol{\varepsilon}^{\mathtt{c}},\, \boldsymbol{\mathsf{q}})) \rangle.
\end{equation}
$$

而势函数 $\Phi^{*}$ 中变量 $\boldsymbol{\mathsf{q}}$ 又可以写成关于 $(\boldsymbol{\varepsilon}^{\mathtt{c}},\, \boldsymbol{\mathsf{f}})$ 的函数，$\hat{\Phi}^{*}(\boldsymbol{\varepsilon}^{\mathtt{c}},\, \boldsymbol{\mathsf{f}})\triangleq \Phi^{*}\big( \boldsymbol{\varepsilon}^{\mathtt{c}},\, \boldsymbol{\mathsf{q}}(\boldsymbol{\varepsilon}^{\mathtt{c}},\, \boldsymbol{\mathsf{f}}) \big)$，继续对 $\mathsf{f}^{(k)}$ 求偏导数，得到

$$
\frac{\partial \hat{\Phi}^{*}}{\partial \mathsf{f}^{(k)}}
= \sum_{k=1}^{M} \frac{\partial \Phi^{*}}{\partial \mathsf{q}^{(l)}} 
\frac{\partial \mathsf{q}^{(l)}}{\partial \mathsf{f}^{(k)} },
$$

结合上述两式，即得到这一节开篇给出的关于 $\boldsymbol{\mathsf{q}}$ 的微分方程 $\eqref{eq:ode_q}$。整理式 $\eqref{eq:ode_q}$ 和 $\eqref{eq:thermo_force}$，得到以下广义标准材料的结构：

$$
\begin{equation}\label{eq:gsm}
\begin{cases}
\boldsymbol{\mathsf{f}}
= -\dfrac{\partial \Pi}{\partial \boldsymbol{\mathsf{q}}}(\boldsymbol{\varepsilon}^{\mathtt{c}},\, \boldsymbol{\mathsf{q}}),&\text{state law},\\
\dot{\boldsymbol{\mathsf{q}}} = \dfrac{\partial \hat{\Phi}^{*}}{\partial \boldsymbol{\mathsf{f}}}(\boldsymbol{\varepsilon}^{\mathtt{c}},\, \boldsymbol{\mathsf{f}}), &\text{evolution law},
\end{cases}
\end{equation}
$$

也可以将上式写成关于广义热力学力的常微分方程。对 $\boldsymbol{\mathsf{f}}$ 求时间导数，并注意到根据式 $\eqref{eq:thermo_force}$，$\partial \mathsf{f}^{(k)}/\partial\mathsf{q}^{(l)}= \partial^{2}\Pi/\partial\mathsf{q}^{(k)}\partial\mathsf{q}^{(l)}$ 的指标 $k,l$ 是对称的，再结合式 $\eqref{eq:temp_b}$，就得到

$$
\begin{equation}\label{eq:ode_force}
\begin{aligned}
\dot{\mathsf{f}}^{(k)}
&= \frac{\partial \mathsf{f}^{(k)}}{\partial \boldsymbol{\varepsilon}^{\mathtt{c}}}
: \dot{\boldsymbol{\varepsilon}}^{\mathtt{c}}
+ \sum_{l=1}^{M} \frac{\partial \mathsf{f}^{(k)}}{\partial \mathsf{q}^{(l)}}
\dot{\mathsf{q}}^{(l)}\\
&= \frac{\partial \mathsf{f}^{(k)}}{\partial \boldsymbol{\varepsilon}^{\mathtt{c}}}
: \dot{\boldsymbol{\varepsilon}}^{\mathtt{c}}
+ \frac{\partial \Phi^{*}}{\partial \mathsf{q}^{(k)}}(\boldsymbol{\varepsilon}^{\mathtt{c}},\, \boldsymbol{\mathsf{q}}(\boldsymbol{\varepsilon}^{\mathtt{c}},\, \boldsymbol{\mathsf{f}})).
\end{aligned}
\end{equation}
$$

### 案例 2：有线性运动硬化的粘塑性本构

带有线性运动硬化的粘塑性本构的状态变量为 $(\boldsymbol{\varepsilon}_{\mathtt{v}},\, p)$，其中 $\boldsymbol{\varepsilon}_{\mathtt{v}}$ 是粘性应变，$p$ 是塑性乘子。对应的热力学力分别为 $\boldsymbol{\tau}_{\mathtt{v}}\triangleq\boldsymbol{\sigma}-\boldsymbol{\beta}$ 和 $\tau_{\mathtt{p}}\triangleq-K(p)$，其中 $\boldsymbol{\beta}$ 是背应力，$K(p)$ 是屈服应力。材料的 Helmholtz 自由能为

$$
\begin{equation}\label{eq:visco_plas_helm}
\begin{aligned}
\omega(\boldsymbol{\varepsilon},\, \boldsymbol{\varepsilon}_{\mathtt{v}},\, p)
&:= \frac{1}{2}(\boldsymbol{\varepsilon}-\boldsymbol{\varepsilon}_{\mathtt{v}})
: \mathbb{L} : (\boldsymbol{\varepsilon}-\boldsymbol{\varepsilon}_{\mathtt{v}})
+ \frac{1}{2} \boldsymbol{\varepsilon}_{\mathtt{v}} : \mathbb{H} : \boldsymbol{\varepsilon}_{\mathtt{v}}\\
&+ \omega_{\mathtt{ih}}(p),\quad
\omega_{\mathtt{ih}}(p) \triangleq 
\int_{0}^{p} R(q)\,\mathrm{d}q,
\end{aligned}
\end{equation}
$$

以及热力学力的势函数为

$$
\begin{equation}\label{eq:visco_plas_dissp}
\varphi^{*}(\boldsymbol{\tau}_{\mathtt{v}},\, \tau_{\mathtt{p}})
:= \frac{\sigma_{0}\dot{\varepsilon}_{0}}{n+1}
\bigg[ \frac{\big( (\boldsymbol{\tau}_{\mathtt{v}})_{\mathtt{eq}}
+ \tau_{\mathtt{p}} \big)^{+}}{\sigma_{0}} \bigg]^{n+1},
\end{equation}
$$

关于状态变量场，假设的变形模式分别为（注意，塑性因子需要单独使用特征函数插值）

$$
\begin{equation}\label{eq:inelas_pattern}
\boldsymbol{\varepsilon}_{\mathtt{v}}(\boldsymbol{y})
= \sum_{k=1}^{M} \mathsf{q}_{\mathtt{v}}^{(k)}\boldsymbol{\mu}^{(k)}(\boldsymbol{y}),\quad
p(\boldsymbol{y})
= \sum_{k=1}^{N} \mathsf{p}^{(k)} \chi^{(k)}(\boldsymbol{y}),
\end{equation}
$$

记状态变量 $\boldsymbol{\mathsf{q}}\triangleq(\boldsymbol{\mathsf{q}}_{\mathtt{v}}^{\top},\, \boldsymbol{\mathsf{p}}^{\top})^{\top}\in\mathbb{R}^{M P}$，将假设的变形模式 $\eqref{eq:inelas_pattern}$，以及材料的 Helmholtz 自由能 $\eqref{eq:visco_plas_helm}$ 代入到宏观势函数 $\Pi$ 的表达式 $\eqref{eq:macro_potential}$ 中，得到

$$
\begin{equation}\label{eq:marco_pi}
\Pi(\boldsymbol{\varepsilon}^{\mathtt{c}},\, \boldsymbol{\mathsf{q}})
= \frac{1}{2}\boldsymbol{\varepsilon}^{\mathtt{c},\top}\mathbb{L}^{\mathtt{c}}\boldsymbol{\varepsilon}^{\mathtt{c}}
- \boldsymbol{\varepsilon}^{\mathtt{c},\top} \boldsymbol{\mathsf{A}}\boldsymbol{\mathsf{q}}_{\mathtt{v}}
+ \frac{1}{2} \boldsymbol{\mathsf{q}}_{\mathtt{v}}^{\top} 
( \boldsymbol{\mathsf{L}} + \boldsymbol{\mathsf{H}} - \boldsymbol{\mathsf{D}} )\boldsymbol{\mathsf{q}}_{\mathtt{v}}
+ \boldsymbol{\mathsf{c}}^{\top} \boldsymbol{\mathsf{w}}_{\mathtt{ih}}(\boldsymbol{\mathsf{p}}),
\end{equation}
$$

其中，

$$
\begin{aligned}
\underbrace{\boldsymbol{\mathsf{H}}}_{M\times M}(k,\, l) 
\triangleq \langle \boldsymbol{\mu}^{(k)} : \mathbb{H} : \boldsymbol{\mu}^{(l)} \rangle,\quad
\underbrace{\boldsymbol{\mathsf{c}}}_{N\times 1}
\triangleq \begin{pmatrix}
c^{(1)} & c^{(2)} & \cdots & c^{(N)}
\end{pmatrix}^{\top},\\
\underbrace{\boldsymbol{\mathsf{w}}_{\mathtt{ih}}}_{N\times 1}
\triangleq \begin{pmatrix}
\omega_{\mathtt{ih}}^{(1)}(\mathsf{p}^{(1)}) & \omega_{\mathtt{ih}}^{(2)}(\mathsf{p}^{(2)}) & \cdots & \omega_{\mathtt{ih}}^{(N)}(\mathsf{p}^{(N)})
\end{pmatrix}^{\top}.
\end{aligned}
$$

热力学力对应为

$$
\begin{aligned}
\boldsymbol{\mathsf{f}}_{\mathtt{v}}
&= -\frac{\partial \Pi}{\partial \boldsymbol{\mathsf{q}}_{\mathtt{v}}}
= \boldsymbol{\mathsf{A}} \boldsymbol{\varepsilon}^{\mathtt{c}}
+ (\boldsymbol{\mathsf{D}} - \boldsymbol{\mathsf{L}} - \boldsymbol{\mathsf{H}})\boldsymbol{\mathsf{q}}_{\mathtt{v}},\\
\boldsymbol{\mathsf{f}}_{\mathtt{p}}
&= -\frac{\partial \Pi}{\partial \boldsymbol{\mathsf{p}}}\\
&= -\begin{pmatrix}
c^{(1)}K^{(1)}(\mathsf{p}^{(1)}) & c^{(2)}K^{(2)}(\mathsf{p}^{(2)}) & \cdots & c^{(N)}K^{(N)}(\mathsf{p}^{(N)})
\end{pmatrix}^{\top}
\end{aligned}
$$

将假设的变形模式 $\eqref{eq:inelas_pattern}$ 代入到热力学力中，得到

$$
\begin{aligned}
\boldsymbol{\tau}_{\mathtt{v}}(\boldsymbol{y},\, \boldsymbol{\varepsilon}^{\mathtt{c}},\, \boldsymbol{\mathsf{q}}_{\mathtt{v}})
&:= \mathbb{L}(\boldsymbol{y}) : \mathbb{E}(\boldsymbol{y}):\boldsymbol{\varepsilon}^{\mathtt{c}}
+ \sum_{k=1}^{M} \mathsf{q}_{\mathtt{v}}^{(k)} \big( \boldsymbol{A}^{(k)}(\boldsymbol{y}) - \mathbb{H}(\boldsymbol{y}):\boldsymbol{\mu}^{(k)}(\boldsymbol{y}) \big),\\
\tau_{\mathtt{p}}(\boldsymbol{y},\, \boldsymbol{\mathsf{p}})
&:= -\sum_{k=1}^{M}\chi^{(k)}(\boldsymbol{y}) K^{(k)}(\mathsf{p}^{(k)})
= \sum_{k=1}^{M} \chi^{(k)}(\boldsymbol{y}) \mathsf{f}_{\mathtt{v}}^{(k)}(\mathsf{p}^{(k)})/c^{(k)}
\end{aligned}
$$

相应的，根据定义式 $\eqref{eq:temp_b}$，并代入具体的材料耗散势 $\eqref{eq:visco_plas_dissp}$，宏观耗散势 $\Phi^{*}$ 可写成

$$
\begin{equation}\label{eq:macro_dissp}
\Phi^{*}(\boldsymbol{\varepsilon}^{\mathtt{c}},\, \boldsymbol{\mathsf{q}})
= \Big\langle \frac{\sigma_{0}\dot{\varepsilon}_{0}}{n+1}
\bigg[ \frac{\big( (\boldsymbol{\tau}_{\mathtt{v}}(\boldsymbol{\varepsilon}^{\mathtt{c}},\, \boldsymbol{\mathsf{q}}_{\mathtt{v}}))_{\mathtt{eq}}
+ \tau_{\mathtt{p}}(\boldsymbol{\mathsf{p}}) \big)^{+}}{\sigma_{0}} \bigg]^{n+1} \Big\rangle,
\end{equation}
$$

然而，观察到材料的耗散势 $\eqref{eq:visco_plas_dissp}$ 中，由（1）von Mises 等效应力运算 $(\bullet)_{\mathtt{eq}}$、（2）Macaulay 括号 $(\bullet)^{+}$ 和（3）指数 $n$ 引入的非线性运算，导致宏观耗散势 $\eqref{eq:macro_dissp}$ 无法像宏观自由能 $\eqref{eq:marco_pi}$ 一样，通过预先计算得到的矩阵表示为关于广义热力学力 $\boldsymbol{\mathsf{f}}$​ 的二次型。根据式 $\eqref{eq:ode_force}$，得到的热力学常微分方程组仍表示为细观场的体积积分形式：

$$
\begin{equation}\label{eq:ode_fv}
\dot{\mathsf{f}}_{\mathtt{v}}^{(k)}
= \frac{\partial \mathsf{f}_{\mathtt{v}}^{(k)}}{\partial \boldsymbol{\varepsilon}^{\mathtt{c}}}
: \dot{\boldsymbol{\varepsilon}}^{\mathtt{c}}
+ \Big\langle \frac{\partial \varphi^{*}}{\partial \boldsymbol{\tau}_{\mathtt{v}}} : \boldsymbol{\mu}^{(k)} \Big\rangle, \quad
k=1,2,\ldots,M.
\end{equation}
$$

根据式 $\eqref{eq:gsm}$，塑性因子 $\boldsymbol{\mathsf{p}}$ 的演化方程可以表示为

$$
\begin{equation}\label{eq:ode_p}
\begin{aligned}
\dot{\mathsf{p}}^{(k)}
&= \Big\langle \frac{\partial \varphi^{*}}{\partial \tau_{\mathtt{p}}} \frac{\partial \tau_{\mathtt{p}}}{\partial \mathsf{f}_{\mathtt{p}}^{(k)}}\Big\rangle
= \frac{1}{c^{(k)}} \Big\langle \frac{\partial \varphi^{*}}{\partial \tau_{\mathtt{p}}} \chi^{(k)}\Big\rangle\\
&=\Big\langle \frac{\partial \varphi^{*}}{\partial \tau_{\mathtt{p}}}\Big\rangle^{(k)},\quad
k=1,2,\ldots,N.
\end{aligned}
\end{equation}
$$

## 对 NTFA 的评论

正面的评价：

1. 相对 TFA 的修正方法，NTFA 的数学结构更加严格。NTFA 建立在变分原理（式 $\eqref{eq:ineq}$）之上，增加变形模式数量总会得到更准确的结果
2. 线上计算量很小。NTFA 线上求解关于系数 $\boldsymbol{\mathsf{q}}$ 的维度为 $M$ 的常微分方程组（式 $\eqref{eq:ode_macro_stat_vars}$），若线下阶段选择恰当的变形模式，这个常微分方程组规模可以很小，同时得到的结果远比同规模无修正 TFA 精确
3. 数学结构很漂亮。NTFA 完美地具像化广义标准材料模型，是力学理论和计算方法结合的优美的例子

负面的评价：

数学过于复杂。NTFA 需要（1）写出耗散势函数 $\Phi$（2）提取模态。其困难如下

1. 首先，对一般的非线性材料，可能就不存在耗散势
2. 其次，假如存在这样的耗散势，但具有强非线性如式 $\eqref{eq:visco_plas_dissp}$，线上阶段的 ODE 还是要对细观物理场作体积积分（见式 $\eqref{eq:ode_fv}$ 和 $\eqref{eq:ode_p}$），因此需要线上存储非弹性变形模式、隐式求解 ODE 又需要在每一次的迭代步中获取式 $\eqref{eq:ode_fv}$ 和 $\eqref{eq:ode_p}$ 的积分项，这会消耗大量的内存和 CPU 时间
3. 然后，假如按照 Suquet 等人提出的对耗散势作 Taylor 展开并使用二阶项近似（在本文档中没有推导这部分的理论），这能将耗散势 ODE 也表示为类似于宏观自由能的结构，但二阶展开的理论也还是复杂了，这也要求 NTFA 对不同组合的非线性材料分别进行讨论
4. 最后，怎么选取合适的模态？从这个文档推导过程可以看到，（1）每一个状态变量一般都要为它设计一组模态（2）如果模态提取得不好就无法得到严格凸的泛函 $\Pi$（式 $\eqref{eq:cond_min_xi}$）。即使 NTFA 数学结构既严谨又优美，最后选取模态又不得不加入一些 ad hoc 或者人为的 criteria 等等