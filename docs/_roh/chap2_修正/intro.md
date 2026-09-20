# 影响张量的修正框架

在上一章可以看到，通过渐进均质化，ROH 可以方便地拓展到三尺度或高阶均质化问题；通过 TFA，ROH 可准确求解热弹性问题，近似给出复杂加载路径下具有材料非线性的单胞的宏观响应，增加分块数量还可进一步提升计算精度；ROH 的实现也可完全嵌入在已有商业软件框架中，例如 MDS。作为一种多尺度计算框架，ROH 虽然具有上述种种优势，但是它在处理材料非线性时也承担了和 TFA 同样的问题：宏观响应过迟和或刚，这也是 Michel 和 Suquet 提出非均匀场变换分析的主要动机。

下面以双相纤维增强型单胞为例，具体介绍 ROH 为什么会有宏观响应过迟或过刚的问题。



在未进行数值离散之前，单胞方程的弱形式和 Green 函数表示的方程解是完全等价的。本论文现在希望证明的是，在数值离散之后，使用单个单元单分块的 ROH 得到的计算结果，和使用常应变单元的 FEM 是相同的。

考虑单胞区域经三角形或四面体单元离散后得到的区域为 $\theta\triangleq\{ \Theta^{(e)} \}$，其单元数量为 $N_{\mathtt{elem}}$。单胞由宏观应变加载，第 $n$ 个增量步取值为 $\boldsymbol{\varepsilon}_{n}^{\mathtt{c}}$，$n=1,2,\ldots$。因为 ROH 和 FEM 分块和单元数量相同，且在单胞域内占据相同空间，之后统一使用 $e$ 表示分块或单元的编号。ROH 分块 $e$ 在第 n 个增量步的应力和应变响应分别记作 $\hat{\boldsymbol{\sigma}}_{n}^{(e)}$ 和 $\hat{\boldsymbol{\varepsilon}}_{n}^{(e)}$，FEM 单元 $e$ 在第 n 个增量步的应力和应变响应分别记作 $\boldsymbol{\sigma}_{n}^{(e)}$ 和 $\boldsymbol{\varepsilon}_{n}^{(e)}$。设 ROH 和 FEM 在第 $n$ 个增量步时所有分块和单元的应力应变响应相同，也即

$$
\begin{equation}\notag
\hat{\boldsymbol{\varepsilon}}_{n}^{(e)}=\boldsymbol{\varepsilon}_{n}^{(e)}, \quad
\hat{\boldsymbol{\sigma}}_{n}^{(e)}=\boldsymbol{\sigma}_{n}^{(e)}, \quad
e=1,2,\ldots,N_{\mathtt{elem}}.
\end{equation}
$$

现证明在第 $n+1$ 个增量步时，ROH 和 FEM 的响应增量也是相同的，也即

$$
\begin{equation}\label{eq:roh_fem_equiv}
\Delta\hat{\boldsymbol{\varepsilon}}_{n+1}^{(e)}=\Delta\boldsymbol{\varepsilon}_{n+1}^{(e)}, \quad
\Delta\hat{\boldsymbol{\sigma}}_{n+1}^{(e)}=\Delta\boldsymbol{\sigma}_{n+1}^{(e)}, \quad
e=1,2,\ldots,N_{\mathtt{elem}}.
\end{equation}
$$

一方面，使用常应变单元的 FEM，单胞方程弱形式的位移场测试函数空间 $H_{\mathtt{per}}^{1}$ 由在单元 $e$ 内线性分布、在单元边界连续、在单胞边界满足周期性条件的子空间 $\mathscr{U}_{\mathtt{per}}^{\theta}$ 近似。空间 $\mathscr{U}_{\mathtt{per}}^{\theta}$ 中的函数在求偏导之后，组成的应变场测试空间记作 $\mathscr{E}_{\mathtt{CSE}}^{\theta}$，其中的函数 $\boldsymbol{\varepsilon}_{\mathtt{CSE}}$ 具有如下性质：

1. 应变场 $\boldsymbol{\varepsilon}_{\mathtt{CSE}}$ 是相容的；
2. 应变场 $\boldsymbol{\varepsilon}_{\mathtt{CSE}}$ 是扰动应变场，也即 $\langle \boldsymbol{\varepsilon}_{\mathtt{CSE}} \rangle=\boldsymbol{0}$；
3. 应变场 $\boldsymbol{\varepsilon}_{\mathtt{CSE}}$ 在每个单元 $e$ 是常值，记作 $\boldsymbol{\varepsilon}_{\mathtt{CSE}}^{(e)}$。

因此，使用常应变有限元方法离散之后的单胞方程弱形式为：寻找 $\Delta\boldsymbol{\varepsilon}_{n+1}\in\mathscr{E}_{\mathtt{CSE}}$，使得

$$
\begin{equation}\label{eq:weak_discret_fem}
\int_{\Theta} 
\Delta\boldsymbol{\sigma}_{n+1}(\Delta\boldsymbol{\varepsilon}_{n+1}) :\delta\boldsymbol{\varepsilon}
\,\mathrm{d}\boldsymbol{y} = 0, \quad
\forall \delta\boldsymbol{\varepsilon} \in \mathscr{E}_{\mathtt{CSE}}.
\end{equation}
$$

另一方面，ROH 给出关于分块应变增量的代数方程组为

$$
\begin{equation}\label{eq:roh_discret}
\Delta\hat{\boldsymbol{\varepsilon}}_{n+1}^{(e)}
= \mathbb{E}^{(e)} : \Delta\boldsymbol{\varepsilon}_{n+1}^{\mathtt{c}}
+ \sum_{a=1}^{N_{\mathtt{elem}}}\mathbb{P}^{(ea)}:\Delta\hat{\boldsymbol{\mu}}_{n+1}^{(a)}, \quad
e=1,2,\ldots,N_{\mathtt{elem}}.
\end{equation}
$$

注意到影响函数 $\mathbb{E}$ 和 $\mathbb{P}^{(a)}$ 通过后处理单胞内弹性方程的应变场获得，上式给出的应变场 $\Delta\hat{\boldsymbol{\varepsilon}}_{n+1}$ 因此也是相容的，所以有

$$
\begin{equation}\label{eq:roh_compatiable}
\Delta\hat{\boldsymbol{\varepsilon}}_{n+1} \in \mathscr{E}_{\mathtt{CSE}}.
\end{equation}
$$

下面将把 ROH 给出的非线性方程组 $\eqref{eq:roh_discret}$ 转换成类似于式 $\eqref{eq:weak_discret_fem}$ 的形式。首先将本征应变和应力的关系式

$$
\begin{equation}\notag
\Delta\hat{\boldsymbol{\mu}}_{n+1}^{(a)}
= \Delta\hat{\boldsymbol{\varepsilon}}_{n+1}^{(a)}
- \mathbb{M}^{(a)} : \Delta\hat{\boldsymbol{\sigma}}_{n+1}^{(a)}
\end{equation}
$$

代入方程组中，得到

$$
\begin{equation}\label{eq:roh_temp}
\begin{aligned}
&\sum_{a=1}^{N_{\mathtt{elem}}}\Big( \mathbb{P}^{(ea)}:\mathbb{M}^{(a)} : \Delta\hat{\boldsymbol{\sigma}}_{n+1}^{(a)}  \Big)\\
=&\, \mathbb{E}^{(e)} : \Delta\boldsymbol{\varepsilon}_{n+1}^{\mathtt{c}}
+ \sum_{a=1}^{N_{\mathtt{elem}}}\Big( \mathbb{P}^{(ea)}:\Delta\hat{\boldsymbol{\varepsilon}}_{n+1}^{(a)} \Big)
- \Delta\hat{\boldsymbol{\varepsilon}}_{n+1}^{(e)}, \quad
e=1,2,\ldots,N_{\mathtt{elem}}.
\end{aligned}
\end{equation}
$$

然后注意到，对任意相容的应变场 $\boldsymbol{\varepsilon}_{\mathtt{cp}}\in\mathscr{E}_{\mathtt{CSE}}$，总满足

$$
\begin{equation}\label{eq:roh_cp_strain_identity}
\boldsymbol{\varepsilon}_{\mathtt{cp}}^{(e)}
= \mathbb{E}^{(e)} : \langle \boldsymbol{\varepsilon}_{\mathtt{cp}} \rangle
+ \sum_{a=1}^{N_{\mathtt{elem}}} \mathbb{P}^{(ea)} : \boldsymbol{\varepsilon}_{\mathtt{cp}}^{(a)},\quad
e = 1,2,\ldots,N_{\mathtt{elem}},
\end{equation}
$$

这一组恒等式可通过如下简单的论证得到：设置本征应变场 $\boldsymbol{\mu}=\boldsymbol{\varepsilon}_{\mathtt{cp}}$，由此得到单胞内应力场 $\boldsymbol{\sigma}=\mathbb{L}:(\boldsymbol{\varepsilon}_{\mathtt{cp}}-\boldsymbol{\mu})\equiv \boldsymbol{0}$，因此总满足平衡方程。由此再结合式 $\eqref{eq:roh_compatiable}$，就得到式 $\eqref{eq:roh_temp}$ 右端项恒等于零，也即

$$
\begin{equation}\label{eq:roh_temp_b}
\sum_{a=1}^{N_{\mathtt{elem}}}\Big( \mathbb{P}^{(ea)}:\mathbb{M}^{(a)} : \Delta\hat{\boldsymbol{\sigma}}_{n+1}^{(a)}  \Big)\\
= \boldsymbol{0}, \quad
e=1,2,\ldots,N_{\mathtt{elem}}.
\end{equation}
$$

再应用互易原理提供的影响张量恒等式

$$
\begin{equation}\notag
c^{(e)} \mathbb{P}^{(ea)} : \mathbb{M}^{(a)}
= \big( c^{(a)} \mathbb{P}^{(ae)} : \mathbb{M}^{(e)} \big)^{\top},
\end{equation}
$$

其中 $\mathbb{M}^{(e)}$ 是单元 $e$ 的弹性柔度张量。代入式 $\eqref{eq:roh_temp_b}$ 中，移项后得到

$$
\begin{equation}\label{eq:roh_temp_c}
\Big( \sum_{a=1}^{N_{\mathtt{elem}}}
c^{(a)} \Delta\hat{\boldsymbol{\sigma}}_{n+1}^{(a)} : \mathbb{P}^{(ae)}  \Big)
: \mathbb{M}^{(e)}
= \boldsymbol{0}, \quad
e=1,2,\ldots,N_{\mathtt{elem}}.
\end{equation}
$$

由于影响函数和应力场在分块上分片常值分布，上式括号项可重新写成

$$
\begin{equation}\notag
\begin{aligned}
\sum_{a=1}^{N_{\mathtt{elem}}}
c^{(a)} \Delta\hat{\boldsymbol{\sigma}}_{n+1}^{(a)} : \mathbb{P}^{(ae)}
&= \frac{1}{|\Theta|}\sum_{a=1}^{N_{\mathtt{elem}}}\int_{\Theta^{(e)}} \Delta\hat{\boldsymbol{\sigma}}_{n+1}^{(a)} : \mathbb{P}^{(ae)} \,\mathrm{d}\boldsymbol{y}\\
&= \frac{1}{|\Theta|}\int_{\Theta} \Delta\hat{\boldsymbol{\sigma}}_{n+1} : \mathbb{P}^{(e)} \,\mathrm{d}\boldsymbol{y},
\end{aligned}
\end{equation}
$$

再根据单元 $e$ 的弹性柔度张量的对称正定性，式 $\eqref{eq:roh_temp_c}$ 可重新写成

$$
\begin{equation}\label{eq:roh_temp_d}
\int_{\Theta} \Delta\hat{\boldsymbol{\sigma}}_{n+1} : \mathbb{P}^{(e)} \,\mathrm{d}\boldsymbol{y}
= \boldsymbol{0}, \quad
e=1,2,\ldots,N_{\mathtt{elem}}.
\end{equation}
$$

选取二阶对称张量的一组线性无关基 $\{\boldsymbol{e}_{i}\}$，$i=1,2,\ldots,N_{\mathtt{dim}}(N_{\mathtt{dim}}+1)/2$。记由影响函数 $\mathbb{P}^{(e)}:\boldsymbol{e}_{i}$ 张成的应变场空间为 $\mathscr{E}_{\mathtt{ROH}}^{\theta}$，由此式 $\eqref{eq:roh_temp_d}$ 可以写成和有限元离散弱形式 $\eqref{eq:weak_discret_fem}$ 类似的表达式：

$$
\begin{equation}\label{eq:weak_discret_roh}
\int_{\Theta} 
\Delta\hat{\boldsymbol{\sigma}}_{n+1}(\Delta\hat{\boldsymbol{\varepsilon}}_{n+1}) :\delta\hat{\boldsymbol{\varepsilon}}
\,\mathrm{d}\boldsymbol{y} = 0, \quad
\forall \delta\hat{\boldsymbol{\varepsilon}} \in \mathscr{E}_{\mathtt{ROH}}^{\theta}.
\end{equation}
$$

由于 ROH 在分块上使用组分自身的本构更新应力，驱动量为分块的平均应变，所以 ROH 和 FEM 在积分点处应力关于应变增量的线性化结果相同。因此证明式 $\eqref{eq:roh_fem_equiv}$ 提出的 ROH 和 FEM 增量响应相同，就等价于证明式 $\eqref{eq:weak_discret_fem}$ 和 $\eqref{eq:weak_discret_roh}$ 的应变场测试空间相同：

$$
\begin{equation}\label{eq:cse_eq_roh}
\mathscr{E}_{\mathtt{CSE}}^{\theta} = \mathscr{E}_{\mathtt{ROH}}^{\theta}.
\end{equation}
$$

一方面，任意由 ROH 影响函数组合得到的应变场 $\hat{\boldsymbol{\varepsilon}}$ 是相容应变场，在单元上分片常值，并且在单胞上体积平均等于零，因此 $\hat{\boldsymbol{\varepsilon}}\in\mathscr{E}_{\mathtt{CSE}}^{\theta}$，所以有

$$
\begin{equation}\label{eq:roh_subset_cse}
\mathscr{E}_{\mathtt{ROH}}^{\theta} \subseteq \mathscr{E}_{\mathtt{CSE}}^{\theta}.
\end{equation}
$$

另一方面，对任意相容的应变扰动场 $\boldsymbol{\varepsilon}\in\mathscr{E}_{\mathtt{CSE}}^{\theta}$，由场变换算子 $\mathcal{D}$ 的性质可得 $\boldsymbol{\varepsilon} = \mathcal{D} * \boldsymbol{\varepsilon}$，而应变场 $\boldsymbol{\varepsilon}$ 在单元内取常值，所以该式又可以写成

$$
\begin{equation}\notag
\boldsymbol{\varepsilon}(\boldsymbol{y})
= \mathcal{D} * \boldsymbol{\varepsilon}
= \sum_{e=1}^{N_{\mathtt{elem}}} \mathbb{P}^{(e)}(\boldsymbol{y}):\boldsymbol{\varepsilon}^{(e)},
\end{equation}
$$

也即应变场 $\boldsymbol{\varepsilon}$ 可表示为影响函数 $\mathbb{P}^{(e)}$ 的线性组合。由 $\boldsymbol{\varepsilon}$ 的任意性可得

$$
\begin{equation}\label{eq:cse_subset_roh}
\mathscr{E}_{\mathtt{CSE}}^{\theta} \subseteq \mathscr{E}_{\mathtt{ROH}}^{\theta} .
\end{equation}
$$

综合式 $\eqref{eq:roh_subset_cse}$ 和式 $\eqref{eq:cse_subset_roh}$ 可得 $\mathscr{E}_{\mathtt{CSE}}^{\theta} = \mathscr{E}_{\mathtt{ROH}}^{\theta}$。再结合之前的论述，公式 $\eqref{eq:roh_fem_equiv}$ 得证。

