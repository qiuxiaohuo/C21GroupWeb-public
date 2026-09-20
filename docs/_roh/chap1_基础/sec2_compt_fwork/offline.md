# 线下求解影响张量



在对方形单胞相对面 $S^{-}$，$S^{+}$  施加周期性边界条件时，方形单胞的相对面节点位移满足如下约束：

$$
\boldsymbol{d}_{a}^{+}-\boldsymbol{d}_{a}^{-} = \boldsymbol{d}_{b}^{+} - \boldsymbol{d}_{b}^{-}, \quad
\forall ~\text{node}~ a,b \in S
$$

### 单胞问题主从节点对的确定

在求解影响张量时，施加在单胞上的宏观应变等于零，这意味着对于方形单胞，8 个角点的位移全部等于零：

$$
\boldsymbol{d}_{a}=0, \quad \text{node}\ a \text{ at corner}
$$

正是因为宏观应变等于零，所以对任意相对面上的节点，位移均相等：

$$
\begin{equation}
\boldsymbol{d}_{a}^{+} = \boldsymbol{d}_{a}^{-}, \quad \forall ~\text{node}~ a,b \in S
\end{equation}
$$

因此，确定约束矩阵 $\boldsymbol{\mathsf{T}}$，就转变为寻找所有在 (1) 边界面 和 (2) 棱线 处的节点对 $\{ a^{+}, a^{-} \}$，$\{a^{++}, a^{+-}, a^{-+}, a^{--}\}$，并在节点对中确定一个节点作为主节点，其它所有节点作为从节点，将从节点的三个自由度全部设置为属性 `slave`，将在之后调整线性方程组时消去该自由度。

### 线性方程组的调整

首先，对于属于位移边界条件的自由度，在组装全局刚度矩阵时直接跳过，由此得到数量为 `NDOF-NdcDOF` 的线性方程组：

$$
\boldsymbol{\mathsf{K}}_{uu} \boldsymbol{\mathsf{d}}_{u} = \boldsymbol{\mathsf{f}}_{u}
$$

将未知量 $\boldsymbol{\mathsf{d}}_{u}$ 分解为 $(\hat{\boldsymbol{\mathsf{d}}}; \boldsymbol{\mathsf{d}}_{s})$，其中 $\hat{\boldsymbol{\mathsf{d}}}$ 表示独立的自由度，维度为 `NiDOF x 1`，$\boldsymbol{\mathsf{d}}_{s}$ 表示从属自由度，维度为 `NscDOF x 1`。对于式 (4) 中的约束关系，可以定义维度为 `NscDOF x NiDOF` 矩阵 $\boldsymbol{\mathsf{A}}$，在第 `sidx` 行处只在列索引为 `s2mMap(sidx)` 处等于1，其它位置均等于零，`s2mMap` 将从属自由度索引映射到对应的主自由度索引。由此，式 (4) 约束关系，以及约束矩阵 $\boldsymbol{\mathsf{T}}$ 可以表示为

$$
\boldsymbol{\mathsf{d}}_{s} = \boldsymbol{\mathsf{A}} \hat{\boldsymbol{\mathsf{d}}}, \quad
\boldsymbol{\mathsf{T}} = \begin{bmatrix}\boldsymbol{\mathsf{I}} \\ \boldsymbol{\mathsf{A}} \end{bmatrix}
$$

所以，添加约束关系后的刚度矩阵和外力列阵等于

$$
\hat{\boldsymbol{\mathsf{K}}} = \begin{bmatrix}\boldsymbol{\mathsf{I}} & \boldsymbol{\mathsf{A}} ^{\top}\end{bmatrix} \boldsymbol{\mathsf{K}} \begin{bmatrix}\boldsymbol{\mathsf{I}} \\ \boldsymbol{\mathsf{A}} \end{bmatrix}, \quad 
\hat{\boldsymbol{\mathsf{f}}} = \begin{bmatrix}\boldsymbol{\mathsf{I}} & \boldsymbol{\mathsf{A}} ^{\top}\end{bmatrix} \boldsymbol{\mathsf{f}}.
$$

此时，$\boldsymbol{\mathsf{A}}^{\top}$ 的作用相当于将索引为 `sidx` 的从属自由度组装在 `s2mMap(sidx)` 位置。

