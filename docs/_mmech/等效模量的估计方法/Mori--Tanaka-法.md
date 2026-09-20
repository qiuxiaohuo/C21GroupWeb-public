# Mori--Tanaka 法

这篇文档先给出 Mori--Tanaka 法一般的张量形式；然后在特殊材料对称性和夹杂形状设置下，给出简化公式；最后按照[简介](../index.md)给出的分析角度评估该估计方法。

## 张量形式

Mori-Tanaka 方法认为，远场应变（应力）等于基体相的应变（应力），**基体相的应变集中张量等于单位张量**：

$$
\mathbb{T}^{(1)} = \mathbb{I}
$$

而等效刚度张量表达式为

$$
\boxed{
\mathbb{L}^{c} = \mathbb{L}^{(1)} 
+ \Big( \sum_{\alpha=2}^{N}
c^{(\alpha)} \big( \mathbb{L}^{(\alpha)} - \mathbb{L}^{(1)} \big) : \mathbb{T}^{(\alpha)} \Big)
: \Big(\sum_{\beta=1}^{N} c^{(\beta)} \mathbb{T}^{(\beta)} \Big)^{-1}
}
$$

其中，集中张量 $\mathbb{T}^{(\alpha)}$ **参考的弹性介质为基体相**。类似的，等效柔度张量表达式为

$$
\boxed{
\mathbb{M}^{c} = \mathbb{M}^{(1)} 
+ \Big( \sum_{\alpha=2}^{N}
c^{(\alpha)} \big( \mathbb{M}^{(\alpha)} - \mathbb{M}^{(1)} \big) : \mathbb{W}^{(\alpha)} \Big)
: \Big(\sum_{\beta=1}^{N} c^{(\beta)} \mathbb{W}^{(\beta)} \Big)^{-1}
}
$$

## 特殊情景

### 各向同性

当夹杂形状为球形，且基体为各向同性材料时，此时 Mori-Tanaka 估计公式退化为 Maxwell 近似公式：

$$
\frac{K^{c} - K^{(1)}}{K^{c} + \frac{4}{3}G^{(1)} }
= \sum_{\alpha=1}^{N} c^{(\alpha)} \frac{K^{(\alpha)} - K^{(1)}}{K^{(\alpha)} + \frac{4}{3}G^{(1)} }, \quad
\frac{G^{c} - G^{(1)}}{G^{c} + H^{(1)} } = \sum_{\alpha=1}^{N} c^{(\alpha)}
\frac{G^{(\alpha)} - G^{(1)}}{G^{(\alpha)} + H^{(1)} }
$$

### 双相纤维增强型单胞

$$
\begin{aligned}
k^{c}&=\frac{c_f k_f\left(k_m+m_m\right)+c_m k_m\left(k_f+m_m\right)}{c_f\left(k_m+m_m\right)+c_m\left(k_f+m_m\right)} \\
m^{c}&=\frac{m_f m_m\left(k_m+2 m_m\right)+k_m m_m\left(c_f m_f+c_m m_m\right)}{k_m m_m+\left(k_m+2 m_m\right)\left(c_f m_m+c_m m_f\right)} \\
p^{c}&=\frac{2 c_f p_f p_m+c_m\left(p_f p_m+p_m^2\right)}{2 c_f p_m+c_m\left(p_f+p_m\right)}\\
l^{c}&=\dfrac{c_f l_f\left(k_m+m_m\right)+c_m l_m\left(k_f+m_m\right)}{c_f\left(k_m+m_m\right)+c_m\left(k_f+m_m\right)} \\
n^{c}&=\left(l-c_f l_f-c_m l_m\right) \frac{\left(l_f-l_m\right)}{\left(k_f-k_m\right)}
+ c_f n_f+c_m n_m
\end{aligned}
$$

## 评估

一般认为，M-T 方法适用于（1）夹杂数量大于等于 2，且夹杂的形状和朝向相同，或（2）夹杂数量等于 2，夹杂朝向随机。夹杂的体积分数一般不超过 20% 或 30%。**对于多相且形状或朝向不一致的夹杂，M-T 估计的等效张量可能不具有主对称性**。然而，对于双相材料，但夹杂可能具有不同的形状或朝向，M-T 估计的等效张量仍具有主对称性。

## 更新日志

### 2026/07/24

1. 创建了文档 ` Mori--Tanaka 法.md`
