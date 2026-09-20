# Gâteaux 导数 与 Fréchet 导数

这篇文章首先给出 Gâteaux 导数 与 Fréchet 导数的定义和符号，然后提供三个简单算例：（1）一维空间中的导数（2）二维空间的方向导数，并给出 Fréchet 导数的图像解释（3）无穷维空间中泛函的变分。

## Gâteaux 导数 与 Fréchet 导数的定义

Gâteaux 导数的定义为

$$
D_{g} F( y;\eta ) \triangleq \lim_{\epsilon\to0} \frac{F( y+\epsilon \eta ) - F( y)}{\epsilon}
= \left.\frac{\mathrm{d} }{\mathrm{d} \epsilon} F( y+\epsilon \eta )\right|_{\epsilon=0}.
$$

由定义式可看到，如果将给定点处函数 $F$ 的 Gâteaux 导数看作是关于 $\eta$ 的新函数，那么 $D_{g}F$ 与 $F$ 属于同一个函数空间。当涉及到定义在无穷维函数空间中的泛函时，需要严格界定泛函定义域所属的函数空间，以及测试函数 $\eta$ 所在的函数空间。此外，根据定义式，Gâteaux 导数是关于函数 $\eta$ 的**一阶齐次函数**，也即

$$
D_{g}F(y;\alpha \eta) = \alpha D_{g}F(y;\eta), \quad \forall \alpha \in \mathbb{R}.
$$

继续推广 Gâteaux 导数的定义到高阶情景，有

$$
D_{g}^{n} F( y;\eta )
= \left.\frac{\mathrm{d}^n }{\mathrm{d} \epsilon^n} F( y+\epsilon \eta )\right|_{\epsilon=0}.
$$

Fréchet 导数是在点 $y$ 处寻找一个**有界线性映射** $D_{F}: X\mapsto Y$，其中 $X$ 是被求导函数定义域所属的空间，$Y$ 是被求导函数值域所属的空间。若关于点 $y$，对任意的 $h\in X$，存在有界线性映射 $D_{F}(y)$，使得

$$
\lim_{\|h\|\to 0} \frac{ \| F(y+h) - F(y) -D_{F}(y)\cdot h \| }{\|h\|} = 0,
$$

那么函数在该点处 Fréchet 可微。

## 一些算例

以下给出的一些算例，是在 Gâteaux 导数存在的假设下按照定义计算得到的。

1. **一维** 函数

$$
y(t) = t^2,
$$

关于 $\mathrm{d}t$ 的 Gâteaux 导数为

$$
D_{g} y(t;\,\mathrm{d} t) = \lim_{\epsilon\to 0} \frac{(t+\epsilon \,\mathrm{d} t)^2 - t^2}{\epsilon} = 2t\,\mathrm{d} t.
$$

同时，常数 $2t$ 自然地定义了点 $t$ 处的有界线性映射，因此也是该函数的 Fréchet 导数。

2. **二维** 函数

$$
y(t, s) = \begin{cases}
\dfrac{t^2 s}{t^4 + s^2}, & \text{if } (t, s) \ne (0, 0), \\
0, & \text{if } (t, s) = (0, 0).
\end{cases}
$$

在点 $(0,0)$ 处，对方向 $\boldsymbol{n}=(n_1,n_2)$ 的 Gâteaux 导数为

$$
D_{g} y(0; \boldsymbol{n}) = \lim_{\epsilon\to0} \frac{n_1^2 n_2}{\epsilon^2 n_1^4 + n_2^2}
= \begin{cases}
n_1^2/n_2, & \text{if } n_2 \neq 0, \\
0, & \text{else}.
\end{cases}
$$

虽然该函数在 $(0,0)$ 处方向导数**处处存在**，但是不存在 Fréchet 导数。假如存在 Fréchet 导数，那么二维空间中的有界线性映射总可以表示为一个二维向量，可以通过两条不同的参数曲线确定；选取方向 $\boldsymbol{n}=(1,0)$ 和 $(0,1)$，就可以确定与有界线性映射等价的向量为 $(0,0)$，对应图中 $z=0$ 平面。但如果选取参数曲线（图中蓝色曲线） $x=t$，$y=t$，在点 $(0,0)$ 处并不与平面 $z=0$  相切。

<figure style="text-align: center;">
  <img 
    src="https://raw.githubusercontent.com/qiuxiaohuo/picgo-bed/main/figsmaxmin.png" 
    style="zoom:50%;" 
  />
  <figcaption>Gâteaux 导数存在，但是 Fréchet 导数不存在的例子</figcaption>
</figure>

3. （来自李少凡 [Introduction to Micromechanics and Nanomechanics](https://www.worldscientific.com/worldscibooks/10.1142/8644)）定义如下**泛函**

$$
\Pi(y) = \int_{0}^{1} p(t) y'^2 + q(t) y^2 + 2f(t)y \,\mathrm{d} t,
$$

泛函 $\Pi( y+\epsilon \eta )$ 可以写作

$$
\begin{aligned}
\Pi( y+\epsilon \eta ) &= \int_{0}^{1} p(t) (y' + \epsilon\eta')^2 + q(t) (y+\epsilon\eta)^2 + 2f(t)(y+\epsilon\eta) \,\mathrm{d} t \\
&= \left\{ \int_{0}^{1} py'^2+qy^2+2fy \,\mathrm{d} t \right\} 
+ 2\epsilon \left\{ \int_{0}^{1} py'\eta'+qy\eta+f\eta \,\mathrm{d} t \right\} 
+ \epsilon^2 \left\{ \int_{0}^{1} p\eta'^2+q\eta^2 \,\mathrm{d} t \right\},
\end{aligned}
$$

所以 Gâteaux 导数等于

$$
D_{g} \Pi( y;\eta ) = 2\int_{0}^{1} \left( py'\eta'+qy\eta+f\eta \right) \,\mathrm{d} t.
$$

而如果泛函定义在 $H^1(0,1)$ 空间中，$p$，$q$ 充分光滑，那么该泛函也存在 Fréchet 导数，并恰好等于 Gâteaux 导数。注意该泛函的 Gâteaux 导数是关于 $\eta$ 的有界线性泛函。

