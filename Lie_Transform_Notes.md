# 李变换笔记

这份笔记是接在生成函数之后的第二份桥接材料。

目标不是一上来就推完整个 Creamer 第 3 节，而是先回答几个最核心的问题：

- 李变换到底是什么
- 为什么它和正则变换有关
- 为什么它特别适合做 `H_3` 的消去
- 它和生成函数是什么关系

## 1. 最短直觉

如果说生成函数像是：

- 一次性写下整张坐标变换说明书

那么李变换更像是：

- 不直接写终点
- 而是写一个“无穷小正则流”
- 然后让变量沿这个流从旧坐标慢慢走到新坐标

所以李变换的关键词不是“整体公式”，而是：

- 生成元
- 流
- 逐阶展开

## 2. 为什么叫“变换”，又为什么叫“流”

在 Lie transform 的语言里，我们不直接写

```text
new variable = some explicit function of old variable
```

而是引入一个辅助参数 `lambda`，让任意量 `A` 随着 `lambda` 演化：

```text
partial_lambda A = {A, W}
```

这里：

- `W` 是生成元
- `{A,W}` 是由 Poisson bracket 产生的无穷小变化
- `lambda=0` 对应原变量
- `lambda=1` 对应变换后的变量

所以它看起来像一条“流线”：

- 从旧变量出发
- 沿着由 `W` 生成的 canonical flow 走到新变量

## 3. 为什么这会保持 canonical 结构

因为 `partial_lambda A = {A,W}`` 这件事本身，就是哈密顿流的形式。

而哈密顿流的一个核心性质就是：

- 它保辛结构
- 因而也保 canonical structure

所以李变换不是“任意流”，而是：

**由某个生成元 `W` 产生的正则流。**

这就是为什么它可以拿来做 canonical transform。

## 4. 为什么它比生成函数更适合 Creamer 第 3 节

Creamer 的目标不是单纯证明“存在某个正则变换”，而是：

- 逐阶看 Hamiltonian 怎么变
- 把三次项 `H_3` 搬走
- 最后还要保留一个便于理解和后续 1D 几何解释的形式

Lie transform 在这方面特别顺手，因为：

- 它天然适合逐阶展开
- Poisson bracket 的阶数计数很清楚
- 很容易直接写出

```text
K = H_2 + H_3 - {H_2,W} + ...
```

然后要求

```text
{H_2,W} = H_3
```

## 5. 它和生成函数是什么关系

当前可以先把两者理解成：

- 生成函数：整体、隐式地定义 canonical transform
- Lie transform：把同一个 canonical transform 改写成由生成元驱动的无穷小流

所以它们不是两套不同目标，而是：

**同一个正则变换问题的两种表达方式。**

## 6. 为什么 Lie transform 对“消项”特别自然

因为 Hamiltonian 在 Lie transform 下也跟着流动。

当只保留到三次阶时，最重要的结构就是：

```text
K = H_2 + H_3 - {H_2,W} + O(4)
```

所以只要我们选 `W` 满足

```text
{H_2,W} = H_3
```

就得到

```text
K = H_2 + O(4)
```

这就是“把 `H_3` 吸收到变换里”。

## 7. 为什么 `W` 要是三次

因为：

- `H_2` 是二次
- 我们希望 `{H_2,W}` 去抵消三次的 `H_3`

而 Poisson bracket 的阶数规则是

```text
deg{F,G} = deg(F) + deg(G) - 2
```

所以必须有

```text
2 + deg(W) - 2 = 3
```

也就是

```text
deg(W) = 3
```

## 8. 这一部分目前最该记住的句子

如果只记一句，就记这个：

**Lie transform 就是用生成元 `W` 定义一个无穷小 canonical flow，并沿着这个 flow 把变量从旧坐标推进到新坐标；对 Creamer 来说，它最重要的作用是让我们能逐阶安排 `\{H_2,W\}=H_3`，从而把 `H_3` 从 Hamiltonian 里搬到变量变换里。**

## 9. 频率分母到底从哪里来

这部分是读 Creamer `(3.4)-(3.5)` 时最容易觉得突兀的地方。

关键点是：

- `H_2` 代表线性波动
- 线性波动的每个 Fourier 模态都有自己的频率 `omega_k`
- `W` 是一个三次泛函，所以它在 Fourier 空间里会由“三个模态的乘积乘上某个待定 kernel”组成

当我们去算

```text
{H_2,W}
```

时，本质上就是让线性 Hamiltonian `H_2` 去作用在 `W` 的三模乘积上。

而 `H_2` 的作用就是给每个模态带上线性振荡频率，所以对一个三模乘积，最后会拉出这样的组合：

```text
+/- omega_1 +/- omega_2 +/- omega_3
```

这和普通时间依赖里

```text
e^{-i omega_1 t} e^{-i omega_2 t} e^{+i omega_3 t}
```

对时间求导会拉出

```text
-i (omega_1 + omega_2 - omega_3)
```

是同一种结构。

所以在 Fourier 空间里，`{H_2,W}`` 的每一项都可以理解成：

```text
(frequency combination) x (kernel of W) x (mode product)
```

接下来我们要求

```text
{H_2,W} = H_3
```

于是就要把 `W` 的 kernel 解出来。这个时候就会出现

```text
kernel of W ~ kernel of H_3 / (frequency combination)
```

所以这些 `+/- omega_1 +/- omega_2 +/- omega_3` 最开始是由 `H_2`
作用在 `W` 上时产生的“系数”，而在反过来解 `W` 的时候，它们就变成了分母。

## 10. 这为什么和 resonance 有关

因为一旦某个 triad 满足

```text
+/- omega_1 +/- omega_2 +/- omega_3 = 0
```

那这个分母就没法除。

这就说明：

- 那部分三次项不是可以被平滑吸收到坐标变换里的非共振项
- 而是真正的共振项

所以可以把判断逻辑记成：

- `delta^2(k_1+k_2+k_3)` 负责波数闭合
- `+/- omega_1 +/- omega_2 +/- omega_3` 负责频率闭合
- 两者同时满足才是真正 resonance

深水重力波的关键事实就是：

- triad 可以因为 `delta` 而出现在 `H_3` 里
- 但由于 `omega = sqrt(g |k|)`，对应频率组合不会真的为零
- 所以 `H_3` 存在，但仍然是非共振的

## 11. Lie transform 不神秘，它只是一个工具

到这一步，最值得明确写下来的理解是：

**Lie transform 不是一种新的物理，也不是额外神秘的对象。它只是一个特别适合做逐阶 canonical 消项的工具。**

对 Creamer 来说，它的作用非常朴素：

- 把整体 canonical transform 改写成无穷小 canonical flow
- 让 Hamiltonian 的变化可以按 Poisson bracket 一阶一阶展开
- 于是可以清楚地安排

```text
{H_2,W} = H_3
```

来消掉三次项

所以如果只保留最实用的理解，可以把 Lie transform 想成：

- 一个把“找正则变换”改写成“解逐阶方程”的工具

## 12. 下一步最自然的问题

这份笔记之后，最自然要继续讲的是：

- 为什么 `partial_lambda A = {A,W}`` 自动就是 canonical flow
- 为什么 Hamiltonian 自己也会按同样规则变换
- 为什么到三次阶时只需要盯住 `H_3 - {H_2,W}`
