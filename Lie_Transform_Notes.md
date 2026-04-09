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

更具体地说，如果 `A(lambda)` 和 `B(lambda)` 都满足

```text
partial_lambda A = {A,W}
partial_lambda B = {B,W}
```

那么它们的 Poisson bracket 也会满足

```text
partial_lambda {A,B} = {{A,B},W}
```

这说明：

- bracket 本身也随着同一个 flow 演化
- Poisson 结构不会被破坏
- 因而 canonical 关系也会被保留下来

所以真正的原因不是“Lie transform 这个名字本身很神秘”，而是：

**它被定义成由 Poisson bracket 生成的 Hamiltonian flow，而 Hamiltonian flow 天生保 Poisson / canonical 结构。**

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

## 10.5 为什么 `(3.1)` 只保留 `\Phi\Phi\Phi` 和 `ZZ\Phi`

这里一个容易误解的点是：

- paper 不是说别的三次组合在数学上绝对写不出来
- 而是说在当前目标下，它们被主动排除了

paper 在 `(3.1)` 后面明确说了两件事：

1. `(3.1)` **不是**最一般的 cubic form
2. 但它是**对 `\phi_s` 奇**的最一般 cubic form

之所以只保留

- `\Phi\Phi\Phi`
- `ZZ\Phi`

是因为作者想要一个：

- near-identity
- 三次精度
- 足以消去原有 `H_3`
- 且不会给 transformed Hamiltonian 引入原 Hamiltonian 里没有的额外结构

在这个标准下：

- `Z^3`
- `Z\Phi^2`

并不是“绝对不能写”，而是会被排除：

- 它们对 `\phi_s` 是偶的
- 并且会向 transformed Hamiltonian 里引入作者不想保留的额外项型

所以 `(3.1)` 应该理解成：

**一个为了当前三次消项任务而选定的最小充分 ansatz，而不是所有可能生成泛函的唯一形式。**

## 10.6 为什么 Hamiltonian 对 `\phi_s` 的奇偶性重要

这点之所以重要，是因为原问题里 Hamiltonian 对 `\phi_s` 是偶的。

从物理上看，这对应的是：

- 动能对速度是二次的
- 势能不依赖速度符号
- 所以把速度势的符号翻转 `\phi_s -> -\phi_s`，Hamiltonian 的数值不变

在表面波问题里，这个性质也可以理解成“时间反演型”的一个影子：

- `\eta` 不变
- `\phi_s` 变号
- 能量不应该因此改号

所以如果你一开始选的生成泛函包含不合适的 `\phi_s` 奇偶结构，
那变换后的 Hamiltonian 就可能产生原问题里没有的 `\phi_s` 奇项。

作者在 `(3.1)` 这里的思路是：

- 既然原 Hamiltonian 在 `\phi_s` 上是偶的
- 那就尽量选一个生成泛函，使 transformed Hamiltonian 继续保持这个性质

这并不是说“奇偶性比别的都神圣”，而是说：

- 这是原系统本来的结构约束
- 保留这个约束可以防止我们人为引入不必要的项型
- 也让变换后的系统更像原系统，而不是被多余自由度污染

所以这点最适合记成一句话：

**关心 `\phi_s` 的奇偶性，不是因为作者在做形式游戏，而是因为原始自由面 Hamiltonian 本来就在 `\phi_s` 上是偶的；如果生成泛函选得不合适，就可能人为制造出原系统没有的奇项结构。**

## 10.7 为什么到了 Lie transform 还是直接用了生成函数的 `B,D`

这点 paper 其实说得很直接，但第一次读时很容易觉得跳步。

关键不是：

- “Lie transform 和 global generating function 完全一样”

而是：

- **在当前只做到三次阶、只关心消去 `H_3` 的层面上，它们可以做到等价。**

paper 在 `(3.10)` 后面明确说：

- 从 Lie flow 得到的新旧变量关系
- 在“第一个非平凡阶”上
- 如果把 `W` 选成 global generating functional `(3.1)` 里三次部分的负号

那么它和 `(3.2)-(3.3)` 得到的关系是一样的。

这就意味着：

- 为了三次阶消项，前面通过 global generating function 求出来的 `B,D`
- 可以直接拿来作为 Lie transform 生成元 `W` 里的 kernels

也就是说，作者的策略是：

1. 先用 global generating function 比较直接地解出 `B,D`
2. 再转到 Lie transform 语言
3. 用同样的 `B,D` 写 `W`
4. 由于两种方法在当前阶数上等价，所以这已经足够保证三次项消失

## 10.8 这样做的理由是什么

主要有三个理由。

### 理由 1：省去重复求解

既然 `B,D` 已经在 global generating function 里通过“消去 cubic term”求出来了，
而 Lie transform 在当前阶数上又和它等价，那就没有必要再重新从零求一遍。

### 理由 2：当前目标只到三次阶

paper 在这里的核心目标是：

- 消去 `H_3`
- 得到一个到这一阶已经足够好的 canonical transform

而不是重新建立一整套完全独立的高阶 Lie kernel 理论。

所以在当前精度下，直接沿用 `B,D` 是最自然的做法。

### 理由 3：真正“特别”的不是 `B,D`，而是后续怎样使用 Lie flow

paper 后面也明确说：

- 如果只看前三阶，两种方法是等价的
- 但更高阶时它们会开始不同
- 经验上 Lie transform 产生的结果更好，尤其在 1D 几何解释、短波骑长波、谱修正这些问题上

所以：

- `B,D` 在起步阶段可以共用
- 真正拉开差距的是高阶结构已经被怎样“预先吸收到变换里”

这也是为什么作者会说 Lie transform seems special, though they had not fully understood why.

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

## 12. 什么是 Jacobi 恒等式

Jacobi 恒等式是 Poisson bracket 最重要的结构恒等式之一：

```text
{A,{B,C}} + {B,{C,A}} + {C,{A,B}} = 0
```

它不是附加假设，而是“一个 bracket 真正像 Poisson bracket 那样工作”必须满足的核心性质。

在这里，它的作用非常具体。

我们已经有

```text
partial_lambda A = {A,W}
partial_lambda B = {B,W}
```

于是

```text
partial_lambda {A,B}
= {{A,W},B} + {A,{B,W}}
```

这时候用 Jacobi 恒等式，就可以把右边改写成

```text
{{A,B},W}
```

所以 Jacobi 恒等式在这里告诉我们：

**“先取 bracket 再沿 flow 演化”和“先各自沿 flow 演化再取 bracket”是兼容的。**

这正是“Lie flow 保 Poisson 结构”的数学核心。

## 13. 下一步最自然的问题

这份笔记之后，最自然要继续讲的是：

- 为什么 `partial_lambda A = {A,W}`` 自动就是 canonical flow
- 为什么 Hamiltonian 自己也会按同样规则变换
- 为什么到三次阶时只需要盯住 `H_3 - {H_2,W}`

## 14. 关于 1D 简化，需要单独记住的一点

进入 section 4 以后，一个最容易被误会的地方是：

- 好像需要先人为指定 `k_1,k_2,k_3` 的大小关系

其实不是。

真正起作用的是：

- 问题退到 `1D`
- 同时 triad 满足

```text
k_1 + k_2 + k_3 = 0
```

在一维里，这已经自动推出：

- 必有一个波数与另外两个符号相反
- 因而那个波数的绝对值等于另外两个绝对值之和

所以 section 3 的一般二维 kernels 一旦代到 1D triad 上，就会自动塌缩：

- `(3.4)` 的分子必有一因子为零，所以 `D=0`
- `(3.5)` 会塌成只含 `|k|` 和 `sgn(k)` 的简单形式

所以这里真正支配 simplification 的，是一维 triad 闭合本身，而不是某个额外指定的大小排序。

## 15. 为什么换成 `chi` 以后，就回到了坐标变换图像

这是 section 4 最值得抓住的转折点。

前面 section 3 做的事情，还是比较抽象：

- 存在一个 canonical transform
- 它能把 `H_3` 吸收到变量变换里

但到 section 4，1D 的 Lie-flow 方程已经被化成 characteristic form。

这时引入 `chi` 的意义就是：

- 不再直接把解看成固定欧拉坐标 `x` 上的 PDE
- 而是把每条 characteristic 用它在初始时刻的位置来编号

所以 `chi` 本质上就是：

- characteristic label
- 初始参数坐标
- 被搬动之前的“点的编号”

而方程

```text
x = chi + lambda ZTilde_0(chi)
```

说的就是：

- 原来标号为 `chi` 的点
- 被水平搬到了当前位置 `x`

所以到了这里，Lie transform 已经不再只是抽象的 Hamiltonian bookkeeping。
它在 1D 里被真正解成了一个具体的坐标重排。

## 16. 为什么 `(4.13)` 可以看成最终的物理表面重建

有了 `chi` 之后，section 4 的逻辑就完整了：

1. 先在 transformed variables 里做简单的线性描述
2. 用 characteristic relation 解出点标签 `chi` 与欧拉位置 `x` 的关系
3. 再把这些线性变量按这个映射重投影回物理表面

所以 `(4.13)` 的本质不是“又写了一个复杂积分公式”，而是：

**把线性变量通过已经解出的 1D 坐标重排，重建成真实物理表面。**

也就是说：

- section 3 告诉你有一个 canonical transform
- section 4 告诉你在 1D 里这个 transform 实际上怎么搬动点
- `(4.13)` 则是把这种点重排重新写回物理表面的最终公式

这也是为什么 `(4.13)` 会自然地产生：

- crest sharpening
- trough flattening
- bound harmonics

因为这些结构不是额外自由波，而是重排后的 Eulerian surface 几何。

## 17. “重建”不是再求动力学，而是坐标翻译

这里值得单独记一句：

**section 4 的 reconstruction 不是在求另一个新的时间演化问题，而是在把参数坐标里的点集重新读成欧拉表面。**

更具体地说：

- `chi` / `y` 是表面点的标签
- `x` 是这些点经过 Lie-flow 重排后出现的欧拉位置
- characteristic relation 给出的是“同一个点从哪里搬到哪里”

所以 `(4.13)` 的作用是：

- 已知点的标签表示
- 已知点如何被水平搬动
- 把它们翻译成固定空间位置 `x` 上看到的真实表面

这就是为什么 reconstruction 会自然产生 bound harmonics：

- 它来自点集重排后的 Eulerian 几何
- 而不是来自新加入的独立自由波

## 18. 关于 section 4 的 `y`

读 section 4 时，一个很容易混淆的点是：`y` 在这里不是二维空间里的“纵向方向轴”。

在 1D Lie-transform / reconstruction 语境里，`y` 更像是：

- 参数坐标
- characteristic label 的一种写法
- 与 `chi` 同类的积分变量

所以当公式从 `chi` 换成 `y` 时，通常不是物理上换了一个新方向，而只是：

- 仍然在描述同一批表面点
- 只是把点标签的记号写成了 `y`

这也是为什么 paper 的 `(4.14)` 写成 `dy`，而很多 practical rewrite 会让人误以为那是欧拉坐标积分。
更稳的理解应该是：

**`dy` 代表对参数点标签积分，`x` 才是最后重建出来的欧拉位置。**
