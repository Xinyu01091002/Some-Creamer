# 1D Creamer Transform 实用笔记

这份笔记只回答一个问题：

**如果不再从 section 2/3 的抽象理论出发，而只想知道 1D 深水里“实用的 Creamer transform 到底怎么用”，应该怎么理解？**

结论先说：

- 1D 的实用 Creamer form 不是一套新的动力学方程
- 它更像一个三步算法：
  1. 给定线性变量
  2. 计算与 Hilbert transform 相关的水平重排
  3. 用这个重排把线性变量重建成物理表面

所以它真正可用的形式，是 section 4 的：

- characteristic remapping
- reconstruction formulas `(4.13)` / `(4.14)`

而不是前面抽象的 `B,D`、`W` 或 Poisson bracket 本身。

## 1. 先分清输入和输出

在 1D 深水 practical use 里，可以把事情分成两类对象：

### 输入：线性变量

这是 transformed / linear variables，一般记成：

- `tilde zeta_0(y)`
- `tilde phi_{s0}(y)`

这里的 `y` 是参数坐标，也可以看成 surface-point label。

这些量满足的是截断后的线性 Hamiltonian 动力学。

### 输出：物理变量

这是你最终真正想要的表面：

- `zeta(x)`
- `phi_s(x)`

这里的 `x` 是欧拉坐标，也就是固定空间位置。

所以 practical Creamer transform 的任务可以压成一句话：

**把参数坐标 `y` 上的线性变量，翻译成欧拉坐标 `x` 上的物理表面。**

## 2. 为什么 1D 会出现 Hilbert transform

在 1D 里，section 3 的一般 kernels 会塌缩成：

- `D = 0`
- `B` 只剩下 `|k|` 和 `sgn(k)` 的简单结构

而 Fourier 空间里：

```text
Hilbert transform  <->  multiplication by  -i sgn(k)
```

所以 1D kernel 一旦只剩 `sgn(k)` 结构，最自然的实空间对象就是 Hilbert transform。

于是作者定义

```text
ZTilde = H[Z]
PhiTilde = H[Phi]
```

这里 `H` 表示 Hilbert transform。

在单频波的线性近似里：

- `cos(kx)` 的 Hilbert transform 是 `sin(kx)`

所以 `ZTilde` 正好可以理解成和表面高度相差 90 度相位的量，也因此有“水平位移”的自然解释。

## 3. characteristic 形式才是 1D practical form 的真正入口

在 `ZTilde` 变量下，1D Lie-flow 方程会变成：

```text
partial_lambda ZTilde   = -(1/2) partial_x (ZTilde^2)
partial_lambda PhiTilde = - ZTilde partial_x PhiTilde
```

第一条可以写成

```text
partial_lambda ZTilde + ZTilde partial_x ZTilde = 0
```

这是标准的一阶 characteristic / transport 方程。

这意味着：

- `ZTilde` 沿 characteristic 保持不变
- characteristic label 可以记成 `chi`

积分以后得到坐标映射：

```text
x = chi + lambda ZTilde_0(chi)
```

这就是实用 1D Creamer form 的第一个核心结果：

**Lie transform 在 1D 里被真正解成了一个水平坐标重排。**

## 4. `chi` 和 `x` 分别是什么

这里一定要分清：

- `chi`：参数坐标，表面点的标签
- `x`：欧拉坐标，点被搬动后出现的空间位置

所以关系

```text
x = chi + lambda ZTilde_0(chi)
```

的意思是：

- 原来编号为 `chi` 的点
- 被水平搬到了 `x`

你可以把这一点看成 1D Creamer transform 的核心几何内容。

## 5. “重建”到底在重建什么

这里最容易误会的是：

- 好像作者又重新求了一个新的表面解

其实不是。

section 4 的 reconstruction 做的是：

**把参数坐标里的同一批表面点，重新按固定欧拉位置来读成物理表面。**

所以 reconstruction 不是：

- 再求一个新的动力学
- 再加入新的 free-wave 模态

而是：

- 已知点标签 `y` / `chi`
- 已知这些点被搬到哪里
- 把它们翻译成 `zeta(x)` 和 `phi_s(x)`

## 6. 实用版算法：1D Creamer transform 怎么用

如果只讲实用层面，可以把它写成一个 3-step algorithm。

### Step 1. 给定线性输入

例如单频输入：

```text
tilde zeta_0(y)   = a sin(k y)
tilde phi_s0(y)   = -(a/k) cos(k y)
```

或者更一般的窄带波包、线性随机海面。

### Step 2. 解水平重排

通过 characteristic relation，给每个参数点 `y`（或 `chi`）计算它对应的欧拉位置：

```text
x = y - tilde zeta_0(y)
```

这里符号写成 `+` 还是 `-` 取决于你用的是 forward 还是 backward form；
paper 在 `(4.13)` 阶段采用的是 backward reconstruction 形式。

这一步给出的是：

- 参数点到欧拉点的映射

### Step 3. 重建物理表面

然后用 `(4.13)` / `(4.14)` 把这些参数点上的线性变量，重新读成：

```text
zeta(x), phi_s(x)
```

这个输出就是 physical surface。

## 7. 为什么这会自动产生 bound harmonics

因为坐标映射本身是非线性的。

在 Fourier 空间里，`(4.14)` 的关键相位会变成

```text
exp[-ik (y - tilde zeta_0(y))]
 = exp(-iky) exp[i k tilde zeta_0(y)]
```

然后展开：

```text
exp[i k tilde zeta_0(y)]
 = 1 + i k tilde zeta_0(y) + (i k tilde zeta_0(y))^2 / 2 + ...
```

这就自动生成：

- 二阶
- 三阶
- 更高阶

的 bound harmonic 结构。

所以：

- 谐波不是因为多加了新的独立自由波
- 而是因为坐标重排把原来简单的相位变成了非线性相位

## 8. 为什么它看起来像 Stokes 波形

因为：

- 线性输入在参数坐标里仍然很简单
- 但一经过水平重排再读回欧拉坐标
- crest 区域会被压缩
- trough 区域会被拉宽

于是自然出现：

- crest sharpening
- trough flattening
- higher harmonics

这就是为什么 section 4 的 1D Lie-transformed waveform 会很接近 Stokes-like 非线性波形。

## 9. 这和“free wave”为什么不同

因为整个过程没有改动线性变量的自由传播结构。

它做的是：

- 在线性变量里保持简单传播
- 在坐标映射和欧拉重建里加入非线性几何

所以额外结构主要应理解为：

- reconstruction-generated bound structure

而不是：

- 新产生的一套独立 free-wave dynamics

## 10. 这一份笔记最该记住的句子

如果只记一句，就记这个：

**1D 深水里 practical Creamer transform 的本质，就是先把线性变量写在参数坐标上，再用 Hilbert-transform 驱动的水平重排 `x(y)`，把同一批表面点重新读成物理欧拉表面 `zeta(x)`。**

## 11. 如果只想实际使用，可以怎么想

以后如果你不是在推理论，而只是想“怎么用 1D Creamer”，就按下面这张心智卡片来：

- 输入：`tilde zeta_0(y), tilde phi_s0(y)`
- 核心几何：`x = y - tilde zeta_0(y)` 这样的 characteristic remapping
- 输出：`zeta(x), phi_s(x)` via `(4.13)` / `(4.14)`
- 非线性来源：不是新动力学，而是欧拉重建
- 谐波性质：主要是 bound harmonics

## 12. 关于 section 4 里的 `y`，一个很容易忘记的提醒

这里的 `y` 不是二维水平坐标里的“另一个方向轴”。

在 1D section 4 里，它本质上和 `chi` 扮演同一类角色：

- 它是参数坐标
- 它是表面点的标签
- 它是重建积分里的积分变量

也就是说，`y` 更接近“Hilbert-transform / characteristic 表示里的点编号”，而不是
`(x,y)` 那种二维平面中的第二个空间方向。

所以当 paper 写出类似 `(4.14)` 的 `dy` 积分时，最好先在脑子里翻译成：

- “我在对参数点集积分”
- 而不是“我在沿几何上的 `y` 方向积分”
