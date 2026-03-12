# Improved linear representation of surface waves. Part 2. Slowly varying bottoms and currents

By JON WRIGHT<sup>1</sup> AND DENNIS B. CREAMER<sup>2</sup> 

$^{1}$ Institute for Nonlinear Science, University of California, San Diego, La Jolla, CA 92093-0402, USA 

$^{2}$ Naval Research Laboratory, Washington, DC 20375-5000, USA 

(Received 15 January 1993 and in revised form 13 August 1993) 

We extend the results of a previous paper to fluids of finite depth. We consider the Hamiltonian theory of waves on the free surface of an incompressible fluid, and derive the canonical transformation that eliminates the leading order of nonlinearity for finite depth. As in the previous paper we propose using the Lie transformation method since it seems to include a nearly correct implementation of short waves interacting with long waves. We show how to use the Eikonal method for slowly varying currents and/or depths in combination with the nonlinear transformation. We note that nonlinear effects are more important in water of finite depth. We note that a nonlinear action conservation law can be derived. 

# 1. Introduction

Accurate descriptions of the sea-surface shape, while a difficult task, is necessary for many applications. Because the dynamical behaviour of the sea surface is difficult to describe exactly, some approximation is usually required. In a previous paper, Creamer et al. (1989, referred to as Paper I herein), we presented a useful approximation scheme that exactly captures the lowest-order nonlinear behaviour of surface waves, does well at higher orders, and also captures the important features of short waves interacting with longer waves. We considered only irrotational motion and ignored the effects of the wind and of surface tension. The method described in that paper has been used to discuss excitations of capillary waves by the wind (Watson & McBride 1992), and to produce a very realistic-looking simulation of an evolving ocean surface. Zakharov (1991) has applied a similar transformation to treat wind wave growth. 

Because of the apparent usefulness of the technique and the obvious importance of finite depth effects we have extended the method of Paper I to treat surface waves in finite depth. We also show how to solve the problem of slowly varying bottom topography and/or slowly varying currents in the Eikonal limit. For a complete discussion of the nonlinear transformation we refer the reader to Paper I. In order to make this paper readable by itself there will be some duplication of material; in particular, some of the notation will be repeated, but the differences caused by the finite depth will be noted. 

The general idea is to replace the usual functions describing the surface, namely the surface potential and the surface elevation, with two new functions and use linear dynamics in the new functions. We require that the solution of the linearized time-evolution equations in the new variable more nearly describes the true solution than does the corresponding linearized solution in terms of surface elevation and velocity potential. 

The Hamiltonian for the system can be expanded in a power series in the wave-slope field. The quadratic term is responsible for describing the usual linear surface waves. Ordinarily the next term in the Hamiltonian is of cubic order and describes the leading nonlinearity. It is this term that we wish to remove so that the leading nonlinearity is quartic in the Hamiltonian or cubic in the equations of motion. This means that corrections to linear theory are a factor of $(\mathrm{slope})^2$ smaller than the linear terms. 

It can be shown that a canonical transformation exists that will remove the third-order terms if there are no resonant interactions at that order. Phillips (1960) showed that there are no three-wave resonances for surface gravity waves of finite depth, so this leading nonlinear term can be removed. As the product of the depth and a typical surface wave vector approaches zero, the nonlinearities became significantly more important for a fixed surface wave slope. Consequently in very shallow water the transformation will not be useful except for small slopes. However, for typical oceanic conditions depending upon both the sea state and the depth, we expect the transformation to be useful for $KD > \frac{1}{2}$ , where $K$ is a dominant wave vector and $D$ is the depth. 

Usually in shallow water the variation of the depth will be important. Consequently we have included a discussion of how to use the nonlinear transformation in conjunction with the Eikonal equations for slowly varying depths. We also note that slowly varying currents can be treated the same way. 

In the next section we establish the notation and in the third section we present the required transformation. In the fourth section we demonstrate the correct method for implementing the Eikonal equations and give some examples to show the applicability of the method. 

# 2. Surface-wave Hamiltonian

Our starting point is the Hamiltonian for fully nonlinear surface waves. In order to express the surface-wave Hamiltonian in a convenient form, as well as for algebraic convenience, we introduce some notation. We wish to describe waves in terms of the canonical variables defined on the water's surface (Zakharov 1968), i.e. the surface elevation $\zeta(x,t)$ and the velocity potential $\phi(x,z,t)$ evaluated at the surface 

$$
\phi_ {s} (\boldsymbol {x}, t) = \phi [ \boldsymbol {x}, \zeta (\boldsymbol {x}, t), t ], \tag {1}
$$

where $\phi(x, t)$ satisfies Laplace's equation in the interior of the water. We write $x$ for $(x, y)$ and $\partial_x$ for $(\partial_x, \partial_y)$ . Three-dimensional partial derivatives of, for example, the velocity potential occur in the Hamiltonian, and these need to be expressed in terms of the velocity potential (or other functions) at the surface. Let $f(x, t)$ be any function at the surface. An interior function $g(x, z, t)$ is defined by solving Laplace's equation 

$$
\nabla^ {2} g = 0,
$$

$$
g (\boldsymbol {x}, \zeta , t) = f (\boldsymbol {x}, t), \tag {2}
$$

and $\pmb {\hat{n}}\cdot \pmb {\nabla}g = 0$ (3) 

on the bottom and sides (if any). Here $\hat{n}$ denotes the normal to the surface. Thus $f$ is the velocity potential at the surface, $g$ is the velocity potential in the interior. The operators $\mathbf{D}_x = (\mathbf{D}_x, \mathbf{D}_y)$ and $\mathbf{D}_z$ are defined by 

$$
\mathrm {D} _ {j} f = \partial_ {j} g | _ {z = \zeta}. \tag {4}
$$

Thus, $\mathbf{D}_j$ is a linear (but non-local) operator from functions of $(x,y)$ to functions of $(x,y)$ . 

The Hamiltonian for surface waves is equal in value to the energy and is (West 1981, p. 33; Henyey et al. 1988) 

$$
H (\zeta , \phi_ {s}) = \frac {1}{2} \int \mathrm {d} ^ {2} x [ \phi_ {s} \left(\mathrm {D} _ {z} - \left(\partial_ {x} \zeta\right) \cdot \mathbf {D} _ {x}\right) \phi_ {s} + g \zeta^ {2} ], \tag {5}
$$

where $\zeta$ and $\phi_s$ are the canonical variables. We have chosen units in which the density $\rho = 1$ . The operators $\mathbf{D}$ are functionals of $\zeta$ , thereby introducing nonlinearities into $H$ . The complexity of surface-wave calculations is entirely in dealing with the D-operators. One method of working with the $\mathbf{D}$ is to expand in a power series in $\zeta$ . Writing $H = H_2 + H_3 + H_4 + \ldots$ we have, for example, 

$$
H _ {2} = \frac {1}{2} \int \mathrm {d} ^ {2} x [ \phi_ {s} (x, t) \theta (x) \phi_ {s} (x, t) + g \zeta^ {2} (x, t) ], \tag {6}
$$

where we have introduced the operator $\theta(x)$ . This operator and the higher-order terms of the Hamiltonian are most conveniently evaluated in Fourier space, where 

$$
\theta (k) = k \tanh  (k h), \quad k = | k |. \tag {7}
$$

The symbol $h$ is the depth. This operator and the expression for $H_{2}$ are the finite-depth versions of the corresponding infinite-depth expressions in Paper I. They were originally introduced by Miles (1977). The expression (7) used in (6) gives the correct linear dispersion relation for finite-depth surface waves. Introducing the Fourier transform of the canonical variables 

$$
\zeta (x) = \int \frac {\mathrm {d} ^ {2} k}{(2 \pi) ^ {2}} \mathrm {e} ^ {\mathrm {i} k \cdot x} \zeta_ {k}, \tag {8}
$$

and $\phi_s(x) = \int \frac{\mathrm{d}^2k}{(2\pi)^2}\mathrm{e}^{\mathrm{i}k\cdot x}\phi_k,$ (9) 

we have 

$$
H _ {3} = \frac {1}{4} (2 \pi) ^ {2} \iiint \frac {\mathrm {d} ^ {2} k _ {1}}{(2 \pi) ^ {2}} \frac {\mathrm {d} ^ {2} k _ {2}}{(2 \pi) ^ {2}} \frac {\mathrm {d} ^ {2} k _ {3}}{(2 \pi) ^ {2}} \delta^ {2} \left(k _ {1} + k _ {2} + k _ {3}\right) \zeta_ {1} \phi_ {2} \phi_ {3} \left(k _ {2} ^ {2} + k _ {3} ^ {2} - k _ {1} ^ {2} - 2 \theta_ {2} \theta_ {3}\right). \tag {10}
$$

We have adopted the notation $\zeta_1 \equiv \zeta_{k_1}$ and similarly for $\phi$ . Again this is the same expression as in Miles (1977). The difference between the above expression for $H_3$ and that in Paper I is that the last term in (10) would just be $k_2 k_3$ . 

The goal is now to find a transformation of variables for which the new $H_{3} = 0$ . 

# 3. The Lie transform of $\zeta$ and $\phi_s$

The description in terms of the Fourier-transformed variables, $\zeta_{k}$ and $\phi_{k}$ turns out to be most convenient. We next introduce two interpolating functions $Z(k,\lambda)$ and $\varPhi(k,\lambda)$ where $\lambda$ varies between zero and one. (The time dependence is temporarily suppressed.) At $\lambda = 0$ , $Z(k,0) = \zeta(k)$ and $\varPhi(k,0) = \phi_s(k)$ . 

The fields $Z(\pmb{k}, \lambda = 1)$ and $\Phi(\pmb{k}, \lambda = 1)$ are to be the new fields. We now introduce a functional $W$ that generates the changes in $Z$ and $\Phi$ as a function of $\lambda$ : 

$$
\begin{array}{l} W = (2 \pi) ^ {2} \int \int \int \frac {\mathrm {d} ^ {2} k _ {1}}{(2 \pi) ^ {2}} \frac {\mathrm {d} ^ {2} k _ {2}}{(2 \pi) ^ {2}} \frac {\mathrm {d} ^ {2} k _ {3}}{(2 \pi) ^ {2}} \delta^ {2} (\pmb {k} _ {1} + \pmb {k} _ {2} + \pmb {k} _ {3}) \\ \times [ D (1, 2, 3) \Phi_ {1} \Phi_ {2} \Phi_ {3} + B (1, 2, 3) Z _ {1} Z _ {2} \phi_ {3} ]. \quad (1 1) \\ \end{array}
$$

Here the notation $\varPhi_{1}$ denotes $\varPhi(k_1,\lambda)$ . The two functions $D(1,2,3) = D(k_1,k_2,k_3)$ and $B$ are to be determined so that the new $H_{3}$ is zero. In Paper I we showed that the new Hamiltonian $K$ , can be written in terms of the old Hamiltonian, $H = H_{2} + H_{3} + H_{4} + \ldots$ as 

$$
K (Z, \Phi) = H _ {2} (Z, \Phi) + H _ {3} (Z, \Phi) - \{H _ {2} (Z, \Phi), W (Z, \Phi) \} + \dots , \tag {12}
$$

where $\{, \}$ is the usual Poisson bracket (see (19)). The idea is now to choose $W$ so that the last two terms cancel. Thus if we ignore terms of fourth order in the slope, $K$ is purely quadratic and the associated equations of motion are linear. It is straightforward to verify that 

$$
\begin{array}{l} \{H _ {2}, W \} = (2 \pi) ^ {2} \iiint \frac {\mathrm {d} ^ {2} k _ {1}}{(2 \pi) ^ {2}} \frac {\mathrm {d} ^ {2} k _ {2}}{(2 \pi) ^ {2}} \frac {\mathrm {d} ^ {2} k _ {3}}{(2 \pi) ^ {2}} \delta^ {2} (\pmb {k} _ {1} + \pmb {k} _ {2} + \pmb {k} _ {3}) \\ \times \{g \zeta_ {3} [ 3 D (1, 2, 3) \phi_ {1} \phi_ {2} + B (1, 2, 3) \zeta_ {1} \zeta_ {2} ] - [ \theta_ {2} B (1, 2, 3) + \theta_ {3} B (3, 1, 2) ] \zeta_ {1} \phi_ {2} \phi_ {3} \}. \quad (1 3) \\ \end{array}
$$

Comparing (10) and (13) we see that if we choose 

$$
D = + \frac {\left[ k _ {1} ^ {2} \left(\theta_ {1} ^ {2} - \left(\theta_ {2} - \theta_ {3}\right) ^ {2}\right) + k _ {2} ^ {2} \left(\theta_ {2} ^ {2} - \left(\theta_ {1} - \theta_ {3}\right) ^ {2}\right) + k _ {3} ^ {2} \left(\theta_ {3} ^ {2} - \left(\theta_ {1} - \theta_ {2}\right) ^ {2}\right) - 2 \theta_ {1} \theta_ {2} \theta_ {3} \left(\theta_ {1} + \theta_ {2} + \theta_ {3}\right) \right]}{(1 2 R g)} \tag {14}
$$

and $B = \frac{[\theta_3(\theta_3(\theta_1 + \theta_2) - \theta_1^2 - \theta_2^2) + \theta_3(k_1^2 + k_2^2 - 2k_3^2) + (\theta_1 - \theta_2)(k_1^2 - k_2^2)]}{(2R)}$ (15) 

with $R = \theta_{1}(\theta_{2} + \theta_{3} - \theta_{1}) + \theta_{2}(\theta_{1} + \theta_{3} - \theta_{2}) + \theta_{3}(\theta_{1} + \theta_{2} - \theta_{3})$ 

then the third-order piece of the new Hamiltonian, $K_{3}$ , is zero. These expressions are similar to those in Paper I, except that some $\theta$ have been replaced by their infinite-depth values. Using the notation 

$$
\bar {\varphi} _ {s} = \Phi (\lambda = 1), \quad \bar {\zeta} = Z (\lambda = 1) \tag {16}
$$

for the new canonical variables we have 

$$
K \left(\bar {\phi} _ {s}, \bar {\zeta}\right) = \frac {1}{2} \int \mathrm {d} ^ {2} x \times \left[ \bar {\phi} _ {s} (x, t) \theta (x) \bar {\phi} _ {s} (x, t) + g \bar {\zeta} ^ {2} (x, t) \right] + O \left(\bar {\zeta} ^ {4}\right). \tag {17}
$$

It remains to implement the transformation between the old and new variables. This identical to what is done in Paper I, but will be repeated here for readability; more detail can be found in the original paper. The functional, $W$ , generates the transformation. Let $A(\lambda)$ be either $\varPhi$ or $Z$ . Then 

$$
\frac {\partial \boldsymbol {A}}{\partial \lambda} = \{A, W \}, \tag {18}
$$

where the Poisson bracket is defined by 

$$
\{A, B \} = \int \mathrm {d} ^ {2} x \left[ \frac {\delta A}{\delta \zeta (x)} \frac {\delta B}{\delta \phi_ {s} (x)} - \frac {\delta A}{\delta \phi_ {s} (x)} \frac {\delta B}{\delta \zeta (x)} \right]. \tag {19}
$$

The boundary conditions are 

$$
\Phi (\lambda = 0) = \phi_ {s}, \quad Z (\lambda = 0) = \zeta . \tag {20}
$$

The equations to be solved are 

$$
\begin{array}{l} \frac {\partial Z _ {1} (\lambda)}{\partial \lambda} = (2 \pi) ^ {2} \iint \frac {\mathrm {d} ^ {2} k _ {2}}{(2 \pi) ^ {2}} \frac {\mathrm {d} ^ {2} k _ {3}}{(2 \pi) ^ {2}} \delta^ {2} (k _ {1} - k _ {2} - k _ {3}) \\ \times [ 3 D (1, 2, 3) \Phi_ {2} (\lambda) \Phi_ {3} (\lambda) + B (2, 3, 1) Z _ {2} (\lambda) Z _ {3} (\lambda) ], \tag {21} \\ \end{array}
$$

$$
\frac {\partial \Phi_ {1} (\lambda)}{\partial \lambda} = - 2 (2 \pi) ^ {2} \iint \frac {\mathrm {d} ^ {2} k _ {2}}{(2 \pi) ^ {2}} \frac {\mathrm {d} ^ {2} k _ {3}}{(2 \pi) ^ {2}} \delta^ {2} \left(k _ {1} - k _ {2} - k _ {3}\right) B (1, 2, 3) Z _ {2} (\lambda) \Phi_ {3} (\lambda). \tag {22}
$$

In our experience, integrating the transformation equations is much easier than integrating the equations of motion, and it is only necessary to perform the transformations at times for which a description of the wave field is desired. 

In order to obtain the velocities it is necessary to calculate $\partial \zeta (x,t) / \partial t$ . This is most easily obtained by differentiating (21) with respect to time and integrating it as a function of $\lambda$ simultaneously with (21) and (22). 

# 4. Eikonal equations

We now consider the possibility that there are slowly varying currents present or that the depth is slowly varying. We wish to demonstrate that the proper action transport equations are given in terms of the new variables. For simplicity we consider only the case of a horizontal current, $\pmb{u}$ . The velocity in the presence of this current is decomposed as usual into 

$$
\boldsymbol {v} = \boldsymbol {u} + \nabla \phi . \tag {23}
$$

An additional term appears in the Hamiltonian: 

$$
H _ {c} = - \int \mathrm {d} \boldsymbol {x} \phi_ {s} (\boldsymbol {x}) \boldsymbol {u} \cdot \nabla \zeta . \tag {24}
$$

Other new terms in the Hamiltonian involve the derivative of the slowly varying current $\pmb{u}$ and are therefore ignored. In order to simplify the discussion, we suppose that the fluid can be approximately divided into three regions. In the first and third regions the bottom and the current do not vary. In the middle region they may both change slowly. 

Now imagine a wave packet that is essentially localized in the first region at an initial time and propagates through region two and into region three so that at sometime it is entirely within region three. If we were solving this problem in the linear wave approximation we would consider the usual ray tracing and action conservation equations as discussed in Whitham (1974). 

We first give a prescription for solving this propagation problem and then we will justify the prescription. We will discuss only the finite-depth problem, but the variable-current solution follows essentially the same steps. 

Step one: transform the initial wave packet using the Lie transform appropriate for the initial depth assuming constant initial depth. 

Step two: propagate from the initial depth to the final depth using the Eikonal equations for a linear system, but using as initial conditions the transformed variables found in step one. 

Step three: transform the final wave packet back to the physical surface variables, now using the Lie transform appropriate for the final depth, assuming constant final depth. 

The argument supporting steps one and three relies heavily on the fact that the Lie 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/44c79985-9f7b-4458-9eab-7f53f03bf0d3/5f5f9fe462e3dceed768c62b1bc4ff7bb544fd2ff4a868d6b56912db489b9716.jpg)



FIGURE 1. Surface waves for a depth of $1\mathrm{m}$ and a wavelength of $2\pi$ . The level of the troughs has been adjusted to better aid the comparison. The dot-dash curve is a linear $\zeta(x) = 0.2\cos x$ , the dashed curve is a Stokes wave, and the solid curve is our nonlinear approximation.


transform of a wave packet remains where it started. It changes its shape, but not its position. No part of the disturbance moves to a new location. This was demonstrated in Paper I. That means that it only knows about the local bottom, which implies that if the local bottom is flat, or nearly so, we can use the Lie transform for a flat bottom appropriate to the local depth and the new wave packet will be the same as if we transformed with the exact transformation. To make this latter point clearer let us review the procedure to find the nonlinear transformation. 

We start with the exact $H$ as given by (5). We then imagine expanding in powers of the surface height. This could be done in principle numerically, but fortunately we do not have to actually do it. Then we imagine finding $W$ in (12) so that the $H_{3}$ term is cancelled. The higher-order terms are neglected and the remaining $H_{2}$ is given by 

$$
H _ {2} = \frac {1}{2} \int d ^ {2} x \left[ \phi_ {s} (x, t) \frac {d}{d z} \phi (x, t) | _ {z = 0} + g \left(\zeta^ {2} (x, t) \right. \right]. \tag {25}
$$

This is identical in form to the usual quadratic Hamiltonian, except that $\phi_s, \zeta$ are the transformed variables. Thus to solve the equations of motion we need only solve linear equations. Now for a slowly varying bottom this problem can be treated via Eikonal methods, Whitham (1974). The important point is that the Lie transform - even the exact form - does not change the form of $H_2$ . It only changes the variables $\phi_s$ and $\zeta$ . Thus the Lie transform for a flat bottom appropriate to the local depth will give the same $H_2$ , and because the transform is local in the sense that it does not move part of the wave packet elsewhere, the exact transform and the locally flat transform will give the same new variables $\phi_s$ and $\zeta$ . 

Thus the prescription is to transform to new variables, propagate via the usual ray equations and action transport equations, and then transform back using (11) for the appropriate depth. Nonlinear effects will not show up in the ray equations, but they will show up in the action equations since it is the action in the new variables that is most nearly conserved. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/44c79985-9f7b-4458-9eab-7f53f03bf0d3/085c7ed1579e5519896cd08dc6d2747738a505b8d390a140a858ef0f9afd5c63.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/44c79985-9f7b-4458-9eab-7f53f03bf0d3/efaa317c7ddc18356eecb3c8b7f89bb785d7562d9e11d1ead1467949bb1a1776.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/44c79985-9f7b-4458-9eab-7f53f03bf0d3/df80b318118c3f4c8788482b2ad5ed7d8b163f373c537846450ad21bbfbac04e.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/44c79985-9f7b-4458-9eab-7f53f03bf0d3/c907559d3a23186f9e03f1024e65f5689686e760ca85233831b5ddea9b4a7d7e.jpg)



FIGURE 2. (a) The solid curve is the initial deep-water wave-height spectrum. The dashed curve is the spectrum after propagation into water $10\mathrm{m}$ deep. Although a wide range of Fourier components were used in the nonlinear transformation, all Fourier components of the wave height for wavenumber greater than 0.2 have been set to zero in order to make it easier to see the effects of propagation. (b) A typical deep-water sea surface generated from the spectrum in (a). The solid line can either be interpreted as the surface in linear approximation or as the transformed surface variables. The dashed line is the actual sea surface predicted by the nonlinear transformation. (c) The surface waves in (b) were transported via the Eikonal equations into water $20\mathrm{m}$ deep. Note that the transformation has raised and narrowed the peaks. (d) The surface waves in (b) were transported via Eikonal equations into water $10\mathrm{m}$ deep.


In order to illustrate the method and to show some of the boundaries of applicability, we have constructed several examples. A sample calculation showing the nonlinear surface which results from a simple sine wave in the linear approximation is shown in figure 1. The input is a linear wave $\zeta(x) = 0.2\cos x$ . The top curve is the linear surface wave and the bottom curve is the Stokes wave. The middle curve is our approximation. The depth is 0.7 and the strength has been adjusted so that all waves have about the same peak-to-peak wave height. This particular situation is very nonlinear. 

In this example with a steep wave, in shallow water, the next-order nonlinear terms 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/44c79985-9f7b-4458-9eab-7f53f03bf0d3/e99da4877509824158e5d022eb687fedbf0ca498b7c23b131e3fe63ef5f9ed1e.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/44c79985-9f7b-4458-9eab-7f53f03bf0d3/318cc821a7dbaec21968b0c2cbc491ab1da585370dae12b3d57d35695fa0e333.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/44c79985-9f7b-4458-9eab-7f53f03bf0d3/4b8d89652f47e0b6f98ab2268349fe0a1fa7267c22996d2cc4b1700f10b357ef.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/44c79985-9f7b-4458-9eab-7f53f03bf0d3/60ac14da6ba876740ae51765141bfd95244a402dbb5e62939b96a3eb43cc90c3.jpg)



FIGURE 3. The conditions for this figure are the same as figure 2 except that the spectrum has been multiplied by a factor of two. At $20\mathrm{m}$ depth the nonlinear transformation has a modest effect as shown in (c), but at $10\mathrm{m}$ depth the transformation has a very large effect and care should be taken in deciding whether it is still meaningful.


are beginning to be important as can be seen by comparing the middle and bottom curves. Nevertheless it is a definite improvement over the linear approximation. 

We now consider the evaluation of a spectrum of waves travelling from deep water to shallow water. Here we will be able to illustrate the calculations for a realistic wave field propagating over a slowly changing bottom. The first spectrum that we used is shown in figure 2(a). The specific form of the one-dimensional wave height spectrum was chosen to be 

$$
\Phi (K) = \operatorname {c o n s t}. \frac {K ^ {2}}{\left(K ^ {2} + K _ {0} ^ {2}\right) ^ {2 . 5}}. \tag {26}
$$

The solid line is the spectrum in deep water, and the dotted line is the spectrum obtained by solving the Eikonal equations for transport across a slowly decreasing depth ending in water $10\mathrm{m}$ deep. Figures 2(b), 2(c) and 2(d) all show a realization of the surface height after propagating to a final depth. The solid curve is the new surface 

height variable and as such does not represent a physical surface. The dotted curve is the actual water surface height. For purposes of plotting, all Fourier components for the wavenumber greater than $0.2\mathrm{rad~m}^{-1}$ have been set equal to zero. Note that in deep water, figure 2(b), the nonlinear transformation does not have much effect. Thus one can interpret all the solid curves as the result of doing linear dynamical calculations. In figure 2(c) the final depth is $20\mathrm{m}$ and we can see that there has been a significant narrowing of the peaks and an increase in height for the higher crests. In figure 2(d) an example is shown at a depth of $10\mathrm{m}$ . Here we begin to see significant nonlinear behaviour. 

In order to see the effect of a more vigorous wave field, we increased the spectrum by a factor of 2, as shown in figure 3. Now there is an increase in nonlinear behaviour - as can be observed in figure $3(d)$ . This suggests that this particular combination of wave height spectrum and depth is nearly too nonlinear to be treated by these methods. 

# 5. Conclusions

In this paper, we have shown how to incorporate the leading-order nonlinearities of shallow-water, irrotational surface waves in a calculational scheme that (i) preserves the canonical, Hamiltonian structure of the theory, (ii) provides a resulting linear theory of surface dynamics, and (iii) is computationally efficient. We showed how to incorporate this scheme into an Eikonal calculation involving a sloping bottom. The inclusion of nonlinearities into the Eikonal method is usually non-trivial (Whitham 1974 devoted a good part of his book to this subject), but we have shown how the leading-order nonlinearities can be trivially incorporated. 

The figures clearly show that typical oceanic conditions can be treated by this method, but that care must be taken if the water becomes too shallow or the height spectrum becomes too large. Because of the complicated dependence on both the energy spectrum and the depth, it is not easy to be precise about the range of validity. As the depth approaches zero, the effective nonlinearity increases and at some depth the fourth-order terms in the transformed Hamiltonian will be comparable to the second-order terms. Other means than the method proposed here are necessary for that situation. Since strong nonlinearities are important in describing solitary shallow-water waves, the proposed method would be of doubtful utility for describing such waves. 

The calculational method developed in this paper has other areas of use. In a recent paper, Ding & Farmer (1993) calculated statistics on breaking wave events by simulations of the sea surface using only linear variables. Given the inadequacy of the linear stochastic model of the ocean surface (visually it does not look at all realistic), Ding & Farmer suggest that something like our method be incorporated into the Monte-Carlo studies. Looking at breaking wave statistics in shallow water would require the extension developed in this paper. 

This work was supported by the Office of Naval Research on the grant N40014-90-J-4135. 

# REFERENCES



CREAMER, D. B., HENYEY, F., SCHULT, R. & WRIGHT, J. 1989 Improved linear representation of ocean surface waves. J. Fluid Mech. 205, 135-161 (referred to herein as Paper I). 





DING, L. & FARMER, D. M. 1993 A Monte-Carlo study on breaking wave statistics and comparison with field observations. Preprint, Institute of Ocean Sciences, Sidney, B.C., Canada. 





HENYEY, F. S., CREAMER, D. B., DYSTHE, K. B., SCHULT, R. L. & WRIGHT, J. A. 1988 The energy and action of small waves riding on large waves. J. Fluid Mech. 189, 443-462. 





MILES, J. W. 1977 On Hamilton's principle for surface waves. J. Fluid Mech. 83, 153-158. 





PHILLIPS, O. M. 1960 On the dynamics of unsteady gravity waves of finite amplitude. Part 1. J. Fluid Mech. 9, 193-217. 





WATSON, K. M. & McBride, J. 1993 Excitation of capillary waves by longer waves. J. Fluid Mech. 250, 103-119. 





WEST, B. J. 1981 Deep Water Gravity Waves, p. 33. Springer. 





WHITHAM, G. B. 1974 Linear and Nonlinear Waves. Wiley. 





ZAKHAROV, V. E. 1968 Stability of periodic waves of finite amplitude on the surface of a deep fluid. J. Appl. Mech. Tech. Phys. 9, 190-194. 





ZAKHAROV, V. E. 1991 Inverse and direct cascade in the wind-driven surface turbulence and wave breaking. Proc. IUTAM Congress on Wave Breaking, Sydney, Australia. 

