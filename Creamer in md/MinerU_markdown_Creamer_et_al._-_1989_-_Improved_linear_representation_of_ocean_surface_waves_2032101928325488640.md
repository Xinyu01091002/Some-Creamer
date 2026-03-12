# Improved linear representation of ocean surface waves

# By DENNIS B. CREAMER, FRANK HENYEY, ROY SCHULT AND JON WRIGHT

Center for Studies of Nonlinear Dynamics, La Jolla Institute, 7855 Fay Avenue, Suite 320, La Jolla, CA 92037, USA 

(Received 16 February 1988 and in revised form 30 January 1989) 

We apply the idea of choosing new variables that are nonlinear functions of the old in order to simplify calculations of irrotational, surface gravity waves. The usual variables consist of the surface elevation and the surface potential, and the transformation to the new variables is a canonical (in Hamilton's sense) one so as to maintain the Hamiltonian structure of the theory. We further consider the approximation of linear dynamics in these new variables. This approximation scheme exactly reproduces the effects of the lowest-order nonlinearities in the usual variables, does well at higher orders, and also captures important features of short waves interacting with longer waves. We provide a physical interpretation of this transformation which is correct in the one-dimensional case, and approximately so in the two-dimensional case. 

# 1. Introduction

The appropriate characterization of the sea surface, while a difficult task, is necessary for a myriad of applications. The dynamical behaviour of the sea surface will never be described exactly, thus requiring some sort of approximation. In this paper we present a useful approximation scheme that exactly captures the lowest-order nonlinear behaviour of surface waves, does well at higher orders, and also captures the important features of short waves interacting with longer waves. We consider only irrotational motion and ignore the effects of the wind and of surface tension. 

The general idea is to replace the usual functions describing the surface, namely the surface potential and the surface elevation, with two new functions and use linear dynamics in the new functions. We require that the solution of the linearized time-evolution equations in the new variables more nearly describes the true solution than does the corresponding linearized solution in terms of surface elevation and velocity potential. In order to determine the new functions we combine a transformation of variables with a perturbation expansion in the wave slope (by variables we mean the functions referred to above). We show that the leading nonlinear terms in the equations of motion (which are of quadratic order) can be entirely removed by this transformation of variables, and hence in the new variables the dynamics is linear plus terms of cubic order in the wave slope. Thus, dynamical effects which are of quadratic order in the wave slope are then contained in this kinematic transformation of variables. If we are concerned with properties of the ocean surface for which it is a good approximation to ignore cubic and higher order nonlinearities, the proposed transformation of variables becomes very useful since it is not then necessary to solve 

partial differential equations for the time evolution of the surface waves. The dynamics of the new variables is given by a linear equation which is identical to the usual linear equations. In addition, one form of the transformation contains a very good approximation to the important dynamics of short waves interacting with long waves and hence is far superior to a perturbation expansion in the slope. 

There is, or course, some cost - the nonlinear transformation of variables must be inverted. The complexity of this nonlinear transformation must be compared to the complexity of solving the dynamical equations. For most applications, one is only interested in the sea surface at discrete times. If the time spacing is not too small, it is more efficient to solve the nonlinear transformation at each time step. For many applications one is not even interested in time evolution as such, but only in its consequences for statistical properties of waves. 

We have not considered the problem of higher-order perturbation theory. For calculations for which effects in higher-order perturbation theory are important, there may be no particular benefit to our change of variables. However, since some of the short-wave long-wave interaction is captured correctly the new variables could be better for higher-order calculations. 

In this paper we do not develop any possible applications of the transformation of variables. However, at this point we shall briefly mention four uses for which the formalism could be useful. There are undoubtedly others. 

(i) Solution of the inverse problem for ocean waves. Suppose one wishes to determine the surface wave field of the ocean from a number of measurements at discrete times. One possibility is to integrate the partial differential equations for the time evolution directly. An alternative is to solve the linear time evolution equations for the new variables and then solve the transformation equations at the discrete times. The length of time between the measurements will determine which method is better by comparing the computational requirements for integrating the evolution equations for one time interval versus the requirements for solving the transformation equations. 

(ii) Short waves riding on long waves. We have discovered that one of the possible transformation of variables leads to a simple, analytic representation of short waves on long waves. It does almost as well as integrating the eikonal equations directly and is much simpler. 

(iii) Statistics of surface gravity waves. Ocean surface waves are a weakly nonlinear system and the statistics of such waves show nearly Gaussian behaviour. In the limit that the nonlinearity goes to zero the statistics would be Gaussian. Since the new variables are more nearly linear (the leading nonlinear behaviour is absent), we expect their statistics to be more nearly Gaussian than the usual variables. The leading deviation from Gaussianity of the usual physical variables should be contained in the transformation. 

(iv) Assess nonlinear wave codes. Wave modelling (Hasselmann 1984) usually includes only the effect of the four-wave resonant interactions to study wave growth. Since our transformation does not affect the lowest-order resonant interactions and removes the effect of three-wave interactions, we view these modelling efforts as describing the dynamics of the transformed variables. Application of wave-modelling results to the real ocean requires the transformation back to the usual physical variables, which, we shall show, does not significantly change the spectrum. Use of the form of the four-wave term in the new variables might give superior results. 

We use a Hamiltonian description of surface waves. Within the framework of 

Hamiltonian mechanics there has evolved a theory of canonical transformations. A description can be found in many text books, such as Goldstein (1980). If the Hamiltonian can be expanded in a perturbation theory in powers of a parameter, and if the leading nonlinearity contains no resonant interactions, there is a theorem that says that there exists a change of variables that will give a new canonical Hamiltonian system for which the leading-order nonlinearity is absent. It is well known that for surface gravity waves there are no three-wave resonant interactions, and so it seems worthwhile to find such transformations. 

The transformation is not unique, and the question naturally arises as to which is in some sense best, and yet still simple. We explicitly exhibit two transformations and show that linear dynamics in one of them captures many important features of surface-wave dynamics beyond that given by the lowest-order nonlinearities. Linearizing in these new variables leads to exactly the same equations of motion as one gets for usual variables; the frequencies are given by the dispersion relation, $\omega = (g k)^{\frac{1}{2}}$ , where $g$ is the gravitational constant and $k$ is the wavenumber. The shape of the surface and the value of the velocity at any point on the surface will, however, be different. 

In the next section we define notation and describe the Hamiltonian for surface gravity waves. Section 3 describes two different canonical transformations from the physical variables to the new variables: one using a global generating function technique and another using Lie transforms. The reader who is concerned only with applications in one dimension can skip to §4 and consider (4.13), relating the physical and new variables, as given. Sections 4, 5, and 6 show various one-dimensional results (using Lie transforms) which can be summarized as follows: 

(i) There exists a practical way to evaluate the Lie transform in terms of the horizontal Lagrangian displacement. 

(ii) The Lie transform provides a good approximation to a Stokes wave of moderate wave slope. 

(iii) In the case of a short-wave packet on a long wave, the transformation correctly incorporates the modulation of the group velocity of the wave packet by the long wave. The position, wavelength, and amplitude of the short wave are correctly given (at least to leading order in the long-wave slope). Because this transformation is time independent, this result occurs as a kinematic effect while it is usually thought of as a dynamic effect due to the advection of the short waves by the long wave. 

(iv) The modification of the surface-height spectrum due to the transformation is proportionately small over the whole range of wavenumbers. 

(v) There is the proper balance of radiation stresses so that there is no change in the average surface elevation. The Stokes drift and the attendant deep-water return current combine properly to ensure this result. 

These results are presented using Lie-transform theory. If we had used some other canonical transformation (such as that from the global generating function in §3) the results would be different. In the case of the global generating function, the approximation to a Stokes wave is still rather good but results 3 and 4 change completely. The modulation of the group velocity of a short-wave packet is incorrectly given, while the modification to the spectrum would be enormous at high wavenumbers, with results similar to the perturbative calculations of Barrick & Weber (1976). Thus the transformed variables in the Lie-transform scheme seem to be special though we have yet to discover why. In §4 we discuss a physical 

interpretation for the one-dimensional case. Finally, in §7 we present some results on applying our formalism to two-dimensional waves. The physical interpretation of the one-dimensional case is extended to two dimensions and the behaviour of short waves on long waves is discussed. In §8 we summarize our results. 

# 2. Surface-wave Hamiltonian

Our starting point is the Hamiltonian for fully nonlinear surface waves. This Hamiltonian, or the related variational principle, has been independently rediscovered by a number of workers. (For a review, see Miles 1981.) Normally, a Legendre transformation carries one from a variational principle to a Hamiltonian. In the case of fluid problems, however, the constructed variational principle is in canonical form, and it is trivial to read off the Hamiltonian. By canonical form for a Lagrangian (Courant-Hilbert 1937) is meant one depending not only on $q$ , but also on $p$ in the following way: 

$$
L \left(p _ {j}, q _ {j}; \dot {p} _ {j}, \dot {q} _ {j}\right) = \sum_ {j} p _ {j} \dot {q} _ {j} - H \left(p _ {j}, q _ {j}\right). \tag {2.1}
$$

Thus, from a canonical form Lagrangian it is trivial to read off the $p, q$ pairs of canonical variables and the Hamiltonian $H(p, q)$ . 

In order to express the surface-wave Hamiltonian in a convenient form, as well as for algebraic convenience, we introduce some notation. We wish to describe waves in terms of the canonical variables defined on the water's surface (Zakharov 1968), i.e. the surface elevation $\zeta (\pmb {x},t)$ and the velocity potential $\phi (\pmb {x},\pmb {z},t)$ evaluated at the surface 

$$
\phi_ {\mathrm {s}} (\boldsymbol {x}, t) = \phi [ \boldsymbol {x}, \zeta (\boldsymbol {x}, t), t ], \tag {2.2}
$$

where $\phi (\pmb {x},z,t)$ satisfies Laplace's equation in the interior of the water. We write $\pmb{x}$ for $(x,y)$ and $\partial_x$ for $(\partial_{x},\partial_{y})$ . Three-dimensional partial derivatives of, for example, the velocity potential occur in the Hamiltonian, and these need to be expressed in terms of the velocity potential (or other functions) at the surface. Let $f(x,t)$ be any function at the surface. An interior function $g(\pmb {x},z,t)$ is defined by solving Laplace's equation 

$$
\nabla^ {2} g = 0,
$$

$$
g (\boldsymbol {x}, \zeta , t) = f (\boldsymbol {x}, t), \tag {2.3}
$$

and $\pmb {\hat{n}}\cdot \pmb {\nabla}g = 0$ (2.4) 

on the bottom and sides (if any). Here $\hat{\pmb{n}}$ denotes the normal to the surface. Thus if $f$ is the velocity potential at the surface, $g$ is the velocity potential in the interior. The operators $\mathbf{D}_x = (\mathbf{D}_x, \mathbf{D}_y)$ and $\mathbf{D}_z$ are defined by 

$$
\mathrm {D} _ {j} f = \partial_ {j} g | _ {z - \zeta}. \tag {2.5}
$$

Thus, $\mathbf{D}_j$ is a linear (but non-local) operator from functions of $(x,y)$ to functions of $(x,y)$ . $\mathbf{D}_x$ differs from $\partial_x$ , since 

$$
\begin{array}{l} \partial_ {x} f (x, t) = \partial_ {x} g [ x, \zeta (x, t) t ] \\ = \mathbf {D} _ {x} f + \left(\partial_ {x} \zeta\right) \mathbf {D} _ {z} f. \tag {2.6} \\ \end{array}
$$

Relations obeyed by the D are 

$$
\mathbf {D} _ {x} + \left(\partial_ {x} \zeta\right) \mathrm {D} _ {z} = \partial_ {x}, \tag {2.7}
$$

and $\mathbf{D}_x\cdot \mathbf{D}_x + \mathbf{D}_z^2 = 0.$ (2.8) 

The Hamiltonian for surface waves is equal in value to the energy and is (West 1981; Henyey et al. 1988) 

$$
H (\zeta , \phi_ {\mathrm {s}}) = \frac {1}{2} \int \mathrm {d} ^ {2} x [ \phi_ {\mathrm {s}} (\mathbf {D} _ {z} - (\partial_ {x} \zeta) \cdot \mathbf {D} _ {x}) \phi_ {\mathrm {s}} + g \zeta^ {2} ], \tag {2.9}
$$

where $\zeta$ and $\phi_{\mathrm{s}}$ are the canonical variables. We have chosen units in which the density $\rho = 1$ . The operators $\mathbf{D}$ are functionals of $\zeta$ , thereby introducing nonlinearities into $H$ . The complexity of surface-wave calculations is entirely in dealing with the $\mathbf{D}$ operators. One method of working with the $\mathbf{D}$ is to expand in a power series in $\zeta$ . Writing $H = H_{2} + H_{3} + H_{4} + \ldots$ we have, for example, 

$$
H _ {2} = \frac {1}{2} \int \mathrm {d} ^ {2} x [ \phi_ {\mathrm {s}} (\pmb {x}, t) \theta (x) \phi_ {\mathrm {s}} (\pmb {x}, t) + g \zeta^ {2} (\pmb {x}, t) ], \tag {2.10}
$$

where we have introduced the operator 

$$
\theta (x) \equiv [ - (\partial_ {x}) ^ {2} ] ^ {\frac {1}{2}}. \tag {2.11}
$$

This operator and the higher-order terms of the Hamiltonian are most conveniently evaluated in Fourier space, where (2.11) becomes just the magnitude of the wave vector, $\theta = |\pmb{k}|$ . Introducing the Fourier transform of the canonical variables 

$$
\zeta (x) = \int \frac {\mathrm {d} ^ {2} k}{(2 \pi) ^ {2}} \mathrm {e} ^ {\mathrm {i} k \cdot x} \zeta_ {k}, \tag {2.12a}
$$

and 

$$
\phi_ {\mathrm {s}} (x) = \int \frac {\mathrm {d} ^ {2} k}{(2 \pi) ^ {2}} \mathrm {e} ^ {\mathrm {i} k \cdot x} \phi_ {k}, \tag {2.12b}
$$

we have 

$$
H _ {3} = \frac {1}{4} (2 \pi) ^ {2} \iiint \frac {\mathrm {d} ^ {2} k _ {1}}{(2 \pi) ^ {2}} \frac {\mathrm {d} ^ {2} k _ {2}}{(2 \pi) ^ {2}} \frac {\mathrm {d} ^ {2} k _ {3}}{(2 \pi) ^ {2}} \delta^ {2} \left(\boldsymbol {k} _ {1} + \boldsymbol {k} _ {2} + \boldsymbol {k} _ {3}\right) \zeta_ {1} \phi_ {2} \phi_ {3} \left(\theta_ {2} ^ {2} + \theta_ {3} ^ {2} - \theta_ {1} ^ {2} - 2 \theta_ {2} \theta_ {3}\right). \tag {2.13}
$$

We have adopted the notation $\zeta_1\equiv \zeta_k$ , and similarly for the $\phi$ 

For surface waves one can, of course, treat the higher-order terms perturbatively. However, one has to be careful about the time development of the wave field in any treatment involving a perturbative parameter, denoted by $\lambda$ . Naive perturbation theory gives results that have secular growth in the form $\lambda t$ , where $t$ is the time, and thus are valid only for small time ( $t \ll 1 / \lambda$ ) and cannot be expected to yield the proper, long-time behaviour of the wave field. Standard texts on classical mechanics suggest developing perturbation theories based on canonical transformations in which the time evolution is properly taken into account. There is a plethora of techniques for accomplishing this transformation, of which we shall consider two: global generating functions and Lie transforms. The common requirement for these different methods is that, after the transformation to the new variables, there is no third-order piece in the Hamiltonian. The fourth-order pieces of the new Hamiltonian will be different for different transformation schemes. However, at this order, the resonant interaction will be the same, while the non-resonant interactions will be different for the various transformations. In fifth and higher orders both resonant and non-resonant interactions are different. Actually it is the difference in how much of the interactions (in higher orders) is already incorporated into the transformation that is responsible for determining which transformation gives the better approximation to exact answers. Our empirical result is that Lie transforms enjoy a privileged position. 

# 3. Lie-transform theory

As mentioned earlier, we shall concentrate on a particular canonical transform known as a Lie transform. In order to elucidate our discussion we first construct a global generating function. We want to perform a transformation on the canonical variables that appear in the Hamiltonian (2.9), $\zeta$ (surface height) and $\phi_{\mathrm{s}}$ (potential on the surface), which preserves the canonical structure and eliminates the cubic terms appearing in the Hamiltonian. This transformation takes $\zeta \rightarrow \bar{\zeta}$ and $\phi_{\mathrm{s}}\rightarrow \bar{\phi}_{\mathrm{s}}$ and can be accomplished by the use of a global, time-independent generating functional $F$ , which depends on one old variable and the other new variable, 

$$
\begin{array}{l} F \left(\phi_ {\mathrm {s}}, \bar {\zeta}\right) = \int \mathrm {d} ^ {2} x \phi_ {\mathrm {s}} (x) \bar {\zeta} (x) - \iiint \frac {\mathrm {d} ^ {2} k _ {1}}{(2 \pi) ^ {2}} \frac {\mathrm {d} ^ {2} k _ {2}}{(2 \pi) ^ {2}} \frac {\mathrm {d} ^ {2} k _ {3}}{(2 \pi) ^ {2}} \\ \times (2 \pi) ^ {2} \delta^ {2} \left(\boldsymbol {k} _ {1} + \boldsymbol {k} _ {2} + \boldsymbol {k} _ {3}\right) [ D (1, 2, 3) \phi_ {1} \phi_ {2} \phi_ {3} + B (1, 2, 3) \bar {\zeta} _ {1} \bar {\zeta} _ {2} \phi_ {3} ], \tag {3.1} \\ \end{array}
$$

where $D(1,2,3)$ and $B(1,2,3)$ are functions which will be chosen so as to eliminate the cubic terms. Notice that (3.1) does not have the most general cubic form. Other possible combinations would add terms to the transformed Hamiltonian which are not in the original Hamiltonian and so their coefficients are set to zero. Note, however, that (3.1) is the most general cubic form odd in $\phi_{\mathrm{s}}$ ; thus the transformed $H$ will remain even in $\phi_{\mathrm{s}}$ . The transformed variables can be determined via 

$$
\begin{array}{l} \zeta (x) = \frac {\delta}{\delta \phi_ {\mathrm {s}} (x)} F (\phi_ {\mathrm {s}}, \bar {\zeta}) \\ = \bar {\zeta} (x) - (2 \pi) ^ {2} \int \int \int \frac {\mathrm {d} ^ {2} k _ {1}}{(2 \pi) ^ {2}} \frac {\mathrm {d} ^ {2} k _ {2}}{(2 \pi) ^ {2}} \frac {\mathrm {d} ^ {2} k _ {3}}{(2 \pi) ^ {2}} \mathrm {e} ^ {- \mathrm {i} k _ {1} \cdot x} \delta^ {2} (\pmb {k} _ {1} + \pmb {k} _ {2} + \pmb {k} _ {3}) \\ \times \left[ 3 D (1, 2, 3) \phi_ {2} \phi_ {3} + B (2, 3, 1) \bar {\zeta} _ {2} \bar {\zeta} _ {3} \right], \tag {3.2} \\ \end{array}
$$

and 

$$
\begin{array}{l} \bar {\phi} _ {s} (x) = \frac {\delta}{\delta \bar {\zeta} (x)} F (\phi_ {s}, \bar {\zeta}) \\ = \phi_ {\mathrm {s}} (x) - 2 (2 \pi) ^ {2} \iiint \frac {\mathrm {d} ^ {2} k _ {1}}{(2 \pi) ^ {2}} \frac {\mathrm {d} ^ {2} k _ {2}}{(2 \pi) ^ {2}} \frac {\mathrm {d} ^ {2} k _ {3}}{(2 \pi) ^ {2}} \mathrm {e} ^ {- \mathrm {i} k _ {1} \cdot x} \delta^ {2} \left(\boldsymbol {k} _ {1} + \boldsymbol {k} _ {2} + \boldsymbol {k} _ {3}\right) B (1, 2, 3) \bar {\zeta} _ {2} \phi_ {3}. \tag {3.3} \\ \end{array}
$$

Notice that, since the generating functional is a mixture of the old and new variables, (3.2) and (3.3) are a pair of implicit equations for the new variables (in terms of the old) or vice versa. Equations (3.2-3.3) are non-trivial to solve. The iterations (3.2) and (3.3) can be used to determine $\zeta$ and $\phi_{\mathrm{s}}$ as expansions in $\bar{\zeta}$ and $\bar{\phi}_{\mathrm{s}}$ which are inserted into the expressions for the Hamiltonian (2.10) and (2.13). After some tedious algebra the new Hamiltonian can be shown to have no cubic terms in the new variables provided that 

$$
D (1, 2, 3) = \frac {\left(\theta_ {1} + \theta_ {2} + \theta_ {3}\right) \left(\theta_ {1} + \theta_ {2} - \theta_ {3}\right) \left(\theta_ {1} - \theta_ {2} - \theta_ {3}\right) \left(\theta_ {1} + \theta_ {3} - \theta_ {2}\right)}{1 2 g \left[ \theta_ {1} \left(\theta_ {2} + \theta_ {3} - \theta_ {1}\right) + \theta_ {2} \left(\theta_ {1} + \theta_ {3} - \theta_ {2}\right) + \theta_ {3} \left(\theta_ {1} + \theta_ {2} - \theta_ {3}\right) \right]}, \tag {3.4}
$$

$$
\text {a n d} \quad B (1, 2, 3) = \frac {\left(\theta_ {1} - \theta_ {2}\right) ^ {2} \left(\theta_ {1} + \theta_ {2}\right) - \theta_ {3} ^ {2} \left(2 \theta_ {3} - \theta_ {1} - \theta_ {2}\right)}{2 \left[ \theta_ {1} \left(\theta_ {2} + \theta_ {3} - \theta_ {1}\right) + \theta_ {2} \left(\theta_ {1} + \theta_ {3} - \theta_ {2}\right) + \theta_ {3} \left(\theta_ {1} + \theta_ {2} - \theta_ {3}\right) \right]}, \tag {3.5}
$$

where the quantity $\theta_{i} = |\pmb{k}_{i}|$ is defined in wavenumber space. These new variables can 

be used to investigate the connection between the physical variables and these variables. The denominator in (3.4) and (3.5) is proportional to 

$$
\left(\omega_ {1} + \omega_ {2} + \omega_ {3}\right) \left(\omega_ {1} + \omega_ {2} - \omega_ {3}\right) \left(\omega_ {1} - \omega_ {2} + \omega_ {3}\right) \left(\omega_ {1} - \omega_ {2} - \omega_ {3}\right),
$$

with the $\omega_{i}$ being the frequencies. Our necessary condition that the three-wave interaction be non-resonant is equivalent to the non-vanishing of this denominator. The practical application of global-generating-function techniques tends to be a bit clumsy since (3.2) and (3.3) must be expanded to the desired order to obtain the functional relations. It is more practical to use a transformation where the generating functional is local, i.e. does not involve both old and new variables. This is most easily accomplished by using a sequence of infinitesimal transformations. This method is commonly known as Lie transforms (for a review see Cary 1981) and is most easily expressed in terms of Poisson brackets. The Poisson bracket of two quantities $A$ and $B$ is denoted by $\{A,B\}$ where 

$$
\{A, B \} = \int \mathrm {d} ^ {2} x \left[ \frac {\delta A}{\delta \zeta (x)} \frac {\delta B}{\delta \phi_ {\mathrm {s}} (x)} - \frac {\delta A}{\delta \phi_ {\mathrm {s}} (x)} \frac {\delta B}{\delta \zeta (x)} \right]. \tag {3.6}
$$

Here $\zeta(x)$ and $\phi_{\mathbf{s}}(x)$ are canonically conjugate variables. 

Lie-transform theories are similar to global canonical transformations; one seeks, by performing the transformation, a simpler Hamiltonian. Instead of considering the global generating function (3.1) and the associated canonical transformation, one considers a family of canonical transformations $\phi_{\mathrm{s}} \rightarrow \Phi(x, \lambda)$ and $\zeta \rightarrow Z(x, \lambda)$ ( $\lambda$ is the parameter describing this family) determined by a local generating function defined at $\lambda$ . There exists a functional $W[\Phi, Z]$ satisfying 

$$
\frac {\partial \mathcal {A}}{\partial \lambda} = \{A, W \}, \tag {3.7}
$$

where $A$ is either $\Phi$ or $Z$ . If $W$ has appropriate properties then $\Phi$ and $Z$ are specified as the unique solutions of (3.7) subject to the appropriate boundary conditions, which we take to be 

$$
\Phi_ {\mathrm {s}} (x, \lambda = 0) = \phi_ {\mathrm {s}} (x), \tag {3.8a}
$$

and 

$$
Z (x, \lambda = 0) = \zeta (x). \tag {3.8b}
$$

The transformed fields or new canonical variables will be taken to be solutions at $\lambda = 1$ : 

$$
\bar {\phi} _ {\mathrm {s}} (x) = \Phi (x, \lambda = 1), \tag {3.9a}
$$

and 

$$
\bar {\zeta} (x) = Z (x, \lambda = 1). \tag {3.9b}
$$

The differential equation in (3.7) can be turned into the integral equations 

$$
Z (x, \lambda) = \zeta (x) + \int_ {0} ^ {\lambda} \mathrm {d} \lambda^ {\prime} \left\{Z \left(x, \lambda^ {\prime}\right), W \right\}, \tag {3.10a}
$$

$$
\Phi (x, \lambda) = \phi_ {\mathrm {s}} (x) + \int_ {0} ^ {\lambda} \mathrm {d} \lambda^ {\prime} \left\{\Phi \left(x, \lambda^ {\prime}\right), W \right\}. \tag {3.10b}
$$

These equations can be iterated to find the new variables as functions of the old variables. To the first non-trivial order, the relations between $\bar{\phi}_{\mathrm{s}},\bar{\zeta}$ and $\phi_{\mathrm{s}},\zeta$ deduced 

from (3.10) are the same as those from (3.2)-(3.3) if $W_{\lambda}$ is chosen to be the negative of the three-wave piece of the global generating functional, (3.1), 

$$
\begin{array}{l} W = (2 \pi) ^ {2} \int \int \int \frac {\mathrm {d} ^ {2} k _ {1}}{(2 \pi) ^ {2}} \frac {\mathrm {d} ^ {2} k _ {2}}{(2 \pi) ^ {2}} \frac {\mathrm {d} ^ {2} k _ {3}}{(2 \pi) ^ {2}} \delta^ {2} (\pmb {k} _ {1} + \pmb {k} _ {2} + \pmb {k} _ {3}) \\ \times [ D (1, 2, 3) \Phi_ {1} \Phi_ {2} \Phi_ {3} + B (1, 2, 3) Z _ {1} Z _ {2} \Phi_ {3} ], \tag {3.11} \\ \end{array}
$$

using notation and definitions of the transforms of $Z(x, \lambda)$ and $\Phi(x, \lambda)$ similar to (2.12). We then have (3.7) explicitly as 

$$
\begin{array}{l} \frac {\partial Z _ {1} (\lambda)}{\partial \lambda} = (2 \pi) ^ {2} \iint \frac {\mathrm {d} ^ {2} k _ {2}}{(2 \pi) ^ {2}} \frac {\mathrm {d} ^ {2} k _ {3}}{(2 \pi) ^ {2}} \delta^ {2} \left(\boldsymbol {k} _ {1} - \boldsymbol {k} _ {2} - \boldsymbol {k} _ {3}\right) \\ \times \left[ 3 D (1, 2, 3) \Phi_ {2} (\lambda) \Phi_ {3} (\lambda) + B (2, 3, 1) Z _ {2} (\lambda) Z _ {3} (\lambda) \right], \tag {3.12a} \\ \end{array}
$$

$$
\frac {\partial \Phi_ {1} (\lambda)}{\partial \lambda} = - 2 (2 \pi) ^ {2} \iint \frac {\mathrm {d} ^ {2} k _ {2}}{(2 \pi) ^ {2}} \frac {\mathrm {d} ^ {2} k _ {3}}{(2 \pi) ^ {2}} \delta^ {2} \left(\boldsymbol {k} _ {1} - \boldsymbol {k} _ {2} - \boldsymbol {k} _ {3}\right) B (1, 2, 3) Z _ {2} (\lambda) \Phi_ {3} (\lambda). \quad (3. 1 2 b)
$$

The equivalence (to third order) of the two different transformations is enough to ensure the absence of three-wave interactions in the Lie-transformed Hamiltonian. To formally see this, recall that the new Hamiltonian $K$ , determines the time evolution of the transformed fields via 

$$
\dot {A} = \{A, K \}, \tag {3.13}
$$

where $A = (\Phi, Z)$ and the overdot refers to the time derivative along a trajectory. Using the old Hamiltonian, $H$ , we have 

$$
\dot {A} = \frac {\partial A}{\partial t} + \left\{A, H \left(\zeta , \phi_ {\mathrm {s}}\right) \right\}. \tag {3.14}
$$

Since the $A$ are functionals of $\zeta$ and $\phi_{\mathbf{s}}$ , and $W$ has no explicit dependence on the time, the first term on the right-hand side is zero. This gives 

$$
K (\Phi , Z) = H \left(\phi_ {\mathrm {s}}, \zeta\right), \tag {3.15}
$$

where $\Phi, Z$ are functionals of $\phi_{\mathrm{s}}, \zeta$ (and vice versa) to be determined from solving (3.12a,b). The physical variables $(\zeta, \phi_{\mathrm{s}})$ are given in terms of the new variables by solving the equations with the boundary conditions specified at $\lambda = 1$ ; $Z(\lambda = 1) \equiv \bar{\zeta}$ and $\Phi(\lambda = 1) \equiv \bar{\phi}_{\mathrm{s}}$ . In order to see the cancellation of the third-order terms in $K$ write (3.12a,b) as the equivalent integral equations (with boundary conditions specified at $\lambda = 1$ ) 

$$
Z (x, \lambda) = \bar {\zeta} (x) - \int_ {\lambda} ^ {1} \mathrm {d} \lambda^ {\prime} \left\{Z \left(x, \lambda^ {\prime}\right), W _ {\lambda^ {\prime}} \right\}, \tag {3.16a}
$$

and $\Phi (x,\lambda) = \bar{\phi}_{\mathbf{s}}(x) - \int_{\lambda}^{1}\mathrm{d}\lambda^{\prime}\{\Phi (x,\lambda^{\prime}),W_{\lambda^{\prime}}\} .$ (3.16b) 

Iterating these equations and using (3.8) the new Hamiltonian (to third order) is then 

$$
K (\bar {\zeta}, \bar {\phi} _ {\mathrm {s}}) = H _ {2} (\bar {\zeta}, \bar {\phi} _ {\mathrm {s}}) + H _ {3} (\bar {\zeta}, \bar {\phi} _ {\mathrm {s}}) - \left\{H _ {2} (\bar {\zeta}, \bar {\phi} _ {\mathrm {s}}), W (\bar {\zeta}, \bar {\phi} _ {\mathrm {s}}) \right\} + \dots . \tag {3.17}
$$

The last two terms cancel since 

$$
\begin{array}{l} \left\{H _ {2}, W _ {1} \right\} = (2 \pi) ^ {2} \iiint \frac {\mathrm {d} ^ {2} k _ {1}}{(2 \pi) ^ {2}} \frac {\mathrm {d} ^ {2} k _ {2}}{(2 \pi) ^ {2}} \frac {\mathrm {d} ^ {2} k _ {3}}{(2 \pi) ^ {2}} \delta^ {2} \left(\boldsymbol {k} _ {1} + \boldsymbol {k} _ {2} + \boldsymbol {k} _ {3}\right) \\ \times \left\{g \bar {\zeta} _ {3} [ 3 D (1, 2, 3) \bar {\phi} _ {1} \bar {\phi} _ {2} + B (1, 2, 3) \bar {\zeta} _ {1}, \zeta_ {2} ] - \left[ \theta_ {2} B (1, 2, 3) + \theta_ {3} B (1, 2, 3) \right] \bar {\zeta} _ {1} \bar {\phi} _ {2} \phi_ {3} \right\} \\ = \frac {1}{4} (2 \pi) ^ {2} \iiint \frac {\mathrm {d} ^ {2} k _ {1}}{(2 \pi) ^ {2}} \frac {\mathrm {d} ^ {2} k _ {2}}{(2 \pi) ^ {2}} \frac {\mathrm {d} ^ {2} k _ {3}}{(2 \pi) ^ {2}} \delta^ {2} \left(\boldsymbol {k} _ {1} + \boldsymbol {k} _ {2} + \boldsymbol {k} _ {3}\right) \bar {\zeta} _ {1} \bar {\phi} _ {2} \bar {\phi} _ {3} \left(\theta_ {2} ^ {2} + \theta_ {3} ^ {2} - \theta_ {1} ^ {2} - 2 \theta_ {2} \theta_ {3}\right), \tag {3.18} \\ \end{array}
$$

which is precisely $H_{3}$ in (2.13). Thus we secure for the Hamiltonian of the new canonical variables 

$$
K (\bar {\zeta}, \bar {\phi} _ {\mathrm {s}}) = \frac {1}{2} \int \mathrm {d} ^ {2} x [ \bar {\phi} _ {\mathrm {s}} (x, t) \theta (x) \bar {\phi} _ {\mathrm {s}} (x, t) + g \bar {\zeta} ^ {2} (x, t) ] + \dots , \tag {3.19}
$$

where ... includes higher than third-order terms. If we truncate $K(\lambda = 1)$ at second order we then have linear equations of motion in the new variables. In any problem, the initial conditions (in time) are specified in terms of the physical variables, so that (3.7), integrated from $\lambda = 0$ to $\lambda = 1$ gives the new variables to be used in (3.19). The time evolution is determined by the quadratic terms in (3.19) and is particularly simple. If one is interested in evolving the wave field to a particular time (i.e. not interested in the intervening behaviour) then this approach is useful because the hard work is done at the initial and final times in transforming between the physical and linear variables. 

# 4. One-dimensional results

We have emphasized that canonical transformations are useful when dealing with a fully evolving (in time) wave field. In order to illustrate the method and to show the quality of the new representation, we discuss some one-dimensional examples. First consider a time-independent (in the appropriate reference frame) wave, e.g. a Stokes wave. Assume that the wave field is described by a Hamiltonian in (3.19), which is only quadratic in the new, linear variables. To describe a Stokes wave, prescribe these new variables to have only one wavenumber component. Then (3.7) integrated from $\lambda = 1$ to $\lambda = 0$ will give the physical variables describing the waveform. In one dimension the functions $B$ and $D$ are particularly simple; $D$ is zero and $B$ is 

$$
B (1, 2, 3) = - \left| k _ {3} \right| \frac {1}{2} \hat {k} _ {1} \hat {k} _ {2}, \tag {4.1}
$$

where, in one dimension, $\hat{k}$ denotes $\operatorname{sgn}(k)$ . Because of the nature of the $B$ and $D$ coefficients the Lie-transform equation can be simply written using the Hilbert transform of $Z$ 

$$
\tilde {Z} (x, \lambda) = \frac {1}{\pi} P \int \frac {\mathrm {d} y Z (y , \lambda)}{x - y}, \tag {4.2}
$$

and similarly for $\Phi$ . Here $P \cap$ is the principal value integral. In Fourier-transformed variables, for example, we have 

$$
\tilde {Z} _ {k} (\lambda) = \int \mathrm {d} x \mathrm {e} ^ {- \mathrm {i} k x} \tilde {Z} (x, \lambda) = - \mathrm {i} \hat {k} Z _ {k} (\lambda). \tag {4.3}
$$

In the linear approximation $\tilde{Z}$ represents the horizontal displacement of the fluid elements in the wave (i.e. the Hilbert transform turns a cosine function into a sine 

function). In the general case, it is still useful to think of $\tilde{Z}$ as being approximately the horizontal displacement. One feature of the Hilbert transform is that the magnitude of the complex function 

$$
f (x) = Z (x, \lambda) + \mathrm {i} \tilde {Z} (x, \lambda) \tag {4.4}
$$

represents the envelope of a wave packet at a given $\lambda$ . For example, if $Z(x)$ is a product of an oscillating wave and a more slowly varying function 

$$
Z (x) = g (x) \cos (k x), \tag {4.5}
$$

then $f(x)\approx g(x)\mathrm{e}^{\mathrm{i}kx},$ (4.6) 

whose magnitude is just the slowly varying function $g$ . In terms of $\tilde{Z}(x, \lambda)$ and $\tilde{\varPhi}(x, \lambda)$ and using (4.1), (3.12) becomes (in real space) 

$$
\frac {\partial}{\partial \lambda} \tilde {Z} (x, \lambda = - \frac {1}{2} \frac {\partial}{\partial x} \tilde {Z} ^ {2} (x, \lambda), \tag {4.7a}
$$

and $\frac{\partial}{\partial\lambda}\tilde{\varPhi} (x,\lambda) = -\tilde{Z} (x,\lambda)\frac{\partial}{\partial x}\tilde{\varPhi} (x,\lambda).$ (4.7 b) 

According to Whitham (1974), (4.7) can be solved using the method of characteristics, yielding 

$$
\tilde {Z} (x, \lambda) = \tilde {Z} _ {0} (\chi), \tag {4.8a}
$$

and $\tilde{\varPhi} (x,\lambda) = \tilde{\varPhi}_{0}(\chi),$ (4.8b) 

where $\chi$ is a shifted horizontal coordinate determined by the equation 

$$
x = \chi + \lambda \tilde {Z} _ {0} (\chi), \tag {4.9a}
$$

or $x = \chi +\lambda \tilde{Z} (x,\lambda),$ (4.9b) 

Here the functions $\tilde{Z}_0(x)$ and $\tilde{\varPhi}_{\mathfrak{o}}(x)$ are the specified boundary conditions at $\lambda = 0$ . Specifying the boundary conditions at $\lambda = 1$ would alter (4.9a) to 

$$
x = \chi - (1 - \lambda) \tilde {Z} _ {0} (\chi), \tag {4.10}
$$

with $\tilde{Z}_0$ and $\tilde{\Phi}_0$ different given functions (at $\lambda = 1$ ). Taking the inverse of (4.2) allows us to express the new, linear variables (i.e. $\lambda = 1$ ) in terms of the physical variables as 

$$
\bar {\zeta} (x) = \frac {1}{\pi} P \int \frac {\mathrm {d} y}{x - y} \tilde {\zeta} (\chi), \tag {4.11a}
$$

and $\bar{\phi}_{\mathrm{s}}(x) = \frac{1}{\pi} P\int \frac{\mathrm{d}y}{x - y}\tilde{\phi}_{\mathrm{s}}(\chi),$ (4.11b) 

where $\chi$ is determined by the implicit equation (4.9a) evaluated at $\lambda = 1$ , $y = \chi + \zeta(\chi)$ . By shifting the integration variable to $\chi$ these can be written 

$$
\bar {\zeta} (x) = \frac {1}{\pi} P \int \frac {\mathrm {d} \chi [ 1 + \tilde {\zeta} ^ {\prime} (\chi) ]}{x - \chi - \tilde {\zeta} (\chi)} \tilde {\zeta} (\chi), \tag {4.12a}
$$

and $\bar{\phi}_{\mathrm{s}}(x) = \frac{1}{\pi} P\int \frac{\mathrm{d}\chi[1 + \tilde{\zeta}'(\chi)]}{x - y - \tilde{\zeta}(\chi)}\tilde{\phi}_{\mathrm{s}}(\chi).$ (4.12b) 

Expressing the physical variables in terms of the linear variables of the truncated 

Hamiltonian (which we now express as $\zeta_0$ and $\phi_{\mathbf{s0}}$ ) in a similar fashion yields (using (4.10) evaluated at $\lambda = 0$ ) 

$$
\zeta (x) = \frac {1}{\pi} P \int \frac {\mathrm {d} y [ 1 - \tilde {\zeta} _ {0} ^ {\prime} (y) ]}{x - y + \tilde {\zeta} _ {0} (y)} \tilde {\zeta} _ {0} (y), \tag {4.13a}
$$

and $\phi_{\mathrm{s}}(x) = \frac{1}{\pi} P\int \frac{\mathrm{d}y[1 - \bar{\zeta}_0(y)]}{x - y + \zeta_0(y)}\tilde{\phi}_{\mathrm{s0}}(y).$ (4.13b) 

In the remainder of this section, and the following two, we shall only be concerned with (4.13), i.e. given the new, linear variables what are the physical variables? 

The Fourier components of the physical variables are easily determined from the transforms of (4.13) as 

$$
\begin{array}{l} \zeta_ {k} = \mathrm {i} \hat {k} \int \mathrm {d} y (1 - \zeta_ {0} ^ {\prime}) \mathrm {e} ^ {- \mathrm {i} k [ y - \tilde {\zeta} _ {0} (y) ]} \tilde {\zeta} _ {0} (y), \\ = \frac {1}{| k |} \int \mathrm {d} y \mathrm {e} ^ {- \mathrm {i} k [ y - \tilde {\zeta} _ {0} (y) ]} \tilde {\zeta} _ {0} ^ {\prime} (y) \\ = \frac {1}{| k |} \int \mathrm {d} y \mathrm {e} ^ {- \mathrm {i} k y} \left[ \mathrm {e} ^ {\mathrm {i} k \bar {\zeta} _ {0} (y)} - 1 \right] \tag {4.14a} \\ \end{array}
$$

and $\phi_{\mathrm{sk}} = \mathrm{i}\hat{k}\int \mathrm{d}y(1 - \tilde{\zeta}_0^{\prime})\mathrm{e}^{-\mathrm{i}k[y - \tilde{\zeta}_0(y)]}\tilde{\phi}_{\mathrm{s0}}(y),$ 

$$
= \frac {1}{| k |} \int \mathrm {d} y \mathrm {e} ^ {- \mathrm {i} k [ y - \tilde {\zeta} _ {0} (y) ]} \tilde {\phi} _ {\mathrm {s} 0} ^ {\prime} (y). \tag {4.14b}
$$

The second forms of both these equations are obtained by noting that the $(1 - \tilde{\zeta}^{\prime})$ can be interpreted as the derivative of the exponential and then integrating by parts. The final form in (4.14a) comes from another integration by parts, noting that there should be no singular contribution at $k = 0$ . As an aside note that any derivatives, such as spatial gradients, of the fields are easily expressed in terms of the corresponding derivatives of the linear variables; multiplying (4.14) by $ik$ and Fourier transforming back to $x$ -space yields 

$$
\zeta^ {\prime} (x) = \frac {1}{\pi} P \int \mathrm {d} y \frac {\tilde {\zeta} _ {0} (y)}{x - y + \tilde {\zeta} _ {0} (y)}, \tag {4.15a}
$$

and $\phi_{\mathrm{s}}^{\prime}(x) = \frac{1}{\pi} P\int \mathrm{d}y\frac{\tilde{\phi}_{\mathrm{s0}}^{\prime}(y)}{x - y + \zeta_{0}(y)}.$ (4.15b) 

Finally we note that the time evolution of the nonlinear, physical fields is easily incorporated into our formula by just letting the linear variables evolve according to their quadratic Hamiltonian. 

Our approximation to a Stokes wave is obtained by $\tilde{\zeta}_0(y) = a\sin (ky - \omega t)$ and $\tilde{\phi}_{\mathrm{so}}(y) = -(a / k)\cos (ky - \omega t)$ . These are the waves of permanent form in a linear theory. The comparison of surface height for an exact Stokes wave and (4.13a) is shown in figure 1. In this figure the dashed line is a sinusoidal, linear wave of height $0.2\mathrm{m}$ and wavelength $2\pi \mathrm{m}$ . The small circles represent our transformed wave (4.13a) while the solid line is a Stokes wave which was calculated using the expansion given by Schwarz (1974). It can be seen that the Lie transform presents a very good approximation to the true Stokes waveform. In figure 2(a,b) we have similar plots for the orbital surface velocities. Figure 2(a) compares the linear (dashed), Lie transformed (dotted), and Stokes (solid) horizontal velocity. The Lie-transformed 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-12/8237f025-efca-4750-807b-3f5bc5cb85c3/09a95c2973f94c73208b59fc2d99b226a5ec3be0bba213245052f456f9466102.jpg)



FIGURE 1. Comparison of linear (dash), Lie-transformed (circles), and Stokes (solid) waves for a wave slope, $k_{0}A$ , of 0.2. The fundamental wavenumber, $k_{0}$ , was chosen to be unity so that in these units the wave height, $A$ , is 0.2.


horizontal velocity, $u$ , is determined from $u \equiv \mathrm{D}_x\phi$ . This expression can be evaluated using the equations of motion derived from (2.9) 

$$
\dot {\zeta} = \left[ \mathbf {D} _ {z} - \left(\partial_ {x} \zeta\right) \mathbf {D} _ {x} \right] \phi \tag {4.16}
$$

and the identity (2.7), yielding 

$$
u = \frac {\left[ \phi_ {s} ^ {\prime} (x) - \zeta^ {\prime} (x) \dot {\zeta} \right]}{1 + \left(\zeta^ {\prime}\right) ^ {2}}, \tag {4.17}
$$

$\dot{\zeta}$ is given by taking the time derivative of $(4.13a)$ , which after some integration by parts becomes 

$$
\dot {\zeta} (x) = \frac {1}{\pi} P \int \mathrm {d} y \frac {\tilde {\zeta} _ {0}}{x - y + \tilde {\zeta} _ {0} (y)}. \tag {4.18}
$$

Figure 2(b) does the same for the vertical velocity, $v \equiv \mathrm{D}_z\phi$ , which again using (4.16) and (2.7) becomes 

$$
v = \frac {\left[ \dot {\zeta} (x) + \zeta^ {\prime} (x) \phi_ {s} ^ {\prime} (x) \right]}{1 + \left(\zeta^ {\prime}\right) ^ {2}}. \tag {4.19}
$$

It should also be noted that in (4.10) evaluated at $\lambda = 0$ , the difference between $x$ and $\chi$ can be interpreted as the horizontal Lagrangian displacement. Since this difference is the Hilbert transform of the vertical height, this shows that the orbital motion of fluid particles is purely circular. For steady, progressive waves this interpretation is consistent with the analysis of Longuet-Higgins (1986). He showed that the Lagrangian displacement follows the motion of a cycloid. The horizontal displacement consists of a piece equivalent to our displacement (to leading order) and another piece which grows secularly in time, being due to the Stokes drift. Our transformation, being time independent, cannot contain such secular growth. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-12/8237f025-efca-4750-807b-3f5bc5cb85c3/ff8000b286d6c08ae2fe5c4be80528e354162ec18c418a79f04c70cb38ebbe53.jpg)



FIGURE 2. Comparison of (a) horizontal and (b) vertical surface velocities for linear (dash), Lie-transformed (circles), and Stokes (solid) waves, for the same parameters as figure 1.


The Lie transformation introduced in the preceding section has implicit higher-order interactions included in it. One may question whether these interactions preserve locality. By this, we mean that a localized (in space) wave packet should not be transformed to a wave packet with non-local pieces to it (i.e. significant pieces outside the packet). In figure 3 we compare a Gaussian wave packet (solid line) to its Lie transform (dotted line) when the total wave slope is approximately one-fourth. If anything, it can be seen that the transformed packet is narrower than the original packet, and so the Lie transform indeed preserves locality. 

Since our transformation induces nonlinearities into a linear wave field one can look at the nature of breaking waves in our theory. One cannot expect to have precisely the same wave breaking conditions as for a Stokes wave since these conditions depend on the form of high-order nonlinearities. While our model incorporates the correct leading order nonlinearity and approximates low-order nonlinearities well, there is no reason to suppose that high-order terms are even approximately correct. In the method of characteristics (cf. (4.8)) wave breaking (or, equivalently, multiple valuedness of the wave height) occurs when the relation between the horizontal coordinate and the shifted coordinate $\chi$ becomes indeterminate (see Whitham 1974 for a discussion of this point). Setting $\lambda = 0$ and differentiating (4.10) with respect to $x$ yields 

$$
1 = \frac {\partial \chi}{\partial x} \left[ 1 - \frac {\partial \tilde {Z} _ {0} (\chi)}{\partial \chi} \right]. \tag {4.20}
$$

Breaking occurs when the term in the brackets vanishes. For a single-component wave, e.g. $Z_0(x) = a\cos kx$ we obtain the condition 

$$
1 = k a \cos k \chi . \tag {4.21}
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-12/8237f025-efca-4750-807b-3f5bc5cb85c3/3889681bc7b09579896723b3540eff55d5f995499967588cf34776f6fed00c7d.jpg)



FIGURE 3. Effects of the Lie transform (dots) on a Gaussian wave packet (solid).


This shows that breaking occurs when the wave slope $(ka)$ equals one, as compared to a true Stokes wave where the wave breaking occurs when the wave slope is about 0.41. 

One important feature of nonlinear surface gravity waves is the appearance of the Stokes-drift current and an attendant return current. This situation has been analysed by Longuet-Higgins & Stewart (1964) who note that the Stokes drift and the return current combine to ensure that there is no change in the surface elevation (i.e. a $k \approx 0$ contribution to the surface wave field), as would be caused by mass transport. One feature of the new representation is that, to lowest order, the return current does balance the Stokes drift and that there is no change of in-surface elevation. This follows from (4.14). Consider a narrow-bandwidth group of waves, such that the nonlinear interactions generate wave components with nearly zero wavenumber. We expand the velocity potential, $\phi$ , in a Fourier series 

$$
\phi (x, z) = \int \frac {\mathrm {d} k}{2 \pi} \mathrm {e} ^ {\mathrm {i} k x} \mathrm {e} ^ {| k | z} \phi_ {k}. \tag {4.22}
$$

The narrow band of waves is centred about some $k_{0}$ with important wave vectors in a band $\Delta k$ about $k_{0}$ , so that the nonlinear interactions generate Fourier components in the range $-\Delta k \lesssim k \lesssim \Delta k$ . It is apparent from (4.22) that the narrow-band component extends to a depth of $1 / k_{0}$ , whereas the small- $k$ terms extend to a depth of $1 / \Delta k$ . We shall show that the currents due to these two terms are equal in magnitude, but opposite in sign. 

In order to obtain the Fourier components of $\phi$ and $\zeta$ in the $\Delta k$ -band, we expand (4.14) for small $k$ 

$$
\zeta_ {k} = \mathrm {i} \hat {k} \int \mathrm {d} y \mathrm {e} ^ {- \mathrm {i} k y} \tilde {\zeta} _ {0} (y) \tilde {\zeta} _ {0} ^ {\prime} (y) \tag {4.23}
$$

and $\phi_{\mathrm{sk}} = \mathrm{i}\hat{k}\int \mathrm{d}y\mathrm{e}^{-\mathrm{i}ky}\tilde{\zeta}_0(y)\tilde{\phi}_{\mathrm{so}}'(y).$ (4.24) 

We observe from (4.3) that the integral of a product of two Hilbert transforms is equal to the integral of the product of the original functions, so that for small $k$ and $k_0 \gg k$ the above equations become 

$$
\zeta_ {k} = \mathrm {i} \hat {k} \int \mathrm {d} y \mathrm {e} ^ {- \mathrm {i} k y} \zeta_ {0} (y) \zeta_ {0} ^ {\prime} (y) \tag {4.25}
$$

and $\phi_{\mathrm{sk}} = \mathrm{i}\hat{k}\int \mathrm{dy} \mathrm{e}^{-\mathrm{i}ky}\zeta_0(y)\phi_{\mathrm{s0}}'(y).$ (4.26) 

The total horizontal current through the vertical plane located at some $x$ is given by 

$$
V (x) = \int_ {- \infty} ^ {\zeta (x)} \frac {\partial \phi (x , z)}{\partial x} d z. \tag {4.27}
$$

For a sufficiently narrow band of waves, there are two contributions to $V(x)$ . The near-surface contribution is readily evaluated from (4.27) to be $\zeta(x) (\partial \phi_{\mathrm{s}} / \partial x)$ , which when averaged over a wavelength $\lambda = 2\pi / \mathbf{k}_0$ gives the familiar Stokes current 

$$
V _ {\text {S t o k e s}} (x) = \left\langle \zeta (x) \frac {\partial \phi_ {\mathrm {s}}}{\partial x} \right\rangle , \tag {4.28}
$$

where $\langle \rangle$ denotes a spatial average over a wavelength, $\lambda$ , about $x$ . To evaluate the deep current we insert (4.22) into (4.27) and keep only leading terms in the bandwidth 

$$
\begin{array}{l} V _ {\mathrm {D e e p}} (x) = \int_ {- \infty} ^ {0} \mathrm {d} z \int_ {- \Delta k} ^ {\Delta k} \frac {\mathrm {d} k}{2 \pi} \mathrm {i} \hat {k} \mathrm {e} ^ {\mathrm {i} k x} \mathrm {e} ^ {| k | z} \phi_ {k} \\ \approx \int \frac {\mathrm {d} k}{2 \pi} \frac {\mathrm {i} k}{| k |} \mathrm {e} ^ {\mathrm {i} k x} \phi_ {k}, \tag {4.29} \\ \end{array}
$$

where the integral is only over $k$ -values such that $k \ll k_0$ . Inserting the expression for $\phi$ , (4.26), we obtain 

$$
\begin{array}{l} V _ {\text {D e e p}} (x) = \mathrm {i} \int \frac {\mathrm {d} k}{2 \pi} \hat {k} \mathrm {e} ^ {\mathrm {i} k x} \left[ \mathrm {i} \hat {k} \int \mathrm {d} y \mathrm {e} ^ {- \mathrm {i} k y} \zeta_ {0} (y) \phi_ {\mathrm {s 0}} ^ {\prime} (y) \right] \\ = - \left\langle \zeta_ {0} (x) \phi_ {\mathbf {s} 0} ^ {\prime} (x) \right\rangle , \tag {4.30} \\ \end{array}
$$

which is minus the near-surface Stokes current (4.28). Thus the deep return current, which arises from nonlinear interactions, correctly cancels the surface Stokes drift. 

Finally, we see from (4.25) that the small- $k$ Fourier components of $\zeta_{k}$ are proportional to $k$ and hence there is no average (i.e. constant) surface elevation. This is in accordance with the analysis of Longuet-Higgins & Stewart (1964). 

# 5. Short waves on long waves

Next we consider the effects that the Lie transform has on the situation of a narrow wave packet riding on a larger, longer wave. With the assumption of good scale separation we write 

$$
\zeta (x) = \zeta_ {\mathrm {L}} (x) + \zeta_ {\mathrm {S}} (x), \tag {5.1}
$$

and a similar expansion for $\tilde{\zeta}(x)$ . Here $\zeta_{\mathrm{L}}$ contains only long-wavelength components while $\zeta_{\mathrm{S}}$ contains only short-wavelength components. It is assumed that $\zeta_{\mathrm{L}}$ is given in terms of $\tilde{\zeta}_{0\mathrm{L}}$ by (4.13). The expansions can be inserted into (4.14a) which, upon keeping only linear terms is $\zeta_{\mathrm{S}}$ and $\zeta_{0\mathrm{S}}$ , becomes (in wavenumber space) 

$$
\zeta_ {\mathrm {S}} (k) = \mathrm {i} \hat {k} \int \mathrm {d} y \mathrm {e} ^ {- \mathrm {i} k [ y - \tilde {\xi} _ {0 L} (y) ]} \tilde {\xi} _ {0 \mathrm {S}} (y). \tag {5.2}
$$

Because of the scale separation we can expand $\zeta_{0\mathrm{L}}(x)$ about the centre of the wave packet, which we take to be $(y_0)$ , 

$$
\tilde {\zeta} _ {0 \mathrm {L}} (y) \approx \tilde {\zeta} _ {0 \mathrm {L}} \left(y _ {0}\right) + \left(y - y _ {0}\right) \tilde {\zeta} _ {0 \mathrm {L}} ^ {\prime} + \dots , \tag {5.3}
$$

where we have discarded wave curvature and higher-order terms. Performing the $y$ -integration yields 

$$
\begin{array}{l} \zeta_ {\mathrm {S}} (k) = \mathrm {e} ^ {- \mathrm {i} k \left(y _ {0} - \tilde {\zeta} _ {0 \mathrm {L}} \left(y _ {0}\right)\right)} \int \mathrm {d} y \mathrm {e} ^ {- \mathrm {i} k \left(y - y _ {0}\right) \left(1 - \zeta_ {0 \mathrm {L}} ^ {\prime} \left(y _ {0}\right) \right.} \int \mathrm {d} p (\hat {k} \hat {p}) \mathrm {e} ^ {\mathrm {i} p \left(y - y _ {0}\right)} \zeta_ {0 \mathrm {S}} (p) \\ = \zeta_ {0 \mathrm {S}} \left[ k \left(1 - \tilde {\zeta} _ {0 \mathrm {L}} ^ {\prime} \left(y _ {0}\right)\right) \right] \mathrm {e} ^ {- \mathrm {i} k \left[ y _ {0} - \tilde {\zeta} _ {0 \mathrm {L}} \left(y _ {0}\right) \right]}. \tag {5.4} \\ \end{array}
$$

Fourier transforming back to coordinate space we secure for the physical wave height 

$$
\zeta_ {\mathrm {S}} (x) = \frac {1}{| 1 - \tilde {\zeta} _ {0 \mathrm {L}} ^ {\prime} \left(y _ {0}\right) |} \zeta_ {0 \mathrm {S}} \left(\frac {x - \tilde {\zeta} _ {0 \mathrm {L}} \left(y _ {0}\right)}{1 - \tilde {\zeta} _ {0 \mathrm {L}} ^ {\prime} \left(y _ {0}\right)}\right), \tag {5.5}
$$

where $\zeta_{\mathbf{os}}(x)$ is the wave function describing the input packet. We immediately note the following points: 

(i) The centre of the wave packet is shifted by an amount $\tilde{\zeta}_{0\mathrm{L}}(y_0)$ , i.e. by the Hilbert transform of the height of the linear, long wave. As an example consider a big wave of the form 

$$
\zeta_ {0 \mathrm {L}} \left(y _ {0}\right) = a \cos k y _ {0}. \tag {5.6}
$$

Then the Hilbert transform is 

$$
\tilde {\zeta} _ {0 \mathrm {L}} \left(y _ {0}\right) = a \sin k y _ {0}. \tag {5.7}
$$

In the linear approximation $\zeta_{0\mathrm{L}}$ is the vertical displacement while $\tilde{\zeta}_{0\mathrm{L}}$ is the horizontal displacement. Thus the shift in (5.5) is precisely the correct shift as given by the long-wave advection. We see that there is no shift of the wave packet at the crest $(y_0 = 0)$ or trough $(y_0 = \pi /k)$ of the big wave. There is a maximal shift on the sides $(y_0 = \pm \pi /2k)$ of the big wave. 

(ii) There is a shift in the short-wave wavenumber due to the term $1 - \tilde{\zeta}_{0\mathrm{L}}$ in the argument of $\tilde{\zeta}_{\mathrm{s}}$ . Since $\tilde{\zeta}_{0\mathrm{L}}'$ is the long-wave straining, this is the correct wavenumber shift. 

(iii) There is a modulation of the width of the packet (and a corresponding change in the height so that the area under the packet remains constant), which is just $1 - ka \cos ky_0$ for the wave considered in (i). Thus the packet becomes narrower and higher at the crest and the opposite at the trough. 

(iv) The form of the wave action is preserved. Of course, the actual value of the action remains constant under any canonical transformation. For short waves on long waves the action is given by an expression similar to that for linear waves with the gravitational constant modified (Garrett & Smith 1976; Henyey et al. 1988) 

$$
A = \frac {1}{2} \int \mathrm {d} x g ^ {\prime} \zeta_ {\mathrm {S}} ^ {2} (x), \tag {5.8}
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-12/8237f025-efca-4750-807b-3f5bc5cb85c3/1c2e97e3c4f9ca056f3d20459f567f36824c53a5f619dd13bc5a50b77dda52fb.jpg)



FIGURE 4. Surface-height field of the long wave (dash) with unit fundamental wavenumber $(k_0)$ and the wave packet (solid) which sits on the long wave.


where the modified gravity is 

$$
g ^ {\prime} = g \left(1 - \tilde {\zeta} _ {0 \mathrm {L}} ^ {\prime}\right), \tag {5.9}
$$

to the leading order in the long-wave wave slope. Inserting (5.5) into (5.8) gives 

$$
A = \frac {1}{2} \int \mathrm {d} x g \zeta_ {\mathrm {0 S}} ^ {2} (x), \tag {5.10}
$$

showing that the value of the action is preserved as well as its form. Preservation of the value of the action by a canonical transformation such as (3.1) requires a modification of the form of the action, (5.8). 

As a specific example consider the linear (i.e. new variable) long wave to have wavelength $2\pi \mathrm{m}$ (i.e. the wavenumber, $k = 1$ ) and the short wave initially to be a Gaussian packet of the form (in wavenumber space) 

$$
\frac {A}{(2 \pi) ^ {\frac {1}{2}} W} \mathrm {e} ^ {- \frac {1}{2} [ (k - k _ {0}) / W ] ^ {2}}, \tag {5.11}
$$

where $W$ denotes the bandwidth, $A$ is the small-wave height, and $k_{0}$ is the central wavenumber. For definiteness, we take $A = 10^{-5} \mathrm{~m}$ , $k_{0} = 72 \mathrm{~m}^{-1}$ and $W = 48 \mathrm{~m}^{-1}$ . This large value of the central wavenumber assures good scale separation while the large bandwidth implies that we can investigate the nature of the interaction between long and short waves as a function of the phase of the long wave. We took the long-wave height to be $0.2 \mathrm{~m}$ implying a wave slope of 0.2. Thus, although interactions of the long wave with itself are important in building up the nonlinear (i.e. Stokes) wave, we want to focus on the strong interactions between waves with widely separated scales. Figure 4 shows the initial wave height $(\zeta_{\mathrm{0L}} + \zeta_{\mathrm{0S}})$ and the 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-12/8237f025-efca-4750-807b-3f5bc5cb85c3/4d9f858f9fdeec955ac81bad1c9647cf086580a1a1d95cf7288e51e3af120835.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-12/8237f025-efca-4750-807b-3f5bc5cb85c3/611bd80184301952fbe92100e098b20a002badc655204851619151ed49779178.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-12/8237f025-efca-4750-807b-3f5bc5cb85c3/5ead78752179debf3a3d332d6a4284d5404c4bebdb560dbf17c660e8a0a3e4b9.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-12/8237f025-efca-4750-807b-3f5bc5cb85c3/1ceb4be889026239a553785d2d475970b87f7b2c774b72e929638bcee3fe1873.jpg)



Position/{large-wave height}



FIGURE 5. Linear (dots) and Lie-transformed (solid) envelope functions squared, $|f(x)|^2$ as a function of the phase: (a) $y_0 = 0$ , (b) $y_0 = \frac{1}{2}\pi$ , (c) $y_0 = \pi$ , (d) $y_0 = -\frac{1}{2}\pi$ . Here $y_0$ is the position of the centre of the wave packet expressed in terms of the long-wave phase ( $y_0 = 0$ is the crest of the long wave).


form of the Gaussian packet $[\zeta_{0\mathrm{S}}(x)]$ for the situation where the packet sits at the crest of the long wave (of the form (5.6)). The total wave height was taken as the initial condition for $(4.13a)$ and the transformed result had the physical $\zeta_{\mathrm{L}}$ subtracted from it. Figure 5 shows the envelope function squared, $|f(x)|^2$ , for the initial $(\zeta_{0\mathrm{S}})$ and Lie-transformed result $(\zeta - \zeta_{\mathrm{L}})$ , at four different values for the phase of the long wave, respectively, $y_0 = 0$ , $\frac{1}{2}\pi$ , $\pi$ , $\frac{3}{2}\pi$ ( $y_0 = 0$ is the crest of the long wave and the wave number, $k$ , is one). It is seen that there is no shift of the wave packet at the crest ( $y_0 = 0$ ) or the trough ( $y_0 = \pi$ ) while the shift at the phase points $y_0 = \frac{1}{2}\pi$ and $y_0 = \frac{3}{2}\pi$ is seen to be given by our formula (5.5). Thus our scale-separated result, (5.5), well approximates the full canonical transformation. 

In this situation the Lie-transformed field is seen to provide better results than some other canonical transformations such as those generated by the global generating function (3.1). As can be seen from figure 5 the Lie-transformed fields give the correct dynamics for a wave packet on a long wave. This is because the dispersion 

relation for the non-physical, linear variables $\bar{\zeta},\bar{\phi}_{\mathrm{s}}$ is unaltered, implying that the wave speed is unchanged. However, it is known that a physical wave packet will speed up or slow down (relative to this constant wave speed) according to the corresponding phase of the long wave at that point. The evolution of our Lie-transformed wave packet has just this behaviour, but it appears as a purely kinematic effect. Global canonical transformation does not give this kinematic behaviour, though the shift of the centre would be approximately correct in the circumstance when the actual shift (which is proportional to the long-wave height) is much less than the width of the packet (which is roughly the short-wave wavenumber). However, when the combination of the shift and width, $k_{\mathrm{S}}\tilde{\zeta}_{0\mathrm{L}}$ is comparable with one, the coherence of the wave packet is destroyed. Our Lie transform seems to properly treat the important interactions (at least in one dimension) between short and long waves to all orders in $k_{\mathrm{S}}\tilde{\zeta}_{0\mathrm{L}}$ . 

# 6. Surface-height spectrum

As discussed in the introduction, since the new variables have linear dynamics it is possible that they have Gaussian statistics. Assuming this behaviour, we now investigate what this implies about the spectrum and statistics of the physical surface height. We shall look at one-dimensional spectra only. 

The surface-wave spectrum is defined as 

$$
\Phi (k) = \int \mathrm {d} x \mathrm {e} ^ {\mathrm {i} k x} \left\langle \zeta \left(\frac {1}{2} x\right) \zeta \left(- \frac {1}{2} x\right) \right\rangle , \tag {6.1}
$$

where $\langle \rangle$ denotes ensemble average. The homogeneity assumed in (6.1) also follows from 

$$
\left\langle \zeta_ {k} \zeta_ {p} \right\rangle \equiv 2 \pi \delta (k + p) \Phi (k). \tag {6.2}
$$

This form is more convenient since we can immediately use (4.14) in (6.2) 

$$
\left\langle \zeta_ {k} \zeta_ {p} \right\rangle = \frac {1}{| k p |} \int \mathrm {d} y _ {1} \mathrm {d} y _ {2} \mathrm {e} ^ {- 1 (k y _ {1} + p y _ {2})} \left\langle \left(\mathrm {e} ^ {\mathrm {i} k \tilde {\zeta} _ {0} (y _ {1})} - 1\right) \left(\mathrm {e} ^ {\mathrm {i} p \tilde {\zeta} _ {0} (y _ {2})} - 1\right) \right\rangle . \tag {6.3}
$$

The $\xi_0$ in (6.3) are assumed to be Gaussian random processes, allowing us to use 

$$
\langle \mathrm {e} ^ {\mathrm {i} a _ {5} ^ {2}} \rangle = \mathrm {e} ^ {- \frac {1}{8} a ^ {2} \langle \xi^ {2} \rangle}, \tag {6.4}
$$

valid for Gaussian variables. Thus 

$$
\langle \zeta_ {k} \zeta_ {p} \rangle = \frac {1}{| k p |} \int \mathrm {d} y _ {1} \mathrm {d} y _ {2} \mathrm {e} ^ {- \mathrm {i} (k y _ {1} + p y _ {2}) \{1 + \mathrm {e} ^ {- \frac {1}{2} \langle [ k \tilde {\zeta} _ {0} (y _ {1}) + p \tilde {\zeta} _ {0} (y _ {2}) ] ^ {2} \rangle} - \mathrm {e} ^ {- \frac {1}{2} k ^ {2} \langle \tilde {\zeta} _ {0} ^ {2} (y _ {1}) \rangle} - \mathrm {e} ^ {- \frac {1}{2} p ^ {2} \langle \tilde {\zeta} _ {0} ^ {2} (y _ {2}) \rangle} \}}. \tag {6.5}
$$

Because of homogeneity $\langle \tilde{\zeta}_0(x)\tilde{\zeta}_0(y)\rangle$ is a function of $x - y$ only, allowing us to explicitly perform one of the integrations. Thus, from (6.2), we find 

$$
\Phi (k) = \frac {1}{k ^ {2}} \int \mathrm {d} y \mathrm {e} ^ {- \mathrm {i} k y} \{1 + \mathrm {e} ^ {- \frac {1}{3} k ^ {2} \langle [ \tilde {\zeta} _ {0} (\frac {1}{2} y) - \tilde {\zeta} _ {0} (- \frac {1}{2} y) ] ^ {2} \rangle} - 2 \mathrm {e} ^ {- \frac {1}{3} k ^ {2} \langle \zeta_ {0} ^ {2} (0) \rangle} \}. \tag {6.6}
$$

Defining the correlation function as the Fourier transform of the spectrum of linear variables, (6.1), 

$$
C _ {0} (y) = \left\langle \zeta_ {0} \left(\frac {1}{2} y\right) \zeta_ {0} \left(- \frac {1}{2} y\right) \right\rangle , \tag {6.7}
$$

we find the correlation function of the Hilbert-transformed variables is the same as for the surface height (by using (4.3)) 

$$
\left\langle \tilde {\xi} _ {0} \left(\frac {1}{2} y\right) \tilde {\xi} _ {0} \left(- \frac {1}{2} y\right) \right\rangle = C _ {0} (y). \tag {6.8}
$$

Using this in (6.6) the spectrum is 

$$
\Phi (k) = \frac {1}{k ^ {2}} \int \mathrm {d} y \mathrm {e} ^ {- \mathrm {i} k y} \left(\mathrm {e} ^ {- k ^ {2} \left[ C _ {0} (0) - C _ {0} (y) \right]} + 1 - 2 \mathrm {e} ^ {- \frac {1}{2} k ^ {2} C _ {0} (0)}\right). \tag {6.9}
$$

Noting that there is no contribution at $k = 0$ (as there should not be if $\Phi_0$ has no contribution at $k = 0$ ) we can rewrite (6.9) to secure 

$$
\Phi (k) = \frac {1}{k ^ {2}} \int \mathrm {d} y \mathrm {e} ^ {- \mathrm {i} k y} \mathrm {e} ^ {- k ^ {2} C _ {0} (0)} \left\{\mathrm {e} ^ {k ^ {2} C _ {0} (y)} - 1 \right\}. \tag {6.10}
$$

In order to interpret this result we perform an expansion of the exponentials in the integrand in powers of the quantity $k^2 C_0$ and write the result as a convolution in wavenumber space 

$$
\Phi (k) = \Phi_ {0} (k) + \frac {1}{2} k ^ {2} \int \mathrm {d} q \left\{\Phi_ {0} (q) \Phi_ {0} (k - q) - 2 \Phi_ {0} (k) \Phi_ {0} (q) \right\}. \tag {6.11}
$$

Notice that our result differs from Barrick & Weber (1976) who obtain a correction which is just the first term in the braces of (6.11). Our result (6.11) includes third-order waves (using their language), but as pointed out above these must be included for consistency and are extremely important. For large values of the spectral parameter $k$ the first term has two regions of the $q$ -integral which contribute, small $q$ and small $k - q$ . Both of these regions give a leading-order contribution of the form $\varPhi_0(k)\varPhi_0(q)$ and the sum cancels the second term in the integral. If we had used a canonical transformation such as (3.1), this cancellation would not have been as complete. That is, some other canonical transformations would have given results where the correction grows as $k^2$ . For a power-law spectrum, the residual correction (for the Lie-transformed field) in the inertial rage is the mean square wave slope times the zeroth-order spectrum. 

Comparison of the modified spectrum is done in figure 6 for a Phillips (1977) colinear spectrum of the form 

$$
\Phi_ {0} (k) = \frac {0 . 0 0 3}{k ^ {3}}, \tag {6.12}
$$

in the inertial range. We have included a roll-off of the spectrum at wavenumbers below about $0.1\mathrm{m}^{-1}$ and above $50\mathrm{m}^{-1}$ . The small circles represent the unperturbed (6.12), while the solid line represents $\varPhi(k)$ of (6.10). For comparative purposes, we show the Barrick & Weber correction term as a dashed line. Note that for large spatial wavenumbers their correction is much larger than the original spectrum and orders of magnitude larger than the consistent correction of (6.11), in line with the above observations on the incomplete cancellation. The dimensionless quantity defining the perturbation in the Barrick & Weber expansion is the short-wave wavenumber (i.e. $k$ ) times the largest-wave height (predominantly from the spectral peak). This quantity is not small as $k$ gets large. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-12/8237f025-efca-4750-807b-3f5bc5cb85c3/c675c9514c5cd0c0d56d24beff03ccc68d98b867113828203127a09cd30f5bd3.jpg)



FIGURE 6. Collinear Phillips spectrum in (6.12) (circles), its Lie transform (solid), and the Barrick & Weber correction (dash).


In order to see explicitly where the difference between Barrick & Weber's and our results comes from, we perform a perturbative expansion of the surface elevation 

$$
\begin{array}{l} \zeta_ {k} = \frac {1}{| k |} \int \mathrm {d} y \tilde {\zeta} _ {0} ^ {\prime} (y) [ 1 + \mathrm {i} k \tilde {\zeta} _ {0} (y) - \frac {1}{2} k ^ {2} \tilde {\zeta} _ {0} ^ {2} (y) + \dots ] \mathrm {e} ^ {- \mathrm {i} k y} \\ = \zeta_ {0 k} - \frac {1}{2} \mathrm {i} | k | \iint \frac {\mathrm {d} k _ {2}}{2 \pi} \delta \left(k - k _ {2} - k _ {3}\right) \hat {k} _ {2} \hat {k} _ {3} \tilde {\zeta} _ {0 k _ {2}} \tilde {\zeta} _ {0 k _ {3}} \\ - \frac {1}{6} \hat {k} k ^ {2} \iint \int \frac {\mathrm {d} k _ {2} \mathrm {d} k _ {3} \mathrm {d} k _ {4}}{(2 \pi) ^ {2}} \delta \left(k - k _ {2} - k _ {3} - k _ {4}\right) \hat {k} _ {2} \hat {k} _ {3} \hat {k} _ {4} \tilde {\zeta} _ {0 k _ {2}} \tilde {\zeta} _ {0 k _ {3}} \tilde {\zeta} _ {0 k _ {4}}. \tag {6.13} \\ \end{array}
$$

The Barrick & Weber result comes from keeping just the first two terms in (6.13) and inserting the answer into the expression for the spectrum (6.2). However, we see that the third term times the first term, which appears as a cross-term in (6.2), is the same order in wave height as the second term squared. Keeping all terms to the same order, as can be done in the Lie transformation, we obtain (6.11) where, as we have noted before, the dimensionless perturbation parameter appears to be the mean-square wave slope. For the large- $k$ portion of the wave spectrum the dominant interactions are the advection and straining by the large-scale waves (small $k$ ) and, as we have discussed in §5, our result, (5.5), consistently includes these interactions; the phase shift and wavenumber shift of the high- $k$ modes are correctly given. Conservation of action is automatically included in our result, ensuring that any modification to the spectrum is small. Other formalisms (such as the global canonical transformation or naive perturbation theory) modify the form of the action and so can lead to large corrections to the spectrum unless done carefully. A consistent 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-12/8237f025-efca-4750-807b-3f5bc5cb85c3/cb33acb59fd29c4a865a7d9417f4eeae3ffa0bb605cede8f7d3ebc52f923e6e8.jpg)



FIGURE 7. Collinear Toba spectrum in (6.14) (circles), its Lie transform (solid), and the Barrick & Weber correction (dash).


expansion should conserve wave action and thus should be an expansion in a small parameter (such as the wave slope). 

In the light of recent experimental and theoretical investigations (Toba 1973; Kitaigorodski 1983; Phillips 1985) of the spatial spectrum we show, in figure 7, comparisons for a collinear spectrum of the form 

$$
\Phi (k) = \frac {0 . 0 0 3}{k ^ {2 . 5}} \tag {6.14}
$$

in the inertial range and with roll-offs similar to those for (6.12). The small circles represent the unperturbed (6.14), while the solid line represents $\Phi(k)$ of (6.11) for this case. Again the modification of the spectrum occurs mainly as an increase in spectral strength at large wavenumbers. 

We also consider the surface-height bispectrum, $\varLambda$ 

$$
\left\langle \zeta_ {k _ {1}} \zeta_ {k _ {2}} \zeta_ {k _ {3}} \right\rangle = \delta \left(k _ {1} + k _ {2} + k _ {3}\right) \Lambda \left(k _ {1}, k _ {2}, k _ {3}\right). \tag {6.15}
$$

For Gaussian variables this quantity is zero. Using (4.14a), the averages in (6.15) are straightforwardly evaluated in terms of $C_0$ or the associated spectrum, $\Phi_0$ . Performing an expansion in $k^2 C_0$ similar to that for (6.11) we secure 

$$
\begin{array}{l} \varLambda (k _ {1}, k _ {2}, k _ {3}) = k _ {3} (\hat {k} _ {1} \hat {k} _ {3}) (\hat {k} _ {2} \hat {k} _ {3}) \varPhi_ {0} (k _ {1}) \varPhi_ {0} (k _ {2}) \\ + k _ {2} \left(\hat {k} _ {1} \hat {k} _ {2}\right) \left(\hat {k} _ {3} \hat {k} _ {2}\right) \Phi_ {0} \left(k _ {1}\right) \Phi_ {0} \left(k _ {3}\right) + k _ {1} \left(\hat {k} _ {2} 1 \hat {k} _ {1}\right) \left(\hat {k} _ {3} \hat {k} _ {1}\right) \Phi_ {0} \left(k _ {2}\right) \Phi_ {0} \left(k _ {3}\right). \tag {6.16} \\ \end{array}
$$

# 7. Two-dimensional results

In order to demonstrate the effect that our Lie transformation has on waves in two dimensions we show some numerical results for the surface height. In figure 8(a) we 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-12/8237f025-efca-4750-807b-3f5bc5cb85c3/d009461214bd0353e4d82d9ab0148888319d4fe44d9ef6680bb86252abc897c0.jpg)



$|k|$ $x$ -position


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-12/8237f025-efca-4750-807b-3f5bc5cb85c3/433c4464c45fe51986e286c2dff7df0711c5b84edd9a118ffd69ce23ca4451bf.jpg)



$|k|$ x-position



FIGURE 8. (a) Contour plot of two-dimensional surface elevation for the sum of two cosine waves, with contours representing the wave slope, $k \times$ wave height, in intervals of 0.02. The two wave vectors are $(k_{1x}, k_{1y}) \approx (0.924, 0.3827)$ and $(k_{2x}, k_{2y}) \approx (0.924, -0.3827)$ . Both wave vectors are of unit magnitude. (b) Contour plot of two-dimensional surface elevation for the Lie transform of the sum of the waves in (a). Again, the contours are in intervals of 0.02.


present a contour plot of the surface-height field for the linear variables. This surface-height field consists of the sum of two cosine waves. The first wave has a wave vector $k_{1x} \approx 0.924 \, \mathrm{m}^{-1}$ and $k_{1y} \approx 0.3827 \, \mathrm{m}^{-1}$ . The second wave has a wave vector $k_{2x} \approx 0.924 \, \mathrm{m}^{-1}$ and $k_{2y} \approx -0.3827 \, \mathrm{m}^{-1}$ . The crest of the sum of the waves occurs at the origin. In figure 8(b) we show the results of inserting this initial height field (and the corresponding velocity potential) into (3.12) and integrating from $\lambda = 1$ to $\lambda = 0$ to obtain the physical surface elevation. The contours in the two figures represent the same values and so we take note of the peaking of the crest at the origin and the flattening of the surrounding troughs. 

Since, in the new variables, the dynamics is linear, a wave packet will propagate in a straight line. Furthermore, since the nonlinear representation maintains the locality, but slightly shifts the position and alters the shape of the packet, the physical wave packet will also propagate approximately in a straight line. In particular, we conclude that a short wave that is being moved around by long waves will, on the average, move with constant velocity and will have constant wavenumber (with the proviso that there are no four-wave resonances present). We tested this conclusion with a ray-tracing calculation. Consider three big waves with wavenumbers $\pmb{k}_1$ (0.05, 0) rad/m, $\pmb{k}_2 = (0.04, 0.05)$ rad/m, and $\pmb{k}_3 = (0, 0.045)$ rad/m, each 

with a $1\mathrm{m}$ amplitude. (Three waves are the minimum needed in two dimensions to prevent the existence of a frame in which the flow is steady.) The ray equations, including intrinsic motion with effective gravity and with advection and straining by the big waves, were integrated for a small wave with initial conditions of $\pmb{q} = (20,0)$ rad/m, for a time of 2000 s, which is comparable with the viscous dissipation time of these waves. 

In order to compare to these ray-tracing results we generalize the effects of long waves on short waves that we discovered in one dimension. In order to do this we must first generalize the notion of the Hilbert transform of a function which was introduced in §4. This is accomplished by an inspection of (4.3). We consider the transform of a function $f(x)$ , which is defined in the Fourier-transformed wave-vector space, $f(\pmb{k}) \equiv \int \mathrm{d}x \mathrm{e}^{\mathrm{i}k \cdot x} f(x)$ , as 

$$
\boldsymbol {F} (\boldsymbol {k}) \equiv - \mathrm {i} \hat {\boldsymbol {k}} f (\boldsymbol {k}), \tag {7.1}
$$

where $\hat{\pmb{k}}$ is a unit vector in the $k$ -direction and is similar to (4.3). Our two-dimensional generalization consists of asserting that the horizontal Lagrangian displacement, which in one dimension is the Hilbert transform of the large-wave height, is given by the transform (7.1) acting on the large wave 

$$
S (\boldsymbol {x}) = - \sum_ {i} h _ {i} \hat {\boldsymbol {k}} _ {i} \sin \left(\boldsymbol {k} _ {i} \cdot \boldsymbol {x} - \omega_ {i} t + \phi_ {i}\right). \tag {7.2}
$$

where $h_i, k_i, \omega_i,$ and $\phi_i$ respectively represent the wave height, wavenumber, frequency, and initial phase for each of the large-wave modes. 

Our interpretation of the one-dimensional result (5.5) or (5.4) can now be extended (cf. the discussion following (5.5)) to two dimensions 

1. The centre of the packet is shifted by an amount $S(x_0)$ given by (7.2). 

2. There is a stretching and a rotation of the short-wave wavenumber given in component form by 

$$
q _ {i} ^ {\prime} = \left[ \delta_ {i j} - \nabla_ {i} S _ {j} \right] q _ {j} \tag {7.3}
$$

with an implied sum over the component index $j$ . 

3. There is a modulation of the width of the packet which is related to the transformation matrix in (7.3) and a corresponding change of its height given by the Jacobian of that transformation. 

For comparison to the ray-tracing calculation we then expect that the packet's wavenumber consists of a time-independent piece (i.e. conserved) and a variation of the form 

$$
\delta \boldsymbol {q} = \sum_ {i} h _ {i} \hat {\boldsymbol {k}} _ {i} \boldsymbol {k} _ {i} \cdot \boldsymbol {q} \cos \left(\boldsymbol {k} _ {i} \cdot \boldsymbol {x} - \omega_ {i} t + \phi_ {i}\right), \tag {7.4}
$$

which is just the gradient of (7.2). The sum is over the long-wave components. Figure 9(a) shows the two components of the deviation 

$$
\Delta \boldsymbol {q} = \boldsymbol {q} _ {\text {r a y}} - [ \delta \boldsymbol {q} + \boldsymbol {q} (t = 0) - \delta \boldsymbol {q} (t = 0) ]. \tag {7.5}
$$

On the average it is seen that $\Delta q_{x}$ is almost zero, while there appears to be an average shift in $\Delta q_{y}$ . 

The expectation for the position of the wave packet involves two effects, the horizontal Lagrangian displacement in (7.2) and a secular piece (not given by the 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-12/8237f025-efca-4750-807b-3f5bc5cb85c3/342c4d16de04a56feadac7705e37f5d9ce1d4d5dbb31e1b203813aac9a135989.jpg)



(w) $x\nabla$


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-12/8237f025-efca-4750-807b-3f5bc5cb85c3/4ae32d90425a45e493f18c878b96a3c24e032b5555d033bf811259ac32271993.jpg)



(w) A


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-12/8237f025-efca-4750-807b-3f5bc5cb85c3/cc9e188d431eaf57dc580bc51dfa7c2c4ad686bac91a0986e68ff0726d84aba2.jpg)



(1-w)²b


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-12/8237f025-efca-4750-807b-3f5bc5cb85c3/9b9a04f1c2c60ae0692cce530adf895777f39c394017b8257199dd7b44f6664c.jpg)



(1-u) $\pmb{b}\pmb{\nabla}$



Hnneepnnnnae



FIGURE 9. (a) The deviations (as a function of time) of the wavenumber of a wave packet (in the eikonal approximation) from that calculated by our improved representation: $\Delta q = q - \delta k$ . Here $\delta k$ represents the stretching and rotation of the wavenumber of the packet due to the interactions with the large waves. (b) The deviations of the position of a wave packet: $\Delta x = x - (S - x_{l})$ . Here $S$ is the Lagrangian displacement discussed in the text and $x_{l}$ is the distance due to Stokes' drift and the group velocity of the packet, $x_{l} = x_{l}$ ( $t = 0$ ) + $T(V_{\mathrm{stokes}} + V_{\mathrm{group}})$ .


linear theory) due to the difference between Lagrangian and Eulerian systems. The latter effect comes from the packet's group velocity 

$$
V _ {q} = \frac {1}{2} \hat {\boldsymbol {q}} (g / q) ^ {\frac {1}{2}} = (0. 3 6 3 6, - 0. 0 1 2 3) \mathrm {m} / \mathrm {s}. \tag {7.6}
$$

and the part of the Stokes-drift velocity that remains for linear, large waves (which is half the total) 

$$
V _ {\text {S t o k e s}} = \sum \frac {1}{2} k \alpha^ {2} \omega = (0. 0 3 3 3, 0. 0 3 4 7) \mathrm {m} / \mathrm {s}, \tag {7.7}
$$

being about $10\%$ of $V_{k}$ . Adding this drift velocity to the group velocity, we predict that the secular term grows at the rate 

$$
V _ {\text {T o t a l}} = (0. 3 9 6 9, 0. 0 2 2 4) \mathrm {m} / \mathrm {s}. \tag {7.8}
$$

In figure 9(b), we plot 

$$
\Delta \boldsymbol {x} = \boldsymbol {x} _ {\text {r a y}} - \left[ \boldsymbol {S} - t \boldsymbol {V} _ {\text {T o t a l}} - \boldsymbol {x} (t = 0) - \boldsymbol {S} (t = 0) \right]. \tag {7.9}
$$

Aside from the oscillating part and a $5\mathrm{m}$ part that would not occur if we have used the average linear $q_{y}$ , there are only $\approx 1\mathrm{m}$ deviations, out of the approximately $800\mathrm{m}$ travelled by the wave. The improved linear representation by itself would be wrong by $(V_{\mathrm{Stokes}})\times 2000\mathrm{s}$ , or about $100\mathrm{m}$ . Although this is about $10\%$ of the total distance, the important conclusion is that the motion of the small wave is predictable, and that there is no chaotic wandering of the wavenumber or position. We have verified that our representation in two dimensions extends the small-wave results presented in one dimension, and that the small-wave advection and straining by the large wave are correctly given. We have also included more long waves in the comparison and obtained similar results. 

Finally, we include the form of the height bispectrum takes in our two-dimensional generalization, 

$$
\begin{array}{l} \Lambda \left(\boldsymbol {k} _ {1}, \boldsymbol {k} _ {2}, \boldsymbol {k} _ {3}\right) = k _ {3} \left(\vec {\boldsymbol {k}} _ {1} \cdot \vec {\boldsymbol {k}} _ {3}\right) \left(\vec {\boldsymbol {k}} _ {2} \cdot \vec {\boldsymbol {k}} _ {3}\right) \Phi_ {0} \left(k _ {1}\right) \Phi_ {0} \left(k _ {2}\right) \\ + k _ {2} \left(\hat {\boldsymbol {k}} _ {1} \cdot \hat {\boldsymbol {k}} _ {2}\right) \left(\hat {\boldsymbol {k}} _ {3} \cdot \hat {\boldsymbol {k}} _ {2}\right) \Phi_ {0} \left(k _ {1}\right) \Phi_ {0} \left(k _ {3}\right) + k _ {1} \left(\hat {\boldsymbol {k}} _ {2} \cdot \hat {\boldsymbol {k}} _ {1}\right) \left(\hat {\boldsymbol {k}} _ {3} \cdot \hat {\boldsymbol {k}} _ {1}\right) \Phi_ {0} \left(k _ {2}\right) \Phi_ {0} \left(k _ {3}\right), \tag {7.10} \\ \end{array}
$$

where $\hat{k}_i$ are unit vectors in the appropriate direction. 

This work was supported by the Defense Advanced Research Projects Agency, and by La Jolla Institute Internal Research Funds. We also wish to thank M. Milder for his suggestion that our model would be useful for studying the interactions of long waves and short waves. 

# REFERENCES



BARRICK, D. E. & WEBER, B. L. 1976 On the nonlinear theory for gravity waves on the ocean's surface - Part II: Interpretation and applications. J. Phys. Oceanogr. 7, 11-21. 





CARY, J. R. 1981 Lie transform perturbation theory for Hamiltonian systems. Phys. Rep. 79, 129-159. 





Courant, R. & Hilbert, D. 1937 Methods of Mathematical Physics. Interscience. 





GARRETT, C. & SMITH, J. 1976 On the interactions between long and short waves. J. Phys. Oceanogr. 6, 925-930. 





GOLDSTEIN, H. 1980 Classical Mechanics. Addison-Wesley. 





HAsSELMann, K. 1984 The science and art of wave prediction - an ode to HO601. In A Celebration in Geophysics and Oceanography - 1982. Scripps Institution of Oceanography Reference Series 84-5, 1984). 





HENYEY, F. S., CREAMER, D. B., DYSTHE, K. B., SCHULT, R. L. & WRIGHT, J. A. 1988 The energy and action of small waves riding on large waves. J. Fluid Mech. 198, 443-462. 





KITAIGORODSKI, S. A. 1983 On the theory of the equilibrium range in the spectrum of wind-generated gravity waves. J. Phys. Oceanogr. 3, 816-827, 





LONGUET-HIGGINs, M. S. 1986 Eulerian and Lagrangian aspects of surface waves. J. Fluid Mech. 173, 683-707. 





LONGUET-HIGGINs, M. S. & STEwART, R. W. 1964 Radiations stresses in water waves; a physical discussion, with applications. Deep-Sea Res. 11, 529-562. 





MILES, J. Hamiltonian formulations for surface waves. Appl. Sci. Res. 37, 103-110. 





PHILLIPS, O. M. 1977 The Dynamics of the Upper Ocean. Cambridge University Press. 





PHILLIPS, O. M. 1985 Spectral and statistical properties of the equilibrium range in wind-generated gravity waves. J. Fluid Mech. 156, 505-531. 





SCHWARZ, L. W. 1974 Computer extension and analytic continuation of Stoke's expansion for gravity waves. J. Fluid Mech. 62, 553-578. 





ToBA, Y. 1973 Local balance in the air-sea boundary processes. III. On the spectrum of wind waves. J. Oceanogr. Soc. Japan 29, 1064-1068. 





WEST, B. J. 1981 Deep Water Gravity Waves, p. 33. Springer. 





WHITHAM, G. B. 1974 Linear and Nonlinear Waves. Wiley. 





ZAKHAROV, V. E. 1968 Stability of periodic waves of finite amplitude on the surface of a deep fluid. J. Appl. Mech. Tech. Phys. 9, 190-194. 

