![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/38aa9c0fe19667a590b9e32a29163a1a80fe36852bccd716395ec29a674de157.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/eed8b48f839e696c0e730d5254b74c17f315e4d6de10f67ab853ff591ea950fe.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/2f4bc32cb361998e286a1c9e519b1f7a3f88a006e89c4f19e29dda87d8f0d6bc.jpg)


# Rough seas – power laws and fractals

- aspects of the geometry of the ocean surface 

Paul H. Taylor 

University of Western Australia 

and Department of Engineering Science University of Oxford 

# What does the sea surface look like?

A simple question for waves on open ocean – assuming uni-directional random waves on deep water 

For linear models of waves $\omega ^ { 2 } = k g$ dispersion equation 

Standard oceanographic models for wave energy spectra : JONSWAP and Pierson-Moskowitz high freq. tail shapes 

$$
s (\omega) \sim \omega^ {- 5} \iff S (k) \sim k ^ {- 3}
$$

If the linear dispersion relationship holds out into the high tails 

$$
\omega^ {2} = k g
$$

period $\scriptstyle ( T = 2 \pi / \omega )$ 1s – 20s wavelength $\scriptstyle ( \lambda = 2 \pi / k )$ 2m – 800m 


Measured frequency power spectrum in Hurricane Camille


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/caa74ef8703dde34fc36f305023966aac29a68e2c7eef975da2c336bbeb5e9a0.jpg)


On open ocean there is a very large range of wavelengths ~ O(103) 

The motion of a FREE wave satisfies the dispersion equation. 

FREE wave components transport wave energy 

BOUND wave components exist because of the nonlinear free-surface boundary conditions and are ‘slaved’ to the free wave components. 

BOUND waves add to the TOTAL wave profile 

- making crests taller and spikier, and troughs less deep and more rounded. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/a845c37b7839095f5cfe27fe39b2cb603ae12ade3db88073afcf41fac471f4a9.jpg)


This work aims to explore the structure of the BOUND wave components 

# Free or bound ?

for free waves, 

a linear high frequency / wavenumber spectra tail can produce difficulties for 

mean square vertical acceleration of the surface 

mean square surface slope in space (also spatial skewness) 

For these are logarithmically singular and k−3 

$$
\overline {{\left(\frac {\partial^ {2} \eta}{\partial t ^ {2}}\right) ^ {2}}} = \int \omega^ {4} S (\omega) d \omega \sim \int \omega^ {- 1} d \omega \sim ?
$$

So 

ARE THE HIGH SPECTRAL TAILS REALLY FREE WAVES ? 

If not, what are they ? 

Or do we blame all difficulties on a short wave cut-off scale? 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/fdace31fc984c36ab25c7fe9adc8d770cc36beb81b8e5fdd0cdf23036b569770.jpg)


# ARE THESE ALL (or mostly) FREE WAVES ?

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/08a4582b815245f5ddfbb2924294bc160d79bc615c5e2fd55aa69a0e431ab749.jpg)


Photos of deep water storm waves in mid-Atlantic, by de Lange 

From ‘Heavy Weather Sailing’ $3 ^ { \mathsf { r d } }$ ed. by K.Adlard Coles 1981 

# ARE THESE ALL (or mostly) FREE WAVES ?

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/6b0d2367b53d2f4889d6d9c033e20815e030280fe92204e8226dd70416e93cc7.jpg)


Even in much less severe seas, some of the waves look very nonlinear. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/93f3e53612d8f88c599bdefdd52426ffb14ef9a5d3a6e2dde62ba1bc6c555c57.jpg)


At small scale in a mild sea, some of the crests look very nonlinear. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/6a7a551fe8458858f2a6d86edc59dd839b3fc7a095e7e961787a9615fc990242.jpg)


What do we know already: 

Linear waves – random, broad-banded – like real ocean ? 

Stokes waves – regular, fully non-linear – 19th century analysis 

2nd order Sharma & Dean – weakly non-linear, limited bandwidth, pairs of linear components 

Numerical solutions of full water wave equations 

∇2 _______________________ free − surface boundary conditions 

Linear 

Nonlinear 

# 2nd order Sharma and Dean (1979), based on Longuet-Higgins and Stewart

– weakly non-linear (2nd order only), limited bandwidth 

– short wave modulation by long wave 

(obtained via generalised Stokes perturbation expansion) 

2 main pieces of physics: 

Long waves accelerate Short waves vertically: variable gravity drives amplitude modulation 

$$
\updownarrow a _ {S M o d} = a _ {S} (1 + a _ {L} k _ {L} C o s [ k _ {L} x ])
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/c39eb5fdfcda9d652f6ca8d72c88b69d70d1afd50d934e3b18282e5c617ae29c.jpg)


Long waves compress and stretch Short waves horizontally: advection so wavelength modulation 

$$
\iff k _ {S M o d} = k _ {S} (1 - a _ {L} k _ {L} C o s [ k _ {L} x ])
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/3e2cc258c41aec6b997f53454afee855df24998e4d5376b0d6a4e120ee51a0fb.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/f4aad8d09c1a5f5638c68911de2c702d9450015ab47a5ba20aeff49b0f114e94.jpg)


For a regular wave, Sharma & Dean $=$ Stokes 2nd order 

# 2nd order Sharma & Dean – weakly non-linear, limited bandwidth

# Creamer to 2nd order can be reduced to Sharma & Dean deep water

Taylor PH (1992) On the kinematics of large ocean waves, Proc. of Behaviour of Offshore Structures Conference, Imperial College, BOSS92, Vol.1, 134. 

# Actually Creamer is better for very long – very short wave pairs

SD has problems with horiz. displacement of short waves by linear motion of long waves 

Creamer D.B., Henyey F., Schult R. and Wright J. (1989) 

J. Fluid Mech. 205, 135-161. 

Improved linear representation of ocean surface waves 

$$
\eta_ {k} = \frac {1}{A b s [ k ]} \int_ {- \infty} ^ {\infty} d x E x p [ - i k x ] (E x p [ i k \eta_ {L H} (x) ] - 1)
$$

where 

is linear horiz. displacement of fluid particles at the surface, $\eta _ { k }$ is the FT of the non-linear vertical surface displacement 

Actually it is based on a Lie-group canonical transformation (change of variables) 

New variables chosen to give best possible linear representation 

- such that inverse transformation captures as much of the (non-resonant) bound wave structure as possible 

# Apply Creamer to a single regular wave

– how much Stokes harmonic structure is generated ? 

Creamer 5th compared to Stokes 5th and linear for AK=0.4 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/14dc1697516995c31f210fc22880d64e6eb097234ee0d00e12fcd48f72669f3e.jpg)


So although Creamer isn’t exact, it does MUCH better than linear 

- exact at $2 ^ { \mathsf { n d } }$ order, ‘best’ approx. at all higher orders 

Exact Stokes regular wave theory 

$$
\begin{array}{l} \eta_ {\text {L i n e a r}} = \varepsilon \cos \theta \\ = \varepsilon \cos (k x - \omega t) \\ \end{array}
$$

$$
\begin{array}{l} \eta^ {\mathrm {S T O K E S}} = \left(\varepsilon - \frac {3}{8} \varepsilon^ {3} - \frac {2 1 1}{1 9 2} \varepsilon^ {5}\right) \cos \theta + \left(\frac {1}{2} \varepsilon^ {2} + \frac {1}{3} \varepsilon^ {4}\right) \cos 2 \theta \\ + \left(\frac {3}{8} \varepsilon^ {3} + \frac {9 9}{1 2 8} \varepsilon^ {5}\right) \cos 3 \theta + \left(\frac {1}{3} \varepsilon^ {4}\right) \cos 4 \theta + \left(\frac {1 2 5}{3 8 4} \varepsilon^ {5}\right) \cos 5 \theta + \dots \quad \omega^ {2} = k g \left(1 + \varepsilon^ {2} k ^ {2} + \dots\right) \\ \end{array}
$$

Exact Stokes theory for regular waves (here to $5 ^ { \mathrm { { t h } } }$ order) 

# Apply Creamer to a single regular wave

$$
\eta_ {\text {L i n e a r}} = \varepsilon \cos \theta
$$

how much exact Stokes harmonic structure is generated ? 

$$
= \varepsilon \cos (k x - \omega t)
$$

$$
\begin{array}{l} \eta^ {\text {S T O K E S}} = \left(\varepsilon - \frac {3}{8} \varepsilon^ {3} - \frac {2 1 1}{1 9 2} \varepsilon^ {5}\right) \cos \theta + \left(\frac {1}{2} \varepsilon^ {2} + \frac {1}{3} \varepsilon^ {4}\right) \cos 2 \\ + \left(\frac {3}{8} \varepsilon^ {3} + \frac {9 9}{1 2 8} \varepsilon^ {5}\right) \cos 3 \theta + \left(\frac {1}{3} \varepsilon^ {4}\right) \cos 4 \theta + \left(\frac {1 2 5}{3 8 4} \varepsilon^ {5}\right) \cos 5 \theta + \dots \quad \omega^ {2} = k g \left(1 + \varepsilon^ {2} k ^ {2} + \dots\right) \\ \end{array}
$$

$$
\begin{array}{l} \eta^ {C R E A M E R} = \left(\varepsilon - \frac {3}{8} \varepsilon^ {3} + \frac {1}{6} \varepsilon^ {5}\right) \cos \theta + \left(\frac {1}{2} \varepsilon^ {2} - \frac {5}{1 2} \varepsilon^ {4}\right) \cos 2 \theta \\ + \left(\frac {3}{8} \varepsilon^ {3} - \frac {6 3}{1 2 8} \varepsilon^ {5}\right) \cos 3 \theta + \left(\frac {1}{3} \varepsilon^ {4}\right) \cos 4 \theta + \left(\frac {1 2 5}{3 8 4} \varepsilon^ {5}\right) \cos 5 \theta + \dots \\ \end{array}
$$

# CREAMER APPROX

$$
\omega^ {2} = k g
$$

Differences in red, Creamer is correct for leading term at each order 

# Surface geometry: regular waves

Creamer approximation for n-th harmonic: 

$$
\eta_ {n} = \left(\frac {1}{n !}\right) \left(\frac {n}{2}\right) ^ {(n - 1)} A ^ {n} k ^ {(n - 1)} \cos n (k x - \omega t) + O (A ^ {n + 2} k ^ {(n + 1)})
$$

$$
N = 1 \quad : \quad 1
$$

$$
N = 2 \quad : \quad 1 / 2
$$

$$
N = 3 \quad : \quad 3 / 8
$$

$$
N = 1 0: 7 1 8 5 / 1 4 5 1 5 2 = 0. 5 3 8 2..
$$

$$
N = 1 0 0: 1. 6 9 0 5 5 \times 1 0 ^ {1 0} \quad \text {w h i c h n o n e o f u s k e w !}
$$

[Note the convergence limit AK < 2/e = 0.73] 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/35fd69a650b2c5c369d2a45256334f83923ef723f602f0ea4e23f668581f51fe.jpg)


# BUT DYNAMICS ARE LINEAR

Stokes – nonlinear dispersion 

$$
\omega^ {2} = k g \left(1 + \underline {{A ^ {2} k ^ {2}}} + \dots\right)
$$

Creamer – linear 

$$
\omega^ {2} = k g
$$

Creamer surface is correct at leading order for all harmonics differences at next term - 2 orders higher 

Lie-group canonical transformation (change of variables) 

Gives best possible linear representation 

– capturing much of the non-resonant bound wave structure 

Creamer transform : 

Linear geometry along a line 

è (Approx. for) Non-linear geometry along a line 

It is a spatial transformation – creating (approx for) the bound harmonics 

- what about random, broad-banded and steep SEA-STATES ? 

# Creamer et al. give

# (Nonlinear) Wavenumber spectrum $S ( k )$ for vertical surface displacement in a random sea-state

$$
S _ {N L} (k) = \frac {1}{k ^ {2}} \int_ {- \infty} ^ {\infty} d x C o s [ k x ] E x p [ - k ^ {2} \sigma^ {2} ] \left(E x p [ k ^ {2} \sigma^ {2} Y _ {0} (x) ] - 1\right)
$$

Surface elevation variance Auto-correlation fnLinear inputs: 

Auto-correlation function $Y _ { o } ( x )$ - related to linear spectrum 

$$
\left\langle \eta_ {L H} \left(x _ {o}\right). \eta_ {L H} \left(x _ {o} + x\right) \right\rangle = \left\langle \eta_ {L} \left(x _ {o}\right). \eta_ {L} \left(x _ {o} + x\right) \right\rangle = \sigma^ {2} Y _ {0} (x)
$$

where $\sigma$ is the standard $\mathsf { d e v } = H _ { s } / 4$ 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/b8daa0c5a16c07151f97592220349adf49ea6224ec5b0518ca2c49c83523f759.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/9647e924dbe424105a71506a239194bd15ec530034b30e0517fc3355352c9c11.jpg)


JONSWAP wavenumber spectrum 

- realistic spectral model 

Auto-correlation $Y _ { o }$ in space also the NewWave shape 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/7c394d2540918b45eba7fb9ee189770a1a240cbbc391397231cf2d43624ddee4.jpg)


Linear random surface elevation profiles 

NSpectral tails $S _ { L I N } { \sim } k ^ { - N }$ with N = 2, 3, 4,,5 with range k = 1,100 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/55ef1d3d1d31a6c7c08a60ef4ec38ec9b96fae63ba9d6b6f3edc60f3916dcdbd.jpg)


These are FRACTAL surfaces - roughness looks the same at all scales - roughness looks the same at all scales 

- as discussed by Mandelbrot (1977) 

What type of linear surface roughness scaling is consistent with water wave hydrodynamics ? 

Surface roughness – input wavenumber spectrum $S _ { L I N } ( k )$ (linear approx.) 

Surface displacement – variance $\sigma ^ { 2 } = \int _ { 0 } ^ { \infty } \ S _ { L I N } ( k ) d k = \overline { { { \eta } ^ { 2 } } }$ 

Zero-crossing wavenumber $k _ { z } = \sqrt { \int _ { 0 } ^ { \infty } k ^ { 2 } S _ { L I N } ( k ) d k \biggl / \int _ { 0 } ^ { \infty } S _ { L I N } ( k ) d k } = \sqrt { \overline { { \eta _ { x } ^ { 2 } } } \biggl / }$ η 2 

Spatial auto-correlation function $Y _ { o } ( x ) = \int _ { 0 } ^ { \infty } C o s ( k x ) S _ { _ { L I N } } ( k ) d k \bigg / \int _ { 0 } ^ { \infty } S _ { _ { L I N } } ( k ) d k$ 

Surface roughness – Creamer output wavenumber spectrum SNL (k) (close to fully nonlinear) 

# What types of spectra are relevant ? Choosing  SLIN (k)

# Linear spectrum

– main components moving according to (close to) linear dispersion 

– phase speed at (close to) linear phase speed, energy at group velocity 

– main upper tail form is generally chosen to be hyperbolic 

$$
S _ {L I N} (k) \sim k ^ {- n}
$$

Upper and lower wavenumber cut-offs ? 

– possible relevance of upper cut-off for structure of whole spectral form, and maybe implications for breaking 

Spectral form only enters Creamer spectrum calculation via $Y _ { 0 } ( x )$ 

So, what is the possible structure of surface spatial auto-correlation? 

Spectral form only enters Creamer spectrum calculation via $Y _ { 0 } ( x )$ 

So, what is the possible structure of surface spatial auto-correlation? 

$$
Y _ {O} (x) = \frac {\int_ {0} ^ {\infty} \cos (k x) S _ {L I N} (k) d k}{\int_ {0} ^ {\infty} S _ {L I N} (k) d k} = 1 - \frac {1}{2} k _ {z} ^ {2} x ^ {2} + \dots .
$$

All linear spectra with an upper cut-off or with a tail decaying faster than $\cdot$ -3 have this local parabolic form 

Y0 (x)___________________________________________________: What is the possible structure of the surface spatial auto-correlation? 


Linear spectrum - hyperbolic


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/6ac122ffdcfabed4974494ac065052da386a3a2f4d3ef5d71351b8b8e86c2345.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/1877e1e97b02727c3fdd2fc1cb03d7821235e9ef77e7930459b2510f8b12b64d.jpg)


$$
Y _ {0} (x) = 1 - \frac {1}{2} k _ {p} ^ {2} x ^ {2} \left(\frac {n - 1}{n - 3}\right) \left(\frac {1 - m ^ {3 - n}}{1 - m ^ {1 - n}}\right) + \dots
$$

$$
Y _ {0} (x) = 1 - \frac {1}{2} k _ {p} ^ {2} x ^ {2} \left(\frac {n + 1}{n - 3}\right) + \dots
$$

Both auto-correlation forms locally parabolic. 

All linear spectra with an upper cut-off or with a tail decaying faster than $\cdot$ -3 have a local parabolic form  Y (x) = 1− 1 $Y _ { 0 } ( x ) = 1 - \frac { 1 } { 2 } k _ { z } ^ { 2 } x ^ { 2 } + . . . . .$ k 2 x 

# Numerical solution for full Creamer nonlinear wavenumber spectrum MATHEMATICA numerical integration

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/9853d1a3038b0cd1fa68647ecf94aff472598d721113986340eb64df64f65f5d.jpg)



Bound difference terms Long waves – set-up and set-down



Bound sum harmonics – locked in phase with linear components


# Set-down /set-up - structure of spectral limit for $\sigma ^ { 2 } k ^ { 2 }$ small $( k \to 0 )$

$$
S _ {N L} (k) \rightarrow \frac {1}{2} k ^ {2} \sigma^ {4} \int_ {- \infty} ^ {\infty} Y _ {0} ^ {2} (x) d x
$$

General result for low tail of wavenumber spectrum 

In terms of the autocorrelation function 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/e8755144431ccd18228d5a2476927cb5739c0a9e88861e0c13b28f51d63bfb8b.jpg)


Are the low tails free waves ? 

LOW TAIL IS ALWAYS k2 

BOUND WAVES 

Depends on overall sea-state rms height 

# Numerical solution for full Creamer nonlinear wavenumber spectrum MATHEMATICA numerical integration

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/2fbb7dab2a7c4b0dd001b3de4ab591b2d6a5fdae7ab339d629584358be13e3a6.jpg)



Bound difference terms Long waves – set-up and set-down



Bound sum harmonics – locked in phase with linear components


Spectral saturation - structure of spectral tail for large2 2 σ k (k → ∞) 

$$
\overline {{S _ {N L} (k) \rightarrow \frac {1}{k ^ {3}} \frac {\operatorname {E x p} \left[ - 1 / \left(2 \sigma^ {2} k _ {z} ^ {2}\right)\right]}{\sigma k _ {z}}}}
$$

Limiting form for high tail of wavenumber spectrum 

- UNIVERSAL result 

Requires $Y _ { 0 } ( x ) \to 1 - \frac { 1 } { 2 } k _ { z } ^ { 2 } x ^ { 2 } + . . . . .$ 

- auto-correlation is locally smooth 

Are the high tails free waves ? 

FAR TAIL IS ALWAYS 3k − 

BOUND WAVES 

Size depends on overall sea-state steepness 

Spectral saturation - structure of spectral tail for large 2 2 σ k (k → ∞) 

$$
S _ {N L} (k) \rightarrow \frac {1}{k ^ {3}} \frac {\operatorname {E x p} \left[ - \left(1 / k _ {z} ^ {2}\right) / \left(2 \sigma^ {2}\right)\right]}{\sigma k _ {z}}
$$

Note appearance of Rayleigh form with x=1/kz 

Exceedance rate for Rayleigh distribution for linear crest elevation 

$$
P (X > x, \sigma) = E x p \Bigl [ - x ^ {2} / (2 \sigma^ {2}) \Bigr ]
$$

This result doesn’t give a limiting crest height but it does link 

the size of the $k ^ { 3 }$ high spectral tail 

to the overall steepness of the sea-state (measured as $\sigma k _ { z } )$ , 

via the probability of an individual crest with an amplitude greater than $1 / k _ { z }$ in a sea-state with an RMS roughness of σ 

SPECTRAL SATURATION IN CREAMER HIGH TAIL 

- strongly amplitude dependent 

$$
\begin{array}{r l} \vdots \\ S _ {N L} (k) & \to \frac {1}{k ^ {3}} \frac {E x p [ - \frac {1 / k _ {z} ^ {2}}{2 \sigma^ {2}} ]}{\sigma k _ {z}} \end{array}
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/3451d12526112fed25e3fb3a03a5ce4aa44b2b545a814d6aa8bac1049ffbd5f2.jpg)


Exponentially small $k ^ { 3 }$ spectral tail for small sea-surface roughness 

$$
\sigma k _ {z} <   \sim 0. 2
$$

Same range as from low order spectral estimates 

Spectral tail depends on global sea-surface roughness only: $\sigma k _ { z }$ 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/63c374d0aea8da09c1924e8c214da539500ab8957bc5d62778c230dd769dbc51.jpg)


So, what do Creamer limiting water waves look like? 

CUSPS occur at random on the surface, ALWAYS PRESENT, 

but number and size depends on wave energy density, each with local wavenumber amplitude spectrum $A ( k ) \sim k ^ { - 3 / 2 }$ 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/ca8d38b8a78de2504188aff3bec1e22da168443ff348e60beeadb1d03e310394.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/de43a233aaece58de191782bbc525e7897af9fc63d5ee864c6b294c843782c2d.jpg)


So what do Creamer limiting water waves look like? 

Cusps occur at random on the surface, 

ALWAYS PRESENT, 

but number and size depending on wave energy density 

with local wavenumber amplitude spectrum $A ( k ) \sim k ^ { - 3 / 2 }$ 

If sea-state is very low, cusps are very small and very rare 

$$
\Rightarrow \sigma k _ {z} <   \sim 0. 2
$$

Probably the occurrence of cusps limits overall surface roughness 

- consistent limit from low order harmonics and high spectral tail 

And what is observed on the open sea 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/455b905ffd43ccbc7969cf84c5d725c59ddc8dc8a66ad92b2502e3626ddf3199.jpg)


# Weak or strong nonlinearity: the vital is:

R. C. T. Rainey 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/21c7174c6b795116b03221db38b19f370f8365bb16de3c08d4dcc9b3ce4e305c.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/9a703c4c48a0c4ad36a92d5729c1e3c2b5d7c263c7bc28e16452b9777a91184d.jpg)



Fig.1 The upper graphs are successive positions of a sheet of particles initially on the zero-pressure surface,for two 1st order waves of the same steepness $k a = 0 . 1 8$ , one twice the length of the other. The lower graphs are the corresponding conventionally defined surfaces (1.1). One wavelength of the longer wave $( 2 \pi )$ is illustrated, and the vertical and horizontal scales are thesame.Theboldlineis theinitial position, when the crestof thelonger wave,and thetroughofthe shorter wave, are both at $\pi$ . Subsequent lines show the surface at successive time-steps equal to $5 \%$ of the time until the conventionally defined crests coincide. The arrow shows where this crest coincidence occurs


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/402ef3b738c0778756d0588c0b1daf9a4bc14466525ac23ed254e5961db327b3.jpg)


An example of a cusp? 

An example of the crest of an extreme standing wave in a tank 

The white circle in the image is a table-tennis ball attracted to the crest 

– this remarkable image is from Prof Stephen Salter of Edinburgh Univ. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/5f407579d07b4d5cdfea6fbec19d836ee92fde7e69e40b68928016698b33f66b.jpg)


An example of a cusp? 

# Full Creamer nonlinear wavenumber spectrum - MATHEMATICA numerical integration

# - closed form low and high wavenumber asymptotes

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/f4e93e2ea7d0c157cdcb4f61a5de6745991c215bf15e27c4058ba2caaac1f1d7.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/201c6542498576bf9107d989cf8309398cf5590427c258f6831cea32db7cbb13.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/87b682d28a9cd17430739c33d4c5514eb2bed1a4241f6e1f9b411e95370301ed.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/4569efd5f910adf9e34d1c36b2a7ed00191be210903cb1445147665c23c9446e.jpg)


Fractals – and ocean waves ? 

Random surface elevation profiles with broad-banded linear $S _ { L I N } ( k ) \sim k ^ { - N }$ 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/b7b0822d32357229973f0275fb6625e825f889a453187cc8460c2173a34bbe35.jpg)


N = 2, <3 NOT compatible 

with Creamer transform 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/23d82cede960795d131ba4c89d17f00de1bef12b86586742196377d5d558bd87.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/b6c3136ee132d2946a1963fa13a6c8c3963561406de0d14f20e20823c6a6387d.jpg)


N >3, 4,,5 compatible up to limiting sea-state steepness ? 

Linear FRACTAL surfaces may be relevant to water wave hydrodynamics - linear part of the spectrum associated with dynamics, transporting energy etc. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/fbce1a408063d3e462a2fe74211b5b3122bd7a8d36a964b9c1bae4e52eb3b854.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/bb3a4234283731bc2d0f4be954128a1a26d6e0039970052b7e318ca44020bd9e.jpg)


# A use for the Creamer transform

Linear initial conditions large error waves 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/42a2e42d4ae47d1b0072a52182bf1d2c19be23e1eacf03c2560719746c625133.jpg)


Creamer-based initial conditions 

~ no error waves 

Wavenumber spectrum via Stokes expansion of Creamer result: 

$$
\frac {1}{k ^ {2}} E x p [ - k ^ {2} \sigma^ {2} ] \left(E x p [ k ^ {2} \sigma^ {2} Y _ {0} (x) ] - 1\right)
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/0f79c12fa58502d7544264f4d6996e16d6ab1bf34be74ce421670adcd8706ce8.jpg)


$$
= \sigma^ {2} Y _ {0} + k ^ {2} \sigma^ {4} \left(- Y _ {0} + \frac {1}{2} Y _ {0} ^ {2}\right) + k ^ {4} \sigma^ {6} \left(\frac {1}{2} Y _ {0} - \frac {1}{2} Y _ {0} ^ {2} + \frac {1}{6} Y _ {0} ^ {3}\right) + \dots
$$

Wave profile contribution 

Linear Quadratic (2nd order) 

Cubic (3rd order) 

If linear spectrum is $S ( k ) \ \sim \ k ^ { n }$ with $n > 3$ 

Leading order expansion of quadratic, cubic,… terms 

ALL decay as $\sim k ^ { n }$ 

$$
\begin{array}{l} S _ {N L} (k) \rightarrow \sigma^ {2} (k / k _ {z}) ^ {- n} \left(1 + k _ {z} ^ {2} \sigma^ {2} (\dots) + k _ {z} ^ {4} \sigma^ {4} (\dots) + k _ {z} ^ {6} \sigma^ {6} (\dots) + \dots\right) \\ f o r k / k _ {z} \rightarrow \infty \\ \end{array}
$$

There is a remarkable degree of internal cancellation going on as naïve scaling suggests that each order should be $k ^ { 2 }$ higher 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/0d053d78d7aa2bdbbc1b09cb61c6ccd6b646c569b78d1e058d1c175a2694de75.jpg)


# Analytical solutions for Stokes-type expansion of Creamer nonlinear wavenumber power spectrum

Linear to Linear+2nd+3rd+4th order 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/7d2561ed2d445726d2da4a8c7e1267f7c36ae15bfa1eb079cd46a1d7d69e9740.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/4cbb659f19d81ea0d6eba0dfaa17528f0e0a24de77fdedb29be148834470e691.jpg)



Moderate steepness $\sigma k _ { p } = 0 . 1 7$


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/6c65cf10154111a9053bfd865f2e180b22876ac30af9b6f169c0e9dbde05f5ee.jpg)



Too steep $\sigma k _ { p } = 0 . 2 8$


So: if linear spectrum is $k ^ { n }$ with $n > 3$ 

Leading order expansion of quadratic, cubic,… terms ALL decay as k-n 

but the (partial) series is still badly behaved if sea is steep enough 

(? relevant for the far tail) 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/3d3663eb651b3d3a759a048a91a8b978cc42321dcd09a629d95b3ebe45e4ab8b.jpg)



SNL(k)/SLIN(kp)


$$
\sigma k _ {p} = 0. 3 3, 0. 3 8, 0. 4 3
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/9e826b1bbc5e84a97be13e22762aa14f383b6086c26209aaaa105cc2191e1ff7.jpg)


Wavenumber spectrum $S _ { N L } ( k )$ 

MUST BE +ve definite 

so this might give a method to estimate limiting sea-state steepness 

$$
\text {L i n e a r} + 2 ^ {\mathrm {n d}} \text {o r d e r :} \sigma k _ {p} <   0. 6 6
$$

Linear+2nd+3rd order : $\sigma k _ { p } < 0 . 3 8$ 

Linear +…+ 11th order: $\sigma k _ { p } < 0 . 1 3$ 

What happens as order → ∞ 

$$
S _ {L I N} (k) \sim k ^ {4} / \left(k _ {p} ^ {2} + k ^ {2}\right) ^ {4}
$$

Fits by Pade approximants 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/eba408aa87000477460ac0399d91326fa6a19b88301cd3c350c4a10417166cd2.jpg)


# ESTIMATING WAVENUMBER HARMONIC GIVING LIMITING STEEPNESS

INPUT LINEAR SPECTRUM 

$$
S _ {L I N} (k) \sim k ^ {4} / \left(k _ {p} ^ {2} + k ^ {2}\right) ^ {4}
$$

Fits by Pade approximants 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/051d5758e0c31e175715982133906025ad15da11c6053c186a81b95f13208e52.jpg)


So, what do Creamer limiting water waves look like? 

CUSPS occur at random on the surface 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/df8031319958d0f81338a2bdb30ac560af1e29c4ea1cf0791c522802c5f6f7c1.jpg)


Comments 

Also look like Rainey’s particle escape trajectories for his linear model for breaking 

So, what do Creamer limiting water waves look like? 

# Cusps :

A: Small but always catastrophic or 

B: Large and dominant 

# Option A: SMALL LIMIT ?

If sea-state is low, cusps are (very) small and (very) rare 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/c13e041d44afcfa863a6f351538fd3b56410baff94145729c298fd878f42c1b9.jpg)


Possibly the occurrence of cusps limits overall surface roughness 

Option B: LARGE LIMIT ? 

If sea-state is very steep (rough), cusps would be common and tall 

# What next

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/e983cddf2ef043b119f9b43b59344eb9b2fbec3eb6f80c8154c35358b83fd4fa.jpg)


We have a model for the tail of wavenumber spectrum S(k) 

• What about tail of frequency spectrum s(ω) ? 

• What about phase speed and frequency dispersion in the tail ? 

No idea how to work out either… 

but I think that all are BOUND so short components move with the long linear crests 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/21135e73d657d6e5fb93363682b48a252445827ac546f8bbb6e6854dee466d69.jpg)


# CONCLUSIONS:

In 1989 Creamer, Henyey, Schult and Wright published a canonical transform for 1-D waves on deep water 

• This gives an excellent approximation for ALL bound waves 

• Linear tail $S _ { L I N } ( k ) \sim k ^ { - N }$ and $N { < } 3$ inconsistent with hydrodynamic 

• Nonlinear high-tail always $\sim k ^ { - 3 }$ 

- BOUND WAVES via long-short interactions 

- UNIVERSAL form, limited and set by overall roughness 

- Spectral form implies cusps in the surface geometry 

• High wavenumber spectral tail limited by saturation of small-scale cusps 

• Maybe whole sea-state steepness limited by appearance of cusp 

Questions ? 

Comments ? 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/da3b0d63bfa0235a8725fd95436058c6d77b32acbf1c6aadac3215bfbb9ddc77.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/1a91d463bcda8421f3abb887d1c7fcce6a1d25cc3fc6d73e0372ff9dbcec7967.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/10f35071e4c8c7d37205e39fbb4e025b9427f89998746a8e749c64dadee213e5.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/06f8542b0c2f83ea9f502a16d32b4c5bc58aae7964d391a7332b2dd83c6b982c.jpg)



Input spectrum – Gaussian


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/899a1a0fa13a312adca35cf4dc5dbf26e09285d1afcfc4e7dfd0ab1efa4bfdd4.jpg)



Auto-correlation


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/ac3bd53d103f8e3818590c3b3793bb97342b91a64fbce127c16e43a0488cadd1.jpg)



Creamer transformed spectrum



(to 4th order)


$$
\begin{array}{l} \mathrm {T R U N C A T E D J O N S W A P S P E C T R U M} \quad k _ {\mathrm {p}} \mathrm {H s} / 4 = 0. 1 5 \\ + \text {C r e a m e r} \\ \end{array}
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/87bdb5068c657c4adfadf83fd2406b5515a9994b5cdfa01eefb1eeb3c9d6d2fa.jpg)


Why the interest – because the non-linear tail is starting to grow for large k 

# PHYSICAL LIMITS

Regular waves – maximum tail for 

$$
\sigma k _ {p} = 1
$$

0.31 

Power-law tail $n { > } 3$ and $m { > } { > } l$ 

$$
\sigma k _ {p} = \sqrt {\frac {n - 3}{n - 1}}
$$

? 

Power-tail tail n=3 – real ocean ? 

$$
\sigma k _ {p} = \sqrt {\frac {m ^ {2} - 1}{2 m ^ {2} \operatorname {L o g} [ m ]}}
$$

? SMALL FOR REAL OCEAN ~0.1？ 

$$
\rightarrow \sqrt {\frac {1}{2 \operatorname {L o g} [ m ]}}
$$

$$
f o r m > > 1
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/cd87df3a03cfa881692b3648c272ba1c25fa905632f328fd875b923e4aa1dbd0.jpg)


Maximum possible spectral tail occurs for _____________________________ ${ \sqrt { a } } \ \sigma \ k _ { p } = 1 { \mathord { \left/ { \vphantom { \sqrt { 2 } } } \right. \kern - delimiterspace } { \sqrt { 2 } } }$ 

taking $S ( k ) \sim k ^ { - n }$ 

Regular waves – maximum tail for $\sigma \ k _ { p } = 1$ 

Power-law tail n>3 and m>>1  n −1 $n { > } 3$ 

Power-tail tail $n { = } 3$ 2 m 2 Log[m] 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/d2d2af5c8ed15bcbd8156ec428dfcf9871fc07a8df4b92cd86c1e71475459144.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/699cbf2e1edee9aa709fb442317c0ac43abeb10fe98a22b31fdc787abd2eadef.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/52afad027eabe5a8a0c83ec4b31a55c81608193b0db611546d959052dc473ff1.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/47837ae2cba10054aae04d85c6f493ee9c8855bac91a98dc8988cbead51177e2.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/a6b0b463494207ac394bd025d0b32e138a7a3d37b74be3fad14147fc56e4261a.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/691a3bc79df5ef545ea7dc1403ed0385b4c1ce2c89372fc74845b7d175a7354e.jpg)



Stormy seas off the Ayrshire coast on Wednesday. Photo by Linda Rayner from Troon.



http://www.bbc.co.uk/news/uk-scotland-30834715


# Regular waves

$$
\eta_ {\text {L i n e a r}} = \varepsilon \cos \theta = \varepsilon \cos (k x - \omega t)
$$

– classical Stokes expansion for waves on deep water (1847) 

$$
\begin{array}{l} \eta^ {\mathrm {S T O K E S}} = \left(\varepsilon - \frac {3}{8} \varepsilon^ {3} - \frac {2 1 1}{1 9 2} \varepsilon^ {5}\right) \cos \theta + \left(\frac {1}{2} \varepsilon^ {2} + \frac {1}{3} \varepsilon^ {4}\right) \cos 2 \theta \\ + \left(\frac {3}{8} \varepsilon^ {3} + \frac {9 9}{1 2 8} \varepsilon^ {5}\right) \cos 3 \theta + \left(\frac {1}{3} \varepsilon^ {4}\right) \cos 4 \theta + \left(\frac {1 2 5}{3 8 4} \varepsilon^ {5}\right) \cos 5 \theta + \dots \\ \end{array}
$$

$$
A n d \qquad \omega^ {2} = k g \left(1 + \varepsilon^ {2} k ^ {2} + \dots\right)
$$

Nonlinear corrections to the surface shape AND to the dispersion equation 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/68165ab3c2d29c1337f5ee251c753ec49f4ae8b2e42b6bd97e367a4b55e16a86.jpg)



Input spectrum – Gaussian


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/d13007ea6fff0fc993fa114ca180c106655747aa6f5e816832a96fa4062da877.jpg)



Auto-correlation


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-13/3ffdc0b7-87b6-49ab-a233-8306a2c44b13/53ec0e489d34f984f0503ad15a1feb6c30800be00cb4c9690eb0200303608dce.jpg)



Creamer transformed spectrum



(to 4th order)
