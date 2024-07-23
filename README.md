# BFSS Model on the Lattice

## The Model

The BFSS model is said to be dual to type II string theory. It is obtained from the reduction of the 9+1 dimensional Supersymmetric Yang-Mills to 0+1 dimensions. 
The complete action of the BFSS model is therefore given by 

$$ S_{M}={N}\int_{0}^{\beta}dt~\mathrm{Tr}\left\lbrace\frac{1}{2}\big(D_{t}X_{M}\big)^{2}-\frac{1}{4}[X_{M},X_{N}]^{2} + i{\bar{\psi}}\gamma^{10}D_{t}\psi\ -{\bar{\psi}}\gamma^{M}[X_{M},\psi]\right\rbrace $$

where $\gamma^M~(M=1, \cdots, 10)$ are the $16\times16$ upper right block matrices of the 10-d Gamma matrices $\Gamma^M$.

Therefore the model is the quantum mechanics of $N\times N$ matrices $X$ and $\psi$ (where $N$ is the number of colors).

<!-- On lattice regularization, the model action becomes  -->

<!-- $$ S_M = S_b + S_f $$ -->

<!-- with -->

<!-- $$ S_{b}\ =\displaystyle \frac{N}{2a}\sum_{t.M}\mathrm{Tr}\left(U X_{M}(t+a)U^{\dagger}-X_{M}(t)\right)^{2}-\frac{N a}{4}\sum_{t.M.N}\mathrm{Tr}[X_{M}(t),X_{N}(t)]^{2} $$ -->

<!-- ![equation](https://latex.codecogs.com/svg.image?&space;S_{b}\=\displaystyle\frac{N}{2a}\sum_{t.M}\mathrm{Tr}\left(U&space;X_{M}(t&plus;a)U^{\dagger}-X_{M}(t)\right)^{2}-\frac{N&space;a}{4}\sum_{t.M.N}\mathrm{Tr}[X_{M}(t),X_{N}(t)]^{2}) -->

<!-- $$\displaystyle S_{f}=i N\sum_{t}\mathrm{Tr}\bar{\psi}(t)\left(\begin{array}{c c}{{0}}&{{D_{+}}}\\ {{D_{-}}}&{{0}}\end{array}\right)\psi(t)-a N\sum_{t,M}\bar{\psi}(t)\gamma^{M}[X_{M}(t),\psi(t)] $$ -->

## The Code
The project was to grasp the essence of lattice field theory, build on the C++ implementation of the BFSS model by Dr. Bergner and write and analyze observables in the simulation. The specific observables written by me are:

- Bosonic Energy: [bfssconfig.h  Line 517](/MCSC-CPPCODE/src/bfssconfig.h#L517)
- Fermionic Energy: [fermionmeasurements.h Line 36](/MCSC-CPPCODE/src/fermionmeasurements.h#L36)
- Gauge Invariant 4-point correlator (not normalised), $\int dt\left \langle ~\mathrm{Tr}(X(t)X(t)) ~ \mathrm{Tr}(X(t + \Delta t)X(t + \Delta t)) ~\right \rangle$: [bfssconfig.h Line 560](/MCSC-CPPCODE/src/bfssconfig.h#L560)

## The Data

Unfortunately, I have lost access to the simulation data for Bosonic and Fermionic energies. Here I present the analysis for the Gauge Invariant 4-point correlator.

(The complete data and the jupyter notebooks of the analysis are in the folder DATA_ANALYSIS)

### [DIFFTEMPS](/DATA_ANALYSIS/CPP/DIFFTEMPS/)

Correlators vs $\Delta t$ for different temperatures, number of lattice sites = 32

Temperatures 0.10 - 0.18   |  Temperatures 0.20 - 0.28
:-------------------------:|:-------------------------:
![](/DATA_ANALYSIS/CPP/DIFFTEMPS/temps_0.18-0.10.png)  |  ![](/DATA_ANALYSIS/CPP/DIFFTEMPS/temps_0.28-0.20.png)

Note that the correlator is symmetric about $\Delta t = 32/2$ due to the periodic boundary conditions. Therefore we have truncated the graph to $\Delta t = 16$. Also the correlator is missing a normalization factor and an overall term to be subtracted from it. The above (and the forthcoming) graphs are only to observe the behavior of the correlator. 

Notice that as the temperature is decreased, the correlator deviates from the expected exponential decay. To verify this behavior, we conducted a few test runs, which are presented here. 

### [SINGLEMIXED](/DATA_ANALYSIS/CPP/SINGLEMIXED/)

Rather than summing over all $X^M$, we consider only two $X^0$ and $X^1$. The corresponding code is at [bfssconfig.h Line 610](/MCSC-CPPCODE/src/bfssconfig.h#L610).

The graph obtained is 
![](/DATA_ANALYSIS/CPP/SINGLEMIXED/singlemixed.png)
