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

## Code
The project was to grasp the essence of lattice field theory, build on the C++ implementation of the BFSS model by Dr. Bergner and write and analyze observables in the simulation. The specific observables written by me are:

- Bosonic Energy: [bfssconfig.h  Line 517](/MCSC-CPPCODE/src/bfssconfig.h#L517)
- Fermionic Energy: [fermionmeasurements.h Line 36](/MCSC-CPPCODE/src/fermionmeasurements.h#L36)
- Gauge Invariant 4-point correlator (not normalised), $\int dt\left \langle ~\mathrm{Tr}(X^M(t)X^N(t)) ~ \mathrm{Tr}(X^M(t + \Delta t)X^N(t + \Delta t)) ~\right \rangle$: [bfssconfig.h Line 560](/MCSC-CPPCODE/src/bfssconfig.h#L560)

Further, in order to verify the simulation data, and scrutinize possible errors, I have also written the codes for 

- $\int dt \langle \mathrm{Tr}(X^0(t)X^0(t))~ \mathrm{Tr}(X^0(t + \Delta t)X^0(t + \Delta t))\rangle$: [bfssconfig.h Line 597](/MCSC-CPPCODE/src/bfssconfig.h#L597).

- $\int dt \langle \mathrm{Tr}(X^0(t)X^1(t))~ \mathrm{Tr}(X^0(t + \Delta t)X^1(t + \Delta t))\rangle$: [bfssconfig.h Line 610](/MCSC-CPPCODE/src/bfssconfig.h#L610).

Further, I also verified the results for the correlators by writing the observables in the FORTRAN code for the same by Dr. Masanori Hanada.

The complete code is in the folder [MCSC-CPPCODE](/MCSC-CPPCODE/).

## The Data

Unfortunately, I have lost access to the simulation data for Bosonic and Fermionic energies. Here I present the analysis for the Gauge Invariant 4-point correlator with only Bosonic action.

(The complete data and the jupyter notebooks of the analysis are in the folder [DATA_ANALYSIS](/DATA_ANALYSIS/))

### [Correlators at Different Temperatures](/DATA_ANALYSIS/CPP/DIFFTEMPS/)

Correlators vs $\Delta t$ for different temperatures

[Number of lattice sites = 32, Number of Colors = 9]

Temperatures 0.10 - 0.18   |  Temperatures 0.20 - 0.28
:-------------------------:|:-------------------------:
![](/DATA_ANALYSIS/CPP/DIFFTEMPS/temps_0.18-0.10.png)  |  ![](/DATA_ANALYSIS/CPP/DIFFTEMPS/temps_0.28-0.20.png)

Note that the correlator is symmetric about $\Delta t = 32/2$ due to the periodic boundary conditions. Therefore we have truncated the graph to $\Delta t = 16$. Also the correlator is missing a normalization factor and an overall term to be subtracted from it. The above (and the forthcoming) graphs are only to observe the behavior of the correlator. 

Notice that as the temperature is decreased, the correlator deviates from the expected exponential decay. 

To verify this behavior, I wrote the correlator observable in the existing FORTRAN code at https://github.com/gbergner/SYM1DMMMT and compared the results. The results are in the folder [FORTRAN/DIFFTEMPS](/DATA_ANALYSIS/FORTRAN/DIFFTEMPS) and are as follows

Temperatures 0.10 - 0.18   |  Temperatures 0.20 - 0.28
:-------------------------:|:-------------------------:
![](/DATA_ANALYSIS/FORTRAN/DIFFTEMPS/temps_0.18-0.10.png)  |  ![](/DATA_ANALYSIS/FORTRAN/DIFFTEMPS/temps_0.28-0.20.png)

We see that the correlator indeed displays the same behavior and the observed behavior is not an artifact of the C++ code. Therefore, to scrutinize the behavior, we conducted the following analysis in an attempt to pinpoint the source of the discrepancy.


### [3 Colors](/DATA_ANALYSIS/CPP/3COLORCORR/) and [6 Colors](/DATA_ANALYSIS/CPP/6COLORCORR/)

To verify that the same effect persists for a lower number of colors too, we considered the case of 6 and 3 colors. 

The data obtained is:

[Number of lattice sites = 16, Temperature = 0.10]

3 Colors   |  6 Colors
:-------------------------:|:-------------------------:
![](/DATA_ANALYSIS/CPP/3COLORCORR/3colorcorr.png)  |  ![](/DATA_ANALYSIS/CPP/6COLORCORR/6colorcorr.png)


### [Single Matrix Correlator](/DATA_ANALYSIS/CPP/SINGLE) and [Single Mixed Matrix Correlator](/DATA_ANALYSIS/CPP/SINGLEMIXED/)

Rather than summing over all $X^M$, we consider two cases, considering only $X^0$, and the other case only the mix term between $X^0$ and $X^1$. The corresponding observables are at [bfssconfig.h Line 597](/MCSC-CPPCODE/src/bfssconfig.h#L597) [bfssconfig.h Line 610](/MCSC-CPPCODE/src/bfssconfig.h#L610).

The graph obtained is:

[Number of lattice sites = 16, Temperature = 0.10, Number of colors = 9]

Single   |  Mixed
:-------------------------:|:-------------------------:
![](/DATA_ANALYSIS/CPP/SINGLE/single.png)  |  ![](/DATA_ANALYSIS/CPP/SINGLEMIXED/singlemixed.png)


### [No Commutator in the Action](/DATA_ANALYSIS/CPP/NOCOMM)

We further investigate the case of turning off the commutator term in the action. 

Turning off the commutator in the action gives the results as follows: 

[Temperature=0.10, Number of colors = 9, Number of lattice sites = 16]

![](/DATA_ANALYSIS/CPP/NOCOMM/nocomm.png)

We see that turning off the commutator term leads to the anomalous behavior turned off and the correlator displaying the expected exponential decay.

### [Correlators With The Fermionic Terms in The Action](/DATA_ANALYSIS/FORTRAN/3COLORWITHFERMIONS/)

To look at the behavior of the correlator under the complete action of the BFSS model, we run the same in the FORTRAN code (since at this stage, it is more efficient and fermionic simulations are exponentially more resource-consuming than bosonic simulations). The results are as follows:

[Temperature=0.10, Number of colors = 3, Number of lattice sites = 16]

![](/DATA_ANALYSIS/FORTRAN/3COLORWITHFERMIONS/fullAction.png)

## Conclusions

From the project, we observe that the 4-point gauge invariant correlator in the pure Bosonic sector of the BFSS model shows deviations from the expected exponential decay. The behavior is not an artifact of the C++ code, but rather the behavior of the model. 

We have further observed that turning off the commutator (interaction) term in the action leads to the expected exponential decay. Therefore we conclude that in the pure Bosonic sector, the presence of the $X-X$ interaction term that is a result of the dimensional reduction is responsible for the anomalous behavior.

Further, in the complete action with the Fermionic terms, the correlator displays the expected exponential decay. The reason for the suppression of the anomalous behavior in the presence of the Fermionic terms is not clear.