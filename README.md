# BFSS Model on the Lattice - Report

## The Model

The BFSS model is said to be dual to type II string theory. It is obtained from the reduction of the 9+1 dimensional Supersymmetric Yang-Mills to 0+1 dimensions. 
The complete action of the BFSS model is therefore given by 

$$ S_{M}={N}\int_{0}^{\beta}dt~\mathrm{Tr}\left\lbrace\frac{1}{2}\big(D_{t}X_{M}\big)^{2}-\frac{1}{4}[X_{M},X_{N}]^{2} + i{\bar{\psi}}\gamma^{10}D_{t}\psi\ -{\bar{\psi}}\gamma^{M}[X_{M},\psi]\right\rbrace $$

where $\gamma^M~(M=1, \cdots, 10)$ are the $16\times16$ upper right block matrices of the 10-d Gamma matrices $\Gamma^M$.

Therefore the model is the quantum mechanics of $N\times N$ matrices $X$ and $\psi$ (where $N$ is the number of colors).

## Project Objectives

The objective of this project was to understand the fundamental principles of lattice field theory, extend the C++ implementation of the BFSS model developed by Dr. Bergner, and formulate and analyze observables within the simulation framework. 

### Lattice Simulations

Lattice formulation is the only known non-perturbative, regularized formalism for Quantum Field Theories. In this, we discretize space-time and form a lattice on which the (finite) path integral is performed.

On a lattice, the gauge fields take up the role of connections between lattice points and, therefore, live on the links between the lattice sites. The action for the gauge fields would be constructed out of plaquette sums as 

$$ S = \frac{\beta}{N} \sum_{n}\sum_{\mu<\nu}\text{Re}~ \text{tr}~(\mathbb{I} - U_{\mu\nu}(n)) $$

As a presentation of the lattice field methods, refer to the file [rudimentary_lattsim.cpp](rudimentary_lattsim.cpp), which presents a rudimentary (not optimized) code for simulating $SU(2)$ gauge theory on a 1+1 dimensional lattice. This code is ugly, inefficient and prone to memory leaks. It is only to display the principles of the lattice simulations. Please note that I have not implemented any observables so far. The only part that has been implemented so far is the lattice action and the updating of the lattice according to the Markov Chain method.

For a more comprehensive introduction to lattice methods, I would point the reader to the book [Quantum Chromodynamics on the Lattice](#ref5)

<!-- On lattice regularization, the model action becomes  -->

<!-- $$ S_M = S_b + S_f $$ -->

<!-- with -->

<!-- $$ S_{b}\ =\displaystyle \frac{N}{2a}\sum_{t.M}\mathrm{Tr}\left(U X_{M}(t+a)U^{\dagger}-X_{M}(t)\right)^{2}-\frac{N a}{4}\sum_{t.M.N}\mathrm{Tr}[X_{M}(t),X_{N}(t)]^{2} $$ -->

<!-- ![equation](https://latex.codecogs.com/svg.image?&space;S_{b}\=\displaystyle\frac{N}{2a}\sum_{t.M}\mathrm{Tr}\left(U&space;X_{M}(t&plus;a)U^{\dagger}-X_{M}(t)\right)^{2}-\frac{N&space;a}{4}\sum_{t.M.N}\mathrm{Tr}[X_{M}(t),X_{N}(t)]^{2}) -->

<!-- $$\displaystyle S_{f}=i N\sum_{t}\mathrm{Tr}\bar{\psi}(t)\left(\begin{array}{c c}{{0}}&{{D_{+}}}\\ {{D_{-}}}&{{0}}\end{array}\right)\psi(t)-a N\sum_{t,M}\bar{\psi}(t)\gamma^{M}[X_{M}(t),\psi(t)] $$ -->

### BFSS Observables

The specific observables I implemented and investigated include:

- Bosonic Energy: [bfssconfig.h  Line 517](/MCSC-CPPCODE/src/bfssconfig.h#L517)
- Fermionic Energy: [fermionmeasurements.h Line 36](/MCSC-CPPCODE/src/fermionmeasurements.h#L36)
- Gauge Invariant 4-point correlator (not normalised), $\int dt\left \langle ~\mathrm{Tr}(X^M(t)X^M(t)) ~ \mathrm{Tr}(X^N(t + \Delta t)X^N(t + \Delta t)) ~\right \rangle$: [bfssconfig.h Line 560](/MCSC-CPPCODE/src/bfssconfig.h#L560)

To validate the simulation data and conduct a thorough examination of potential errors, I have developed supplementary code for

- $\int dt \langle \mathrm{Tr}(X^0(t)X^0(t))~ \mathrm{Tr}(X^0(t + \Delta t)X^0(t + \Delta t))\rangle$: [bfssconfig.h Line 597](/MCSC-CPPCODE/src/bfssconfig.h#L597).

- $\int dt \langle \mathrm{Tr}(X^0(t)X^1(t))~ \mathrm{Tr}(X^0(t + \Delta t)X^1(t + \Delta t))\rangle$: [bfssconfig.h Line 610](/MCSC-CPPCODE/src/bfssconfig.h#L610).

Additionally, to validate the results obtained for the correlators, I implemented the corresponding observables in the FORTRAN code developed by Dr. Masanori Hanada, which served as an independent verification of our findings.

The complete code is in the folder [MCSC-CPPCODE](/MCSC-CPPCODE/).

## Results

Due to unforeseen circumstances, the simulation data for Bosonic and Fermionic energies are no longer available. In this section, we present an analysis of the Gauge Invariant 4-point correlator, focusing solely on the Bosonic sector.

(The complete data and the jupyter notebooks of the analysis are in the folder [DATA_ANALYSIS](/DATA_ANALYSIS/))

### [Correlators at Different Temperatures](/DATA_ANALYSIS/CPP/DIFFTEMPS/)

Correlators vs $\Delta t$ for different temperatures

[Number of lattice sites = 32, Number of Colors = 9]

Temperatures 0.10 - 0.18   |  Temperatures 0.20 - 0.28
:-------------------------:|:-------------------------:
![](/DATA_ANALYSIS/CPP/DIFFTEMPS/temps_0.18-0.10.png)  |  ![](/DATA_ANALYSIS/CPP/DIFFTEMPS/temps_0.28-0.20.png)

It is important to note that the correlator exhibits symmetry about $\Delta t = (\text{Number of lattice sites})/2$ due to the periodic boundary conditions imposed. Consequently, we have truncated the graphical representation to $\Delta t = 16$. Furthermore, it should be noted that the correlator lacks a normalization factor and an overall term that should be subtracted from it. The graphs presented here, as well as those that follow, are intended solely to illustrate the behavior of the correlator.

Upon examination, it is evident that as the temperature decreases, the correlator's behavior deviates from the anticipated exponential decay.

To validate this observed phenomenon, we implemented the correlator observable in the existing FORTRAN codebase, which can be found at https://github.com/gbergner/SYM1DMMMT. Subsequently, we conducted a comparative analysis of the results. The outcomes of this investigation are available in the directory [FORTRAN/DIFFTEMPS](/DATA_ANALYSIS/FORTRAN/DIFFTEMPS) and are presented as follows

Temperatures 0.10 - 0.18   |  Temperatures 0.20 - 0.28
:-------------------------:|:-------------------------:
![](/DATA_ANALYSIS/FORTRAN/DIFFTEMPS/temps_0.18-0.10.png)  |  ![](/DATA_ANALYSIS/FORTRAN/DIFFTEMPS/temps_0.28-0.20.png)

The correlator exhibits consistent behavior across implementations, indicating that the observed phenomenon is not an artifact of the C++ codebase. Therefore, to scrutinize the behavior, we conducted the following analysis in an attempt to pinpoint the source of the discrepancy.


### [3 Colors](/DATA_ANALYSIS/CPP/3COLORCORR/) and [6 Colors](/DATA_ANALYSIS/CPP/6COLORCORR/)

To verify that the same effect persists for a lower number of colors too, we considered the case of 6 and 3 colors. 

The data obtained is:

[Number of lattice sites = 16, Temperature = 0.10]

3 Colors   |  6 Colors
:-------------------------:|:-------------------------:
![](/DATA_ANALYSIS/CPP/3COLORCORR/3colorcorr.png)  |  ![](/DATA_ANALYSIS/CPP/6COLORCORR/6colorcorr.png)


### [Single Matrix Correlator](/DATA_ANALYSIS/CPP/SINGLE) and [Single Mixed Matrix Correlator](/DATA_ANALYSIS/CPP/SINGLEMIXED/)

Rather than summing over all $X^M$, we consider two cases

- considering only $X^0$, for the correlator
- considering a case of mixed terms between $X^0$ and $X^1$ to ensure that the behavior persists even in the mixed case.
<!-- The corresponding observables are at [bfssconfig.h Line 597](/MCSC-CPPCODE/src/bfssconfig.h#L597) [bfssconfig.h Line 610](/MCSC-CPPCODE/src/bfssconfig.h#L610). -->

The graph obtained is:

[Number of lattice sites = 16, Temperature = 0.10, Number of colors = 9]

Single   |  Mixed
:-------------------------:|:-------------------------:
![](/DATA_ANALYSIS/CPP/SINGLE/single.png)  |  ![](/DATA_ANALYSIS/CPP/SINGLEMIXED/singlemixed.png)

The correlator for the mixed case too, although being negligible in magnitude, displays the same behavior as the correlator for the single case.


### [No Commutator in the Action](/DATA_ANALYSIS/CPP/NOCOMM)

We further investigate the case of turning off the commutator term in the action. 

Turning off the commutator in the action gives the results as follows: 

[Temperature=0.10, Number of colors = 9, Number of lattice sites = 16]

![](/DATA_ANALYSIS/CPP/NOCOMM/nocomm.png)

We observe that turning off the commutator term results in the suppression of the anomalous behavior, with the correlator exhibiting the expected exponential decay.

### [Correlators With The Fermionic Terms in The Action](/DATA_ANALYSIS/FORTRAN/3COLORWITHFERMIONS/)

To investigate the behavior of the correlator under the complete action of the BFSS model, we employ the FORTRAN implementation of the code. This choice is motivated by the superior computational efficiency of the FORTRAN code as compared to the C++ code which is currently under development, particularly for fermionic simulations, which are significantly more resource-intensive compared to their bosonic counterparts.

The results are as follows:

[Temperature=0.10, Number of colors = 3, Number of lattice sites = 16]

![](/DATA_ANALYSIS/FORTRAN/3COLORWITHFERMIONS/fullAction.png)

## Conclusions

In our investigation of the BFSS model, we have observed that the 4-point gauge invariant correlator in the purely Bosonic sector exhibits deviations from the expected exponential decay. This phenomenon is not attributable to any artifacts in the C++ implementation but rather reflects an intrinsic property of the model itself.

Further analysis reveals that when the commutator (interaction) term in the action is nullified, the correlator reverts to the expected exponential decay behavior. This leads us to conclude that within the purely Bosonic sector, the presence of the X-X interaction term is responsible for the observed anomalous behavior.

Interestingly, when considering the complete action inclusive of the Fermionic terms, the correlator demonstrates the expected exponential decay. However, the underlying mechanism by which the Fermionic terms suppress the anomalous behavior remains unclear and requires further investigation.

Some similar studies indicate that such an oscillatory behavior is a sign of the model exhibiting confinement [[6](#ref6)-[7](#ref7)]. This translates to the statement that the glueballs in the Bosonic sector of the BFSS model are confined. However, such a statement here would require further theoretical and numerical verification to be stated concretely and can be only stated now as a speculation and not a conclusive statement.

## References

[1] M. Hanada, “BFSS code manual.” 

[2] V. G. Filev and D. O’Connor, “The BFSS model on the lattice,” J. High Energ. Phys., vol. 2016, no. 5, p. 167, May 2016, doi: 10.1007/JHEP05(2016)167.

[3] M. Hanada, J. Nishimura, Y. Sekino, and T. Yoneya, “Direct test of the gauge-gravity correspondence for Matrix theory correlation functions,” J. High Energ. Phys., vol. 2011, no. 12, p. 20, Dec. 2011, doi: 10.1007/JHEP12(2011)020.

[4] M. Hanada, “What lattice theorists can do for superstring/M-theory,” Int. J. Mod. Phys. A, vol. 31, no. 22, p. 1643006, Aug. 2016, doi: 10.1142/S0217751X16430065.

<a id="ref5"></a>

[5] C. Gattringer and C. B. Lang, Quantum Chromodynamics on the Lattice: An Introductory Presentation, vol. 788. in Lecture Notes in Physics, vol. 788. Berlin, Heidelberg: Springer Berlin Heidelberg, 2010. doi: 10.1007/978-3-642-01850-3.

<a id="ref6"></a>

[6] O. Oliveira, D. Dudal, and P. J. Silva, “Glueball spectral densities from the lattice,” Oct. 29, 2012, arXiv: arXiv:1210.7794. Accessed: Jul. 24, 2024. [Online]. Available: http://arxiv.org/abs/1210.7794

<a id="ref7"></a>

[7] L. C. Loveridge, O. Oliveira, and P. J. Silva, “Schwinger function, confinement, and positivity violation in pure gauge QED,” Phys. Rev. D, vol. 106, no. 1, p. L011502, Jul. 2022, doi: 10.1103/PhysRevD.106.L011502.


## Appendix

Statistical analysis functions used in the analysis of the correlator:

```python
def average(list):
    return sum(list)/len(list)

def standarddeviation(list):
    averageofthelist = average(list)
    newlist = [ (elementofthelist - averageofthelist)**2 for elementofthelist in list ]
    sigmasquared  = average(newlist)
    sigma = sigmasquared**0.5
    return sigma

def autocorrelation(list, lag):
    xi = list[:-lag]
    xiplust = list[lag:]
    xitimesxiplust = [ i*iplust for (i, iplust) in zip(xi, xiplust) ]
    averageofxitimexiplust = average(xitimesxiplust)
    averageofxi = average(xi)
    averageofxiplust = average(xiplust)
    autocorrelationvalue = averageofxitimexiplust - (averageofxi*averageofxiplust)
    return autocorrelationvalue

def integratedautocorrelation(list):
    selfcorrelation = standarddeviation(list)
    integratedautocorrelationvalues = [1/2]
    for lag in range(1, len(list)):
        normalisedautocorrelation = autocorrelation(list, lag)/selfcorrelation
        integratedautocorrelationvalues.append(integratedautocorrelationvalues[-1] + normalisedautocorrelation)

    return integratedautocorrelationvalues

def error(list):
    baseerror = standarddeviation(list)/((len(list))**0.5)
    integratedautocorrelationvalues = integratedautocorrelation(list)
    for i in range(1,len(integratedautocorrelationvalues)):
        if integratedautocorrelationvalues[i-1] > integratedautocorrelationvalues[i]:
            newerror = ((2*integratedautocorrelationvalues[i])**0.5)*baseerror
            break
    else:
        newerror = ((2*integratedautocorrelationvalues[-1])**0.5)*baseerror
    
    return newerror
```
