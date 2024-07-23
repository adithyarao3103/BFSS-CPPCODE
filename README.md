# BFSS Model on the Lattice

## The Model

The BFSS model is said to be dual to type II string theory. It is obtained from the reduction of the 9+1 dimensional Supersymmetric Yang-Mills to 0+1 dimensions. 
The complete action of the BFSS model is therefore given by 

$$ S_{M}={N}\int_{0}^{\beta}dt~\mathrm{Tr}\left\lbrace\frac{1}{2}\big(D_{t}X_{M}\big)^{2}-\frac{1}{4}[X_{M},X_{N}]^{2} + i{\bar{\psi}}\gamma^{10}D_{t}\psi\ -{\bar{\psi}}\gamma^{M}[X_{M},\psi]\right\rbrace $$

where $\gamma^M~(M=1, \cdots, 10)$ are the $16\times16$ upper right block matrices of the 10-d Gamma matrices $\Gamma^M$.

Therefore the model is the quantum mechanics of $N\times N$ matrices $X$ and $\psi$ (where $N$ is the number of colors).

On lattice regularization, the model action becomes 

$$ S_M = S_b + S_f $$

with

$$ S_{b}\ =\ \frac{N}{2a}\sum_{t.M}\mathrm{Tr}\left(U X_{M}(t+a)U^{\dagger}-X_{M}(t)\right)^{2}-\frac{N a}{4}\sum_{t.M.N}\mathrm{Tr}[X_{M}(t),X_{N}(t)]^{2} $$

$$S_{f}=i N\sum_{t}\mathrm{Tr}\bar{\psi}(t)\left(\begin{array}{c c}{{0}}&{{D_{+}}}\\ {{D_{-}}}&{{0}}\end{array}\right)\psi(t)-a N\sum_{t,M}\bar{\psi}(t)\gamma^{M}[X_{M}(t),\psi(t)], $$

## The Project
