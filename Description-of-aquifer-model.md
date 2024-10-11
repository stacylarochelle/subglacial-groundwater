### Equations
*Equations in this section are taken from Chapter 7 of Wang (2000) Theory of Linear Poroelasticity*

This poroelastic model solves for pore pressure $p$ and displacement components $u$ and $v$ inside an aquifer cross-section in a state of plane strain, i.e., meaning that the displacement components are independent of $z$, such that strains $\varepsilon_{zz} = \varepsilon_{xz} = \varepsilon_{yz}= 0$. All quantities considered here describe changes with respect to a reference state (e.g., hydro-static equilibrium). The aquifer has spatially uniform hydro-mechanical properties. 

To solve for our 3 unknown $p$, $u$ and $v$, we have to solve 2 mechanical equations and 1 fluid pressure diffusion equation at every gridpoint. The two mechanical equations to be solved are obtained by substituting the constitutive equation relating stress ($\sigma_{ij}$), strains ($\varepsilon_{ij}$) and pore pressure ($p$): 
$$\sigma_{ij} = 2G \varepsilon_{ij} + 2G\frac{\nu}{1-2\nu}\varepsilon_{kk}\delta_{ij}-\alpha p \delta_{ij} \tag{1} $$
and the definition of infinitesimal strain: 
$$\varepsilon_{ij} = \frac{1}{2}\left(\frac{\partial u_i}{\partial x_j}+\frac{\partial u_j}{\partial x_i}\right)\tag{2} $$
into the mechanical equilibrium equations in plane strain with zero body forces: 
$$\frac{\partial \sigma_{xx}}{\partial x}+\frac{\partial \sigma_{xy}}{\partial y}=0\tag{3} $$
$$\frac{\partial \sigma_{xy}}{\partial x}+\frac{\partial \sigma_{yy}}{\partial y}=0\tag{4} $$
such that: 
$$G\left(\frac{\partial^2 u}{\partial x^2}+\frac{\partial^2 u}{\partial y^2}\right)+ \frac{G}{1-2\nu}\left(\frac{\partial^2 u}{\partial x^2}+\frac{\partial^2 v}{\partial x \partial y}\right)=\alpha \frac{\partial p}{\partial x}\tag{5} $$
$$G\left(\frac{\partial^2 v}{\partial x^2}+\frac{\partial^2 v}{\partial y^2}\right)+ \frac{G}{1-2\nu}\left(\frac{\partial^2 u}{\partial x \partial y}+\frac{\partial^2 v}{\partial y^2}\right)=\alpha \frac{\partial p}{\partial y}\tag{6} $$

Here $G$ is the elastic shear modulus, $\nu$ is the Poisson ratio, $\alpha$ is the poroelastic Biot-Willis coefficient and $\varepsilon_{kk} = \varepsilon_{xx}+\varepsilon_{yy}$ is the volumetric strain.

The plane-strain fluid diffusion equation in terms of displacements is as follows:
$$\alpha \frac{\partial}{\partial t}\left(\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}\right)+S_{\varepsilon}\frac{\partial p}{\partial t} = \frac{k}{\mu}\left(\frac{\partial^2 p}{\partial x^2}+\frac{\partial^2 p}{\partial y^2}\right)\tag{7} $$
where $\alpha$ is, again, the poroelastic Biot-Willis coefficient, $k$ is permeability and $\mu$ is water viscosity.  $S_\varepsilon$ is the specific storage coefficient at constant strain and is equal to: 
$$S_\varepsilon=\frac{\alpha^2(1-2\nu_u)(1-2\nu)}{2G(\nu_u-\nu)}\tag{8} $$
where $\nu_u$ is the undrained Poisson ratio. 

### Fully-implicit finite-difference implementation
*Inspired from https://www-udc.ig.utexas.edu/external/becker/Geodynamics557.pdf*

We first discretize Eqs. (5), (6) and (7) using second-order, centered finite differences in space and forward first-order differences in time: 
$$ \frac{2-2\nu}{1-2\nu}\left(\frac{u^{n+1}_{i+1,j}-2u^{n+1}_{i,j}+u^{n+1}_{i-1,j}}{(\Delta x)^2}\right)+\left(\frac{u^{n+1}_{i,j+1}-2u^{n+1}_{i,j}+u^{n+1}_{i,j-1}}{(\Delta y)^2}\right)+\\
\frac{1}{1-2\nu}\left(\frac{v^{n+1}_{i+1,j+1}-v^{n+1}_{i+1,j-1}-v^{n+1}_{i-1,j+1}+v^{n+1}_{i-1,j-1}}{4\Delta x\Delta y}\right) =\frac{\alpha}{G}\frac{p^{n+1}_{i+1,j}-p^{n+1}_{i-1,j}}{2\Delta x}\tag{9}$$

$$ \frac{2-2\nu}{1-2\nu}\left(\frac{v^{n+1}_{i,j+1}-2v^{n+1}_{i,j}+v^{n+1}_{i,j-1}}{(\Delta y)^2}\right)+\left(\frac{v^{n+1}_{i+1,j}-2v^{n+1}_{i,j}+v^{n+1}_{i-1,j}}{(\Delta x)^2}\right)+\\
\frac{1}{1-2\nu}\left(\frac{u^{n+1}_{i+1,j+1}-u^{n+1}_{i+1,j-1}-u^{n+1}_{i-1,j+1}+u^{n+1}_{i-1,j-1}}{4\Delta x\Delta y}\right) =\frac{\alpha}{G}\frac{p^{n+1}_{i,j+1}-p^{n+1}_{i,j-1}}{2\Delta y}\tag{10}$$

$$\frac{\alpha}{\Delta t}\left(\frac{u^{n+1}_{i+1,j}-u^{n+1}_{i-1,j}+u^{n}_{i+1,j}-u^{n}_{i-1,j}}{2\Delta x}+\frac{v^{n+1}_{i,j+1}-v^{n+1}_{i,j-1}+v^{n}_{i,j+1}-v^{n}_{i,j-1}}{2\Delta y}\right)+S_\varepsilon\left(\frac{p^{n+1}_{i,j}-p^n_{i,j}}{\Delta t}\right)=\\
\frac{k}{\mu}\left(\frac{p^{n+1}_{i+1,j}-2p^{n+1}_{i,j}+p^{n+1}_{i-1,j}}{2\Delta x}+\frac{p^{n+1}_{i,j+1}-2p^{n+1}_{i,j}+p^{n+1}_{i,j-1}}{2\Delta y}\right)\tag{11}$$

Where the superscript indicates the time step and the subscripts the spatial positions. We can rearrange the equations to have all unknown quantities (those at time $n+1$) are on the LHS and known quantities (those at time $n$) on the RHS: 

$$ \frac{2-2\nu}{1-2\nu}\left(\frac{u^{n+1}_{i+1,j}-2u^{n+1}_{i,j}+u^{n+1}_{i-1,j}}{(\Delta x)^2}\right)+\left(\frac{u^{n+1}_{i,j+1}-2u^{n+1}_{i,j}+u^{n+1}_{i,j-1}}{(\Delta y)^2}\right)+\\
\frac{1}{1-2\nu}\left(\frac{v^{n+1}_{i+1,j+1}-v^{n+1}_{i+1,j-1}-v^{n+1}_{i-1,j+1}+v^{n+1}_{i-1,j-1}}{4\Delta x\Delta y}\right) -\frac{\alpha}{G}\frac{p^{n+1}_{i+1,j}-p^{n+1}_{i-1,j}}{2\Delta x}=0 \tag{12}$$

$$ \frac{2-2\nu}{1-2\nu}\left(\frac{v^{n+1}_{i,j+1}-2v^{n+1}_{i,j}+v^{n+1}_{i,j-1}}{(\Delta y)^2}\right)+\left(\frac{v^{n+1}_{i+1,j}-2v^{n+1}_{i,j}+v^{n+1}_{i-1,j}}{(\Delta x)^2}\right)+\\
\frac{1}{1-2\nu}\left(\frac{u^{n+1}_{i+1,j+1}-u^{n+1}_{i+1,j-1}-u^{n+1}_{i-1,j+1}+u^{n+1}_{i-1,j-1}}{4\Delta x\Delta y}\right) -\frac{\alpha}{G}\frac{p^{n+1}_{i,j+1}-p^{n+1}_{i,j-1}}{2\Delta y}=0 \tag{13}$$

$$\alpha\left(\frac{u^{n+1}_{i+1,j}-u^{n+1}_{i-1,j}}{2\Delta x}+\frac{v^{n+1}_{i,j+1}-v^{n+1}_{i,j-1}}{2\Delta y}\right)-\frac{k\Delta t}{\mu}\left(\frac{p^{n+1}_{i+1,j}+p^{n+1}_{i-1,j}}{2\Delta x}+\frac{p^{n+1}_{i,j+1}+p^{n+1}_{i,j-1}}{2\Delta y}\right)\\+\left(\frac{2k\Delta t}{\mu}\left(\frac{1}{2\Delta x}+\frac{1}{2\Delta y}\right)+S_\varepsilon\right)p^{n+1}_{i,j}= -\alpha\left(\frac{u^{n}_{i+1,j}-u^{n}_{i-1,j}}{2\Delta x}+\frac{v^{n}_{i,j+1}-v^{n}_{i,j-1}}{2\Delta y}\right)+S_\varepsilon p^n_{i,j} \tag{14}$$

such that the equations form a system of equation of type $\mathbf{Ax} = \mathbf{b}$ where $m\times m$ matrix $\mathbf{A}$ stores all the constant coefficients, $m \times 1$ vector $\mathbf{x}$ contains the unknown quantities $u$, $v$ and $p$ at time $n+1$ and $m\times 1$ vector $\mathbf{b}$ contains the RHS of Eqs. (12-14), which are either $\mathbf{0}$ or functions of known quantities at time $n$.

- Note that to make the solution second-order accurate in time, we could, fairly easily, implement a Crank-Nicolson scheme.

### Boundary conditions
In addition to the three equations above, we also have to add 3 additional equations at each boundary to $\mathbf{A}$ to implement the boundary conditions. 

*Mechanical boundary conditions:*
- Fixed displacement boundary conditions: e.g., $u^{n+1}_{1,j}=0$ to prescribe zero horizontal displacement on the left boundary. 
- Fixed displacement gradient boundary condition: e.g., $v^{n+1}_{1,j}-v^{n+1}_{3,j}=0$ to prescribe zero gradient in vertical displacements on the left boundary. 
- Fixed stress boundary conditions are implemented through the constitutive relationships:
	- e.g., to set $\sigma_{yy}(x,0) = \sigma_{ice}$: 
$$ \sigma_{yy} = 2G\varepsilon_{yy}+2G\frac{\nu}{1-2\nu}(\epsilon_{xx}+\varepsilon_{yy})-\alpha p$$
$$ \sigma_{ice} = 2G\frac{\partial v}{\partial y}+2G\frac{\nu}{1-2\nu}\left(\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}\right)-\alpha p$$
$$2G\left(\frac{v^{n+1}_{i,3}-v^{n+1}_{i,1}}{2\Delta y}\right)+2G\frac{\nu}{1-2\nu}\left(\frac{u^{n+1}_{i+1,2}-v^{n+1}_{i-1,2}}{2\Delta x}+\frac{v^{n+1}_{i,3}-v^{n+1}_{i,1}}{2\Delta y}\right)-\alpha p^{n+1}_{i,2} = \sigma_{ice}$$
- e.g., to set $\sigma_{xy}(x,0) = \tau_{ice}$: 
$$ \sigma_{xy} = 2G\varepsilon_{xy} $$
$$\tau_{ice}= G\left(\frac{\partial u}{\partial y}+\frac{\partial v}{\partial x}\right)$$ 
$$ G\left(\frac{u^{n+1}_{i,3}-u^{n+1}_{i,1}}{2\Delta y}+\frac{v^{n+1}_{i+1,2}-v^{n+1}_{i-1,2}}{2\Delta x}\right)=\tau_{ice}$$

Again, all unknown quantities at $n+1$ to be solved for are on the LHS of the equations and the known quantities -- $\sigma_{ice}$ and $\tau_{ice}$ -- are on the RHS. 

*Fluid boundary conditions:*
- Fixed pore pressure boundary conditions: e.g., $p^{n+1}_{i,1} =0$ or $p^{n+1}_{i,1} =p_{subglacial}$ to prescribe pore pressure at the top boundary. 
- No-flow boundary conditions: e.g., $p^{n+1}_{i,m+1}-p^{n+1}_{i,m-1}=0$ to set a no-flow boundary condition at the aquifer-bedrock interface.

### Solving for $\textbf{x}$:
- The distributions of $u$, $v$ and $p$ at time $n+1$ stored in $\textbf{x}$ can then be solved for by taking the inverse of $\textbf{A}$ and multiplying by $\textbf{b}$:
$$\textbf{x} = \textbf{A}^{-1}\textbf{b}$$
- If the boundary conditions remain of the same type at all times, then $\textbf{A}$ is constant through time and, therefore, only needs to be inverted once.
- If the type of boundary condition evolves over time, however, $\textbf{A}$ has to be inverted for at each timestep, making the whole computation much more time intensive. 
	-  In the problem set up we are first considering where ice acts as a confining unit to the aquifer below (e.g., $\partial p/\partial y = 0$) and seawater sets a constant pressure at the top of the aquifer (e.g., $p = 0$ for constant sea level), the fluid boundary condition at the top of the aquifer will change from a no-flow boundary condition to a fixed pore pressure one. This means that $\textbf{A}$ and, thus, its inverse will have to be updated at each time step. 
	-  For the mechanical boundary conditions, only the RHS of the equations will change with grounding line migration (e.g., from $\sigma_{yy} = \sigma_{ice}$ to $\sigma_{yy} = 0$), so the rows of $\textbf{A}$ storing these equations do not need to be updated.
	- The good news is that once we add subglacial hydrology to the model and prescribe $p = p_{subglacial}$ at the top of the aquifer, $\textbf{A}$ will be constant through time. 
- Since we are only changing a few rows of $\textbf{A}$ at a time, it may be possible to significantly speed up the computation of its inverse by making use of the inverse from the previous time step through the Woodbury matrix identity. 
- To speed up inverse computation, equations should be arranged such that the matrix $\mathbf{A}$ is as "diagonal" as possible. 
