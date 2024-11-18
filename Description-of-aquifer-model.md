### Equations
*Equations in this section are taken from Chapter 7 of Wang (2000) Theory of Linear Poroelasticity*

This poroelastic model solves for pore pressure $p$ and displacement components $u$ and $v$ inside an aquifer cross-section in a state of plane strain, meaning that the displacement components are independent of $z$, such that strains $\varepsilon_{zz} = \varepsilon_{xz} = \varepsilon_{yz}= 0$. All quantities considered here describe changes with respect to a reference state (e.g., hydrostatic equilibrium). The aquifer has spatially uniform hydro-mechanical properties. 

To solve for our 3 unknown $p$ (pore pressure), $u$ (horizontal displacement) and $v$ (vertical displacement), we have to solve 2 mechanical equations and 1 fluid pressure diffusion equation at every gridpoint. The two mechanical equations to be solved are obtained by substituting the constitutive equation relating stress ($\sigma_{ij}$), strains ($\varepsilon_{ij}$) and pore pressure ($p$): 

$$\sigma_{ij} = 2G \varepsilon_{ij} + 2G\frac{\nu}{1-2\nu}\varepsilon_{kk}\delta_{ij}-\alpha p \delta_{ij}$$

and the definition of infinitesimal strain: 

$$\varepsilon_{ij} = \frac{1}{2}\left(\frac{\partial u_i}{\partial x_j}+\frac{\partial u_j}{\partial x_i}\right)$$

into the mechanical equilibrium equations in plane strain with zero body forces: 

$$\frac{\partial \sigma_{xx}}{\partial x}+\frac{\partial \sigma_{xy}}{\partial y}=0$$

$$\frac{\partial \sigma_{xy}}{\partial x}+\frac{\partial \sigma_{yy}}{\partial y}=0$$

such that: 

$$G\left(\frac{\partial^2 u}{\partial x^2}+\frac{\partial^2 u}{\partial y^2}\right)+ \frac{G}{1-2\nu}\left(\frac{\partial^2 u}{\partial x^2}+\frac{\partial^2 v}{\partial x \partial y}\right)=\alpha \frac{\partial p}{\partial x}$$

$$G\left(\frac{\partial^2 v}{\partial x^2}+\frac{\partial^2 v}{\partial y^2}\right)+ \frac{G}{1-2\nu}\left(\frac{\partial^2 u}{\partial x \partial y}+\frac{\partial^2 v}{\partial y^2}\right)=\alpha \frac{\partial p}{\partial y}$$

Here $G$ is the elastic shear modulus, $\nu$ is the Poisson ratio, $\alpha$ is the poroelastic Biot-Willis coefficient and $\varepsilon_{kk} = \varepsilon_{xx}+\varepsilon_{yy}$ is the volumetric strain.

The plane-strain fluid diffusion equation in terms of displacements is as follows:

$$\alpha \frac{\partial}{\partial t}\left(\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}\right)+S_{\varepsilon}\frac{\partial p}{\partial t} = \frac{k}{\mu}\left(\frac{\partial^2 p}{\partial x^2}+\frac{\partial^2 p}{\partial y^2}\right)$$

where $\alpha$ is, again, the poroelastic Biot-Willis coefficient, $k$ is permeability and $\mu$ is water viscosity.  $S_\varepsilon$ is the specific storage coefficient at constant strain and is equal to: 

$$S_\varepsilon=\frac{\alpha^2(1-2\nu_u)(1-2\nu)}{2G(\nu_u-\nu)}$$

where $\nu_u$ is the undrained Poisson ratio. 

### Fully-implicit finite-difference implementation

We first discretize our three equations using second-order, centered finite differences in space and forward first-order differences in time: 

$$\frac{2-2\nu}{1-2\nu}\left(\frac{u_{i+1,j}^{n+1}-2u_{i,j}^{n+1}+u_{i-1,j}^{n+1}}{(\Delta x)^2}\right)+\left(\frac{u_{i,j+1}^{n+1}-2u_{i,j}^{n+1}+u_{i,j-1}^{n+1}}{(\Delta y)^2}\right)+
\frac{1}{1-2\nu}\left(\frac{v_{i+1,j+1}^{n+1}-v_{i+1,j-1}^{n+1}-v_{i-1,j+1}^{n+1}+v_{i-1,j-1}^{n+1}}{4\Delta x\Delta y}\right) =\frac{\alpha}{G}\frac{p_{i+1,j}^{n+1}-p_{i-1,j}^{n+1}}{2\Delta x}$$

$$\frac{2-2\nu}{1-2\nu}\left(\frac{v_{i,j+1}^{n+1}-2v_{i,j}^{n+1}+v_{i,j-1}^{n+1}}{(\Delta y)^2}\right)+\left(\frac{v_{i+1,j}^{n+1}-2v_{i,j}^{n+1}+v_{i-1,j}^{n+1}}{(\Delta x)^2}\right)+
\frac{1}{1-2\nu}\left(\frac{u_{i+1,j+1}^{n+1}-u_{i+1,j-1}^{n+1}-u_{i-1,j+1}^{n+1}+u_{i-1,j-1}^{n+1}}{4\Delta x\Delta y}\right) =\frac{\alpha}{G}\frac{p_{i,j+1}^{n+1}-p_{i,j-1}^{n+1}}{2\Delta y}$$

$$\frac{\alpha}{\Delta t}\left(\frac{u_{i+1,j}^{n+1}-u_{i-1,j}^{n+1}+u_{i+1,j}^{n}-u_{i-1,j}^{n}}{2\Delta x}+\frac{v_{i,j+1}^{n+1}-v_{i,j-1}^{n+1}+v_{i,j+1}^{n}-v_{i,j-1}^{n}}{2\Delta y}\right)+S_\varepsilon\left(\frac{p_{i,j}^{n+1}-p_{i,j}^n}{\Delta t}\right)=\frac{k}{\mu}\left(\frac{p_{i+1,j}^{n+1}-2p_{i,j}^{n+1}+p_{i-1,j}^{n+1}}{\Delta x^2}+\frac{p_{i,j+1}^{n+1}-2p_{i,j}^{n+1}+p_{i,j-1}^{n+1}}{\Delta y^2}\right)$$

Where the superscript indicates the time step and the subscripts the spatial positions. We can rearrange the equations such that all unknown quantities (those at time $n+1$) are on the LHS and known quantities (those at time $n$) on the RHS: 

$$\frac{2-2\nu}{1-2\nu}\left(\frac{u_{i+1,j}^{n+1}+u_{i-1,j}^{n+1}}{(\Delta x)^2}\right)+\left(\frac{u_{i,j+1}^{n+1}+u_{i,j-1}^{n+1}}{(\Delta y)^2}\right)-2\left(\frac{2-2\nu}{1-2\nu}\frac{1}{(\Delta x)^2}+\frac{1}{(\Delta y)^2}\right)u_{i,j}^{n+1}\\
+\frac{1}{1-2\nu}\left(\frac{v_{i+1,j+1}^{n+1}-v_{i+1,j-1}^{n+1}-v_{i-1,j+1}^{n+1}+v_{i-1,j-1}^{n+1}}{4\Delta x\Delta y}\right) -\frac{\alpha}{G}\frac{p_{i+1,j}^{n+1}-p_{i-1,j}^{n+1}}{2\Delta x}=0$$

$$\frac{2-2\nu}{1-2\nu}\left(\frac{v_{i,j+1}^{n+1}+v_{i,j-1}^{n+1}}{(\Delta y)^2}\right)+\left(\frac{v_{i+1,j}^{n+1}+v_{i-1,j}^{n+1}}{(\Delta x)^2}\right)-2\left(\frac{2-2\nu}{1-2\nu}\frac{1}{(\Delta y)^2}+\frac{1}{(\Delta x)^2}\right)v_{i,j}^{n+1}\\
+\frac{1}{1-2\nu}\left(\frac{u_{i+1,j+1}^{n+1}-u_{i+1,j-1}^{n+1}-u_{i-1,j+1}^{n+1}+u_{i-1,j-1}^{n+1}}{4\Delta x\Delta y}\right) -\frac{\alpha}{G}\frac{p_{i,j+1}^{n+1}-p_{i,j-1}^{n+1}}{2\Delta y}=0$$

$$\frac{1}{\Delta t}\left(\frac{u_{i+1,j}^{n+1}-u_{i-1,j}^{n+1}}{2\Delta x}+\frac{v_{i,j+1}^{n+1}-v_{i,j-1}^{n+1}}{2\Delta y}\right)-\frac{k}{\alpha \mu}\left(\frac{p_{i+1,j}^{n+1}+p_{i-1,j}^{n+1}}{\Delta x^2}+\frac{p_{i,j+1}^{n+1}+p_{i,j-1}^{n+1}}{\Delta y^2}\right)+\left(\frac{2k}{\alpha\mu}\left(\frac{1}{\Delta x^2}+\frac{1}{\Delta y^2}\right)+\frac{S_\varepsilon}{\alpha \Delta t}\right)p_{i,j}^{n+1}= -\frac{1}{\Delta t}\left(\frac{u_{i+1,j}^{n}-u_{i-1,j}^{n}}{2\Delta x}+\frac{v_{i,j+1}^{n}-v_{i,j-1}^{n}}{2\Delta y}\right)+\frac{S_\varepsilon}{\alpha \Delta t} p_{i,j}^n$$

such that the equations form a system of equation of type $\mathbf{Ax} = \mathbf{b}$ where $m\times m$ matrix $\mathbf{A}$ stores all the constant coefficients, $m \times 1$ vector $\mathbf{x}$ contains the unknown quantities $u$, $v$ and $p$ at time $n+1$ and $m\times 1$ vector $\mathbf{b}$ contains the RHS of the equations, which are either $\mathbf{0}$ or functions of known quantities at time $n$.

- Note that to make the solution second-order accurate in time, we would have to implement a Crank-Nicolson scheme.

### Boundary conditions
Model 1.1 has the following boundary conditions: 

**Left boundary ($i==0$)**: $u = 0$, $\frac{\partial v}{\partial x} = 0$, and $\frac{\partial p}{\partial x} = 0$ such that the problem is horizontally symmetric

These are implemented by:
- Replacing the first mechanical equation by $u_{i,j}^{n+1} = 0$ and setting $u_{i-1,j}^{n+1}=0$ in $\mathbf{A}$ and $u_{i-1,j}^{n}=0$ in the $\mathbf{b}$ vector;
- Replacing $v_{i-1,j}^{n+1}$ with $v_{i+1,j}^{n+1}$ in $\mathbf{A}$ since $v_{i+1,j}^{n+1} - v_{i-1,j}^{n+1} = 0$;
- Replacing $p_{i-1,j}^{n+1}$ with $p_{i+1,j}^{n+1}$ in $\mathbf{A}$ since $p_{i+1,j}^{n+1} - p_{i-1,j}^{n+1} = 0$.

**Right boundary ($i==nx-1$)**: $u = 0$, $v = 0$, and $\frac{\partial p}{\partial x} = 0$

These are implemented by:
- Replacing the first mechanical equation by $u_{i,j}^{n+1} = 0$ and setting $u_{i+1,j}^{n+1}=0$ in $\mathbf{A}$ and $u_{i+1,j}^{n}=0$ in the $\mathbf{b}$ vector;
- Replacing the first mechanical equation by $v_{i,j}^{n+1} = 0$ and setting $v_{i+1,j}^{n+1}=0$ in $\mathbf{A}$;
- Replacing $p_{i+1,j}^{n+1}$ with $p_{i-1,j}^{n+1}$ in $\mathbf{A}$ since $p_{i+1,j}^{n+1} - p_{i-1,j}^{n+1} = 0$.

**Bottom boundary ($j==0$)**: $u = 0$, $v = 0$, and $\frac{\partial p}{\partial y} = 0$

These are implemented by:
- Replacing the first mechanical equation by $u_{i,j}^{n+1} = 0$ and setting $u_{i,j-1}^{n+1}=0$ in $\mathbf{A}$;
- Replacing the first mechanical equation by $v_{i,j}^{n+1} = 0$ and setting $v_{i,j-1}^{n+1}=0$ in $\mathbf{A}$ and $v_{i,j-1}^{n}=0$ in $\mathbf{b}$;
- Replacing $p_{i,j-1}^{n+1}$ with $p_{i,j+1}^{n+1}$ in $\mathbf{A}$ since $p_{i,j+1}^{n+1} - p_{i,j-1}^{n+1} = 0$.

**Top boundary**: $\sigma_{yy} = \sigma_{ice}$, $\sigma_{xy} = \tau_{ice}$

Fixed stress boundary conditions are implemented through the constitutive relationships.

- For $\sigma_{yy}(x,0) = \sigma_{ice}$:
 
$$\sigma_{yy} = 2G\varepsilon_{yy}+2G\frac{\nu}{1-2\nu}(\epsilon_{xx}+\varepsilon_{yy})-\alpha p$$

$$\sigma_{ice} = 2G\frac{\partial v}{\partial y}+2G\frac{\nu}{1-2\nu}\left(\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}\right)-\alpha p$$

$$\frac{v_{i,j+1}^{n+1}-v_{i,j-1}^{n+1}}{2\Delta y} = \left(\frac{1-2\nu}{1-\nu}\right)\frac{1}{2G}(\sigma_{ice,i,j}^{n+1}+\alpha p_{i,j}^{n+1})-\left(\frac{\nu}{1-\nu}\right)\left(\frac{u_{i+1,j}^{n+1}-u_{i-1,j}^{n+1}}{2 \Delta x}\right)$$

- For $\sigma_{xy}(x,0) = \tau_{ice}$:

$$ \sigma_{xy} = 2G\varepsilon_{xy} $$

$$\tau_{ice}= G\left(\frac{\partial u}{\partial y}+\frac{\partial v}{\partial x}\right)$$ 

$$\frac{u_{i,j+1}^{n+1}-u_{i,j-1}^{n+1}}{2\Delta y}= -\frac{v_{i+1,j}^{n+1}-v_{i-1,j}^{n+1}}{2\Delta x}+\frac{\tau_{ice,i,j}^{n+1}}{G}$$

These expressions can then be used to substitute the boundary values in the fluid diffusion equation and the 2 mechanical equations, such that the three equations become: 

$$\frac{2-2\nu}{1-2\nu}\left(\frac{u_{i+1,j}^{n+1}+u_{i-1,j}^{n+1}}{\Delta x ^2}\right)+\frac{2u_{i,j-1}^{n+1}}{\Delta y^2}-\left(\frac{v_{i+1,j}^{n+1}-v_{i+1,j}^{n-1}}{\Delta x \Delta y}\right)-2\left(\frac{2-2\nu}{1-2\nu(\Delta x)^2}+\frac{1}{\Delta y^2}-\frac{\nu}{(1-2\nu)(1-\nu)4(\Delta x)^2}\right)u_{i,j}^{n+1} -\frac{\nu}{(1-2\nu)(1-\nu)}\left(\frac{u_{i+2,j}^{n+1}+u_{i-2,j}^{n+1}}{4(\Delta x)^2}\right)+\left(\frac{2\nu-1}{2-2\nu}\right)\frac{\alpha}{G}\left(\frac{p_{i+1,j}^{n+1}-p_{i-1,j}^{n+1}}{2\Delta x}\right) = -\frac{2\tau_{ice,i,j}^{n+1}}{G\Delta y}-\frac{1}{2G(1-\nu)}\left(\frac{\sigma_{ice,i+1,j}^{n+1}-\sigma_{ice,i-1,j}^{n+1}}{2\Delta x}\right)$$

$$\frac{2-2\nu}{1-2\nu}\left(\frac{2v_{i,j-1}^{n+1}}{(\Delta y)^2}\right)+\left(\frac{v_{i+1,j}^{n+1}+v_{i-1,j}^{n+1}}{(\Delta x)^2}\right)-\frac{1}{1-2\nu}\left(\frac{v_{i+2,j}^{n+1}+v_{i-2,j}^{n+1}}{4(\Delta x)^2}\right)-2\left(\frac{2-2\nu}{1-2\nu}\frac{1}{(\Delta y)^2}+\frac{1}{(\Delta x)^2}-\frac{1}{4(1-2\nu)(\Delta x)^2}\right)v_{i,j}^{n+1}\\
-\frac{2\nu}{1-2\nu}\left(\frac{u_{i+1,j}^{n+1}-u_{i-1,j}^{n+1}}{\Delta x\Delta y}\right) -\frac{\alpha}{G}\frac{p_{i,j+1}^{n+1}-p_{i,j-1}^{n+1}}{2\Delta y} +\frac{2\alpha}{G\Delta y}p_{i,j}^{n+1}=-\frac{2\sigma_{ice,i,j}^{n+1}}{G\Delta y}-\frac{1}{G(1-2\nu)}\left(\frac{\tau_{ice,i+1,j}^{n+1}-\tau_{ice,i-1,j}^{n+1}}{2\Delta x}\right)$$

$$\frac{1}{\Delta t}\left(\frac{1-2\nu}{1-\nu}\right)\left(\frac{u_{i+1,j}^{n+1}-u_{i-1,j}^{n+1}}{2\Delta x}\right)-\frac{k}{\alpha \mu}\left(\frac{p_{i+1,j}^{n+1}+p_{i-1,j}^{n+1}}{\Delta x^2}+\frac{p_{i,j+1}^{n+1}+p_{i,j-1}^{n+1}}{\Delta y^2}\right)+\left(\frac{2k}{\alpha\mu}\left(\frac{1}{\Delta x^2}+\frac{1}{\Delta y^2}\right)+\frac{S_\varepsilon}{\alpha \Delta t}+\left(\frac{1-2\nu}{1-\nu}\right)\frac{\alpha}{G\Delta t} \right)p_{i,j}^{n+1}= -\frac{1}{\Delta t}\left(\frac{1-2\nu}{1-\nu}\right)\left(\frac{u_{i+1,j}^{n}-u_{i-1,j}^{n}}{2\Delta x}\right)+\left(\frac{S_\varepsilon}{\alpha \Delta t} -\left(\frac{1-2\nu}{1-\nu}\right)\frac{\alpha}{G\Delta t}\right)p_{i,j}^n-\left(\frac{1-2\nu}{1-\nu}\right)\frac{1}{G\Delta t}(\sigma_{ice,i,j}^{n+1}+\sigma_{ice,i,j}^{n})$$

$(p_{i,j+1} - p_{i,j-1})/2\Delta y$ is adjusted based on the fluid boundary condition -- e.g., $p = 0$, $p = p_{subglacial}$ or $\partial p/ \partial y = 0$.

Note that the substitutions introduce terms with subscripts $i-2$ and $i+2$ that need to be handled carefully at the horizontal boundaries. 

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

