## horizontal_diffusion.py
# This function calculates and plots the horizontal 1D diffusion of pore pressure for a
# zero-flux boundary condition on the left and a fixed p = 0 boundary condition on the right.

from matplotlib import rcParams
rcParams['font.family'] = 'Helvetica'
rcParams['font.size'] = 13
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as mcolors
from matplotlib.cm import ScalarMappable
import numpy as np


def horizontal_diffusion(nx,dx,x,times,dt,nt,gamma,c,sigma_ice,k,mu,B,nu,nu_u,G,alpha):

    # Set up inverse matrix
    A = np.zeros((nx, nx))

    p = np.zeros(nx)
    q = np.zeros(nt-1)
    p_all = np.zeros((nx,nt-1))      
    eps_zz = np.zeros((nx,nt-1))  
    eps_xx = np.zeros((nx,nt-1))  
    u = np.zeros((nx,nt-1))

    factor = (k/mu)*(2*G/(alpha*(1-nu)))

    for j in range(nx):
        if j==nx-1: # Fixed pore pressure boundary condition on the right
            A[j,j] = 1 
        elif j==0: # Zero-flux boundary condition on the left (for symmetry)
            A[j,j+1] = -2*factor/dx**2
            A[j,j] = 2*factor/dx**2+(3/(B*(1+nu_u)))/dt
        else:
            A[j,j+1] = -factor/dx**2
            A[j,j-1] = -factor/dx**2
            A[j,j] = 2*factor/dx**2+(3/(B*(1+nu_u)))/dt

    # Main time loop
    for n in range(nt-1):
        b = np.zeros(nx)
        for j in range(nx):
            if j == nx - 1: 
                b[j] = 0 # Zero pressure boundary condition on the right
            else:
                b[j] = (3/(B*(1+nu_u)))/(dt)*p[j]-(sigma_ice[j,n+1]-sigma_ice[j,n])/dt
        
        p = np.linalg.solve(A,b) # Solve for pore pressure
        
        p_all[:,n] = p
        eps_zz[:,n] = 1/(2*G)*((1-nu)*sigma_ice[:,n+1]+(1-2*nu)*alpha*p_all[:,n]) # Compute vertical strain
        eps_xx[:,n] = 1/(2*G)*(-nu*sigma_ice[:,n+1]+(1-2*nu)*alpha*p_all[:,n]) # Compute horizontal strain
        q[n] = -(k/mu)*(p[nx-1]-p[nx-2])/dx # Compute discharge 

    # Integrate depth-constant eps_zz to get vertical displacement
    v = eps_zz*1000 # Assume 1 km thick aquifer

    # Integrate eps_xx to get horizontal displacement
    for n in range(nt-1):
        u[:,n] = np.cumsum(eps_xx[:,n]*dx)

    # Plots
    colors = cm.get_cmap('YlOrRd', nt-1)
    plt.figure(figsize=(9,5.5))

    # Plot pore pressure profiles over time
    plt.subplot(2,2,1)
    for n in range(nt-1):
        if n == 0 or n % 10 == 0: 
            plt.plot(x/1e3,p_all[:,n]/1e6,'.-', markersize=3,color=colors(n))
    plt.xlabel('x [km]')
    plt.ylabel('Pore pressure, $p$ [MPa]')
    plt.xlim(-5,5)

    # Plot vertical displacement profiles over time
    plt.subplot(2,2,2)
    for n in range(nt-1):
        if n == 0 or n % 10 == 0: 
            plt.plot(x/1e3,v[:,n]*1e3, '.-', markersize=3,color=colors(n))
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.ylabel('Vertical displacement, $v$ [mm]')
    plt.xlim(-5,5)
    plt.xlabel('x [km]')

    # Plot horizontal discharge rate over time
    plt.subplot(2,2,3)
    for n in range(nt-1):
        plt.plot(times[n]/(3600*24),q[n]*(3600*24),'.',markersize=3,color=colors(n))
    plt.xlabel('Time [yr]')
    plt.ylabel('Horizontal discharge rate, $q$ [m/day]')
    plt.xlim(0,150)

    # Plot horizontal displacement profiles over time
    plt.subplot(2,2,4)
    for n in range(nt-1):
        if n == 0 or n % 10 == 0: 
            plt.plot(x/1e3,u[:,n]*1e3, '.-', markersize=3,color=colors(n))
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.ylabel('Horizontal displacement, $u$ [mm]')
    plt.xlabel('x [km]')
    plt.xlim(-5,5)

    plt.tight_layout()
    plt.show()
    