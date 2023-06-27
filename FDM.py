import numpy as np

platelength = 50
platewidth  = 50
maxiter     = 700

alpha       = 2
deltax      = 1

deltat      = (deltax**2)/(4*alpha)
gamma       = (alpha * deltat)/(deltax**2)

#   Boundary conditions
utop    = 100
uleft   = 0.0
ubottom = 0.0
uright  = 0.0

#   Initial condition
uinit   = 0

def initializeu(maxiter, ni=platelength, nj=platewidth):
    
    #   Initialize solution: the grid of u(k,i,j)
    u   = np.full((maxiter,ni,nj),uinit)

    #   Set the boundary conditions
    u[:,0,:]    = utop
    u[:,:,0]    = utop
    u[:,-1,:]    = utop
    u[:,:,-1]    = utop

    return u

#   Original code to compute heatmap at each space and time step

def calculate(u):
    nk,ni,nj = u.shape
    for k in range (0,nk-1):
        for i in range (0,ni-1):
            for j in range (i,nj-1):
                u[k+1,i,j]  =   gamma * (u[k][i+1][j] + 
                                         i[k][i-1][j] + 
                                         u[k][i][j+1] +
                                         u[k][i][j-1] -
                                         4*u[k][i][j]) + u[k][i][j]
    return u

                                        
# calculate(initializeu(maxiter,ni=platelength,nj=platewidth))
# test