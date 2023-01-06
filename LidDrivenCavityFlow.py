import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm


# Initialise all the constants needed inside 
N_POINTS = 41
DOMAIN_SIZE = 1.0
N_ITERATIONS = 5000
DELTA_T=1e-3
NU=0.000001
RHO=1.125 
U_TOP=1.0
N_PRESSURE_POISSON_ITERATIONS=50
CFL_number = 0.5

Re = (U_TOP*DOMAIN_SIZE)/NU
print("Reynolds Number = "+str(Re))

def main():
    delta_X = DOMAIN_SIZE/(N_POINTS-1)  # Since spacing is uniform, delta defined in only one direction

    # Uniform Grid on both x and y directions
    x = np.linspace(0.0,DOMAIN_SIZE,N_POINTS)
    y = np.linspace(0.0,DOMAIN_SIZE,N_POINTS)
    
    X,Y = np.meshgrid(x,y)   # Creating Meshgrid 

    # Defining Empty Arrays for Velocities with Zeros Initialised at all points

    u_prev = np.zeros_like(X)
    v_prev = np.zeros_like(X)
    p_prev = np.zeros_like(X)

    u_time = np.zeros(shape=(N_POINTS,N_POINTS,N_ITERATIONS))
    v_time = np.zeros_like(u_time)
    p_time = np.zeros_like(u_time)


    def central_difference_x(f,delta_X):
        diff = np.zeros_like(f)
        diff[1:-1,1:-1] = (f[1:-1,2: ] - f[1:-1,0:-2])/(2*delta_X)
        return diff

    def central_difference_y(f,delta_X):
        diff=np.zeros_like(f)
        diff[1:-1,1:-1] = (f[2: ,1:-1] - f[0:-2,1:-1])/(2*delta_X)
        return diff

    def laplace(f,delta_X):
        diff= np.zeros_like(f)
        diff[1:-1,1:-1] = (f[1:-1,0:-2] + f[0:-2,1:-1] - 4*f[1:-1,1:-1] + f[1:-1,2: ]+f[2:,1:-1])/delta_X**2
        return diff

    Delta_T_MAX = (
        0.5 * DELTA_T**2 / NU
    )
    if DELTA_T > CFL_number * Delta_T_MAX:
        print("Delta_T = "+ str(DELTA_T)+"Delta_T_max=",str(Delta_T_MAX))
        raise RuntimeError("Stability is not guarenteed")
        
    ii =0
    
    for __ in tqdm(range(N_ITERATIONS)):
        du_prev_dx = central_difference_x(u_prev,delta_X)
        du_prev_dy = central_difference_y(u_prev,delta_X)
        dv_prev_dx = central_difference_x(v_prev,delta_X)
        dv_prev_dy = central_difference_y(v_prev,delta_X)
        laplace_u_prev = laplace(u_prev,delta_X)
        laplace_v_prev = laplace(v_prev,delta_X)

        #Calculating Tentative Velocities without pressure gradients
        #by rearranging the momentum equation

        u_tent = DELTA_T*(-1*(u_prev * du_prev_dx + v_prev*du_prev_dy ) + NU * laplace_u_prev )

        v_tent = DELTA_T*(-1*(u_prev * dv_prev_dx +  v_prev*dv_prev_dy) + NU * laplace_v_prev )

        # Applying Velocity Boundary Conditions on all four sides of the cavity

        u_tent[0,:] = 0.0    
        u_tent[:,0] = 0.0    
        u_tent[:,-1] =0.0
        u_tent[-1,:] = U_TOP        
        v_tent[0,:] = 0.0
        v_tent[:,0] = 0.0
        v_tent[-1,:] =0.0
        v_tent[:,-1] = 0.0


        d_u_tent_d_x = central_difference_x(u_tent,delta_X)
        d_v_tent_d_y = central_difference_y(v_tent,delta_X)

        # Pressure Correction by solving the pressure poisson equation

        rhs = RHO/DELTA_T *(d_u_tent_d_x + d_v_tent_d_y)

        for __ in range(N_PRESSURE_POISSON_ITERATIONS):
            p_next = np.zeros_like(p_prev)
            p_next[1:-1,1:-1] = 0.25 * ( p_prev[1:-1,0:-2] + p_prev[0:-2,1:-1] + p_prev[1:-1,2:] + p_prev[2:,1:-1] - (delta_X**2) * rhs[1:-1,1:-1] )

            # Pressure Boundary Conditions - Homogenous Neumann Boundary Conditions 
            # Except for the top, where it is a Dirichlet BC

            p_next[:,-1] = p_next[:,-2]
            p_next[0,:]  = p_next[1,:]
            p_next[:,0]  = p_next[:, 1]
            p_next[-1,:] = 0.0

            p_prev = p_next
            
    
        d_p_next_d_x = central_difference_x(p_next,delta_X)
        d_p_next_d_y = central_difference_y(p_next,delta_X)

        
    
        # Adding Velocity Corrections 

        u_next = u_tent  -  (DELTA_T/RHO) * d_p_next_d_x

        v_next = v_tent  -  (DELTA_T/RHO) * d_p_next_d_y

   

    
        u_next[0,:] = 0.0    
        u_next[:,0] = 0.0    
        u_next[:,-1] =0.0
        u_next[-1,:] = U_TOP        
        v_next[0,:] = 0.0
        v_next[:,0] = 0.0
        v_next[-1,:] =0.0
        v_next[:,-1] = 0.0
    
        # Advance in time
        u_prev = u_next
        v_prev = v_next
        p_prev = p_next
        
        u_time[:,:,ii] = u_prev
        v_time[:,:,ii] = v_prev
        p_time[:,:,ii] = p_prev
        ii=ii+1

        

        


    # Plot Section
    plt.figure()
    plt.contourf(X,Y,p_next)
    plt.streamplot(X,Y,u_next,v_next,color="black")
    plt.colorbar
    plt.show()

    plt.figure()
    plt.quiver(X,Y,u_next,v_next,color="black")
    plt.streamplot(X,Y,u_next,v_next,color="black")
    plt.show()

    print(ii)
    plt.figure()
    for ii in range(0,50,1):
        plt.contourf(X,Y,p_time[:,:,ii])
        plt.pause(0.1)  
        plt.show()  






if __name__== "__main__":
    main()