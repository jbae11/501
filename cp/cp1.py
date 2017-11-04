import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import scipy.special as special


def analytical_ss (r_list):
    temp = (T_0 / 2) * (1 - np.cos(np.pi * r_list / r_list[-1]))
    return temp


def numerical (r_grid, t_max):
    # GRID AND STUFF 
    R = 3
    k = 0.15
    rho = 8000. / 1000000.
    cp = 500
    alpha = k / (rho*cp)
    T_0 = 1

    t_grid = 1000

    r_list = np.linspace(0, R, r_grid) 
    dr = r_list[1] - r_list[0]
    print(r_list)
    t_list = np.linspace(0, t_max, t_grid)
    dt = t_list[1] - t_list[0]
    print(t_list)

    if dt/dr > 1e8:
        raise ValueError('Too high of dt/dr dude')  

    t = np.zeros((len(r_list), len(t_list)), dtype=float)

    # Apply Initial Condition
    t[:,0] = (T_0 / 2) * (1 - np.cos(np.pi * r_list / R))
    print(t[:,0])

    x = dt/dr
    
    plt.plot(r_list, t[:,0], label='t = 0 secs')

    for timestep in range(1, len(t_list)):
        t[0,timestep] = t[0,timestep-1] + (alpha * (x/dr) * (2*t[1,timestep-1] - 2*t[0,timestep-1]))
        t[-1,timestep] = t[-1, timestep-1] + alpha * (x/dr) * (2*t[-2,timestep-1] - 2*t[-1,timestep-1])
        for space in range(1, len(r_list)-1):
            first_term = alpha * (x/r_list[space]) * (t[space+1,timestep-1] - t[space-1, timestep-1])
            second_term = alpha * (x/dr) * (t[space+1, timestep-1] - (2*t[space, timestep-1]) + t[space-1, timestep-1])
            t[space, timestep] = t[space, timestep-1] + first_term + second_term

        if timestep%100 ==0:
            plt.plot(r_list, t[:,timestep], label='t = %s secs' %str(t_list[timestep])[:5])

    print(t)
    plt.ylabel(r'$\frac{Temperature}{T_0}$')
    plt.xlabel('r [cm]')
    plt.legend(loc=(1.0,0))
    plt.title('Numerical Solution of Time-evolution of Heat Profile')
    plt.savefig('Numerical', format='png')
    plt.show()

def analytical (r_grid, t_max, t_grid, n):
    R = 0.03
    k = 15
    rho = 8000
    cp = 500
    alpha = 3.75e-6
    T_0 = 1
    r_list = np.linspace(0, R, r_grid)
    dr = r_list[1] - r_list[0]
    t_list = np.linspace(0, t_max, t_grid)
    dt = t_list[1] - t_list[0]

    # create list of betas for n:
    betas = []
    for i in range(1, n):
        betas.append(i* np.pi / R)

    t = [0] * r_grid
    t_compiled_betas = []
    t_tot = np.ndarray(shape=(r_grid, t_grid), dtype=float)

    for time in range(0, len(t_list)):
        for b in betas:
            for space in range(1, len(r_list)-1):
                integral = integrate.quad(integrand1, 0, R, args=(T_0, b, R))
                integral2 = integrate.quad(integrand2, 0, R, args=(b))
                #print(integral)
                #print(integral2)
                spatial_term = integral[0] / (r_list[space] * integral2[0])
                print(spatial_term)
                #integrand = (T_0 * r_list[space]**3 * np.sin(b*r_list[space]) / 2) * (1-np.cos(np.pi*r_list[space]/R))
                #integrand2 = r_list[space]**2 * (np.sin(b*r_list[space]))**2
                #spatial_term = integrate.quad(integrand,0,R) / (r_list[space] * integrate.quad( r_list[space] * integrand2, 0, R ))
                temporal_term = np.exp(-alpha * b**2 * t_list[time]) 
                print(temporal_term)
                print(space)
                print(r_list[space])
                t[space] = spatial_term * temporal_term * np.sin(b*r_list[space])
                print(t)
            t_compiled_betas = [x+y for x,y in zip(t_compiled_betas, t)]
        t_tot = np.vstack((t_tot, t_compiled_betas))

    print(t_tot)

def integrand1(r,b,R,T_0):
    return T_0 * r**3 * np.sin(b*r) / 2 * (1-np.cos(np.pi * r/R))
def integrand2(r,b):
    return r**2 * (np.sin(b*r)**2)

numerical(25, 70)
