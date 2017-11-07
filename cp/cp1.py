import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import scipy.special as special
import scipy.optimize


def steady_state (r_grid):
    r_list = np.linspace(0, 3, r_grid)
    temp = (3/2) * ((np.pi**2) + 6) / (3 * (np.pi**2))
    t = [temp] * r_grid

    plt.plot(r_list, t, label = 'Steady State Temperature')
    plt.ylabel(r'$\frac{Temperature}{T_0}$')
    plt.xlabel('r [cm]')
    plt.legend()
    plt.title('Steady State Temperature of the Sphere')
    plt.savefig('steady_state.png', format='png')
    plt.show()

def numerical (r_grid, t_max):
    # GRID AND STUFF 
    R = 3
    k = 0.15
    rho = 8000. / 1000000.
    cp = 500
    alpha = k / (rho*cp)
    T_0 = 1

    t_grid = 10000

    r_list = np.linspace(0, R, r_grid) 
    dr = r_list[1] - r_list[0]
    t_list = np.linspace(0, t_max, t_grid)
    dt = t_list[1] - t_list[0]

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
    plt.legend()
    plt.title('Numerical Solution of Time-evolution of Heat Profile')
    plt.savefig('Numerical.png', format='png')
    plt.show()

def beta_equation(beta):
    y = 3*beta - np.tan(beta * 3)
    return y 

def analytical (r_grid, t_max, t_grid, n):
    R = 3
    k = 0.15
    rho = 8000. / 1000000.
    cp = 500
    alpha = k / (rho*cp)
    T_0 = 1

    #########################
    # finding betas
    points = np.arange(-0.1, 1e3, 1e-4)
    pot_betas = []
    for beta in points:
        if beta_equation(beta) > 0:
            flag = 0
        elif beta_equation(beta) < 0 and flag == 0:
            pot_betas.append(beta)
            flag = 1

    # pick only first n betas
    betas = pot_betas[:n]
    print(betas)
    ###########################

    T_0 = 1
    r_list = np.linspace(0, R, r_grid)
    dr = r_list[1] - r_list[0]
    t_list = np.linspace(0, t_max, t_grid)
    dt = t_list[1] - t_list[0]

    t = [0] * r_grid
    t_compiled_betas = [0] * r_grid
    t_tot = (T_0 / 2) * (1 - np.cos(np.pi * r_list / R))
    t_ss = 3*T_0 /2 * ((np.pi**2 + 6)/(3*np.pi**2))
    print(t_ss)
    steady_state = [t_ss] * r_grid
    # plot initial condition
    plt.plot(r_list, t_tot, label = 't = 0 secs')
    plt.plot(r_list, steady_state, label = 'Steady State')

    def integrand1(r, b, T_0, R, t_ss):
        return ((r*T_0*0.5*(1-np.cos(np.pi*r / R)) - t_ss * r)*np.sin(b * r))
    def integrand2(r, b):
        return ((np.sin(b*r))**2) 

    print(betas)
    for time in range(1, len(t_list)):
        t_compiled_betas = [0] * r_grid
        for b in betas[1:]:
            t = [0] * r_grid
            for space in range(1, len(r_list)):
                r = r_list[space]
                integral1 = integrate.quad(integrand1, 0, R, args=(b, T_0, R, t_ss))
                integral2 = integrate.quad(integrand2, 0, R, args=(b))
                a = integral1[0] / integral2[0]
                spatial = np.sin(b*r) / r
                temporal = np.exp(-1 * b**2 * alpha * t_list[time])
                t[space] = a * spatial * temporal
                t[0] = t[1]
            t_compiled_betas = [x+y for x,y in zip(t_compiled_betas, t)]

        t_compiled_betas = [x+t_ss for x in t_compiled_betas]    
        plt.plot(r_list, t_compiled_betas, label='t = %s secs' %str(t_list[time])[:5])
        t_tot = np.vstack((t_tot, t_compiled_betas))

    print(t_tot)
    plt.ylabel(r'$\frac{Temperature}{T_0}$')
    plt.xlabel('r [cm]')
    plt.legend()
    plt.title('Analytical Solution of Time-evolution of Heat Profile')
    plt.savefig('Analytical.png', format='png')
    plt.show()



def anal_n_diff (r_grid, time, n):
    R = 3
    k = 0.15
    rho = 8000. / 1000000.
    cp = 500
    alpha = k / (rho*cp)
    T_0 = 1

    #########################
    # finding betas
    points = np.arange(-0.1, 1e3, 1e-4)
    pot_betas = []
    count = 0
    for beta in points:
        if beta_equation(beta) > 0:
            flag = 0
        elif beta_equation(beta) < 0 and flag == 0:
            pot_betas.append(beta)
            count = count + 1
            flag = 1
        if count > n:
            break
    betas = pot_betas
    ###########################

    T_0 = 1
    r_list = np.linspace(0, R, r_grid)
    dr = r_list[1] - r_list[0]

    t = [0] * r_grid
    t_compiled_betas = [0] * r_grid
    t_ss = 3*T_0 /2 * ((np.pi**2 + 6)/(3*np.pi**2))
    def integrand1(r, b, T_0, R, t_ss):
        return ((r*T_0*0.5*(1-np.cos(np.pi*r / R)) - t_ss * r)*np.sin(b * r))
    def integrand2(r, b):
        return ((np.sin(b*r))**2) 

    t_compiled_betas = [0] * r_grid
    t_old = [0] * r_grid

    for i in range(1,len(betas)):
        t = [0] * r_grid
        for space in range(1, len(r_list)):
            r = r_list[space]
            integral1 = integrate.quad(integrand1, 0, R, args=(betas[i], T_0, R, t_ss))
            integral2 = integrate.quad(integrand2, 0, R, args=(betas[i]))
            a = integral1[0] / integral2[0]
            spatial = np.sin(betas[i]*r) / r
            temporal = np.exp(-1 * betas[i]**2 * alpha * time)
            t[space] = a * spatial * temporal
            t[0] = t[1]
        t_compiled_betas = [x+y for x,y in zip(t_compiled_betas, t)]
        t_compiled_betas = [x+t_ss for x in t_compiled_betas]
        
        if sum([abs(x-y) for x, y in zip(t_compiled_betas, t_old)]) < 1e-3:
            print('Error Less than 1e-4 in n = %s' %str(i))    
        plt.plot(r_list, t_compiled_betas, label='n = %s terms' %str(i))
        t_old = np.copy(t_compiled_betas)
        t_compiled_betas = [x-t_ss for x in t_compiled_betas]
        

    plt.ylabel(r'$\frac{Temperature}{T_0}$')
    plt.xlabel('r [cm]')
    plt.legend()
    plt.title('Terms in the Series Solution for Convergence. t = %s secs' %str(time))
    plt.savefig('terms_%s.png' %str(time), format='png')
    plt.show()


####################################

#steady_state(100)
#analytical(100, 70, 10, 20)
#numerical(100, 70)
anal_n_diff (100, 4, 15)