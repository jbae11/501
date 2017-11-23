import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import scipy.special as special
import scipy.optimize
from iapws import IAPWS97

def fit_poly_cp(order_to):
    """ fits a polynomial for Temp and C_p of Water"""
    temp = np.linspace(290+273, 320+273, 100)
    cp =[]
    for i in temp:
        cp.append(IAPWS97(T=i, P=15.17).cp)
    plt.plot(temp, cp, 'ro', label='Data from IAPWS97' )

    for i in range(1, order_to+1):
        eq = np.polyfit(temp, cp, i)
        print(eq)
        eq = list(reversed(eq))
        print(eq)
        x = np.polynomial.polynomial.Polynomial(eq)
        integral = integrate.quad(x, 290+273, 325+273)
        x = list(x)
        x = list(reversed(x))
        print(x)
        y = np.polyval(x, temp)
        
        plt.plot(temp, y, label='Order %s'%str(i))
        plt.xlim(290+273, 325+273)
        print('Result of the Integral for order %s is:' %str(i))
        print(integral[0])
        print('C = ')
        print(132.446*integral[0] / (366*2*.47**2))
        print('\n \n')

    plt.xlabel('T [K]')
    plt.ylabel('C_p [J/gK]')
    plt.legend()
    plt.title('C_p value vs fit')
    plt.savefig('cp_plot.png', format='png')
    plt.show()


def root_solver():
    """ fits a polynomial for Temp and C_p of Water"""
    z = np.linspace(0, 366, 10600)
    t_c = []
    for i in z:
        root = np.roots([-1.47e-6, 0.0025, -1.41, 276-np.cos(np.pi * i / 366)])
        print(root)
        filtered_root = float(str(root[0])[1:-4])
        t_c.append(filtered_root)
    plt.plot(t_c, z)
    plt.show()


def fuel_rod_temp():
    """ temperature distribution in the fuel rod T_f(r,z)"""

    # k of UOX at 300C:
    k = 5.7 
    change = 10000
    zee = np.linspace(0,366, 10)
    dz = zee[1]-zee[0]
    q_vol = 164.1*np.sin(np.pi * zee / 366)
    ar = np.linspace(0, 0.47, 50)
    dr = ar[1]-ar[0]
    t = np.ndarray(shape=(len(zee),len(ar)), dtype=float)
    # insulated bc at z=0 and z=L
    for r in range(0,len(ar)):
        t[0][r] = 563
        t[1][r] = 563
        t[-1][r] = 598
        t[-2][r] = 598


    trold = np.zeros(len(ar))
    while zchange > 1e-3:
        for z in range(2,len(zee)-1):
            while change > 1e-3:
                for r in range(2,len(ar)-2):
                    one = - (t[z-1][r] + t[z+1][r])/(dz**2)
                    two = t[z][r-1]/(ar[r] * dr)
                    three = - (t[z][r-1] + t[z][r+1]) / (dr**2)
                    four = - q_vol[z] / k
                    denom = -2/(dz**2) -2/(dr**2) + 1/(ar[r]*dr)
                    t[z][r] = (one + two + three + four) / denom

                # BC at r=0
                t[z][0] = t[z][1]
                tr = t[z][:]
                print('iteration')
                change = max(tr-trold)
                trold = tr
            one = - (t[z-1][r] + t[z+1][r])/(dz**2)
            two = t[z][r-1]/(ar[r] * dr)
            three = - (t[z][r-1] + t[z][r+1]) / (dr**2)
            four = - q_vol[z] / k
            denom = -2/(dz**2) -2/(dr**2) + 1/(ar[r]*dr)
            t[z][r] = (one + two + three + four) / denom    
            print('ziteration')
    print(t)
def numerical (r_grid, t_max, t_grid):
    # GRID AND STUFF 
    R = 3
    k = 0.15
    rho = 8000. / 1000000.
    cp = 500
    alpha = k / (rho*cp)
    T_0 = 1

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

        if timestep%1000 ==0:
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

def gimme_a_t (r_grid, t_list):
    R = 3
    k = 0.15
    rho = 8000. / 1000000.
    cp = 500
    alpha = k / (rho*cp)
    T_0 = 1
    n = 30
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

    for time in t_list:
        t_compiled_betas = [0] * r_grid
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
        plt.plot(r_list, t_compiled_betas, label='t = %s secs' %str(time)[:5])
        

    plt.ylabel(r'$\frac{Temperature}{T_0}$')
    plt.xlabel('r [cm]')
    plt.legend()
    plt.title('converged T(r,t) for all Seven Cases')
    plt.savefig('converged.png', format='png')
    plt.show()

def anal_champ (r_grid, t_list):
    R = 3
    k = 0.15
    rho = 8000. / 1000000.
    cp = 500
    alpha = k / (rho*cp)
    T_0 = 1
    n = 20
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

    for time in t_list:
        t_compiled_betas = [0] * r_grid
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
        plt.plot(r_list, t_compiled_betas, label='Analytical Solution' )

####################################


def gridd(r_grid_list, t_max, time):
    # GRID AND STUFF 
    R = 3
    k = 0.15
    rho = 8000. / 1000000.
    cp = 500
    alpha = k / (rho*cp)
    T_0 = 1

    t_grid = 10000
    
    for r_grid in r_grid_list:
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

        for timestep in range(1, len(t_list)):
            t[0,timestep] = t[0,timestep-1] + (alpha * (x/dr) * (2*t[1,timestep-1] - 2*t[0,timestep-1]))
            t[-1,timestep] = t[-1, timestep-1] + alpha * (x/dr) * (2*t[-2,timestep-1] - 2*t[-1,timestep-1])
            for space in range(1, len(r_list)-1):
                first_term = alpha * (x/r_list[space]) * (t[space+1,timestep-1] - t[space-1, timestep-1])
                second_term = alpha * (x/dr) * (t[space+1, timestep-1] - (2*t[space, timestep-1]) + t[space-1, timestep-1])
                t[space, timestep] = t[space, timestep-1] + first_term + second_term

            if abs(t_list[timestep] - time) < 1e-2:
                plt.plot(r_list, t[:,timestep], label='grid = '+ str(r_grid) + ' (dr = %s cm)' %str(dr)[:5])
                break


def num_champ(r_grid_list, t_max, time):
    # GRID AND STUFF 
    R = 3
    k = 0.15
    rho = 8000. / 1000000.
    cp = 500
    alpha = k / (rho*cp)
    T_0 = 1

    t_grid = 10000
    
    for r_grid in r_grid_list:
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

        for timestep in range(1, len(t_list)):
            t[0,timestep] = t[0,timestep-1] + (alpha * (x/dr) * (2*t[1,timestep-1] - 2*t[0,timestep-1]))
            t[-1,timestep] = t[-1, timestep-1] + alpha * (x/dr) * (2*t[-2,timestep-1] - 2*t[-1,timestep-1])
            for space in range(1, len(r_list)-1):
                first_term = alpha * (x/r_list[space]) * (t[space+1,timestep-1] - t[space-1, timestep-1])
                second_term = alpha * (x/dr) * (t[space+1, timestep-1] - (2*t[space, timestep-1]) + t[space-1, timestep-1])
                t[space, timestep] = t[space, timestep-1] + first_term + second_term

            if abs(t_list[timestep] - time) < 1e-2:
                plt.plot(r_list, t[:,timestep], label='Numerical Solution')
                break

def hagridd (t_grid, t_max=70):
    # GRID AND STUFF 
    R = 3
    k = 0.15
    rho = 8000. / 1000000.
    cp = 500
    alpha = k / (rho*cp)
    T_0 = 1

    t_list = np.linspace(0, t_max, t_grid)
    dt = t_list[1] - t_list[0]
    r_list = np.arange(0, R, (7/10)*dt) 
    dr = r_list[1] - r_list[0]


    if dt/dr > 1e8:
        raise ValueError('Too high of dt/dr dude')  

    t = np.zeros((len(r_list), len(t_list)), dtype=float)

    # Apply Initial Condition
    t[:,0] = (T_0 / 2) * (1 - np.cos(np.pi * r_list / R))
    x = dt/dr

    for timestep in range(1, len(t_list)):
        t[0,timestep] = t[0,timestep-1] + (alpha * (x/dr) * (2*t[1,timestep-1] - 2*t[0,timestep-1]))
        t[-1,timestep] = t[-1, timestep-1] + alpha * (x/dr) * (2*t[-2,timestep-1] - 2*t[-1,timestep-1])
        for space in range(1, len(r_list)-1):
            first_term = alpha * (x/r_list[space]) * (t[space+1,timestep-1] - t[space-1, timestep-1])
            second_term = alpha * (x/dr) * (t[space+1, timestep-1] - (2*t[space, timestep-1]) + t[space-1, timestep-1])
            t[space, timestep] = t[space, timestep-1] + first_term + second_term
        if timestep % max([(timestep/10),1]) == 0:
            plt.plot(r_list, t[:,timestep], label='t = %s secs' %str(t_list[timestep])[:5])

    plt.ylabel(r'$\frac{Temperature}{T_0}$')
    plt.xlabel('r [cm]')
    plt.legend()
    plt.title('Temporal Grid Refinement Study at t_grid = ' + str(t_grid) + ' (dt = %s)' %str(dt)[:4])
    plt.savefig('t_grid_%s.png' %t_grid, format='png')
    plt.show()


fuel_rod_temp()