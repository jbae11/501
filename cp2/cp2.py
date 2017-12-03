import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import scipy.special as special
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
        print(x)
        integral = integrate.quad(x, a=290+273, b=325+273)
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
    z = np.linspace(0, 366, 100)
    t_c = []
    for i in z:
        root = np.roots([-1.47e-6, 0.0025, -1.41, 265-np.cos(np.pi * i / 366)])
        print(root)
        print(type(root))
        filtered_root = float(str(root[0])[1:7])
        t_c.append(filtered_root)

    plt.plot(t_c, z)
    plt.show()

def low_key():
    z = np.linspace(0, 366, 1000)
    t_c = []
    for i in z:
        temp = 18.02*(1-np.cos(np.pi*i/366)) + 563
        t_c.append(temp)

    plt.plot(z, t_c)
    plt.show()

def shit():
    z = np.linspace(0, 366, 1000)
    t_c = []
    rho_c = []
    for i in z:
        kewl = 581.02-18.011*np.cos(np.pi*i / 366)
        t_c.append(kewl)
        rhoc = IAPWS97(T=kewl, P=15.17).rho /1000
        rho_c.append(rhoc)

    plt.plot(rho_c, z)
    plt.xlabel('Density of Water [g/cm^3]')
    plt.ylabel('z [cm]')
    plt.title('Water Density vs z')
    plt.savefig('rho_c_z.png', format='png')
    plt.show() 

"""
    plt.plot(t_c, z)
    plt.xlabel('Coolant Temperature [K]')
    plt.ylabel('z [cm]')
    plt.title('Coolant temperature vs z')
    plt.savefig('t_c_z.png', format='png')
    plt.show()
"""


def find_tc(z_grid):
    z = np.linspace(0, 366, z_grid)
    t_c = []
    rho_c = []
    for i in z:
        kewl = 581.02-18.011*np.cos(np.pi*i / 366)
        t_c.append(kewl)

    return t_c


def find_h(z_grid):
    ri = 0.47
    ro = 0.625
    v = 350
    z = np.linspace(0, 366, z_grid)
    t_c = []
    h_list = []
    for i in z:
        temp = 581.02-18.011*np.cos(np.pi*i /366)
        t_c.append(temp)
    for i in t_c:
        L = 4 * np.pi*(ro**2 - ri**2) / (2*np.pi*ro + 2*np.pi*ri)
        mu = IAPWS97(T=i, P=15.17).mu * 10
        rho = IAPWS97(T=i, P=15.17).rho /1000
        cp = IAPWS97(T=i, P=15.17).cp
        k = IAPWS97(T=i, P=15.17).k /100
        Re = rho*v*L / mu
        Pr = mu*cp / k
        Nu = 0.023 * Re**(0.8) * Pr**(0.4)
        h = Nu * k / L
        h_list.append(h)


    return h_list

def q_vol_z():
    z = np.linspace(0, 366, 1000)
    q = 164.1 * np.sin(np.pi*z / 366)
    plt.plot(z, q)
    plt.xlabel('Height [cm]')
    plt.ylabel('Volumetric Heat Generation [W/m^3]')
    plt.title('Volumetric Heat Generation with Respect to Height')
    plt.savefig('q_vol_z.png', format='png')
    plt.show()
    plt.close()

def four():
    R = 0.47
    k = 0.057
    L = 366
    q_vol_list = [116.03, 164.1, 116.03]
    h_list = [3.591, 3.575, 3.567]
    t_list = [568.22, 581.04, 593.78]
    L_list = [L/4, L/2, L*3/4]
    labels = ['L/4', 'L/2', '3L/4']
    for i in range(0, len(h_list)):
        c_2 = (q_vol_list[i]*R/2) + q_vol_list[i] * h_list[i] * R**2 /(4*k) + h_list[i] * t_list[i]
        r = np.linspace(0, 0.47, 100)
        t_r = (-q_vol_list[i]*r**2)/(4*k) + c_2
        plt.plot(r, t_r, label='z = ' + str(L_list[i])[:5] + ' cm (%s)' %labels[i])

    plt.legend()
    plt.xlabel('Radial Distance [cm]')
    plt.ylabel('Fuel Temperature [K]')
    plt.title('Fuel Temperature vs Radial Distance')
    plt.savefig('t_f_r.png', format='png')
    plt.show()
    plt.close()

def fuel_rod_temp(z_grid):
    """ temperature distribution in the fuel rod T_f(r,z)"""


    tol = 1e-3

    # k of UOX at 300C:
    k = 0.057 
    change = 10000
    z_list = np.linspace(0,366, z_grid)
    dz = z_list[1]-z_list[0]
    q_vol = 164.1*np.sin(np.pi * z_list / 366)
    r_list = np.linspace(0, 0.47, 50)
    dr = r_list[1]-r_list[0]
    h = find_h(z_grid)
    t_c = find_tc(z_grid)

    # 2d matrix with len(z) columns and len(r) rows
    t = np.zeros(shape=(len(z_list),len(r_list)), dtype=float)
    print(len(t))
    print(len(t[0]))
    # insulated bc at z=0 and z=L
    # convective bc at r = 0.47
    # r=0 bc

    trold = np.zeros(len(r_list))
    while True:
        told = np.copy(t)
        for z in range(1,len(z_list)-1):
            while True:
                tr = t[z][:]
                trold = tr.copy()
                for r in range(1,len(r_list)-1):
                    # z BC at the end
                    if z == len(z_list):
                        one = (2*t[z-1][r])/(dz**2)
                        two = (t[z][r+1] - t[z][r-1]) / (r_list[r] * 2*dr)
                        three = (t[z][r-1] + t[z][r+1]) / (dr**2)
                        four = q_vol[z] / k
                        denom = 2/(dz**2) + 2/(dr**2) 
                    # z BC at the beginning
                    elif z == 0:
                        one = (2*t[z+1][r])/(dz**2)
                        two = (t[z][r+1] - t[z][r-1]) / (r_list[r] * 2*dr)
                        three = (t[z][r-1] + t[z][r+1]) / (dr**2)
                        four = q_vol[z] / k
                        denom = 2/(dz**2) + 2/(dr**2) 
                    # all other scenarios
                    else:
                        one = (2*t[z+1][r])/(dz**2)
                        two = (t[z][r+1] - t[z][r-1]) / (r_list[r] * 2*dr)
                        three = (t[z][r-1] + t[z][r+1]) / (dr**2)
                        four = q_vol[z] / k
                        denom = 2/(dz**2) + 2/(dr**2) 

                    # do the addition to find t[z][r]
                    t[z][r] = (one + two + three + four) / denom
                # BC at r=0
                t[z][0] = t[z][1]
                #convective BC:
                t[z][-1] = (t[z][-2] + dr* (h[z]/k) * t_c[z]) / (1+ dr*h[z]/k)
                # (-h[z] * t_c[z] - (k*t[z][-2]/dr)) / (k/dr - h[z])
                print('iteration')
                conv = abs(tr-trold)
                if max(conv) < tol:
                    break

        conver = np.abs(t-told)
        if np.max(conver) < tol:
            break
        print('ziteration')

    plt.plot(r_list, t[int(z_grid/4)], label='L/4')

    plt.plot(r_list, t[int(z_grid/2)], label='L/2')

    plt.plot(r_list, t[int(3*z_grid/4)], label='3L/4')
    plt.legend()
    plt.show()


########################################################

########################################################

########################################################

########################################################

########################################################

########################################################
fuel_rod_temp(10)