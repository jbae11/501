import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import scipy.special as special
import scipy.optimize
from iapws import IAPWS97
from sympy import *

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
    z = np.linspace(0, 366, 100)
    t_c = []
    for i in z:
        root = np.roots([-1.47e-6, 0.0025, -1.41, 276.23-np.cos(np.pi * i / 366)])
        print(root)
        print(type(root))
        filtered_root = float(str(root[0])[1:7])
        t_c.append(filtered_root)

    plt.plot(t_c, z)
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

def find_tc():
    z = np.linspace(0, 366, 1000)
    t_c = []
    rho_c = []
    for i in z:
        kewl = 581.02-18.011*np.cos(np.pi*i / 366)
        t_c.append(kewl)

    return t_c


def find_h():
    ri = 0.47
    ro = 0.625
    v = 350
    z = np.linspace(0, 366, 1000)
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

    plt.plot(h_list, z)
    plt.xlabel('Heat Transfer Coefficient [W/cm^2 k]')
    plt.ylabel('z [cm]')
    plt.title('Heat Transfer Coefficient vs z')
    plt.savefig('h_z.png', format='png')
    plt.show() 

    return h_list

def find_c2():
    q_vol_list = [116.03, 164.1, 116.03]
    h_list = [3.591, 3.575, 3.567]
    t_list = [568.22, 581.04, 593.78]
    L_list = [366/4, 366/2, 366*3/4]
    for i in range(0, len(h_list)):
        c_2 = q_vol_list[i] * (4.122/h_list[i] + .96) + t_list[i]
        r = np.linspace(0, 0.47, 100)
        t_r = (-q_vol_list[i]*r**2)/(4*0.057) + c_2
        plt.plot(r, t_r, label='z = %s cm'%str(L_list[i])[:5])

    plt.legend()
    plt.xlabel('Radial Distance [cm]')
    plt.ylabel('Fuel Temperature [K]')
    plt.title('Fuel Temperature vsRadial Distance')
    plt.savefig('t_f_r.png', format='png')
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


########################################################

########################################################

########################################################

########################################################

########################################################

########################################################

find_h()