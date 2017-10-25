import numpy as np
import matplotlib.pyplot as plt
# constants

def hw_6_2(start, end, step):
    #HW 6 Problem 2

    #Steady State Equation
    maximum = 100
    c_val = 0
    for c in np.linspace(start, end, step):
        a= 5.
        l= 0.5
        t_0= 20.
        x=np.linspace(0,l,100)
        k = 1.28
        t_ss=-(c*np.exp(-a*x)/((a**2)*k))+(x/l)*(t_0 - (c*(1-np.exp(-a*l)))/((a**2)*k))+(c/((a**2)*k)) +273
        diff = abs(40 +273 - max(t_ss))
        plt.plot(x, t_ss, label=str(c))
        if diff < maximum:
            maximum = diff
            c_val = c
    print(c_val)
    plt.plot(x, np.linspace(313, 313, 100))
    max_temp_x = np.log(-(a*k)/(c*l)*(t_0-(c*(1-np.exp(-a*l)))/(a**2*k)))/(-a)
    print('Maximum Temperature Occurs at : %s' %str(max_temp_x))
    return c_val

def hw_6_2_b(c_val):
    l = 0.5
    a = 5
    k = 1.28
    t_0 = 20
    x=np.linspace(0,l,100)
    print('The c constant that makes the maximum temperature closest to 40C (313K) is:')
    print(c_val)
    c = c_val
    t_ss=-(c*np.exp(-a*x)/((a**2)*k))+(x/l)*(t_0 - (c*(1-np.exp(-a*l)))/((a**2)*k))+(c/((a**2)*k)) + 273
    plt.plot(x, t_ss, label='Analytical Solution')
    plt.ylabel('Temperature [K]')
    plt.xlabel('x [m]')
    plt.legend(loc=(1.0,0))


def hw_6_3(x_grid, c):

    l = 0.5
    a = 5
    k = 1.28
    t_0 = 20

    # space grid definition
    x = np.linspace(0,l,x_grid)
    dx = x[1]-x[0]
    # heat generation term
    q = c * np.exp(-a * x)

    heat_gen = [(qx*(dx**2))/(2*k) for qx in q]
    array = np.zeros(x_grid)
    oldarray = np.zeros(x_grid)
        
    # boundary conditions
    array[0] = 273
    array[-1] = t_0 + 273

    eps =1000
    iterations = 0
    error = 100
    # run loop while eps is small (converges or 10000 iterations)
    while True:
        print('Iteration ' + str(iterations))
        iterations += 1
        print('Error ' + str(error))
        print(x_grid)
        for space in range(1, x_grid-1):
            array[space] = ((array[space - 1] + array[space + 1])/2) + heat_gen[space]
        print(array)
        error = max(abs(array - oldarray))
        print(error)
#        error = max([abs(x-y) for x,y in zip(oldarray, array)])
        if error > 1e-3:
            oldarray = np.copy(array)
        else:
            break

    print('Number of Iterations it took:' + str(iterations))
    print('Epsilon: ' + str(eps))
    print('Convergence: %s' %str(error))
    print('The Temperature Profile is:')
    print(array)
    label = '%s max_temp = %s (%s)' %(str(dx), str(max(array)), str(abs(313-max(array))/max(array)) )
    plt.plot(x, array, label=label)



#################################
#################################
#################################
#################################

#hw_6_2(3700, 4000, 1000)
#plt.legend()
#axes = plt.gca()
#plt.ylabel('Temperature [K]')
#plt.xlabel('x [m]')
#plt.show()
hw_6_2_b(3767.5675)

for x_grid in [10, 30, 100, 300]:
   hw_6_3(x_grid, 3767.5675)

plt.legend()
plt.ylabel('Temperature [K]')
plt.xlabel('x [m]')
plt.show()
