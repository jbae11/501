import numpy as np
import matplotlib.pyplot as plt


#GLOBAL VAIRABLE CONSTANTS
# m, J, s, kg
R = 0.03
k = 15
rho = 8000
cp = 500
alpha = 3.75e-6
T_0 = 1

def analytical_ss (r_list):
    temp = (T_0 / 2) * (1 - np.cos(np.pi * r_list / r_list[-1]))
    return temp


def numerical (r_grid, t_max, t_grid):
    # GRID AND STUFF 
    r_list = np.linspace(0, R, r_grid)
    dr = r_list[1] - r_list[0]
    t_list = np.linspace(0, t_max, t_grid)
    dt = t_list[1] - t_list[0]

    t = np.ndarray(shape=(r_grid, t_grid), dtype=float)
    told = np.ndarray(shape=(r_grid, t_grid), dtype=float)

    # Apply Initial Condition
    t[0] = (T_0 / 2) * (1 - np.cos(np.pi * r_list / r_list[-1]))

    change = 1000
    iterations = 0
    while change > 1e-3:
        print('Iteration %i' %iterations)
        print('Change %f' %change)
        iterations += 1
        for timestep in range(1, t_grid-1):
            for radius in range(1, r_grid-1):
                array[radius, timestep] = (1 / alpha)*((t[radius, timestep+1] - t[radius, timestep]) / dt )










analytical_temp = analytical_ss(r_list)
plt.plot(r_list, analytical_temp, label='analytical')
plt.show()
