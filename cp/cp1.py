import numpy as np
import matplotlib.pyplot as plt


def analytical(r_list, T_0):
    length = len(r_list)
    temp = (T_0 / 2) * (1 - np.cos(np.pi * r_list / r_list[-1]))
    return temp

r_list = np.linspace(0, 1, 100)
T_0 = 20
analytical_temp = analytical(r_list, T_0)
plt.plot(r_list, analytical_temp, label='analytical')
plt.show()