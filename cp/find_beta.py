import numpy as np


def beta_equation(beta):
    y = 3*beta - np.tan(beta * 3)
    return y 

R = 3
k = 0.15
rho = 8000. / 1000000.
cp = 500
alpha = k / (rho*cp)
T_0 = 1

points = np.arange(0, 1e2, 1e-5)
pot_betas = []
count = 0
flag = 1
for beta in points:
    if beta_equation(beta) > 0:
        flag = 0
    elif beta_equation(beta) < 0 and flag == 0:
        pot_betas.append(beta)
        count = count + 1
        flag = 1
betas = pot_betas



file = open('betas', 'w')
count = 0
for b in betas[:-1]:
    count = count + 1
    print(b) 
    file.write(str(count) + ' & ' + str(b) + '\\\ \n')

