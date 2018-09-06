import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0, 100, 1000)

H = 50
delta = 2.5

T_base = 2
T_star = 1

T = T_star + 0.5*(T_base - T_star)*(1.0 + np.tanh((x - (H) - delta)/(0.5*delta)))

plt.plot(x, T)
plt.plot([H, H], [0, T_base])
#plt.plot([H-delta, H-delta], [0, T_base])
plt.savefig("profile2.png")

