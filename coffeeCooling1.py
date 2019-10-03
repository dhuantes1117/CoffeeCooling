import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
c = 4181    # J kg-1 degC
h = 15      # Wdeg C^-1 m^2
A = 5.E-2   # m^2
m = 0.35    # kg
T_env = 22  # degC
T0    = 90  # degC

r = (h * A) / (m * c) # degC s^-1

def coolingLaw(T, t):
    dTdt = -r * (T[0] - T_env)
    return np.array([dTdt])

times = np.linspace(0, 3600, num=601)
solution = odeint(coolingLaw, np.array([T0]), times)

T = solution[:,0]

plt.plot(times / 60, T)
plt.title("Coffee temperature over time")
plt.xlabel("Time $(min)$")
plt.ylabel("Temperature $(^oC)$")
plt.savefig("Challenge1Figure.png")
plt.show()
plt.close()
