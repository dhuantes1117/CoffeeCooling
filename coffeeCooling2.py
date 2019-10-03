import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
c = 4181    # J kg-1 degC
h = 15      # Wdeg C^-1 m^2
A = 5.E-2   # m^2
m = 0.35    # kg
T_env = 73.4  # degF
T0    = 174   # degF

A = 0.00448347
m = 0.353817

T_env = (T_env - 32) * (5 / 9.0)
T0 = (T0 - 32) * (5 / 9.0)

r = (h * A) / (m * c) # degC s^-1

def coolingLaw(T, t):
    dTdt = -r * (T[0] - T_env)
    return np.array([dTdt])

times = np.linspace(0, 3600, num=601)
solution = odeint(coolingLaw, np.array([T0]), times)

T = solution[:,0]

plt.plot(times / 60, (T * (9 / 5.0)) + 32)
plt.title("Coffee temperature over time")
plt.xlabel("Time $(min)$")
plt.ylabel("Temperature $(^oF)$")
plt.ylim(70, 180)
plt.savefig("Challenge2Figure.png")
plt.show()
plt.close()
