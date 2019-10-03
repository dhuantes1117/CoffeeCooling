import numpy as np
from numpy import pi as PI
from numpy.polynomial.polynomial import Polynomial as P
import matplotlib.pyplot as plt
from scipy.integrate import odeint
#------------------------------------------------------------------------METHOD DEFINITIONS-------------------------------------------------
"""
 newtCool - single term of newton's law of cooling in it's differential form. Returns the rate of change of temperature with respect to time for the given conditions specified by the arguments.

 Arguments:
 r  - the r constant for the object in question
 T  - the temperature of whatever object heat transfer is occuring
 To - the temperature of the object
"""
def newtCool(r, T_self, T_other):
    dTdt = r * (T_other - T_self)
    return dTdt
""" 
 dT_ldt - newtons cooling law with T_l, temperature of liquid (coffee, in our
 case). created in order to determine at provided values of temperature what
 the combined rate of loss of heat of the mug to the environment (first term)
 and to the coffee cup (second term) are.
 
 Arguments:
 should be sent to a tuple...
 r_env - r constant for area of liquid in contact with air (environment)
 T_l   - current temperature of liquid inside of cup
 T_env - temperature of surrounding environment (should not change)
 r_c   - r constant for area of liquid in contact with mug
 T_c   - current temperature of the cup
"""
def dT_ldt(r_env, T_l, T_env, r_c, T_c):
    #        heat lost to environment         heat lost from cup
    dT_ldt = newtCool(r_env, T_l, T_env) + newtCool(r_c, T_l, T_c)
    return dT_ldt
""" 
 dT_cdt - newtons cooling law with T_c, temperature of a cup containing some
 liquid with temperature T_l (coffee, in our case). created in order to
 determine at provided values of temperature what the combined rate of gain of
 heat of the liquid (assumef to be greater than the cup) to the cup (first
 term) and the rate of heat lost to the environment (second term) are.

 Arguments:
 should be sent to a tuple...
 r_ce  - r constant for area of cup in contact with air (environment)
 T_l   - current temperature of liquid inside of the cup
 T_env - temperature of surrounding environment (should not change)
 r_c   - r constant for area of liquid in contact with mug
 T_c - current temperature of the cup
"""
def dT_cdt(r_ce, T_l, T_env, r_c, T_c):
    #        heat gained from liquid      heat lost to sides
    dT_cdt = newtCool(r_c, T_c, T_l) + newtCool(r_ce, T_c, T_env)
    return dT_cdt

def cool(T_l, T_c, T_env, r_env, r_c, r_ce, t_list):
    liquidTemp  = T_l
    cupTemp     = T_c
    envTemp     = T_env
    liquidTempList, cupTempList, envTempList = [liquidTemp], [cupTemp], [envTemp]
    d_liquidTempList, d_cupTempList, d_envTempList = [], [], []
    for t in t_list:
        liquidTempChange    = dT_ldt(r_env, liquidTemp, envTemp, r_c, cupTemp)
        cupTempChange       = dT_cdt(r_ce, liquidTemp, envTemp, r_c, cupTemp)
        liquidTemp  += liquidTempChange
        cupTemp     += cupTempChange
        print("t = %g, liquidTemp = %f, cupTemp = %f, envTemp = %f\n cupTempChange = %g, liquidTempChange = %g" %(t, liquidTemp, cupTemp, envTemp, cupTempChange, liquidTempChange))
        liquidTempList.append(liquidTemp)
        cupTempList.append(cupTemp)
        envTempList.append(envTemp)
        d_liquidTempList.append(liquidTempChange)
        d_cupTempList.append(cupTempChange)
        d_envTempList.append(0)
    return (liquidTempList, cupTempList, envTempList)

"""
 multiple interface cooling
 ARGS: T: initial temperature (list)
 t: current time
 r_tuple: list of r values for each coffee-surface interface
"""
def MICooling(T, t, *Rgs):
    dTdt_list = np.array([coolingLaw(T, t, r) for r in Rgs])
    total = np.sum(np.ndarray.transpose(dTdt_list), 0) 

def getR(cupProperties, Volume):
    H, r_a, r_b = cupProperties
    return (((H * r_a**3) - (r_a * Volume) + (r_b * Volume)) / H)**(1 / 3.0)

def geth(cupProperties):
    H, r_a, r_b = cupProperties
    return H / ((r_b / r_a) - 1)

#Calculated by Patrick Cook
def SurfaceAreafromVolume(cupProperties, Volume):
    H, r_a, r_b = cupProperties
    R = getR(cupProperties, Volume)
    h = geth(cupProperties)
    return (PI**2 * R**4 + 9 * ((Volume + (PI / 3) * r_a**2 * h) / R)**2)**(1 / 2.0) - PI * r_a * (r_a**2 + h**2)**(1 / 2.0) + PI * r_a**2

def frustumVolume(h, r1, r2):
    return (1 / 3.0) * PI * h * (r1**2 + r2**2 + (r1 * r2))

def cupVolume(cupProperties, volum, thick):
    H, r_a, r_b = cupProperties
    return frustumVolume(geth(cupProperties) + thick, r_a + thick, r_b + thick) - volum

def r(h, A, p, v, c):
    return (h * A) / (p * v * c)
# ------------------------------------------------------VARIABLE DECLARATIONS---------------------------------------------------------------
                   # coffee cup parameters
max_l = 11.27      # cm
max_SA = 263       # cm^2
max_height = 11.34 # quick mafs
rA = 3             # cm
rA_ext = 3.3       # cm
rB = 3.86          # cm
rB_ext = 4.2       # cm
thickness = (rB_ext - rB + rA_ext - rA) / 2 # cm

# given parameters
c_water = 4181     # J kg-1 degC
c_ceramic = 900
h_water = 15       # Wdeg C^-1 m^-2
h_pork  = h_water * 0.606 / 1.5 
A = 5.E-2          # m^2
m = 0.35           # kg

T_env = 25         # degC
T_cup = T_env      # degC
T_coffee = 100     # degC

p_water = 997
p_pork  = 2403

# ------------------------------------------------------UNIT CONVERSIONS---------------------------------------------------------------
max_l /= 100      # m
max_height /= 100 # quick mafs
rA /= 100         # m
rB /= 100         # m
rB_ext /= 100     # m
rA_ext /= 100     # m
thickness /= 100  # m

cup = (max_height, rA, rB)
ext = (max_height + thickness, rA_ext, rB_ext)

#calculating physical parameters of our cup/coffee system
cup          = (max_height, rA, rB)
ext          = (max_height + thickness, rA_ext, rB_ext)
vol_coffee   = 8 # oz
vol_coffee  *= 2.95735e-5 #oz to m^3 conversion
vol_container= frustumVolume(vol_coffee / (PI * rA ** 2), rA_ext, getR(cup, vol_coffee) + thickness)
vol_cup      = vol_container - vol_coffee
m_coffee     = p_water * vol_coffee #kg conversion
sa_coffee    = SurfaceAreafromVolume(cup, vol_coffee) # m^2
sa_ext       = SurfaceAreafromVolume((max_height + thickness, rA_ext, rB_ext), vol_container)
#calculate r constants for different interfaces
r_env = r(h_water, PI * getR(cup, vol_coffee)**2, p_water, vol_coffee, c_water)
r_cup = r(h_water, sa_coffee, p_water, vol_coffee, c_water)
r_ce  = r(h_pork, sa_ext, p_pork, vol_cup, c_ceramic)

times = np.linspace(0, 3600, num=36e3 + 1, endpoint=True)

lTemps, cTemps, eTemps = cool(T_coffee, T_env, T_env , r_env, r_cup, r_ce, times)

times /= 60

plt.plot(times, lTemps[0:-1])
plt.plot(times, cTemps[0:-1])
plt.plot(times, eTemps[0:-1])
plt.legend(("liquid temperature", "cup temperature", "environment temperature"))
plt.savefig("coolingfig.png")
plt.xlabel("time (min)")
plt.ylabel("temperature ($^o$C)")
plt.title("Heat loss of coffee in a mug over time")
plt.show()
plt.close()
