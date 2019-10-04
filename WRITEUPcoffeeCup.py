import numpy as np
from numpy import pi as PI
from numpy import e 
from numpy.polynomial.polynomial import Polynomial as P
import matplotlib.pyplot as plt
from scipy.integrate import odeint
#------------------------------------------------------------------------METHOD DEFINITIONS-------------------------------------------------
def newtCool(r, T_self, T_other):
    dTdt = r * (T_other - T_self)
    return dTdt

def dT_ldt(r_env, T_l, T_env, r_c, T_c):
    #        heat lost to environment         heat lost from cup
    dT_ldt = newtCool(r_env, T_l, T_env) + newtCool(r_c, T_l, T_c)
    return dT_ldt

def dT_cdt(r_ce, T_l, T_env, r_c, T_c):
    #        heat gained from liquid      heat lost to sides
    dT_cdt = newtCool(r_c, T_c, T_l) + newtCool(r_ce, T_c, T_env)
    return dT_cdt

def cool(T_l, T_c, T_env, r_env, r_c, r_ce, t_list):
    liquidTemp  = T_l
    cupTemp     = T_c
    envTemp     = T_env
    liquidTempList, cupTempList, envTempList = [liquidTemp], [cupTemp], [envTemp]
    DeltaT = t_list[1] - t_list[0]
    for t in t_list:
        liquidTempChange    = dT_ldt(r_env, liquidTemp, envTemp, r_c, cupTemp)
        cupTempChange       = dT_cdt(r_ce, liquidTemp, envTemp, r_c, cupTemp)
        liquidTemp  += liquidTempChange * DeltaT
        cupTemp     += cupTempChange * DeltaT
        liquidTempList.append(liquidTemp)
        cupTempList.append(cupTemp)
        envTempList.append(envTemp)
    return (liquidTempList, cupTempList, envTempList)

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

def r(h, A, p, v, c):
    return (h * A) / (p * v * c)

def expFit(x, y):
    logFit = np.polyfit(np.array(x), np.log(np.array(y)), 1, w=np.sqrt(y))
    return np.array([e**logFit[1], logFit[0]])
# ------------------------------------------------------VARIABLE DECLARATIONS---------------------------------------------------------------
                   # coffee cup parameters
max_l = 11.27      # cm
max_height = 11.34 # cm 
rA = 3             # cm
rA_ext = 3.3       # cm
rB = 3.86          # cm
rB_ext = 4.2       # cm
thickness = (rB_ext - rB + rA_ext - rA) / 2 # cm

c_water = 4181     # J kg-1 degC
c_ceramic = 1085   # J kg-1 degC
h_water = 15       # Wdeg C^-1 m^-2
k_pork  = 1.5      # Wdeg C^-1 m^-2
p_water = 997      # kg / m^3
p_pork  = 2403     # kg / m^3

T_env = 23         # degC
T_cup = T_env      # degC
T_coffee = 74      # degC
# ------------------------------------------------------UNIT CONVERSIONS---------------------------------------------------------------
max_l /= 100      # m
max_height /= 100 # m
rA /= 100         # m
rB /= 100         # m
rB_ext /= 100     # m
rA_ext /= 100     # m
thickness /= 100  # m

#calculating physical parameters of our cup/coffee system
cup          = (max_height, rA, rB)
ext          = (max_height + thickness, rA_ext, rB_ext)
vol_coffee   = 12 # oz
vol_coffee  *= 2.95735e-5 #oz to m^3 conversion
vol_container= frustumVolume(vol_coffee / (PI * rA ** 2), rA_ext, getR(cup, vol_coffee) + thickness)
vol_cup      = vol_container - vol_coffee
m_coffee     = p_water * vol_coffee #kg conversion
sa_coffee    = SurfaceAreafromVolume(cup, vol_coffee) # m^2
sa_ext       = SurfaceAreafromVolume(ext, vol_container) # m^2



# for the h value for r_env I used the given value of 15 W m^-2 K^-1
#r_env = r(h_water, PI * getR(cup, vol_coffee)**2, p_water, vol_coffee, c_water)
r_env = r(5, PI * getR(cup, vol_coffee)**2, p_water, vol_coffee, c_water)
# for the h value of r_cup I used it's thermal conductivity / the thickness of the wall, using the heat transfer of a pipe as a (strong) approximation
#r_cup = r(k_pork / thickness, sa_coffee, p_water, vol_coffee, c_water)
r_cup = r(12.87, sa_coffee, p_water, vol_coffee, c_water)
# for the h value of r_ce I multiplied it's linear thermal conductivity by the surface area and divided by its volume as a rough estimate of it's thermal conductivity per unit area
#r_ce  = r(k_pork * sa_ext / vol_container, sa_ext, p_pork, vol_cup, c_ceramic)
r_ce  = r(1.528, sa_ext, p_pork, vol_cup, c_ceramic)

print("h_water: %g" %(h_water))
print("h_pork: %g" %(k_pork / thickness))
print("h_pork to Env:  %g"  %(k_pork * sa_ext / vol_container))

print("r_env: %g" %(r_env))
print("r_cup: %g" %(r_cup))
print("r_ce:  %g"  %(r_ce))

print("r_env: %g" %(r_env))
print("r_cup: %g" %(r_cup))
print("r_ce:  %g"  %(r_ce))

print("Coffee Area: %g" %(PI * getR(cup, vol_coffee)**2))
print("Coffee mass: %g" %(vol_coffee * p_water))

times = np.arange(0, 3600, 0.1)

lTemps, cTemps, eTemps = cool(T_coffee, T_env, T_env , r_env, r_cup, r_ce, times)
times = times / 60


experimentalX   = np.array([1, 2, 3, 4, 5, 6, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 600])
experimentalY   = np.array([165, 163, 161, 159, 156, 154, 146, 139, 132, 126 , 121, 116, 112, 108, 106 , 102 , 100])
minErr = 9e9
minErrnumi = -1
minErrnumj = -1
minErrnumk = -1
for i in np.arange(1, 1, 0.1):#water range
    for j in np.arange(-1, 1, 0.1):#air range
        for k in np.arange(-1, 1, 0.1):#pork to air range
            h_waterenv = ((i)**-1 + (j)**-1)**-1
            h_watercup = ((i)**-1 + (k_pork / thickness)**-1)**-1
            h_cupenvir = ((k_pork / thickness)**-1 + (k)**-1)**-1
            r_env = r(i, PI * getR(cup, vol_coffee)**2, p_water, vol_coffee, c_water)
            r_cup = r(j, sa_coffee, p_water, vol_coffee, c_water)
            lTemps, cTemps, eTemps = cool(T_coffee, T_env, T_env , r_env, r_cup, r_ce, times)
            lTemps = np.array(lTemps[0:-1])
            cTemps = np.array(cTemps[0:-1])
            eTemps = np.array(eTemps[0:-1])
            coeff           = expFit(times, lTemps)
            expectedY       = coeff[0] * e**(coeff[1] * experimentalX)
            chisq           = np.sum((experimentalY - expectedY)**2)
            if minErr > chisq:
                minErr = chisq
                minErrnumi = i
                minErrnumj = j
                minErrnumk = k
print("A minimum error of %g occured at i = %g, j = %g, k = %g" %(minErr, minErrnumi, minErrnumj, minErrnumk))
print("bye bye")
exit()
lTemps = np.array(lTemps[0:-1])
cTemps = np.array(cTemps[0:-1])
eTemps = np.array(eTemps[0:-1])

coeff           = expFit(times, lTemps)
experimentalX   = np.array([1, 2, 3, 4, 5, 6, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 600])
experimentalY   = np.array([165, 163, 161, 159, 156, 154, 146, 139, 132, 126 , 121, 116, 112, 108, 106 , 102 , 100])
expectedY       = coeff[0] * e**(coeff[1] * experimentalX)
chisq           = np.sum((experimentalY - expectedY)**2 / expectedY)
print("Chi Squared: %g" %(chisq))
plt.plot(times, lTemps * (9 / 5.0) + 32)
plt.plot(times, cTemps * (9 / 5.0) + 32)
plt.plot(times, eTemps * (9 / 5.0) + 32)
plt.legend(("Liquid temperature", "Cup temperature", "Environment Temperature"))
plt.savefig("coolingfig2.png")
plt.xlabel("Time (min)")
plt.ylabel("Temperature ($^o$F)")
plt.ylim(70, 180.01)
plt.title("Heat loss of coffee in a mug over time")
plt.show()
plt.close()
