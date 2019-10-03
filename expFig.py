import numpy as np
from numpy import e
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import sys

def readvalues(filename):
    values = {}
    cat = []
    filestring  = open(filename).read().strip()
    lines = [line.split(", ") for line in filestring.split("\n") if not line[0] == "#"]
    cat = lines[0]
    values = [[float(d) for d in l] for l in lines[1:]]
    x = [l[0] for l in values]
    y = [l[1] for l in values]
    return (cat, x, y)

cat, x, y = readvalues("realExperiment.csv")
logFit = np.polyfit(np.array(x), np.log(np.array(y)), 1, w=np.sqrt(y))
xFit = np.linspace(0, 60, 3600)
txt = "y = %f e$^{%fx}$" %(e**logFit[1], logFit[0])
yFit = e**logFit[1] * e**(xFit * logFit[0])
plt.figtext(0.5, 0.8, txt)
plt.plot(xFit, yFit, "g-")
plt.plot(x, y, "bo")
plt.xlabel(cat[0])
plt.ylabel(cat[1])
plt.ylim(70, 180)
plt.title("Coffee cooling experimental results")
plt.savefig("experimentalResults.png")
plt.show()
plt.close()
