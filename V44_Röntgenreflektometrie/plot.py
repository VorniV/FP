import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
from scipy import stats
from uncertainties import ufloat
from scipy.optimize import curve_fit


detektor_phi, detektor_i = np.genfromtxt("data/Detektorscan.txt", unpack=True)

# fit
def gauss(x, amp, mu, sigma):
    return amp/(sigma*np.sqrt(2*np.pi)) * np.exp(-(x- mu)**2/(2*sigma**2))

# fit and errors
detektor_params, detektor_pcov = curve_fit(gauss, detektor_phi, detektor_i)
detektor_err = np.sqrt(np.diag(detektor_pcov))

# Plot
detektor_phi_new = np.linspace(detektor_phi[0]-0.05, detektor_phi[-1]+0.05, 10000)
plt.figure()
plt.xlabel(r"$\theta$ / \si{\degree}")
plt.ylabel(r"Anzahl Ev $\cdot 10^{6}$")
plt.plot(detektor_phi, detektor_i*1e-6, ".", label="Datenpunkte")
plt.plot(detektor_phi_new, gauss(detektor_phi_new, *detektor_params)*1e-6, label="Ausgleichskurve")
plt.legend(loc="best") 
plt.grid()
plt.tight_layout()
plt.savefig("build/detektorscan.pdf")
##############
messung_phi, messung_i = np.genfromtxt("data/refl1.txt", unpack=True)

diffus_phi, diffus_i = np.genfromtxt("data/refl1untergrund.txt", unpack=True)

#relative data
rel_i = messung_i - diffus_i
rel_phi = diffus_phi

#Plot
plt.figure()
plt.xlabel(r"$\theta$ / \si{\degree}")
plt.yscale("log")
plt.ylabel(r"Anzahl Events")
plt.plot(messung_phi, messung_i, label="Messwerte")
plt.plot(diffus_phi, diffus_i, label="Diffuser Scan")
plt.plot(rel_phi, rel_i, label="Korrigierte Messwerte")
plt.legend(loc="best") 
plt.grid()
plt.tight_layout()
plt.savefig("build/messwerte_relativ.pdf")

