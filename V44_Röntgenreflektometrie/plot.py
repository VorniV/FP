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

detektor_params, detektor_pcov = curve_fit(gauss, detektor_phi, detektor_i)
detektor_err = np.sqrt(np.diag(detektor_pcov))

# Plot
detektor_phi_new = np.linspace(detektor_phi[0]-0.05, detektor_phi[-1]+0.05, 10000)
plt.figure()
plt.xlabel(r"$\theta$ / \si{\degree}")
plt.ylabel(r"Anzahl Events $\cdot 10^{6}$")
plt.plot(detektor_phi, detektor_i*1e-6, ".", label="Datenpunkte")
plt.plot(detektor_phi_new, gauss(detektor_phi_new, *detektor_params)*1e-6, label="Ausgleichskurve")
plt.legend(loc="best") 
plt.grid()
plt.tight_layout()
plt.savefig("build/detektorscan.pdf")

amp = ufloat(detektor_params[0], detektor_err[0])
mu = ufloat(detektor_params[1], detektor_err[1])
sigma = ufloat(detektor_params[2], detektor_err[2])
print(f"amp= {amp:.4f}")
print(f"mu= {mu:.4f}")
print(f"sigma= {sigma:.4f}")

amp1 = amp/(sigma*np.sqrt(2*np.pi))
halb = sigma*2*np.sqrt(2*np.log(2))
print(f"MaxI= {amp1:.4f}")
print(f"Halbwertsbreite= {halb:.4f}")
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

# Normierung und Silizium Theoriekurve

alpha_crit = 0.315
alpha_si = 0.223
r_lambda = 1.54e-10
k = 2*np.pi / r_lambda
n = 1 - 7.6e-6 + 1.54e-8j*141/(4*np.pi)

# Ideale Kurve
def theorie(alpha):
    return (np.abs((k * np.sin(alpha)- k*np.sqrt(n**2-np.cos(alpha)**2))/(k * np.sin(alpha)+ k*np.sqrt(n**2-np.cos(alpha)**2))))**2

# Normierung
rel_i = rel_i/rel_i[rel_phi == alpha_crit]

# Plot
plt.figure()
plt.xlabel(r"$\theta$ / \si{\degree}")
plt.yscale("log")
plt.ylabel(r"Normierte Anzahl Events")
plt.plot(rel_phi, rel_i, label="Normierte Werte")
plt.plot(rel_phi, theorie(np.deg2rad(rel_phi)), label="Theoriekurve für ideal glatte Siliziumoberfläche")
plt.legend(loc="best") 
plt.grid()
plt.tight_layout()
plt.savefig("build/silizium.pdf")

# z_scan

z, z_i = np.genfromtxt("data/zscan1.txt", unpack=True)

plt.figure()
plt.xlabel(r"$z$ / \si{\milli\meter}")
plt.ylabel(r"Anzahl Events $\cdot 10^{5}$")
#plt.axvline(x = -0.36, linestyle="--", color="r", label=r"Strahlbreite $d$")
#plt.axvline(x = -0.12, linestyle="--", color="r")
plt.plot(z, z_i*1e-5, ".", label="Messwerte")
plt.legend(loc="best") 
plt.grid()
plt.tight_layout()
plt.savefig("build/zscan.pdf")

#rocking
rocking_phi, rocking_i = np.genfromtxt("data/rockingscan0_genauer_absorber=100.txt", unpack=True)

# Plot
plt.figure()
plt.xlabel(r"$\theta$ / \si{\degree}")
#plt.yscale("log")
plt.ylabel(r"Anzahl Events")
plt.plot(rocking_phi, rocking_i, ".", label="Messwerte")
plt.legend(loc="best") 
plt.grid()
plt.tight_layout()
plt.savefig("build/rocking.pdf")


# G-Faktor
D = 0.0215
a_g = 0.40
d_0 = np.sin(np.deg2rad(a_g))*D
G_faktor = D/d_0

# Liste der jeweiligen G-Faktoren
faktor = []
for i in rel_phi:
    if i < a_g:
        faktor = np.append(faktor, G_faktor*np.sin(np.deg2rad(i)))
    else:
        faktor = np.append(faktor, 1)

# Anwendung G-Faktor und erneute Normalisierung
faktor[0]= 0.001
rel_i = rel_i/faktor
rel_i = rel_i/rel_i[rel_phi == alpha_crit]

# Bestimmung der Dicke 

#np.savetxt("data/korrigiert.txt",np.array([rel_phi, rel_i]).T, fmt="%.3e")
#num, min_phi, min_i = np.genfromtxt("data/minima.txt", unpack=True)
