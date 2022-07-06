import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
from scipy import stats
from uncertainties import ufloat
from scipy.optimize import curve_fit
import scipy.constants as const


lambda_, deg1_1, min1_1, deg1_2, min1_2,  deg2_1, min2_1, deg2_2, min2_2, deg3_1, min3_1, deg3_2, min3_2 = np.genfromtxt("data/daten.txt", unpack=True)

z,B=np.genfromtxt("data/bfeld.txt", unpack=True)

plt.figure()
plt.xlabel(r"z / $\si{\milli\meter}$")
plt.ylabel(r"B(z) / $\si{\milli\tesla}$")
plt.plot(z, B,".", label="B(z)")
plt.legend(loc="best")
plt.grid()
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/bfeld.pdf')


# umrechnung
deg1_1 = deg1_1 + min1_1/60 
deg1_2 = deg1_2 + min1_2/60
deg2_1 = deg2_1 + min2_1/60
deg2_2 = deg2_2 + min2_2/60
deg3_1 = deg3_1 + min3_1/60
deg3_2 = deg3_2 + min3_2/60

winkelprobe1 = 0.5*(deg1_1 - deg1_2)
winkelprobe2 = 0.5*(deg2_1 - deg2_2) 
winkelrein = 0.5*(deg3_1 - deg3_2) 

winkelprobe1 = np.deg2rad(winkelprobe1)/1.36e-3
winkelprobe2 = np.deg2rad(winkelprobe2)/1.296e-3
winkelrein = np.deg2rad(winkelrein)/5.11e-3
####

np.savetxt("data/winkelprobe1.txt", winkelprobe1)
np.savetxt("data/winkelprobe2.txt", winkelprobe2)
np.savetxt("data/winkelrein.txt", winkelrein)



plt.figure()
plt.xlabel(r"$\lambda^2$ / $(\si{\micro\meter})^2$")
plt.ylabel(r"$\theta_{KR}$ / $\frac{\mathrm{rad}}{\mathrm{m}}$")
plt.plot((lambda_**2)*1e12, winkelprobe1,".", mew=3, label="Probe 1")
plt.plot((lambda_**2)*1e12, winkelprobe2,".", mew=3, label="Probe 2")
plt.plot((lambda_**2)*1e12, winkelrein,".", mew=3, label="Reine Probe")
plt.legend(loc="best")
plt.grid()
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/winkel.pdf')

diff1 = np.abs(winkelprobe1-winkelrein) 
diff2 = np.abs(winkelprobe2-winkelrein)

np.savetxt("data/differenz1.txt", diff1)
np.savetxt("data/differenz2.txt", diff2)


N1 = 1.2e24
N2 = 2.8e24
n = 3.57
B = 421e-3


def mass(x, a):
    return a * x 

params1, cov1 = curve_fit(mass, lambda_**2, diff1, p0=[12e12] )

lambda_new = np.linspace(min(lambda_)-2e-6, max(lambda_)+0.5e-7, 100)
plt.figure()
plt.xlabel(r"$\lambda^2$ / $(\si{\micro\meter})^2$")
plt.ylabel(r"$\theta_\mathrm{frei,1}$ / $\frac{\mathrm{rad}}{\mathrm{m}}$")
plt.plot((lambda_**2)*1e12, diff1, ".", mew=3,  label="Daten")
plt.plot((lambda_new**2)*1e12, mass(lambda_new**2, *params1), label="Ausgleichsgerade")
plt.legend(loc="best")
plt.grid()
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig("build/Probe1.pdf")

err1 = np.diag(np.sqrt(cov1))
a1 = ufloat(params1[0], err1[0])

mass1 = unp.sqrt(const.elementary_charge**3 *N1 *B / (8* np.pi**2 *const.epsilon_0 * const.speed_of_light**3 * n * a1))

params2, cov2 = curve_fit(mass, lambda_**2, diff2, p0=[12e12] )

plt.figure()
plt.xlabel(r"$\lambda^2$ / $(\si{\micro\meter})^2$")
plt.ylabel(r"$\theta_\mathrm{frei,2}$ / $\frac{\mathrm{rad}}{\mathrm{m}}$")
plt.plot((lambda_**2)*1e12, diff2, ".", mew=3,  label="Daten")
plt.plot((lambda_new**2)*1e12, mass(lambda_new**2, *params2),  label="Ausgleichsgerade")
plt.legend(loc="best")
plt.grid()
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig("build/Probe2.pdf")

err2 = np.diag(np.sqrt(cov2))
a2 = ufloat(params2[0], err2[0])

mass2 = unp.sqrt(const.elementary_charge**3 *N2 *B / (8* np.pi**2 *const.epsilon_0 * const.speed_of_light**3 * n * a2))

f=  open("build/results.txt", "w")
f.write(f"Steigung 1: {a1}, Masse 1: {mass1/const.electron_mass}\n")
f.write(f"Steigung 2: {a2}, Masse 2: {mass2/const.electron_mass}")
f.close()