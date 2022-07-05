import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
from scipy import stats
from uncertainties import ufloat
from scipy.optimize import curve_fit
import scipy.constants as const

<<<<<<< HEAD
||||||| constructed merge base
x = np.linspace(0, 10, 1000)
y = x ** np.sin(x)
=======
plt.rcParams.update({'font.size': 24})
>>>>>>> backup

<<<<<<< HEAD
||||||| constructed merge base
plt.subplot(1, 2, 1)
plt.plot(x, y, label='Kurve')
plt.xlabel(r'$\alpha \:/\: \si{\ohm}$')
plt.ylabel(r'$y \:/\: \si{\micro\joule}$')
plt.legend(loc='best')

plt.subplot(1, 2, 2)
plt.plot(x, y, label='Kurve')
plt.xlabel(r'$\alpha \:/\: \si{\ohm}$')
plt.ylabel(r'$y \:/\: \si{\micro\joule}$')
plt.legend(loc='best')

# in matplotlibrc leider (noch) nicht möglich
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot.pdf')
=======
lambda_, deg1_1, min1_1, deg1_2, min1_2,  deg2_1, min2_1, deg2_2, min2_2, deg3_1, min3_1, deg3_2, min3_2 = np.genfromtxt("data/daten.txt", unpack=True)


# umrechnung
deg1_1 = deg1_1 + min1_1/60 
deg1_2 = deg1_2 + min1_2/60
deg2_1 = deg2_1 + min2_1/60
deg2_2 = deg2_2 + min2_2/60
deg3_1 = deg3_1 + min3_1/60
deg3_2 = deg3_2 + min3_2/60

winkelprobe1 = 0.5*(deg1_1 - deg1_2)
winkelprobe2 = 0.5*(deg2_1 - deg2_2) 
winkelprobe3 = 0.5*(deg3_1 - deg3_2) 

winkelprobe1 = np.deg2rad(winkelprobe1)/5.11e-3
winkelprobe2 = np.deg2rad(winkelprobe2)/1.296e-3
winkelprobe3 = np.deg2rad(winkelprobe3)/1.36e-3
####

np.savetxt("data/winkelprobe1.txt", winkelprobe1)
np.savetxt("data/winkelprobe2.txt", winkelprobe2)
np.savetxt("data/winkelprobe3.txt", winkelprobe3)



plt.figure(figsize=(15,8))
plt.xlabel(r"$\lambda^2$ / $(\mathrm{\mu m})^2$")
plt.ylabel(r"$\frac{\theta}{\mathrm{d}}$ / $\frac{\mathrm{rad}}{\mathrm{m}}$")
plt.plot((lambda_**2)*1e12, winkelprobe1, "x", mew=3, markersize=14, label="Reine Probe")
plt.plot((lambda_**2)*1e12, winkelprobe2, "x", mew=3, markersize=14, label="Probe 1")
plt.plot((lambda_**2)*1e12, winkelprobe3, "x", mew=3, markersize=14, label="Probe 2")
plt.legend(loc="best")
plt.grid()
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot.pdf')

diff21 = np.abs(winkelprobe2-winkelprobe1) 
diff31 = np.abs(winkelprobe3-winkelprobe1)

np.savetxt("data/differenz1.txt", diff21)
np.savetxt("data/differenz2.txt", diff31)
>>>>>>> backup
