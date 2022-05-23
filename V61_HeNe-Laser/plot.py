
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from scipy.optimize import curve_fit
from uncertainties import ufloat
import uncertainties.unumpy as unp

# TEM00 -------------------------------------------------------------------------------------------

r, I = np.genfromtxt('TEM00.txt', unpack=True)
I = I * 1e3

def exp1(r, I0, r0, w):
    return I0 * np.exp(-(r - r0)**2/(2*w**2)) + 0.2

par, cov = optimize.curve_fit(exp1, r, I)#, p0=[0.9, 90])
        #I0 = ufloat(par[0], np.sqrt(cov[0][0]))
        #r0 = ufloat(par[1], np.sqrt(cov[1][1]))
        #w = ufloat(par[2], np.sqrt(cov[2][2]))
        #print(I0, r0, w)
I0 = par[0]
r0 = par[1]
w  = par[2]

print(I0, r0, w)

rx = np.linspace(np.min(r), np.max(r), 1000)

plt.plot(r, I, 'x', color='r', label='Daten TEM$_{00}$-Mode')
plt.plot(rx, exp1(rx, I0, r0, w), '-', color='b', label='Ausgleichsrechnung')
plt.xlim(np.min(r)-3, np.max(r)+3)
plt.xlabel(r'Abstand $r$ (mm)')
plt.ylabel(r'Intensität $I(r)$ ($\mu$W)')
#plt.ylim(-5e4, Imax+5e4)
plt.legend(loc='best')

plt.tight_layout()
plt.savefig('TEM00.pdf')
plt.close()


# TEM10 -------------------------------------------------------------------------------------------

r, I = np.genfromtxt('TEM01.txt', unpack=True)
#print(np.min(I))
#I = I - np.min(I)

def exp2(r, I0, r0, r1, w):
    return I0 * 4* (r-r1)**2 * np.exp(-2*(r - r0)**2/(2*w**2)) + 0.2

par, cov = optimize.curve_fit(exp2, r, I)#, p0=[0.151, 0.5, 8])
        #I0 = ufloat(par[0], np.sqrt(cov[0][0]))
        #r0 = ufloat(par[1], np.sqrt(cov[1][1]))
        #w = ufloat(par[2], np.sqrt(cov[2][2]))
I0 = par[0]
r0 = par[1]
r1 = par[2]
w  = par[3]

print(I0, r0, r1, w)

rx = np.linspace(np.min(r)-3, np.max(r)+3, 1000)

plt.plot(r, I, 'x', color='r', label='Daten TEM$_{01}$-Mode')
plt.plot(rx, exp2(rx, I0, r0, w), '-', color='b', label='Ausgleichsrechnung')
plt.xlim(np.min(r)-3, np.max(r)+3)
plt.xlabel(r'Abstand $r$ (mm)')
plt.ylabel(r'Intensität $I(r)$ ($\mu$W)')
#plt.ylim(-5e4, Imax+5e4)
plt.legend(loc='best')

plt.tight_layout()
plt.savefig('TEM01.pdf')
plt.close()


# Polarisation ----------------------------------------------------------------------------------

phi, I = np.genfromtxt('polarisation.txt', unpack=True)

def cos(phi, I0, phi0):
    return I0 * (np.cos((phi + phi0)*(2*np.pi/360)))**2 

par, cov = optimize.curve_fit(cos, phi, I, p0=[0.9, 90])
        #I0 = ufloat(par[0], np.sqrt(cov[0][0]))
        #phi0 = ufloat(par[1], np.sqrt(cov[1][1]))
I0   = par[0]
phi0 = par[1]
print(I0, phi0)

phix = np.linspace(np.min(phi), np.max(phi), 1000)

plt.plot(phi, I, 'x', color='r', label='Daten')
plt.plot(phix, cos(phix, I0, phi0), '-', color='b', label='Ausgleichsrechnung')
plt.xlim(np.min(phi)-10, np.max(phi)+10)
plt.xlabel(r'Polarisationswinkel $\phi$ (°)')
plt.ylabel(r'Intensität $I(\phi)$ (mW)')
#plt.ylim(-5e4, Imax+5e4)
plt.legend(loc='best')

plt.tight_layout()
plt.savefig('pol.pdf')
plt.close()


# Modendifferenz ------------------------------------------------------------------------------------

Per Hand

# Doppler-Verschiebung ------------------------------------------------------------------------------

Per Hand

# Bestimmung der Wellenlänge -------------------------------------------------------------------------

n1, d1 = np.genfromtxt('gitter1.txt', unpack=True)
n2, d2 = np.genfromtxt('gitter2.txt', unpack=True)
n3, d3 = np.genfromtxt('gitter2.txt', unpack=True)
n4, d4 = np.genfromtxt('gitter2.txt', unpack=True)

g1 = 80000   #Gitterkonstante (Passende Einheit!)
g2 = 100000
g3 = 600000
g4 = 1200000
L1 = 0.5   #Abstand Gitter-Schirm (Passende Einheit!)
L2 = 0.3

def lam(n, d, L, g):
    return unp.sin( unp.tan(d/L) ) / (g * np.sqrt(n**2) )

lam1 = np.mean( lam(n1, d1, L1, g1) )
lam2 = np.mean( lam(n2, d2, L1, g2) )
lam3 = np.mean( lam(n3, d3, L2, g3) )
lam4 = np.mean( lam(n4, d4, L2, g4) )
f1 = (np.std (lam(n1, d1, L1, g1)) )**2
f2 = (np.std (lam(n2, d2, L1, g2)) )**2
f3 = (np.std (lam(n3, d3, L2, g3)) )**2
f4 = (np.std (lam(n4, d4, L2, g4)) )**2

lamges = (lam1 + lam2 +lam3 +lam4) / 4
fges = f1 + f2 + f3 +f4

lamges2 = np.mean( lam(n1, d1, L1, g1), lam(n2, d2, L1, g2), lam(n3, d3, L2, g3), lam(n4, d4, L2, g4) )
fges2 =(np.std( lam(n1, d1, L1, g1), lam(n2, d2, L1, g2), lam(n3, d3, L2, g3), lam(n4, d4, L2, g4)) )**2

print(lam1, f1, lam2, f2, lam3, f3, lam4, f4, lamges, fges, lamges2, fges2)
#np.savetxt('auswertung.dat', np.column_stack([lam1, lam2, lamges]), header='l1, l2, lges')

