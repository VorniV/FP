
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
from scipy.optimize import curve_fit
from uncertainties import ufloat
import uncertainties.unumpy as unp

#--------- Magnetfeld Eichung ------------------------------------------------------

I = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8])
B0 = np.array([7.6  , 54.2 , 102.5, 151.6, 203.1, 255.1, 297.4, 338.3, 381.9, 418.7, 454.3, 482.1, 509.1, 530.0 , 547.8, 563.3, 577.5])

B = B0 * 10**(-3)

def Mag(x, a, b, c):
    return a*x**2 + b*x + c

par, cov = curve_fit(Mag, I, B)

errpar = np.sqrt(np.diag(cov))

print(par[0],par[1], par[2] ,errpar[0], errpar[1], errpar[2])

x = np.linspace(0, 9)

plt.plot(I, B, 'x', color='r', label='Messwerte B')
plt.plot(x, Mag(x, *par), '-', color='b', label='Ausgleichskurve')
plt.xlim(0,9)
plt.xlabel(r'Stromstärke I [A]')
plt.ylabel(r'Magnetfeldstärke B [T]')
plt.ylim(0, 0.6)
plt.legend(loc='best')

plt.tight_layout()
plt.savefig('Magnetfeld.pdf')
plt.close()

           #------- Magnetfeld Berechnung -----------------

BRot = Mag(8,par[0],par[1], par[2])                                             #Magnetfeld Rot
FBRot = np.sqrt( (8**2 * errpar[0])**2  +  (8 *errpar[1])**2  +  (errpar[2])**2 )        #Fehler Magnetfeld Rot

BBlau1 = Mag(3.75,par[0],par[1], par[2])                                           #Magnetfeld Blau Sigma
FBBlau1 = np.sqrt( (3.75**2 * errpar[0])**2  +  (3.75 *errpar[1])**2  +  (errpar[2])**2 )      #Fehler

BBlau2 = Mag(8,par[0],par[1], par[2])                                           #Magnetfeld Blau Lambda
FBBlau2 = np.sqrt( (8**2 * errpar[0])**2  +  (8 *errpar[1])**2  +  (errpar[2])**2 )      #Fehler

print('Magnetfeld Rot: ', BRot, FBRot)
print('Magnetfeld Blau Sigma: ', BBlau1, FBBlau1)
print('Magnetfeld Blau Lambda: ', BBlau2, FBBlau2)

#-------- Spektrallinien -----------------------------------------------

h = 6.626 * 10**(-34)
c = 299792458
mu = 9.274 * 10**(-24)
Lambda1 = 643.8 * 10**(-9)          #Rote Wellenlänge
Lambda2 = 480 * 10**(-9)            #Blaue Wellenlänge


def Lam(e,E,F):                                                     #Wellenlängenverschiebung
    return (1/2) * (e / E) * F 

        #------------------ Rot ----------------------------------------

s = np.array([92.0 , 98 , 103, 102, 106, 109, 113, 117, 119, 122])   #klein delta s
S = np.array([198.0, 204, 207, 214, 218, 225, 234, 235, 255, 261])   #Groß Delta S

L = 4.894 *10**(-11)                                           #Dispersionsgebiet Rot (Lambda_D)


LambdaRot = np.mean(Lam(s,S,L))                                    #Wellenlängenverschiebung Rot
FehlerRot = np.std(Lam(s,S,L)) / np.sqrt( len(s))
                        

gRot = LambdaRot * ( (h * c) / (mu * BRot * Lambda1**2))            #Lande-Faktor Rot
FgRot= ((h* c) / (mu * Lambda1**2))  *  np.sqrt( ( (1/BRot) * FehlerRot)**2 + ((LambdaRot / BRot**2) * FBRot)**2 )

print('Lambda Rot: ',LambdaRot, FehlerRot)  
print('Lande-Faktor Rot: ', gRot, FgRot)

        #----------------- Blau (Sigma) ------------------------------

s1 = np.array([72, 74, 82, 80, 84, 78, 80, 82, 78, 86])             #klein delta s
S1 = np.array([152, 156, 158, 156, 162, 158, 154, 160, 168, 156])   #Groß Delta S

L1 = 2.695 * 10**(-11)                                              #Dispersionsgebiet Blau Sigma (Lambda_D)


LambdaBlau1 = np.mean(Lam(s1,S1,L1))                                    #Wellenlängenverschiebung Blau Sigma
FehlerBlau1 = np.std(Lam(s1,S1,L1)) / np.sqrt( len(s1))
                        

gBlau1 = LambdaBlau1 * ( (h * c) / (mu * BBlau1 * Lambda2**2))            #Lande-Faktor Blau Sigma
FgBlau1= ((h* c) / (mu * Lambda2**2))  *  np.sqrt( ( (1/BBlau1) * FehlerBlau1)**2 + ((LambdaBlau1 / BBlau1**2) * FBBlau1)**2 )

print('Lambda Blau Sigma: ',LambdaBlau1, FehlerBlau1)  
print('Lande-Faktor Blau Sigma: ', gBlau1, FgBlau1)

        #-------------- Blau (Lambda) ------------------------------
#
#s2 = np.array([])   #klein delta s
#S2 = np.array([])   #Groß Delta S
#
#L1 = 2.695 * 10**(-11)                                              #Dispersionsgebiet Blau Lambda (Lambda_D)
#
#
#LambdaBlau2 = np.mean(Lam(s2,S2,L2))                                    #Wellenlängenverschiebung Blau Lambda
#FehlerBlau2 = np.std(Lam(s2,S2,L2)) / np.sqrt( len(s2))
#                        
#
#gBlau2 = LambdaBlau2 * ( (h * c) / (mu * BBlau2 * Lambda2**2))            #Lande-Faktor Blau Lambda
#FgBlau2= ((h* c) / (mu * Lambda2**2))  *  np.sqrt( ( (1/BBlau2) * FehlerBlau2)**2 + ((LambdaBlau2 / BBlau2**2) * FBBlau2)**2 )
#
#print('Lambda Blau Lambda: ',LambdaBlau2, FehlerBlau2)  
#print('Lande-Faktor Blau Lambda: ', gBlau2, FgBlau2)








