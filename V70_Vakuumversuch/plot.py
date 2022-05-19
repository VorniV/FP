import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import uncertainties as unc
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy import stats
import pandas as pd
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
from tabulate import tabulate

#--------------Daten----------------------
#Turbopumpe:
VT = ufloat(33.0, 3.3)

# Anfangsdruck p0 = 5 * 10^-3 mbar
p0T = ufloat(5 * 10**(-3), 0.3*5 * 10**(-3))

# Enddruck pe = 1.5 * 10^-5 mbar
peT = ufloat(1.5 * 10**(-5), 0.3*1.5 * 10**(-5))

# Daten für die p(t)-Kurve der Turbopumpe
t1, p_1T, p_2T, p_3T = np.genfromtxt('data/evakuierung_turbo.txt',comments='#',delimiter=',',unpack=True)

# Daten für die Leckrate der Turbopumpe
t, p_1TL1, p_2TL1, p_3TL1 = np.genfromtxt('data/leck_turbo_pg=1e-4.txt',comments='#',delimiter=',',unpack=True)
t, p_1TL2, p_2TL2, p_3TL2 = np.genfromtxt('data/leck_turbo_pg=2e-4.txt',comments='#',delimiter=',',unpack=True)
t, p_1TL5, p_2TL5, p_3TL5 = np.genfromtxt('data/leck_turbo_pg=5e-5.txt',comments='#',delimiter=',',unpack=True)
t, p_1TL7, p_2TL7, p_3TL7 = np.genfromtxt('data/leck_turbo_pg=7e-5.txt',comments='#',delimiter=',',unpack=True)


#Drehschieberpumpe:
VD = ufloat(10.0, 0.8)

# Anfangsdruck p0 = 998mbar
p0D = ufloat(998, 0.003*1200)

# Enddruck pe = 0.012 mbar
peD = ufloat(0.012, 0.0012)

# Daten für die p(t)-Kurve der Drehschieberpumpe
t2, p_1D, p_2D, p_3D = np.genfromtxt('data/evakuierung_dreh.txt',comments='#',delimiter=',',unpack=True)

# Daten für die Leckrate der Drehschieberpumpe
t3, p_1DL1, p_2DL1, p_3DL1 = np.genfromtxt('data/leck_dreh_pg=0.5.txt',comments='#',delimiter=',',unpack=True)
t3, p_1DL9, p_2DL9, p_3DL9 = np.genfromtxt('data/leck_dreh_pg=9.7.txt',comments='#',delimiter=',',unpack=True)
t3, p_1DL50, p_2DL50, p_3DL50 = np.genfromtxt('data/leck_dreh_pg=50.txt',comments='#',delimiter=',',unpack=True)
t3, p_1DL100, p_2DL100, p_3DL100 = np.genfromtxt('data/leck_dreh_pg=100.txt',comments='#',delimiter=',',unpack=True)

#----------------Berechnungen--------------
def mittel(a,b,c):
    if (np.size(a) != np.size(b) or np.size(b) != np.size(c)):
        print("NICHT ALLE WERTE EINGETRAGEN")
        if(np.size(a)>np.size(b)):
            print("FEHLER b)")
        else:
            if(np.size(b)>np.size(c)):
                print("FEHLER c)")
            else:
                print("FEHLERa)")
    arr = unp.uarray(np.zeros(np.size(a)),np.zeros(np.size(a)))
    for i in range(0,np.size(a)):
        arr[i] = unc.ufloat(np.mean([a[i], b[i], c[i]]), np.std([a[i], b[i], c[i]])/ np.sqrt(3))
    return arr

def logprinter(t,a,b,c,name):
   
    d=mittel(a,b,c)
    
    if(name == "dreh"):
        d1=fehler_dreh(d)
        table1 ={'Zeit': t, 'log': logungdreh(d1,0,np.size(d1))}
        headers = ["t [s]","$\ln(\frac{p-p_{E}}{p_{0}-p_{E}})$"]
        print("\n", tabulate(table1, headers,tablefmt = "latex_raw"))  
    else:
        d1=fehler_turbo(d)
        table1 ={'Zeit': t, 'log': logung(d1,0,np.size(d1))}
        headers = ["t [s]","$\ln(\frac{p-p_{E}}{p_{0}-p_{E}})$"]
        print("\n", tabulate(table1, headers,tablefmt = "latex_raw"))

def turboprinter(t,a,b,c):
    a1=err_turbo(a)
    b1=err_turbo(b)
    c1=err_turbo(c)
    d=mittel(a,b,c)
    d1=fehler_turbo(d)
    table1 ={'Zeit': t,'Messreihe 1': a1, 'Messreihe 2': b1,  'Messreihe 3': c1, 'gemittelte Messwerte': d1}
    headers = ["t [s]", "p1 [mbar]", "p2 [mbar]", "p3 [mbar]", "\bar{p} [mbar]"]
    print("\n", tabulate(table1, headers,tablefmt = "latex_raw"))    
      

def drehprinter(t,a,b,c):
    a1=err_dreh(a)
    b1=err_dreh(b)
    c1=err_dreh(c)
    d=mittel(a,b,c)
    d1=fehler_dreh(d)
    table1 ={'Zeit': t,'Messreihe 1': a1, 'Messreihe 2': b1,  'Messreihe 3': c1, 'gemittelte Messwerte': d1}
    headers = ["t [s]", "p1 [mbar]", "p2 [mbar]", "p3 [mbar]", " $\bar{p}$ [mbar]"]
    print("\n", tabulate(table1, headers,tablefmt = "latex_raw"))    
    

def err_turbo(a):
    a = noms(a)
    err = np.zeros(np.size(a))
    for i in range(0,np.size(a)):
        if(a[i] <= 100):
            err[i] = 0.3 * a[i]
        else:
            err[i] = 0.5 * a[i]
    b =unp.uarray(a,err) 
    #b =unp.uarray(a,err/3)
    return b

def fehler_turbo(a):
    a = noms(a)
    err = stds(a)
    for i in range(0,np.size(a)):
        if(a[i] <= 100):
            err[i] = err[i]+0.3 * a[i]
        else:
            err[i] = err[i]+0.5 * a[i]
    b =unp.uarray(a,err) 
    #b =unp.uarray(a,err/3)
    return b

def err_dreh(a):                    
    a = noms(a)
    err = np.zeros(np.size(a))
    for i in range(0,np.size(a)):
        if(a[i] >= 10):
            err[i] = 0.003*1200
        else:
            if(a[i] >= 2 *10**(-3)):
                err[i] = a[i] * 0.1
            else:
                err[i] = a[i] * 2
    b =unp.uarray(a,err)    
    #b =unp.uarray(a,err/3)     
    return b

def fehler_dreh(a):                    
    a = noms(a)
    err = stds(a)
    for i in range(0,np.size(a)):
        if(a[i] >= 10):
            err[i] = err[i]+0.003*1200
        else:
            if(a[i] >= 2 *10**(-3)):
                err[i] = err[i]+a[i] * 0.1
            else:
                err[i] = err[i]+a[i] * 2
    b =unp.uarray(a,err)    
    #b =unp.uarray(a,err/3)     
    return b

def err_mess(a):
    if (a[0] >= 900):
        return err_dreh(a)
    else:
        return err_turbo(a)

def logung(a, grenz_1, grenz_2):
    return unp.log((a[grenz_1:grenz_2]-peT)/(p0T-peT))

def logungdreh(a, grenz_1, grenz_2):
    return unp.log((a[grenz_1:grenz_2]-peD)/(p0D-peD))

def lin(t,a,b):
    return a * t + b

size_label = 17

def saugung(t,a,b,c,name):
 
    t1=t
    mean = mittel(a,b,c)
    x=fehler_turbo(mean)
    grenz_1 = 5
    grenz_2 = 10

    V = unc.ufloat(33, 3.3)
    #theo = 77
    plt.clf()
    # für den ersten Bereich

    params_1, cov_1 = curve_fit(lin, t1[:grenz_1], noms(logung(x, 0, grenz_1)))

    cov_1 = np.sqrt(np.diag(cov_1))
    print("1. Fit:\n",f"m = {params_1[0]:.4f} \pm {cov_1[0]:.5f} \t n = {params_1[1]:.4f} \pm {cov_1[0]:.5f}")
    uparams_1 = unp.uarray(params_1, cov_1)
    S_1 = -1 * V * uparams_1[0]
    print(f"Saugvermögen: {noms(S_1):.4f} \pm {stds(S_1):.5f}")
    #rel_abw(theo,S_1)
    plt.errorbar(t1[:grenz_1], noms(logung(x, 0, grenz_1)), yerr= stds(logung(x, 0, grenz_1)), fmt='r.',label= "Messdaten Bereich 1")
    plt.plot(t1[:grenz_1], lin(t1[:grenz_1],*params_1), color='red',label="Fit Bereich 1")
       
    # für den zweiten Bereich

    params_2, cov_2 = curve_fit(lin, t1[grenz_1:grenz_2], noms(logung(x, grenz_1, grenz_2)))

    cov_2 = np.sqrt(np.diag(cov_2))
    print("2. Fit:\n",f"m = {params_2[0]:.4f} \pm {cov_2[0]:.5f} 1/s \t n = {params_2[1]:.4f} \pm {cov_2[0]:.5f}")
    uparams_2 = unp.uarray(params_2, cov_2)
    S_2 = -1 * V * uparams_2[0]
    print(f"Saugvermögen: {noms(S_2):.4f} \pm {stds(S_2):.5f}")
    #rel_abw(theo,S_2)
    plt.errorbar(t1[grenz_1:grenz_2], noms(logung(x, grenz_1,grenz_2)), yerr= stds(logung(x, grenz_1, grenz_2)), fmt='g.',label= "Messdaten Bereich 2")
    plt.plot(t1[grenz_1-1:grenz_2], lin(t1[grenz_1-1:grenz_2],*params_2), color='green',label="Fit Bereich 2")
    
    # für den dritten Bereich

    params_3, cov_3 = curve_fit(lin, t1[grenz_2:np.size(mean)], noms(logung(x, grenz_2, np.size(mean))))

    cov_3 = np.sqrt(np.diag(cov_3))
    print("3. Fit:\n",f"m = {params_3[0]:.4f} \pm {cov_3[0]:.5f} \t n = {params_3[1]:.4f} \pm {cov_3[0]:.5f}")
    uparams_3 = unp.uarray(params_3, cov_3)
    S_3 = -1 * V * uparams_3[0]
    print(f"Saugvermögen: {noms(S_3):.4f} \pm {stds(S_3):.5f}")
    #rel_abw(theo,S_3)
    plt.errorbar(t1[grenz_2:np.size(mean)], noms(logung(x, grenz_2,np.size(mean))), yerr= stds(logung(x, grenz_2,np.size(mean))), fmt='b.',label= "Messdaten Bereich 3")
    plt.plot(t1[grenz_2-1:np.size(mean)], lin(t1[grenz_2-1:np.size(mean)],*params_3), color='blue',label="Fit Bereich 3")

    plt.xlabel(r"$t$ $ [s]$")
    plt.ylabel(r"$\ln(\frac{p-p_{E}}{p_{0}-p_{E}})$")
    plt.tight_layout()
    plt.legend(loc = 'best')
    plt.savefig("build/plots/plot_" + name + ".pdf")
    #err_mean =mean
    #mean = err_mess(mean)
    #print(stds(err_mean - mean))
    
    saug = [S_1,S_2,S_3]
    return saug

def drehsaugung(t,a,b,c,name):
    grenz_1 = 21 
    grenz_2 = 35 #Arrayindizes 
    t1=t
    
    V = unc.ufloat(34,3.4)
    #theo = 1.1
    mean = mittel(a,b,c)
    x=fehler_dreh(mean)
    plt.clf()
    # für den ersten Bereich

    params_1, cov_1 = curve_fit(lin, t1[:grenz_1], noms(logungdreh(x, 0, grenz_1)))

    cov_1 = np.sqrt(np.diag(cov_1))
    print("1. Fit:\n",f"m = {params_1[0]:.4f} \pm {cov_1[0]:.5f} \t n = {params_1[1]:.4f} \pm {cov_1[0]:.5f}")
    uparams_1 = unp.uarray(params_1, cov_1)
    S_1 = -1 * V * uparams_1[0]
    print(f"Saugvermögen: {noms(S_1):.4f} \pm {stds(S_1):.5f}")
    #rel_abw(theo,S_1)
    plt.errorbar(t1[:grenz_1], noms(logungdreh(x, 0, grenz_1)), yerr= stds(logungdreh(x, 0, grenz_1)), fmt='k,',label= "Messdaten")
    plt.plot(t1[:grenz_1], lin(t1[:grenz_1],*params_1), color='red',label="Fit Bereich 1")
       
    # für den zweiten Bereich

    params_2, cov_2 = curve_fit(lin, t1[grenz_1:grenz_2], noms(logungdreh(x, grenz_1, grenz_2)))

    cov_2 = np.sqrt(np.diag(cov_2))
    print("2. Fit:\n",f"m = {params_2[0]:.4f} \pm {cov_2[0]:.5f} 1/s \t n = {params_2[1]:.4f} \pm {cov_2[0]:.5f}")
    uparams_2 = unp.uarray(params_2, cov_2)
    S_2 = -1 * V * uparams_2[0]
    print(f"Saugvermögen: {noms(S_2):.4f} \pm {stds(S_2):.5f}")
    #rel_abw(theo,S_2)
    plt.errorbar(t1[grenz_1:grenz_2], noms(logungdreh(x, grenz_1,grenz_2)), yerr= stds(logungdreh(x, grenz_1, grenz_2)), fmt='k,')
    plt.plot(t1[grenz_1-1:grenz_2], lin(t1[grenz_1-1:grenz_2],*params_2), color='green',label="Fit Bereich 2")
    
    # für den dritten Bereich

    params_3, cov_3 = curve_fit(lin, t1[grenz_2:np.size(mean)], noms(logungdreh(x, grenz_2, np.size(mean))))

    cov_3 = np.sqrt(np.diag(cov_3))
    print("3. Fit:\n",f"m = {params_3[0]:.4f} \pm {cov_3[0]:.5f} \t n = {params_3[1]:.4f} \pm {cov_3[0]:.5f}")
    uparams_3 = unp.uarray(params_3, cov_3)
    S_3 = -1 * V * uparams_3[0]
    print(f"Saugvermögen: {noms(S_3):.4f} \pm {stds(S_3):.5f}")
    #rel_abw(theo,S_3)
    plt.errorbar(t1[grenz_2:np.size(mean)], noms(logungdreh(x, grenz_2,np.size(mean))), yerr= stds(logungdreh(x, grenz_2,np.size(mean))), fmt='k,')
    plt.plot(t1[grenz_2-1:np.size(mean)], lin(t1[grenz_2-1:np.size(mean)],*params_3), color='blue',label="Fit Bereich 3")

    plt.xlabel(r"$t$ $ [s]$")
    plt.ylabel(r"$\ln(\frac{p-p_{E}}{p_{0}-p_{E}})$")
    plt.tight_layout()
    plt.legend(loc = 'best')
    plt.savefig("build/plots/plot_" + name + ".pdf")
    
    saug = [S_1,S_2,S_3]
    return saug

def leckrate(t,a,b,c,name,p):

    mean = mittel(a,b,c)
    
    if(name == "turbo"):
        V = unc.ufloat(33,3.3)
        theo = 77
        x=fehler_turbo(mean)
    else:
        V = unc.ufloat(34,3.4)
        theo = 1.1
        x=fehler_dreh(mean)
    
    params_1, cov_1 = curve_fit(lin, t, noms(mean))
    cov_1 = np.sqrt(np.diag(cov_1))
    print(f"Leck ", name , p ,f"mbar Fit:\n",f"m = {params_1[0]:.8f} \pm {cov_1[0]:.8f} \t n = {params_1[1]:.8f} \pm {cov_1[0]:.8f}")
    uparams = unp.uarray(params_1, cov_1)
    S = V / p *uparams[0]
    print(f"Saugvermögen {name}_{p}: {noms(S):.4f} \pm {stds(S):.5f}")
    #rel_abw(theo,S)
    plt.clf()
    plt.errorbar(t, noms(x), yerr= stds(x), fmt='r.',label= "Messdaten")
    plt.plot(t, lin(t,*params_1), color='green',label="Ausgleichsgerade")

    plt.xlabel(r"$t$ $ [s]$")
    plt.ylabel(r"$p$ $ [mbar]$")
    plt.tight_layout()
    plt.legend(loc = 'best')
    plt.savefig(f"build/plots/leck_{name}_{p}.pdf")

    return S

def saugdreh(log,a,b,c,d,s1,s2,s3,s4,slog):
    
    t_1 =np.linspace(0.1,1000,1000)
    dreh_theo = np.ones(1000) * 1.1

    plt.figure()
    plt.plot(t_1,dreh_theo, label = "Theoriewert" )
    plt.errorbar((a[0]+a[-1])/2, noms(s1), yerr=stds(s1), fmt='x', label = "Leck 0.5 mbar")
    plt.errorbar((b[0]+b[-1])/2, noms(s2), yerr=stds(s2), fmt='x', label = "Leck 9.7 mbar")
    plt.errorbar((c[0]+c[-1])/2, noms(s3), yerr=stds(s3), fmt='x', label = "Leck 50 mbar")
    plt.errorbar((d[0]+d[-1])/2, noms(s4), yerr=stds(s4), fmt='x', label = "Leck 100 mbar")
    plt.errorbar((log[0]+log[21])/2, noms(slog[0]),yerr=stds(slog[0]), fmt='x', label = "Evakuierung Druckbereich 1")
    plt.errorbar((log[21]+log[35])/2, noms(slog[1]),yerr=stds(slog[1]), fmt='x', label = "Evakuierung Druckbereich 2")
    plt.errorbar((log[35]+log[-1])/2, noms(slog[2]),yerr=stds(slog[2]), fmt='x', label = "Evakuierung Druckbereich 3")
    plt.xscale('log')
    plt.rc('axes', labelsize= size_label)
    plt.xlabel(r"$p$ $ [mbar]$")
    plt.ylabel(r"$S$ $ [\frac{l}{s}]$")
    plt.tight_layout()
    plt.legend(loc = 'best')
    plt.savefig("build/plots/saug_dreh.pdf")


def saugturbo(log,a,b,c,d,s1,s2,s3,s4,slog):
    
    t_1 =np.linspace(10**(-3),0.01,1000)
    theo = np.ones(1000) * 77

    plt.figure()
    plt.plot(t_1,theo, label = "Theoriewert" )
    plt.errorbar((a[0]+a[-1])/2, noms(s1), yerr=stds(s1), fmt='x', label = "Leck 1e-4 mbar")
    plt.errorbar((b[0]+b[-1])/2, noms(s2), yerr=stds(s2), fmt='x', label = "Leck 2e-4 mbar")
    plt.errorbar((c[0]+c[-1])/2, noms(s3), yerr=stds(s3), fmt='x', label = "Leck 5e-5 mbar")
    plt.errorbar((d[0]+d[-1])/2, noms(s4), yerr=stds(s4), fmt='x', label = "Leck 7e-5 mbar")
    plt.errorbar((log[0]+log[5])/2, noms(slog[0]),yerr=stds(slog[0]), fmt='x', label = "Evakuierung Druckbereich 1")
    plt.errorbar((log[5]+log[10])/2, noms(slog[1]),yerr=stds(slog[1]), fmt='x', label = "Evakuierung Druckbereich 2")
    plt.errorbar((log[10]+log[-1])/2, noms(slog[2]),yerr=stds(slog[2]), fmt='x', label = "Evakuierung Druckbereich 3")
    plt.xscale('log')
    plt.rc('axes', labelsize= size_label)
    plt.xlabel(r"$p$ $ [mbar]$")
    plt.ylabel(r"$S$ $ [\frac{l}{s}]$")
    plt.tight_layout()
    plt.legend(loc = 'best')
    plt.savefig("build/plots/saug_turbo.pdf")


#------------------Ausgabe----------------------------------------

#turboprinter(t1,p_1T,p_2T,p_3T) 
#turboprinter(t, p_1TL1, p_2TL1, p_3TL1)
#turboprinter(t, p_1TL2, p_2TL2, p_3TL2)
#turboprinter(t, p_1TL5, p_2TL5, p_3TL5)
#turboprinter(t, p_1TL7, p_2TL7, p_3TL7)
#drehprinter(t2, p_1D, p_2D, p_3D)
#drehprinter(t3, p_1DL1, p_2DL1, p_3DL1)
#drehprinter(t3, p_1DL9, p_2DL9, p_3DL9)
#drehprinter(t3, p_1DL50, p_2DL50, p_3DL50)
#drehprinter(t3, p_1DL100, p_2DL100, p_3DL100)
#logprinter(t2, p_1D, p_2D, p_3D,'dreh')
#logprinter(t1,p_1T,p_2T,p_3T,'turbo')

#saug_turbo_ev=saugung(t1,p_1T,p_2T,p_3T,'ev_turbo')
#saug_dreh_ev=drehsaugung(t2, p_1D, p_2D, p_3D,'ev_dreh')

#saug_turbo_1 = leckrate(t, p_1TL1, p_2TL1, p_3TL1,'turbo',1*10**(-4))
#saug_turbo_2 = leckrate(t, p_1TL2, p_2TL2, p_3TL2,'turbo',2*10**(-4))
#saug_turbo_5 = leckrate(t, p_1TL5, p_2TL5, p_3TL5,'turbo',5*10**(-5))
#saug_turbo_7 = leckrate(t, p_1TL7, p_2TL7, p_3TL7,'turbo',7*10**(-5))

#saug_dreh_1 = leckrate(t3, p_1DL1, p_2DL1, p_3DL1,'dreh',0.5)
#saug_dreh_9 = leckrate(t3, p_1DL9, p_2DL9, p_3DL9,'dreh',9.7)
#saug_dreh_50 = leckrate(t3, p_1DL50, p_2DL50, p_3DL50,'dreh',50)
#saug_dreh_100 = leckrate(t3, p_1DL100, p_2DL100, p_3DL100,'dreh',100)

#saugdreh(p_1D,p_1DL1,p_1DL9,p_1DL50,p_1DL100,saug_dreh_1,saug_dreh_9,saug_dreh_50,saug_dreh_100,saug_dreh_ev)
#saugturbo(p_1T,p_1TL1,p_1TL2,p_1TL5,p_1TL7,saug_turbo_1,saug_turbo_2,saug_turbo_5,saug_turbo_7,saug_turbo_ev)