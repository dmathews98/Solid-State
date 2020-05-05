'''
Fitting Iridium data to Murnaghan and Vinet EoS
'''

from __future__ import division
import numpy as np
import lmfit
import matplotlib.pyplot as plt
import math as m

    # x = V
def PV_MurnaghanEoS(x, Vo, B, Bp):
    P = (B/Bp)*((x/Vo)**(-1.0*Bp) - 1.0)
    return P

def PV_VinetEoS(x, Vo, B, Bp):
    z = (x/Vo)**(1.0/3.0)
    P = 3.0*B*((1.0-z)/(z**2.0))*np.exp(1.5*(Bp - 1.0)*(1.0 - z))
    return P

    # x = P
def VP_MurnaghanEoS(x, Vo, B, Bp):
    V = Vo*((((Bp/B)*x) + 1.0)**(-1.0/Bp))
    return V

def fit_method(method, xvals, yvals, Vo_g, B_g, Bp_g):
    model = lmfit.Model(method)

    parameters = model.make_params(Vo=Vo_g, B=B_g, Bp=Bp_g)
    parameters['Vo'].vary=False

    result = model.fit(yvals, x=xvals, params=parameters)

    print(result.fit_report())
    print('\n')

    return model, result

def main():
    file = 'Iridium-data.txt'

    # read in data
    xvals = [] # Volume / Angstroms^3
    yvals = [] # Pressure /GPa
    f = open(file, 'r')
    for line in f.readlines():
        vals  = line.split(' , ')
        xvals.append(float(vals[1]))
        yvals.append(float(vals[0]))
    f.close()

    # a~3.839 angstrom so UC vol roughly 56.579 ang^3
    # experiment results give it as 55.059 at 0 GPa
    # bulk mod unknown so guess Bp = 4 GPa and B = 5 GPa
    Vo_g = 55.059
    B_g = 300.0
    Bp_g = 5.0

    Method = PV_MurnaghanEoS
    M_model, M_result = fit_method(Method, xvals, yvals, Vo_g, B_g, Bp_g)

    Method = PV_VinetEoS
    V_model, V_result = fit_method(Method, xvals, yvals, Vo_g, B_g, Bp_g)

    Method = VP_MurnaghanEoS
    Mvp_model, Mvp_result = fit_method(Method, yvals, xvals, Vo_g, B_g, Bp_g)

    plt.figure(1)
    plt.plot(xvals, yvals, 'bo')

    xplot = np.linspace(47.7, 55.2, 1000)

    plt.plot(xplot, M_model.eval(M_result.params, x=xplot), 'r--')

    plt.plot(xplot, V_model.eval(V_result.params, x=xplot), 'g--')

    plt.title('Fitting of Murnaghan and Vinet P(V) EoS to experimental data')
    plt.legend(['Data', 'Murnaghan Fit', 'Vinet Fit'], loc='upper right')
    plt.xlabel(r'Volume $\AA$')
    plt.ylabel('Pressure GPa')
    plt.savefig('PVfits.pdf')

    plt.figure(2)
    plt.plot(yvals, xvals, 'bo')

    yplot = np.linspace(-0.1, 72.6, 1000)
    plt.plot(yplot, Mvp_model.eval(Mvp_result.params, x=yplot), 'r--')

    plt.title('Fitting of Murnaghan V(P) EoS to experimental data')
    plt.legend(['Data', 'Murnaghan Fit'], loc='upper right')
    plt.ylabel(r'Volume $\AA$')
    plt.xlabel('Pressure GPa')
    plt.savefig('VPfit.pdf')

    #plt.show()

main()
