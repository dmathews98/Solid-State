

from __future__ import division
import sys
import numpy as np
import matplotlib.pyplot as plt
import lmfit

def linear(x, m, b):
    return m*x + b

def fit_method(method, xvals, yvals, m_g, b_g):
    model = lmfit.Model(method)

    parameters = model.make_params(m=m_g, b=b_g)

    result = model.fit(yvals, x=xvals, params=parameters)

    print('')
    print(result.fit_report())
    print('')

    return model, result

if __name__ == "__main__":
    # I/O

    if len(sys.argv) != 3:
        print('Usage: python det_eg.py filename')
        quit()
    else:
        [prog, file, Tc] = sys.argv

    Tc = float(Tc)

    data = np.loadtxt(str(file))   # T^2 [K^2]  C/T [mJ /mol /K^2]

    T = np.sqrt(data[:,0])   # [K]

    # (a)
    # [mJ /mol /K]
    C_s = data[:,1] * T

    # (b)
    beta = 0.0568

    gamma = 0.596

    # [mJ /mol /K]
    C_ph = beta * T**3.0

    # [mJ /mol /K]
    C_se = C_s - C_ph

    # (c)
    # [mJ /mol /K]
    C_n = gamma * T + C_ph

    # (d)
    plt.figure(1)
    plt.plot(T, C_se, 'g-o', label=r'$C_{se}$', markersize=2.5)
    plt.plot(T, C_n, 'b-o', label=r'$C_n$', markersize=2.5)
    plt.plot(T, C_s, 'r-o', label=r'$C_s$', markersize=2.5)
    plt.xlabel('T (K)')
    plt.ylabel('Heat Capacity (mJ /mol /K)')
    plt.title('Heat Capacities of Galium v Temperature')
    plt.legend(loc='best')

    d_out = np.transpose(np.vstack([np.transpose(T), np.transpose(C_n), np.transpose(C_s), np.transpose(C_se)]))

    # np.savetxt('Ga_data.txt', d_out)

    print('')
    print('      T          C_n       C_s        C_se')
    print(d_out)
    print('')

    # (e)
    # [mJ /mol /K]
    C_ne_Tc = gamma * Tc

    print('C_ne_Tc = '+str(C_ne_Tc)+' mJ /mol /K')

    # (f)
    plt.figure(2)
    plt.grid()
    plt.plot(Tc/T, C_se/C_ne_Tc, 'r-o', markersize=2.5)
    plt.yscale('log')
    plt.xlabel(r'$\frac{Tc}{T}$')
    plt.ylabel(r'$\frac{C_{se}}{C_{ne}(T_c)}$')

    # (g)


    # (h)
    indices = np.where(T < Tc/2.0)
    max = max(indices[0])
    min = min(indices[0])+1

    # grad = (C_se[min]/C_ne_Tc - C_se[max]/C_ne_Tc)/(Tc/T[min] - Tc/T[max])

    m_g = 2.0
    b_g = 2.0

    model, res = fit_method(linear, Tc/T[min:max], C_se[min:max]/C_ne_Tc, m_g, b_g)

    grad = float(res.params['m'])
    print('Gradient '+str(grad))

    k = 1.38e-23
    e = 1.6e-19

    Eg = -grad*k/e * 1e3
    d0 = Eg/2.0

    print('Eg '+str(Eg)+' meV')

    print('Delta(0) '+str(d0)+' meV')

    twogap = (3.5 * k * Tc / e) * 1e3
    gap = twogap / 2.0
    print('BCS theory: Delta(0) '+str(gap)+' meV')

    plt.show()
