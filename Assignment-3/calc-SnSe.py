

import numpy as np
import matplotlib.pyplot as plt

def intrinsic_carrier_conc(kbT, m_n, m_p, E_g):
    t1 = ((kbT)/(2.0 * np.pi * hbar**2.0))**(3.0/2.0)
    t2 = (m_n*m_p)**(3.0/4.0)
    t3 = np.exp(-E_g/(2.0*kbT))
    return (2.0 * t1 * t2 * t3)

def effective_conc(kbT, m_x):
    return (2.0 * ( (m_x * kbT)/(2.0 * np.pi * hbar**2.0) )**(3.0/2.0))

def extrinsic_carrier_conc(kbT, m_x, N_X, E_x):
    return (2.0 * N_X * ( 1.0 + ( 1.0 + (4.0 * N_X/effective_conc(kbT, m_x) * np.exp(E_x/(kbT))) )**(1.0/2.0) )**(-1.0))

def main():
    global e
    e = 1.6e-19
    global c
    c = 3e8   # ms^-1
    global ev
    eV = 1.6e-19
    global m_e
    # m_e = 0.511e6 # eV c^-2
    m_e = 9.10938356e-31  # kg
    global k_b
    # k_b = 8.617e-5   # eV K^-1
    k_b = 1.38064852e-23  # J K^-1
    global hbar
    # hbar = 6.582e-16   # eV s
    hbar = 1.0545718e-34  # J s

    #   SnSe
    E_g = 0.8 * eV # eV
    m_n = 0.76 * m_e   # eV
    m_p = 1.5 * m_e   # eV
    N_A = 1e17 * 1e6   # m^-3
    E_a = 0.1 * eV  # eV
    T = 300.0   # K
    mu_p = 200  # cm^2 V^-1 s^-1

    T_range = np.linspace(100,1000,100000)
    carrier_concs = np.array([])

    for T in T_range:
        kbT = k_b*T
        pi = intrinsic_carrier_conc(kbT, m_n, m_p, E_g) * 1e-6
        pe = extrinsic_carrier_conc(kbT, m_p, N_A, E_a) * 1e-6

        total_carrier_conc = 2*pi + pe
        carrier_concs = np.append(carrier_concs, total_carrier_conc)

    plt.plot(1.0/T_range, carrier_concs, 'black')
    plt.vlines(1.0/100, 0, 1e20, label='100K', color='red', linestyle='dotted')
    plt.vlines(1.0/300, 0, 1e20, label='300K', color='green', linestyle='dotted')
    plt.vlines(1.0/1000, 0, 1e20, label='1000K', color='blue', linestyle='dotted')
    plt.yscale('log')

    plt.title('Total Carrier Concentration of SnSe vs Reciprocal Temperature')
    plt.xlabel(r'T$^{-1}$ /K$^{-1}$')
    plt.ylabel(r'Total Carrier Concentration /cm$^{-3}$')
    plt.legend(loc='best')

    plt.ylim([0,1e19])
    plt.show()

main()
