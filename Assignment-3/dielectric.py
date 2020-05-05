import numpy as np
import matplotlib.pyplot as plt

def wl_plasma(n, m):
    omega_p = np.sqrt(n * e**2.0 / (e_inf * e_0 * m))
    wl_p = 2.0*np.pi*c/omega_p

    return wl_p

def epsilon_wl(wl, wl_p):
    epsilon = e_inf*(1.0 - (wl/wl_p)**2.0)

    return epsilon

def real_reflectance(epsilon):
    n = np.sqrt(epsilon)
    R = (n-1)**2.0/((n+1)**2.0)

    return R

def main():
    global e_inf
    e_inf = 1.0
    global e
    e = 1.6e-19
    global e_0
    e_0 = 8.85e-12
    global m_e
    m_e = 9.11e-31
    global c
    c = 3e8

    m_star = 1.0
    m = m_star * m_e

    n_800 = 1.741e27
    n_900 = 1.376e27   # m^-3
    n_1000 = 1.114e27
    n_1100 = 9.21e26

    wl_max = 801e-9
    wl_min = 350e-9

    wl_range = np.linspace(wl_min, wl_max, 3000)

    wl_p_800 = wl_plasma(n_800, m)
    wl_p_900 = wl_plasma(n_900, m)
    wl_p_1000 = wl_plasma(n_1000, m)
    wl_p_1100 = wl_plasma(n_1100, m)

    epsilon_wl_l_800 = epsilon_wl(wl_range, wl_p_800)
    reflectance_l_800 = real_reflectance(epsilon_wl_l_800)

    epsilon_wl_l_900 = epsilon_wl(wl_range, wl_p_900)
    reflectance_l_900 = real_reflectance(epsilon_wl_l_900)

    epsilon_wl_l_1000 = epsilon_wl(wl_range, wl_p_1000)
    reflectance_l_1000 = real_reflectance(epsilon_wl_l_1000)

    epsilon_wl_l_1100 = epsilon_wl(wl_range, wl_p_1100)
    reflectance_l_1100 = real_reflectance(epsilon_wl_l_1100)

    plt.plot(wl_range, reflectance_l_800, 'y', label='1.741e21 cm^-3')
    plt.plot(wl_range, reflectance_l_900, 'r', label='1.376e21 cm^-3')
    plt.plot(wl_range, reflectance_l_1000, 'b', label='1.114e21 cm^-3')
    plt.plot(wl_range, reflectance_l_1100, 'g', label='9.210e20 cm^-3')
    plt.xlabel('Wavelength (m)')
    plt.ylabel('Reflectance')
    plt.legend(loc='best')
    plt.show()

main()
