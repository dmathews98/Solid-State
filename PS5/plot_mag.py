import numpy as np
import matplotlib.pyplot as plt

def Curie(C, T):
    return (float(C)/T)

def twolevel(N, mu, B, T):
    kb = 1.0
    return (N*kb*np.tanh(mu*B/(kb*T)))

def main():
    N = 1
    mu = 1
    B = 1
    C = 1

    t = np.arange(0.0001, 50, 0.001)

    curie_dat = Curie(C, t)
    twolev_dat = twolevel(N, mu, B, t)

    #t = 1/t

    plt.plot(t, curie_dat)
    plt.plot(t, twolev_dat)
    plt.yscale('log')
    plt.show()


main()
