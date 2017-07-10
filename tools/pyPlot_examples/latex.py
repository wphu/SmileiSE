# -*- coding: utf-8 -*-

from __future__ import unicode_literals
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['text.usetex'] = True
#matplotlib.rcParams['text.latex.unicode'] = True



plt.figure(1, figsize=(10, 6))
ax = plt.axes([0.1, 0.1, 0.8, 0.7])
t = np.arange(0.0, 1.0 + 0.01, 0.01)
s = np.cos(2*2*np.pi*t) + 2
plt.plot(t, s)

plt.xlabel(r'x label (mm): ' r'$\mu=100\ \mu s, \sigma=15^{10}\ m^{-3}$', fontsize=16)
plt.ylabel(r'Velocity $(m/s)$', fontsize=16)
#plt.title(r'\TeX\ is Number $\displaystyle\sum_{n=1}^\infty'
#          r'\frac{-e^{i\pi}}{2^n}$!', fontsize=16, color='r')
plt.title(r"$\phi\ (kT_e)$",fontsize=18,fontweight='bold')
plt.grid(True)
#plt.savefig("latex.png", dpi=300)
plt.savefig("latex.pdf", dpi=300)
#plt.show()
