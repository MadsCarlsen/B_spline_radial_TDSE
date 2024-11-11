import numpy as np 
import matplotlib.pyplot as plt
from scipy.io import FortranFile 
from matplotlib.colors import LogNorm

# Loads a binary fortran data file 
def load_fortran_binary(path, dim=1): 
    f = FortranFile(path, 'r')
    data = f.read_reals(float)
    dim1 = int(len(data)/dim) 
    data = data.reshape((dim1, dim), order="F")  # Reshaping for 2D arrays 
    f.close()
    return data

E_list = np.loadtxt('spectrum_data/window_energies.txt')

l_max = 25

tot_spec = np.zeros_like(E_list)

fig, ax = plt.subplots()
for li in range(l_max+1): 
    spec_l = 0.001**4 * np.loadtxt(f'spectrum_data/window_partial_spectrum_{li}.txt')
    tot_spec += spec_l
    ax.plot(E_list, spec_l)
ax.set_yscale('log')

fig, ax = plt.subplots()
ax.plot(E_list, tot_spec)
ax.set_yscale('log')

#plt.show()


# PLOT THE PMD! 
theta_list = np.loadtxt('spectrum_data/theta_list.txt')
E_pos = E_list[E_list > 0]
k_list = np.sqrt(2*E_pos)

THETA, K = np.meshgrid(theta_list, E_pos, indexing='ij')
Nt = len(theta_list)
Nk = len(k_list)

X = K * np.sin(THETA)
Z = K * np.cos(THETA)

PMD = load_fortran_binary('spectrum_data/PMD.dat', Nk)
PMD /= K  # Normalize to momentum scale 
print(np.max(PMD))
#PMD /= np.max(PMD)

fig, ax = plt.subplots()
ax.axis('equal')
log_order = 3
im = ax.pcolormesh(Z, X, PMD, norm=LogNorm(vmax=10**9, vmin=1*10**(-log_order)), cmap='viridis', shading='gouraud')#, shading='gouraud')#, shading='gouraud')

cbar = fig.colorbar(im, ax=ax)#, aspect=10, fraction=0.039, pad=0.15)
ax.set_ylim(0, np.max(X))
ax.set_xlim(np.min(Z), np.max(Z))

plt.tight_layout()
plt.show()