import numpy as np 
import matplotlib.pyplot as plt 
from scipy.interpolate import CubicSpline 
from scipy.special import sph_harm, gamma
from matplotlib.colors import LogNorm

def load_spectrum(l_nr):  
    spec = np.loadtxt('spectrum_data/partial_spectrum_'+ str(l_nr) + '.txt')
    energy = np.loadtxt('eigenstate_data/energies_continuum_' + str(l_nr) + '.txt')
    return spec[:,0] + 1j*spec[:,1], energy 

l_max = 10

energy_l0 = np.loadtxt('eigenstate_data/energies_continuum_0.txt')

spline_list = []
spline_list_r = []
spline_list_i = []

fig, ax = plt.subplots()

for i in range(l_max+1): 
    spec_l, energy_l = load_spectrum(i)
    k_l = np.sqrt(2*energy_l)
    #spec_l *= np.sqrt(k_l)  # Scale with sqrt(k) to change from energy to momentum normalized
    
    # Add Coulomb phase shift and momentum normalization
    phase = np.angle(gamma(i + 1 - 1j/k_l))
    spec_l *= (-1j)**i * np.exp(1j * phase) #/ k_l  

    # Spline it to interpolate between different k values later
    spline = CubicSpline(k_l, spec_l, extrapolate=False, bc_type='natural')
    spline_r = CubicSpline(k_l, np.real(spec_l), extrapolate=False, bc_type='natural')
    spline_i = CubicSpline(k_l, np.imag(spec_l), extrapolate=False, bc_type='natural')
    spline_list.append(spline)
    spline_list_r.append(spline_r)
    spline_list_i.append(spline_i)
    
    if i < 1: 
        plt.plot(k_l, np.real(spec_l), 'o', c=f'C{i}')
        plt.plot(k_l, np.imag(spec_l), 'x', c=f'C{i}')

E_list = np.linspace(min(energy_l0), max(energy_l0), 2000)

plt.plot(np.sqrt(2*E_list), spline_list_r[0](np.sqrt(2*E_list)))


k_max = 1#k_list[-1] 
k_min = np.sqrt(2*E_list[0])

print(k_max, k_min)

# For now Cartesian? 
kz_list = np.linspace(-k_max, k_max, 300)
kx_list = np.linspace(-1, 1, 300) 

KZ, KX = np.meshgrid(kz_list, kx_list, indexing='ij')

PMD = np.zeros_like(KX, dtype=complex)

for i, kz in enumerate(kz_list): 
    for j, kx in enumerate(kx_list): 
        # First determine angular coordinates 
        k = np.sqrt(kx**2 + kz**2)
        theta = np.arccos(kz / k)

        # Check if within range 
        if k <= k_min: 
            PMD[i,j] = 0.
            continue 

        # Add up the projektion on the Coulomn WF 
        PMD_i = 0.
        for li in range(l_max+1): 
            spline_i = spline_list_r[li](k) + 1j*spline_list_i[li](k)
            if spline_i != spline_i: 
                continue 
            PMD_i +=  spline_i * sph_harm(0, li, 0., theta)  
        PMD[i,j] = PMD_i / np.sqrt(k)

# Now plot the stuffs? 
PMD = np.abs(PMD)**2 
np.nan_to_num(PMD, copy=False)
print('Max val : ', np.max(PMD))
#PMD = PMD / np.max(PMD)

log_order = 6
fig, ax = plt.subplots()
ax.axis('equal')
im = ax.pcolormesh(KZ, KX, PMD, norm=LogNorm(vmax=1, vmin=1*10**(-log_order)), cmap='turbo', shading='gouraud')

plt.colorbar(im)
plt.show()


