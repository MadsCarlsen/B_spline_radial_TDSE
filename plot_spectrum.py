import numpy as np 
import matplotlib.pyplot as plt 
from scipy.interpolate import CubicSpline 

def load_spectrum(l_nr):  
    spec = np.loadtxt('spectrum_data/partial_spectrum_'+ str(l_nr) + '.txt')
    energy = np.loadtxt('eigenstate_data/energies_continuum_' + str(l_nr) + '.txt')
    return np.abs(spec[:,0] + 1j*spec[:,1])**2, energy 

l_max = 10

energy_l0 = np.loadtxt('eigenstate_data/energies_continuum_0.txt')

E_list = np.linspace(min(energy_l0), max(energy_l0), 2000)
tot_spec = np.zeros_like(E_list)

fig, ax = plt.subplots()

for i in range(l_max+1): 
    # Make spline of each partial spectrum and add them together 
    spec_l, energy_l = load_spectrum(i)
    spline = CubicSpline(energy_l, spec_l, extrapolate=False, bc_type='natural')

    temp = spline(E_list)
    np.nan_to_num(temp, copy=False)
    tot_spec += temp 

    if i in [0,1,4]: 
        #ax.plot(energy_l, spec_l, c=f'C{i+1}')
        ax.plot(E_list, temp, 'o', c=f'C{i+1}')
        ax.plot()

ax.plot(E_list, tot_spec, c='C0')
ax.set_yscale('log')

# Plot the total spectrum 
lambd = 532 # nm
I0 = 0.2
E_max = np.sqrt(I0 * 1.0e14 / 3.50945e16)  # Maximum electric field amplitude in a.u.

#omega = 0.5 - 0.125 
omega = 2.*np.pi * 137.036 / (lambd * 1.0e-9 / 5.29177e-11)
Up = E_max**2 / (4. * omega**2)  # Ponderomotive energy in a.u
Ip = 0.5 
Nc = 2
E_conv = 12*Up

laser_duration = 2*np.pi*Nc / omega 

r_max = 0.55 * laser_duration * np.sqrt(2*E_conv)
delta_r = 1/3 * 2*np.pi / np.sqrt(2*E_conv)
dt = 2*np.pi/omega / 50 #np.sqrt(12*0.05*omega / (E_conv**3))


for i in range(100): 
    if i*omega > Ip + Up: 
        min_photon = i 
        break 

print('Up : ', Up)
print('Min. photons : ', min_photon)
print('Laser duration : ', laser_duration, laser_duration * 24.1 *1e-3)
print('r_max : ', r_max)
print('delta_r : ', delta_r)
print('dt : ', dt)

fig, ax = plt.subplots()
ax.plot(E_list , tot_spec) 
ax.set_yscale('log')

plt.axvline(2*Up, ls='--', c='k')
plt.axvline(10*Up, ls='--', c='k')

max_val_index = np.argmax(tot_spec)



for i in range(min_photon, min_photon + 11): 
    #plt.axvline((E_list[max_val_index] + omega*i), ls='--', c='gray', alpha=0.5)
    plt.axvline((-Ip - Up + omega*i), ls='--', c='gray', alpha=0.5)

plt.show()