import numpy as np

lambd = 800  # nm
I0 = 1. 
E_max = np.sqrt(I0 * 1.0e14 / 3.50945e16)  # Maximum electric field amplitude in a.u.

omega = 2.*np.pi * 137.036 / (lambd * 1.0e-9 / 5.29177e-11)
Up = E_max**2 / (4. * omega**2)  # Ponderomotive energy in a.u
Ip = 0.5 
Nc = 2
E_conv = 12*Up

laser_duration = 2*np.pi*Nc / omega 

r_max = 0.55 * laser_duration * np.sqrt(2*E_conv)
delta_r = 1/3 * 2*np.pi / np.sqrt(2*E_conv)
dt = 2*np.pi/omega / 50 #np.sqrt(12*0.05*omega / (E_conv**3))

min_photon = 0
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