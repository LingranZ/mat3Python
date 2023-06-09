import numpy as np
import matplotlib.pyplot as plt

r=0.3
h=15

s=4*np.pi*r**2

sigma=5.67*1e-8 # W/(m*m*K)

T= np.linspace(20,300, 500)
Tair=20

PHI_conv=h*(T-Tair)*s
PHI_ray= sigma * (  (T+273)**4-(Tair+273)**4  )*s
PHI_ray_approx= sigma * (T+273)**4 *s


plt.figure(figsize=(8, 6)) 

plt.plot(T-Tair, PHI_conv, 'r-',label='Convection',lw=2)    
plt.plot(T-Tair, PHI_ray, 'g-',label='Rayonnement',lw=2)    
plt.plot(T-Tair, PHI_ray_approx, 'b-',label='Rayonnement (approximation)',lw=2)    
    
plt.grid(ls='--')
plt.legend(loc=0)    
plt.xlabel('T-Tair ($^\circ C$)')
plt.ylabel(r'Flux $\Phi $ (w)')
#plt.axis([50,70,0,4000])
