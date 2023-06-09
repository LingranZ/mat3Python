import numpy as np
import matplotlib.pyplot as plt



def chaleursphere(r,rho,Cp,DT):
    V=4/3*np.pi*r**3
    Q=rho*V*Cp*DT
    return (Q)
    
def chaleurplaque(e,D, rho,Cp,DT):
    V=e*D*D
    Q=rho*V*Cp*DT
    return (Q)

z=np.linspace(0.1,0.5, 100)
DT=140




# Acier	
Cp= 435; rho= 7500; 
plt.plot(z,chaleursphere(z,rho,Cp, DT),'-g*',label='Sphère acier')


# Béton	
Cp= 840; rho= 2350; D=0.8;
plt.plot(z,chaleurplaque(z,D,rho,Cp, DT),'-cs',label='plaque béton')




plt.grid()
plt.legend(loc=0)    
plt.xlabel('rayon sphère ou épaisseur plaque (m)')
plt.ylabel('Chaleur (J)')

plt.title("Chaleur en joule emmagasinée \n dans une sphère en acier ou une plaque en béton")

plt.axis([0.3025,0.3050, 0.535e8, 0.54e8])