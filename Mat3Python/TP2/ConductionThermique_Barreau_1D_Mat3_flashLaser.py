
# -------------------------------------------------------------------------
#                   Programme Python
#      L'équation de la chaleur en conduction dans un barreau
#                  dT/dt - D d²T/dx² = 0
# -------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm


# -------------------------------------------------------------------------
#                 Choix des paramètres
# -------------------------------------------------------------------------

L=0.01  # La longueur du barreau

ne=100 # Nombre d'élement

dx = L/ne   #La longueur de chaque élément

x = np.linspace(0, L-dx, ne)+0.5*dx # Maillage pour l'espace : Coordonnées des points


tfin = 0.01 # Temps final
dt = 1e-4 # Pas de temps
Nt=int(tfin/dt)
t = np.linspace(0, tfin, int(tfin/dt)+1)    # Maillage pour le temps



# -------------------------------------------------------------------------
#                 Choix des matériaux
# -------------------------------------------------------------------------

# cp: Capacité thermique massique (J K−1 kg−1) 
# rho: Masse volumique ( kg/m3 ) 
# conduct_thermique : Conductivité thermique ( W m−1 K−1 ) 

# Fer
cp= 444; rho= 7860; conduct_thermique= 80 

# Acier	
#cp= 435; rho= 7500; conduct_thermique=50

# Aluminium		
#cp= 897; rho=2700; conduct_thermique=220

# Cuivre	
#cp= 385; rho= 8960; conduct_thermique=368

# Bois		
#cp= 2000; rho= 800; conduct_thermique=0.04


D= conduct_thermique/(rho*cp) # diffusivité thermique (m2/s)
F = D*dt/(dx**2)



# -------------------------------------------------------------------------
#                  Les conditions initiales
# -------------------------------------------------------------------------

# Conditions initiales:
Tini=np.ones(ne)*20; Tini[0]=60

plot_Tinit=0
if plot_Tinit:
    plt.figure()    
    plt.plot(x, Tini, '-*b')
    plt.grid()
    plt.xlabel('Distance (m)')
    plt.ylabel(r'Température ($^\circ C$)')
    plt.title('Température initiale')
    #plt.savefig('Température initiale.png',dpi=50,format='png',transparent=False)
   

# -------------------------------------------------------------------------
#         L'evolution temporelle de la température dans le barreau
# -------------------------------------------------------------------------


#print ('Tini=',Tini)


    
Tavant=Tini
Tapres=np.zeros(ne)


tableau_temps_Tx=np.zeros((Nt+1, ne))
tableau_temps_Tx[0,:]=Tini

for i in range(0, Nt):  
    # Calculer Tapres aux points de maille externes : point 0
    Tapres[0] = Tavant[0] + F*(Tavant[1] - Tavant[0])
    
    # Calculer Tapres aux points de maille internes : points 1-> ne-1
    for j in range(1, ne-1):
        Tapres[j] = Tavant[j] + F*(Tavant[j-1] - 2*Tavant[j] + Tavant[j+1])

    # Calculer Tapres aux points de maille externes : point ne
    
    Tapres[ne-1] = Tavant[ne-1] - F*(Tavant[ne-1] - Tavant[ne-2])   
    
    # Préciser les conditions limites (?)
    #if t[i]<0.2:  
        #Tapres[0] = 60; 
    
    # Mettre à jour Tavant (avant le prochain pas de temps)
    Tavant= Tapres
    tableau_temps_Tx[i+1,:]=Tapres
   
    

#print ('Tapres=',Tapres)

# -------------------------------------------------------------------------
#             Traitement des résultats 
# -------------------------------------------------------------------------


t_index=np.array([0,0.04,0.2,0.25,0.5,1])*Nt
t_index = t_index.astype(int)

color = iter(cm.rainbow(np.linspace(0, 1,len(t_index))) )


plt.figure() 


for k in t_index:
    plt.plot(x*1000, tableau_temps_Tx[k,:], '-',color=next(color),markersize=5,label='t='+str(t[k])+'s')  # Tracé de la température   
     
plt.grid()
plt.legend(loc=0)    
plt.xlabel('Distance (mm)')
plt.ylabel(r'Température ($^\circ C$)')
plt.title('Evolution de la température')





plt.figure() 


plt.plot(t, tableau_temps_Tx[:,0], '-g*',markersize=5, label='Le premier élément')  # Tracé de la température   
plt.plot(t, tableau_temps_Tx[:,-1], '-r*',markersize=5, label='Le dernier élément')  # Tracé de la température   
     
plt.grid()
plt.legend(loc=0)    
plt.xlabel('Temps (s)')
plt.ylabel(r'Température ($^\circ C$)')
plt.title('Les températures au sein des éléments')

values=tableau_temps_Tx[:,1]

oo=np.argmax(max(values))
