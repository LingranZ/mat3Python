
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

L=0.2  # La longueur du barreau
ne=101 # Nombre d'élement

dx = L/ne   #La longueur de chaque élément
x = np.linspace(0, L-dx, ne)+0.5*dx # Maillage pour l'espace : Coordonnées des points



tfin = 600  # Temps final
dt = 0.010 # Pas de temps
Nt=int(tfin/dt) # Nombre total de pas
t = np.linspace(0, tfin, Nt)    # Maillage pour le temps

Tchaud=40 # la température la plus chaude dans le barreau

# -------------------------------------------------------------------------
#                 Choix des matériaux
# -------------------------------------------------------------------------

# Lambda: Conductivité thermique ( W m−1 K−1 ) 
# cp: Capacité thermique massique (J K−1 kg−1) 
# rho: Masse volumique ( kg/m3 ) 


# Fer
cp= 444; rho= 7860; Lambda= 80 

# Acier	
#cp= 435; rho= 7500; Lambda=50

# Aluminium		
#cp= 897; rho=2700; Lambda=220

# Cuivre	
#cp= 385; rho= 8960; Lambda=368

# Bois		
#cp= 2000; rho= 800; Lambda=0.04


D= Lambda/(rho*cp) # diffusivité thermique (m2/s)
F = D*dt/(dx**2)



# -------------------------------------------------------------------------
#                  Les conditions initiales
# -------------------------------------------------------------------------

# Conditions initiales:
Tini=np.ones(ne)*20; 
Tini[int(ne/2)]=Tchaud

# -------------------------------------------------------------------------
#         L'evolution temporelle de la température dans le barreau
# -------------------------------------------------------------------------

    
Tavant=Tini
Tapres=np.zeros(ne)


tableau_temps_Tx=np.zeros((Nt, ne))
tableau_temps_Tx[0,:]=Tini

for i in range(0, Nt-1):  
    
    # Calculer Tapres aux points de maille internes
    for j in range(1, ne-1):
        Tapres[j] = Tavant[j] + F*(Tavant[j-1] - 2*Tavant[j] + Tavant[j+1])

    # Préciser les conditions limites
    #Tapres[int(ne/2)] = Tchaud; Tapres[-1] = 20;  Tapres[0] = 20#
    
    idm=int(ne/2)
    n_chaud=5
    Tapres[idm-n_chaud:idm+n_chaud+1] = Tchaud; Tapres[-1] = 20;  Tapres[0] = 20#
    
    # Mettre à jour Tavant (avant le prochain pas de temps)
    Tavant= Tapres
    tableau_temps_Tx[i+1,:]=Tapres
    
    

# -------------------------------------------------------------------------
#             Traitement des résultats 
# -------------------------------------------------------------------------

#Sélectionner plusieurs moments du temps de calcul 
t_index=np.array([0,0.01,0.05,0.1,1])*(Nt-1) 
t_index = t_index.astype(int)

#Générer plusieurs couleurs d'arc-en-ciel
color = iter(cm.rainbow(np.linspace(0, 1,len(t_index))) )


plt.figure() 


for k in t_index:
    # Tracé la courbe de la température, pour chaque moment sélectionné
    Label='t='+str( round(t[k],2)  )+'s'
    plt.plot(x, tableau_temps_Tx[k,:], '-',color=next(color),markersize=5,label=Label)    
       
plt.grid()
plt.legend(loc=0)    
plt.xlabel('Distance (m)')
plt.ylabel(r'Température ($^\circ C$)')
plt.title('Evolution de la température')

