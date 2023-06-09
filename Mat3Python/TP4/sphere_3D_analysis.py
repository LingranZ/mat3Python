# -------------------------------------------------------------------------
#                      Programme Python
#                 l'équation de la chaleur dans une sphere
# -------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm


def diff_analysis(n,dt,mat,volumeVersion):
    # -------------------------------------------------------------------------
    #                 Choix des paramètres
    # -------------------------------------------------------------------------
    
    #n=10 # Nombre total de couches
    #dt=0.01 # pas de temps
    #mat='Fer'
    #volumeVersion=3 # 1, 2,3
    
    R=0.010 # le rayon maximum
    
    dr=R/(n-0.5) # L'épaisseur de chaque couche ; # dr/2+(n-1)*dr=R
    #print ('dr=',dr)
    
    couche = np.linspace(0, n-1, n)    # index de chaque couche: [0,1, 2, ...,n-1]
    r_couche=dr/2 + couche*dr # rayon de chaque couche
    
    
    
    tfin=4000 # Temps final
    
    duree_chaud=tfin*0.05
    
    Nt=int(tfin/dt) # Nombre total de pas
    t = np.linspace(0, tfin, Nt)    # Maillage pour le temps
    
    
    
    
    # -------------------------------------------------------------------------
    #                 Choix des matériaux
    # -------------------------------------------------------------------------
    
    # Lambda: Conductivité thermique ( W m−1 K−1 ) 
    # cp: Capacité thermique massique (J K−1 kg−1) 
    # rho: Masse volumique ( kg/m3 ) 
    
    if mat=='Fer':
        # Fer
        cp= 444; rho= 7860; Lambda= 80; mat='Fer'
        
        # Acier	
        #cp= 435; rho= 7500; Lambda=50; mat='Acier'
        
        # Aluminium		
        #cp= 897; rho=2700; Lambda=220; mat='Aluminium'
        
        # Cuivre	
        #cp= 385; rho= 8960; Lambda=368; mat='Cuivre'
    if mat=='Bois':
        # Bois		
        cp= 2000; rho= 800; Lambda=0.04; mat='Bois'
    
    
    h= 8 # W m-2 K-1
    To=400
    Tair=20
    
    
    D= Lambda/(rho*cp) # diffusivité thermique (m2/s)
    M = h/(rho*cp)
    
    
    # Nombre de Biot
    Bi=h*(R/3)/Lambda
    Tau=rho*cp*(R/3)/h
    F=D*dt/(dr**2)
    
    # -------------------------------------------------------------------------
    #                  Les conditions initiales
    # -------------------------------------------------------------------------
    
    # Conditions initiales:
    Tini=np.ones(n)*To 
    
    # -------------------------------------------------------------------------
    #         L'evolution temporelle de la température dans la sphere
    # -------------------------------------------------------------------------
    
    Tavant=Tini
    Tapres=np.zeros(n)
    
    
    tableau_temps_Tx=np.zeros((Nt, n))
    tableau_temps_Tx[0,:]=Tini
    
    for i in range(0, Nt-1): 
        # pour la couche 0 :
        Vo=4/3*np.pi*r_couche[0]**3
        So=4*np.pi*r_couche[0]**2
        Tapres[0] = Tavant[0] + So/Vo*D*dt/dr*(Tavant[1]-Tavant[0])
        
        # pour la couche i, i=[1,2,...,n-2]
        for j in range(1, n-1):
            if volumeVersion==1:
                Vi=4*np.pi*r_couche[j-1]**2*dr
            elif volumeVersion==2:
                Vi=4*np.pi*r_couche[j]**2*dr
            else:
                Vi=4/3*np.pi*(r_couche[j]**3-r_couche[j-1]**3)
            Sinti=4*np.pi*r_couche[j-1]**2
            Sexti=4*np.pi*r_couche[j]**2
        
            
            Tapres[j] = Tavant[j] - Sinti/Vi*D*dt/dr*(Tavant[j]-Tavant[j-1]) + Sexti/Vi*D*dt/dr*(Tavant[j+1]-Tavant[j])
        
        # pour la couche n-1   
        if volumeVersion==1:
            Vfin=4*np.pi*r_couche[n-2]**2*dr
        elif volumeVersion==2:
            Vfin=4*np.pi*r_couche[n-1]**2*dr
        else:
            Vfin=4/3*np.pi*(r_couche[n-1]**3-r_couche[n-2]**3)
                
        
        Sintfin=4*np.pi*r_couche[n-2]**2
        Sextfin=4*np.pi*r_couche[n-1]**2
            
        
        Tapres[n-1] = Tavant[n-1]  -Sintfin/Vfin*D*dt/dr* (Tavant[n-1]-Tavant[n-2])- Sextfin/Vfin*M*dt *(Tavant[n-1]-Tair)   
    
        
        # Mettre à jour Tavant (avant le prochain pas de temps)
        Tavant= Tapres
        tableau_temps_Tx[i+1,:]=Tapres
     
        
    # -------------------------------------------------------------------------
    #             Traitement des résultats 
    # -------------------------------------------------------------------------
    
    #Sélectionner plusieurs moments du temps de calcul 
    t_index=np.array([0,0.1,0.2,0.5,1])*(Nt-1) 
    t_index = t_index.astype(int)
    
    
    
    #Générer plusieurs couleurs d'arc-en-ciel
    color = iter(cm.rainbow(np.linspace(0, 1,len(t_index))) )
    plt.figure(figsize=(8, 6)) 
    titre=mat+', VolumeVersion='+str(volumeVersion)+', dt='+str(dt)+', n='+str(n)+', F='+str(round(F,4))+', Bi='+str(round(Bi,5))
    
    plt.subplot(2,1,1)
    plt.title(titre)
    for k in t_index:
        # Tracé la courbe de la température, pour chaque moment sélectionné
        Label='t='+str( round(t[k],2)  )+'s'
        plt.plot(r_couche*1e3, tableau_temps_Tx[k,:], '-',color=next(color),markersize=5,label=Label,lw=2)    
           
    plt.grid(ls='--')
    plt.legend(loc=0)    
    plt.xlabel('Rayon (mm)')
    plt.ylabel(r'Température ($^\circ C$)')
    plt.axis([0,10,0,500])
    #plt.title('Evolution de la température, matériaux='+mat)
    
    
    plt.subplot(2,1,2)
    plt.plot(t, (tableau_temps_Tx[:,-1]-Tair)/(To-Tair), '-m',label='(T-Tair)/(To-Tair)',lw=2)
    plt.plot(t, np.exp(-t/Tau), '-c',label='exp(-t/tau)',lw=2) 
    
    plt.grid(ls='--')
    plt.legend(loc=0)    
    plt.xlabel('Temps (s)')
    plt.ylabel(r' ')
    plt.axis([0,4000,0,1])
    
    plt.tight_layout()
    plt.savefig(titre+'.png',dpi=50,format='png',transparent=False)
    
    A=(tableau_temps_Tx[:,-1]-Tair)/(To-Tair)
    B=np.exp(-t/Tau)
    
    diff=max(abs(A-B))
    return (F,Bi,diff)


'''

# -------------------------------------------------------------------------
#             dt change, volumeVersion=3 
# -------------------------------------------------------------------------
    
#n=10 # Nombre total de couches
#dt=0.01 # pas de temps
#mat='Fer'
#volumeVersion=3 # 1, 2,3

DT=np.linspace(0.001,0.01,10)

n=10 # Nombre total de couches
dt=0.01 # pas de temps
mat='Fer'
volumeVersion=3 # 1, 2,3

FF=[]
BBI=[]
DIFF=[]

for i in DT:
    print (i)
    out=diff_analysis(n,i,mat,volumeVersion) #(F,Bi,diff)
    FF.append(out[0])
    BBI.append(out[1])
    DIFF.append(out[2])



titre='analysis, dtChange, '+mat+', VolumeVersion='+str(volumeVersion)+', n='+str(n)


plt.figure(figsize=(8, 6)) 
plt.plot(DT,DIFF, '-sr',label='Diff')

plt.plot(DT,FF, '-*b',label='F')
plt.plot(DT,BBI, '-<c',label='Bi')

plt.xlabel('dt (s)')

plt.legend(loc=0)
plt.title(titre)

plt.tight_layout()
plt.grid(ls='--')
plt.savefig(titre+'.png',dpi=50,format='png',transparent=False)

'''


'''

# -------------------------------------------------------------------------
#             dt change, volumeVersion=1,2,3,n=10
# -------------------------------------------------------------------------
    
#n=10 # Nombre total de couches
#dt=0.01 # pas de temps
#mat='Fer'
#volumeVersion=3 # 1, 2,3

DT=np.linspace(0.001,0.01,10)

n=10 # Nombre total de couches
mat='Fer'


FF_v1=[]
BBI_v1=[]
DIFF_v1=[]


FF_v2=[]
BBI_v2=[]
DIFF_v2=[]

FF_v3=[]
BBI_v3=[]
DIFF_v3=[]

for i in DT:
    print (i)
    volumeVersion=1
    out_v1=diff_analysis(n,i,mat,volumeVersion) #(F,Bi,diff)
    FF_v1.append(out_v1[0])
    BBI_v1.append(out_v1[1])
    DIFF_v1.append(out_v1[2])
    
    volumeVersion=2
    out_v2=diff_analysis(n,i,mat,volumeVersion) #(F,Bi,diff)
    FF_v2.append(out_v2[0])
    BBI_v2.append(out_v2[1])
    DIFF_v2.append(out_v2[2])
    
    volumeVersion=3
    out_v3=diff_analysis(n,i,mat,volumeVersion) #(F,Bi,diff)
    FF_v3.append(out_v3[0])
    BBI_v3.append(out_v3[1])
    DIFF_v3.append(out_v3[2])
    
    



titre='analysis, dtChange, '+mat+', VolumeVersion=1,2,3, n='+str(n)+', Bi=0.00033'


plt.figure(figsize=(8, 6)) 
plt.plot(DT,DIFF_v1, '--sr',label='Diff_v1')
plt.plot(DT,DIFF_v2, '-og',label='Diff_v2')
plt.plot(DT,DIFF_v3, '-.<b',label='Diff_v3')

plt.plot(DT,FF_v2, '*c',label='F')


plt.xlabel('dt (s)')

plt.legend(loc=0)
plt.title(titre)

plt.tight_layout()
plt.grid(ls='--')
plt.savefig(titre+'.png',dpi=50,format='png',transparent=False)

'''




# -------------------------------------------------------------------------
#             n change, volumeVersion=1,2,3,dt=0.001
# -------------------------------------------------------------------------
    
#n=10 # Nombre total de couches
#dt=0.01 # pas de temps
#mat='Fer'
#volumeVersion=3 # 1, 2,3

dt=0.001
N=np.linspace(10,100,10)

mat='Fer'


FF_v1=[]
BBI_v1=[]
DIFF_v1=[]


FF_v2=[]
BBI_v2=[]
DIFF_v2=[]

FF_v3=[]
BBI_v3=[]
DIFF_v3=[]

for i in N:
    i=int(i)
    print (i)
    volumeVersion=1
    out_v1=diff_analysis(i,dt,mat,volumeVersion) #(F,Bi,diff)
    FF_v1.append(out_v1[0])
    BBI_v1.append(out_v1[1])
    DIFF_v1.append(out_v1[2])
    print ('volumeVersion=',volumeVersion)
    
    volumeVersion=2
    out_v2=diff_analysis(i,dt,mat,volumeVersion) #(F,Bi,diff)
    FF_v2.append(out_v2[0])
    BBI_v2.append(out_v2[1])
    DIFF_v2.append(out_v2[2])
    print ('volumeVersion=',volumeVersion)
    
    volumeVersion=3
    out_v3=diff_analysis(i,dt,mat,volumeVersion) #(F,Bi,diff)
    FF_v3.append(out_v3[0])
    BBI_v3.append(out_v3[1])
    DIFF_v3.append(out_v3[2])
    print ('volumeVersion=',volumeVersion)
    
    



titre='analysis, NChange, '+mat+', VolumeVersion=1,2,3, dt=0.001'


plt.figure(figsize=(8, 6)) 
plt.plot(N,DIFF_v1, '--sr',label='Diff_v1')
plt.plot(N,DIFF_v2, '-og',label='Diff_v2')
plt.plot(N,DIFF_v3, '-.<b',label='Diff_v3')

plt.plot(N,FF_v2, '*c',label='F')


plt.xlabel('dt (s)')

plt.legend(loc=0)
plt.title(titre)

plt.tight_layout()
plt.grid(ls='--')
plt.savefig(titre+'.png',dpi=50,format='png',transparent=False)






 
